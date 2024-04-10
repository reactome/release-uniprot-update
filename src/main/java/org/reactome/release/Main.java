package org.reactome.release;

import org.gk.model.GKInstance;
import org.gk.model.InstanceDisplayNameGenerator;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.reactome.release.reports.DuplicateAccessionReport;
import org.reactome.release.reports.Reportable;
import org.reactome.release.reports.TrEMBLAccessionReport;
import org.reactome.util.general.DBUtils;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.time.ZonedDateTime;
import java.time.format.DateTimeFormatter;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import static org.reactome.release.Utils.emptyListIfNull;
import static org.reactome.release.Utils.isTrEMBLId;
import static org.reactome.util.general.DBUtils.getCuratorDbAdaptor;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 *         Created 7/31/2023
 */
public class Main {
    private Path uniprotUpdateDirectoryPath;

    public static void main(String[] args) throws Exception {
        Main main = new Main();

        String configFilePathAsString = args.length > 0 ? args[0] : getDefaultConfigFilePath().toString();
        Properties configProperties = getConfigProperties(configFilePathAsString);

        main.run(configProperties);
    }

    @SuppressWarnings("unchecked")
    private void run(Properties configProperties) throws Exception {
        MySQLAdaptor dba = getCuratorDbAdaptor(configProperties);

        List<String> skipList = getSkipList();

        this.uniprotUpdateDirectoryPath = Paths.get(configProperties.getProperty("uniprotUpdateDirectory"));

        GKInstance uniProtReferenceDatabase = getUniProtReferenceDatabase(dba);
        GKInstance instanceEdit = getInstanceEdit(dba, "UniProt Update on " + getTodaysDate());

        Map<Integer, String> taxonIdToSpeciesName = getTaxonIdToSpeciesName();

        // Counters
        int totalNumberOfDbInstances;
        int numberOfInstancesInSwissProtFile = 0;
        int numberOfObsoleteInstancesWithNoEWAS = 0;
        int numberOfNewSwissProtInstances = 0;

        Map<Long,String> duplicateDbIdToReferenceGeneProductAccession = new HashMap<>();

        System.out.println("Populating rgp accession to db id...");
        Map<String, Long> rgpAccessionToDbId = getRGPAccessionToDbIdMap(dba);
        totalNumberOfDbInstances = rgpAccessionToDbId.size();
        System.out.println("Populating isoform accession to db id...");
        Map<String, Long> isoformAccessionToDbId = getIsoformAccessionToDbIdMap(dba);
        System.out.println("Populating rds identifier to db id...");
        Map<String, Long> rdsIdentifierToDbId = getRDSIdentifierToDbIdMap(dba);

        Map<String, List<String>> secondaryAccessionToPrimaryAccessionList = new HashMap<>();
        Map<String, String> misMatchedIsoformAccessionToRGPAccession = new HashMap<>();

        BufferedWriter sequenceReportWriter = Files.newBufferedWriter(
            getUniprotUpdateDirectoryPath().resolve("sequence_uniprot_report.txt"));
        BufferedWriter referenceDNASequenceReportWriter = Files.newBufferedWriter(
            getUniprotUpdateDirectoryPath().resolve("reference_DNA_sequence_report.txt"));

        String line;
        StringBuilder entryBuilder = new StringBuilder();

        int recordCounter = 0;

        SwissProtFileProcessor swissProtFileProcessor = new SwissProtFileProcessor(getUniprotUpdateDirectoryPath());
        BufferedReader swissProtFileReader = swissProtFileProcessor.getFileReader();
        while ((line = swissProtFileReader.readLine()) != null) {
            entryBuilder.append(line);

            if (line.contains("</entry>")) {
                String entry = entryBuilder.toString();
                entryBuilder = new StringBuilder();

                if (recordCounter % 1000 == 0) {
                    if (recordCounter != 0) {
                        dba.commit();
                        System.out.println(String.format("%d records processed and committed", recordCounter));
                    }
                    dba.startTransaction();
                }
                recordCounter += 1;

                List<String> accessions = matchMultipleValues(entry, "<accession>(.*?)</accession>");
                String primaryAccession = accessions.remove(0);
                for (String secondaryAccession : accessions) {
                    secondaryAccessionToPrimaryAccessionList.computeIfAbsent(
                        secondaryAccession, k -> new ArrayList<>()).add(primaryAccession);
                }

                String organismName = matchSingleValue(entry, "<name type=\"scientific\">(.*?)</name>");
                String taxon = "";
                Map<String, GKInstance> speciesNameToInstanceCache = new HashMap<>();
                GKInstance speciesInstance = null;
                for (String speciesName : taxonIdToSpeciesName.values()) {
                    if (organismName.contains(speciesName)) {
                        taxon = speciesName;
                        speciesInstance = getSpeciesInstance(dba, taxon, speciesNameToInstanceCache);
                    }
                }

                if (taxon.length() < 2 && !rgpAccessionToDbId.containsKey(primaryAccession)) {
                    continue;
                }

                numberOfInstancesInSwissProtFile += 1;

                String id = matchSingleValue(entry, "<name>([A-Za-z0-9_]*)</name>");

                accessions.add(0, id);
                String description = matchSingleValue(entry, "<protein(.*)</protein>");

                String fullName = matchSingleValue(description, "<recommendedName>\\s+<fullName>(.*?)</fullName>");
                if (fullName.isEmpty()) {
                    fullName = matchSingleValue(description,
                        "<recommendedName ref=\"\\d+\">\\s+<fullName>(.*)</fullName>");
                }
                String recommendedName = !fullName.isEmpty() ? fullName : "No name";

                description =
                    description
                    .replaceAll("</fullName>","")
                    .replaceAll("<fullName>","")
                    .replaceAll("</recommendedName>","")
                    .replaceAll("<recommendedName>"," recommendedName: ")
                    .replaceAll("</alternativeName>","")
                    .replaceAll("</shortName>","")
                    .replaceAll("<alternativeName>"," alternativeName: ")
                    .replaceAll("<shortName>"," shortName: ")
                    .replaceAll("<recommendedName ref=\"\\d+\"","")
                    .replaceAll("<ecNumber>", " ecNumber: ")
                    .replaceAll("</ecNumber>", "")
                    .replaceAll(" +"," ")
                    .replaceAll("\\n","\t")
                    .replaceAll("\\t"," ")
                    .replaceAll(">","")
                    .replaceAll("<","")
                    .replaceAll("type=\"fragments?\"","")
                    .replaceAll("^\\s+","")
                    .replaceAll("\\s+$","");

                Integer sequenceLength = Integer.parseInt(
                    matchSingleValue(entry, "<sequence.*length=\"(\\d+)\""));

                String checksum = matchSingleValue(entry, "<sequence.*checksum=\"([0-9A-F]+)\"");

                List<String> geneNames = matchMultipleValues(entry, "<gene>(.*?)</gene>").stream().flatMap(
                    names -> Arrays.stream(names.trim().split("\\s{2,}")).map(geneName ->
                        geneName
                        .replaceAll("</name>","")
                        .replaceAll("<name.*?>","")
                        .replaceAll(" {2}", "")
                    )
                ).distinct().collect(Collectors.toList());


                String name = !geneNames.isEmpty() &&!geneNames.get(0).isEmpty() ?
                    geneNames.get(0) :
                    recommendedName;

                List<GKInstance> referenceDNASequences = new ArrayList<>();
                if (taxon.contains("Homo sapiens")) {
                    String typeValueRegex = "<property type=\"gene ID\" value=\"(ENSG.*?)\"";
                    Set<String> uniqueEnsEMBLGeneIds = new HashSet<>();
                    uniqueEnsEMBLGeneIds.addAll(matchMultipleValues(entry, typeValueRegex));
                    uniqueEnsEMBLGeneIds = uniqueEnsEMBLGeneIds
                        .stream()
                        .distinct()
                        .map(this::removeVersionNumber)
                        .collect(Collectors.toSet());

                    if (uniqueEnsEMBLGeneIds.size() > 1) {
                        referenceDNASequenceReportWriter.write("Multiple gene ids -- " +
                            String.join("\t", primaryAccession, name, uniqueEnsEMBLGeneIds.toString()) + "\n");
                    }
                    GKInstance humanEnsEMBLGeneReferenceDatabase = getHumanEnsEMBLGeneReferenceDatabase(dba);

                    for (String ensEMBLGeneId : uniqueEnsEMBLGeneIds) {
                        GKInstance referenceDNASequence;

                        if (rdsIdentifierToDbId.containsKey(ensEMBLGeneId)) {
                            referenceDNASequenceReportWriter.write("Checking existing reference DNA sequence for " +
                                ensEMBLGeneId + " with db_id " + rdsIdentifierToDbId.get(ensEMBLGeneId) + "\n");

                            long rdsDbId = rdsIdentifierToDbId.get(ensEMBLGeneId);
                            referenceDNASequence = fetchReferenceDNASequenceByDbId(dba, rdsDbId);

                            GKInstance existingRDSReferenceDatabase = (GKInstance)
                                referenceDNASequence.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
                            boolean isUpdateToReferenceDNASequence = false;
                            if (existingRDSReferenceDatabase == null ||
                                !sameDbId(existingRDSReferenceDatabase, humanEnsEMBLGeneReferenceDatabase)) {
                                referenceDNASequence.addAttributeValue(
                                    ReactomeJavaConstants.referenceDatabase, humanEnsEMBLGeneReferenceDatabase);
                                isUpdateToReferenceDNASequence = true;
                            }

                            List<String> existingGeneNames = (List<String>)
                                referenceDNASequence.getAttributeValuesList(ReactomeJavaConstants.geneName);
                            if (existingGeneNames == null || areDifferentLists(existingGeneNames, geneNames)) {
                                referenceDNASequence.setAttributeValue(ReactomeJavaConstants.geneName, geneNames);
                                isUpdateToReferenceDNASequence = true;
                            }
                            GKInstance existingSpeciesInstance = (GKInstance)
                                referenceDNASequence.getAttributeValue(ReactomeJavaConstants.species);
                            if (existingSpeciesInstance == null ||
                                (speciesInstance != null &&
                                !existingSpeciesInstance.getDBID().equals(speciesInstance.getDBID()))) {

                                referenceDNASequence.addAttributeValue(ReactomeJavaConstants.species, speciesInstance);
                                isUpdateToReferenceDNASequence = true;
                            }

                            String existingIdentifier = (String)
                                referenceDNASequence.getAttributeValue(ReactomeJavaConstants.identifier);
                            if (existingIdentifier == null || !existingIdentifier.equals(ensEMBLGeneId)) {
                                referenceDNASequence.setAttributeValue(
                                    ReactomeJavaConstants.identifier, ensEMBLGeneId);
                                isUpdateToReferenceDNASequence = true;
                            }

                            if (isUpdateToReferenceDNASequence) {
                                referenceDNASequenceReportWriter.write(
                                    "Updating existing reference DNA sequence for " + ensEMBLGeneId + " with db_id " +
                                    rdsIdentifierToDbId.get(ensEMBLGeneId) + "\n"
                                );
                                dba.loadInstanceAttributeValues(referenceDNASequence);
                                InstanceDisplayNameGenerator.setDisplayName(referenceDNASequence);
                                referenceDNASequence.addAttributeValue(ReactomeJavaConstants.modified, instanceEdit);
                                dba.updateInstance(referenceDNASequence);
                            }
                        } else {
                            if (uniqueEnsEMBLGeneIds.size() > 1 && !onEnsEMBLPrimaryAssembly(ensEMBLGeneId)) {
                                // Reference DNA Sequences to be created only for primary gene ids for a UniProt entry
                                // When there is only one gene id for a UniProt entry, it is assumed to be the primary
                                // id
                                referenceDNASequenceReportWriter.write(ensEMBLGeneId + " is not a primary/canonical " +
                                    "gene -- skipping creation of ReferenceDNASequence\n"
                                );
                                continue;
                            }

                            referenceDNASequence = new GKInstance(
                                dba.getSchema().getClassByName(ReactomeJavaConstants.ReferenceDNASequence));
                            referenceDNASequence.setDbAdaptor(dba);
                            referenceDNASequence.setAttributeValue(
                                ReactomeJavaConstants.referenceDatabase, humanEnsEMBLGeneReferenceDatabase);
                            referenceDNASequence.setAttributeValue(
                                ReactomeJavaConstants.identifier, ensEMBLGeneId
                            );
                            referenceDNASequence.setAttributeValue(
                                ReactomeJavaConstants.created, instanceEdit
                            );
                            referenceDNASequence.setAttributeValue(
                                ReactomeJavaConstants.geneName, geneNames
                            );
                            referenceDNASequence.setAttributeValue(
                                ReactomeJavaConstants.species, speciesInstance
                            );
                            InstanceDisplayNameGenerator.setDisplayName(referenceDNASequence);

                            long referenceDNASequenceDbId = dba.storeInstance(referenceDNASequence);
                            referenceDNASequenceReportWriter.write("Reference DNA sequence with db_id " +
                                referenceDNASequenceDbId + " created for " + ensEMBLGeneId + "\n");
                            rdsIdentifierToDbId.put(ensEMBLGeneId, referenceDNASequenceDbId);
                        }
                        referenceDNASequences.add(referenceDNASequence);
                    }
                }
                List<String> keywords = matchMultipleValues(entry, "<keyword id=\".*?\">(.*?)</keyword>");

                String comments = parseComments(entry);

                List<String> isoformIds = matchMultipleValues(entry, "<isoform>\\s*<id>([A-Z0-9-]*)");

                List<String> chains = parseChains(entry);

                Map<String,List<?>> values = new HashMap<>();
                values.put(ReactomeJavaConstants.secondaryIdentifier, accessions);
                values.put(ReactomeJavaConstants.description, Collections.singletonList(description));
                values.put(ReactomeJavaConstants.sequenceLength, Collections.singletonList(sequenceLength));
                values.put(ReactomeJavaConstants.species, Collections.singletonList(speciesInstance));
                values.put(ReactomeJavaConstants.checksum, Collections.singletonList(checksum));
                values.put(ReactomeJavaConstants.name, Collections.singletonList(name));
                values.put(ReactomeJavaConstants.geneName, geneNames);
                values.put(ReactomeJavaConstants.comment, Collections.singletonList(comments));
                values.put(ReactomeJavaConstants.keyword, keywords);
                values.put(ReactomeJavaConstants.chain, chains);
                if (taxon.contains("Homo sapiens")) {
                    values.put(ReactomeJavaConstants.referenceGene, referenceDNASequences);
                }
                if (!rgpAccessionToDbId.containsKey(primaryAccession)) {
                    numberOfNewSwissProtInstances += 1;

                    GKInstance newReferenceGeneProductInstance =
                        new GKInstance(dba.getSchema().getClassByName(ReactomeJavaConstants.ReferenceGeneProduct));
                    newReferenceGeneProductInstance.setDbAdaptor(dba);
                    newReferenceGeneProductInstance.setAttributeValue(
                        ReactomeJavaConstants.referenceDatabase, uniProtReferenceDatabase);
                    newReferenceGeneProductInstance.setAttributeValue(
                        ReactomeJavaConstants.identifier, primaryAccession);
                    newReferenceGeneProductInstance.setAttributeValue(ReactomeJavaConstants.created, instanceEdit);
                    long newRGPDbId = dba.storeInstance(newReferenceGeneProductInstance);

                    System.out.println(String.format("New UniProt:%s\t%d", primaryAccession, newRGPDbId));
                    updateInstance(newReferenceGeneProductInstance, values, sequenceReportWriter);
                    for (String isoformId : isoformIds) {
                        if (!isoformId.contains(primaryAccession)) {
                            misMatchedIsoformAccessionToRGPAccession.put(isoformId, primaryAccession);
                        }

                        GKInstance newIsoformInstance = new GKInstance(
                            dba.getSchema().getClassByName(ReactomeJavaConstants.ReferenceIsoform));
                        newIsoformInstance.setDbAdaptor(dba);
                        newIsoformInstance.setAttributeValue(
                            ReactomeJavaConstants.referenceDatabase, uniProtReferenceDatabase);
                        newIsoformInstance.setAttributeValue(ReactomeJavaConstants.identifier, primaryAccession);
                        newIsoformInstance.setAttributeValue(
                            ReactomeJavaConstants.isoformParent, newReferenceGeneProductInstance);
                        newIsoformInstance.setAttributeValue(ReactomeJavaConstants.created, instanceEdit);
                        newIsoformInstance.setAttributeValue(ReactomeJavaConstants.variantIdentifier, isoformId);

                        updateInstance(newIsoformInstance, values, sequenceReportWriter);
                    }
                } else {
                    Collection<GKInstance> existingReferenceGeneProductInstances = dba.fetchInstanceByAttribute(
                        ReactomeJavaConstants.ReferenceGeneProduct,
                        ReactomeJavaConstants.identifier,
                        "=",
                        primaryAccession
                    );
                    boolean duplicateFlag = false;
                    for (GKInstance existingReferenceGeneProductInstance : existingReferenceGeneProductInstances) {
                        if (isAReferenceIsoform(existingReferenceGeneProductInstance)) {
                            continue;
                        }

                        if (duplicateFlag) {
                            duplicateDbIdToReferenceGeneProductAccession.put(
                                existingReferenceGeneProductInstance.getDBID(), primaryAccession);
                            continue;
                        }

                        System.out.println(String.format("Updating master sequence...%d\t%s",
                            existingReferenceGeneProductInstance.getDBID(), primaryAccession));

                        dba.loadInstanceAttributeValues(existingReferenceGeneProductInstance);
                        existingReferenceGeneProductInstance.addAttributeValue(
                            ReactomeJavaConstants.modified, instanceEdit);

                        updateInstance(existingReferenceGeneProductInstance, values, sequenceReportWriter);

                        duplicateFlag = true;

                        if (values.get(ReactomeJavaConstants.species).isEmpty()) {
                            values.put(ReactomeJavaConstants.species, Collections.singletonList((GKInstance)
                                existingReferenceGeneProductInstance.getAttributeValue(ReactomeJavaConstants.species))
                            );
                        }
                        for (String isoformId : isoformIds) {
                            if (isoformId.contains(primaryAccession)) {
                                Collection<GKInstance> isoformInstances = dba.fetchInstanceByAttribute(
                                    ReactomeJavaConstants.ReferenceIsoform,
                                    ReactomeJavaConstants.variantIdentifier,
                                    "=",
                                    isoformId
                                );
                                if (!isoformInstances.isEmpty()) {
                                    for (GKInstance isoformInstance : isoformInstances) {
                                        String isoformAccession = (String) isoformInstance.getAttributeValue(
                                            ReactomeJavaConstants.variantIdentifier);
                                        if (!isoformAccession.contains(primaryAccession)) {
                                            continue;
                                        }
                                        System.out.println(String.format("Existing isoform update: %s\tMaster: %d",
                                            isoformAccession, existingReferenceGeneProductInstance.getDBID()));

                                        dba.loadInstanceAttributeValues(isoformInstance);
                                        isoformInstance.setAttributeValue(ReactomeJavaConstants.isoformParent,
                                            existingReferenceGeneProductInstance);
                                        isoformInstance.addAttributeValue(ReactomeJavaConstants.modified,
                                            instanceEdit);

                                        updateInstance(isoformInstance, values, sequenceReportWriter);

                                        isoformAccessionToDbId.remove(isoformId);
                                    }
                                } else {
                                    GKInstance isoformInstance = new GKInstance(
                                        dba.getSchema().getClassByName(ReactomeJavaConstants.ReferenceIsoform)
                                    );
                                    isoformInstance.setDbAdaptor(dba);
                                    isoformInstance.setAttributeValue(ReactomeJavaConstants.identifier,
                                        primaryAccession);
                                    isoformInstance.setAttributeValue(ReactomeJavaConstants.isoformParent,
                                        existingReferenceGeneProductInstance);
                                    isoformInstance.setAttributeValue(ReactomeJavaConstants.created,
                                        instanceEdit);
                                    isoformInstance.setAttributeValue(ReactomeJavaConstants.variantIdentifier,
                                        isoformId);
                                    long isoformDbId = dba.storeInstance(isoformInstance);

                                    System.out.println(String.format("New isoform: %s\t%d\tMaster: %d",
                                        isoformId, isoformDbId, existingReferenceGeneProductInstance.getDBID()));

                                    updateInstance(isoformInstance, values, sequenceReportWriter);
                                }
                            } else {
                                misMatchedIsoformAccessionToRGPAccession.put(isoformId, primaryAccession);
                            }
                        }
                        rgpAccessionToDbId.remove(primaryAccession);
                    }
                }
            }
        }

        referenceDNASequenceReportWriter.close();
        sequenceReportWriter.close();

        System.out.println(recordCounter + " records processed and committed");
        System.out.println("All records in " + swissProtFileProcessor.getSwissProtFilePath() + " processed");

        System.out.println("Starting clean-up tasks after processing UniProt XML");
        dba.startTransaction();

        System.out.println("Updating mis-matched isoforms");

        for (String misMatchedIsoformAccession : misMatchedIsoformAccessionToRGPAccession.keySet()) {
            List<GKInstance> isoformParents = new ArrayList<>();

            Iterator<GKInstance> isoformInstanceIterator = (dba.fetchInstanceByAttribute(
                ReactomeJavaConstants.ReferenceIsoform,
                ReactomeJavaConstants.variantIdentifier,
                "=",
                misMatchedIsoformAccession
            )).iterator();

            GKInstance isoformInstance = isoformInstanceIterator.hasNext() ? isoformInstanceIterator.next() : null;
            if (isoformInstance != null) {
                GKInstance isoformParent =
                    (GKInstance) isoformInstance.getAttributeValue(ReactomeJavaConstants.isoformParent);
                if (isoformParent == null) {
                    continue;
                }
                isoformParents.add(isoformParent);
            }

            Iterator<GKInstance> mismatchedParentIterator = (dba.fetchInstanceByAttribute(
                ReactomeJavaConstants.ReferenceGeneProduct,
                ReactomeJavaConstants.identifier,
                "=",
                misMatchedIsoformAccessionToRGPAccession.get(misMatchedIsoformAccession)
            )).iterator();

            GKInstance mismatchedParent = mismatchedParentIterator.hasNext() ? mismatchedParentIterator.next() : null;
            if (mismatchedParent != null && isoformInstance != null) {
                isoformParents.add(mismatchedParent);
                long isoformInstanceDbId = isoformInstance.getDBID();
                System.out.println(String.format("Mismatched parent: %s(%d)\t%s\n",
                    misMatchedIsoformAccession,
                    isoformInstanceDbId,
                    misMatchedIsoformAccessionToRGPAccession.get(misMatchedIsoformAccession))
                );

                isoformInstance.setAttributeValue(ReactomeJavaConstants.isoformParent, isoformParents);
                dba.updateInstance(isoformInstance);
            }
        }

        System.out.println("Mis-matched isoform updates complete");

        System.out.println("Updating display names...");

        updateDisplayNames(dba, ReactomeJavaConstants.ReferenceGeneProduct);
        updateDisplayNames(dba, ReactomeJavaConstants.ReferenceIsoform);

        System.out.println("Done");

        System.out.println("Remaining instances:" + rgpAccessionToDbId.keySet().size());

        System.out.println("Deleting obsolete instances with no referrers...");

        Iterator<String> rgpAccessionsIterator = rgpAccessionToDbId.keySet().iterator();
        List<String> tremblAccessions = new ArrayList<>();
        while (rgpAccessionsIterator.hasNext()) {
            String rgpAccession = rgpAccessionsIterator.next();

            if (isTrEMBLId(rgpAccession)) {
                tremblAccessions.add(rgpAccession);
                rgpAccessionsIterator.remove();
            } else {
                Collection<GKInstance> obsoleteReferenceGeneProductInstances = dba.fetchInstanceByAttribute(
                    ReactomeJavaConstants.ReferenceGeneProduct,
                    ReactomeJavaConstants.identifier,
                    "=",
                    rgpAccession
                );

                boolean isObsoleteRGPDeleted = false;
                for (GKInstance obsoleteReferenceGeneProductInstance : obsoleteReferenceGeneProductInstances) {
                    String variantIdentifier = null;
                    if (isAReferenceIsoform(obsoleteReferenceGeneProductInstance)) {
                        variantIdentifier = (String) obsoleteReferenceGeneProductInstance.getAttributeValue(
                            ReactomeJavaConstants.variantIdentifier);
                    }

                    if (variantIdentifier != null) {
                        continue;
                    }

                    long obsoleteRGPDbId = obsoleteReferenceGeneProductInstance.getDBID();
                    Collection<GKInstance> referrers = getRGPReferrers(obsoleteReferenceGeneProductInstance);
                    if (referrers == null || referrers.isEmpty()) {
                        System.out.println("Deleting " + obsoleteRGPDbId + "...");
                        dba.deleteByDBID(obsoleteRGPDbId);
                        numberOfObsoleteInstancesWithNoEWAS += 1;
                        isObsoleteRGPDeleted = true;
                    }
                }
                if (isObsoleteRGPDeleted) {
                    rgpAccessionsIterator.remove();
                }
            }
        }
        Reportable trEMBLAccessionReport = new TrEMBLAccessionReport(getUniprotUpdateDirectoryPath(), tremblAccessions);
        trEMBLAccessionReport.writeReport();

        List<Long> dbIdsToSkip = new ArrayList<>();
        Iterator<String> isoformAccessionIterator = isoformAccessionToDbId.keySet().iterator();
        while (isoformAccessionIterator.hasNext()) {
            String isoformAccession = isoformAccessionIterator.next();
            Iterator<GKInstance> isoformInstanceIterator = (dba.fetchInstanceByAttribute(
                ReactomeJavaConstants.ReferenceIsoform,
                ReactomeJavaConstants.variantIdentifier,
                "=",
                isoformAccession
            )).iterator();

            GKInstance isoformInstance = isoformInstanceIterator.hasNext() ? isoformInstanceIterator.next() : null;
            if (isoformInstance == null) {
                System.out.println(isoformAccession + " is not a variant identifier for any ReferenceIsoform");
                continue;
            }

            long obsoleteIsoformDbId = isoformInstance.getDBID();
            GKInstance isoformParent =
                (GKInstance) isoformInstance.getAttributeValue(ReactomeJavaConstants.isoformParent);
            if (isoformParent == null) {
                System.out.println(isoformInstance.getDBID());
                dbIdsToSkip.add(obsoleteIsoformDbId);
                continue;
            }

            String isoformParentIdentifier = (String)
                isoformParent.getAttributeValue(ReactomeJavaConstants.identifier);
            if (isoformParentIdentifier == null || isoformParentIdentifier.isEmpty()) {
                continue;
            }

            Collection<GKInstance> referrers = getRGPReferrers(isoformInstance);
            if (referrers == null || referrers.isEmpty()) {
                System.out.println("Deleting " + obsoleteIsoformDbId + "...");
                dba.deleteByDBID(obsoleteIsoformDbId);
                numberOfObsoleteInstancesWithNoEWAS += 1;
                isoformAccessionIterator.remove();
            }
        }
        System.out.println("Done.");

        System.out.println("Preparing reports...");
        Set<Long> noReferrerDbIds = new HashSet<>();

        Reportable duplicateAccessionReport = new DuplicateAccessionReport(
            getUniprotUpdateDirectoryPath(), duplicateDbIdToReferenceGeneProductAccession);
        duplicateAccessionReport.writeReport();

        List<String> skipReplaceableReportLines = new ArrayList<>();
        List<String> skipNoReplacementReportLines = new ArrayList<>();

        BufferedWriter wikiWriter = Files.newBufferedWriter(getUniprotUpdateDirectoryPath().resolve("uniprot.wiki"));

        wikiWriter.write(
        "{| class=\"wikitable\"\n" +
            "|+ Obsolete UniProt Instances (with replacement UniProt)\n" +
            "|-\n" +
            "! Replacement UniProt\n" +
            "! Obsolete UniProt\n" +
            "! Reactome instances with obsolete UniProt\n" +
            "! EWAS associated with obsolete UniProt\n" +
            "! Species\n" +
            "|-\n"
        );

        rgpAccessionsIterator = rgpAccessionToDbId.keySet().iterator();
        while (rgpAccessionsIterator.hasNext()) {
            String rgpAccession = rgpAccessionsIterator.next();

            boolean isSecondaryAccession = false;

            if (secondaryAccessionToPrimaryAccessionList.containsKey(rgpAccession)) {
                List<String> alternateAccessions = secondaryAccessionToPrimaryAccessionList.get(rgpAccession);
                if (alternateAccessions == null) {
                    System.err.println("Zero alternate accessions for " + rgpAccession + ": " +
                        alternateAccessions);
                    continue;
                }
                isSecondaryAccession = true;

                long obsoleteDbId;
                List<Long> referrerDbIds = new ArrayList<>();
                String speciesName;

                Collection<GKInstance> obsoleteRGPInstances = dba.fetchInstanceByAttribute(
                    ReactomeJavaConstants.ReferenceGeneProduct,
                    ReactomeJavaConstants.identifier,
                    "=",
                    rgpAccession
                );
                for (GKInstance obsoleteRGPInstance : obsoleteRGPInstances) {
                    String variantIdentifier = null;
                    if (isAReferenceIsoform(obsoleteRGPInstance)) {
                        variantIdentifier = (String) obsoleteRGPInstance.getAttributeValue(
                            ReactomeJavaConstants.variantIdentifier);
                    }

                    if (variantIdentifier != null) {
                        continue;
                    }
                    obsoleteDbId = obsoleteRGPInstance.getDBID();
                    speciesName = getSpeciesName(obsoleteRGPInstance);

                    List<GKInstance> referrers = emptyListIfNull((List<GKInstance>) obsoleteRGPInstance.getReferers(
                        ReactomeJavaConstants.referenceEntity));
                    for (GKInstance referrer : referrers) {
                        if (referrer.getSchemClass().isa(ReactomeJavaConstants.EntityWithAccessionedSequence)) {
                            referrerDbIds.add(referrer.getDBID());
                        }
                    }

                    if (!referrerDbIds.isEmpty()) {
                        StringBuilder reportLineBuilder = new StringBuilder();
                        reportLineBuilder.append("||");
                        reportLineBuilder.append(String.join("|",
                            alternateAccessions
                                .stream()
                                .map(
                                    alternateAccession ->
                                    String.format("[https://www.uniprot.org/uniprot/%s %s]",
                                        alternateAccession, alternateAccession))
                                .collect(Collectors.toList())
                          ));
                        reportLineBuilder.append("\n");
                        for (String alternateAccession : alternateAccessions) {
                            System.out.println(String.format("%s\t%s\t%s",
                                rgpAccession, alternateAccession, obsoleteDbId));
                        }
                        reportLineBuilder.append(String.format("|%s\n", rgpAccession));
                        reportLineBuilder.append(String.format(
                            "|[https://curator.reactome.org/cgi-bin/instancebrowser?DB=%s&ID=%d& %d]\n",
                            dba.getDBName(), obsoleteDbId, obsoleteDbId
                            ));
                        reportLineBuilder.append(String.format("||%s\n", String.join(
                            "|", referrerDbIds.stream().map(Object::toString).collect(Collectors.toList())
                        )));
                        reportLineBuilder.append("|" + speciesName + "\n");
                        reportLineBuilder.append("|-\n");

                        String reportLine = reportLineBuilder.toString();
                        if (skipList.stream().anyMatch(accession -> accession.equals(rgpAccession))) {
                            skipReplaceableReportLines.add(reportLine);
                        } else {
                            wikiWriter.write(reportLine);
                        }
                    } else {
                        noReferrerDbIds.add(obsoleteDbId);
                    }
                }
            }
            if (isSecondaryAccession) {
                rgpAccessionsIterator.remove();
            }
        }

        wikiWriter.write("|}\n\n-----\n");

        wikiWriter.write(
            "{| class=\"wikitable\"\n" +
                "|+ Obsolete UniProt Instances (deleted forever, no replacement)\n" +
                "|-\n" +
                "! Obsolete UniProt\n" +
                "! Reactome instances with obsolete UniProt\n" +
                "! EWAS associated with obsolete UniProt\n" +
                "! Species\n" +
                "|-\n"
        );

        for (String rgpAccession : rgpAccessionToDbId.keySet()) {
            long obsoleteDbId = -1L;
            List<String> referrerIds = new ArrayList<>();
            String speciesName = "";

            Collection<GKInstance> obsoleteRGPInstances = emptyListIfNull(dba.fetchInstanceByAttribute(
                ReactomeJavaConstants.ReferenceGeneProduct,
                ReactomeJavaConstants.identifier,
                "=",
                rgpAccession
            ));

            for (GKInstance obsoleteRGPInstance : obsoleteRGPInstances) {
                String variantIdentifier = null;
                if (isAReferenceIsoform(obsoleteRGPInstance)) {
                    variantIdentifier = (String) obsoleteRGPInstance.getAttributeValue(
                        ReactomeJavaConstants.variantIdentifier);
                }

                if (variantIdentifier != null) {
                    continue;
                }
                obsoleteDbId = obsoleteRGPInstance.getDBID();
                speciesName = getSpeciesName(obsoleteRGPInstance);

                List<GKInstance> referrers = emptyListIfNull(
                    (List<GKInstance>) obsoleteRGPInstance.getReferers(ReactomeJavaConstants.referenceEntity));
                for (GKInstance referrer : referrers) {
                    if (referrer.getSchemClass().isa(ReactomeJavaConstants.EntityWithAccessionedSequence)) {
                        GKInstance referrerStableIdInstance =
                            (GKInstance) referrer.getAttributeValue(ReactomeJavaConstants.stableIdentifier);
                        if (referrerStableIdInstance != null) {
                            String referrerStableId =
                                (String) referrerStableIdInstance.getAttributeValue(ReactomeJavaConstants.identifier);
                            referrerIds.add(referrerStableId);
                        } else {
                            referrerIds.add(referrer.getDBID().toString());
                        }
                    }
                }
            }

            System.out.println(rgpAccession);
            if (!referrerIds.isEmpty()) {
                StringBuilder reportLineBuilder = new StringBuilder();
                //reportLineBuilder.append("|\n");
                reportLineBuilder.append(String.format("||%s\n", rgpAccession));
                reportLineBuilder.append(String.format(
                    "|[https://curator.reactome.org/cgi-bin/instancebrowser?DB=%s&ID=%d& %d]\n",
                    dba.getDBName(), obsoleteDbId, obsoleteDbId
                ));
                reportLineBuilder.append(String.format("||%s\n", String.join(
                    "|", referrerIds.stream().map(Object::toString).collect(Collectors.toList())
                )));
                reportLineBuilder.append(String.format("|%s\n", speciesName));
                reportLineBuilder.append("|-\n");

                String reportLine = reportLineBuilder.toString();
                if (skipList.stream().anyMatch(accession -> accession.equals(rgpAccession))) {
                    skipNoReplacementReportLines.add(reportLine);
                } else {
                    wikiWriter.write(reportLine);
                }
            } else {
                noReferrerDbIds.add(obsoleteDbId);
            }
        }

        for (String isoformAccession : isoformAccessionToDbId.keySet()) {
            Collection<GKInstance> isoformInstances = emptyListIfNull(dba.fetchInstanceByAttribute(
                ReactomeJavaConstants.ReferenceIsoform,
                ReactomeJavaConstants.variantIdentifier,
                "=",
                isoformAccession
            ));
            String speciesName;
            for (GKInstance isoformInstance : isoformInstances) {
                List<String> referrerIds = new ArrayList<>();
                long isoformInstanceDbId = isoformInstance.getDBID();
                speciesName = getSpeciesName(isoformInstance);

                List<GKInstance> referrers = emptyListIfNull(
                    (List<GKInstance>) isoformInstance.getReferers(ReactomeJavaConstants.referenceEntity));
                for (GKInstance referrer : referrers) {
                    if (referrer.getSchemClass().isa(ReactomeJavaConstants.EntityWithAccessionedSequence)) {
                        GKInstance referrerStableIdInstance =
                            (GKInstance) referrer.getAttributeValue(ReactomeJavaConstants.stableIdentifier);
                        if (referrerStableIdInstance != null) {
                            String referrerStableId =
                                (String) referrerStableIdInstance.getAttributeValue(ReactomeJavaConstants.identifier);
                            referrerIds.add(referrerStableId);
                        } else {
                            referrerIds.add(referrer.getDBID().toString());
                        }
                    }
                }

                if (!referrerIds.isEmpty()) {
                    StringBuilder reportLineBuilder = new StringBuilder();
                    //reportLineBuilder.append("|\n");
                    reportLineBuilder.append(String.format("||%s\n", isoformAccession));
                    reportLineBuilder.append(String.format(
                        "|[https://curator.reactome.org/cgi-bin/instancebrowser?DB=%s&ID=%d& %d]\n",
                        dba.getDBName(), isoformInstanceDbId, isoformInstanceDbId
                    ));
                    reportLineBuilder.append(String.format("||%s\n", String.join(
                        "|", referrerIds.stream().map(Object::toString).collect(Collectors.toList())
                    )));
                    reportLineBuilder.append(String.format("|%s\n", speciesName));
                    reportLineBuilder.append("|-\n");

                    String reportLine = reportLineBuilder.toString();
                    if (skipList.stream().anyMatch(accession -> accession.equals(isoformAccession))) {
                        skipNoReplacementReportLines.add(reportLine);
                    } else {
                        wikiWriter.write(reportLine);
                    }
                } else {
                    noReferrerDbIds.add(isoformInstanceDbId);
                }
            }
        }

        wikiWriter.write("|}\n-----\n");

        wikiWriter.write(
            "{| class=\"wikitable\"\n" +
                "|+ SKIPLIST Obsolete UniProt Instances (with replacement UniProt)\n" +
                "|-\n" +
                "! Replacement UniProt\n" +
                "! Obsolete UniProt\n" +
                "! Reactome instances with obsolete UniProt\n" +
                "! EWAS associated with obsolete UniProt\n" +
                "! Species\n" +
                "|-\n"
        );

        for (String skipReplaceableReportLine : skipReplaceableReportLines) {
            wikiWriter.append(skipReplaceableReportLine);
        }

        wikiWriter.write("|}\n-----\n");

        wikiWriter.write(
            "{| class=\"wikitable\"\n" +
                "|+ SKIPLIST Obsolete UniProt Instances (deleted forever, no replacement)\n" +
                "|-\n" +
                "! Obsolete UniProt\n" +
                "! Reactome instances with obsolete UniProt\n" +
                "! EWAS associated with obsolete UniProt\n" +
                "! Species\n" +
                "|-\n"
        );

        for (String skipNoReplacementLine : skipNoReplacementReportLines) {
            wikiWriter.append(skipNoReplacementLine);
        }
        wikiWriter.write("|}\n");
        wikiWriter.close();

        System.out.println("\nDeleting DBID with obsolete UniProt and no referrers (2nd round during wiki report)...");

        NEXT:for (long noReferrerDbId : noReferrerDbIds) {
            for (long dbIdToSkip : dbIdsToSkip) {
                if (noReferrerDbId == dbIdToSkip) {
                    continue NEXT;
                }

                dba.deleteByDBID(noReferrerDbId);
                System.out.println("Deleting DBID: " + noReferrerDbId);
            }
        }

        System.out.println("Checking for duplicate isoform instances...");

        List<GKInstance> referenceIsoformUniProtInstances = emptyListIfNull((List<GKInstance>) (
            dba.fetchInstancesByClass(ReactomeJavaConstants.ReferenceIsoform))
            .stream()
            .filter(isoform -> hasUniProtReferenceDatabase((GKInstance) isoform))
            .collect(Collectors.toList()));

        Map<String,List<Long>> variantIdentifierToDbId = new HashMap<>();
        for (GKInstance referenceIsoformUniProtInstance : referenceIsoformUniProtInstances) {
            String variantIdentifier =
                (String) referenceIsoformUniProtInstance.getAttributeValue(ReactomeJavaConstants.variantIdentifier);
            long isoformDbId = referenceIsoformUniProtInstance.getDBID();

            if (variantIdentifier == null || variantIdentifier.isEmpty()) {
                System.out.println(String.format("ReferenceIsoform %s has no variant identifier", isoformDbId));
                continue;
            }

            if (variantIdentifierToDbId.containsKey(variantIdentifier)) {
                variantIdentifierToDbId.get(variantIdentifier).add(isoformDbId);
                System.out.println(String.format("Multiple instance for %s:\t%s",
                    variantIdentifier, variantIdentifierToDbId.get(variantIdentifier)));
            } else {
                variantIdentifierToDbId.computeIfAbsent(variantIdentifier, k -> new ArrayList<>()).add(isoformDbId);
            }
        }

        dba.commit();
        System.out.println("UniProt Update has completed");
        System.out.println("Total db instances: " + totalNumberOfDbInstances);
        System.out.println("Total SwissProt instances in file: " + numberOfInstancesInSwissProtFile);
        System.out.println("Obsolete instances with no referrers: " + numberOfObsoleteInstancesWithNoEWAS);
        System.out.println("Number of new SwissProt instances: " + numberOfNewSwissProtInstances);
    }

    private static Path getDefaultConfigFilePath() throws URISyntaxException {
        return Paths.get(Main.class.getClassLoader().getResource("config.properties").toURI());
    }

    private static Properties getConfigProperties(String configFilePathAsString) throws IOException {
        Properties configProperties = new Properties();
        configProperties.load(Files.newInputStream(Paths.get(configFilePathAsString)));
        return configProperties;
    }

    private Path getUniprotUpdateDirectoryPath() {
        return this.uniprotUpdateDirectoryPath;
    }

    @SuppressWarnings("unchecked")
    private GKInstance getUniProtReferenceDatabase(MySQLAdaptor dba) throws Exception {
        Collection<GKInstance> uniProtReferenceDatabaseInstances = dba.fetchInstanceByAttribute(
            ReactomeJavaConstants.ReferenceDatabase, ReactomeJavaConstants.name, "=", "UniProt");

        if (uniProtReferenceDatabaseInstances == null || uniProtReferenceDatabaseInstances.isEmpty()) {
            throw new RuntimeException("Could not find UniProt Reference Database in " + dba);
        }
        return uniProtReferenceDatabaseInstances.iterator().next();
    }

    private GKInstance getInstanceEdit(MySQLAdaptor dba, String note) throws Exception {
        GKInstance instanceEdit = new GKInstance(dba.getSchema().getClassByName(ReactomeJavaConstants.InstanceEdit));
        instanceEdit.setDbAdaptor(dba);
        instanceEdit.setAttributeValue(
            ReactomeJavaConstants.author, getOrCreatePersonInstance(dba));
        instanceEdit.setAttributeValue(ReactomeJavaConstants.note, note);
        instanceEdit.setAttributeValue(ReactomeJavaConstants.dateTime, getCurrentDateTime());
        InstanceDisplayNameGenerator.setDisplayName(instanceEdit);
        dba.storeInstance(instanceEdit);
        return instanceEdit;
    }

    @SuppressWarnings("unchecked")
    private GKInstance getOrCreatePersonInstance(MySQLAdaptor dba) throws Exception {
        final String personSurname = "Weiser";
        final String personInitials = "JD";

        Collection<GKInstance> personInstances = dba.fetchInstancesByClass(ReactomeJavaConstants.Person);
        List<GKInstance> matchedPersonInstances = new ArrayList<>();
        for (GKInstance personInstance : personInstances) {
            String personInstanceSurname = (String) personInstance.getAttributeValue(ReactomeJavaConstants.surname);
            String personInstanceInitials = (String) personInstance.getAttributeValue(ReactomeJavaConstants.initial);
            if (personInstanceSurname != null && personInstanceSurname.equals(personSurname) &&
                personInstanceInitials != null && personInstanceInitials.equals(personInitials)) {
                matchedPersonInstances.add(personInstance);
            }
        }

        if (!matchedPersonInstances.isEmpty()) {
            return matchedPersonInstances.get(0);
        } else {
            GKInstance personInstance =
                new GKInstance(dba.getSchema().getClassByName(ReactomeJavaConstants.Person));
            personInstance.setAttributeValue(ReactomeJavaConstants.surname, personSurname);
            personInstance.setAttributeValue(ReactomeJavaConstants.initial, personInitials);
            InstanceDisplayNameGenerator.setDisplayName(personInstance);
            dba.storeInstance(personInstance);
            return personInstance;
        }
    }

    private String getCurrentDateTime() {
        return ZonedDateTime.now().format(DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss"));
    }

    private String getCurrentDate() {
        return ZonedDateTime.now().format(DateTimeFormatter.ofPattern("EEE MMM dd YYYY"));
    }

    private String getTodaysDate() {
        return ZonedDateTime.now().format(DateTimeFormatter.ofPattern("yyyyMMdd"));
    }

    private Map<Integer, String> getTaxonIdToSpeciesName() {
        Map<Integer, String> taxonIdToSpeciesName = new HashMap<>();
        taxonIdToSpeciesName.put(9606, "Homo sapiens");
        taxonIdToSpeciesName.put(10090, "Mus musculus");
        taxonIdToSpeciesName.put(10116, "Rattus norvegicus");
        taxonIdToSpeciesName.put(9913, "Bos taurus");
        taxonIdToSpeciesName.put(9031, "Gallus gallus");
        taxonIdToSpeciesName.put(7227, "Drosophila melanogaster");
        taxonIdToSpeciesName.put(6239, "Caenorhabditis elegans");
        taxonIdToSpeciesName.put(4932, "Saccharomyces cerevisiae");
        taxonIdToSpeciesName.put(4896, "Schizosaccharomyces pombe");
        taxonIdToSpeciesName.put(11695, "Human immunodeficiency virus type 1");
        taxonIdToSpeciesName.put(11718, "Human immunodeficiency virus type 2");
        taxonIdToSpeciesName.put(132504, "Influenza A virus");
        return taxonIdToSpeciesName;
    }

    @SuppressWarnings("unchecked")
    private Map<String, Long> getRGPAccessionToDbIdMap(MySQLAdaptor dba) throws Exception {
        Collection<GKInstance> instances = dba.fetchInstanceByAttribute(
            ReactomeJavaConstants.ReferenceGeneProduct,
            ReactomeJavaConstants.referenceDatabase,
            "=",
            getUniProtReferenceDatabase(dba)
        );

        Map<String, Long> identifierToDbId = new HashMap<>();
        for (GKInstance instance : instances) {
            String identifier = (String) instance.getAttributeValue(ReactomeJavaConstants.identifier);

            if (identifier != null && !identifier.isEmpty()) {
                identifierToDbId.put(identifier, instance.getDBID());
            }
        }
        return identifierToDbId;
    }

    private Map<String, Long> getIsoformAccessionToDbIdMap(MySQLAdaptor dba) throws Exception {
        return getIdentifierToDbIdMap(
            dba, ReactomeJavaConstants.ReferenceIsoform, ReactomeJavaConstants.variantIdentifier
        );
    }

    private Map<String, Long> getRDSIdentifierToDbIdMap(MySQLAdaptor dba) throws Exception {
        return getIdentifierToDbIdMap(
            dba, ReactomeJavaConstants.ReferenceDNASequence, ReactomeJavaConstants.identifier
        );
    }

    @SuppressWarnings("unchecked")
    private Map<String, Long> getIdentifierToDbIdMap(MySQLAdaptor dba, String className, String identifierAttribute) throws Exception {
        Collection<GKInstance> instances = dba.fetchInstancesByClass(className);

        Map<String, Long> identifierToDbId = new HashMap<>();
        for (GKInstance instance : instances) {
            String identifier = (String) instance.getAttributeValue(identifierAttribute);

            if (identifier != null && !identifier.isEmpty()) {
                identifierToDbId.put(identifier, instance.getDBID());
            }
        }
        return identifierToDbId;
    }

    private List<String> getSkipList() throws IOException, URISyntaxException {
        List<String> skipListIds = new ArrayList<>();

        Path resourcesDirectory = Paths.get(this.getClass().getClassLoader().getResource(".").toURI());
        List<Path> skipListFilePaths = Files.list(resourcesDirectory)
            .filter(path -> path.getFileName().toString().startsWith("skiplist"))
            .collect(Collectors.toList());

        for (Path skipListFilePath :  skipListFilePaths) {
            for (String uniProtIdToSkip : Files.readAllLines(skipListFilePath)) {
                if (isValidUniProtId(uniProtIdToSkip)) {
                    skipListIds.add(uniProtIdToSkip);
                }
            }
        }

        return skipListIds;
    }

    private boolean isValidUniProtId(String potentialUniProtId) {
        final List<Integer> validUniProtIdLengths = Arrays.asList(6, 10);
        return validUniProtIdLengths.contains(potentialUniProtId.length());
    }

    private boolean isAReferenceIsoform(GKInstance rgpInstance) {
        return rgpInstance.getSchemClass().isa(ReactomeJavaConstants.ReferenceIsoform);
    }

    private GKInstance getSpeciesInstance(MySQLAdaptor dba, String speciesName, Map<String, GKInstance> speciesCache)
        throws Exception {

        if (speciesCache != null && speciesCache.get(speciesName) != null) {
            return speciesCache.get(speciesName);
        }

        GKInstance speciesInstance = getExistingSpeciesInstance(dba, speciesName);
        if (speciesInstance != null) {
            return speciesInstance;
        } else {
            speciesInstance = createNewSpeciesInstance(dba, speciesName);
            dba.storeInstance(speciesInstance);
            return speciesInstance;
        }
    }

    @SuppressWarnings("unchecked")
    private GKInstance getHumanEnsEMBLGeneReferenceDatabase(MySQLAdaptor dba) throws Exception {
        Collection<GKInstance> ensEMBLHumanReferenceDatabaseInstances =
            dba.fetchInstanceByAttribute(
                ReactomeJavaConstants.ReferenceDatabase,
                ReactomeJavaConstants.name,
                "=",
                "ENSEMBL"
            );

        if (ensEMBLHumanReferenceDatabaseInstances == null || ensEMBLHumanReferenceDatabaseInstances.isEmpty()) {
            throw new RuntimeException("Could not get EnsEMBL human gene reference database from " + dba);
        }

        return ensEMBLHumanReferenceDatabaseInstances.iterator().next();
    }

    private GKInstance fetchReferenceDNASequenceByDbId(MySQLAdaptor dba, long referenceDNASequenceDbId) throws Exception {
        return dba.fetchInstance(referenceDNASequenceDbId);
    }

    private boolean sameDbId(GKInstance instance1, GKInstance instance2) {
        return instance1.getDBID().equals(instance2.getDBID());
    }

    private boolean areDifferentLists(List<?> list1, List<?> list2) {
        if (list1 == list2) {
            return false;
        } else if (list1 == null || list2 == null) {
            return true;
        }

        if (list1.size() != list2.size()) {
            return true;
        }

        for (Object element1 : list1) {
            if (!list2.contains(element1)) {
                return true;
            }
        }
        return false;
    }

    private boolean onEnsEMBLPrimaryAssembly(String ensEMBLGeneId) throws InterruptedException {
        final List<String> primaryAssemblyRegions = Arrays.asList(
            "1","2","3","4","5","6","7","8","9","10",
            "11","12","13","14","15","16","17","18","19","20",
            "21","22","X","Y","MT"
        );
        final int maxQueryAttempts = 5;

        int queryAttempts = 0;
        String ensEMBLIdData = null;
        while (queryAttempts < maxQueryAttempts && ensEMBLIdData == null) {
            try {
                queryAttempts += 1;
                ensEMBLIdData = queryEnsEMBLRESTAPI(ensEMBLGeneId);
            } catch (IOException e) {
                System.out.println("IOException when querying  " + ensEMBLGeneId + ": " + e.getMessage());
                // Sleep progressively longer for each query attempt to allow the server more time to respond
                Thread.sleep(queryAttempts * 500);
            }
        }
        Pattern seqRegionPattern = Pattern.compile("\"seq_region_name\":(\".*?\")");
        Matcher seqRegionMatcher = seqRegionPattern.matcher(ensEMBLIdData);

        if (seqRegionMatcher.find()) {
            String seqRegion = seqRegionMatcher.group(1);
            return primaryAssemblyRegions.contains(seqRegion);
        }
        return false;
    }

    private String queryEnsEMBLRESTAPI(String ensEMBLGeneId) throws IOException, InterruptedException {
        URL ensemblLookupURL = new URL(
            "https://rest.ensembl.org/lookup/id/" + ensEMBLGeneId + "?content-type=application/json"
        );

        HttpURLConnection httpURLConnection = (HttpURLConnection) ensemblLookupURL.openConnection();
        if (httpURLConnection.getResponseCode() == HttpURLConnection.HTTP_BAD_REQUEST) {
            if (getError(httpURLConnection).contains("not found")) {
                return "";
            } else {
                System.out.println(String.format("Bad request for %s:  Sleeping for 5 seconds and retrying", ensemblLookupURL));
                Thread.sleep(5000);
            }

            return queryEnsEMBLRESTAPI(ensEMBLGeneId);
        }
        BufferedReader ensEMBLInputReader =
            new BufferedReader(new InputStreamReader(httpURLConnection.getInputStream()));
        String inputLine;
        StringBuilder content = new StringBuilder();
        while ((inputLine = ensEMBLInputReader.readLine()) != null) {
            content.append(inputLine);
        }
        ensEMBLInputReader.close();
        httpURLConnection.disconnect();

        return content.toString();
    }

    private String getError(HttpURLConnection urlConnection) throws IOException {
        BufferedReader errorInputReader =
            new BufferedReader(new InputStreamReader(urlConnection.getErrorStream()));
        String inputLine;
        StringBuilder content = new StringBuilder();
        while ((inputLine = errorInputReader.readLine()) != null) {
            content.append(inputLine);
        }
        errorInputReader.close();
        return content.toString();
    }

    private void updateInstance(
        GKInstance instance, Map<String, List<?>> values, BufferedWriter sequenceReportWriter) throws Exception {

        boolean isInstanceChanged = false;
        for (String attributeName : values.keySet()) {
            List<?> newValuesForAttribute =
                values.get(attributeName).stream().filter(Objects::nonNull).collect(Collectors.toList());
            if (newValuesForAttribute.size() == 0) {
                System.out.println("WARNING: No new values for " + attributeName + " on " + instance.getDBID() +
                    " skipping attribute update");
                continue;
            }

            if (attributeName.toLowerCase().equals(ReactomeJavaConstants.checksum)) {
                Boolean oldSequenceChangedValue =
                    (Boolean) instance.getAttributeValue("isSequenceChanged");

                Boolean newSequenceChangedValue =
                    getNewIsSequenceChangedAttributeValue(instance, newValuesForAttribute);

                if (oldSequenceChangedValue == null || !oldSequenceChangedValue.equals(newSequenceChangedValue)) {
                    instance.setAttributeValue("isSequenceChanged", newSequenceChangedValue);
                    System.out.println(String.format("%s (%d) has a new is_sequence_changed value",
                        instance.getDisplayName(), instance.getDBID()));
                    isInstanceChanged = true;
                }
            }

            if (attributeName.toLowerCase().equals(ReactomeJavaConstants.chain)) {
                boolean chainChangeLogUpdated =
                    updateChainLog(instance, (List<String>) newValuesForAttribute, sequenceReportWriter);
                if (hasChains(instance) && chainChangeLogUpdated) {
                    List<GKInstance> ewasInstances = getAllEwasInstances(instance);

                    for (GKInstance ewasInstance : ewasInstances) {
                        reportChangedChainForEWASInstance(instance, ewasInstance);
                    }
                    /*
                    Collection<GKInstance> allEwasInstances = new ArrayList<>();
                    Collection<GKInstance> referenceEntityEwasInstances =
                        instance.getReferers(ReactomeJavaConstants.referenceEntity);
                    if (referenceEntityEwasInstances != null) {
                        allEwasInstances.addAll(referenceEntityEwasInstances);
                    }
                    Collection<GKInstance> hasModifiedResidueInstances = new ArrayList<>();
                    Collection<GKInstance> referenceSequenceModifiedResidues =
                        instance.getReferers(ReactomeJavaConstants.referenceSequence);
                    if (referenceSequenceModifiedResidues != null) {
                        hasModifiedResidueInstances.addAll(referenceSequenceModifiedResidues);
                    }

                    Collection<GKInstance> secondReferenceSequenceModifiedResidues =
                        instance.getReferers(ReactomeJavaConstants.secondReferenceSequence);
                    if (secondReferenceSequenceModifiedResidues != null) {
                        hasModifiedResidueInstances.addAll(secondReferenceSequenceModifiedResidues);
                    }

                    for (GKInstance hasModifiedResidueInstance : hasModifiedResidueInstances) {
                        Collection<GKInstance> hasModifiedEwasInstances =
                            hasModifiedResidueInstance.getReferers(ReactomeJavaConstants.hasModifiedResidue);
                        if (hasModifiedEwasInstances != null) {
                            allEwasInstances.addAll(hasModifiedEwasInstances);
                        }
                    }
*/
                }
            }

            if (valuesChanged(instance, attributeName, newValuesForAttribute)) {
                if (isSingleAttribute(attributeName)) {
                    instance.setAttributeValue(attributeName, newValuesForAttribute.get(0));
                } else {
                    instance.setAttributeValue(attributeName, newValuesForAttribute);
                }

                isInstanceChanged = true;
            }
        }

        if (isInstanceChanged) {
            MySQLAdaptor dba = (MySQLAdaptor) instance.getDbAdaptor();
            if (instance.getDBID() == null) {
                dba.storeInstance(instance);
            } else {
                dba.updateInstance(instance);
            }
        }
    }

    private Boolean getNewIsSequenceChangedAttributeValue(GKInstance instance, List<?> newValues) throws Exception {
        String oldChecksum = (String) instance.getAttributeValue(ReactomeJavaConstants.checksum);
        String newChecksum = newValues.get(0).toString();

        return isSequenceChanged(oldChecksum, newChecksum);
    }

    private boolean isSequenceChanged(String oldChecksum, String newChecksum) {
        return oldChecksum != null && newChecksum != null && !oldChecksum.equals(newChecksum);
    }

    @SuppressWarnings("unchecked")
    private boolean updateChainLog(GKInstance instance, List<String> newChainValues, BufferedWriter sequenceReportWriter)
        throws Exception {
        boolean chainLogChanged = false;

        List<String> oldChainValues = instance.getAttributeValuesList(ReactomeJavaConstants.chain);
        String date = getCurrentDate();

        String referenceGeneProductDescription = getReferenceGeneProductDescription(instance);

        for (String oldChainValue : oldChainValues) {
            if (!newChainValues.contains(oldChainValue)) {
                String logEntry = String.format("%s for %d removed on %s", oldChainValue, instance.getDBID(), date);
                sequenceReportWriter.write(logEntry + " for " + referenceGeneProductDescription + "\n");

                String existingLog = (String) instance.getAttributeValue("_chainChangeLog");
                String fullLog =
                    existingLog != null ?
                    existingLog + ";" + logEntry :
                    logEntry;

                instance.addAttributeValue("_chainChangeLog", fullLog);
                System.out.println("old chain removed for " + instance.getDBID());
                chainLogChanged = true;
            }
        }

        for (String newChainValue : newChainValues) {
            if (!oldChainValues.contains(newChainValue)) {
                String logEntry = String.format("%s for %d added on %s", newChainValue, instance.getDBID(), date);
                sequenceReportWriter.write(logEntry + " for " + referenceGeneProductDescription + "\n");


                String existingLog = (String) instance.getAttributeValue("_chainChangeLog");
                String fullLog =
                    existingLog != null ?
                        existingLog + ";" + logEntry :
                        logEntry;

                instance.addAttributeValue("_chainChangeLog", fullLog);
                System.out.println("new chain added for " + instance.getDBID());
                chainLogChanged = true;
            }
        }
        return chainLogChanged;
    }

    private boolean hasChains(GKInstance instance) throws Exception {
        List<String> chainValues = instance.getAttributeValuesList(ReactomeJavaConstants.chain);
        return chainValues != null && !chainValues.isEmpty();
    }

    private String getReferenceGeneProductDescription(GKInstance rgpInstance) throws Exception {
        String referenceGeneProductDescription = rgpInstance.getDBID() != null ? rgpInstance.getDBID().toString() : "";

        String rgpName = (String) rgpInstance.getAttributeValue(ReactomeJavaConstants.name);
        if (rgpName != null && !rgpName.isEmpty()) {
            referenceGeneProductDescription += " - " + rgpName;
        }

        GKInstance speciesInstance = (GKInstance) rgpInstance.getAttributeValue(ReactomeJavaConstants.species);
        if (speciesInstance != null) {
            referenceGeneProductDescription += " (" + speciesInstance.getDisplayName() + ")";
        }

        return referenceGeneProductDescription;
    }

    private boolean valuesChanged(GKInstance instance, String attributeName, List<?> newValues) throws Exception {
        List<?> currentValuesToCompare;
        List<?> newValuesToCompare;

        List<?> currentValues = instance.getAttributeValuesList(attributeName);

        if (instance.getSchemClass().getAttribute(attributeName).isInstanceTypeAttribute()) {
            currentValuesToCompare = currentValues.stream()
                .map(inst -> ((GKInstance) inst).getDBID()).collect(Collectors.toList());
            newValuesToCompare = newValues.stream()
                .map(inst -> ((GKInstance) inst).getDBID()).collect(Collectors.toList());
        } else {
            currentValuesToCompare = currentValues.stream().map(Object::toString).collect(Collectors.toList());
            newValuesToCompare = newValues.stream().map(Object::toString).collect(Collectors.toList());
        }

        if (!areDifferentLists(currentValuesToCompare, newValuesToCompare)) {
            return false;
        }

        System.out.println(String.format("%s changed for instance %d", attributeName, instance.getDBID()));
        System.out.println(String.format("old attribute values - %s",
            currentValuesToCompare.stream().map(Object::toString).sorted().collect(Collectors.joining(","))));
        System.out.println(String.format("new attribute values - %s",
            newValuesToCompare.stream().map(Object::toString).sorted().collect(Collectors.joining(","))));

        return true;
    }

    private boolean isSingleAttribute(String attributeName) {
        return Arrays.asList(
            ReactomeJavaConstants.species,
            ReactomeJavaConstants.sequenceLength,
            ReactomeJavaConstants.checksum,
            ReactomeJavaConstants.comment
        ).contains(attributeName);
    }

    private String getSpeciesName(GKInstance instance) throws Exception {
        GKInstance species = (GKInstance) instance.getAttributeValue(ReactomeJavaConstants.species);
        if (species != null) {
            return species.getDisplayName();
        }
        return "";
    }

    private boolean hasUniProtReferenceDatabase(GKInstance instance) {
        GKInstance referenceDatabase;
        try {
            referenceDatabase = (GKInstance) instance.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
        } catch (Exception e) {
            throw new RuntimeException("Unable to fetch reference databases from " + instance);
        }

        return referenceDatabase != null && referenceDatabase.getDisplayName().toLowerCase().contains("uniprot");
    }

    @SuppressWarnings("unchecked")
    private GKInstance getExistingSpeciesInstance(MySQLAdaptor dba, String speciesName) throws Exception {
        Collection<GKInstance> speciesInstances = dba.fetchInstanceByAttribute(
            ReactomeJavaConstants.Species,
            ReactomeJavaConstants.name,
            "=",
            speciesName
        );

        if (speciesInstances == null) {
            return null;
        }
        return speciesInstances.iterator().next();
    }

    private GKInstance createNewSpeciesInstance(MySQLAdaptor dba, String speciesName) throws Exception {
        GKInstance speciesInstance = new GKInstance(dba.getSchema().getClassByName(ReactomeJavaConstants.Species));
        speciesInstance.setDbAdaptor(dba);
        speciesInstance.setAttributeValue(ReactomeJavaConstants.name, speciesName);
        InstanceDisplayNameGenerator.setDisplayName(speciesInstance);
        return speciesInstance;
    }

    @SuppressWarnings("unchecked")
    private Collection<GKInstance> getRGPReferrers(GKInstance rgpInstance) throws Exception {
        Collection<GKInstance> referrers = new ArrayList<>();

        final List<String> reverseAttributes = Arrays.asList(
            ReactomeJavaConstants.referenceEntity,
            ReactomeJavaConstants.referenceSequence,
            ReactomeJavaConstants.secondReferenceSequence,
            ReactomeJavaConstants.isoformParent
        );

        for (String reverseAttribute : reverseAttributes) {
            Collection<GKInstance> reverseAttributeReferrers = rgpInstance.getReferers(reverseAttribute);
            if (reverseAttributeReferrers != null) {
                referrers.addAll(reverseAttributeReferrers);
            }
        }
        return referrers;
    }


    private String matchSingleValue(String entry, String patternString) {
        List<String> values = matchMultipleValues(entry, patternString);
        return !values.isEmpty() ? values.get(0) : "";
    }

    private List<String> matchMultipleValues(String entry, String patternString) {
        Pattern pattern = Pattern.compile(patternString, Pattern.MULTILINE);
        Matcher matcher = pattern.matcher(entry);

        List<String> multipleValues = new ArrayList<>();
        while (matcher.find()) {
            multipleValues.add(matcher.group(1));
        }
        return multipleValues;
    }

    private String removeVersionNumber(String identifier) {
        Pattern identifierWithVersionNumberPattern = Pattern.compile("(.*)\\.\\d+$");
        Matcher identifierWithVersionNumberMatcher = identifierWithVersionNumberPattern.matcher(identifier);

        String identifierWithoutVersionNumber;

        if (identifierWithVersionNumberMatcher.find()) {
            identifierWithoutVersionNumber = identifierWithVersionNumberMatcher.group(1);
        } else {
            identifierWithoutVersionNumber = identifier;
        }
        return identifierWithoutVersionNumber;
    }

    private String parseComments(String entry) {
        Pattern commentsPattern = Pattern.compile(
            "<comment type=\"([A-Za-z ]*?)\".*?\\s+<text.*?>(.*?)</text>", Pattern.MULTILINE);
        Matcher commentsMatcher = commentsPattern.matcher(entry);

        StringBuilder comments = new StringBuilder();
        while (commentsMatcher.find()) {
            String commentType = commentsMatcher.group(1).toUpperCase();
            String commentText = commentsMatcher.group(2);

            comments.append(commentType).append(" ").append(commentText);
        }
        return comments.toString();
    }

    private List<String> parseChains(String entry) {
        List<String> featuresTypes = Arrays.asList(
            "initiator methionine",
            "chain",
            "peptide",
            "propeptide",
            "signal peptide",
            "transit peptide"
        );
        Pattern chainsPattern =
            Pattern.compile( "<feature.*?type=\"(" + String.join("|", featuresTypes) + ")\"(.*?)</feature>");
        Matcher chainsMatcher = chainsPattern.matcher(entry);

        List<String> chains = new ArrayList<>();
        while (chainsMatcher.find()) {
            String featureType = chainsMatcher.group(1);
            String featureContent = chainsMatcher.group(2);

            if (featureType.equals("initiator methionine")) {
                chains.add(parseInitiatorMethionineChain(featureContent));
            } else {
                chains.add(parseGenericChain(featureType, featureContent));
            }
        }
        return chains;
    }

    private String parseInitiatorMethionineChain(String featureContent) {
        String position = matchSingleValue(featureContent, "<location>\\n\\s+<position position=\"(\\d+)\"");
        return "initiator methionine:" + position;
    }

    private String parseGenericChain(String featureType, String featureContent) {
        String chainStart = matchSingleValue(featureContent, "<begin position=\"(\\d+)\"");
        String chainEnd = matchSingleValue(featureContent, "<end position=\"(\\d+)");

        return featureType + ":" + chainStart + "-" + chainEnd;
    }

    @SuppressWarnings("unchecked")
    private void updateDisplayNames(MySQLAdaptor dba, String className) throws Exception {
        Collection<GKInstance> instances = dba.fetchInstancesByClass(className);
        for (GKInstance instance : instances) {
            InstanceDisplayNameGenerator.setDisplayName(instance);
            dba.updateInstanceAttribute(instance, ReactomeJavaConstants._displayName);
        }
    }

    private List<GKInstance> getAllEwasInstances(GKInstance referenceGeneProduct) throws Exception {
        List<GKInstance> allEwasInstances = new ArrayList<>();

        Collection<GKInstance> referenceEntityEwasInstances =
            referenceGeneProduct.getReferers(ReactomeJavaConstants.referenceEntity);
        if (referenceEntityEwasInstances != null) {
            allEwasInstances.addAll(referenceEntityEwasInstances);
        }
        Collection<GKInstance> hasModifiedResidueInstances = new ArrayList<>();
        Collection<GKInstance> referenceSequenceModifiedResidues =
            referenceGeneProduct.getReferers(ReactomeJavaConstants.referenceSequence);
        if (referenceSequenceModifiedResidues != null) {
            hasModifiedResidueInstances.addAll(referenceSequenceModifiedResidues);
        }

        Collection<GKInstance> secondReferenceSequenceModifiedResidues =
            referenceGeneProduct.getReferers(ReactomeJavaConstants.secondReferenceSequence);
        if (secondReferenceSequenceModifiedResidues != null) {
            hasModifiedResidueInstances.addAll(secondReferenceSequenceModifiedResidues);
        }

        for (GKInstance hasModifiedResidueInstance : hasModifiedResidueInstances) {
            Collection<GKInstance> hasModifiedEwasInstances =
                hasModifiedResidueInstance.getReferers(ReactomeJavaConstants.hasModifiedResidue);
            if (hasModifiedEwasInstances != null) {
                allEwasInstances.addAll(hasModifiedEwasInstances);
            }
        }

        return allEwasInstances;
    }

    private void reportChangedChainForEWASInstance(GKInstance referenceGeneProduct, GKInstance ewas)
        throws Exception {
        //"RGP db id\tRGP Accession\tEWAS db id\tEWAS name\tEWAS author (created or last modified)\t"
        String reportLine = String.join("\t",
            referenceGeneProduct.getDBID().toString(),
            (String) referenceGeneProduct.getAttributeValue(ReactomeJavaConstants.identifier),
            ewas.getDBID().toString(),
            ewas.getDisplayName(),
            getAuthor(ewas)
        ).concat(System.lineSeparator());

        Files.write(
            getUniprotUpdateDirectoryPath().resolve("ewasCoordinatesReport.txt"),
            reportLine.getBytes(),
            StandardOpenOption.CREATE, StandardOpenOption.APPEND
        );
    }

    private String getAuthor(GKInstance ewas) throws Exception {
        GKInstance ewasCreatedInstanceEdit = (GKInstance) ewas.getAttributeValue(ReactomeJavaConstants.created);
        if (ewasCreatedInstanceEdit != null) {
            return getAuthorFromInstanceEdit(ewasCreatedInstanceEdit);
        } else {
            return getLastModifiedAuthor(ewas);
        }
    }

    private String getLastModifiedAuthor(GKInstance ewas) throws Exception {
        List<GKInstance> ewasModifiedInstanceEdits = ewas.getAttributeValuesList(ReactomeJavaConstants.modified);
        if (ewasModifiedInstanceEdits != null && !ewasModifiedInstanceEdits.isEmpty()) {
            GKInstance ewasMostRecentModifiedInstanceEdit = ewasModifiedInstanceEdits.get(0);
            return getAuthorFromInstanceEdit(ewasMostRecentModifiedInstanceEdit);
        } else {
            return "Unknown author";
        }
    }

    private String getAuthorFromInstanceEdit(GKInstance instanceEdit) throws Exception {
        GKInstance instanceEditAuthor = (GKInstance) instanceEdit.getAttributeValue(ReactomeJavaConstants.author);
        if (instanceEditAuthor != null) {
            return instanceEditAuthor.getDisplayName();
        } else {
            return "Unknown author";
        }
    }
}
