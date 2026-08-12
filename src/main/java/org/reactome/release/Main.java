package org.reactome.release;

import org.gk.model.ReactomeJavaConstants;
import org.reactome.curation.model.SimpleInstance;
import org.reactome.release.reports.DuplicateAccessionReport;
import org.reactome.release.reports.Reportable;
import org.reactome.release.reports.TrEMBLAccessionReport;
import org.reactome.release.utils.CuratorToolAPI;

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

import static org.reactome.release.utils.Utils.emptyListIfNull;
import static org.reactome.release.utils.Utils.isTrEMBLId;

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
        CuratorToolAPI curatorToolAPI = new CuratorToolAPI(Long.parseLong(configProperties.getProperty("personId")));

        List<String> skipList = getSkipList();

        this.uniprotUpdateDirectoryPath = Paths.get(configProperties.getProperty("uniprotUpdateDirectory"));

        SimpleInstance uniProtReferenceDatabase = curatorToolAPI.fetchUniProtReferenceDatabase();

        Map<Integer, String> taxonIdToSpeciesName = getTaxonIdToSpeciesName();

        // Counters
        int totalNumberOfDbInstances;
        int numberOfInstancesInSwissProtFile = 0;
        int numberOfObsoleteInstancesWithNoEWAS = 0;
        int numberOfNewSwissProtInstances = 0;

        Map<Long,String> duplicateDbIdToReferenceGeneProductAccession = new HashMap<>();

        System.out.println("Populating rgp accession to db id...");
        Map<String, Long> rgpAccessionToDbId = curatorToolAPI.getRGPAccessionToDbIdMap();
        totalNumberOfDbInstances = rgpAccessionToDbId.size();
        System.out.println("Populating isoform accession to db id...");
        Map<String, Long> isoformAccessionToDbId = curatorToolAPI.getIsoformAccessionToDbIdMap();
        System.out.println("Populating rds identifier to db id...");
        Map<String, Long> rdsIdentifierToDbId = curatorToolAPI.getRDSIdentifierToDbIdMap();

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
                        System.out.println(String.format("%d records processed and committed", recordCounter));
                    }
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
                SimpleInstance speciesInstance = null;
                for (String speciesName : taxonIdToSpeciesName.values()) {
                    if (organismName.contains(speciesName)) {
                        taxon = speciesName;
                        speciesInstance = curatorToolAPI.getSpeciesInstance(taxon);
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

                Integer sequenceLength = parseSequenceLength(entry, primaryAccession);

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

                List<SimpleInstance> referenceDNASequences = new ArrayList<>();
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
                    SimpleInstance humanEnsEMBLGeneReferenceDatabase = curatorToolAPI.getHumanEnsEMBLGeneReferenceDatabase();

                    for (String ensEMBLGeneId : uniqueEnsEMBLGeneIds) {
                        SimpleInstance referenceDNASequence;

                        if (rdsIdentifierToDbId.containsKey(ensEMBLGeneId)) {
                            referenceDNASequenceReportWriter.write("Checking existing reference DNA sequence for " +
                                ensEMBLGeneId + " with db_id " + rdsIdentifierToDbId.get(ensEMBLGeneId) + "\n");

                            long rdsDbId = rdsIdentifierToDbId.get(ensEMBLGeneId);
                            referenceDNASequence = fetchReferenceDNASequenceByDbId(curatorToolAPI, rdsDbId);

                            SimpleInstance existingRDSReferenceDatabase = (SimpleInstance)
                                referenceDNASequence.getAttribute(ReactomeJavaConstants.referenceDatabase);
                            boolean isUpdateToReferenceDNASequence = false;
                            if (existingRDSReferenceDatabase == null ||
                                !sameDbId(existingRDSReferenceDatabase, humanEnsEMBLGeneReferenceDatabase)) {
                                referenceDNASequence.setAttribute(
                                    ReactomeJavaConstants.referenceDatabase, humanEnsEMBLGeneReferenceDatabase);
                                isUpdateToReferenceDNASequence = true;
                            }

                            List<String> existingGeneNames = (List<String>)
                                referenceDNASequence.getAttribute(ReactomeJavaConstants.geneName);
                            if (existingGeneNames == null || areDifferentLists(existingGeneNames, geneNames)) {
                                referenceDNASequence.setAttribute(ReactomeJavaConstants.geneName, geneNames);
                                isUpdateToReferenceDNASequence = true;
                            }
                            SimpleInstance existingSpeciesInstance = (SimpleInstance)
                                referenceDNASequence.getAttribute(ReactomeJavaConstants.species);
                            if (existingSpeciesInstance == null ||
                                (speciesInstance != null &&
                                !existingSpeciesInstance.getDbId().equals(speciesInstance.getDbId()))) {

                                referenceDNASequence.setAttribute(ReactomeJavaConstants.species, speciesInstance);
                                isUpdateToReferenceDNASequence = true;
                            }

                            String existingIdentifier = (String)
                                referenceDNASequence.getAttribute(ReactomeJavaConstants.identifier);
                            if (existingIdentifier == null || !existingIdentifier.equals(ensEMBLGeneId)) {
                                referenceDNASequence.setAttribute(ReactomeJavaConstants.identifier, ensEMBLGeneId);
                                isUpdateToReferenceDNASequence = true;
                            }

                            if (isUpdateToReferenceDNASequence) {
                                referenceDNASequenceReportWriter.write(
                                    "Updating existing reference DNA sequence for " + ensEMBLGeneId + " with db_id " +
                                    rdsIdentifierToDbId.get(ensEMBLGeneId) + "\n"
                                );
                                referenceDNASequence.setDisplayName(curatorToolAPI.getReferenceSequenceDisplayName(referenceDNASequence));
                                curatorToolAPI.commit(referenceDNASequence);
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

                            referenceDNASequence = new SimpleInstance();
                            referenceDNASequence.setSchemaClassName(ReactomeJavaConstants.ReferenceDNASequence);
                            referenceDNASequence.setAttribute(
                                ReactomeJavaConstants.referenceDatabase, humanEnsEMBLGeneReferenceDatabase);
                            referenceDNASequence.setAttribute(ReactomeJavaConstants.identifier, ensEMBLGeneId);
                            referenceDNASequence.setAttribute(ReactomeJavaConstants.geneName, geneNames);
                            referenceDNASequence.setAttribute(ReactomeJavaConstants.species, speciesInstance);

                            referenceDNASequence.setDisplayName(curatorToolAPI.getReferenceSequenceDisplayName(referenceDNASequence));

                            long referenceDNASequenceDbId = curatorToolAPI.commit(referenceDNASequence).getDbId();
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

                    SimpleInstance newReferenceGeneProductInstance = new SimpleInstance();
                    newReferenceGeneProductInstance.setSchemaClassName(ReactomeJavaConstants.ReferenceGeneProduct);
                    newReferenceGeneProductInstance.setAttribute(
                        ReactomeJavaConstants.referenceDatabase, uniProtReferenceDatabase);
                    newReferenceGeneProductInstance.setAttribute(ReactomeJavaConstants.identifier, primaryAccession);
                    long newRGPDbId = curatorToolAPI.commit(newReferenceGeneProductInstance).getDbId();

                    System.out.println(String.format("New UniProt:%s\t%d", primaryAccession, newRGPDbId));
                    updateInstance(curatorToolAPI, newReferenceGeneProductInstance, values, sequenceReportWriter);
                    for (String isoformId : isoformIds) {
                        if (!isoformId.contains(primaryAccession)) {
                            misMatchedIsoformAccessionToRGPAccession.put(isoformId, primaryAccession);
                        }

                        SimpleInstance newIsoformInstance = new SimpleInstance();
                        newIsoformInstance.setSchemaClassName(ReactomeJavaConstants.ReferenceIsoform);
                        newIsoformInstance.setAttribute(
                            ReactomeJavaConstants.referenceDatabase, uniProtReferenceDatabase);
                        newIsoformInstance.setAttribute(ReactomeJavaConstants.identifier, primaryAccession);
                        newIsoformInstance.setAttribute(
                            ReactomeJavaConstants.isoformParent, newReferenceGeneProductInstance);
                        newIsoformInstance.setAttribute(ReactomeJavaConstants.variantIdentifier, isoformId);

                        updateInstance(curatorToolAPI, newIsoformInstance, values, sequenceReportWriter);
                    }
                } else {
                    Collection<SimpleInstance> existingReferenceGeneProductInstances = curatorToolAPI.getReferenceGeneProductsByIdentifier(primaryAccession);
                    boolean duplicateFlag = false;
                    for (SimpleInstance existingReferenceGeneProductInstance : existingReferenceGeneProductInstances) {
                        if (isAReferenceIsoform(existingReferenceGeneProductInstance)) {
                            continue;
                        }

                        if (duplicateFlag) {
                            duplicateDbIdToReferenceGeneProductAccession.put(
                                existingReferenceGeneProductInstance.getDbId(), primaryAccession);
                            continue;
                        }

                        System.out.println(String.format("Updating master sequence...%d\t%s",
                            existingReferenceGeneProductInstance.getDbId(), primaryAccession));

                        updateInstance(curatorToolAPI, existingReferenceGeneProductInstance, values, sequenceReportWriter);

                        duplicateFlag = true;

                        if (values.get(ReactomeJavaConstants.species).isEmpty()) {
                            values.put(ReactomeJavaConstants.species, Collections.singletonList((SimpleInstance)
                                existingReferenceGeneProductInstance.getAttribute(ReactomeJavaConstants.species))
                            );
                        }
                        for (String isoformId : isoformIds) {
                            if (isoformId.contains(primaryAccession)) {
                                List<SimpleInstance> isoformInstances = curatorToolAPI.getReferenceIsoformByVariantIdentifier(isoformId);
                                if (!isoformInstances.isEmpty()) {
                                    for (SimpleInstance isoformInstance : isoformInstances) {
                                        String isoformAccession = (String) isoformInstance.getAttribute(
                                            ReactomeJavaConstants.variantIdentifier);
                                        if (!isoformAccession.contains(primaryAccession)) {
                                            continue;
                                        }
                                        System.out.println(String.format("Existing isoform update: %s\tMaster: %d",
                                            isoformAccession, existingReferenceGeneProductInstance.getDbId()));

                                        isoformInstance.setAttribute(ReactomeJavaConstants.isoformParent,
                                            existingReferenceGeneProductInstance);

                                        updateInstance(curatorToolAPI, isoformInstance, values, sequenceReportWriter);

                                        isoformAccessionToDbId.remove(isoformId);
                                    }
                                } else {
                                    SimpleInstance isoformInstance = new SimpleInstance();
                                    isoformInstance.setSchemaClassName(ReactomeJavaConstants.ReferenceIsoform);
                                    isoformInstance.setAttribute(ReactomeJavaConstants.identifier,
                                        primaryAccession);
                                    isoformInstance.setAttribute(ReactomeJavaConstants.isoformParent,
                                        existingReferenceGeneProductInstance);
                                    isoformInstance.setAttribute(ReactomeJavaConstants.variantIdentifier,
                                        isoformId);
                                    long isoformDbId = curatorToolAPI.commit(isoformInstance).getDbId();

                                    System.out.println(String.format("New isoform: %s\t%d\tMaster: %d",
                                        isoformId, isoformDbId, existingReferenceGeneProductInstance.getDbId()));

                                    updateInstance(curatorToolAPI, isoformInstance, values, sequenceReportWriter);
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

        System.out.println("Updating mis-matched isoforms");

        for (String misMatchedIsoformAccession : misMatchedIsoformAccessionToRGPAccession.keySet()) {
            List<SimpleInstance> isoformParents = new ArrayList<>();

            List<SimpleInstance> isoformInstances = curatorToolAPI.getReferenceIsoformByVariantIdentifier(misMatchedIsoformAccession);

            SimpleInstance isoformInstance = !isoformInstances.isEmpty() ? isoformInstances.get(0) : null;
            if (isoformInstance != null) {
                SimpleInstance isoformParent =
                    (SimpleInstance) isoformInstance.getAttribute(ReactomeJavaConstants.isoformParent);
                if (isoformParent == null) {
                    continue;
                }
                isoformParents.add(isoformParent);
            }

            List<SimpleInstance> mismatchedParents = curatorToolAPI.getReferenceGeneProductsByIdentifier(misMatchedIsoformAccession);

            SimpleInstance mismatchedParent = !mismatchedParents.isEmpty() ? mismatchedParents.get(0) : null;
            if (mismatchedParent != null && isoformInstance != null) {
                isoformParents.add(mismatchedParent);
                long isoformInstanceDbId = isoformInstance.getDbId();
                System.out.println(String.format("Mismatched parent: %s(%d)\t%s\n",
                    misMatchedIsoformAccession,
                    isoformInstanceDbId,
                    misMatchedIsoformAccessionToRGPAccession.get(misMatchedIsoformAccession))
                );

                isoformInstance.setAttribute(ReactomeJavaConstants.isoformParent, isoformParents);
                curatorToolAPI.commit(isoformInstance);
            }
        }

        System.out.println("Mis-matched isoform updates complete");

        System.out.println("Updating display names...");

        curatorToolAPI.updateReferenceGeneProductDisplayNames();
        //curatorToolAPI.updateReferenceIsoformDisplayNames();

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
                List<SimpleInstance> obsoleteReferenceGeneProductInstances =
                    curatorToolAPI.getReferenceGeneProductsByIdentifier(rgpAccession);

                boolean isObsoleteRGPDeleted = false;
                for (SimpleInstance obsoleteReferenceGeneProductInstance : obsoleteReferenceGeneProductInstances) {
                    String variantIdentifier = null;
                    if (isAReferenceIsoform(obsoleteReferenceGeneProductInstance)) {
                        variantIdentifier = (String) obsoleteReferenceGeneProductInstance.getAttribute(
                            ReactomeJavaConstants.variantIdentifier);
                    }

                    if (variantIdentifier != null) {
                        continue;
                    }

                    long obsoleteRGPDbId = obsoleteReferenceGeneProductInstance.getDbId();
                    List<SimpleInstance> referrers = getRGPReferrers(curatorToolAPI, obsoleteReferenceGeneProductInstance);
                    if (referrers == null || referrers.isEmpty()) {
                        System.out.println("Deleting " + obsoleteRGPDbId + "...");
                        curatorToolAPI.deleteInstance(obsoleteReferenceGeneProductInstance);
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
            List<SimpleInstance> isoformInstances = curatorToolAPI.getReferenceIsoformByVariantIdentifier(isoformAccession);

            SimpleInstance isoformInstance = !isoformInstances.isEmpty() ? isoformInstances.get(0) : null;
            if (isoformInstance == null) {
                System.out.println(isoformAccession + " is not a variant identifier for any ReferenceIsoform");
                continue;
            }

            long obsoleteIsoformDbId = isoformInstance.getDbId();
            SimpleInstance isoformParent =
                (SimpleInstance) isoformInstance.getAttribute(ReactomeJavaConstants.isoformParent);
            if (isoformParent == null) {
                System.out.println(isoformInstance.getDbId());
                dbIdsToSkip.add(obsoleteIsoformDbId);
                continue;
            }

            String isoformParentIdentifier = (String)
                isoformParent.getAttribute(ReactomeJavaConstants.identifier);
            if (isoformParentIdentifier == null || isoformParentIdentifier.isEmpty()) {
                continue;
            }

            List<SimpleInstance> referrers = getRGPReferrers(curatorToolAPI, isoformInstance);
            if (referrers == null || referrers.isEmpty()) {
                System.out.println("Deleting " + obsoleteIsoformDbId + "...");
                curatorToolAPI.deleteInstance(isoformInstance);
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

        List<String> plantReplaceableReportLines = new ArrayList<>();
        List<String> plantNoReplacementReportLines = new ArrayList<>();

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

                List<SimpleInstance> obsoleteRGPInstances = curatorToolAPI.getReferenceGeneProductsByIdentifier(rgpAccession);
                for (SimpleInstance obsoleteRGPInstance : obsoleteRGPInstances) {
                    String variantIdentifier = null;
                    if (isAReferenceIsoform(obsoleteRGPInstance)) {
                        variantIdentifier = (String) obsoleteRGPInstance.getAttribute(
                            ReactomeJavaConstants.variantIdentifier);
                    }

                    if (variantIdentifier != null) {
                        continue;
                    }
                    obsoleteDbId = obsoleteRGPInstance.getDbId();
                    speciesName = getSpeciesName(obsoleteRGPInstance);

                    List<SimpleInstance> referrers = emptyListIfNull(
                        curatorToolAPI.getReferrers(obsoleteRGPInstance, ReactomeJavaConstants.referenceEntity)
                    );
                    for (SimpleInstance referrer : referrers) {
                        if (referrer.getSchemaClassName().equals(ReactomeJavaConstants.EntityWithAccessionedSequence)) {
                            referrerDbIds.add(referrer.getDbId());
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
                            "|[https://newcurator.reactome.org/curatorgraph/dataSchema/DatabaseObject/instance/%d %d]\n",
                            obsoleteDbId, obsoleteDbId
                            ));
                        reportLineBuilder.append(String.format("||%s\n", String.join(
                            "|", referrerDbIds.stream().map(Object::toString).collect(Collectors.toList())
                        )));
                        reportLineBuilder.append("|" + speciesName + "\n");
                        reportLineBuilder.append("|-\n");

                        String reportLine = reportLineBuilder.toString();
                        if (skipList.stream().anyMatch(accession -> accession.equals(rgpAccession))) {
                            skipReplaceableReportLines.add(reportLine);
                        } else if (isPlantReportLine(reportLine)) {
                            plantReplaceableReportLines.add(reportLine);
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

            List<SimpleInstance> obsoleteRGPInstances = curatorToolAPI.getReferenceGeneProductsByIdentifier(rgpAccession);

            for (SimpleInstance obsoleteRGPInstance : obsoleteRGPInstances) {
                String variantIdentifier = null;
                if (isAReferenceIsoform(obsoleteRGPInstance)) {
                    variantIdentifier = (String) obsoleteRGPInstance.getAttribute(
                        ReactomeJavaConstants.variantIdentifier);
                }

                if (variantIdentifier != null) {
                    continue;
                }
                obsoleteDbId = obsoleteRGPInstance.getDbId();
                speciesName = getSpeciesName(obsoleteRGPInstance);

                List<SimpleInstance> referrers =
                    curatorToolAPI.getReferrers(obsoleteRGPInstance, ReactomeJavaConstants.referenceEntity);
                for (SimpleInstance referrer : referrers) {
                    if (referrer.getSchemaClassName().equals(ReactomeJavaConstants.EntityWithAccessionedSequence)) {
                        SimpleInstance referrerStableIdInstance =
                            (SimpleInstance) referrer.getAttribute(ReactomeJavaConstants.stableIdentifier);
                        if (referrerStableIdInstance != null) {
                            String referrerStableId =
                                (String) referrerStableIdInstance.getAttribute(ReactomeJavaConstants.identifier);
                            referrerIds.add(referrerStableId);
                        } else {
                            referrerIds.add(referrer.getDbId().toString());
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
                    "|[https://newcurator.reactome.org/curatorgraph/dataSchema/DatabaseObject/instance/%d %d]\n",
                    obsoleteDbId, obsoleteDbId
                ));
                reportLineBuilder.append(String.format("||%s\n", String.join(
                    "|", referrerIds.stream().map(Object::toString).collect(Collectors.toList())
                )));
                reportLineBuilder.append(String.format("|%s\n", speciesName));
                reportLineBuilder.append("|-\n");

                String reportLine = reportLineBuilder.toString();
                if (skipList.stream().anyMatch(accession -> accession.equals(rgpAccession))) {
                    skipNoReplacementReportLines.add(reportLine);
                } else if (isPlantReportLine(reportLine)) {
                    plantNoReplacementReportLines.add(reportLine);
                } else {
                    wikiWriter.write(reportLine);
                }
            } else {
                noReferrerDbIds.add(obsoleteDbId);
            }
        }

        for (String isoformAccession : isoformAccessionToDbId.keySet()) {
            List<SimpleInstance> isoformInstances = curatorToolAPI.getReferenceIsoformByVariantIdentifier(isoformAccession);
            String speciesName;
            for (SimpleInstance isoformInstance : isoformInstances) {
                List<String> referrerIds = new ArrayList<>();
                long isoformInstanceDbId = isoformInstance.getDbId();
                speciesName = getSpeciesName(isoformInstance);

                List<SimpleInstance> referrers =
                    curatorToolAPI.getReferrers(isoformInstance, ReactomeJavaConstants.referenceEntity);
                for (SimpleInstance referrer : referrers) {
                    if (referrer.getSchemaClassName().equals(ReactomeJavaConstants.EntityWithAccessionedSequence)) {
                        SimpleInstance referrerStableIdInstance =
                            (SimpleInstance) referrer.getAttribute(ReactomeJavaConstants.stableIdentifier);
                        if (referrerStableIdInstance != null) {
                            String referrerStableId =
                                (String) referrerStableIdInstance.getAttribute(ReactomeJavaConstants.identifier);
                            referrerIds.add(referrerStableId);
                        } else {
                            referrerIds.add(referrer.getDbId().toString());
                        }
                    }
                }

                if (!referrerIds.isEmpty()) {
                    StringBuilder reportLineBuilder = new StringBuilder();
                    //reportLineBuilder.append("|\n");
                    reportLineBuilder.append(String.format("||%s\n", isoformAccession));
                    reportLineBuilder.append(String.format(
                        "|[https://newcurator.reactome.org/curatorgraph/dataSchema/DatabaseObject/instance/%d %d]\n",
                        isoformInstanceDbId, isoformInstanceDbId
                    ));
                    reportLineBuilder.append(String.format("||%s\n", String.join(
                        "|", referrerIds.stream().map(Object::toString).collect(Collectors.toList())
                    )));
                    reportLineBuilder.append(String.format("|%s\n", speciesName));
                    reportLineBuilder.append("|-\n");

                    String reportLine = reportLineBuilder.toString();
                    if (skipList.stream().anyMatch(accession -> accession.equals(isoformAccession))) {
                        skipNoReplacementReportLines.add(reportLine);
                    } else if (isPlantReportLine(reportLine)) {
                        plantNoReplacementReportLines.add(reportLine);
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
                "|+ PLANT Obsolete UniProt Instances (with replacement UniProt)\n" +
                "|-\n" +
                "! Replacement UniProt\n" +
                "! Obsolete UniProt\n" +
                "! Reactome instances with obsolete UniProt\n" +
                "! EWAS associated with obsolete UniProt\n" +
                "! Species\n" +
                "|-\n"
        );

        for (String plantReplaceableReportLine : plantReplaceableReportLines) {
            wikiWriter.append(plantReplaceableReportLine);
        }

        wikiWriter.write("|}\n-----\n");

        wikiWriter.write(
            "{| class=\"wikitable\"\n" +
                "|+ PLANT Obsolete UniProt Instances (deleted forever, no replacement)\n" +
                "|-\n" +
                "! Obsolete UniProt\n" +
                "! Reactome instances with obsolete UniProt\n" +
                "! EWAS associated with obsolete UniProt\n" +
                "! Species\n" +
                "|-\n"
        );

        for (String plantNoReplacementReportLine : plantNoReplacementReportLines) {
            wikiWriter.append(plantNoReplacementReportLine);
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

                //dba.deleteByDBID(noReferrerDbId);
                System.out.println("Deleting DBID: " + noReferrerDbId);
                curatorToolAPI.deleteByDbId(noReferrerDbId);
            }
        }

        System.out.println("Checking for duplicate isoform instances...");

        List<SimpleInstance> referenceIsoformUniProtInstances =
            curatorToolAPI.fetchUniProtReferenceIsoformInstances();

        Map<String,List<Long>> variantIdentifierToDbId = new HashMap<>();
        for (SimpleInstance referenceIsoformUniProtInstance : referenceIsoformUniProtInstances) {
            String variantIdentifier =
                (String) referenceIsoformUniProtInstance.getAttribute(ReactomeJavaConstants.variantIdentifier);
            long isoformDbId = referenceIsoformUniProtInstance.getDbId();

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

    private String getCurrentDate() {
        return ZonedDateTime.now().format(DateTimeFormatter.ofPattern("EEE MMM dd YYYY"));
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

    private List<String> getSkipList() {
        final BufferedReader skipListWithNoReplacement = getSkipListFileBufferedReader("skiplist_no_replacement.txt");
        final BufferedReader skipListWithReplacement = getSkipListFileBufferedReader("skiplist_with_replacement.txt");

        List<String> skipListIds = new ArrayList<>();
        skipListIds.addAll(skipListWithNoReplacement.lines().filter(this::isValidUniProtId).collect(Collectors.toList()));
        skipListIds.addAll(skipListWithReplacement.lines().filter(this::isValidUniProtId).collect(Collectors.toList()));
        return skipListIds;
    }

    private BufferedReader getSkipListFileBufferedReader(String skipListFileName) {
        return new BufferedReader(new InputStreamReader(
            this.getClass().getClassLoader().getResourceAsStream(skipListFileName)
        ));
    }

    private boolean isValidUniProtId(String potentialUniProtId) {
        final List<Integer> validUniProtIdLengths = Arrays.asList(6, 10);
        return validUniProtIdLengths.contains(potentialUniProtId.length());
    }

    private boolean isAReferenceIsoform(SimpleInstance rgpInstance) {
        return rgpInstance.getSchemaClassName().equals(ReactomeJavaConstants.ReferenceIsoform);
    }

    private SimpleInstance fetchReferenceDNASequenceByDbId(CuratorToolAPI curatorToolAPI, long referenceDNASequenceDbId) {
        return curatorToolAPI.findByDbId(referenceDNASequenceDbId);
    }

    private boolean sameDbId(SimpleInstance instance1, SimpleInstance instance2) {
        return instance1.getDbId().equals(instance2.getDbId());
    }

    private boolean areDifferentLists(List<?> list1, List<?> list2) {
        if (list1 == list2) {
            return false;
        } else if (list1 == null || list2 == null) {
            return true;
        }

        return !list1.equals(list2);
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
        CuratorToolAPI curatorToolAPI,
        SimpleInstance instance,
        Map<String, List<?>> values,
        BufferedWriter sequenceReportWriter
    ) throws Exception {

        boolean isInstanceChanged = false;
        for (String attributeName : values.keySet()) {
            List<?> newValuesForAttribute =
                values.get(attributeName).stream().filter(Objects::nonNull).collect(Collectors.toList());
            if (newValuesForAttribute.isEmpty()) {
                System.out.println("WARNING: No new values for " + attributeName + " on " + instance.getDbId() +
                    " skipping attribute update");
                continue;
            }

            if (attributeName.toLowerCase().equals(ReactomeJavaConstants.checksum)) {
                Boolean oldSequenceChangedValue = (Boolean) instance.getAttribute("isSequenceChanged");

                Boolean newSequenceChangedValue =
                    getNewIsSequenceChangedAttributeValue(instance, newValuesForAttribute);

                if (oldSequenceChangedValue == null || !oldSequenceChangedValue.equals(newSequenceChangedValue)) {
                    instance.setAttribute("isSequenceChanged", newSequenceChangedValue);
                    System.out.println(String.format("%s (%d) has a new is_sequence_changed value",
                        instance.getDisplayName(), instance.getDbId()));
                    isInstanceChanged = true;
                }
            }

            if (attributeName.toLowerCase().equals(ReactomeJavaConstants.chain)) {
                boolean chainChangeLogUpdated =
                    updateChainLog(instance, (List<String>) newValuesForAttribute, sequenceReportWriter);
                if (hasChains(instance) && chainChangeLogUpdated) {
                    List<SimpleInstance> ewasInstances = getAllEwasInstances(curatorToolAPI, instance);

                    for (SimpleInstance ewasInstance : ewasInstances) {
                        reportChangedChainForEWASInstance(instance, ewasInstance);
                    }
                }
            }

            if (valuesChanged(instance, attributeName, newValuesForAttribute)) {
                if (isSingleAttribute(attributeName)) {
                    instance.setAttribute(attributeName, newValuesForAttribute.get(0));
                } else {
                    instance.setAttribute(attributeName, newValuesForAttribute);
                }

                isInstanceChanged = true;
            }
        }

        if (isInstanceChanged) {
            curatorToolAPI.commit(instance);
        }
    }

    private Boolean getNewIsSequenceChangedAttributeValue(SimpleInstance instance, List<?> newValues) {
        String oldChecksum = (String) instance.getAttribute(ReactomeJavaConstants.checksum);
        String newChecksum = newValues.get(0).toString();

        return isSequenceChanged(oldChecksum, newChecksum);
    }

    private boolean isSequenceChanged(String oldChecksum, String newChecksum) {
        return oldChecksum != null && newChecksum != null && !oldChecksum.equals(newChecksum);
    }

    @SuppressWarnings("unchecked")
    private boolean updateChainLog(SimpleInstance instance, List<String> newChainValues, BufferedWriter sequenceReportWriter)
        throws Exception {
        boolean chainLogChanged = false;

        List<String> oldChainValues = (List<String>) instance.getAttribute(ReactomeJavaConstants.chain);
        String date = getCurrentDate();

        String referenceGeneProductDescription = getReferenceGeneProductDescription(instance);

        for (String oldChainValue : oldChainValues) {
            if (!newChainValues.contains(oldChainValue)) {
                String logEntry = String.format("%s for %d removed on %s", oldChainValue, instance.getDbId(), date);
                sequenceReportWriter.write(logEntry + " for " + referenceGeneProductDescription + "\n");

                String existingLog = (String) instance.getAttribute("_chainChangeLog");
                String fullLog =
                    existingLog != null ?
                    existingLog + ";" + logEntry :
                    logEntry;

                instance.setAttribute("_chainChangeLog", fullLog);
                System.out.println("old chain removed for " + instance.getDbId());
                chainLogChanged = true;
            }
        }

        for (String newChainValue : newChainValues) {
            if (!oldChainValues.contains(newChainValue)) {
                String logEntry = String.format("%s for %d added on %s", newChainValue, instance.getDbId(), date);
                sequenceReportWriter.write(logEntry + " for " + referenceGeneProductDescription + "\n");


                String existingLog = (String) instance.getAttribute("_chainChangeLog");
                String fullLog =
                    existingLog != null ?
                        existingLog + ";" + logEntry :
                        logEntry;

                instance.setAttribute("_chainChangeLog", fullLog);
                System.out.println("new chain added for " + instance.getDbId());
                chainLogChanged = true;
            }
        }
        return chainLogChanged;
    }

    private boolean hasChains(SimpleInstance instance) {
        List<String> chainValues = (List<String>) instance.getAttribute(ReactomeJavaConstants.chain);
        return chainValues != null && !chainValues.isEmpty();
    }

    private String getReferenceGeneProductDescription(SimpleInstance rgpInstance) {
        String referenceGeneProductDescription = rgpInstance.getDbId() != null ? rgpInstance.getDbId().toString() : "";

        String rgpName = (String) rgpInstance.getAttribute(ReactomeJavaConstants.name);
        if (rgpName != null && !rgpName.isEmpty()) {
            referenceGeneProductDescription += " - " + rgpName;
        }

        SimpleInstance speciesInstance = (SimpleInstance) rgpInstance.getAttribute(ReactomeJavaConstants.species);
        if (speciesInstance != null) {
            referenceGeneProductDescription += " (" + speciesInstance.getDisplayName() + ")";
        }

        return referenceGeneProductDescription;
    }

    private boolean valuesChanged(SimpleInstance instance, String attributeName, List<?> newValues) {
        List<String> currentValuesToCompare = toComparableValues(getAttributeValues(instance, attributeName));
        List<String> newValuesToCompare = toComparableValues(newValues);

        if (!areDifferentLists(currentValuesToCompare, newValuesToCompare)) {
            return false;
        }

        System.out.println(String.format("%s changed for instance %d", attributeName, instance.getDbId()));
        System.out.println(String.format("old attribute values - %s", String.join(",", currentValuesToCompare)));
        System.out.println(String.format("new attribute values - %s", String.join(",", newValuesToCompare)));

        return true;
    }

    /**
     * SimpleInstance holds a single-valued attribute as the value itself, a multi-valued attribute as a List, and
     * returns null for an attribute with no value, so every read is normalized to a list here.
     */
    private List<Object> getAttributeValues(SimpleInstance instance, String attributeName) {
        Object value = instance.getAttribute(attributeName);
        if (value == null) {
            return Collections.emptyList();
        }
        return value instanceof List ? new ArrayList<>((List<?>) value) : Collections.singletonList(value);
    }

    private List<String> toComparableValues(List<?> values) {
        return values.stream().map(this::toComparableValue).collect(Collectors.toList());
    }

    /**
     * SimpleInstance carries no schema, so there is nothing to ask whether the attribute is instance-typed; it is
     * decided per value instead. An instance value is compared by dbId and anything else by its string form, which
     * also keeps an Integer and a Long holding the same number from reading as a change.
     */
    private String toComparableValue(Object value) {
        return value instanceof SimpleInstance ?
            String.valueOf(((SimpleInstance) value).getDbId()) :
            value.toString();
    }

    private boolean isSingleAttribute(String attributeName) {
        return Arrays.asList(
            ReactomeJavaConstants.species,
            ReactomeJavaConstants.sequenceLength,
            ReactomeJavaConstants.checksum,
            ReactomeJavaConstants.comment
        ).contains(attributeName);
    }

    private String getSpeciesName(SimpleInstance instance) {
        SimpleInstance species = (SimpleInstance) instance.getAttribute(ReactomeJavaConstants.species);
        if (species != null) {
            return species.getDisplayName();
        }
        return "";
    }

    @SuppressWarnings("unchecked")
    private List<SimpleInstance> getRGPReferrers(CuratorToolAPI curatorToolAPI, SimpleInstance rgpInstance) throws Exception {
        List<SimpleInstance> referrers = new ArrayList<>();

        final List<String> reverseAttributes = Arrays.asList(
            ReactomeJavaConstants.referenceEntity,
            ReactomeJavaConstants.referenceSequence,
            ReactomeJavaConstants.secondReferenceSequence,
            ReactomeJavaConstants.isoformParent
        );

        for (String reverseAttribute : reverseAttributes) {
            List<SimpleInstance> reverseAttributeReferrers =
                curatorToolAPI.getReferrers(rgpInstance, reverseAttribute);
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

    /**
     * Returns null when the entry has no parsable sequence length. updateInstance skips a null attribute value with a
     * warning, so one malformed entry costs its sequence length rather than failing the whole run.
     *
     * @param entry - the UniProt XML entry to parse.
     * @param accession - the entry's primary accession, for the warning message.
     * @return the sequence length, or null if the entry has none.
     */
    private Integer parseSequenceLength(String entry, String accession) {
        String sequenceLength = matchSingleValue(entry, "<sequence.*length=\"(\\d+)\"");
        if (sequenceLength.isEmpty()) {
            System.out.println("WARNING: No sequence length found for " + accession);
            return null;
        }
        return Integer.valueOf(sequenceLength);
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
        String position = matchSingleValue(featureContent, "<location.*?>\\s+<position position=\"(\\d+)\"");
        return "initiator methionine:" + position;
    }

    private String parseGenericChain(String featureType, String featureContent) {
        String chainStart = matchSingleValue(featureContent, "<begin position=\"(\\d+)\"");
        String chainEnd = matchSingleValue(featureContent, "<end position=\"(\\d+)");

        return featureType + ":" + chainStart + "-" + chainEnd;
    }



    private List<SimpleInstance> getAllEwasInstances(CuratorToolAPI curatorToolAPI, SimpleInstance referenceGeneProduct) throws Exception {
        List<SimpleInstance> allEwasInstances = new ArrayList<>();

        List<SimpleInstance> referenceEntityEwasInstances =
            curatorToolAPI.getReferrers(referenceGeneProduct, ReactomeJavaConstants.referenceEntity);
        if (referenceEntityEwasInstances != null) {
            allEwasInstances.addAll(referenceEntityEwasInstances);
        }
        List<SimpleInstance> hasModifiedResidueInstances = new ArrayList<>();
        List<SimpleInstance> referenceSequenceModifiedResidues =
            curatorToolAPI.getReferrers(referenceGeneProduct, ReactomeJavaConstants.referenceSequence);

        if (referenceSequenceModifiedResidues != null) {
            hasModifiedResidueInstances.addAll(referenceSequenceModifiedResidues);
        }

        List<SimpleInstance> secondReferenceSequenceModifiedResidues =
            curatorToolAPI.getReferrers(referenceGeneProduct, ReactomeJavaConstants.secondReferenceSequence);
        if (secondReferenceSequenceModifiedResidues != null) {
            hasModifiedResidueInstances.addAll(secondReferenceSequenceModifiedResidues);
        }

        for (SimpleInstance hasModifiedResidueInstance : hasModifiedResidueInstances) {
            List<SimpleInstance> hasModifiedEwasInstances =
                curatorToolAPI.getReferrers(hasModifiedResidueInstance, ReactomeJavaConstants.hasModifiedResidue);
            if (hasModifiedEwasInstances != null) {
                allEwasInstances.addAll(hasModifiedEwasInstances);
            }
        }

        return allEwasInstances;
    }

    private void reportChangedChainForEWASInstance(SimpleInstance referenceGeneProduct, SimpleInstance ewas)
        throws Exception {
        //"RGP db id\tRGP Accession\tEWAS db id\tEWAS name\tEWAS author (created or last modified)\t"
        String reportLine = String.join("\t",
            referenceGeneProduct.getDbId().toString(),
            (String) referenceGeneProduct.getAttribute(ReactomeJavaConstants.identifier),
            ewas.getDbId().toString(),
            ewas.getDisplayName(),
            getAuthor(ewas)
        ).concat(System.lineSeparator());

        Files.write(
            getUniprotUpdateDirectoryPath().resolve("ewasCoordinatesReport.txt"),
            reportLine.getBytes(),
            StandardOpenOption.CREATE, StandardOpenOption.APPEND
        );
    }

    private String getAuthor(SimpleInstance ewas) throws Exception {
        SimpleInstance ewasCreatedInstanceEdit = (SimpleInstance) ewas.getAttribute(ReactomeJavaConstants.created);
        if (ewasCreatedInstanceEdit != null) {
            return getAuthorFromInstanceEdit(ewasCreatedInstanceEdit);
        } else {
            return getLastModifiedAuthor(ewas);
        }
    }

    private String getLastModifiedAuthor(SimpleInstance ewas) throws Exception {
        List<SimpleInstance> ewasModifiedInstanceEdits = (List<SimpleInstance>) ewas.getAttribute(ReactomeJavaConstants.modified);
        if (ewasModifiedInstanceEdits != null && !ewasModifiedInstanceEdits.isEmpty()) {
            SimpleInstance ewasMostRecentModifiedInstanceEdit = ewasModifiedInstanceEdits.get(0);
            return getAuthorFromInstanceEdit(ewasMostRecentModifiedInstanceEdit);
        } else {
            return "Unknown author";
        }
    }

    private String getAuthorFromInstanceEdit(SimpleInstance instanceEdit) throws Exception {
        SimpleInstance instanceEditAuthor = (SimpleInstance) instanceEdit.getAttribute(ReactomeJavaConstants.author);
        if (instanceEditAuthor != null) {
            return instanceEditAuthor.getDisplayName();
        } else {
            return "Unknown author";
        }
    }

    private boolean isPlantReportLine(String reportLine) {
        final List<String> plantSpeciesNames = Arrays.asList("Arabidopsis thaliana", "Oryza sativa");
        return plantSpeciesNames.stream().anyMatch(reportLine::contains);
    }
}
