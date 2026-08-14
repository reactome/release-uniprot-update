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
import static org.reactome.release.utils.Utils.getAttributeValues;
import static org.reactome.release.utils.Utils.getTrEMBLIds;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 *         Created 7/31/2023
 */
public class Main {
    private Path uniprotUpdateDirectoryPath;

    // Held as fields, rather than as locals of run, so that close can release them whether or not the run completed:
    // an unclosed report writer loses whatever is still sitting in its buffer, and the curator tool API holds a Spring
    // context whose threads keep the JVM alive.
    private CuratorToolAPI curatorToolAPI;
    private BufferedReader swissProtFileReader;
    private BufferedWriter sequenceReportWriter;
    private BufferedWriter referenceDNASequenceReportWriter;
    private BufferedWriter wikiWriter;

    private SimpleInstance humanEnsEMBLGeneReferenceDatabase;

    public static void main(String[] args) throws Exception {
        Main main = new Main();

        String configFilePathAsString = args.length > 0 ? args[0] : getDefaultConfigFilePath().toString();
        Properties configProperties = getConfigProperties(configFilePathAsString);

        try {
            main.run(configProperties);
        } finally {
            main.close();
        }
    }

    /**
     * Closes everything the run opened, whether or not it completed. Each resource is closed independently so that one
     * failure to close does not keep the others open or replace the exception that ended the run.
     */
    private void close() {
        closeQuietly(swissProtFileReader, "SwissProt file reader");
        closeQuietly(sequenceReportWriter, "sequence report writer");
        closeQuietly(referenceDNASequenceReportWriter, "reference DNA sequence report writer");
        closeQuietly(wikiWriter, "wiki report writer");

        if (curatorToolAPI != null) {
            curatorToolAPI.close();
        }
    }

    private void closeQuietly(Closeable resource, String resourceDescription) {
        if (resource == null) {
            return;
        }

        try {
            resource.close();
        } catch (IOException e) {
            System.err.println("Unable to close the " + resourceDescription + ": " + e.getMessage());
        }
    }

    @SuppressWarnings("unchecked")
    private void run(Properties configProperties) throws Exception {
        curatorToolAPI = new CuratorToolAPI(Long.parseLong(configProperties.getProperty("personId")));

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
        System.out.println("Populating rgp identifier to instance...");
        Map<String, List<SimpleInstance>> rgpIdentifierToInstances = curatorToolAPI.getRGPIdentifierToInstancesMap();
        System.out.println("Populating isoform accession to instance...");
        Map<String, List<SimpleInstance>> isoformAccessionToInstances =
            curatorToolAPI.getIsoformAccessionToInstancesMap();
        // Held apart from the index above, which is what the run looks isoforms up in and so must keep every isoform
        // of the database: this is the set of accessions still to be accounted for, which the run empties as the
        // SwissProt file turns out to carry them.
        Set<String> remainingIsoformAccessions =
            curatorToolAPI.getUniProtIsoformAccessions(isoformAccessionToInstances);
        System.out.println("Populating rds identifier to instance...");
        Map<String, SimpleInstance> rdsIdentifierToInstance = curatorToolAPI.getRDSIdentifierToInstanceMap();

        Map<String, List<String>> secondaryAccessionToPrimaryAccessionList = new HashMap<>();
        Map<String, String> misMatchedIsoformAccessionToRGPAccession = new HashMap<>();

        sequenceReportWriter = Files.newBufferedWriter(
            getUniprotUpdateDirectoryPath().resolve("sequence_uniprot_report.txt"));
        referenceDNASequenceReportWriter = Files.newBufferedWriter(
            getUniprotUpdateDirectoryPath().resolve("reference_DNA_sequence_report.txt"));

        String line;
        StringBuilder entryBuilder = new StringBuilder();

        int recordCounter = 0;

        SwissProtFileProcessor swissProtFileProcessor = new SwissProtFileProcessor(getUniprotUpdateDirectoryPath());
        swissProtFileReader = swissProtFileProcessor.getFileReader();
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
                if (accessions.isEmpty()) {
                    System.out.println("WARNING: No accession in record " + recordCounter + " -- skipping it");
                    continue;
                }
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
                        // The first matching species name wins; without stopping here a later match would replace it
                        // and each match would cost another species query.
                        break;
                    }
                }

                if (speciesInstance == null && !rgpAccessionToDbId.containsKey(primaryAccession)) {
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

                // [^>]* keeps the match inside the one <sequence> tag: an entry is read as a single line, so a greedy
                // .* here would run past the tag and take a checksum from anywhere later in the entry.
                String checksum = matchSingleValue(entry, "<sequence[^>]*checksum=\"([0-9A-F]+)\"");

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

                    for (String ensEMBLGeneId : uniqueEnsEMBLGeneIds) {
                        SimpleInstance referenceDNASequence = refreshIndexedInstance(
                            rdsIdentifierToInstance, ensEMBLGeneId);

                        if (referenceDNASequence != null) {
                            referenceDNASequenceReportWriter.write("Checking existing reference DNA sequence for " +
                                ensEMBLGeneId + " with db_id " + referenceDNASequence.getDbId() + "\n");

                            SimpleInstance existingRDSReferenceDatabase = (SimpleInstance)
                                referenceDNASequence.getAttribute(ReactomeJavaConstants.referenceDatabase);
                            boolean isUpdateToReferenceDNASequence = false;
                            if (existingRDSReferenceDatabase == null ||
                                !sameDbId(existingRDSReferenceDatabase, getHumanEnsEMBLGeneReferenceDatabase())) {
                                referenceDNASequence.setAttribute(
                                    ReactomeJavaConstants.referenceDatabase, getHumanEnsEMBLGeneReferenceDatabase());
                                isUpdateToReferenceDNASequence = true;
                            }

                            List<Object> existingGeneNames =
                                getAttributeValues(referenceDNASequence, ReactomeJavaConstants.geneName);
                            if (areDifferentLists(existingGeneNames, geneNames)) {
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
                                    referenceDNASequence.getDbId() + "\n"
                                );
                                referenceDNASequence.setDisplayName(
                                    curatorToolAPI.getReferenceSequenceDisplayName(referenceDNASequence));
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

                            referenceDNASequence.setDisplayName(
                                curatorToolAPI.getReferenceSequenceDisplayName(referenceDNASequence));

                            long referenceDNASequenceDbId = curatorToolAPI.commit(referenceDNASequence).getDbId();
                            referenceDNASequenceReportWriter.write("Reference DNA sequence with db_id " +
                                referenceDNASequenceDbId + " created for " + ensEMBLGeneId + "\n");
                            rdsIdentifierToInstance.put(ensEMBLGeneId, referenceDNASequence);
                        }
                        referenceDNASequences.add(referenceDNASequence);
                    }
                }
                List<String> keywords = matchMultipleValues(entry, "<keyword id=\".*?\">(.*?)</keyword>");

                List<String> comments = parseComments(entry);

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
                values.put(ReactomeJavaConstants.comment, comments);
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
                    // updateInstance commits this instance a second time, and the commit above has left the copy
                    // held here out of date with what is stored (its InstanceEdits in particular), so it is re-read
                    // before being written again.
                    newReferenceGeneProductInstance = curatorToolAPI.refresh(newReferenceGeneProductInstance);

                    System.out.println(String.format("New UniProt:%s\t%d", primaryAccession, newRGPDbId));
                    updateInstance(curatorToolAPI, newReferenceGeneProductInstance, values, sequenceReportWriter);
                    for (String isoformId : isoformIds) {
                        if (!isoformId.contains(primaryAccession)) {
                            // Recorded for the mis-matched isoform clean-up below and not created here, matching how
                            // the branch for an existing ReferenceGeneProduct handles a mis-matched isoform.
                            misMatchedIsoformAccessionToRGPAccession.put(isoformId, primaryAccession);
                            continue;
                        }

                        SimpleInstance newIsoformInstance = new SimpleInstance();
                        newIsoformInstance.setSchemaClassName(ReactomeJavaConstants.ReferenceIsoform);
                        newIsoformInstance.setAttribute(
                            ReactomeJavaConstants.referenceDatabase, uniProtReferenceDatabase);
                        newIsoformInstance.setAttribute(ReactomeJavaConstants.identifier, primaryAccession);
                        newIsoformInstance.setAttribute(ReactomeJavaConstants.isoformParent,
                            Collections.singletonList(newReferenceGeneProductInstance));
                        newIsoformInstance.setAttribute(ReactomeJavaConstants.variantIdentifier, isoformId);

                        updateInstance(curatorToolAPI, newIsoformInstance, values, sequenceReportWriter);
                        indexIsoform(isoformAccessionToInstances, isoformId, newIsoformInstance);
                    }
                } else {
                    Collection<SimpleInstance> existingReferenceGeneProductInstances =
                        refreshIndexedInstances(rgpIdentifierToInstances, primaryAccession);
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

                        updateInstance(
                            curatorToolAPI, existingReferenceGeneProductInstance, values, sequenceReportWriter);

                        duplicateFlag = true;

                        // The species values are a singleton list holding whatever the file gave, so the list is never
                        // empty -- it holds a null when no species name matched. That is the case where the existing
                        // instance's species is carried over to the isoform updates below.
                        if (speciesInstance == null) {
                            values.put(ReactomeJavaConstants.species, Collections.singletonList((SimpleInstance)
                                existingReferenceGeneProductInstance.getAttribute(ReactomeJavaConstants.species))
                            );
                        }
                        for (String isoformId : isoformIds) {
                            if (isoformId.contains(primaryAccession)) {
                                List<SimpleInstance> isoformInstances =
                                    refreshIndexedInstances(isoformAccessionToInstances, isoformId);
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
                                            Collections.singletonList(existingReferenceGeneProductInstance));

                                        updateInstance(curatorToolAPI, isoformInstance, values, sequenceReportWriter);

                                        remainingIsoformAccessions.remove(isoformId);
                                    }
                                } else {
                                    SimpleInstance isoformInstance = new SimpleInstance();
                                    isoformInstance.setSchemaClassName(ReactomeJavaConstants.ReferenceIsoform);
                                    isoformInstance.setAttribute(ReactomeJavaConstants.identifier,
                                        primaryAccession);
                                    isoformInstance.setAttribute(ReactomeJavaConstants.isoformParent,
                                        Collections.singletonList(existingReferenceGeneProductInstance));
                                    isoformInstance.setAttribute(ReactomeJavaConstants.variantIdentifier,
                                        isoformId);
                                    long isoformDbId = curatorToolAPI.commit(isoformInstance).getDbId();
                                    // As with a new ReferenceGeneProduct above, updateInstance commits this instance
                                    // a second time, so it is re-read first.
                                    isoformInstance = curatorToolAPI.refresh(isoformInstance);

                                    System.out.println(String.format("New isoform: %s\t%d\tMaster: %d",
                                        isoformId, isoformDbId, existingReferenceGeneProductInstance.getDbId()));

                                    updateInstance(curatorToolAPI, isoformInstance, values, sequenceReportWriter);
                                    indexIsoform(isoformAccessionToInstances, isoformId, isoformInstance);
                                }
                            } else {
                                misMatchedIsoformAccessionToRGPAccession.put(isoformId, primaryAccession);
                            }
                        }
                    }
                    // Removed once the entry has been processed, not inside the loop above: the loop skips
                    // ReferenceIsoforms, so an accession whose only instances are isoforms never reached this and was
                    // then reported as obsolete despite the file still carrying it.
                    rgpAccessionToDbId.remove(primaryAccession);
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

            List<SimpleInstance> isoformInstances =
                refreshIndexedInstances(isoformAccessionToInstances, misMatchedIsoformAccession);

            SimpleInstance isoformInstance = !isoformInstances.isEmpty() ? isoformInstances.get(0) : null;
            if (isoformInstance != null) {
                List<Object> existingParents =
                    getAttributeValues(isoformInstance, ReactomeJavaConstants.isoformParent);
                if (existingParents.isEmpty()) {
                    continue;
                }
                // All of the existing parents are kept: the commit below replaces the attribute's values, so any
                // parent left out here would be dropped from the instance.
                existingParents.forEach(existingParent -> isoformParents.add((SimpleInstance) existingParent));
            }

            List<SimpleInstance> mismatchedParents =
                refreshIndexedInstances(rgpIdentifierToInstances, misMatchedIsoformAccession);

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
        curatorToolAPI.updateReferenceIsoformDisplayNames();

        System.out.println("Done");

        System.out.println("Remaining instances:" + rgpAccessionToDbId.keySet().size());

        System.out.println("Deleting obsolete instances with no referrers...");

        // Looked up for every remaining accession at once, rather than accession by accession inside the loop: each
        // look-up is a request to UniProt, and one request answers for a hundred accessions.
        Set<String> tremblAccessionSet = getTrEMBLIds(rgpAccessionToDbId.keySet());

        Iterator<String> rgpAccessionsIterator = rgpAccessionToDbId.keySet().iterator();
        List<String> tremblAccessions = new ArrayList<>();
        while (rgpAccessionsIterator.hasNext()) {
            String rgpAccession = rgpAccessionsIterator.next();

            if (tremblAccessionSet.contains(rgpAccession)) {
                tremblAccessions.add(rgpAccession);
                rgpAccessionsIterator.remove();
            } else {
                List<SimpleInstance> obsoleteReferenceGeneProductInstances =
                    refreshIndexedInstances(rgpIdentifierToInstances, rgpAccession);

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
                    if (hasNoRGPReferrers(curatorToolAPI, obsoleteReferenceGeneProductInstance)) {
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
        Reportable trEMBLAccessionReport =
            new TrEMBLAccessionReport(getUniprotUpdateDirectoryPath(), tremblAccessions);
        trEMBLAccessionReport.writeReport();

        List<Long> dbIdsToSkip = new ArrayList<>();
        Iterator<String> isoformAccessionIterator = remainingIsoformAccessions.iterator();
        while (isoformAccessionIterator.hasNext()) {
            String isoformAccession = isoformAccessionIterator.next();
            List<SimpleInstance> isoformInstances =
                refreshIndexedInstances(isoformAccessionToInstances, isoformAccession);

            SimpleInstance isoformInstance = !isoformInstances.isEmpty() ? isoformInstances.get(0) : null;
            if (isoformInstance == null) {
                System.out.println(isoformAccession + " is not a variant identifier for any ReferenceIsoform");
                continue;
            }

            long obsoleteIsoformDbId = isoformInstance.getDbId();
            List<Object> isoformParents =
                getAttributeValues(isoformInstance, ReactomeJavaConstants.isoformParent);
            if (isoformParents.isEmpty()) {
                System.out.println(isoformInstance.getDbId());
                dbIdsToSkip.add(obsoleteIsoformDbId);
                continue;
            }

            SimpleInstance isoformParent = (SimpleInstance) isoformParents.get(0);
            String isoformParentIdentifier = (String)
                isoformParent.getAttribute(ReactomeJavaConstants.identifier);
            if (isoformParentIdentifier == null || isoformParentIdentifier.isEmpty()) {
                continue;
            }

            if (hasNoRGPReferrers(curatorToolAPI, isoformInstance)) {
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

        wikiWriter = Files.newBufferedWriter(getUniprotUpdateDirectoryPath().resolve("uniprot.wiki"));

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
                isSecondaryAccession = true;

                List<SimpleInstance> obsoleteRGPInstances =
                    refreshIndexedInstances(rgpIdentifierToInstances, rgpAccession);
                for (SimpleInstance obsoleteRGPInstance : obsoleteRGPInstances) {
                    String variantIdentifier = null;
                    if (isAReferenceIsoform(obsoleteRGPInstance)) {
                        variantIdentifier = (String) obsoleteRGPInstance.getAttribute(
                            ReactomeJavaConstants.variantIdentifier);
                    }

                    if (variantIdentifier != null) {
                        continue;
                    }
                    // Declared per instance: two instances can share an accession, and each one gets its own report
                    // line, so referrers must not carry over from the instance before it.
                    long obsoleteDbId = obsoleteRGPInstance.getDbId();
                    List<Long> referrerDbIds = new ArrayList<>();
                    String speciesName = getSpeciesName(obsoleteRGPInstance);

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
                    } else if (hasNoRGPReferrers(curatorToolAPI, obsoleteRGPInstance)) {
                        // The report lists EWAS referrers only, but the instances collected here are deleted further
                        // down, so an instance is only added once nothing at all refers to it -- the same check the
                        // first round of deletions makes.
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
            System.out.println(rgpAccession);

            List<SimpleInstance> obsoleteRGPInstances =
                refreshIndexedInstances(rgpIdentifierToInstances, rgpAccession);

            for (SimpleInstance obsoleteRGPInstance : obsoleteRGPInstances) {
                String variantIdentifier = null;
                if (isAReferenceIsoform(obsoleteRGPInstance)) {
                    variantIdentifier = (String) obsoleteRGPInstance.getAttribute(
                        ReactomeJavaConstants.variantIdentifier);
                }

                if (variantIdentifier != null) {
                    continue;
                }
                // Declared per instance: two instances can share an accession, and each one gets its own report
                // line, so referrers must not carry over from the instance before it. Declaring them here also
                // means no report line and no deletion for an accession whose instances are all isoforms, where
                // there is no obsolete dbId to report or delete.
                long obsoleteDbId = obsoleteRGPInstance.getDbId();
                List<String> referrerIds = new ArrayList<>();
                String speciesName = getSpeciesName(obsoleteRGPInstance);

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
                } else if (hasNoRGPReferrers(curatorToolAPI, obsoleteRGPInstance)) {
                    // The report lists EWAS referrers only, but the instances collected here are deleted further
                    // down, so an instance is only added once nothing at all refers to it -- the same check the
                    // first round of deletions makes.
                    noReferrerDbIds.add(obsoleteDbId);
                }
            }
        }

        for (String isoformAccession : remainingIsoformAccessions) {
            List<SimpleInstance> isoformInstances =
                refreshIndexedInstances(isoformAccessionToInstances, isoformAccession);
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
                } else if (hasNoRGPReferrers(curatorToolAPI, isoformInstance)) {
                    // The report lists EWAS referrers only, but the instances collected here are deleted further
                    // down, so an instance is only added once nothing at all refers to it -- the same check the
                    // first round of deletions makes.
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

        // The skip check is a condition on the delete, not a loop around it: nesting the delete inside a loop over
        // dbIdsToSkip deleted nothing at all when that list was empty, and deleted each dbId once per non-matching
        // entry when it was not.
        for (Long noReferrerDbId : noReferrerDbIds) {
            if (dbIdsToSkip.contains(noReferrerDbId)) {
                continue;
            }

            System.out.println("Deleting DBID: " + noReferrerDbId);
            curatorToolAPI.deleteByDbId(noReferrerDbId);
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

    private SimpleInstance getHumanEnsEMBLGeneReferenceDatabase() {
        if (humanEnsEMBLGeneReferenceDatabase == null) {
            humanEnsEMBLGeneReferenceDatabase = curatorToolAPI.getHumanEnsEMBLGeneReferenceDatabase();
        }
        return humanEnsEMBLGeneReferenceDatabase;
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
        return ZonedDateTime.now().format(DateTimeFormatter.ofPattern("EEE MMM dd yyyy"));
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

    private List<String> getSkipList() throws IOException {
        List<String> skipListIds = new ArrayList<>();

        try (
            BufferedReader skipListWithNoReplacement = getSkipListFileBufferedReader("skiplist_no_replacement.txt");
            BufferedReader skipListWithReplacement = getSkipListFileBufferedReader("skiplist_with_replacement.txt")
        ) {
            skipListIds.addAll(
                skipListWithNoReplacement.lines().filter(this::isValidUniProtId).collect(Collectors.toList()));
            skipListIds.addAll(
                skipListWithReplacement.lines().filter(this::isValidUniProtId).collect(Collectors.toList()));
        }
        return skipListIds;
    }

    private BufferedReader getSkipListFileBufferedReader(String skipListFileName) {
        InputStream skipListFileInputStream = this.getClass().getClassLoader().getResourceAsStream(skipListFileName);
        // Reported here rather than as the NullPointerException an absent resource would otherwise cause inside
        // InputStreamReader.
        if (skipListFileInputStream == null) {
            throw new IllegalStateException("The skip list file " + skipListFileName + " is not on the class path");
        }

        return new BufferedReader(new InputStreamReader(skipListFileInputStream));
    }

    private boolean isValidUniProtId(String potentialUniProtId) {
        final List<Integer> validUniProtIdLengths = Arrays.asList(6, 10);
        return validUniProtIdLengths.contains(potentialUniProtId.length());
    }

    private boolean isAReferenceIsoform(SimpleInstance rgpInstance) {
        return rgpInstance.getSchemaClassName().equals(ReactomeJavaConstants.ReferenceIsoform);
    }

    /**
     * Returns the indexed instance for the key, up to date with the database, or null if the index holds none for the
     * key or this run has deleted the one it held. The up-to-date copy takes the place of the indexed one, so that
     * the next look-up of the same key does not have to read it back again.
     *
     * @param index - instances of the database, indexed by the identifier they are looked up by.
     * @param key - the identifier to look up.
     * @return the instance for the key, or null if there is none.
     */
    private SimpleInstance refreshIndexedInstance(Map<String, SimpleInstance> index, String key) {
        SimpleInstance indexedInstance = index.get(key);
        if (indexedInstance == null) {
            return null;
        }

        SimpleInstance refreshedInstance = curatorToolAPI.refresh(indexedInstance);
        if (refreshedInstance == null) {
            index.remove(key);
        } else {
            index.put(key, refreshedInstance);
        }
        return refreshedInstance;
    }

    /**
     * As refreshIndexedInstance, for an index whose key can hold more than one instance. Instances this run has
     * deleted are dropped, so the list returned is empty rather than null where the key has nothing left behind it.
     *
     * @param index - instances of the database, indexed by the identifier they are looked up by.
     * @param key - the identifier to look up.
     * @return the instances for the key, empty if there are none.
     */
    private List<SimpleInstance> refreshIndexedInstances(Map<String, List<SimpleInstance>> index, String key) {
        List<SimpleInstance> indexedInstances = index.get(key);
        if (indexedInstances == null) {
            return Collections.emptyList();
        }

        indexedInstances.replaceAll(curatorToolAPI::refresh);
        indexedInstances.removeIf(Objects::isNull);
        return indexedInstances;
    }

    /**
     * Adds an isoform created by this run to the index, so that it is found by a later look-up of its variant
     * identifier -- the mis-matched isoform clean-up looks up isoform ids belonging to other entries, which the entry
     * they belong to may have created by then.
     *
     * An isoform with no dbId was never stored, so it is left out: the index is what the run treats as the contents
     * of the database, and committing such an instance from a later look-up would store it a second time.
     *
     * @param index - the isoforms of the database, indexed by variant identifier.
     * @param variantIdentifier - the variant identifier of the created isoform.
     * @param isoformInstance - the created isoform.
     */
    private void indexIsoform(
        Map<String, List<SimpleInstance>> index, String variantIdentifier, SimpleInstance isoformInstance) {
        if (isoformInstance.getDbId() == null) {
            return;
        }

        index.computeIfAbsent(variantIdentifier, k -> new ArrayList<>()).add(isoformInstance);
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
        if (ensEMBLIdData == null) {
            System.err.println("Unable to query EnsEMBL for " + ensEMBLGeneId + " after " + maxQueryAttempts +
                " attempts -- treating it as not on the primary assembly");
            return false;
        }

        // The capture group is inside the quotes so that the region name is compared with the unquoted values in
        // primaryAssemblyRegions.
        Pattern seqRegionPattern = Pattern.compile("\"seq_region_name\":\"(.*?)\"");
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
                System.out.println(String.format(
                    "Bad request for %s:  Sleeping for 5 seconds and retrying", ensemblLookupURL));
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

    /**
     * Records each chain value added to, or removed from, the instance in the sequence report and reports whether
     * there was anything to record.
     *
     * The entries used to also be appended to the instance's own "_chainChangeLog" attribute. No class in the
     * graph-core model hierarchy has that slot, so curator-tool-ws could not find a set method for it and dropped the
     * value on every commit (logging "Cannot find a set method for _chainChangeLog" as it did so): the log has only
     * ever reached the sequence report. Persisting it again means adding the attribute to graph-core first.
     *
     * @param instance - the instance whose chain values are being replaced.
     * @param newChainValues - the chain values parsed from the UniProt entry.
     * @param sequenceReportWriter - the writer for the sequence report.
     * @return true if any chain value was added or removed, false if the chain values are unchanged.
     */
    @SuppressWarnings("unchecked")
    private boolean updateChainLog(
        SimpleInstance instance, List<String> newChainValues, BufferedWriter sequenceReportWriter
    ) throws Exception {
        boolean chainLogChanged = false;

        // getAttribute returns null for an instance with no chain values (e.g. one created by this run), so the
        // list is defaulted here rather than dereferenced below.
        List<String> oldChainValues =
            emptyListIfNull((List<String>) instance.getAttribute(ReactomeJavaConstants.chain));
        String date = getCurrentDate();

        String referenceGeneProductDescription = getReferenceGeneProductDescription(instance);

        for (String oldChainValue : oldChainValues) {
            if (!newChainValues.contains(oldChainValue)) {
                String logEntry = String.format("%s for %d removed on %s", oldChainValue, instance.getDbId(), date);
                sequenceReportWriter.write(logEntry + " for " + referenceGeneProductDescription + "\n");

                System.out.println("old chain removed for " + instance.getDbId());
                chainLogChanged = true;
            }
        }

        for (String newChainValue : newChainValues) {
            if (!oldChainValues.contains(newChainValue)) {
                String logEntry = String.format("%s for %d added on %s", newChainValue, instance.getDbId(), date);
                sequenceReportWriter.write(logEntry + " for " + referenceGeneProductDescription + "\n");

                System.out.println("new chain added for " + instance.getDbId());
                chainLogChanged = true;
            }
        }
        return chainLogChanged;
    }

    private boolean hasChains(SimpleInstance instance) {
        return !getAttributeValues(instance, ReactomeJavaConstants.chain).isEmpty();
    }

    private String getReferenceGeneProductDescription(SimpleInstance rgpInstance) {
        String referenceGeneProductDescription = rgpInstance.getDbId() != null ? rgpInstance.getDbId().toString() : "";

        List<Object> rgpNames = getAttributeValues(rgpInstance, ReactomeJavaConstants.name);
        if (!rgpNames.isEmpty()) {
            referenceGeneProductDescription += " - " + rgpNames.get(0);
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

    /**
     * The attributes whose value is set as the value itself rather than as a List. These are the single-valued
     * attributes of the graph model: setting a multi-valued one to a bare value leaves curator-tool-ws unable to find
     * a set method for it (its setter takes a List), so it logs the miss and drops the value silently on commit.
     *
     * @param attributeName - the name of the attribute being set.
     * @return true if the attribute holds a single value, false if it holds a List.
     */
    private boolean isSingleAttribute(String attributeName) {
        return Arrays.asList(
            ReactomeJavaConstants.species,
            ReactomeJavaConstants.sequenceLength,
            ReactomeJavaConstants.checksum
        ).contains(attributeName);
    }

    private String getSpeciesName(SimpleInstance instance) {
        SimpleInstance species = (SimpleInstance) instance.getAttribute(ReactomeJavaConstants.species);
        if (species != null) {
            return species.getDisplayName();
        }
        return "";
    }

    /**
     * Reports whether anything in the database refers to the instance through any of the attributes a reference
     * gene product is referred to by. Every caller asks only whether there are referrers at all, which is why this
     * answers that rather than handing the referrers themselves back.
     *
     * The referrers are fetched once for all of the attributes. The single-attribute form of getReferrers fetches
     * every referrer of the instance and then keeps one attribute's worth of them, so calling it for each attribute
     * in turn fetched the same referrers over again for each one.
     *
     * @param curatorToolAPI - the API to fetch the referrers through.
     * @param rgpInstance - the instance to look for referrers of.
     * @return true if nothing refers to the instance, false if anything does.
     */
    private boolean hasNoRGPReferrers(CuratorToolAPI curatorToolAPI, SimpleInstance rgpInstance) throws Exception {
        final List<String> reverseAttributes = Arrays.asList(
            ReactomeJavaConstants.referenceEntity,
            ReactomeJavaConstants.referenceSequence,
            ReactomeJavaConstants.secondReferenceSequence,
            ReactomeJavaConstants.isoformParent
        );

        return emptyListIfNull(curatorToolAPI.getReferrers(rgpInstance))
            .stream()
            .filter(namedReferrers -> reverseAttributes.contains(namedReferrers.getAttributeName()))
            .allMatch(namedReferrers -> emptyListIfNull(namedReferrers.getReferrers()).isEmpty());
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
        // [^>]* rather than .*, for the reason given where the checksum is parsed.
        String sequenceLength = matchSingleValue(entry, "<sequence[^>]*length=\"(\\d+)\"");
        if (sequenceLength.isEmpty()) {
            System.out.println("WARNING: No sequence length found for " + accession);
            return null;
        }
        return Integer.valueOf(sequenceLength);
    }

    private List<String> parseComments(String entry) {
        // DOTALL so that a <text> body wrapped over several lines is captured whole rather than skipped. It also lets
        // the gap before <text> span lines, so that gap is guarded against </?comment: a comment type that owns no
        // <text> of its own (interaction, alternative products) would otherwise reach forward and take the text of a
        // later comment, filing it under the wrong type.
        Pattern commentsPattern = Pattern.compile(
            "<comment type=\"([A-Za-z ]*?)\"(?:(?!</?comment).)*?<text[^>]*>(.*?)</text>",
            Pattern.MULTILINE | Pattern.DOTALL);
        Matcher commentsMatcher = commentsPattern.matcher(entry);

        // Joined with a space rather than appended straight onto each other, which ran the end of one comment into the
        // type of the next ("FUNCTION ...SUBUNIT ...").
        List<String> comments = new ArrayList<>();
        while (commentsMatcher.find()) {
            String commentType = commentsMatcher.group(1).toUpperCase();
            String commentText = commentsMatcher.group(2);

            comments.add(commentType + " " + commentText);
        }
        return comments;
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



    private List<SimpleInstance> getAllEwasInstances(
        CuratorToolAPI curatorToolAPI, SimpleInstance referenceGeneProduct) throws Exception {
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
        List<Object> ewasModifiedInstanceEdits = getAttributeValues(ewas, ReactomeJavaConstants.modified);
        if (!ewasModifiedInstanceEdits.isEmpty()) {
            // The modified InstanceEdits are held oldest first, so the most recent one is the last, not the first.
            SimpleInstance ewasMostRecentModifiedInstanceEdit =
                (SimpleInstance) ewasModifiedInstanceEdits.get(ewasModifiedInstanceEdits.size() - 1);
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
