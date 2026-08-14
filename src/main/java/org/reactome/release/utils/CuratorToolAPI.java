package org.reactome.release.utils;

import org.gk.model.InstanceDisplayNameGenerator;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.reactome.curation.CuratorToolWsApplication;
import org.reactome.curation.controller.CurationController;
import org.reactome.curation.model.InstanceList;
import org.reactome.curation.model.NamedReferrerList;
import org.reactome.curation.model.SimpleInstance;
import org.reactome.server.graph.domain.model.DatabaseObject;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.WebApplicationType;
import org.springframework.boot.builder.SpringApplicationBuilder;
import org.springframework.context.ConfigurableApplicationContext;

import java.util.*;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;

import static org.reactome.release.utils.Utils.getFirstAttributeValueAsString;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 * Created 7/5/2026
 */
public class CuratorToolAPI {

    private static final Logger logger = LoggerFactory.getLogger(CuratorToolAPI.class);

    private static final int PAGE_SIZE = 500;

    private static CurationController controller;
    private long personId;

    // The GO instances are loaded once, up-front, and kept in memory for the whole run. The server rejects a
    // commit whose "modified" InstanceEdit does not match the stored one (optimistic locking) and a commit
    // replaces *all* of the stored instance's attributes, so an in-memory copy must be re-read once this run has
    // written to it. These sets record what this run has written to and what it has removed.
    private final Set<Long> committedDbIds = new HashSet<>();
    private final Set<Long> deletedDbIds = new HashSet<>();

    // Source of the negative placeholder dbIds given to instances that are not in the database yet.
    private final AtomicLong placeholderDbIdCounter = new AtomicLong();

    // A species instance is wanted once per SwissProt entry, from the same dozen species names for the whole run, so
    // the instance found (or created) for a name is kept rather than queried again for every entry of that species.
    private final Map<String, SimpleInstance> speciesNameToInstance = new HashMap<>();

    // The instances behind the indexes the run looks reference gene products and isoforms up in; see
    // fetchInflatedRGPInstances.
    private List<SimpleInstance> inflatedRGPInstances;

    private ConfigurableApplicationContext applicationContext;

    public CuratorToolAPI(long personId) {
        if (controller == null) {
            controller = this.initController();
            if (controller == null) {
                throw new IllegalStateException("Failed to initialize CuratorToolAPI: controller is null");
            }
        }
        this.personId = personId;
    }

    // The following code is copied directly from the slicing tool project.
    private CurationController initController() {
        try {
            // curator-tool-ws's bundled application.properties forces DEBUG for these loggers and binds the
            // HTTP connector to 9090. System properties outrank a classpath application.properties in Spring
            // Boot's precedence order, so these settings take hold for the batch run without editing
            // curator-tool-ws. (SpringApplicationBuilder.properties(...) are default/lowest precedence and
            // would NOT override application.properties.)
            System.setProperty("logging.level.org.springframework.data.neo4j", "WARN");
            System.setProperty("logging.level.org.springframework.security", "WARN");
            // Disable the HTTP server; the full servlet context is kept for correct AspectJ wiring.
            System.setProperty("server.port", "-1");

            applicationContext = new SpringApplicationBuilder(CuratorToolWsApplication.class)
                .web(WebApplicationType.SERVLET)
                .run();
            return applicationContext.getBean(CurationController.class);
        }
        catch (Exception e) {
            logger.error("GraphDBInstanceManager.initController(): " + e.getMessage(), e);
        }
        return null;
    }

    public SimpleInstance commit(SimpleInstance simpleInstance) {
        if (simpleInstance.getDefaultPersonId() == null) {
            simpleInstance.setDefaultPersonId(getPersonId());
        }

        boolean isNewInstance = simpleInstance.getDbId() == null || simpleInstance.getDbId() < 0;
        if (simpleInstance.getDbId() == null) {
            simpleInstance.setDbId(nextPlaceholderDbId());
        }

        SimpleInstance committedInstance = controller.commit(simpleInstance);
        // Nothing was stored, so callers cannot go on to use the instance -- for a new one, the dbId it would have
        // been given is the very thing that is missing. Reported here rather than as a NullPointerException at
        // whichever caller dereferences the result first.
        if (committedInstance == null) {
            throw new IllegalStateException("Commit of " + simpleInstance + " returned no instance");
        }

        if (isNewInstance) {
            // The response carries the dbId the database assigned in place of the placeholder, and it is the only
            // place it is reported, so it is copied back onto the instance the caller holds.
            simpleInstance.setDbId(committedInstance.getDbId());
        }

        Long dbId = simpleInstance.getDbId();
        if (dbId != null && dbId > 0) {
            committedDbIds.add(dbId);
            deletedDbIds.remove(dbId);
        }
        return committedInstance;
    }

    /**
     * Returns a placeholder dbId for an instance that is not in the database yet. curator-tool-ws identifies a new
     * instance by a NEGATIVE dbId: it is what makes a commit store the instance (recording the author in its
     * "created" slot rather than in "modified") and replace the placeholder with a real dbId. A null dbId is not a
     * substitute -- curator-tool-ws unboxes the dbId without a null check while working out what to store
     * (DatabaseObjectInstanceConverter.convert and CurationService.grepNewInstances), so committing an instance
     * with a null dbId fails with a NullPointerException.
     *
     * @return a negative dbId, unused by any other instance created during this run.
     */
    private long nextPlaceholderDbId() {
        return -placeholderDbIdCounter.incrementAndGet();
    }

    /**
     * Returns a copy of the instance that is up to date with the database if this run has already committed it,
     * and the instance itself otherwise. Returns null if this run has deleted the instance, in which case the
     * caller must stop using it -- committing it would re-create the deleted instance.
     *
     * An instance that this run has committed is stale in two ways: its "modified" InstanceEdit no longer matches
     * the stored one, so a further commit is rejected with an InstanceChangedException, and any attribute written
     * by that commit (or by a commit of a fresh copy of the same instance elsewhere in the run) still holds its
     * pre-commit value, which the next commit would write back over the stored value.
     *
     * @param instance - the instance to refresh.
     * @return an up-to-date instance, or null if the instance has been deleted by this run.
     */
    public SimpleInstance refresh(SimpleInstance instance) {
        Long dbId = instance.getDbId();
        if (dbId == null) {
            return instance; // Never committed, so there is nothing stored to be out of date with.
        }
        if (deletedDbIds.contains(dbId)) {
            return null;
        }
        return committedDbIds.contains(dbId) ? inflate(instance) : instance;
    }

    public SimpleInstance findByDbId(long dbId) {
        DatabaseObject databaseObject = controller.findByDdId(dbId);
        if (databaseObject == null) {
            return null;
        }

        try {
            return controller.getConverter().convert(databaseObject);
        } catch (Exception e) {
            throw new RuntimeException("Unable to convert DatabaseObject " + databaseObject + " to SimpleInstance", e);
        }
    }

    public SimpleInstance findByDisplayName(String className, String displayName) {
        // NB: the controller takes the display name first and a comma-separated list of class names second.
        return controller.findByDisplayName(displayName, className);
    }

    public SimpleInstance fetchUniProtReferenceDatabase() {
        SimpleInstance uniProtReferenceDatabase = findByDisplayName(ReactomeJavaConstants.ReferenceDatabase, "UniProt");
        if (uniProtReferenceDatabase == null) {
            throw new IllegalStateException("No " + ReactomeJavaConstants.ReferenceDatabase +
                " instance with the display name 'UniProt' exists in the database");
        }
        return uniProtReferenceDatabase;
    }

    public Map<String, Long> getRGPAccessionToDbIdMap() {
        Map<String, Long> identifierToDbId = new HashMap<>();
        for (SimpleInstance rgp : fetchUniProtRGPInstances()) {
            String identifier = getIdentifierFromDisplayName(rgp);
            if (identifier != null && !identifier.isEmpty()) {
                identifierToDbId.put(identifier, rgp.getDbId());
            }
        }
        return identifierToDbId;
    }

    /**
     * Returns every reference gene product of the database, inflated, indexed by its identifier. The instances
     * themselves are kept rather than just their dbIds because inflating them is the cost of building this index in
     * the first place: with them in hand, the run has each instance's attributes without searching for the same
     * accession again, entry by entry, as it goes through the SwissProt file.
     *
     * An identifier maps to a list because the database can hold more than one instance for it: the duplicate master
     * sequences the run reports on, and the isoforms of a master sequence, which carry their parent's accession as
     * their own identifier.
     *
     * No reference database filter is applied, matching getReferenceGeneProductsByIdentifier, which this index
     * stands in for.
     *
     * An indexed instance goes stale once this run commits to it, so read it back through refresh.
     *
     * @return the reference gene products of the database, indexed by identifier.
     */
    public Map<String, List<SimpleInstance>> getRGPIdentifierToInstancesMap() {
        Map<String, List<SimpleInstance>> rgpIdentifierToInstances = new HashMap<>();
        for (SimpleInstance referenceGeneProduct : fetchInflatedRGPInstances()) {
            String identifier = (String) referenceGeneProduct.getAttribute(ReactomeJavaConstants.identifier);
            if (identifier != null && !identifier.isEmpty()) {
                rgpIdentifierToInstances
                    .computeIfAbsent(identifier, k -> new ArrayList<>())
                    .add(referenceGeneProduct);
            }
        }
        return rgpIdentifierToInstances;
    }

    /**
     * Returns every ReferenceIsoform of the database, inflated, indexed by its variant identifier -- kept rather
     * than discarded for its dbId for the reason given on getRGPIdentifierToInstancesMap, and drawn from the same
     * instances, since a ReferenceIsoform inherits from ReferenceGeneProduct and so is one of them.
     *
     * A variant identifier maps to a list because the database can hold more than one isoform for it -- the
     * duplicates the run reports on.
     *
     * @return the isoforms of the database, indexed by variant identifier.
     */
    public Map<String, List<SimpleInstance>> getIsoformAccessionToInstancesMap() {
        Map<String, List<SimpleInstance>> isoformIdentifierToInstances = new HashMap<>();
        for (SimpleInstance referenceGeneProduct : fetchInflatedRGPInstances()) {
            if (!ReactomeJavaConstants.ReferenceIsoform.equals(referenceGeneProduct.getSchemaClassName())) {
                continue;
            }

            String variantIdentifier =
                (String) referenceGeneProduct.getAttribute(ReactomeJavaConstants.variantIdentifier);
            if (variantIdentifier != null && !variantIdentifier.isEmpty()) {
                isoformIdentifierToInstances
                    .computeIfAbsent(variantIdentifier, k -> new ArrayList<>())
                    .add(referenceGeneProduct);
            }
        }
        return isoformIdentifierToInstances;
    }

    /**
     * Returns the variant identifiers of the UniProt isoforms of the database -- the isoforms the run accounts for
     * against the SwissProt file, and deletes or reports on where the file no longer carries them. The index itself
     * is wider than this, holding every isoform whatever its reference database, because it stands in for a search
     * that applied no such filter.
     *
     * Which isoforms are the UniProt ones is settled by the dbIds a UniProt-filtered search returns, rather than by
     * reading a reference database off the instances and deciding here what counts as UniProt. That search returns
     * shells, which is all this needs, so it inflates nothing.
     *
     * @param isoformAccessionToInstances - the isoforms of the database, indexed by variant identifier.
     * @return the variant identifiers of the UniProt isoforms.
     */
    public Set<String> getUniProtIsoformAccessions(Map<String, List<SimpleInstance>> isoformAccessionToInstances) {
        Set<Long> uniProtIsoformDbIds =
            fetchInstancesForClass(ReactomeJavaConstants.ReferenceIsoform, "UniProt")
                .stream()
                .map(SimpleInstance::getDbId)
                .collect(Collectors.toSet());

        Set<String> uniProtIsoformAccessions = new HashSet<>();
        for (Map.Entry<String, List<SimpleInstance>> indexedIsoforms : isoformAccessionToInstances.entrySet()) {
            boolean anyIsAUniProtIsoform = indexedIsoforms.getValue()
                .stream()
                .anyMatch(isoform -> uniProtIsoformDbIds.contains(isoform.getDbId()));

            if (anyIsAUniProtIsoform) {
                uniProtIsoformAccessions.add(indexedIsoforms.getKey());
            }
        }
        return uniProtIsoformAccessions;
    }

    /**
     * Returns every ReferenceDNASequence in the database, inflated, indexed by its identifier -- kept rather than
     * discarded for its dbId for the reason given on getRGPIdentifierToInstancesMap.
     *
     * @return the reference DNA sequences of the database, indexed by identifier.
     */
    public Map<String, SimpleInstance> getRDSIdentifierToInstanceMap() {
        Map<String, SimpleInstance> rdsIdentifierToInstance = new HashMap<>();
        for (SimpleInstance referenceDNASequence : fetchRDSInstances()) {
            String rdsIdentifier = (String) referenceDNASequence.getAttribute(ReactomeJavaConstants.identifier);
            if (rdsIdentifier != null && !rdsIdentifier.isEmpty()) {
                rdsIdentifierToInstance.put(rdsIdentifier, referenceDNASequence);
            }
        }
        return rdsIdentifierToInstance;
    }

    public SimpleInstance getSpeciesInstance(String speciesName) throws Exception {
        SimpleInstance cachedSpeciesInstance = speciesNameToInstance.get(speciesName);
        if (cachedSpeciesInstance != null) {
            return cachedSpeciesInstance;
        }

        SimpleInstance speciesInstance = fetchSpecies(speciesName);
        if (speciesInstance == null) {
            speciesInstance = createNewSpeciesInstance(speciesName);
            commit(speciesInstance);
        }
        speciesNameToInstance.put(speciesName, speciesInstance);
        return speciesInstance;
    }

    public SimpleInstance getHumanEnsEMBLGeneReferenceDatabase() {
        InstanceList ensEMBLHumanReferenceDatabaseInstances = controller.searchInstances(
            ReactomeJavaConstants.ReferenceDatabase,
            0,
            1,
            Optional.of(ReactomeJavaConstants.name),
            Optional.of("equal"),
            Optional.of("ENSEMBL")
        );

        if (ensEMBLHumanReferenceDatabaseInstances == null || ensEMBLHumanReferenceDatabaseInstances.isEmpty()) {
            throw new RuntimeException("Could not get EnsEMBL human gene reference database");
        }

        return ensEMBLHumanReferenceDatabaseInstances.getInstances().get(0);
    }

    public List<SimpleInstance> getReferenceGeneProductsByIdentifier(String identifier) {
        return getInstancesByAttribute(
            ReactomeJavaConstants.ReferenceGeneProduct, ReactomeJavaConstants.identifier, identifier);
    }

    public List<SimpleInstance> getReferenceIsoformByVariantIdentifier(String variantIdentifier) {
        return getInstancesByAttribute(
            ReactomeJavaConstants.ReferenceIsoform, ReactomeJavaConstants.variantIdentifier, variantIdentifier
        );
    }

    public void deleteInstance(SimpleInstance instance) {
        if (instance.getDefaultPersonId() == null) {
            instance.setDefaultPersonId(getPersonId());
        }

        controller.delete(instance);
        if (instance.getDbId() != null) {
            deletedDbIds.add(instance.getDbId());
            committedDbIds.remove(instance.getDbId());
        }
    }

    public void deleteByDbId(long noReferrerDbId) {
        SimpleInstance instance = controller.findByDdIdInInstance(noReferrerDbId);
        // Nothing stored under that db id -- there is nothing to delete, and passing the null on would fail inside
        // curator-tool-ws instead of saying which db id was asked for.
        if (instance == null) {
            logger.warn("No instance found for db id " + noReferrerDbId + " -- nothing to delete");
            return;
        }
        deleteInstance(instance);
    }

    public void close() {
        // Null when this instance did not create the context: the controller is static and only the first instance
        // built one.
        if (applicationContext != null) {
            applicationContext.close();
        }
    }

    public SimpleInstance inflate(SimpleInstance shellInstance) {
        return controller.findByDdIdInInstance(shellInstance.getDbId());
    }

    public List<SimpleInstance> getReferrers(SimpleInstance instance, String referrerAttributeName) throws Exception {
        return controller.getReferrers(instance.getDbId())
            .stream()
            .filter(g -> referrerAttributeName.equals(g.getAttributeName()))
            .findFirst()
            .map(NamedReferrerList::getReferrers)
            .orElse(Collections.emptyList());
    }

    public Collection<NamedReferrerList> getReferrers(SimpleInstance instance) throws Exception {
        return controller.getReferrers(instance.getDbId());
    }

    public long getPersonId() {
        return this.personId;
    }

    private List<SimpleInstance> fetchUniProtRGPInstances() {
        return fetchInstancesForClass(ReactomeJavaConstants.ReferenceGeneProduct, "UniProt");
    }

    /**
     * Returns every instance carrying the ReferenceGeneProduct label -- the master sequences and, since a
     * ReferenceIsoform inherits from ReferenceGeneProduct, the isoforms too -- inflated. Held once the first index is
     * built from it, so that the indexes built afterwards share the one set of instances rather than inflating the
     * same instance again for each of them.
     *
     * @return the reference gene products of the database, inflated.
     */
    private List<SimpleInstance> fetchInflatedRGPInstances() {
        if (inflatedRGPInstances == null) {
            inflatedRGPInstances = fetchInstancesForClass(ReactomeJavaConstants.ReferenceGeneProduct)
                .parallelStream()
                .map(this::inflate)
                .filter(Objects::nonNull)
                .collect(Collectors.toList());
        }
        return inflatedRGPInstances;
    }

    public List<SimpleInstance> fetchUniProtReferenceIsoformInstances() {
        return fetchInstancesForClass(ReactomeJavaConstants.ReferenceIsoform, "UniProt")
            .parallelStream()
            .map(this::inflate)
            .collect(Collectors.toList());
    }

    public void updateReferenceGeneProductDisplayNames() {
        for (SimpleInstance rgpInstance : fetchUniProtRGPInstances()) {
            updateDisplayName(rgpInstance);
        }
    }

    /**
     * A ReferenceIsoform node carries the label of every class it inherits from, so a search for
     * ReferenceGeneProduct returns the isoforms too and updateReferenceGeneProductDisplayNames covers them -- which is
     * why Main filters isoforms out of such a search's results where it wants the master sequences alone. This pass
     * exists so that the isoforms are covered even if that ceases to hold; where it does hold, the display names it
     * computes are the ones already stored and nothing is committed.
     */
    public void updateReferenceIsoformDisplayNames() {
        // The shell instances of the search are enough here, since updateDisplayName inflates what it is given --
        // fetchUniProtReferenceIsoformInstances would inflate every isoform a second time.
        for (SimpleInstance referenceIsoformInstance :
             fetchInstancesForClass(ReactomeJavaConstants.ReferenceIsoform, "UniProt")) {

            updateDisplayName(referenceIsoformInstance);
        }
    }

    private void updateDisplayName(SimpleInstance referenceSequence) {
        // A search returns shell instances -- dbId, displayName and schemaClass only -- so the instance has to be
        // inflated before there are any attributes to build a display name out of. Committing a shell would also
        // clear every attribute of the stored instance, since a commit replaces all of them.
        SimpleInstance inflatedReferenceSequence = inflate(referenceSequence);
        if (inflatedReferenceSequence == null) {
            logger.warn("No instance found for db id " + referenceSequence.getDbId() +
                " -- skipping its display name update");
            return;
        }

        String currentDisplayName = inflatedReferenceSequence.getDisplayName();
        String newDisplayName = getReferenceSequenceDisplayName(inflatedReferenceSequence);

        if (!newDisplayName.equals(currentDisplayName)) {
            inflatedReferenceSequence.setDisplayName(newDisplayName);
            commit(inflatedReferenceSequence);
        }
    }

    public String getReferenceSequenceDisplayName(SimpleInstance referenceSequence) {
        String dbName = null;
        SimpleInstance refDB =
            (SimpleInstance) referenceSequence.getAttribute(ReactomeJavaConstants.referenceDatabase);
        if (refDB != null) {
            dbName = refDB.getDisplayName();
        }
        if (dbName == null) {
            dbName = "Unknown";
        }

        String identifier = getFirstAttributeValueAsString(referenceSequence, ReactomeJavaConstants.variantIdentifier);

        if (identifier == null) {
            identifier = getFirstAttributeValueAsString(referenceSequence, ReactomeJavaConstants.identifier);
        }
        if (identifier == null) {
            identifier = "Unknown";
        }

        // geneName and name are multi-valued attributes, so their values are Lists rather than the Strings this used
        // to cast them to; the first value is the one that belongs in the display name.
        String name = getFirstAttributeValueAsString(referenceSequence, ReactomeJavaConstants.geneName);
        if (name == null) {
            name = getFirstAttributeValueAsString(referenceSequence, ReactomeJavaConstants.name);
        }
        if (name == null) {
            name = "Unknown";
        }
        return dbName + ":" + identifier + " " + name;
    }

    private List<SimpleInstance> fetchRDSInstances() {
        return fetchInstancesForClass(ReactomeJavaConstants.ReferenceDNASequence)
            .parallelStream()
            .map(this::inflate)
            .collect(Collectors.toList());
    }

    private List<SimpleInstance> fetchInstancesForClass(String className) {
        return fetchInstancesForClass(className, null);
    }

    private List<SimpleInstance> fetchInstancesForClass(String className, String referenceDatabaseName) {
        // curator-tool-ws only skips the reference database filter when the query parameters are ABSENT.
        // An empty search key is still a present Optional and matches a reference database whose display
        // name is "" -- i.e. nothing at all.
        boolean filterByReferenceDatabase = referenceDatabaseName != null && !referenceDatabaseName.isEmpty();
        Optional<String> attribute = filterByReferenceDatabase
            ? Optional.of(ReactomeJavaConstants.referenceDatabase)
            : Optional.empty();
        Optional<String> operand = filterByReferenceDatabase ? Optional.of("equal") : Optional.empty();
        Optional<String> searchKey = filterByReferenceDatabase ? Optional.of(referenceDatabaseName) : Optional.empty();

        List<SimpleInstance> instances = new ArrayList<>();

        int skip = 0;
        Integer total = null;
        do {
            InstanceList page = controller.searchInstances(className, skip, PAGE_SIZE, attribute, operand, searchKey);

            if (total == null) {
                total = page.getTotalCount() != null ? page.getTotalCount() : 0;
            }
            instances.addAll(page.getInstances());
            skip += PAGE_SIZE;
        } while (skip < total);

        return instances;
    }

    private List<SimpleInstance> getInstancesByAttribute(String className, String attributeName, String attributeValue) {
        List<SimpleInstance> instances = new ArrayList<>();

        int skip = 0;
        Integer total = null;
        do {
            InstanceList page = controller.searchInstances(
                className,
                skip,
                PAGE_SIZE,
                Optional.of(attributeName),
                Optional.of("equal"),
                Optional.of(attributeValue)
            );

            if (total == null) {
                total = page.getTotalCount() != null ? page.getTotalCount() : 0;
            }
            instances.addAll(page.getInstances());
            skip += PAGE_SIZE;
        } while (skip < total);

        return instances.parallelStream().map(this::inflate).collect(Collectors.toList());
    }

    private String getIdentifierFromDisplayName(SimpleInstance referenceSequence) {
        String displayName = referenceSequence.getDisplayName();
        int colonIndex = displayName.indexOf(':');
        int spaceIndex = displayName.indexOf(' ', colonIndex + 1);
        return displayName.substring(colonIndex + 1, spaceIndex < 0 ? displayName.length() : spaceIndex);
    }

    public SimpleInstance fetchSpecies(String speciesName) {
        InstanceList speciesInstances = controller.searchInstances(
            ReactomeJavaConstants.Species,
            0,
            1,
            Optional.of(ReactomeJavaConstants.name),
            Optional.of("equal"),
            Optional.of(speciesName)
        );

        return !speciesInstances.isEmpty() ? speciesInstances.getInstances().get(0) : null;
    }

    private SimpleInstance createNewSpeciesInstance(String speciesName) throws Exception {
        SimpleInstance speciesInstance = new SimpleInstance();
        speciesInstance.setSchemaClassName(ReactomeJavaConstants.Species);
        speciesInstance.setAttribute(ReactomeJavaConstants.name, Collections.singletonList(speciesName));
        speciesInstance.setDisplayName(speciesName);
        return speciesInstance;
    }

}
