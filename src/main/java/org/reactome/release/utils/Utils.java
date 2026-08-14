package org.reactome.release.utils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.reactome.curation.model.SimpleInstance;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 *         Created 8/28/2023
 */
public class Utils {
    private static final String UNIPROT_ACCESSIONS_URL = "https://rest.uniprot.org/uniprotkb/accessions";

    // The most accessions UniProt takes in one request to its accessions end point.
    private static final int MAX_ACCESSIONS_PER_REQUEST = 100;

    // The accession formats UniProt assigns; taken from the pattern UniProt documents for them. Anything else is
    // held by no UniProt entry, reviewed or not, so it is answered for without a request being made at all -- which
    // also keeps a malformed accession from failing the request of the batch it would otherwise be sent in.
    private static final Pattern UNIPROT_ACCESSION_PATTERN = Pattern.compile(
        "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}");

    public static String getUpdateDirectory() throws URISyntaxException {
        return Paths.get(Utils.class.getClassLoader().getResource(".").toURI()).toString();
    }

    /**
     * Returns which of the accessions UniProt holds as TrEMBL (i.e. unreviewed) entries.
     *
     * The accessions are asked about a hundred at a time rather than one at a time, and only their reviewed status
     * is asked for rather than their whole entries, because this is asked for every accession of the database the
     * SwissProt file no longer carries and each request is a round trip to UniProt.
     *
     * @param potentialTrEMBLIds - the accessions to look up.
     * @return the accessions that are TrEMBL accessions.
     */
    public static Set<String> getTrEMBLIds(Collection<String> potentialTrEMBLIds) {
        List<String> accessions = potentialTrEMBLIds.stream().distinct().collect(Collectors.toList());

        // Only accessions of a format UniProt assigns go in a batch: its accessions end point rejects a request
        // carrying one of any other format, and rejects it whole, so a single malformed accession would cost the
        // answer for the hundred it was sent with. The rest fall through to being looked up on their own below.
        List<String> batchableAccessions =
            accessions.stream().filter(Utils::isAUniProtAccession).collect(Collectors.toList());

        Map<String, Boolean> accessionToIsTrEMBL = new HashMap<>();
        for (int batchStart = 0; batchStart < batchableAccessions.size(); batchStart += MAX_ACCESSIONS_PER_REQUEST) {
            List<String> batch = batchableAccessions.subList(
                batchStart, Math.min(batchStart + MAX_ACCESSIONS_PER_REQUEST, batchableAccessions.size()));

            accessionToIsTrEMBL.putAll(fetchAccessionToIsTrEMBL(batch));
        }

        Set<String> tremblIds = new HashSet<>();
        for (String accession : accessions) {
            // An accession no batch answered for is looked up on its own, so that the answer for it is the one a
            // look-up of it on its own has always given: it was held back from the batches for its format, or the
            // batch carrying it failed, or it is now a secondary accession of some entry and was answered for under
            // that entry's primary accession rather than under itself.
            boolean isTrEMBL = accessionToIsTrEMBL.containsKey(accession)
                ? accessionToIsTrEMBL.get(accession)
                : isTrEMBLId(accession);

            if (isTrEMBL) {
                tremblIds.add(accession);
            }
        }
        return tremblIds;
    }

    public static boolean isTrEMBLId(String potentialTrEMBLId) {
        final String uniProtAccessionURLAsString = "https://rest.uniprot.org/uniprotkb/" + potentialTrEMBLId + ".txt";

        HttpURLConnection uniProtAccessionHttpURLConnection = null;
        // Read to the end and closed, rather than stopped at the first line of interest and left open: a response
        // left undrained keeps its connection out of the keep-alive pool, so each accession looked up after it pays
        // to open a new one.
        try {
            URL uniProtAccessionURL = new URL(uniProtAccessionURLAsString);
            uniProtAccessionHttpURLConnection = (HttpURLConnection) uniProtAccessionURL.openConnection();
            uniProtAccessionHttpURLConnection.setRequestMethod("GET");

            try (BufferedReader uniProtAccessionReader = new BufferedReader(
                new InputStreamReader(uniProtAccessionHttpURLConnection.getInputStream()))) {

                List<String> entryLines = uniProtAccessionReader.lines().collect(Collectors.toList());
                return entryLines.stream().anyMatch(line -> line.matches("^.*(Unreviewed|TrEMBL).*$"));
            }
        } catch (IOException e) {
            if (serverUnavailable(uniProtAccessionHttpURLConnection)) {
                throw new RuntimeException("Unable to connect to UniProt RESTful server ", e);
            } else {
                System.err.println("Unable to get content from " + uniProtAccessionURLAsString + ": " + e);
                return false;
            }
        }
    }

    /**
     * Asks UniProt for the reviewed status of a batch of accessions, as the accession and whether it is a TrEMBL one.
     *
     * Only the accessions UniProt answers for are in the map returned: an accession it holds no entry under is left
     * out rather than reported as not being a TrEMBL accession, so that the caller can tell the two apart.
     *
     * @param accessions - the accessions of one request's worth.
     * @return each answered accession, and whether it is a TrEMBL accession.
     */
    private static Map<String, Boolean> fetchAccessionToIsTrEMBL(List<String> accessions) {
        final String uniProtAccessionsURLAsString = UNIPROT_ACCESSIONS_URL +
            "?accessions=" + String.join(",", accessions) + "&fields=accession,reviewed&format=tsv";

        HttpURLConnection uniProtAccessionsHttpURLConnection = null;
        try {
            URL uniProtAccessionsURL = new URL(uniProtAccessionsURLAsString);
            uniProtAccessionsHttpURLConnection = (HttpURLConnection) uniProtAccessionsURL.openConnection();
            uniProtAccessionsHttpURLConnection.setRequestMethod("GET");

            try (BufferedReader uniProtAccessionsReader = new BufferedReader(
                new InputStreamReader(uniProtAccessionsHttpURLConnection.getInputStream()))) {

                return parseAccessionToIsTrEMBL(uniProtAccessionsReader.lines().collect(Collectors.toList()));
            }
        } catch (IOException e) {
            if (serverUnavailable(uniProtAccessionsHttpURLConnection)) {
                throw new RuntimeException("Unable to connect to UniProt RESTful server ", e);
            }

            // Left empty rather than answered for here, so that each accession of the batch falls back to being
            // looked up on its own and one failed request does not decide for a hundred accessions at once.
            System.err.println("Unable to get content from " + uniProtAccessionsURLAsString + ": " + e);
            return Collections.emptyMap();
        }
    }

    /**
     * Parses the tab separated accession and reviewed status of each row of a UniProt response. The header row, and
     * any row without both columns, is passed over.
     *
     * @param responseLines - the lines of the response.
     * @return each accession of the response, and whether it is a TrEMBL accession.
     */
    private static Map<String, Boolean> parseAccessionToIsTrEMBL(List<String> responseLines) {
        final String unreviewedStatus = "unreviewed";
        final String reviewedStatus = "reviewed";

        Map<String, Boolean> accessionToIsTrEMBL = new HashMap<>();
        for (String responseLine : responseLines) {
            String[] columns = responseLine.split("\t");
            if (columns.length < 2) {
                continue;
            }

            String accession = columns[0];
            String status = columns[1].toLowerCase();
            // Anything that is neither status is not an entry's row -- the header row above all -- and saying so is
            // left to the accession being looked up on its own rather than guessed at from an unrecognised row.
            if (status.equals(unreviewedStatus) || status.equals(reviewedStatus)) {
                accessionToIsTrEMBL.put(accession, status.equals(unreviewedStatus));
            }
        }
        return accessionToIsTrEMBL;
    }

    private static boolean isAUniProtAccession(String potentialUniProtAccession) {
        return potentialUniProtAccession != null &&
            UNIPROT_ACCESSION_PATTERN.matcher(potentialUniProtAccession).matches();
    }

    /**
     * SimpleInstance holds a single-valued attribute as the value itself, a multi-valued attribute as a List, and
     * returns null for an attribute with no value, so every read is normalized to a list here.
     *
     * @param instance - the instance to read the attribute from.
     * @param attributeName - the name of the attribute to read.
     * @return the attribute's values, empty if it has none.
     */
    public static List<Object> getAttributeValues(SimpleInstance instance, String attributeName) {
        Object value = instance.getAttribute(attributeName);
        if (value == null) {
            return Collections.emptyList();
        }
        return value instanceof List ? new ArrayList<>((List<?>) value) : Collections.singletonList(value);
    }

    /**
     * Returns the attribute's first value as a String, for a single-valued String attribute or where only the first
     * value of a multi-valued attribute is wanted (e.g. the gene name used in a display name).
     *
     * @param instance - the instance to read the attribute from.
     * @param attributeName - the name of the attribute to read.
     * @return the attribute's first value as a String, or null if it has no values.
     */
    public static String getFirstAttributeValueAsString(SimpleInstance instance, String attributeName) {
        List<Object> values = getAttributeValues(instance, attributeName);
        return values.isEmpty() ? null : values.get(0).toString();
    }

    public static <E> List<E> emptyListIfNull(List<E> list) {
        return list != null ? list : new ArrayList<>();
    }

    public static <E> List<E> emptyListIfNull(Collection<E> collection) {
        return collection != null ? new ArrayList<>(collection) : new ArrayList<>();
    }

    public static void writeAndPrint(String line) throws IOException {
        Files.write(
            Paths.get("test.txt"),
            line.getBytes(),
            StandardOpenOption.CREATE, StandardOpenOption.APPEND
        );
        System.out.println(line);
    }

    private static boolean serverUnavailable(HttpURLConnection uniProtAccessionHttpURLConnection) {
        try {
            return uniProtAccessionHttpURLConnection.getResponseCode() >= 500;
        } catch (IOException e) {
            return true;
        }
    }
}
