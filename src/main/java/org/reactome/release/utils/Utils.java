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
import java.util.List;

import org.reactome.curation.model.SimpleInstance;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 *         Created 8/28/2023
 */
public class Utils {
    public static String getUpdateDirectory() throws URISyntaxException {
        return Paths.get(Utils.class.getClassLoader().getResource(".").toURI()).toString();
    }

    public static boolean isTrEMBLId(String potentialTrEMBLId) {
        final String uniProtAccessionURLAsString = "https://rest.uniprot.org/uniprotkb/" + potentialTrEMBLId + ".txt";

        HttpURLConnection uniProtAccessionHttpURLConnection = null;
        try {
            URL uniProtAccessionURL = new URL(uniProtAccessionURLAsString);
            uniProtAccessionHttpURLConnection = (HttpURLConnection) uniProtAccessionURL.openConnection();
            uniProtAccessionHttpURLConnection.setRequestMethod("GET");
            BufferedReader uniProtAccessionReader = new BufferedReader(
                new InputStreamReader(uniProtAccessionHttpURLConnection.getInputStream()));

            return uniProtAccessionReader.lines().anyMatch(
                line -> line.matches("^.*(Unreviewed|TrEMBL).*$")
            );
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
