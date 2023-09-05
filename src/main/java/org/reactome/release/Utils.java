package org.reactome.release;

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
import java.util.List;

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
