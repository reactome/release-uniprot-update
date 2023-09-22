package org.reactome.release;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.zip.GZIPInputStream;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 *         Created 9/19/2023
 */
public class SwissProtFileProcessor {
    private Path updateDirectoryPath;
    private Path swissProtFilePath;

    public SwissProtFileProcessor(Path updateDirectoryPath) {
        this.updateDirectoryPath = updateDirectoryPath;
    }

    public Path getSwissProtFilePath() throws IOException {
        if (this.swissProtFilePath == null) {
            this.swissProtFilePath = Files.list(getUpdateDirectoryPath())
                .filter(path -> path.toString().contains("uniprot_sprot.xml"))
                .findFirst()
                .orElseThrow(() -> new RuntimeException("Can't find SwissProt file uniprot_sprot.xml[.gz]"));
        }
        return this.swissProtFilePath;
    }

    public BufferedReader getFileReader() throws IOException {
        gunzipSwissProtFileIfZipped();
        return Files.newBufferedReader(getSwissProtFilePath());
    }

    /**
     * Gunzips SwissProt XML file is gzipped.
     * @return <code>true</code> if file is gunzipped;<code>false</code> otherwise if unchanged
     */
    private boolean gunzipSwissProtFileIfZipped() throws IOException {
        if (getSwissProtFilePath().endsWith(".gz")) {
            System.out.println("Found SwissProt file with .gz extension - unzipping");
            gunzipOrThrow(getSwissProtFilePath());
            this.swissProtFilePath = Paths.get(getSwissProtFilePath().toString().replace(".gz",""));
            return true;
        }
        return false;
    }

    private void gunzipOrThrow(Path filePath) throws IOException {
        try (GZIPInputStream gzipInputStream = new GZIPInputStream(new FileInputStream(filePath.toFile()));
             FileOutputStream unzipOutputStream = new FileOutputStream(
                 filePath.toString().replace(".gz",""))) {

            byte[] buffer = new byte[1024];
            int len;
            while ((len = gzipInputStream.read(buffer)) > 0) {
                unzipOutputStream.write(buffer, 0, len);
            }
        }
    }

    private Path getUpdateDirectoryPath() {
        return this.updateDirectoryPath;
    }
}
