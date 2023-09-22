package org.reactome.release.reports;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 *         Created 9/19/2023
 */
public interface Reportable {

    default void writeReport() throws IOException {
        writeHeader();
        writeBody();
    }

    default void writeHeader() throws IOException {
        Files.write(
            getFilePath(),
            getHeader().getBytes(),
            StandardOpenOption.CREATE, StandardOpenOption.APPEND
        );
    }

    void writeBody() throws IOException;

    Path getFilePath();

    String getHeader();
}
