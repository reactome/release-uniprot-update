package org.reactome.release.reports;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.List;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 *         Created 9/19/2023
 */
public class TrEMBLAccessionReport implements Reportable {
    private Path outputDirectoryPath;
    private List<String> trEMBLAccessions;

    public TrEMBLAccessionReport(Path outputDirectoryPath, List<String> trEMBLAccessions) {
        this.outputDirectoryPath = outputDirectoryPath;
        this.trEMBLAccessions = trEMBLAccessions;
    }

    @Override
    public void writeBody() throws IOException {
        for (String trEMBLAccession : getTrEMBLAccessions()) {
            Files.write(
                getFilePath(),
                trEMBLAccession.concat(System.lineSeparator()).getBytes(),
                StandardOpenOption.APPEND
            );
        }
    }

    @Override
    public Path getFilePath() {
        return this.outputDirectoryPath.resolve("trembl_to_update.acc");
    }

    @Override
    public String getHeader() {
        return "TrEMBL_Accessions";
    }

    public List<String> getTrEMBLAccessions() {
        return trEMBLAccessions;
    }
}
