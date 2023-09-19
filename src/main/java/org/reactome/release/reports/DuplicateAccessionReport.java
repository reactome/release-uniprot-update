package org.reactome.release.reports;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;

import static org.reactome.release.Utils.isTrEMBLId;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 *         Created 9/19/2023
 */
public class DuplicateAccessionReport implements Reportable {
    private Path outputFilePath;
    private Map<Long, String> duplicateDbIdToUniProtAccession;

    public DuplicateAccessionReport(Path outputFilePath, Map<Long, String> duplicateDbIdToUniProtAccession) {
        this.outputFilePath = outputFilePath;
        this.duplicateDbIdToUniProtAccession = duplicateDbIdToUniProtAccession;
    }

    @Override
    public void writeBody() throws IOException {
        BufferedWriter duplicateDbIDWriter = Files.newBufferedWriter(getFilePath());
        for (long duplicatedDbId : getDuplicateDbIdToUniProtAccession().keySet()) {
            String rgpDuplicateAccession = getDuplicateDbIdToUniProtAccession().get(duplicatedDbId);
            if (isTrEMBLId(rgpDuplicateAccession)) {
                continue;
            }
            duplicateDbIDWriter.write(String.format("%s\t%s%n",rgpDuplicateAccession, duplicatedDbId));
        }
        duplicateDbIDWriter.close();
    }

    @Override
    public Path getFilePath() {
        return this.outputFilePath.resolve("duplicated_db_id.txt");
    }

    @Override
    public String getHeader() {
        final String[] headerColumns = {"ReferenceGeneProduct_Db_Id", "Duplicate_UniProt_Accession"};
        return String.join("\t", headerColumns).concat(System.lineSeparator());
    }

    public Map<Long, String> getDuplicateDbIdToUniProtAccession() {
        return this.duplicateDbIdToUniProtAccession;
    }
}
