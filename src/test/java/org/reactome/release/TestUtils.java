package org.reactome.release;

import org.junit.Test;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.equalTo;
import static org.hamcrest.core.Is.is;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 *         Created 8/29/2023
 */
public class TestUtils {

    @Test
    public void swissProtAccessionGivesFalseForIsTrEMBLID() {
        final String swissProtAccession = "P12345";

        assertThat(
            Utils.isTrEMBLId(swissProtAccession),
            is(equalTo(false))
        );
    }

    @Test
    public void tremblAccessionGivesTrueForIsTrEMBLID() {
        final String trEMBLAccession = "A0A024QZQ1";

        assertThat(
            Utils.isTrEMBLId(trEMBLAccession),
            is(equalTo(true))
        );
    }

    @Test
    public void nonUniProtAccessionGivesFalseForIsTrEMBLID() {
        final String nonUniProtAccession = "A123";

        assertThat(
            Utils.isTrEMBLId(nonUniProtAccession),
            is(equalTo(false))
        );
    }
}
