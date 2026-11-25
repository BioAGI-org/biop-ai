# Copyright 2008-2016 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Testing online code for fetching NCBI qblast.

Uses Bio.Blast.qblast() to run some online blast queries, get XML blast results
back, and then checks Bio.Blast.read() can read them. One smoke test does the
same with the deprecated Bio.Blast.NCBIWWW and Bio.Blast.NCBIXML modules, always
against a mock file, so it does not repeat an online search.

Goals:
    - Make sure that all retrieval is working as expected.
    - Make sure we can parse the latest XML format being used by the NCBI.

If an internet connection is available and run_tests.py has not given the argument
'--offline', these tests will run online (and can take a long time). Otherwise, the
module runs offline using a mocked version of urllib.request.urlopen, which returns
file contents instead of internet responses.

IMPORTANT:
If you add new tests (or change existing tests) you must provide new 'mock' xml
files which fulfill the requirements of the respective tests. Add them to
MOCK_RESULTS under the name of the test that uses them.

"""

import unittest
import warnings
from io import BufferedReader
from io import BytesIO
from unittest import mock
from urllib.error import HTTPError

import requires_internet
import support

from Bio import Blast
from Bio import BiopythonDeprecationWarning
from Bio import MissingExternalDependencyError

with warnings.catch_warnings():
    # Still tested by test_deprecated_ncbiwww_ncbixml while they ship.
    warnings.simplefilter("ignore", BiopythonDeprecationWarning)
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML

Blast.email = "biopython@biopython.org"
NCBIWWW.email = "biopython@biopython.org"


try:
    requires_internet.check()
    # this will set 'require_internet.check.available' (if not already done before)
    # it is False if we don't have internet or if we run the tests with --offline
except MissingExternalDependencyError:
    # need to catch this error, because we want to run the tests even offline
    pass  # 'requires_internet_check.available' is now False

# When offline, each qblast call gets the 'wait' page and then the next of
# these result files, mimicking an NCBI qblast response.
MOCK_RESULTS = {
    "test_blastp_nr_actin": ["mock_actin.xml"],
    "test_deprecated_ncbiwww_ncbixml": ["mock_actin.xml"],
    "test_discomegablast": ["mock_disco.xml"],
    "test_orchid_est": ["mock_orchid.xml"],
    "test_pcr_primers": ["mock_pcr.xml"],
    "test_short_query": ["mock_short_empty.xml", "mock_short_result.xml"],
}

# Tests that use their mock files even when online.
ALWAYS_MOCKED = {"test_deprecated_ncbiwww_ncbixml"}


#####################################################################

# List of qblast requests stored as a tuple of parameters:
# - program
# - database
# - query identifier or sequence
# - expectation value threshold
# - Entrez filter string (or None)
# - list of hit identifiers expected to be found (or None if expect 0)


class TestQblast(unittest.TestCase):
    def setUp(self):
        if (
            not requires_internet.check.available
            or self._testMethodName in ALWAYS_MOCKED
        ):
            # if the tests are running with --offline or if no internet is
            # available, we mock urlopen, so that it returns files from the
            # Blast folder. Bio.Blast.qblast peeks at the stream, so each
            # response is a BufferedReader.
            blast = support.DATA / "Blast"
            wait = (blast / "mock_wait.html").read_bytes()
            responses = iter(
                [
                    BufferedReader(BytesIO(data))
                    for name in MOCK_RESULTS.get(self._testMethodName, [])
                    for data in (wait, (blast / name).read_bytes())
                ]
            )
            for patcher in (
                mock.patch.object(Blast, "urlopen", side_effect=responses),
                mock.patch.object(NCBIWWW, "urlopen", side_effect=responses),
                # We don't want to wait 20 sec for each test... Blast.time is
                # the time module itself, so sleep is mocked only while a test
                # here runs, not for every other test in the process.
                mock.patch.object(Blast.time, "sleep"),
            ):
                patcher.start()
                self.addCleanup(patcher.stop)

    def test_blastp_nr_actin(self):
        # Simple protein blast filtered for rat only, using protein
        # GI:160837788 aka NP_075631.2
        # the actin related protein 2/3 complex, subunit 1B [Mus musculus]
        self.run_qblast(
            "blastp",
            "nr",
            "NP_075631.2",
            0.001,
            "rat [ORGN]",
            {"megablast": "FALSE"},
            [
                "NP_112408.1",
                "AAH59131.1",
                "EDM14357.1",
                "NP_001008766.1",
                "NP_001102411.1",
                "EDL80109.1",
                "EDL80106.1",
                "NP_001100434.1",
                "AAI67084.1",
            ],
        )

    def test_pcr_primers(self):
        # This next example finds PCR primer matches in Chimpanzees, e.g. BRCA1:
        self.run_qblast(
            "blastn",
            "nr",
            "GTACCTTGATTTCGTATTC" + ("N" * 30) + "GACTCTACTACCTTTACCC",
            10,
            "pan [ORGN]",
            {"megablast": "FALSE"},
            [
                "XM_034941187.1",
                "XM_034941186.1",
                "XM_034941185.1",
                "XM_034941184.1",
                "XM_034941183.1",
                "XM_034941182.1",
                "XM_034941180.1",
                "XM_034941179.1",
                "XM_034941178.1",
                "XM_034941177.1",
            ],
        )

    def test_orchid_est(self):
        # Try an orchid EST (nucleotide) sequence against NR using BLASTX
        self.run_qblast(
            "blastx",
            "nr",
            """>gi|116660609|gb|EG558220.1|EG558220 CR02019H04 Leaf CR02 cDNA library Catharanthus roseus cDNA clone CR02019H04 5', mRNA sequence
               CTCCATTCCCTCTCTATTTTCAGTCTAATCAAATTAGAGCTTAAAAGAATGAGATTTTTAACAAATAAAA
               AAACATAGGGGAGATTTCATAAAAGTTATATTAGTGATTTGAAGAATATTTTAGTCTATTTTTTTTTTTT
               TCTTTTTTTGATGAAGAAAGGGTATATAAAATCAAGAATCTGGGGTGTTTGTGTTGACTTGGGTCGGGTG
               TGTATAATTCTTGATTTTTTCAGGTAGTTGAAAAGGTAGGGAGAAAAGTGGAGAAGCCTAAGCTGATATT
               GAAATTCATATGGATGGAAAAGAACATTGGTTTAGGATTGGATCAAAAAATAGGTGGACATGGAACTGTA
               CCACTACGTCCTTACTATTTTTGGCCGAGGAAAGATGCTTGGGAAGAACTTAAAACAGTTTTAGAAAGCA
               AGCCATGGATTTCTCAGAAGAAAATGATTATACTTCTTAATCAGGCAACTGATATTATCAATTTATGGCA
               GCAGAGTGGTGGCTCCTTGTCCCAGCAGCAGTAATTACTTTTTTTTCTCTTTTTGTTTCCAAATTAAGAA
               ACATTAGTATCATATGGCTATTTGCTCAATTGCAGATTTCTTTCTTTTGTGAATG""",
            0.0000001,
            None,
            {"megablast": "FALSE"},
            [
                "XP_021665344.1",
                "XP_021615158.1",
                "XP_017223689.1",
                "OMP06800.1",
                "XP_021634873.1",
                "XP_021299673.1",
                "XP_002311451.2",
                "XP_021976565.1",
                "OMO90244.1",
            ],
        )

    def test_discomegablast(self):
        self.run_qblast(
            "blastn",
            "nr",
            """>some sequence
               ATGAAGATCTTCCAGATCCAGTGCAGCAGCTTCAAGGAGAGCAGGTGGCAGAAGAGCAAGTGCGACAACT
               GCCTGAAGTTCCACATCGACATCAACAACAACAGCAAGACCAGCAACACCGACACCGACTTCGACGCCAA
               CACCAACATCAACAGCAACATCAACAGCAACATCAACAGCAACATCAACATCAACAACAGCGGCAACAAC
               AACAAGAACAGCAACAACATCGAGATCACCGAGAACATCGACAACAAGGCCAAGATCATCAACAAGCACA
               TCAAGACCATCACCAACAGCAAGCCCATCCCCATCCCCATCCCCACCCCCACCCCCATCAGCATCAAGGA
               GAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGATG
               AAGAGCACCATCAACCTGAGGAGCGAGGACACCACCAGCAACAAGAGCACCATCGTGTTCACCGAGTGCC
               TGGAGTACAAGGGCCACCAGTGGAGGCCCAACATCTGCGTGACCTGCTTCAGCCCCAAGAACAAGCACAA
               GAACGTGCTGCCCGAGACCAGCACCCCCCTGATCAGCCAGAGCAGCCAGACCAGCACCATCACCCCCAGC
               AGCAGCAGCACCAGCACCAGCACCAGCAGCATCAGCACCCACAAGACCGCCAACAACAAGACCGTGATCA
               CCTACATCAGCAGCACCACCACCACCACCACCACCAGCAGCAGCAGCAGCAGCCCCCCCAGCAGCAGCAT
               CGCCGGCATCACCAACCCCACCAGCAGGAGCAGCAGCCCCATCCTGAAGAGCGTGCCCCCCAGCGCCTAC
               AGCAACGTGGTGATCCCCATCAACAACATCAACAACAGCAACAGCAACAGCAGCAGCGGCGGCGGCAACA
               ACAACAACAAGAGCATCAGCACCCCCAGCAGCCCCATCATCAGCAGGCCCATCACCAACAAGATCAACAA
               CAACAACAACAACAACCAGCCCCAGCTGCACTACAACCAGCCCCAGAGCAGCAGCGTGAGCACCACCAGC
               AGCCCCATCATCAGGCCCGTGCTGAGGAGGCAGTTCCAGAGCTTCCCCAGCAACCCCAAGATCAGCAAGG
               CCATCCTGGAGCAGTGCAACATCATCAACAACAACAGCAACAGCAACAACAGCAACAACAAGGACCCCGT
               GATCCTGTGCAAGTACACCATCGAGAGCCAGCCCAAGAGCAACATCAGCGTGCTGAAGCCCACCCTGGTG
               GAGTTCATCAACCAGCCCGACAGCAAGGACGACGAGAGCAGCGTGAAGAGCCCCCCCCTGCCCGTGGAGA
               GCCAGCCCATCTTCAACAGCAAGCAGAGCGCCACCATGGACGGCATCACCACCCACAAGAGCGTGAGCAT
               CACCATCAGCACCAGCACCAGCCCCAGCAGCACCACCACCACCACCAGCACCACCACCAGCATCATCGCC
               GAGGAGCCCAGCAGCCCCATCCTGCCCACCGCCAGCCCCAGCAGCAGCAGCAGCAGCATCATCACCACCG
               CCACCGCCAGCACCATCCCCATGAGCCCCAGCCTGCCCAGCATCCCCTTCCACGAGTTCGAGACCATGGA
               GAGCAGCACCACCACCACCCTGCTGAGCGAGAACAACGGCGGCGGCGGCGGCAGCAGCTGCAACGACAAC
               AGCAGGAGGAACAGCCTGAACATCCTGCCCCTGAGGCTGAAGAGCTTCAGCTTCAGCGCCCCCCAGAGCG
               ACAGCATGATCGAGCAGCCCGAGGACGACCCCTTCTTCGACTTCGAGGACCTGAGCGACGACGACGACAG
               CAACGACAACGACGACGAGGAGCTGAAGGAGATCAACGGCGAGAAGATCATCCAGCAGAACGACCTGACC
               CCCACCACCACCATCACCAGCACCACCACCATCCTGCAGAGCCCCACCCTGGAGAAGACCCTGAGCACCA
               CCACCACCACCACCATCCCCAGCCCCAGCACCAACAGCAGGAGCATCTGCAACACCCTGATGGACAGCAC
               CGACAGCATCAACAACACCAACACCAACACCAACACCAACACCAACACCAACACCAACACCAACACCAAC
               ACCAACACCAACACCAACACCAACGCCAACATCAACAACAAGGTGAGCACCACCACCACCACCACCACCA
               CCAAGAGGAGGAGCCTGAAGATGGACCAGTTCAAGGAGAAGGAGGACGAGTGGGACCAGGGCGTGGACCT
               GACCAGCTTCCTGAAGAGGAAGCCCACCCTGCAGAGGGACTTCAGCTACTGCAACAACAAGGTGATGGAG
               ATCAGCAGCGTGAAGGAGGAGGCCAAGAGGCTGCACGGCGGCACCGGCTACATCCACCAGTTCGCCTTCG
               AGGCCTTCAAGGACATCCTGGAGGCCAAGCAGACCCAGATCAACAGGGCCTTCTGCAGCCAGAAGATCGA
               CGCCCCCGACTGCGAGATGCTGATCAACGAGATCAACACCGCCAAGAAGCTGCTGGAGGACCTGCTGGAG
               CTGAACAGCAACAGCAGCGGCAGCGGCAACAACAGCAACGACAACAGCGGCAGCAGCAGCCCCAGCAGCA
               GCAAGACCAACACCCTGAACCAGCAGAGCATCTGCATCAAGAGCGAGATCCAACGATACGTTGAAATTCG
               CTTGTGTGCCACTGGTAAATCCACCCCCCCTAAGCCTCTAATAGGGAGACCTTAG""",
            0.0000001,
            None,
            {"template_type": 0, "template_length": 18, "megablast": "on"},  # noqa 231
            ["XM_635681.1", "XM_008496783.1"],
        )

    def run_qblast(
        self,
        program,
        database,
        query,
        e_value,
        entrez_filter,
        additional_args,
        expected_hits,
    ):
        """Do qblast searches with given parameters and analyze results."""
        try:
            handle = Blast.qblast(
                program,
                database,
                query,
                alignments=10,
                descriptions=10,
                hitlist_size=10,
                entrez_query=entrez_filter,
                expect=e_value,
                **additional_args,
            )
        except HTTPError:
            # e.g. a proxy error
            raise MissingExternalDependencyError("internet connection failed") from None

        with handle:
            record = Blast.read(handle)

        # Check the recorded input parameters agree with those requested
        self.assertEqual(float(record.param["expect"]), e_value)
        self.assertEqual(record.program.lower(), program)
        self.assertLessEqual(len(record), 10)

        # Check the expected result(s) are found in the hits
        self.assertGreater(len(record), 0)  # Expected some hits!
        hit_ids = [hit.target.id for hit in record]
        self.assertTrue(
            any(
                expected_hit in hit_id.split("|")
                for expected_hit in expected_hits
                for hit_id in hit_ids
            ),
            f"Missing all expected hits ({', '.join(expected_hits)}), "
            f"instead have: {', '.join(hit_ids)}",
        )

    def test_deprecated_ncbiwww_ncbixml(self):
        """Smoke test of the deprecated NCBIWWW.qblast and NCBIXML.read.

        This always runs against mock_actin.xml (see ALWAYS_MOCKED), as
        test_blastp_nr_actin already runs this search online.
        """
        handle = NCBIWWW.qblast(
            "blastp",
            "nr",
            "NP_075631.2",
            hitlist_size=10,
            entrez_query="rat [ORGN]",
            expect=0.001,
        )
        record = NCBIXML.read(handle)
        handle.close()
        self.assertEqual(record.application.lower(), "blastp")
        self.assertEqual(float(record.expect), 0.001)
        self.assertIn(
            "NP_112408.1",
            [part for a in record.alignments for part in a.hit_id.split("|")],
        )

    def test_parse_qblast_ref_page(self):
        with open(support.DATA / "Blast" / "html_msgid_29_blastx_001.html", "rb") as f:
            handle = BytesIO(f.read())
        self.assertRaises(ValueError, Blast._parse_qblast_ref_page, handle)

    def test_short_query(self):
        """Test SHORT_QUERY_ADJUST parameter."""
        # Should give no hits:
        with Blast.qblast("blastp", "nr", "ICWENRMP", hitlist_size=5) as my_search:
            my_hits = Blast.read(my_search)
        self.assertEqual(len(my_hits), 0)

        # Should give hits:
        with Blast.qblast(
            "blastp", "nr", "ICWENRMP", hitlist_size=5, short_query=True
        ) as my_search:
            my_hits = Blast.read(my_search)
        self.assertNotEqual(len(my_hits), 0)

        # Disabled for 1.81 release
        # Query does not seem to be returning anything so
        # maybe the test needs to be updated?
        #
        # with warnings.catch_warnings(record=True) as w:
        #     # Cause all warnings to always be triggered.
        #     warnings.simplefilter("always")
        #     # Trigger a warning.
        #     my_search = Blast.qblast(
        #         "blastn", "nt", "ATGTCAACTTCAGAA", hitlist_size=5, short_query=True
        #     )
        #     # Verify some things
        #     self.assertEqual(w[0].category, BiopythonWarning)
        #     self.assertIn("blastn", str(w[0].message))
        #     my_hits = Blast.read(my_search)
        #     my_search.close()
        #     self.assertNotEqual(len(my_hits), 0)

    def test_error_conditions(self):
        """Test if exceptions were properly handled."""
        self.assertRaises(
            ValueError,
            Blast.qblast,
            "megablast",
            "nt",
            "ATGCGTACGCAGCTAAAGTAAACCTATCGCGTCTCCT",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
