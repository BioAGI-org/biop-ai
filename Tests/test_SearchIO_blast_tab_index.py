# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO blast-tab indexing."""

import unittest
from unittest import mock

import support
from search_tests_common import CheckIndex
from search_tests_common import CheckRaw
from search_tests_common import EndOfFileGuard

from Bio import SearchIO


class BlastTabRawCases(CheckRaw):
    """Check BLAST tabular get_raw method."""

    fmt = "blast-tab"

    def test_multi_block_bgzf(self):
        """Test get_raw on a BGZF copy with many blocks, as on the plain file."""
        self.check_raw_multi_block_bgzf(
            support.DATA / "Blast" / "tab_2226_tblastn_001.txt"
        )
        self.check_raw_multi_block_bgzf(
            support.DATA / "Blast" / "tab_2226_tblastn_005.txt", comments=True
        )

    def test_blasttab_2226_multiple_first(self):
        """Test blast-tab raw string retrieval, BLAST 2.2.26+, multiple queries, first (tab_2226_tblastn_001.txt)."""
        filename = support.DATA / "Blast" / "tab_2226_tblastn_001.txt"
        raw = """gi|16080617|ref|NP_391444.1|	gi|145479850|ref|XM_001425911.1|	34.88	43	28	0	31	73	1744	1872	1e-05	34.7
gi|16080617|ref|NP_391444.1|	gi|72012412|ref|XM_777959.1|	33.90	59	31	1	44	94	1057	1233	1e-04	31.6
gi|16080617|ref|NP_391444.1|	gi|115975252|ref|XM_001180111.1|	33.90	59	31	1	44	94	1057	1233	1e-04	31.6
"""
        self.check_raw(filename, "gi|16080617|ref|NP_391444.1|", raw)

    def test_blasttab_2226_multiple_last(self):
        """Test blast-tab raw string retrieval, BLAST 2.2.26+, multiple queries, last (tab_2226_tblastn_001.txt)."""
        filename = support.DATA / "Blast" / "tab_2226_tblastn_001.txt"
        raw = """gi|11464971:4-101	gi|350596019|ref|XM_003360601.2|	95.92	98	4	0	1	98	95	388	2e-67	 199
gi|11464971:4-101	gi|350596019|ref|XM_003360601.2|	29.58	71	46	2	30	96	542	754	4e-05	32.7
gi|11464971:4-101	gi|301779869|ref|XM_002925302.1|	97.96	98	2	0	1	98	78	371	2e-67	 202
gi|11464971:4-101	gi|301779869|ref|XM_002925302.1|	30.00	100	64	2	3	96	804	1103	3e-09	45.1
gi|11464971:4-101	gi|296223671|ref|XM_002757683.1|	97.96	98	2	0	1	98	161	454	4e-67	 202
gi|11464971:4-101	gi|296223671|ref|XM_002757683.1|	30.00	100	64	2	3	96	866	1165	3e-09	45.1
gi|11464971:4-101	gi|338714227|ref|XM_001492113.3|	97.96	98	2	0	1	98	173	466	2e-66	 202
gi|11464971:4-101	gi|338714227|ref|XM_001492113.3|	31.00	100	63	2	3	96	899	1198	1e-09	46.6
gi|11464971:4-101	gi|365982352|ref|XM_003667962.1|	30.77	52	27	1	12	54	3181	3336	1.7	19.6
"""
        self.check_raw(filename, "gi|11464971:4-101", raw)

    def test_blasttab_2226_single(self):
        """Test blast-tab raw string retrieval, BLAST 2.2.26+, single query (tab_2226_tblastn_004.txt)."""
        filename = support.DATA / "Blast" / "tab_2226_tblastn_004.txt"
        raw = """gi|11464971:4-101	gi|350596019|ref|XM_003360601.2|	95.92	98	4	0	1	98	95	388	2e-67	 199
gi|11464971:4-101	gi|350596019|ref|XM_003360601.2|	29.58	71	46	2	30	96	542	754	4e-05	32.7
gi|11464971:4-101	gi|301779869|ref|XM_002925302.1|	97.96	98	2	0	1	98	78	371	2e-67	 202
gi|11464971:4-101	gi|301779869|ref|XM_002925302.1|	30.00	100	64	2	3	96	804	1103	3e-09	45.1
gi|11464971:4-101	gi|296223671|ref|XM_002757683.1|	97.96	98	2	0	1	98	161	454	4e-67	 202
gi|11464971:4-101	gi|296223671|ref|XM_002757683.1|	30.00	100	64	2	3	96	866	1165	3e-09	45.1
gi|11464971:4-101	gi|338714227|ref|XM_001492113.3|	97.96	98	2	0	1	98	173	466	2e-66	 202
gi|11464971:4-101	gi|338714227|ref|XM_001492113.3|	31.00	100	63	2	3	96	899	1198	1e-09	46.6
gi|11464971:4-101	gi|365982352|ref|XM_003667962.1|	30.77	52	27	1	12	54	3181	3336	1.7	19.6
"""
        self.check_raw(filename, "gi|11464971:4-101", raw)

    def test_blasttab_2226_multiple_first_commented(self):
        """Test blast-tab raw string retrieval, BLAST 2.2.26+, multiple queries, first, commented (tab_2226_tblastn_005.txt)."""
        filename = support.DATA / "Blast" / "tab_2226_tblastn_005.txt"
        raw = """# TBLASTN 2.2.26+
# Query: random_s00
# Database: db/minirefseq_mrna
# 0 hits found
"""
        self.check_raw(filename, "random_s00", raw, comments=True)

    def test_blasttab_2226_multiple_middle_commented(self):
        """Test blast-tab raw string retrieval, BLAST 2.2.26+, multiple queries, middle, commented (tab_2226_tblastn_005.txt)."""
        filename = support.DATA / "Blast" / "tab_2226_tblastn_005.txt"
        raw = """# TBLASTN 2.2.26+
# Query: gi|16080617|ref|NP_391444.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168]
# Database: db/minirefseq_mrna
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 3 hits found
gi|16080617|ref|NP_391444.1|	gi|145479850|ref|XM_001425911.1|	34.88	43	28	0	31	73	1744	1872	1e-05	34.7
gi|16080617|ref|NP_391444.1|	gi|72012412|ref|XM_777959.1|	33.90	59	31	1	44	94	1057	1233	1e-04	31.6
gi|16080617|ref|NP_391444.1|	gi|115975252|ref|XM_001180111.1|	33.90	59	31	1	44	94	1057	1233	1e-04	31.6
"""
        self.check_raw(filename, "gi|16080617|ref|NP_391444.1|", raw, comments=True)

    def test_blasttab_2226_multiple_last_commented(self):
        """Test blast-tab raw string retrieval, BLAST 2.2.26+, multiple queries, last, commented (tab_2226_tblastn_005.txt)."""
        filename = support.DATA / "Blast" / "tab_2226_tblastn_005.txt"
        raw = """# TBLASTN 2.2.26+
# Query: gi|11464971:4-101 pleckstrin [Mus musculus]
# Database: db/minirefseq_mrna
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 9 hits found
gi|11464971:4-101	gi|350596019|ref|XM_003360601.2|	95.92	98	4	0	1	98	95	388	2e-67	 199
gi|11464971:4-101	gi|350596019|ref|XM_003360601.2|	29.58	71	46	2	30	96	542	754	4e-05	32.7
gi|11464971:4-101	gi|301779869|ref|XM_002925302.1|	97.96	98	2	0	1	98	78	371	2e-67	 202
gi|11464971:4-101	gi|301779869|ref|XM_002925302.1|	30.00	100	64	2	3	96	804	1103	3e-09	45.1
gi|11464971:4-101	gi|296223671|ref|XM_002757683.1|	97.96	98	2	0	1	98	161	454	4e-67	 202
gi|11464971:4-101	gi|296223671|ref|XM_002757683.1|	30.00	100	64	2	3	96	866	1165	3e-09	45.1
gi|11464971:4-101	gi|338714227|ref|XM_001492113.3|	97.96	98	2	0	1	98	173	466	2e-66	 202
gi|11464971:4-101	gi|338714227|ref|XM_001492113.3|	31.00	100	63	2	3	96	899	1198	1e-09	46.6
gi|11464971:4-101	gi|365982352|ref|XM_003667962.1|	30.77	52	27	1	12	54	3181	3336	1.7	19.6
"""
        self.check_raw(filename, "gi|11464971:4-101", raw, comments=True)

    def test_blasttab_2226_single_commented(self):
        """Test blast-tab raw string retrieval, BLAST 2.2.26+, single query, commented (tab_2226_tblastn_008.txt)."""
        filename = support.DATA / "Blast" / "tab_2226_tblastn_008.txt"
        raw = """# TBLASTN 2.2.26+
# Query: gi|11464971:4-101 pleckstrin [Mus musculus]
# Database: db/minirefseq_mrna
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 9 hits found
gi|11464971:4-101	gi|350596019|ref|XM_003360601.2|	95.92	98	4	0	1	98	95	388	2e-67	 199
gi|11464971:4-101	gi|350596019|ref|XM_003360601.2|	29.58	71	46	2	30	96	542	754	4e-05	32.7
gi|11464971:4-101	gi|301779869|ref|XM_002925302.1|	97.96	98	2	0	1	98	78	371	2e-67	 202
gi|11464971:4-101	gi|301779869|ref|XM_002925302.1|	30.00	100	64	2	3	96	804	1103	3e-09	45.1
gi|11464971:4-101	gi|296223671|ref|XM_002757683.1|	97.96	98	2	0	1	98	161	454	4e-67	 202
gi|11464971:4-101	gi|296223671|ref|XM_002757683.1|	30.00	100	64	2	3	96	866	1165	3e-09	45.1
gi|11464971:4-101	gi|338714227|ref|XM_001492113.3|	97.96	98	2	0	1	98	173	466	2e-66	 202
gi|11464971:4-101	gi|338714227|ref|XM_001492113.3|	31.00	100	63	2	3	96	899	1198	1e-09	46.6
gi|11464971:4-101	gi|365982352|ref|XM_003667962.1|	30.77	52	27	1	12	54	3181	3336	1.7	19.6
"""
        self.check_raw(filename, "gi|11464971:4-101", raw, comments=True)

    def test_blasttab_2226_many_hits(self):
        """Test blast-tab raw string retrieval, BLAST 2.2.26+, query with many hits (tab_2226_tblastn_001.txt)."""
        filename = support.DATA / "Blast" / "tab_2226_tblastn_001.txt"
        # The query's first hit line, 40000 times
        self.check_raw_repeated(
            filename,
            "gi|16080617|ref|NP_391444.1|",
            b"gi|16080617|ref|NP_391444.1|\tgi|145479850|",
            b"\n",
            40000,
        )

    def test_blasttab_2226_many_hits_commented(self):
        """Test blast-tab raw string retrieval, BLAST 2.2.26+, query with many hits, commented (tab_2226_tblastn_005.txt)."""
        filename = support.DATA / "Blast" / "tab_2226_tblastn_005.txt"
        # The query's first hit line, 40000 times
        self.check_raw_repeated(
            filename,
            "gi|16080617|ref|NP_391444.1|",
            b"gi|16080617|ref|NP_391444.1|\tgi|145479850|",
            b"\n",
            40000,
            comments=True,
        )


class BlastTabIndexCases(CheckIndex):
    fmt = "blast-tab"

    def test_blasttab_2226_tblastn_001(self):
        """Test blast-tab indexing, BLAST 2.2.26+, multiple queries."""
        filename = support.DATA / "Blast" / "tab_2226_tblastn_001.txt"
        self.check_index(filename, self.fmt)

    def test_blasttab_2226_tblastn_002(self):
        """Test blast-tab indexing, BLAST 2.2.26+, single query, no hits."""
        filename = support.DATA / "Blast" / "tab_2226_tblastn_002.txt"
        self.check_index(filename, self.fmt)

    def test_blasttab_2226_tblastn_004(self):
        """Test blast-tab indexing, BLAST 2.2.26+, single query, multiple hits."""
        filename = support.DATA / "Blast" / "tab_2226_tblastn_004.txt"
        self.check_index(filename, self.fmt)

    def test_blasttab_2226_tblastn_005(self):
        """Test blast-tab indexing, BLAST 2.2.26+, multiple queries, commented."""
        filename = support.DATA / "Blast" / "tab_2226_tblastn_005.txt"
        self.check_index(filename, self.fmt, comments=True)

    def test_blasttab_2226_tblastn_006(self):
        """Test blast-tab indexing, BLAST 2.2.26+, single query, no hits, commented."""
        filename = support.DATA / "Blast" / "tab_2226_tblastn_006.txt"
        self.check_index(filename, self.fmt, comments=True)

    def test_blasttab_comment_sing(self):
        """Test blast-tab indexing, BLAST 2.2.26+, single query, multiple hits, commented."""
        filename = support.DATA / "Blast" / "tab_2226_tblastn_008.txt"
        self.check_index(filename, self.fmt, comments=True)

    def test_blasttab_2226_tblastn_009(self):
        """Test blast-tab indexing, BLAST 2.2.26+, custom columns."""
        filename = support.DATA / "Blast" / "tab_2226_tblastn_009.txt"
        self.check_index(filename, self.fmt, fields=["qseqid", "sseqid"])

    def test_blasttab_2226_tblastn_010(self):
        """Test blast-tab indexing, BLAST 2.2.26+, custom columns, commented."""
        filename = support.DATA / "Blast" / "tab_2226_tblastn_010.txt"
        self.check_index(filename, self.fmt, comments=True)

    def test_blasttab_2226_tblastn_011(self):
        """Test blast-tab indexing, BLAST 2.2.26+, all columns, commented."""
        filename = support.DATA / "Blast" / "tab_2226_tblastn_011.txt"
        self.check_index(filename, self.fmt, comments=True)


class TruncatedCases(unittest.TestCase):
    """A file cut short after indexing is read back without looping."""

    def index(self, data):
        """Index data as blast-tab through EndOfFileGuard; return both."""
        handle = EndOfFileGuard(data)
        with mock.patch(
            "Bio.SearchIO._index._open_for_random_access", return_value=handle
        ):
            index = SearchIO.index("test.txt", "blast-tab")
        self.addCleanup(index.close)
        return index, handle

    def test_query_after_end_of_file(self):
        """Test reading back a query that starts after the new end of file."""
        # A file can shrink under its index: truncated in place while a
        # SearchIO.index dict is open, or replaced after index_db wrote its
        # offsets to disk.
        filename = support.DATA / "Blast" / "tab_2226_tblastn_001.txt"
        with open(filename, "rb") as handle:
            data = handle.read()
        index, handle = self.index(data)
        key = "gi|11464971:4-101"
        self.assertEqual(list(index), ["gi|16080617|ref|NP_391444.1|", key])
        offset = data.index(b"\n" + key.encode() + b"\t") + 1
        # Cut the file after its first line, well before the last query.
        handle.truncate(data.index(b"\n") + 1)
        for name, get in (("get_raw", index.get_raw), ("getitem", index.__getitem__)):
            with self.subTest(name):
                with self.assertRaises(ValueError) as cm:
                    get(key)
                self.assertEqual(
                    str(cm.exception), f"File ended before the query at offset {offset}"
                )

    def test_empty_query_id_at_end_of_file(self):
        """Test reading back a query with an empty ID that ends the file."""
        # The first query's ID column is blanked, and the file is cut after
        # indexing where that query ends. At end of file the ID column of
        # the empty "line" is empty too, so the reader must stop there, not
        # take it as more of this query.
        filename = support.DATA / "Blast" / "tab_2226_tblastn_001.txt"
        with open(filename, "rb") as handle:
            data = handle.read()
        data = data.replace(b"gi|16080617|ref|NP_391444.1|\t", b"\t")
        index, handle = self.index(data)
        key = "gi|11464971:4-101"
        self.assertEqual(list(index), ["", key])
        # Cut the file where the query with the empty ID ends.
        end = data.index(b"\n" + key.encode() + b"\t") + 1
        handle.truncate(end)
        self.assertEqual(index.get_raw(""), data[:end])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
