# Copyright 2005 by Iddo Friedberg.  All rights reserved.
# Revisions copyright 2006-2013,2017 by Peter Cock. All rights reserved.
# Revisions copyright 2008 by Frank Kauff. All rights reserved.
# Revisions copyright 2009 by Michiel de Hoon. All rights reserved.
# Revisions copyright 2015 by Joe Cora. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for Nexus module."""

import json
import os.path
import platform
import subprocess
import sys
import tempfile
import unittest
from io import StringIO
from unittest import mock

import support

import Bio
from Bio import Phylo
from Bio import SeqIO
from Bio.Align import Alignment
from Bio.Align import MultipleSeqAlignment
from Bio.Align.nexus import AlignmentIterator
from Bio.Align.nexus import AlignmentWriter
from Bio.AlignIO.NexusIO import NexusIterator
from Bio.AlignIO.NexusIO import NexusWriter
from Bio.Nexus import Nexus
from Bio.Nexus import Trees
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

try:
    from Bio.Nexus import cnexus
except ImportError:
    cnexus = None


class OldSelfTests(unittest.TestCase):
    """Test cases originally in Nexus.py via __main__."""

    def test_trees_and_taxa_block(self):
        """Basic tree file with TREES and TAXA block."""
        nexus1 = Nexus.Nexus()
        nexus1.read(support.DATA / "Nexus" / "bats.nex")

    def test_data_and_codons_block(self):
        """Simple sequence data file with DATA and CODONS block."""
        nexus2 = Nexus.Nexus()
        nexus2.read(support.DATA / "Nexus" / "codonposset.nex")

    def test_data_sets_trees_unknown_block(self):
        """Sequence data file with DATA, SETS, TREES and an unknown block."""
        nexus3 = Nexus.Nexus()
        nexus3.read(support.DATA / "Nexus" / "test_Nexus_input.nex")

    def test_taxa_and_characters_block(self):
        """Taxa and characters multi-state block."""
        nexus4 = Nexus.Nexus()
        nexus4.read(support.DATA / "Nexus" / "vSysLab_Ganaspidium_multistate.nex")

    def test_taxa_and_characters_with_many_codings_one_without_state(self):
        """Taxa and chr blocks, over 9 codings, 1 character without states."""
        nexus5 = Nexus.Nexus()
        nexus5.read(
            support.DATA / "Nexus" / "vSysLab_Heptascelio_no-states_10+chars.nex"
        )

    def test_taxa_and_characters_with_many_codings_two_without_state(self):
        """Taxa and chr blocks, over 9 codings, 2 character without states."""
        nexus6 = Nexus.Nexus()
        # TODO: Implement continuous datatype:
        # Bio.Nexus.StandardData.NexusError: Unsupported datatype: continuous
        self.assertRaises(
            Nexus.NexusError,
            nexus6.read,
            support.DATA / "Nexus" / "vSysLab_Oreiscelio_discrete+continuous.nex",
        )

    def test_standard_data_error_is_nexus_error(self):
        """A bad standard-datatype coding raises Nexus.NexusError."""
        nexus = Nexus.Nexus()
        with self.assertRaisesRegex(Nexus.NexusError, "Improper character"):
            nexus.read(
                "#NEXUS\nbegin data;\ndimensions ntax=2 nchar=3;\n"
                'format datatype=standard symbols="01";\n'
                "matrix\na 01)\nb 010\n;\nend;\n"
            )


class NexusTest1(unittest.TestCase):
    def setUp(self):
        self.testfile_dir = support.DATA / "Nexus"
        self.handle = open(os.path.join(self.testfile_dir, "test_Nexus_input.nex"))

    def tearDown(self):
        self.handle.close()

    def test_WriteToFileName(self):
        """Test writing to a given filename."""
        filename = support.DATA / "Nexus" / "test_temp.nex"
        if os.path.isfile(filename):
            os.remove(filename)
        n = Nexus.Nexus(self.handle)
        n.write_nexus_data(filename)
        self.assertTrue(os.path.isfile(filename))
        os.remove(filename)

    def test_write_with_dups(self):
        # see issue: biopython/Bio/Nexus/Nexus.py _unique_label() eval error #633
        records = [
            SeqRecord(
                Seq("ATGCTGCTGAT"), id="foo", annotations={"molecule_type": "DNA"}
            )
            for _ in range(4)
        ]
        out_file = StringIO()
        self.assertEqual(4, SeqIO.write(records, out_file, "nexus"))

    def test_NexusTest1(self):
        """Test Nexus module."""
        # check data of main nexus file
        n = Nexus.Nexus(self.handle)
        self.assertEqual(
            os.path.normpath(n.filename),
            os.path.normpath(support.DATA / "Nexus" / "test_Nexus_input.nex"),
        )
        self.assertEqual(n.ntax, 9)
        self.assertEqual(n.nchar, 48)
        self.assertEqual(n.datatype, "dna")
        self.assertTrue(n.interleave)
        self.assertEqual(n.missing, "?")
        self.assertEqual(n.gap, "-")
        self.assertEqual(
            n.taxlabels,
            [
                "t1",
                "t2 the name",
                "isn'that [a] strange name?",
                "one should be punished, for (that)!",
                "t5",
                "t6",
                "t7",
                "t8",
                "t9",
            ],
        )
        self.assertEqual(
            n.charlabels,
            {
                0: "a",
                1: "b",
                2: "c",
                4: "f",
                9: "A",
                10: "B",
                22: "x",
                23: "y",
                29: "1,2,3 can't decide for a name?!",
                47: "final",
            },
        )
        self.assertEqual(
            n.charsets,
            {
                "big": [0, 2, 4, 6],
                "bigchunk": [
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    16,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                    24,
                    25,
                    26,
                    27,
                    28,
                    29,
                    30,
                    31,
                    32,
                    33,
                    34,
                    35,
                    36,
                    37,
                    38,
                    39,
                    40,
                    41,
                    42,
                    43,
                    44,
                    45,
                    46,
                ],
                "byname": [
                    0,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    16,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                    24,
                    25,
                    26,
                    27,
                    28,
                    29,
                ],
                "c1": [0, 1, 2, 3, 4, 5, 6, 7],
                "c2": [8, 9, 10, 11, 12, 13, 14, 15],
                "c3": [16, 17, 18, 19, 20, 21, 22, 23],
                "firsthalf": [
                    0,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    16,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                ],
                "mix": [0, 1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46],
                "mux": [
                    0,
                    1,
                    4,
                    7,
                    8,
                    10,
                    13,
                    16,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                    25,
                    28,
                    31,
                    34,
                    37,
                    40,
                    43,
                    46,
                ],
                "pos1": [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45],
                "pos2": [1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46],
                "pos3": [2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47],
                "secondhalf": [
                    24,
                    25,
                    26,
                    27,
                    28,
                    29,
                    30,
                    31,
                    32,
                    33,
                    34,
                    35,
                    36,
                    37,
                    38,
                    39,
                    40,
                    41,
                    42,
                    43,
                    44,
                    45,
                    46,
                    47,
                ],
            },
        )
        self.assertEqual(
            n.taxsets,
            {
                "normal": [
                    "isn'that [a] strange name?",
                    "one should be punished, for (that)!",
                    "t1",
                    "t5",
                    "t6",
                    "t8",
                ],
                "reference": [
                    "isn'that [a] strange name?",
                    "one should be punished, for (that)!",
                    "t1",
                    "t2 the name",
                    "t5",
                    "t6",
                ],
                "tbyname1": [
                    "isn'that [a] strange name?",
                    "one should be punished, for (that)!",
                    "t1",
                    "t2 the name",
                    "t5",
                    "t6",
                ],
                "tbyname2": [
                    "isn'that [a] strange name?",
                    "one should be punished, for (that)!",
                    "t2 the name",
                    "t5",
                    "t6",
                    "t7",
                ],
                "tbyname3": ["t1", "t2 the name"],
            },
        )
        self.assertEqual(len(n.charpartitions), 2)
        self.assertIn("codons", n.charpartitions)
        self.assertIn("part", n.charpartitions)
        self.assertEqual(
            n.charpartitions["codons"],
            {
                "a": [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45],
                "b": [1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46],
                "c": [2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47],
            },
        )
        self.assertEqual(
            n.charpartitions["part"],
            {
                "one": [0, 1, 2, 3, 4, 5, 6, 7],
                "three": [16, 17, 18, 19, 20, 21, 22, 23],
                "two": [8, 9, 10, 11, 12, 13, 14, 15],
            },
        )
        self.assertEqual(list(n.taxpartitions), ["taxpart"])
        self.assertEqual(
            n.taxpartitions["taxpart"],
            {
                "badnames": [
                    "isn'that [a] strange name?",
                    "one should be punished, for (that)!",
                    "t2 the name",
                ],
                "goodnames": ["t1", "t5", "t6", "t7", "t8", "t9"],
            },
        )

        # now we check excluding characters, deleting taxa,
        # and exporting adjusted sets
        f1 = tempfile.NamedTemporaryFile("w+")
        n.write_nexus_data(f1, delete=["t1", "t7"], exclude=n.invert(n.charsets["big"]))
        f1.seek(0)
        nf1 = Nexus.Nexus(f1)
        self.assertEqual(os.path.normpath(nf1.filename), os.path.normpath(f1.name))
        self.assertEqual(nf1.ntax, 7)
        self.assertEqual(nf1.nchar, 4)
        self.assertEqual(nf1.datatype, "dna")
        self.assertFalse(nf1.interleave)
        self.assertEqual(nf1.missing, "?")
        self.assertEqual(nf1.gap, "-")
        self.assertEqual(
            nf1.taxlabels,
            [
                "t2 the name",
                "isn'that [a] strange name?",
                "one should be punished, for (that)!",
                "t5",
                "t6",
                "t8",
                "t9",
            ],
        )
        self.assertEqual(nf1.charlabels, {0: "a", 1: "c", 2: "f"})
        self.assertEqual(
            nf1.charsets,
            {
                "big": [0, 1, 2, 3],
                "bigchunk": [1, 2, 3],
                "byname": [0, 2, 3],
                "c1": [0, 1, 2, 3],
                "firsthalf": [0, 1, 2, 3],
                "mix": [0, 2],
                "mux": [0, 2],
                "pos1": [0, 3],
                "pos2": [2],
                "pos3": [1],
            },
        )
        self.assertEqual(
            nf1.taxsets,
            {
                "normal": [
                    "isn'that [a] strange name?",
                    "one should be punished, for (that)!",
                    "t5",
                    "t6",
                    "t8",
                ],
                "reference": [
                    "isn'that [a] strange name?",
                    "one should be punished, for (that)!",
                    "t2 the name",
                    "t5",
                    "t6",
                ],
                "tbyname1": [
                    "isn'that [a] strange name?",
                    "one should be punished, for (that)!",
                    "t2 the name",
                    "t5",
                    "t6",
                ],
                "tbyname2": [
                    "isn'that [a] strange name?",
                    "one should be punished, for (that)!",
                    "t2 the name",
                    "t5",
                    "t6",
                ],
                "tbyname3": ["t2 the name"],
            },
        )
        self.assertEqual(len(nf1.charpartitions), 2)
        self.assertIn("codons", nf1.charpartitions)
        self.assertIn("part", nf1.charpartitions)
        self.assertEqual(
            nf1.charpartitions["codons"], {"a": [0, 3], "b": [2], "c": [1]}
        )
        self.assertEqual(nf1.charpartitions["part"], {"one": [0, 1, 2, 3]})

        self.assertEqual(list(nf1.taxpartitions), ["taxpart"])
        self.assertEqual(
            nf1.taxpartitions["taxpart"],
            {
                "badnames": [
                    "isn'that [a] strange name?",
                    "one should be punished, for (that)!",
                    "t2 the name",
                ],
                "goodnames": ["t5", "t6", "t8", "t9"],
            },
        )

        f2 = tempfile.NamedTemporaryFile("w+")
        n.write_nexus_data(f2, delete=["t2_the_name"], exclude=list(range(3, 40, 4)))
        f2.seek(0)
        nf2 = Nexus.Nexus(f2)
        self.assertEqual(os.path.normpath(nf2.filename), os.path.normpath(f2.name))
        self.assertEqual(nf2.ntax, 9)
        self.assertEqual(nf2.nchar, 38)
        self.assertEqual(nf2.datatype, "dna")
        self.assertFalse(nf2.interleave)
        self.assertEqual(nf2.missing, "?")
        self.assertEqual(nf2.gap, "-")
        self.assertEqual(
            nf2.taxlabels,
            [
                "t1",
                "t2 the name",
                "isn'that [a] strange name?",
                "one should be punished, for (that)!",
                "t5",
                "t6",
                "t7",
                "t8",
                "t9",
            ],
        )
        self.assertEqual(
            nf2.charlabels,
            {
                0: "a",
                1: "b",
                2: "c",
                3: "f",
                7: "A",
                8: "B",
                17: "x",
                22: "1,2,3 can't decide for a name?!",
                37: "final",
            },
        )
        self.assertEqual(
            nf2.charsets,
            {
                "big": [0, 2, 3, 5],
                "bigchunk": [
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    16,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                    24,
                    25,
                    26,
                    27,
                    28,
                    29,
                    30,
                    31,
                    32,
                    33,
                    34,
                    35,
                    36,
                ],
                "byname": [
                    0,
                    3,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    16,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                ],
                "c1": [0, 1, 2, 3, 4, 5],
                "c2": [6, 7, 8, 9, 10, 11],
                "c3": [12, 13, 14, 15, 16, 17],
                "firsthalf": [
                    0,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    16,
                    17,
                ],
                "mix": [0, 1, 3, 8, 10, 12, 17, 19, 21, 26, 28, 30, 33, 36],
                "mux": [
                    0,
                    1,
                    3,
                    6,
                    8,
                    10,
                    12,
                    13,
                    14,
                    15,
                    16,
                    17,
                    19,
                    21,
                    26,
                    28,
                    30,
                    33,
                    36,
                ],
                "pos1": [0, 5, 7, 9, 14, 16, 18, 23, 25, 27, 32, 35],
                "pos2": [1, 3, 8, 10, 12, 17, 19, 21, 26, 28, 30, 33, 36],
                "pos3": [2, 4, 6, 11, 13, 15, 20, 22, 24, 29, 31, 34, 37],
                "secondhalf": [
                    18,
                    19,
                    20,
                    21,
                    22,
                    23,
                    24,
                    25,
                    26,
                    27,
                    28,
                    29,
                    30,
                    31,
                    32,
                    33,
                    34,
                    35,
                    36,
                    37,
                ],
            },
        )

        self.assertEqual(
            nf2.taxsets,
            {
                "normal": [
                    "isn'that [a] strange name?",
                    "one should be punished, for (that)!",
                    "t1",
                    "t5",
                    "t6",
                    "t8",
                ],
                "reference": [
                    "isn'that [a] strange name?",
                    "one should be punished, for (that)!",
                    "t1",
                    "t2 the name",
                    "t5",
                    "t6",
                ],
                "tbyname1": [
                    "isn'that [a] strange name?",
                    "one should be punished, for (that)!",
                    "t1",
                    "t2 the name",
                    "t5",
                    "t6",
                ],
                "tbyname2": [
                    "isn'that [a] strange name?",
                    "one should be punished, for (that)!",
                    "t2 the name",
                    "t5",
                    "t6",
                    "t7",
                ],
                "tbyname3": ["t1", "t2 the name"],
            },
        )
        self.assertEqual(len(nf2.charpartitions), 2)
        self.assertIn("codons", nf2.charpartitions)
        self.assertIn("part", nf2.charpartitions)
        self.assertEqual(
            nf2.charpartitions["codons"],
            {
                "a": [0, 5, 7, 9, 14, 16, 18, 23, 25, 27, 32, 35],
                "b": [1, 3, 8, 10, 12, 17, 19, 21, 26, 28, 30, 33, 36],
                "c": [2, 4, 6, 11, 13, 15, 20, 22, 24, 29, 31, 34, 37],
            },
        )
        self.assertEqual(
            nf2.charpartitions["part"],
            {
                "one": [0, 1, 2, 3, 4, 5],
                "three": [12, 13, 14, 15, 16, 17],
                "two": [6, 7, 8, 9, 10, 11],
            },
        )
        self.assertEqual(list(nf2.taxpartitions), ["taxpart"])
        self.assertEqual(
            nf2.taxpartitions["taxpart"],
            {
                "badnames": [
                    "isn'that [a] strange name?",
                    "one should be punished, for (that)!",
                    "t2 the name",
                ],
                "goodnames": ["t1", "t5", "t6", "t7", "t8", "t9"],
            },
        )
        # check the stepmatrix
        self.assertEqual(
            n.weighted_stepmatrix(name="matrix_test"),
            """\
usertype matrix_test stepmatrix=5
        A        C        G        T        -
[A]     .       2.40     2.57     2.43     2.43     
[C]    2.40      .       2.28     2.12     2.14     
[G]    2.57     2.28      .       2.31     2.31     
[T]    2.43     2.12     2.31      .       2.14     
[-]    2.43     2.14     2.31     2.14      .       
;
""",  # noqa : W291
        )

    def test_write_alignment(self):
        # Default causes no interleave (columns <= 1000)
        records = [
            SeqRecord(
                Seq("ATGCTGCTGA" * 90), id=_id, annotations={"molecule_type": "DNA"}
            )
            for _id in ["foo", "bar", "baz"]
        ]
        a = Alignment(records)

        handle = StringIO()
        AlignmentWriter(handle).write([a])
        handle.seek(0)
        data = handle.read()
        self.assertEqual(
            data,
            """\
#NEXUS
begin data;
dimensions ntax=3 nchar=900;
format datatype=dna missing=? gap=-;
matrix
foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
;
end;
""",
        )
        # Default causes interleave (columns > 1000)
        records = [
            SeqRecord(
                Seq("ATGCTGCTGA" * 110), id=_id, annotations={"molecule_type": "DNA"}
            )
            for _id in ["foo", "bar", "baz"]
        ]
        a = Alignment(records)
        handle = StringIO()
        AlignmentWriter(handle).write([a])
        handle.seek(0)
        data = handle.read()
        self.assertEqual(
            data,
            """\
#NEXUS
begin data;
dimensions ntax=3 nchar=1100;
format datatype=dna missing=? gap=- interleave;
matrix
foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA

foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA

foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA

foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA

foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA

foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA

foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA

foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA

foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA

foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA

foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA

foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA

foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA

foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA

foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA

foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA

;
end;
""",
        )

        # Override interleave: True
        records = [
            SeqRecord(
                Seq("ATGCTGCTGA" * 9), id=_id, annotations={"molecule_type": "DNA"}
            )
            for _id in ["foo", "bar", "baz"]
        ]
        a = Alignment(records)
        handle = StringIO()
        AlignmentWriter(handle, interleave=True).write([a])
        handle.seek(0)
        data = handle.read()
        self.assertEqual(
            data,
            """\
#NEXUS
begin data;
dimensions ntax=3 nchar=90;
format datatype=dna missing=? gap=- interleave;
matrix
foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA

foo ATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGA

;
end;
""",
        )

        # Override interleave: False
        records = [
            SeqRecord(
                Seq("ATGCTGCTGA" * 110), id=_id, annotations={"molecule_type": "DNA"}
            )
            for _id in ["foo", "bar", "baz"]
        ]
        a = Alignment(records)
        handle = StringIO()
        AlignmentWriter(handle, interleave=False).write([a])
        handle.seek(0)
        data = handle.read()
        self.assertEqual(
            data,
            """\
#NEXUS
begin data;
dimensions ntax=3 nchar=1100;
format datatype=dna missing=? gap=-;
matrix
foo ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
bar ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
baz ATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGAATGCTGCTGA
;
end;
""",
        )

    def test_write_alignment_msa(self):
        # Default causes no interleave (columns <= 1000)
        records = [
            SeqRecord(
                Seq("ATGCTGCTGA" * 90), id=_id, annotations={"molecule_type": "DNA"}
            )
            for _id in ["foo", "bar", "baz"]
        ]
        a = MultipleSeqAlignment(records)

        handle = StringIO()
        NexusWriter(handle).write_alignment(a)
        handle.seek(0)
        data = handle.read()
        self.assertIn("ATGCTGCTGA" * 90, data)

        # Default causes interleave (columns > 1000)
        records = [
            SeqRecord(
                Seq("ATGCTGCTGA" * 110), id=_id, annotations={"molecule_type": "DNA"}
            )
            for _id in ["foo", "bar", "baz"]
        ]
        a = MultipleSeqAlignment(records)
        handle = StringIO()
        NexusWriter(handle).write_alignment(a)
        handle.seek(0)
        data = handle.read()
        self.assertNotIn("ATGCTGCTGA" * 90, data)
        self.assertIn("ATGCTGCTGA" * 7, data)

        # Override interleave: True
        records = [
            SeqRecord(
                Seq("ATGCTGCTGA" * 9), id=_id, annotations={"molecule_type": "DNA"}
            )
            for _id in ["foo", "bar", "baz"]
        ]
        a = MultipleSeqAlignment(records)
        handle = StringIO()
        NexusWriter(handle).write_alignment(a, interleave=True)
        handle.seek(0)
        data = handle.read()
        self.assertNotIn("ATGCTGCTGA" * 9, data)
        self.assertIn("ATGCTGCTGA" * 7, data)

        # Override interleave: False
        records = [
            SeqRecord(
                Seq("ATGCTGCTGA" * 110), id=_id, annotations={"molecule_type": "DNA"}
            )
            for _id in ["foo", "bar", "baz"]
        ]
        a = MultipleSeqAlignment(records)
        handle = StringIO()
        NexusWriter(handle).write_alignment(a, interleave=False)
        handle.seek(0)
        data = handle.read()
        self.assertIn("ATGCTGCTGA" * 110, data)

    def test_TreeTest1(self):
        """Test Tree module."""
        n = Nexus.Nexus(self.handle)
        t3 = n.trees[2]
        t2 = n.trees[2]
        t3.root_with_outgroup(["t1", "t5"])
        self.assertEqual(
            str(t3),
            "tree tree1 = (((((('one should be punished, for (that)!','isn''that [a] strange name?'),'t2 the name'),t8,t9),t6),t7),(t5,t1));",
        )
        self.assertEqual(t3.is_monophyletic(["t8", "t9", "t6", "t7"]), -1)
        self.assertEqual(t3.is_monophyletic(["t1", "t5"]), 13)
        t3.split(parent_id=t3.search_taxon("t9"))
        stdout = sys.stdout
        try:
            sys.stdout = StringIO()
            t3.display()
            output = sys.stdout.getvalue()
        finally:
            sys.stdout = stdout
        expected = """\
  #                            taxon            prev            succ    brlen blen (sum)  support              comment
  1    'isn''that [a] strange name?'               2              []   100.00     119.84    10.00                    -
  2                                -               4          [3, 1]     0.40      19.84     0.30                    -
  3 'one should be punished, for (that)!'               2              []     0.50      20.34        -                    -
  4                                -               6          [2, 5]     4.00      19.44     3.00                    -
  5                    't2 the name'               4              []     0.30      19.74        -                    -
  6                                -               9       [4, 7, 8]     2.00      15.44     1.00                    -
  7                               t8               6              []     1.20      16.64        -                    -
  8                               t9               6        [17, 18]     3.40      18.84        -                    -
  9                                -              11         [6, 10]     0.44      13.44    33.00                    -
 10                               t6               9              []     1.00      14.44        -                    -
 11                                -              16         [9, 12]    13.00      13.00    12.00                    -
 12                               t7              11              []    99.90     112.90        -                    -
 13                                -              16        [14, 15]     0.00       0.00     0.00                    -
 14                               t5              13              []    99.00      99.00        -                    -
 15                               t1              13              []     0.98       0.98        -                    -
 16                                -            None        [11, 13]     0.00       0.00        -                    -
 17                              t90               8              []     1.00      19.84        -                    -
 18                              t91               8              []     1.00      19.84        -                    -

Root:  16
"""
        self.assertEqual(len(output.split("\n")), len(expected.split("\n")))
        for l1, l2 in zip(output.split("\n"), expected.split("\n")):
            self.assertEqual(l1, l2)
        self.assertEqual(output, expected)
        self.assertEqual(t3.is_compatible(t2, threshold=0.3), [])

    def test_TreeTest2(self):
        """Handle text labels on internal nodes."""
        ts1b = (
            "(Cephalotaxus:125.000000,(Taxus:100.000000,Torreya:100.000000)"
            "TT1:25.000000)Taxaceae:90.000000;"
        )
        tree = Trees.Tree(ts1b)
        self.assertEqual(
            self._get_flat_nodes(tree),
            [
                ("Taxaceae", 90.0, None, None),
                ("Cephalotaxus", 125.0, None, None),
                ("TT1", 25.0, None, None),
                ("Taxus", 100.0, None, None),
                ("Torreya", 100.0, None, None),
            ],
        )
        tree.prune("Torreya")
        self.assertEqual(tree.all_ids(), [0, 1, 3])
        ts1c = (
            "(Cephalotaxus:125.000000,(Taxus:100.000000,Torreya:100.000000)"
            "25.000000)90.000000;"
        )
        tree = Trees.Tree(ts1c)
        self.assertEqual(
            self._get_flat_nodes(tree),
            [
                (None, 90.0, None, None),
                ("Cephalotaxus", 125.0, None, None),
                (None, 25.0, None, None),
                ("Taxus", 100.0, None, None),
                ("Torreya", 100.0, None, None),
            ],
        )
        self.assertFalse(tree.has_support())
        with self.assertRaises(Exception) as context:
            tree.randomize()
        self.assertIn(
            "Either number of taxa or list of taxa must be specified.",
            str(context.exception),
        )
        tree_rand = Trees.Tree(ts1c)
        tree_rand.randomize(ntax=4)
        self.assertEqual(
            sorted(tree_rand.get_taxa()), ["taxon1", "taxon2", "taxon3", "taxon4"]
        )
        tree.branchlength2support()
        tree.convert_absolute_support(2)
        self.assertEqual(
            self._get_flat_nodes(tree),
            [
                (None, 0.0, 90.0, None),
                ("Cephalotaxus", 0.0, 62.5, None),
                (None, 0.0, 12.5, None),
                ("Taxus", 0.0, 50.0, None),
                ("Torreya", 0.0, 50.0, None),
            ],
        )

        ts2 = (
            "(((t9:0.385832, (t8:0.445135,t4:0.41401)C:0.024032)B:0.041436,"
            "t6:0.392496)A:0.0291131, t2:0.497673, ((t0:0.301171,"
            "t7:0.482152)E:0.0268148, ((t5:0.0984167,t3:0.488578)G:0.0349662,"
            "t1:0.130208)F:0.0318288)D:0.0273876);"
        )
        tree = Trees.Tree(ts2)
        tree.branchlength2support()
        supports = []
        for i in tree.all_ids():
            node = tree.node(i)
            data = node.get_data()
            supports.append(data.support)
        self.assertEqual(
            supports,
            [
                0.0,
                0.0291131,
                0.041436,
                0.385832,
                0.024032,
                0.445135,
                0.41401,
                0.392496,
                0.497673,
                0.0273876,
                0.0268148,
                0.301171,
                0.482152,
                0.0318288,
                0.0349662,
                0.0984167,
                0.488578,
                0.130208,
            ],
        )
        ts3 = (
            "(((B 9:0.385832, (C 8:0.445135, C4:0.41401)C:0.024032)B:0.041436,"
            "A 6:0.392496)A:0.0291131, t2:0.497673, ((E 0:0.301171,"
            "E 7:0.482152)E:0.0268148, ((G 5:0.0984167,G 3:0.488578)G:0.0349662,"
            "F 1:0.130208)F:0.0318288)D:0.0273876);"
        )
        self.assertFalse(tree.is_identical(Trees.Tree(ts3)))
        tree = Trees.Tree(ts3)
        self.assertTrue(tree.is_bifurcating())
        self.assertTrue(tree.is_bifurcating(1))
        self.assertEqual(
            [tree.distance(0, n) for n in tree.all_ids()],
            [
                0.0,
                0.0291131,
                0.0705491,
                0.4563811,
                0.0945811,
                0.5397161,
                0.5085911,
                0.4216091,
                0.497673,
                0.0273876,
                0.0542024,
                0.3553734,
                0.5363544,
                0.0592164,
                0.09418259999999999,
                0.1925993,
                0.5827606,
                0.1894244,
            ],
        )

        subtree = tree.set_subtree(10)
        self.assertEqual(sorted(subtree), ["E 0", "E 7"])
        tree.collapse_genera()
        self.assertEqual(tree.all_ids(), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 17])

    def test_merge_with_support(self):
        """Test merge_with_support and consensus method."""
        ts1 = (
            "(((B 9:0.385832, (C 8:0.445135, C 4:0.41401)C:0.024032)B:0.041436,"
            "A 6:0.392496)A:0.0291131, t2:0.497673, ((E 0:0.301171,"
            "E 7:0.482152)E:0.0268148, ((G 5:0.0984167,G 3:0.488578)G:0.0349662,"
            "F 1:0.130208)F:0.0318288)D:0.0273876);"
        )

        tbs1 = (
            "(((B 9:0.385832, (C 8:0.445135, C 4:0.41401)C:0.024032)B:0.041436,"
            "A 6:0.392496)A:0.0291131, t2:0.497673, ((G 5:0.0984167,"
            "G 3:0.488578)E:0.0268148, ((E 0:0.301171, E 7:0.482152)G:0.0349662,"
            "F 1:0.130208)F:0.0318288)D:0.0273876);"
        )

        tbs2 = (
            "(((B 9:0.385832,A 6:0.392496 C:0.024032)B:0.041436, (C 8:0.445135,"
            "C 4:0.41401))A:0.0291131, t2:0.497673, ((E 0:0.301171, E 7:0.482152)"
            "E:0.0268148, ((G 5:0.0984167,G 3:0.488578)G:0.0349662,F 1:0.130208)"
            "F:0.0318288)D:0.0273876);"
        )

        t1 = Trees.Tree(ts1)
        tb1 = Trees.Tree(tbs1)
        tb2 = Trees.Tree(tbs2)

        t1.branchlength2support()
        tb1.branchlength2support()
        tb2.branchlength2support()

        t1.merge_with_support(bstrees=[tb1, tb2], threshold=0.2)

        supports = []
        for i in t1.all_ids():
            node = t1.node(i)
            data = node.get_data()
            supports.append(data.support)
        self.assertTrue(
            supports,
            [
                0.0,
                1.0,
                0.04,
                1.0,
                0.5,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                0.5,
                1.0,
                1.0,
                0.5,
                1.0,
                1.0,
                1.0,
                1.0,
            ],
        )

    def test_to_string(self):
        """Test to_string method."""
        t = Trees.Tree(
            "(((B 9:0.385832, (C 8:0.445135, C 4:0.41401)C:0.024032)B:0.041436,"
            "A 6:0.392496)A:0.0291131, t2:0.497673, ((E 0:0.301171,"
            "E 7:0.482152)E:0.0268148, ((G 5:0.0984167,G 3:0.488578)G:0.0349662,"
            "F 1:0.130208)F:0.0318288)D:0.0273876);"
        )
        tree_a_val = "((A 6,(B 9,(C 8,C 4))),t2,((E 0,E 7),(F 1,(G 5,G 3))));"
        self.assertEqual(t.to_string(ladderize="LEFT"), "tree a_tree = " + tree_a_val)

    def test_two_string_taxonomies(self):
        """Two non-numeric values on one node raise TreeError, not assert."""
        with self.assertRaisesRegex(
            Trees.TreeError, "Two string taxonomies 'x' and 'y'"
        ):
            Trees.Tree("((A,B)x:y:0.1,C);")

    def test_large_newick(self):
        with open(
            os.path.join(self.testfile_dir, "int_node_labels.nwk")
        ) as large_ex_handle:
            tree = Trees.Tree(large_ex_handle.read())

    def _get_flat_nodes(self, tree):
        cur_nodes = [tree.node(tree.root)]
        nodedata = []
        while len(cur_nodes) > 0:
            new_nodes = []
            for cur_node in cur_nodes:
                nodedata.append(
                    (
                        cur_node.data.taxon,
                        cur_node.data.branchlength,
                        cur_node.data.support,
                        cur_node.data.comment,
                    )
                )
                new_nodes.extend([tree.node(nid) for nid in cur_node.get_succ()])
            cur_nodes = new_nodes
        return nodedata

    def test_NexusComments(self):
        """Test the ability to parse nexus comments at internal and leaf nodes."""
        # A tree with simple comments throughout the tree.
        ts1b = "((12:0.13,19[&comment1]:0.13)[&comment2]:0.1,(20:0.171,11:0.171):0.13)[&comment3];"
        tree = Trees.Tree(ts1b)
        self.assertEqual(
            self._get_flat_nodes(tree),
            [
                (None, 0.0, None, "[&comment3]"),
                (None, 0.1, None, "[&comment2]"),
                (None, 0.13, None, None),
                ("12", 0.13, None, None),
                ("19", 0.13, None, "[&comment1]"),
                ("20", 0.171, None, None),
                ("11", 0.171, None, None),
            ],
        )

        # A tree with more complex comments throughout the tree.
        # This is typical of the MCC trees produced by `treeannotator` in the beast-mcmc suite of phylogenetic tools
        # The key difference being tested here is the ability to parse internal node comments that include ','.
        ts1b = "(((9[&rate_range={1.3E-5,0.10958320752991428},height_95%_HPD={0.309132419999969,0.3091324199999691},length_range={3.513906814545109E-4,0.4381986285528381},height_median=0.309132419999969,length_95%_HPD={0.003011577063374571,0.08041621647998398}]:0.055354097721950546,5[&rate_range={1.3E-5,0.10958320752991428},height_95%_HPD={0.309132419999969,0.3091324199999691},length_range={3.865051168833178E-5,0.4391594442572986},height_median=0.309132419999969,length_95%_HPD={0.003011577063374571,0.08041621647998398}]:0.055354097721950546)[&height_95%_HPD={0.3110921040545068,0.38690865205576275},length_range={0.09675588357303178,0.4332959544380489},length_95%_HPD={0.16680375169879613,0.36500804261814374}]:0.20039426358269385)[&height_95%_HPD={0.5289500597932948,0.6973881165460601},length_range={0.02586430194846201,0.29509451958008265},length_95%_HPD={0.0840287249314221,0.2411078625957056}]:0.23042678598484334)[&height_95%_HPD={0.7527502510685965,0.821862094763501},height_median=0.8014438411766163,height=0.795965080422763,posterior=1.0,height_range={0.49863013698599995,0.821862094763501},length=0.0];"
        tree = Trees.Tree(ts1b)
        self.assertEqual(
            self._get_flat_nodes(tree),
            [
                (
                    None,
                    0.0,
                    None,
                    "[&height_95%_HPD={0.7527502510685965,0.821862094763501},height_median=0.8014438411766163,height=0.795965080422763,posterior=1.0,height_range={0.49863013698599995,0.821862094763501},length=0.0]",
                ),
                (
                    None,
                    0.23042678598484334,
                    None,
                    "[&height_95%_HPD={0.5289500597932948,0.6973881165460601},length_range={0.02586430194846201,0.29509451958008265},length_95%_HPD={0.0840287249314221,0.2411078625957056}]",
                ),
                (
                    None,
                    0.20039426358269385,
                    None,
                    "[&height_95%_HPD={0.3110921040545068,0.38690865205576275},length_range={0.09675588357303178,0.4332959544380489},length_95%_HPD={0.16680375169879613,0.36500804261814374}]",
                ),
                (
                    "9",
                    0.055354097721950546,
                    None,
                    "[&rate_range={1.3E-5,0.10958320752991428},height_95%_HPD={0.309132419999969,0.3091324199999691},length_range={3.513906814545109E-4,0.4381986285528381},height_median=0.309132419999969,length_95%_HPD={0.003011577063374571,0.08041621647998398}]",
                ),
                (
                    "5",
                    0.055354097721950546,
                    None,
                    "[&rate_range={1.3E-5,0.10958320752991428},height_95%_HPD={0.309132419999969,0.3091324199999691},length_range={3.865051168833178E-5,0.4391594442572986},height_median=0.309132419999969,length_95%_HPD={0.003011577063374571,0.08041621647998398}]",
                ),
            ],
        )

    def test_NexusComments_two_per_node(self):
        """Strip both comments of a node and keep the first (issue 3456)."""
        # MrBayes writes one comment before the branch length (the clade
        # posterior) and one after it (the branch length statistics).
        tree = Trees.Tree(
            "((A[&p=1]:0.1[&l=2],B[&p=1]:0.2[&l=3])[&p=0.9]:0.3[&l=1],C:0.5);"
        )
        self.assertEqual(
            self._get_flat_nodes(tree),
            [
                (None, 0.0, None, None),
                (None, 0.3, None, "[&p=0.9]"),
                ("C", 0.5, None, None),
                ("A", 0.1, None, "[&p=1]"),
                ("B", 0.2, None, "[&p=1]"),
            ],
        )

    def test_NexusComments_labelled_node(self):
        """Keep an internal node's label as its taxon when it has comments."""
        for text in ("((A,B)x[&p=1]:0.1[&l=2],C);", "((A,B)x[&p=1]:0.1,C);"):
            with self.subTest(text=text):
                tree = Trees.Tree(text)
                data = tree.node(tree.node(tree.root).succ[0]).data
                self.assertEqual(data.taxon, "x")
                self.assertIsNone(data.support)
                self.assertEqual(data.branchlength, 0.1)
                self.assertEqual(data.comment, "[&p=1]")
                # a string support made this raise ValueError
                tree.to_string(plain=False)

    def test_NexusComments_to_string(self):
        """Write node comments when ignore_comments is False."""
        tree = Trees.Tree("((A[&rate=1]:0.1,B:0.2)[&rate=2]:0.3,C:0.4);")
        self.assertNotIn("[&rate", tree.to_string(plain=False))
        text = tree.to_string(plain=False, plain_newick=True, ignore_comments=False)
        # the old code read a missing attribute and dropped every comment
        self.assertEqual(
            text, "((A[&rate=1]:0.10000,B:0.20000)[&rate=2]0.00000:0.30000,C:0.40000)"
        )
        again = Trees.Tree(text)
        for old, new in zip(tree.chain.values(), again.chain.values()):
            self.assertEqual(new.data.taxon, old.data.taxon)
            self.assertEqual(new.data.comment, old.data.comment)

    def test_NexusComments_to_string_nodecomment(self):
        """A custom data class with nodecomment, not comment, still writes it."""

        class Data:
            def __init__(
                self, taxon=None, branchlength=0.0, support=None, comment=None
            ):
                self.taxon = taxon
                self.branchlength = branchlength
                self.support = support
                self.nodecomment = None

        tree = Trees.Tree("(A:0.1,B:0.2);", data=Data)
        tree.node(tree.search_taxon("A")).data.nodecomment = "[&rate=1]"
        text = tree.to_string(plain=False, plain_newick=True, ignore_comments=False)
        self.assertEqual(text, "(A[&rate=1]:0.10000,B:0.20000)")

    def test_NexusComments_to_string_ignore_comments(self):
        """With ignore_comments, to_string does not read the comments."""

        class Data(Trees.NodeData):
            @property
            def nodecomment(self):
                raise RuntimeError("nodecomment read")

        tree = Trees.Tree("(A:0.1,B:0.2);", data=Data)
        text = tree.to_string(plain=False, plain_newick=True)
        self.assertEqual(text, "(A:0.10000,B:0.20000)")

    def test_NexusComments_second_unmatched(self):
        """An unclosed second comment raises TreeError."""
        with self.assertRaisesRegex(Trees.TreeError, r"Found \[& without matching \]"):
            Trees.Tree("(A[&p=1]:0.1[&l=2,B);")

    def test_NexusComments_bracket_before_comment(self):
        """Find a comment's ] after its [&, past a ] in a quoted label."""
        # Searching for ] from the start of the text finds the label's ]
        # first. The old code then lost the comment; the loop that removes
        # every comment would never end.
        tree = Trees.Tree("((A,B)'x]y'[&p]:0.1,C);")
        data = tree.node(tree.node(tree.root).succ[0]).data
        self.assertEqual(data.comment, "[&p]")
        self.assertEqual(data.branchlength, 0.1)

    def test_empty_tree_string(self):
        """An empty tree, with or without ';' or spaces, has only its root."""
        for text in ("", ";", "   ", " ; "):
            with self.subTest(text=text):
                tree = Trees.Tree(text)
                self.assertEqual(len(tree.chain), 1)
                self.assertEqual(tree.node(tree.root).succ, [])

    def test_NexusComments_parentheses(self):
        """Accept ( and ) inside a comment on a leaf, a clade or the root."""
        A = ("A", 0.1, None, "[&x=(1)]")
        root = (None, 0.0, None, None)
        for text, nodes in [
            ("(A:0.1[&x=(1)],B);", [root, A, ("B", 0.0, None, None)]),
            ("(A[&x=(1)]:0.1,B);", [root, A, ("B", 0.0, None, None)]),
            (
                "((A,B)[&x=(1)]:0.1,C);",
                [
                    root,
                    (None, 0.1, None, "[&x=(1)]"),
                    ("C", 0.0, None, None),
                    ("A", 0.0, None, None),
                    ("B", 0.0, None, None),
                ],
            ),
            (
                "(A,B)[&x=(1)];",
                [
                    (None, 0.0, None, "[&x=(1)]"),
                    ("A", 0.0, None, None),
                    ("B", 0.0, None, None),
                ],
            ),
            # a tree of one leaf is read as "(A)"; it raised IndexError
            ("A:0.1[&x=(1)];", [root, A]),
            ("A;", [root, ("A", 0.0, None, None)]),
        ]:
            with self.subTest(text=text):
                self.assertEqual(self._get_flat_nodes(Trees.Tree(text)), nodes)
        with self.assertRaisesRegex(Trees.TreeError, "Parentheses do not match"):
            Trees.Tree("(A[&x=(1)](B),C);")

    def test_NexusComments_many_per_node(self):
        """Strip a million comments from one node in linear time.

        Tree._get_values rebuilt the node's text once for each comment it
        removed, so this 4 MB tree took minutes to parse. A stray [ before
        a comment, which that rebuilding joined to the text after the
        comment to make a new one, is now an error.
        """
        code = r"""
import json
import sys
from Bio.Nexus import Trees
tree = Trees.Tree("(A:1" + "[&x]" * int(sys.argv[1]) + ",B:2);")
data = tree.node(tree.node(tree.root).succ[0]).data
print(json.dumps([data.taxon, data.branchlength, data.comment]))
"""
        self.assertEqual(parse_in_child(self, code, 1_000_000), ["A", 1.0, "[&x]"])
        with self.assertRaisesRegex(Trees.TreeError, r"\[ before a comment"):
            Trees.Tree("(A:1[[&a]&b],B:2);")


class TestSelf(unittest.TestCase):
    def test_repeated_names_no_taxa(self):
        handle = StringIO(
            """#NEXUS
        [TITLE: NoName]
        begin data;
        dimensions ntax=4 nchar=50;
        format interleave datatype=protein   gap=- symbols="FSTNKEYVQMCLAWPHDRIG";
        matrix
        CYS1_DICDI          -----MKVIL LFVLAVFTVF VSS------- --------RG IPPEEQ----
        ALEU_HORVU          MAHARVLLLA LAVLATAAVA VASSSSFADS NPIRPVTDRA ASTLESAVLG
        CATH_HUMAN          ------MWAT LPLLCAGAWL LGV------- -PVCGAAELS VNSLEK----
        CYS1_DICDI          -----MKVIL LFVLAVFTVF VSS------- --------RG IPPEEQ---X
        ;
        end;
        """
        )
        alignments = AlignmentIterator(handle)
        alignment = next(alignments)
        self.assertEqual(alignment.shape, (4, 50))
        self.assertEqual(
            str(alignment),
            """\
CYS1_DICD         0 -----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ---- 26
ALEU_HORV         0 MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLG 50
CATH_HUMA         0 ------MWATLPLLCAGAWLLGV--------PVCGAAELSVNSLEK---- 32
CYS1_DICD         0 -----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ---X 27
""",
        )
        self.assertEqual(
            alignment[0], "-----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ----"
        )
        self.assertEqual(alignment.sequences[0].id, "CYS1_DICDI")
        self.assertEqual(
            alignment[1], "MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLG"
        )
        self.assertEqual(alignment.sequences[1].id, "ALEU_HORVU")
        self.assertEqual(
            alignment[2], "------MWATLPLLCAGAWLLGV--------PVCGAAELSVNSLEK----"
        )
        self.assertEqual(alignment.sequences[2].id, "CATH_HUMAN")
        self.assertEqual(
            alignment[3], "-----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ---X"
        )
        self.assertEqual(alignment.sequences[3].id, "CYS1_DICDI")
        self.assertRaises(StopIteration, next, alignments)

    def test_repeated_names_no_taxa_msa(self):
        handle = StringIO(
            """#NEXUS
        [TITLE: NoName]
        begin data;
        dimensions ntax=4 nchar=50;
        format interleave datatype=protein   gap=- symbols="FSTNKEYVQMCLAWPHDRIG";
        matrix
        CYS1_DICDI          -----MKVIL LFVLAVFTVF VSS------- --------RG IPPEEQ----
        ALEU_HORVU          MAHARVLLLA LAVLATAAVA VASSSSFADS NPIRPVTDRA ASTLESAVLG
        CATH_HUMAN          ------MWAT LPLLCAGAWL LGV------- -PVCGAAELS VNSLEK----
        CYS1_DICDI          -----MKVIL LFVLAVFTVF VSS------- --------RG IPPEEQ---X
        ;
        end;
        """
        )
        alignments = NexusIterator(handle)
        alignment = next(alignments)
        self.assertEqual(len(alignment), 4)
        self.assertEqual(
            str(alignment),
            """\
Alignment with 4 rows and 50 columns
-----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ---- CYS1_DICDI
MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLG ALEU_HORVU
------MWATLPLLCAGAWLLGV--------PVCGAAELSVNSLEK---- CATH_HUMAN
-----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ---X CYS1_DICDI.copy""",
        )
        records = iter(alignment)
        record = next(records)
        self.assertEqual(
            record.seq, Seq("-----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ----")
        )
        self.assertEqual(record.name, "CYS1_DICDI")
        self.assertEqual(record.id, "CYS1_DICDI")
        record = next(records)
        self.assertEqual(
            record.seq, Seq("MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLG")
        )
        self.assertEqual(record.name, "ALEU_HORVU")
        self.assertEqual(record.id, "ALEU_HORVU")
        record = next(records)
        self.assertEqual(
            record.seq, Seq("------MWATLPLLCAGAWLLGV--------PVCGAAELSVNSLEK----")
        )
        self.assertEqual(record.name, "CATH_HUMAN")
        self.assertEqual(record.id, "CATH_HUMAN")
        record = next(records)
        self.assertEqual(
            record.seq, Seq("-----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ---X")
        )
        self.assertEqual(record.name, "CYS1_DICDI")
        self.assertEqual(record.id, "CYS1_DICDI.copy")
        self.assertRaises(StopIteration, next, records)
        self.assertRaises(StopIteration, next, alignments)

    def test_repeated_names_with_taxa(self):
        handle = StringIO(
            """#NEXUS
        [TITLE: NoName]
        begin taxa
        CYS1_DICDI
        ALEU_HORVU
        CATH_HUMAN
        CYS1_DICDI;
        end;
        begin data;
        dimensions ntax=4 nchar=50;
        format interleave datatype=protein   gap=- symbols="FSTNKEYVQMCLAWPHDRIG";
        matrix
        CYS1_DICDI          -----MKVIL LFVLAVFTVF VSS------- --------RG IPPEEQ----
        ALEU_HORVU          MAHARVLLLA LAVLATAAVA VASSSSFADS NPIRPVTDRA ASTLESAVLG
        CATH_HUMAN          ------MWAT LPLLCAGAWL LGV------- -PVCGAAELS VNSLEK----
        CYS1_DICDI          -----MKVIL LFVLAVFTVF VSS------- --------RG IPPEEQ---X
        ;
        end;
        """
        )
        alignments = AlignmentIterator(handle)
        alignment = next(alignments)
        self.assertEqual(len(alignment), 4)
        self.assertEqual(
            str(alignment),
            """\
CYS1_DICD         0 -----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ---- 26
ALEU_HORV         0 MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLG 50
CATH_HUMA         0 ------MWATLPLLCAGAWLLGV--------PVCGAAELSVNSLEK---- 32
CYS1_DICD         0 -----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ---X 27
""",
        )
        self.assertEqual(
            alignment[0], "-----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ----"
        )
        self.assertEqual(alignment.sequences[0].id, "CYS1_DICDI")
        self.assertEqual(
            alignment[1], "MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLG"
        )
        self.assertEqual(alignment.sequences[1].id, "ALEU_HORVU")
        self.assertEqual(
            alignment[2], "------MWATLPLLCAGAWLLGV--------PVCGAAELSVNSLEK----"
        )
        self.assertEqual(alignment.sequences[2].id, "CATH_HUMAN")
        self.assertEqual(
            alignment[3], "-----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ---X"
        )
        self.assertEqual(alignment.sequences[3].id, "CYS1_DICDI")
        self.assertRaises(StopIteration, next, alignments)

    def test_repeated_names_with_taxa_msa(self):
        handle = StringIO(
            """#NEXUS
        [TITLE: NoName]
        begin taxa
        CYS1_DICDI
        ALEU_HORVU
        CATH_HUMAN
        CYS1_DICDI;
        end;
        begin data;
        dimensions ntax=4 nchar=50;
        format interleave datatype=protein   gap=- symbols="FSTNKEYVQMCLAWPHDRIG";
        matrix
        CYS1_DICDI          -----MKVIL LFVLAVFTVF VSS------- --------RG IPPEEQ----
        ALEU_HORVU          MAHARVLLLA LAVLATAAVA VASSSSFADS NPIRPVTDRA ASTLESAVLG
        CATH_HUMAN          ------MWAT LPLLCAGAWL LGV------- -PVCGAAELS VNSLEK----
        CYS1_DICDI          -----MKVIL LFVLAVFTVF VSS------- --------RG IPPEEQ---X
        ;
        end;
        """
        )
        alignments = NexusIterator(handle)
        alignment = next(alignments)
        self.assertEqual(len(alignment), 4)
        self.assertEqual(
            str(alignment),
            """\
Alignment with 4 rows and 50 columns
-----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ---- CYS1_DICDI
MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLG ALEU_HORVU
------MWATLPLLCAGAWLLGV--------PVCGAAELSVNSLEK---- CATH_HUMAN
-----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ---X CYS1_DICDI.copy""",
        )
        records = iter(alignment)
        record = next(records)
        self.assertEqual(
            record.seq, Seq("-----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ----")
        )
        self.assertEqual(record.name, "CYS1_DICDI")
        self.assertEqual(record.id, "CYS1_DICDI")
        record = next(records)
        self.assertEqual(
            record.seq, Seq("MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTLESAVLG")
        )
        self.assertEqual(record.name, "ALEU_HORVU")
        self.assertEqual(record.id, "ALEU_HORVU")
        record = next(records)
        self.assertEqual(
            record.seq, Seq("------MWATLPLLCAGAWLLGV--------PVCGAAELSVNSLEK----")
        )
        self.assertEqual(record.name, "CATH_HUMAN")
        self.assertEqual(record.id, "CATH_HUMAN")
        record = next(records)
        self.assertEqual(
            record.seq, Seq("-----MKVILLFVLAVFTVFVSS---------------RGIPPEEQ---X")
        )
        self.assertEqual(record.name, "CYS1_DICDI")
        self.assertEqual(record.id, "CYS1_DICDI.copy")
        self.assertRaises(StopIteration, next, records)
        self.assertRaises(StopIteration, next, alignments)

    def test_empty_file_read(self):
        with self.assertRaises(ValueError) as cm:
            AlignmentIterator(StringIO())
        self.assertEqual(str(cm.exception), "Empty file.")

    def test_empty_file_read_msa(self):
        self.assertEqual([], list(NexusIterator(StringIO())))

    def test_multiple_output(self):
        records = [
            SeqRecord(
                Seq("ATGCTGCTGAT"), id="foo", annotations={"molecule_type": "DNA"}
            ),
            SeqRecord(
                Seq("ATGCTGCAGAT"), id="bar", annotations={"molecule_type": "DNA"}
            ),
            SeqRecord(
                Seq("ATGCTGCGGAT"), id="baz", annotations={"molecule_type": "DNA"}
            ),
        ]
        a = Alignment(records)

        handle = StringIO()
        AlignmentWriter(handle).write([a])
        handle.seek(0)
        data = handle.read()
        self.assertEqual(
            data,
            """\
#NEXUS
begin data;
dimensions ntax=3 nchar=11;
format datatype=dna missing=? gap=-;
matrix
foo ATGCTGCTGAT
bar ATGCTGCAGAT
baz ATGCTGCGGAT
;
end;
""",
        )
        handle = StringIO()
        with self.assertRaises(ValueError):
            AlignmentWriter(handle).write([a, a])

    def test_multiple_output_msa(self):
        records = [
            SeqRecord(
                Seq("ATGCTGCTGAT"), id="foo", annotations={"molecule_type": "DNA"}
            ),
            SeqRecord(
                Seq("ATGCTGCAGAT"), id="bar", annotations={"molecule_type": "DNA"}
            ),
            SeqRecord(
                Seq("ATGCTGCGGAT"), id="baz", annotations={"molecule_type": "DNA"}
            ),
        ]
        a = MultipleSeqAlignment(records)

        handle = StringIO()
        NexusWriter(handle).write_file([a])
        handle.seek(0)
        data = handle.read()
        self.assertTrue(data.startswith("#NEXUS\nbegin data;\n"), data)
        self.assertTrue(data.endswith("end;\n"), data)

        handle = StringIO()
        with self.assertRaises(ValueError):
            NexusWriter(handle).write_file([a, a])

    def test_translate_dashed_taxon_names(self):
        # Unquoted taxon labels containing a hyphen must be preserved in the
        # TRANSLATE block instead of raising "Missing ','"
        # (biopython/biopython#1022).
        handle = StringIO(
            """#NEXUS
        begin trees;
        translate
        1 dog-cat,
        2 mouse-rat,
        3 foo-bar
        ;
        tree t1 = [&U] (1,(2,3));
        end;
        """
        )
        nexus = Nexus.Nexus(handle)
        self.assertEqual(nexus.translate, {1: "dog-cat", 2: "mouse-rat", 3: "foo-bar"})
        # Bio.Phylo, the path the issue reported, goes through the same parser.
        # NexusIO passes each name through Nexus.safename, which quotes it.
        handle.seek(0)
        tree = Phylo.read(handle, "nexus")
        self.assertEqual(
            [clade.name for clade in tree.get_terminals()],
            ["'dog-cat'", "'mouse-rat'", "'foo-bar'"],
        )

    def test_translate_format_error(self):
        # A bad or missing taxon number is a format error, but an unrelated
        # exception must not be disguised as one (BioAGI-org/biop-ai#552).
        nexus = Nexus.Nexus()
        for options in ("1 a, x b", ""):
            with self.assertRaisesRegex(Nexus.NexusError, "Format error"):
                nexus._translate(options)
        bug = AttributeError("bug in quotestrip")
        with mock.patch.object(Nexus, "quotestrip", side_effect=bug):
            with self.assertRaises(AttributeError) as cm:
                nexus._translate("1 a, 2 b")
        self.assertIs(cm.exception, bug)
        # an error from reading the number is not a format error either
        with self.assertRaisesRegex(Nexus.NexusError, "^Missing closing quote "):
            nexus._translate("'1")
        # A missing label was stored as None, and the tree using it then
        # failed with AttributeError.
        text = "#NEXUS\nbegin trees;\ntranslate 1;\ntree t = (1,1);\nend;\n"
        with self.assertRaisesRegex(Nexus.NexusError, "^Missing label for 1 "):
            Nexus.Nexus(text)

    def test_charpartition_quoted_subset_names_round_trip(self):
        # append_sets must quote charpartition subset names that contain
        # spaces or commas, so the written file reads back the same
        # (BioAGI-org/biop-ai#547).
        data = """#NEXUS
begin data;
dimensions ntax=2 nchar=6;
format datatype=dna missing=? gap=-;
matrix
a ACGTAC
b ACGTAC
;
end;
"""
        nexus = Nexus.Nexus(
            data
            + "begin sets;\n"
            + "charpartition part = 'gene one': 1-3, 'gene, two': 4-6;\n"
            + "end;\n"
        )
        expected = {"part": {"gene one": [0, 1, 2], "gene, two": [3, 4, 5]}}
        self.assertEqual(nexus.charpartitions, expected)
        again = Nexus.Nexus(data + nexus.append_sets())
        self.assertEqual(again.charpartitions, expected)


class TestStandardDatatype(unittest.TestCase):
    """Matrices with datatype=standard, such as morphological characters.

    See biopython/biopython#2814.
    """

    # Mesquite writes SYMBOLS with a space between (and before) the symbols.
    # The NEXUS standard says this is the same as SYMBOLS="012".
    template = """\
#NEXUS
begin taxa;
    dimensions ntax=3;
    taxlabels alpha beta gamma;
end;
begin characters;
    dimensions nchar=5;
    format datatype = standard gap = - missing = ? %s;
    matrix
    alpha 01201
    beta  %s
    gamma 22100
;
end;
"""

    def test_commandline_keeps_quoted_values_whole(self):
        """A quoted option value with spaces or "=" in it is one value."""
        line = Nexus.Commandline(
            "format datatype=standard symbols = \"0 1 2\" equate='R=(02) Y=(1)'",
            "characters",
        )
        self.assertEqual(
            line.options,
            {
                "datatype": "standard",
                "symbols": '"0 1 2"',
                "equate": "'R=(02) Y=(1)'",
            },
        )
        # A doubled quote stands for the quote itself, and text straight after
        # the closing quote stays in the value, as it did before
        line = Nexus.Commandline("format equate='it''s a' x='ab'cd", "characters")
        self.assertEqual(line.options, {"equate": "'it''s a'", "x": "'ab'cd"})

    def test_symbols_with_and_without_whitespace(self):
        """SYMBOLS = "  0 1 2" reads the same as SYMBOLS="012"."""
        for symbols in ('symbols = "  0 1 2"', 'symbols="012"', "symbols='0 1 2'"):
            with self.subTest(symbols=symbols):
                n = Nexus.Nexus(self.template % (symbols, "1?-21"))
                self.assertEqual("".join(n.symbols), "012")
                self.assertEqual(
                    {taxon: str(row) for taxon, row in n.matrix.items()},
                    {"alpha": "01201", "beta": "1?-21", "gamma": "22100"},
                )

    def test_phylo_reads_trees(self):
        """Bio.Phylo reads the trees from a file with such a matrix."""
        text = self.template % ('symbols = "  0 1 2"', "1?-21")
        text += "begin trees;\n    tree one = ((alpha,beta),gamma);\nend;\n"
        (tree,) = Phylo.parse(StringIO(text), "nexus")
        self.assertEqual(
            [clade.name for clade in tree.get_terminals()], ["alpha", "beta", "gamma"]
        )

    def test_read_msa(self):
        """Bio.AlignIO and Bio.SeqIO read a standard matrix as one letter per character."""
        text = self.template % ('symbols = "  0 1 2"', "1?-21")
        alignment = next(NexusIterator(StringIO(text)))
        self.assertEqual(
            [record.id for record in alignment], ["alpha", "beta", "gamma"]
        )
        self.assertEqual(
            [record.seq for record in alignment],
            [Seq("01201"), Seq("1?-21"), Seq("22100")],
        )
        records = SeqIO.parse(StringIO(text), "nexus")
        self.assertEqual(
            [str(record.seq) for record in records], ["01201", "1?-21", "22100"]
        )

    def test_read(self):
        """Bio.Align reads a standard matrix as one letter per character."""
        text = self.template % ('symbols = "  0 1 2"', "1?-21")
        alignment = next(AlignmentIterator(StringIO(text)))
        self.assertEqual(alignment.shape, (3, 5))
        self.assertEqual(
            [record.id for record in alignment.sequences], ["alpha", "beta", "gamma"]
        )
        self.assertEqual(list(alignment), ["01201", "1?-21", "22100"])

    def test_coding_of_one_state(self):
        """A polymorphic (1) or uncertain {1} coding of one state is that state."""
        text = self.template % ('symbols="012"', "1(1)-{2}1")
        alignment = next(NexusIterator(StringIO(text)))
        self.assertEqual(alignment[1].seq, Seq("11-21"))
        alignment = next(AlignmentIterator(StringIO(text)))
        self.assertEqual(alignment[1], "11-21")

    def test_polymorphic_or_uncertain_coding_raises(self):
        """A coding of several states cannot be one letter, so reading raises ValueError."""
        for row, message in (
            ("1(01)210", r"Taxon beta: character 2 has the polymorphic coding \(01\)"),
            ("12{02}11", r"Taxon beta: character 3 has the uncertain coding \{02\}"),
        ):
            text = self.template % ('symbols="012"', row)
            # Bio.Nexus itself keeps the coding
            n = Nexus.Nexus(text)
            self.assertEqual(str(n.matrix["beta"]), row)
            with self.subTest(row=row, reader="Bio.AlignIO"):
                with self.assertRaisesRegex(ValueError, message):
                    next(NexusIterator(StringIO(text)))
            with self.subTest(row=row, reader="Bio.Align"):
                with self.assertRaisesRegex(ValueError, message):
                    next(AlignmentIterator(StringIO(text)))


class TestCommandDispatch(unittest.TestCase):
    def test_error_in_known_command_is_not_unknown_command(self):
        """An AttributeError inside a known command's handler reaches the caller."""
        text = "#NEXUS\nbegin data;\n%s;\nend;\n"
        bug = AttributeError("bug in _title")
        with mock.patch.object(Nexus.Nexus, "_title", side_effect=bug):
            with self.assertRaises(AttributeError) as cm:
                Nexus.Nexus(text % "title x")
        self.assertIs(cm.exception, bug)
        with self.assertRaisesRegex(Nexus.NexusError, "^Unknown command: nosuch "):
            Nexus.Nexus(text % "nosuch")


class TestUnclosedQuote(unittest.TestCase):
    """A quoted word that runs to the end of the text raises NexusError.

    The file may close the quote, but not within the part of it that one
    parser reads, such as a line of the matrix.
    """

    def setUp(self):
        # Before the fix next_word kept reading the empty buffer, adding
        # None to the word until memory ran out. Stop such a loop early so
        # that the test fails instead.
        next_char = Nexus.CharBuffer.__next__

        def bounded_next(buffer):
            if not buffer.buffer:
                buffer.reads_at_end = getattr(buffer, "reads_at_end", 0) + 1
                if buffer.reads_at_end > 1000:
                    raise AssertionError("CharBuffer read past its end 1000 times")
            return next_char(buffer)

        patcher = mock.patch.object(Nexus.CharBuffer, "__next__", bounded_next)
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_next_word(self):
        """CharBuffer.next_word raises NexusError on a word left open."""
        for text, word in (("'t1", "'t1"), ('"a b', '"a b'), ("'it''", "'it'")):
            with self.subTest(text=text):
                with self.assertRaises(Nexus.NexusError) as cm:
                    Nexus.CharBuffer(text).next_word()
                self.assertEqual(str(cm.exception), f"Missing closing quote in: {word}")
        # Closed quotes, and a doubled quote inside them, read as before
        self.assertEqual(Nexus.CharBuffer("'it''s' x").next_word(), "'it's'")
        self.assertEqual(Nexus.CharBuffer('"a b" x').next_word(), '"a b"')

    def test_read(self):
        """Reading such a file raises NexusError."""
        matrix = (
            "begin data;\ndimensions ntax=2 nchar=4;\nformat datatype=dna;\nmatrix\n"
        )
        for text, word in (
            # a taxon label split over two lines of the matrix
            (matrix + "'t1\nt2' ACGT\nt3 ACGT\n;\nend;\n", "'t1"),
        ):
            with self.subTest(text=text):
                with self.assertRaisesRegex(
                    Nexus.NexusError, f"^Missing closing quote in: {word}$"
                ):
                    Nexus.Nexus("#NEXUS\n" + text)


class TestCharBufferWhitespace(unittest.TestCase):
    def test_peek_and_next_skip_the_same_whitespace(self):
        """next_nonwhitespace returns the character peek_nonwhitespace saw."""
        # Before the fix next_nonwhitespace skipped only " \t\n" and so
        # returned the \v, \f, no-break space or \r itself.
        for text in ("\v=", "\f(", "\xa0-", "\r*"):
            with self.subTest(text=text):
                buffer = Nexus.CharBuffer(text)
                self.assertEqual(buffer.peek_nonwhitespace(), text[-1])
                self.assertEqual(buffer.next_nonwhitespace(), text[-1])
                self.assertIsNone(buffer.next_nonwhitespace())
        self.assertEqual(
            Nexus._split_tree_command("\f* t1 = (a,b)"), ("t1", 1.0, False, "(a,b)")
        )


class TestPartitionQuotes(unittest.TestCase):
    """A comma inside a quoted name does not split a partition."""

    def test_comma_in_quoted_name(self):
        """Single and double quotes both keep a comma in a name."""
        for quote, other in (("'", '"'), ('"', "'")):
            with self.subTest(quote=quote):
                # A comma in a double-quoted name split the partition there,
                # raising "Missing closing quote in: \"a". The other quote
                # character inside a name must not close the quote.
                nexus = Nexus.Nexus(
                    "#NEXUS\nbegin data;\ndimensions ntax=4 nchar=6;\n"
                    "format datatype=dna;\nmatrix\n'a,b' ACGTAC\n"
                    f"{quote}a{other}b,c{quote} ACGTAC\nc ACGTAC\nd ACGTAC\n"
                    ";\nend;\nbegin sets;\ntaxpartition tp = "
                    f"one: {quote}a,b{quote}, two: {quote}a{other}b,c{quote} c, "
                    "three: d;\ncharpartition cp = "
                    f"{quote}x,y{quote}: 1-2, {quote}x{other}y,z{quote}: 3-4, "
                    "z: 5-6;\nend;\n"
                )
                self.assertEqual(
                    nexus.taxpartitions,
                    {
                        "tp": {
                            "one": ["a,b"],
                            "two": [f"a{other}b,c", "c"],
                            "three": ["d"],
                        }
                    },
                )
                self.assertEqual(
                    nexus.charpartitions,
                    {"cp": {"x,y": [0, 1], f"x{other}y,z": [2, 3], "z": [4, 5]}},
                )


class TestMatchchar(unittest.TestCase):
    """A matchchar in the matrix stands for the first taxon's character."""

    def setUp(self):
        # Before the fix _matrix made a new Seq for every matchchar it
        # replaced, and replaced a matchchar with itself forever when the
        # first taxon had the matchchar at that position. Stop such a loop
        # early so that the test fails instead.
        made = 0

        def bounded_seq(*args, **kwargs):
            nonlocal made
            made += 1
            if made > 1000:
                raise AssertionError("_matrix made more than 1000 Seq objects")
            return Seq(*args, **kwargs)

        patcher = mock.patch.object(Nexus, "Seq", bounded_seq)
        patcher.start()
        self.addCleanup(patcher.stop)

    def parse(self, matchchar, *rows, interleave=""):
        """Return the matrix of a DNA data block as a dict of strings."""
        nexus = Nexus.Nexus(
            "#NEXUS\nbegin data;\ndimensions ntax=2 nchar=4;\n"
            f"format datatype=dna matchchar={matchchar}{interleave};\n"
            "matrix\n" + "\n".join(rows) + "\n;\nend;\n"
        )
        return {taxon: str(seq) for taxon, seq in nexus.matrix.items()}

    def test_substitution(self):
        """Each matchchar is replaced by the first taxon's character."""
        self.assertEqual(
            self.parse(".", "a ACGT", "b ..-T"), {"a": "ACGT", "b": "AC-T"}
        )
        # in an interleaved matrix, by the first taxon's in the same block
        self.assertEqual(
            self.parse(".", "a AC", "b .-", "a GT", "b A.", interleave=" interleave"),
            {"a": "ACGT", "b": "A-AT"},
        )

    def test_matchchar_in_first_taxon(self):
        """A matchchar that is also a gap, missing or state symbol parses."""
        for matchchar, first, second, expected in (
            ("-", "ACG-", "AC--", "ACG-"),
            ("?", "ACG?", "ACG?", "ACG?"),
            ("A", "ACGT", "ACGA", "ACGT"),
        ):
            with self.subTest(matchchar=matchchar):
                self.assertEqual(
                    self.parse(matchchar, "a " + first, "b " + second),
                    {"a": first, "b": expected},
                )

    def test_matchchar_past_first_taxon(self):
        """A matchchar past the end of the first taxon raises NexusError."""
        with self.assertRaisesRegex(
            Nexus.NexusError,
            r"^Taxon b: matchchar \. past the end of the first taxon's sequence ",
        ):
            self.parse(".", "a ACGT", "b ACGTA.")

    def test_short_rows_after_long_reference(self):
        """Short rows after a long first row take time in their own length."""
        # The first taxon's sequence was decoded again for every later row
        # of the block, so 100,000 one-character rows under a 4,000,000
        # character one took minutes before the length check rejected them.
        code = r"""
import json
import sys
from Bio.Nexus import Nexus
n, length = int(sys.argv[1]), int(sys.argv[2])
try:
    Nexus.Nexus(
        f"#NEXUS\nbegin data;\ndimensions ntax={n} nchar={length};\n"
        "format datatype=dna matchchar=. interleave;\nmatrix\n"
        f"t0 {'A' * length}\n" + "".join(f"t{i} .\n" for i in range(1, n))
        + ";\nend;\n"
    )
except Nexus.NexusError as error:
    print(json.dumps(str(error)))
"""
        self.assertEqual(
            parse_in_child(self, code, 100_000, 4_000_000),
            "Matrix Nchar 4000000 does not match data length (1) for taxon t1",
        )


class TestIndexRange(unittest.TestCase):
    """Character and taxon numbers must lie between 1 and nchar or ntax."""

    def parse(self, commands, charlabels=""):
        """Return a Nexus object of one taxon and two characters."""
        return Nexus.Nexus(
            "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=2;\nformat datatype=dna;\n"
            f"{charlabels}matrix\nt1 AC\n;\nend;\nbegin sets;\n{commands}end;\n"
        )

    def test_valid_numbers(self):
        """Numbers from 1 to nchar or ntax read as before."""
        nexus = self.parse(
            "charset c = 1 - 2;\ncharset d = 2 1 - 2\\2;\n"
            "charpartition p = a: 1, b: 2;\ntaxset t = 1;\n",
            charlabels="charlabels 1 x, 2 y;\n",
        )
        self.assertEqual(nexus.charsets, {"c": [0, 1], "d": [0, 1]})
        self.assertEqual(nexus.charpartitions, {"p": {"a": [0], "b": [1]}})
        self.assertEqual(nexus.taxsets, {"t": ["t1"]})
        self.assertEqual(nexus.charlabels, {0: "x", 1: "y"})

    def test_character_number_out_of_range(self):
        """A character number below 1 or above nchar raises NexusError."""
        for commands, charlabels, number in (
            # a quoted negative number was accepted, so '-200000000' - 0 was
            # a list of 200 million numbers
            ("charset c = '-5' - 2;\n", "", -5),
            ("charset c = 0;\n", "", 0),
            ("charset c = 1 - 3;\n", "", 3),
            ("charpartition p = a: '-5' - 1, b: 2;\n", "", -5),
            ("codonposset p = N: '-5' - 1, 1: 2;\n", "", -5),
            ("", "charlabels '-3' x;\n", -3),
            ("", "charlabels 3 x;\n", 3),
        ):
            with self.subTest(commands=commands, charlabels=charlabels):
                with self.assertRaisesRegex(
                    Nexus.NexusError,
                    rf"^Illegal character identifier: {number}, "
                    r"not between 1 and nchar \(=2\)\.$",
                ):
                    self.parse(commands, charlabels)

    def test_taxon_number_out_of_range(self):
        """A taxon number below 1 or above ntax raises NexusError."""
        for commands, number in (
            ("taxset t = 0;\n", 0),
            ("taxset t = 2;\n", 2),
            ("taxpartition p = a: '-1';\n", -1),
        ):
            with self.subTest(commands=commands):
                with self.assertRaisesRegex(
                    Nexus.NexusError,
                    rf"^Illegal taxon identifier: {number}, "
                    r"not between 1 and ntax \(=1\)\.$",
                ):
                    self.parse(commands)

    def test_taxon_number_without_label(self):
        """A taxon number with no taxon label yet raises NexusError."""
        # Without a matrix, or before it, taxlabels is shorter than ntax, and
        # such a number raised IndexError from Nexus.Nexus().
        for nexus in (
            "begin taxa;\ndimensions ntax=2;\ntaxlabels a b;\nend;\n"
            "begin sets;\ntaxset t = 1;\nend;\n",
            "begin taxa;\ndimensions ntax=2;\nend;\n"
            "begin sets;\ntaxpartition p = a: 2;\nend;\n",
        ):
            with self.subTest(nexus=nexus):
                with self.assertRaisesRegex(
                    Nexus.NexusError, r"^Unknown taxon identifier: [12]$"
                ):
                    Nexus.Nexus("#NEXUS\n" + nexus)

    def test_only_format_errors_are_formatting_errors(self):
        """A badly formed list is a formatting error; MemoryError is not."""
        nexus = self.parse("charset c1 = 1;\n")
        for line in ("x = 1-", "x = 1-2\\0", "x = 1-2\\q", "x = c1-2"):
            with self.subTest(line=line):
                with self.assertRaisesRegex(Nexus.NexusError, "^Formatting error"):
                    nexus._charset(line)
        with self.assertRaisesRegex(Nexus.NexusError, "^Formatting error"):
            nexus._taxset("x = t1-")
        # Every exception but NexusError was reported as "Formatting error".
        with mock.patch.object(nexus, "_resolve", side_effect=MemoryError):
            with self.assertRaises(MemoryError):
                nexus._charset("x = 1-2")

    def test_standard_and_vector_qualifiers(self):
        """A (STANDARD) qualifier after a set name is accepted; VECTOR is not."""
        # (STANDARD) was rejected as "Formatting error in line: ...".
        nexus = self.parse(
            "charset c (STANDARD) = 1-2;\ntaxset t (Standard) = t1;\n"
            "charpartition p (standard) = a: 1, b: 2;\n"
        )
        self.assertEqual(nexus.charsets, {"c": [0, 1]})
        self.assertEqual(nexus.taxsets, {"t": ["t1"]})
        self.assertEqual(nexus.charpartitions, {"p": {"a": [0], "b": [1]}})
        for line, message in (
            ("c (VECTOR) = 11", "^Unsupported VECTOR format in line: c "),
            ("c (FOO) = 1", "^Unknown qualifier in line: c "),
            ("c (STANDARD = 1", "^Formatting error in line: c "),
        ):
            with self.subTest(line=line):
                with self.assertRaisesRegex(Nexus.NexusError, message):
                    nexus._charset(line)


class TestNcharBound(unittest.TestCase):
    """A stated nchar must fit in the input, as the matrix needs that many."""

    def test_nchar_larger_than_input(self):
        """An nchar that no matrix in the input could hold raises NexusError."""
        # with no matrix, nchar was only a number from the file, and a
        # charset range up to it was expanded in full: nchar=200000000 in a
        # 101-byte file needed about 18 GB
        no_matrix = (
            "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=1000;\nend;\n"
            "begin sets;\ncharset c = 1 - 1000;\nend;\n"
        )
        # a later dimensions command after a real matrix did the same
        after_matrix = (
            "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=4;\n"
            "format datatype=dna;\nmatrix\nt1 ACGT\n;\nend;\n"
            "begin characters;\ndimensions nchar=1000;\nend;\n"
            "begin sets;\ncharset c = 1 - 1000;\nend;\n"
        )
        for text in (no_matrix, after_matrix):
            for source in (text, StringIO(text)):
                with self.subTest(text=text, source=type(source).__name__):
                    with self.assertRaisesRegex(
                        Nexus.NexusError,
                        r"^Dimensions nchar=1000 is larger than the whole "
                        rf"input \({len(text)} characters\)\.$",
                    ):
                        Nexus.Nexus(source)

    def test_sets_before_matrix(self):
        """Charset and charlabels before the matrix, as Mesquite writes them."""
        nexus = Nexus.Nexus(
            "#NEXUS\nbegin characters;\ndimensions ntax=2 nchar=4;\n"
            "format datatype=dna;\ncharset c = 2 - 4;\ncharlabels 1 a, 4 d;\n"
            "matrix\nt1 ACGT\nt2 ACGA\n;\nend;\n"
        )
        self.assertEqual(nexus.nchar, 4)
        self.assertEqual(nexus.charsets, {"c": [1, 2, 3]})
        self.assertEqual(nexus.charlabels, {0: "a", 3: "d"})
        self.assertEqual(str(nexus.matrix["t2"]), "ACGA")

    @unittest.skipIf(platform.python_implementation() == "PyPy", "no tracemalloc")
    def test_repeated_range_and_name(self):
        """A range or set name repeated in a charset is stored only once."""
        # the repeats were all expanded into one list before duplicates were
        # removed, so memory grew with repeats times nchar: about 16 MB and
        # 4 MB for these 5 KB and 2 KB inputs (a 36 KB file needed 766 MiB)
        import tracemalloc

        head = (
            "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=1000;\n"
            f"format datatype=dna;\nmatrix\nt1 {'A' * 1000}\n;\nend;\n"
            "begin sets;\n"
        )
        for sets, charsets in (
            ("charset c = " + "1-1000 " * 500 + ";\n", {"c": list(range(1000))}),
            (
                "charset c = 1-1000;\ncharset d = " + "c " * 500 + ";\n",
                {"c": list(range(1000)), "d": list(range(1000))},
            ),
        ):
            with self.subTest(sets=sets[:20]):
                was_tracing = tracemalloc.is_tracing()
                if not was_tracing:
                    tracemalloc.start()
                try:
                    tracemalloc.reset_peak()
                    before = tracemalloc.get_traced_memory()[0]
                    nexus = Nexus.Nexus(head + sets + "end;\n")
                    peak = tracemalloc.get_traced_memory()[1] - before
                finally:
                    if not was_tracing:
                        tracemalloc.stop()
                self.assertEqual(nexus.charsets, charsets)
                self.assertLess(peak, 1_000_000)


def parse_in_child(testcase, code, *args):
    """Run code in a child Python that imports this checkout's Bio.

    Return what the child prints, as JSON. The tests of quadratic and cubic
    parsing below use this: the old code took minutes on their input, and
    the child is stopped after 20 seconds, so on that code the test fails
    instead of only running slowly. The fixed code takes a few seconds at most,
    even on a busy machine, so the limit is kept low: old code that takes 7
    minutes here must still fail on a machine many times faster.
    """
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        filter(
            None,
            [os.path.dirname(os.path.dirname(Bio.__file__)), env.get("PYTHONPATH")],
        )
    )
    try:
        result = subprocess.run(
            [sys.executable, "-W", "ignore", "-c", code, *map(str, args)],
            capture_output=True,
            text=True,
            env=env,
            timeout=20,
        )
    except subprocess.TimeoutExpired:
        testcase.fail("parsing did not finish within 20 seconds")
    testcase.assertEqual(result.returncode, 0, result.stderr)
    return json.loads(result.stdout)


class TestLongInput(unittest.TestCase):
    """Long commands, and many commands, parse in linear time.

    CharBuffer, which reads the words of a command, removed the first item
    of a list of the command's characters for each character it read, and
    joined all the characters left for each item of a set; the block reader
    removed the first item of the list of commands for each command. So the
    time grew with the square of the input: each input below took several
    minutes to parse, and now takes a second or two.
    """

    def test_translate(self):
        """A TRANSLATE table of 100,000 taxa."""
        n = 100_000
        code = r"""
import json
import sys
from Bio.Nexus import Nexus
n = int(sys.argv[1])
table = ",\n".join(f"{i} taxon_{i}" for i in range(1, n + 1))
nexus = Nexus.Nexus(
    f"#NEXUS\nbegin trees;\ntranslate\n{table};\ntree t = (1,{n});\nend;\n"
)
print(json.dumps([list(nexus.translate.items()), nexus.trees[0].get_taxa()]))
"""
        translate, taxa = parse_in_child(self, code, n)
        self.assertEqual(translate, [[i, f"taxon_{i}"] for i in range(1, n + 1)])
        self.assertEqual(taxa, ["taxon_1", f"taxon_{n}"])

    def test_charset(self):
        """A CHARSET of 160,000 items, half of them ranges."""
        code = r"""
import json
import sys
from Bio.Nexus import Nexus
n = int(sys.argv[1])
nexus = Nexus.Nexus(
    "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=4;\nformat datatype=dna;\n"
    "matrix\nt1 ACGT\n;\nend;\nbegin sets;\n"
    f"charset c = {'2 3-4 ' * n};\nend;\n"
)
print(json.dumps(nexus.charsets))
"""
        self.assertEqual(parse_in_child(self, code, 80_000), {"c": [1, 2, 3]})

    def test_many_commands(self):
        """A block of 1,600,000 commands."""
        n = 1_600_000
        code = r"""
import json
import sys
from Bio.Nexus import Nexus
n = int(sys.argv[1])
nexus = Nexus.Nexus(
    "#NEXUS\nbegin foo;\n" + "x;" * n + "\nend;\nbegin data;\n"
    "dimensions ntax=1 nchar=4;\nformat datatype=dna;\nmatrix\nt1 ACGT\n;\nend;\n"
)
print(
    json.dumps(
        [
            [[block.title, block.commandlines] for block in nexus.unknown_blocks],
            str(nexus.matrix["t1"]),
        ]
    )
)
"""
        blocks, row = parse_in_child(self, code, n)
        self.assertEqual(blocks, [["foo", [["x"] * n]]])
        self.assertEqual(row, "ACGT")


class TestRepeatedLabels(unittest.TestCase):
    """Repeated taxon labels in a matrix are renamed in linear time.

    The repeats of a label are renamed "a.copy", "a.copy1", "a.copy2" and so
    on. The name for the k-th repeat was found by walking past the k-1 names
    given before it, checking each against a list of all the labels so far,
    so the time grew with the cube of the number of repeats: a matrix of
    5,000 rows named "a", a 20 KB file, took several minutes to parse. Each
    row after the first block of an interleaved matrix was also looked up in
    a dict of all the taxa built again for that row.
    """

    def test_negative_ntax_interleaved(self):
        """An interleaved matrix with a negative ntax raises NexusError."""
        # Each row starts a new block, before the first block has any taxa.
        for ntax in (-1, -2):
            text = (
                "#NEXUS\nbegin data;\n"
                f"dimensions ntax={ntax} nchar=4;\n"
                "format datatype=dna interleave;\n"
                "matrix\nA ACGT\nB ACGT\n;\nend;\n"
            )
            with self.subTest(ntax=ntax):
                with self.assertRaises(Nexus.NexusError) as cm:
                    Nexus.Nexus(text)
                self.assertIn("not in first block", str(cm.exception))

    @staticmethod
    def renamed(n):
        """Return the names n rows named "a" are given."""
        return ["a", "a.copy"] + [f"a.copy{i}" for i in range(1, n - 1)]

    def test_matrix(self):
        """A matrix of 5,000 rows named "a"."""
        n = 5000
        code = r"""
import json
import sys
from Bio.Nexus import Nexus
n = int(sys.argv[1])
nexus = Nexus.Nexus(
    f"#NEXUS\nbegin data;\ndimensions ntax={n} nchar=1;\n"
    "format datatype=dna;\nmatrix\n" + "a A\n" * n + ";\nend;\n"
)
print(json.dumps([nexus.unaltered_taxlabels, nexus.taxlabels, list(nexus.matrix)]))
"""
        unaltered, taxlabels, matrix = parse_in_child(self, code, n)
        self.assertEqual(unaltered, ["a"] * n)
        self.assertEqual(taxlabels, self.renamed(n))
        self.assertEqual(matrix, self.renamed(n))

    def test_interleaved_matrix(self):
        """An interleaved matrix of two blocks of 5,000 rows named "a"."""
        n = 5000
        code = r"""
import json
import sys
from Bio.Nexus import Nexus
n = int(sys.argv[1])
nexus = Nexus.Nexus(
    f"#NEXUS\nbegin data;\ndimensions ntax={n} nchar=2;\n"
    "format datatype=dna interleave;\nmatrix\n" + "a A\n" * n + "a C\n" * n
    + ";\nend;\n"
)
print(
    json.dumps(
        [nexus.taxlabels, {label: str(seq) for label, seq in nexus.matrix.items()}]
    )
)
"""
        taxlabels, matrix = parse_in_child(self, code, n)
        self.assertEqual(taxlabels, self.renamed(n))
        self.assertEqual(matrix, dict.fromkeys(self.renamed(n), "AC"))

    def test_labels_like_copies(self):
        """Labels that already look like copies are renamed as before."""
        read = ["a", "a", "a.copy", "a", "a.copy1", "a.copy", "a", "b.copy", "b"]
        renamed = [
            "a",
            "a.copy",
            "a.copy1",
            "a.copy2",
            "a.copy3",
            "a.copy4",
            "a.copy5",
            "b.copy",
            "b",
        ]
        n = len(read)
        rows = "".join(f"{label} A\n" for label in read)
        for interleave, nchar, matrix in (("", 1, rows), ("interleave", 2, rows * 2)):
            with self.subTest(interleave=interleave):
                nexus = Nexus.Nexus(
                    f"#NEXUS\nbegin data;\ndimensions ntax={n} nchar={nchar};\n"
                    f"format datatype=dna {interleave};\nmatrix\n{matrix};\nend;\n"
                )
                self.assertEqual(nexus.unaltered_taxlabels, read)
                self.assertEqual(nexus.taxlabels, renamed)

    def test_copy_without_number(self):
        """Repeats of a label like "seq.copy_a" or "copycat" get ".copy" added."""
        # The text after "copy" was read as a copy number, so these raised
        # ValueError from int().
        for label in ("seq.copy_a", "copycat", "x.copyx", "x.copy²"):
            renamed = [label, f"{label}.copy", f"{label}.copy1"]
            with self.subTest(label=label, source="matrix"):
                nexus = Nexus.Nexus(
                    "#NEXUS\nbegin data;\ndimensions ntax=3 nchar=1;\n"
                    f"format datatype=dna;\nmatrix\n{label} A\n{label} A\n{label} A\n"
                    ";\nend;\n"
                )
                self.assertEqual(nexus.taxlabels, renamed)
            with self.subTest(label=label, source="add_sequence"):
                nexus = Nexus.Nexus()
                for _ in range(3):
                    nexus.add_sequence(label, "A")
                self.assertEqual(nexus.taxlabels, renamed)

    def test_explicit_copies(self):
        """Rows "a.copy1" to "a.copy15000", then the same rows again."""
        # Each repeat of the second half started its own walk up the chain
        # "a.copy1", "a.copy2", ... past all the names the first half and
        # the repeats before it had taken, which took minutes.
        n = 15_000
        code = r"""
import json
import sys
from Bio.Nexus import Nexus
n = int(sys.argv[1])
rows = "".join(f"a.copy{i} A\n" for i in range(1, n + 1))
nexus = Nexus.Nexus(
    f"#NEXUS\nbegin data;\ndimensions ntax={2 * n} nchar=1;\n"
    f"format datatype=dna;\nmatrix\n{rows * 2};\nend;\n"
)
print(json.dumps(nexus.taxlabels))
"""
        self.assertEqual(
            parse_in_child(self, code, n),
            [f"a.copy{i}" for i in range(1, 2 * n + 1)],
        )

    def test_repeats_meeting_in_one_chain(self):
        """Repeats of different labels in one chain of copies are named as before."""
        read = [
            "a.copy1",
            "a.copy3",
            "a",
            "a.copy",
            "a.copy2",
            "a.copy1",
            "a",
            "a.copy3",
            "a.copy",
            "a",
            "b",
            "a.copy1.copy",
            "a.copy1",
            "b",
        ]
        renamed = [
            "a.copy1",
            "a.copy3",
            "a",
            "a.copy",
            "a.copy2",
            "a.copy4",
            "a.copy5",
            "a.copy6",
            "a.copy7",
            "a.copy8",
            "b",
            "a.copy1.copy",
            "a.copy9",
            "b.copy",
        ]
        n = len(read)
        rows = "".join(f"{label} A\n" for label in read)
        for interleave, nchar, matrix in (("", 1, rows), ("interleave", 2, rows * 2)):
            with self.subTest(interleave=interleave):
                nexus = Nexus.Nexus(
                    f"#NEXUS\nbegin data;\ndimensions ntax={n} nchar={nchar};\n"
                    f"format datatype=dna {interleave};\nmatrix\n{matrix};\nend;\n"
                )
                self.assertEqual(nexus.unaltered_taxlabels, read)
                self.assertEqual(nexus.taxlabels, renamed)
                self.assertEqual(list(nexus.matrix), renamed)
        nexus = Nexus.Nexus()
        for label in read:
            nexus.add_sequence(label, "A")
        self.assertEqual(nexus.taxlabels, renamed)

    def test_add_sequence(self):
        """Adding 4,000 sequences named "a", as the NEXUS writers do."""
        # add_sequence walked the chain of copies from "a" for each repeat,
        # checking each name against the list of all the labels, which took
        # minutes. Writing an alignment as NEXUS calls it for each record.
        n = 4000
        code = r"""
import json
import sys
from Bio.Nexus import Nexus
n = int(sys.argv[1])
nexus = Nexus.Nexus()
for _ in range(n):
    nexus.add_sequence("a", "A")
print(json.dumps([nexus.unaltered_taxlabels, nexus.taxlabels, list(nexus.matrix)]))
"""
        unaltered, taxlabels, matrix = parse_in_child(self, code, n)
        self.assertEqual(unaltered, ["a"] * n)
        self.assertEqual(taxlabels, self.renamed(n))
        self.assertEqual(matrix, self.renamed(n))

    def test_add_sequence_freed_name(self):
        """add_sequence gives again a copy name the caller has removed."""
        nexus = Nexus.Nexus()
        for _ in range(3):
            nexus.add_sequence("a", "A")
        self.assertEqual(nexus.taxlabels, ["a", "a.copy", "a.copy1"])
        nexus.taxlabels.remove("a.copy")
        del nexus.matrix["a.copy"]
        nexus.add_sequence("a", "A")
        self.assertEqual(nexus.taxlabels, ["a", "a.copy1", "a.copy"])

    def test_add_sequence_empty_first(self):
        """add_sequence takes an empty sequence as the first one."""
        # With no matrix yet and nothing to pad, this raised TypeError.
        nexus = Nexus.Nexus()
        nexus.add_sequence("a", "")
        self.assertEqual(nexus.matrix, {"a": Seq("")})
        self.assertEqual(nexus.ntax, 1)
        self.assertEqual(nexus.nchar, 0)
        self.assertEqual(nexus.taxlabels, ["a"])
        nexus.add_sequence("b", "AC")
        self.assertEqual(nexus.matrix, {"a": Seq("--"), "b": Seq("AC")})


class TestWideClade(unittest.TestCase):
    """A clade with many children parses in linear time.

    At each comma between two children of a clade, Tree._parse added up the
    lengths of all the pieces of the clade's text before that comma, so the
    time grew with the square of the number of children: the tree below, one
    clade of 50,000 leaves, took several minutes to parse, and now takes
    under a second.
    """

    def test_star_tree(self):
        """A NEXUS tree of one clade with 50,000 leaves."""
        n = 50_000
        # the commas inside the comments do not separate children
        code = r"""
import json
import sys
from Bio.Nexus import Nexus
n = int(sys.argv[1])
newick = ",".join(f"t{i}[&r={{{i},{i}}}]:{float(i)}" for i in range(n))
nexus = Nexus.Nexus(f"#NEXUS\nbegin trees;\ntree star = ({newick});\nend;\n")
tree = nexus.trees[0]
children = [tree.node(i).data for i in tree.node(tree.root).succ]
print(json.dumps([[data.taxon, data.branchlength, data.comment] for data in children]))
"""
        self.assertEqual(
            parse_in_child(self, code, n),
            [[f"t{i}", float(i), f"[&r={{{i},{i}}}]"] for i in range(n)],
        )

    def test_deep_tree(self):
        """A wide clade nested 900 deep parses in linear time.

        Tree._parse scanned the text of each clade again at every level of
        nesting above it, so the time grew with the tree length times the
        depth: this 1.5 MB tree took about 100 seconds to parse.
        """
        code = r"""
import json
import sys
from Bio.Nexus import Trees
n, depth = map(int, sys.argv[1:])
star = "(" + ",".join(f"t{i}" for i in range(n)) + ")"
tree = Trees.Tree("(" * depth + star + ",x)" * depth + ";")
node, siblings = tree.node(tree.root), []
for level in range(depth):
    inner, outer = (tree.node(i) for i in node.succ)
    siblings.append(outer.data.taxon)
    node = inner
print(json.dumps([siblings, [tree.node(i).data.taxon for i in node.succ]]))
"""
        siblings, leaves = parse_in_child(self, code, 200_000, 900)
        self.assertEqual(siblings, ["x"] * 900)
        self.assertEqual(leaves, [f"t{i}" for i in range(200_000)])


@unittest.skipIf(cnexus is None, "C extension Bio.Nexus.cnexus not available")
@unittest.skipUnless(
    sys.platform.startswith("linux") and sys.maxsize > 2**32,
    "needs RLIMIT_AS to be enforced (Linux, 64-bit)",
)
@unittest.skipUnless(
    platform.python_implementation() == "CPython",
    "PyPy's cpyext copies the input string, so it runs out of memory "
    "before scanfile's own allocation is reached",
)
class TestCnexusAllocationFailure(unittest.TestCase):
    """cnexus.scanfile must raise MemoryError when it cannot allocate."""

    def test_scanfile_out_of_memory(self):
        # Regression test: scanfile set MemoryError when its working buffer
        # could not be allocated, but then carried on and wrote through the
        # NULL pointer (segfault). Cap the child's address space so that a
        # 64 MB input fits but a second 64 MB buffer does not.
        code = """
import resource
import sys
from Bio.Nexus import cnexus
text = "A" * (64 * 1024 * 1024)
cnexus.scanfile("warm up;")
try:
    with open("/proc/self/status") as handle:
        vmsize = next(
            int(line.split()[1]) * 1024
            for line in handle
            if line.startswith("VmSize:")
        )
    limit = vmsize + 16 * 1024 * 1024
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
except (OSError, StopIteration, ValueError) as err:
    print("SKIP", repr(err))
    sys.exit()
try:
    cnexus.scanfile(text)
except MemoryError:
    print("MemoryError")
"""
        # CI runs the suite with AddressSanitizer preloaded, whose allocator
        # aborts on failure by default; ask it to return NULL like malloc.
        env = os.environ.copy()
        asan_options = env.get("ASAN_OPTIONS")
        env["ASAN_OPTIONS"] = ":".join(
            filter(None, [asan_options, "allocator_may_return_null=1"])
        )
        result = subprocess.run(
            [sys.executable, "-W", "ignore", "-c", code],
            capture_output=True,
            text=True,
            env=env,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        if result.stdout.startswith("SKIP"):
            self.skipTest(f"cannot cap the address space: {result.stdout.strip()}")
        self.assertEqual(result.stdout.strip(), "MemoryError")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
