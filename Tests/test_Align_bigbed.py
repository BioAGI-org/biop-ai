# Copyright 2022 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Align.bigbed module."""
import os
import struct
import subprocess
import sys
import tempfile
import unittest
import zlib
from io import BytesIO
from io import StringIO
from unittest import mock

import support

import Bio
from Bio import Align
from Bio import SeqIO
from Bio.Align import Alignments
from Bio.Align import bigbed
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.bigbed."
    ) from None


for i, argument in enumerate(sys.argv):
    if argument == "--big":
        big = True
        sys.argv.pop(i)
        break
else:
    big = False


class TestAlign_dna_rna(unittest.TestCase):
    # The bigBed file dna_rna.bb was generated using the commands
    # sort -k1,1 -k2,2n dna_rna.bed > dna_rna.sorted.bed
    # twoBitInfo hg38.2bit hg38.chrom.sizes
    # bedToBigBed dna_rna.sorted.bed hg38.chrom.sizes dna_rna.bb

    path = support.DATA / "Blat" / "dna_rna.bb"

    def setUp(self):
        data = {}
        with SeqIO.parse(support.DATA / "Blat" / "dna.fa", "fasta") as records:
            for record in records:
                name, start_end = record.id.split(":")
                assert name == "chr3"
                start, end = start_end.split("-")
                start = int(start)
                end = int(end)
                sequence = str(record.seq)
                assert len(sequence) == end - start
                data[start] = sequence
        self.dna = Seq(data, length=198295559)  # hg38 chr3
        with SeqIO.parse(support.DATA / "Blat" / "rna.fa", "fasta") as records:
            self.rna = {record.id: record.seq for record in records}

    def test_reading(self):
        """Test parsing dna_rna.bb."""
        path = support.DATA / "Blat" / "dna_rna.bb"
        with Align.parse(path, "bigbed") as alignments:
            self.check_alignments(alignments)
            alignments = iter(alignments)
            self.check_alignments(alignments)
        with Align.parse(path, "bigbed") as alignments:
            self.check_alignments(alignments)
        with self.assertRaises(AttributeError):
            alignments._stream
        with Align.parse(path, "bigbed") as alignments:
            pass
        with self.assertRaises(AttributeError):
            alignments._stream

    def test_writing(self):
        """Test writing dna_rna.bb."""
        with tempfile.TemporaryFile() as output:
            with Align.parse(self.path, "bigbed") as alignments:
                Align.write(alignments, output, "bigbed")
            output.flush()
            output.seek(0)
            alignments = Align.parse(output, "bigbed")
            self.check_alignments(alignments)

    def check_alignments(self, alignments):
        self.assertEqual(
            str(alignments.declaration),
            """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
""",
        )
        self.assertEqual(len(alignments.targets), 1)
        self.assertEqual(alignments.targets[0].id, "chr3")
        self.assertEqual(len(alignments.targets[0]), 198295559)
        self.assertEqual(len(alignments), 4)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 1711))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_046654.1")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[42530895, 42530958, 42532020,
                           42532095, 42532563, 42532606],
                          [     181,      118,      118,
                                 43,       43,        0]])
                # fmt: on
            )
        )
        alignment.target.seq = self.dna
        alignment.query.seq = self.rna[alignment.query.id]
        self.assertTrue(
            np.array_equal(
                alignment.substitutions,
                # fmt: off
            np.array([[36.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0., 40.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0., 57.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0., 42.,  0.,  0.,  0.,  0.],
                      [ 2.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  3.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                     ])
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGTacgt")
        # The modified RNAs have gaps in their sequence. As this information is
        # not stored in a BED file, we cannot calculate the substitution matrix.
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (181 aligned letters; 175 identities; 6 mismatches; 1530 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 181:
        identities = 175,
        mismatches = 6.
    gaps = 1530:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1530:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1530:
                open_internal_deletions = 2,
                extend_internal_deletions = 1528;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 1530)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 1530)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 1530)
        self.assertEqual(counts.gaps, 1530)
        self.assertEqual(counts.aligned, 181)
        self.assertEqual(counts.identities, 175)
        self.assertEqual(counts.mismatches, 6)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 978)
        self.assertEqual(alignment.shape, (2, 1711))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_046654.1_modified")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[42530895, 42530922, 42530958, 42532020, 42532037,
                           42532039, 42532095, 42532563, 42532606],
                          [     179,      152,      116,      116,       99,
                                 99,       43,       43,        0]])
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (179 aligned letters; 0 identities; 0 mismatches; 1532 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 179:
        identities = 0,
        mismatches = 0.
    gaps = 1532:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 1532:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 1532:
                open_internal_deletions = 3,
                extend_internal_deletions = 1529;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 1532)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 1532)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 1532)
        self.assertEqual(counts.gaps, 1532)
        self.assertEqual(counts.aligned, 179)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 5407))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_111921.1")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[48663767, 48663813, 48665640,
                           48665722, 48669098, 48669174],
                          [       0,       46,       46,
                                128,      128,      204]])
                # fmt: on
            )
        )
        alignment.target.seq = self.dna
        alignment.query.seq = self.rna[alignment.query.id]
        self.assertTrue(
            np.array_equal(
                alignment.substitutions,
                # fmt: off
            np.array([[53.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0., 35.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0., 50.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0., 27.,  0.,  0.,  0.,  0.],
                      [ 9.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  7.,  0.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0., 16.,  0.,  0.,  0.,  0.,  0.],
                      [ 0.,  0.,  0.,  7.,  0.,  0.,  0.,  0.],
                     ])
            )
        )
        self.assertEqual(alignment.substitutions.alphabet, "ACGTacgt")
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (204 aligned letters; 165 identities; 39 mismatches; 5203 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 204:
        identities = 165,
        mismatches = 39.
    gaps = 5203:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 5203:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 5203:
                open_internal_deletions = 2,
                extend_internal_deletions = 5201;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 5203)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 5203)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 5203)
        self.assertEqual(counts.gaps, 5203)
        self.assertEqual(counts.aligned, 204)
        self.assertEqual(counts.identities, 165)
        self.assertEqual(counts.mismatches, 39)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 972)
        self.assertEqual(alignment.shape, (2, 5407))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "NR_111921.1_modified")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[48663767, 48663795, 48663796, 48663813, 48665640,
                           48665716, 48665722, 48669098, 48669174],
                          [       0,       28,       28,       45,       45,
                                121,      127,      127,      203]])
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (203 aligned letters; 0 identities; 0 mismatches; 5204 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 203:
        identities = 0,
        mismatches = 0.
    gaps = 5204:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 5204:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 5204:
                open_internal_deletions = 3,
                extend_internal_deletions = 5201;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 5204)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 5204)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 5204)
        self.assertEqual(counts.gaps, 5204)
        self.assertEqual(counts.aligned, 203)
        self.assertRaises(StopIteration, next, alignments)


class TestAlign_dna(unittest.TestCase):
    # The bigBed file psl_34_001.bb was generated using the commands
    # sort -k1,1 -k2,2n psl_34_001.bed > psl_34_001.sorted.bed
    # twoBitInfo hg19.2bit hg19.chrom.sizes
    # bedToBigBed psl_34_001.sorted.bed hg19.chrom.sizes psl_34_001.bb

    # The bigBed file psl_34_003.bb was generated using the commands
    # sort -k1,1 -k2,2n psl_34_003.bed > psl_34_003.sorted.bed
    # twoBitInfo hg19.2bit hg19.chrom.sizes
    # bedToBigBed psl_34_003.sorted.bed hg19.chrom.sizes psl_34_003.bb

    # The bigBed file psl_34_004.bb was generated using the commands
    # sort -k1,1 -k2,2n psl_34_004.bed > psl_34_004.sorted.bed
    # twoBitInfo hg19.2bit hg19.chrom.sizes
    # bedToBigBed psl_34_004.sorted.bed hg19.chrom.sizes psl_34_004.bb

    # The bigBed file psl_34_005.bb was generated using the commands
    # sort -k1,1 -k2,2n psl_34_005.bed > psl_34_005.sorted.bed
    # twoBitInfo hg19.2bit hg19.chrom.sizes
    # bedToBigBed psl_34_005.sorted.bed hg19.chrom.sizes psl_34_005.bb

    def check_alignments_psl_34_001(self, alignments):
        self.assertEqual(
            str(alignments.declaration),
            """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
""",
        )
        self.assertEqual(len(alignments.targets), 10)
        self.assertEqual(alignments.targets[0].id, "chr1")
        self.assertEqual(len(alignments.targets[0]), 249250621)
        self.assertEqual(alignments.targets[1].id, "chr10")
        self.assertEqual(len(alignments.targets[1]), 135534747)
        self.assertEqual(alignments.targets[2].id, "chr13")
        self.assertEqual(len(alignments.targets[2]), 115169878)
        self.assertEqual(alignments.targets[3].id, "chr18")
        self.assertEqual(len(alignments.targets[3]), 78077248)
        self.assertEqual(alignments.targets[4].id, "chr19")
        self.assertEqual(len(alignments.targets[4]), 59128983)
        self.assertEqual(alignments.targets[5].id, "chr2")
        self.assertEqual(len(alignments.targets[5]), 243199373)
        self.assertEqual(alignments.targets[6].id, "chr22")
        self.assertEqual(len(alignments.targets[6]), 51304566)
        self.assertEqual(alignments.targets[7].id, "chr4")
        self.assertEqual(len(alignments.targets[7]), 191154276)
        self.assertEqual(alignments.targets[8].id, "chr8")
        self.assertEqual(len(alignments.targets[8]), 146364022)
        self.assertEqual(alignments.targets[9].id, "chr9")
        self.assertEqual(len(alignments.targets[9]), 141213431)
        self.assertEqual(len(alignments), 22)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 50))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[1207056, 1207106],
                          [      0,      50]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (50 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 50:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 50)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 33))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[10271783, 10271816],
                          [       0,       33]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (33 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 33:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 33)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 946)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[39368490, 39368526],
                          [      36,        0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 36)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 824)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[61700837, 61700871],
                          [       0,       34]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (34 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 34:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 34)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 942)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[220325687, 220325721],
                          [       34,         0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (34 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 34:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 34)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 834)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[99388555, 99388591],
                          [      36,        0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 36)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 920)
        self.assertEqual(alignment.shape, (2, 25))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[112178171, 112178196],
                          [       25,         0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (25 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 25:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 25)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 912)
        self.assertEqual(alignment.shape, (2, 51))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[52759147, 52759154, 52759160, 52759198],
                          [       0,        7,        7,       45]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (45 aligned letters; 0 identities; 0 mismatches; 6 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 45:
        identities = 0,
        mismatches = 0.
    gaps = 6:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 6:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 6:
                open_internal_deletions = 1,
                extend_internal_deletions = 5;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 6)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 6)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 6)
        self.assertEqual(counts.gaps, 6)
        self.assertEqual(counts.aligned, 45)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[23891310, 23891349],
                          [       0,       39]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (39 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 39:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 39)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 930)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[43252217, 43252245],
                          [       0,       28]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (28 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 28:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 28)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 848)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[553742, 553781],
                          [    39,      0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (39 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 39:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 39)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 890)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[35483340, 35483365, 35483499, 35483510],
                          [       0,       25,       25,       36]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 134 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 134:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 134:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 134:
                open_internal_deletions = 1,
                extend_internal_deletions = 133;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 134)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 134)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 134)
        self.assertEqual(counts.gaps, 134)
        self.assertEqual(counts.aligned, 36)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[54017130, 54017169],
                          [      39,        0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (39 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 39:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 39)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[53575980, 53575997],
                          [      17,        0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 17)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 946)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[120641740, 120641776],
                          [       36,         0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 36)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 682)
        self.assertEqual(alignment.shape, (2, 44))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[183925984, 183925990, 183926028],
                          [        0,         6,        44]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (44 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 44:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 44)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 834)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[42144400, 42144436],
                          [       0,       36]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 36)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 892)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[48997405, 48997442],
                          [      37,        0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (37 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 37:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 37)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 572)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[37558157, 37558167, 37558173, 37558191],
                          [      28,       18,       18,        0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (28 aligned letters; 0 identities; 0 mismatches; 6 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 28:
        identities = 0,
        mismatches = 0.
    gaps = 6:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 6:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 6:
                open_internal_deletions = 1,
                extend_internal_deletions = 5;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 6)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 6)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 6)
        self.assertEqual(counts.gaps, 6)
        self.assertEqual(counts.aligned, 28)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[61646095, 61646111],
                          [       0,       16]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (16 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 16:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 16)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[95160479, 95160520],
                          [       0,       41]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (41 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 41:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 41)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 854)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[85737865, 85737906],
                          [       0,       41]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (41 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 41:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 41)
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_001(self):
        """Test reading psl_34_001.bb."""
        path = support.DATA / "Blat" / "psl_34_001.bb"
        with Align.parse(path, "bigbed") as alignments:
            self.check_alignments_psl_34_001(alignments)

    def test_writing_psl_34_001(self):
        """Test writing psl_34_001.bb."""
        path = support.DATA / "Blat" / "psl_34_001.bb"
        with tempfile.TemporaryFile() as output:
            with Align.parse(path, "bigbed") as alignments:
                Align.write(alignments, output, "bigbed")
            output.flush()
            output.seek(0)
            alignments = Align.parse(output, "bigbed")
            self.check_alignments_psl_34_001(alignments)

    def check_alignments_psl_34_003(self, alignments):
        self.assertEqual(
            str(alignments.declaration),
            """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
""",
        )
        self.assertEqual(len(alignments.targets), 3)
        self.assertEqual(alignments.targets[0].id, "chr1")
        self.assertEqual(len(alignments.targets[0]), 249250621)
        self.assertEqual(alignments.targets[1].id, "chr2")
        self.assertEqual(len(alignments.targets[1]), 243199373)
        self.assertEqual(alignments.targets[2].id, "chr4")
        self.assertEqual(len(alignments.targets[2]), 191154276)
        self.assertEqual(len(alignments), 3)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 33))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[10271783, 10271816],
                          [       0,       33]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (33 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 33:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 33)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[53575980, 53575997],
                          [      17,        0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 17)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[61646095, 61646111],
                          [       0,       16]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (16 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 16:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 16)
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_003(self):
        """Test reading psl_34_003.bb."""
        path = support.DATA / "Blat" / "psl_34_003.bb"
        with Align.parse(path, "bigbed") as alignments:
            self.check_alignments_psl_34_003(alignments)

    def test_writing_psl_34_003(self):
        """Test writing psl_34_003.bb."""
        path = support.DATA / "Blat" / "psl_34_003.bb"
        with tempfile.TemporaryFile() as output:
            with Align.parse(path, "bigbed") as alignments:
                Align.write(alignments, output, "bigbed")
            output.flush()
            output.seek(0)
            alignments = Align.parse(output, "bigbed")
            self.check_alignments_psl_34_003(alignments)

    def check_alignments_psl_34_004(self, alignments):
        self.assertEqual(
            str(alignments.declaration),
            """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
""",
        )
        self.assertEqual(len(alignments.targets), 10)
        self.assertEqual(alignments.targets[0].id, "chr1")
        self.assertEqual(len(alignments.targets[0]), 249250621)
        self.assertEqual(alignments.targets[1].id, "chr10")
        self.assertEqual(len(alignments.targets[1]), 135534747)
        self.assertEqual(alignments.targets[2].id, "chr13")
        self.assertEqual(len(alignments.targets[2]), 115169878)
        self.assertEqual(alignments.targets[3].id, "chr18")
        self.assertEqual(len(alignments.targets[3]), 78077248)
        self.assertEqual(alignments.targets[4].id, "chr19")
        self.assertEqual(len(alignments.targets[4]), 59128983)
        self.assertEqual(alignments.targets[5].id, "chr2")
        self.assertEqual(len(alignments.targets[5]), 243199373)
        self.assertEqual(alignments.targets[6].id, "chr22")
        self.assertEqual(len(alignments.targets[6]), 51304566)
        self.assertEqual(alignments.targets[7].id, "chr4")
        self.assertEqual(len(alignments.targets[7]), 191154276)
        self.assertEqual(alignments.targets[8].id, "chr8")
        self.assertEqual(len(alignments.targets[8]), 146364022)
        self.assertEqual(alignments.targets[9].id, "chr9")
        self.assertEqual(len(alignments.targets[9]), 141213431)
        self.assertEqual(len(alignments), 19)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 50))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[1207056, 1207106],
                          [      0,      50]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (50 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 50:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 50)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 946)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[39368490, 39368526],
                          [      36,        0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 36)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 824)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[61700837, 61700871],
                          [       0,       34]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (34 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 34:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 34)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 942)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[220325687, 220325721],
                          [       34,         0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (34 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 34:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 34)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 834)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[99388555, 99388591],
                          [      36,        0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 36)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 920)
        self.assertEqual(alignment.shape, (2, 25))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[112178171, 112178196],
                          [       25,         0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (25 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 25:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 25)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 912)
        self.assertEqual(alignment.shape, (2, 51))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[52759147, 52759154, 52759160, 52759198],
                          [       0,        7,        7,       45]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (45 aligned letters; 0 identities; 0 mismatches; 6 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 45:
        identities = 0,
        mismatches = 0.
    gaps = 6:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 6:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 6:
                open_internal_deletions = 1,
                extend_internal_deletions = 5;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 6)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 6)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 6)
        self.assertEqual(counts.gaps, 6)
        self.assertEqual(counts.aligned, 45)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[23891310, 23891349],
                          [       0,       39]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (39 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 39:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 39)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 930)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[43252217, 43252245],
                          [       0,       28]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (28 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 28:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 28)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 848)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[553742, 553781],
                          [    39,      0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (39 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 39:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 39)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 890)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[35483340, 35483365, 35483499, 35483510],
                          [       0,       25,       25,       36]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 134 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 134:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 134:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 134:
                open_internal_deletions = 1,
                extend_internal_deletions = 133;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 134)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 134)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 134)
        self.assertEqual(counts.gaps, 134)
        self.assertEqual(counts.aligned, 36)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[54017130, 54017169],
                          [      39,        0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (39 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 39:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 39)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 946)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[120641740, 120641776],
                          [       36,         0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 36)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 682)
        self.assertEqual(alignment.shape, (2, 44))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[183925984, 183925990, 183926028],
                          [        0,         6,        44]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (44 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 44:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 44)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 834)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[42144400, 42144436],
                          [       0,       36]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 36)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 892)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[48997405, 48997442],
                          [      37,        0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (37 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 37:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 37)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 572)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[37558157, 37558167, 37558173, 37558191],
                          [      28,       18,       18,        0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (28 aligned letters; 0 identities; 0 mismatches; 6 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 28:
        identities = 0,
        mismatches = 0.
    gaps = 6:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 6:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 6:
                open_internal_deletions = 1,
                extend_internal_deletions = 5;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 6)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 6)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 6)
        self.assertEqual(counts.gaps, 6)
        self.assertEqual(counts.aligned, 28)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[95160479, 95160520],
                          [       0,       41]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (41 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 41:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 41)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 854)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[85737865, 85737906],
                          [       0,       41]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (41 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 41:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 41)
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_004(self):
        """Test reading psl_34_004.bb."""
        path = support.DATA / "Blat" / "psl_34_004.bb"
        with Align.parse(path, "bigbed") as alignments:
            self.check_alignments_psl_34_004(alignments)

    def test_writing_psl_34_004(self):
        """Test writing psl_34_004.bb."""
        path = support.DATA / "Blat" / "psl_34_004.bb"
        with tempfile.TemporaryFile() as output:
            with Align.parse(path, "bigbed") as alignments:
                Align.write(alignments, output, "bigbed")
            output.flush()
            output.seek(0)
            alignments = Align.parse(output, "bigbed")
            self.check_alignments_psl_34_004(alignments)

    def check_alignments_psl_34_005(self, alignments):
        self.assertEqual(
            str(alignments.declaration),
            """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
""",
        )
        self.assertEqual(len(alignments.targets), 10)
        self.assertEqual(alignments.targets[0].id, "chr1")
        self.assertEqual(len(alignments.targets[0]), 249250621)
        self.assertEqual(alignments.targets[1].id, "chr10")
        self.assertEqual(len(alignments.targets[1]), 135534747)
        self.assertEqual(alignments.targets[2].id, "chr13")
        self.assertEqual(len(alignments.targets[2]), 115169878)
        self.assertEqual(alignments.targets[3].id, "chr18")
        self.assertEqual(len(alignments.targets[3]), 78077248)
        self.assertEqual(alignments.targets[4].id, "chr19")
        self.assertEqual(len(alignments.targets[4]), 59128983)
        self.assertEqual(alignments.targets[5].id, "chr2")
        self.assertEqual(len(alignments.targets[5]), 243199373)
        self.assertEqual(alignments.targets[6].id, "chr22")
        self.assertEqual(len(alignments.targets[6]), 51304566)
        self.assertEqual(alignments.targets[7].id, "chr4")
        self.assertEqual(len(alignments.targets[7]), 191154276)
        self.assertEqual(alignments.targets[8].id, "chr8")
        self.assertEqual(len(alignments.targets[8]), 146364022)
        self.assertEqual(alignments.targets[9].id, "chr9")
        self.assertEqual(len(alignments.targets[9]), 141213431)
        self.assertEqual(len(alignments), 22)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 50))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[1207056, 1207106],
                          [      0,      50]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (50 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 50:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 50)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 33))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[10271783, 10271816],
                          [       0,       33]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (33 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 33:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 33)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 946)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[39368490, 39368526],
                          [      36,        0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 36)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 824)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[61700837, 61700871],
                          [       0,       34]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (34 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 34:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 34)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 942)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[220325687, 220325721],
                          [       34,         0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (34 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 34:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 34)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 834)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[99388555, 99388591],
                          [      36,        0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 36)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 920)
        self.assertEqual(alignment.shape, (2, 25))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr10")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[112178171, 112178196],
                          [       25,         0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (25 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 25:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 25)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 912)
        self.assertEqual(alignment.shape, (2, 51))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[52759147, 52759154, 52759160, 52759198],
                          [       0,        7,        7,       45]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (45 aligned letters; 0 identities; 0 mismatches; 6 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 45:
        identities = 0,
        mismatches = 0.
    gaps = 6:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 6:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 6:
                open_internal_deletions = 1,
                extend_internal_deletions = 5;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 6)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 6)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 6)
        self.assertEqual(counts.gaps, 6)
        self.assertEqual(counts.aligned, 45)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[23891310, 23891349],
                          [       0,       39]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (39 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 39:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 39)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 930)
        self.assertEqual(alignment.shape, (2, 28))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr18")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[43252217, 43252245],
                          [       0,       28]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (28 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 28:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 28)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 848)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[553742, 553781],
                          [    39,      0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (39 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 39:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 39)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 890)
        self.assertEqual(alignment.shape, (2, 170))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[35483340, 35483365, 35483499, 35483510],
                          [       0,       25,       25,       36]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 134 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 134:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 134:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 134:
                open_internal_deletions = 1,
                extend_internal_deletions = 133;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 134)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 134)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 134)
        self.assertEqual(counts.gaps, 134)
        self.assertEqual(counts.aligned, 36)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 39))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr19")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[54017130, 54017169],
                          [      39,        0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (39 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 39:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 39)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 17))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[53575980, 53575997],
                          [      17,        0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (17 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 17:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 17)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 946)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[120641740, 120641776],
                          [       36,         0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 36)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 682)
        self.assertEqual(alignment.shape, (2, 44))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[183925984, 183925990, 183926028],
                          [        0,         6,        44]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (44 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 44:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 44)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 834)
        self.assertEqual(alignment.shape, (2, 36))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[42144400, 42144436],
                          [       0,       36]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (36 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 36:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 36)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 892)
        self.assertEqual(alignment.shape, (2, 37))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr22")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[48997405, 48997442],
                          [      37,        0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (37 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 37:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 37)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 572)
        self.assertEqual(alignment.shape, (2, 34))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[37558157, 37558167, 37558173, 37558191],
                          [      28,       18,       18,        0]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (28 aligned letters; 0 identities; 0 mismatches; 6 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 28:
        identities = 0,
        mismatches = 0.
    gaps = 6:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 6:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 6:
                open_internal_deletions = 1,
                extend_internal_deletions = 5;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 6)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 6)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 6)
        self.assertEqual(counts.gaps, 6)
        self.assertEqual(counts.aligned, 28)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 16))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "hg18_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[61646095, 61646111],
                          [       0,       16]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (16 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 16:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 16)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr8")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[95160479, 95160520],
                          [       0,       41]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (41 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 41:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 41)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 854)
        self.assertEqual(alignment.shape, (2, 41))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr9")
        self.assertEqual(alignment.query.id, "hg19_dna")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[85737865, 85737906],
                          [       0,       41]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (41 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 41:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 41)
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_34_005(self):
        """Test reading psl_34_005.bb."""
        path = support.DATA / "Blat" / "psl_34_005.bb"
        with Align.parse(path, "bigbed") as alignments:
            self.check_alignments_psl_34_005(alignments)

    def test_writing_psl_34_005(self):
        """Test writing psl_34_005.bb."""
        path = support.DATA / "Blat" / "psl_34_005.bb"
        with tempfile.TemporaryFile() as output:
            with Align.parse(path, "bigbed") as alignments:
                Align.write(alignments, output, "bigbed")
            output.flush()
            output.seek(0)
            alignments = Align.parse(output, "bigbed")
            self.check_alignments_psl_34_005(alignments)


class TestAlign_dnax_prot(unittest.TestCase):
    # The bigBed file psl_35_001.bb was generated using the commands
    # sort -k1,1 -k2,2n psl_35_001.bed > psl_35_001.sorted.bed
    # twoBitInfo hg38.2bit hg38.chrom.sizes
    # bedToBigBed psl_35_001.sorted.bed hg38.chrom.sizes psl_35_001.bb

    path = support.DATA / "Blat" / "psl_35_001.bb"

    def check_alignments(self, alignments):
        self.assertEqual(
            str(alignments.declaration),
            """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
""",
        )
        self.assertEqual(len(alignments.targets), 2)
        self.assertEqual(alignments.targets[0].id, "chr13")
        self.assertEqual(len(alignments.targets[0]), 114364328)
        self.assertEqual(alignments.targets[1].id, "chr4")
        self.assertEqual(len(alignments.targets[1]), 190214555)
        self.assertEqual(len(alignments), 8)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 986)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[75549820, 75549865, 75567225, 75567312],
                          [       0,       45,       45,      132]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (132 aligned letters; 0 identities; 0 mismatches; 17360 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 132:
        identities = 0,
        mismatches = 0.
    gaps = 17360:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 17360:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 17360:
                open_internal_deletions = 1,
                extend_internal_deletions = 17359;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 17360)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 17360)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 17360)
        self.assertEqual(counts.gaps, 17360)
        self.assertEqual(counts.aligned, 132)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[75560749, 75560881],
                          [       0,      132]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (132 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 132:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 132)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[75566694, 75566850],
                          [       0,      156]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (156 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 156:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 156)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[75569459, 75569507],
                          [       0,       48]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (48 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 48:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 48)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[75594914, 75594989],
                          [       0,       75]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (75 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 75:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 75)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1000)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr13")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[75604767, 75604827, 75605728, 75605809],
                          [       0,       60,       60,      141]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (141 aligned letters; 0 identities; 0 mismatches; 901 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 141:
        identities = 0,
        mismatches = 0.
    gaps = 901:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 901:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 901:
                open_internal_deletions = 1,
                extend_internal_deletions = 900;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 901)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 901)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 901)
        self.assertEqual(counts.gaps, 901)
        self.assertEqual(counts.aligned, 141)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 166)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[41257605, 41257731, 41263227, 41263290],
                          [       0,      126,      126,      189]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (189 aligned letters; 0 identities; 0 mismatches; 5496 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 189:
        identities = 0,
        mismatches = 0.
    gaps = 5496:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 5496:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 5496:
                open_internal_deletions = 1,
                extend_internal_deletions = 5495;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 5496)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 5496)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 5496)
        self.assertEqual(counts.gaps, 5496)
        self.assertEqual(counts.aligned, 189)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 530)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr4")
        self.assertEqual(alignment.query.id, "CAG33136.1")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates,
                # fmt: off
                np.array([[41260685, 41260787],
                          [       0,      102]]),
                # fmt: on
            )
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (102 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 102:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 102)
        self.assertRaises(StopIteration, next, alignments)

    def test_reading_psl_35_001(self):
        """Test parsing psl_35_001.bb."""
        with Align.parse(self.path, "bigbed") as alignments:
            self.check_alignments(alignments)

    def test_writing_psl_35_001(self):
        """Test writing psl_35_001.bb."""
        with tempfile.TemporaryFile() as output:
            with Align.parse(self.path, "bigbed") as alignments:
                Align.write(alignments, output, "bigbed")
            output.flush()
            output.seek(0)
            alignments = Align.parse(output, "bigbed")
            self.check_alignments(alignments)


class TestAlign_bed12(unittest.TestCase):
    # The bigBed files were generated using the commands
    # twoBitInfo hg19.2bit hg19.chrom.sizes
    # bedToBigBed bed3.bed hg19.chrom.sizes bed3.bb
    # bedToBigBed bed4.bed hg19.chrom.sizes bed4.bb
    # bedToBigBed bed5.bed hg19.chrom.sizes bed5.bb
    # bedToBigBed bed6.bed hg19.chrom.sizes bed6.bb
    # bedToBigBed bed7.bed hg19.chrom.sizes bed7.bb
    # bedToBigBed bed8.bed hg19.chrom.sizes bed8.bb
    # bedToBigBed bed9.bed hg19.chrom.sizes bed9.bb
    # bedToBigBed bed12.bed hg19.chrom.sizes bed12.bb

    def check_autosql(self, declaration, bedN, msg):
        if bedN == 3:
            self.assertEqual(
                str(declaration),
                """\
table bed
"Browser Extensible Data"
(
   string chrom;         "Reference sequence chromosome or scaffold"
   uint   chromStart;    "Start position in chromosome"
   uint   chromEnd;      "End position in chromosome"
)
""",
                msg=msg,
            )
        elif bedN == 4:
            self.assertEqual(
                str(declaration),
                """\
table bed
"Browser Extensible Data"
(
   string chrom;         "Reference sequence chromosome or scaffold"
   uint   chromStart;    "Start position in chromosome"
   uint   chromEnd;      "End position in chromosome"
   string name;          "Name of item."
)
""",
                msg=msg,
            )
        elif bedN == 5:
            self.assertEqual(
                str(declaration),
                """\
table bed
"Browser Extensible Data"
(
   string chrom;         "Reference sequence chromosome or scaffold"
   uint   chromStart;    "Start position in chromosome"
   uint   chromEnd;      "End position in chromosome"
   string name;          "Name of item."
   uint   score;         "Score (0-1000)"
)
""",
                msg=msg,
            )
        elif bedN == 6:
            self.assertEqual(
                str(declaration),
                """\
table bed
"Browser Extensible Data"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position in chromosome"
   uint    chromEnd;      "End position in chromosome"
   string  name;          "Name of item."
   uint    score;         "Score (0-1000)"
   char[1] strand;        "+ or - for strand"
)
""",
                msg=msg,
            )
        elif bedN == 7:
            self.assertEqual(
                str(declaration),
                """\
table bed
"Browser Extensible Data"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position in chromosome"
   uint    chromEnd;      "End position in chromosome"
   string  name;          "Name of item."
   uint    score;         "Score (0-1000)"
   char[1] strand;        "+ or - for strand"
   uint    thickStart;    "Start of where display should be thick (start codon)"
)
""",
                msg=msg,
            )
        elif bedN == 8:
            self.assertEqual(
                str(declaration),
                """\
table bed
"Browser Extensible Data"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position in chromosome"
   uint    chromEnd;      "End position in chromosome"
   string  name;          "Name of item."
   uint    score;         "Score (0-1000)"
   char[1] strand;        "+ or - for strand"
   uint    thickStart;    "Start of where display should be thick (start codon)"
   uint    thickEnd;      "End of where display should be thick (stop codon)"
)
""",
                msg=msg,
            )
        elif bedN == 9:
            self.assertEqual(
                str(declaration),
                """\
table bed
"Browser Extensible Data"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position in chromosome"
   uint    chromEnd;      "End position in chromosome"
   string  name;          "Name of item."
   uint    score;         "Score (0-1000)"
   char[1] strand;        "+ or - for strand"
   uint    thickStart;    "Start of where display should be thick (start codon)"
   uint    thickEnd;      "End of where display should be thick (stop codon)"
   uint    reserved;      "Used as itemRgb as of 2004-11-22"
)
""",
                msg=msg,
            )
        elif bedN == 12:
            self.assertEqual(
                str(declaration),
                """\
table bed
"Browser Extensible Data"
(
   string          chrom;          "Reference sequence chromosome or scaffold"
   uint            chromStart;     "Start position in chromosome"
   uint            chromEnd;       "End position in chromosome"
   string          name;           "Name of item."
   uint            score;          "Score (0-1000)"
   char[1]         strand;         "+ or - for strand"
   uint            thickStart;     "Start of where display should be thick (start codon)"
   uint            thickEnd;       "End of where display should be thick (stop codon)"
   uint            reserved;       "Used as itemRgb as of 2004-11-22"
   int             blockCount;     "Number of blocks"
   int[blockCount] blockSizes;     "Comma separated list of block sizes"
   int[blockCount] chromStarts;    "Start positions relative to chromStart"
)
""",
                msg=msg,
            )

    def check_alignments(self, alignments, bedN, msg):
        self.assertEqual(len(alignments), 2)
        alignment = next(alignments)
        if bedN >= 5:
            self.assertEqual(alignment.score, 960, msg=msg)
        self.assertEqual(alignment.shape, (2, 4000), msg=msg)
        self.assertLess(
            alignment.coordinates[0, 0], alignment.coordinates[0, -1], msg=msg
        )
        self.assertLess(
            alignment.coordinates[1, 0], alignment.coordinates[1, -1], msg=msg
        )
        self.assertEqual(len(alignment), 2, msg=msg)
        self.assertIs(alignment.sequences[0], alignment.target, msg=msg)
        self.assertIs(alignment.sequences[1], alignment.query, msg=msg)
        self.assertEqual(alignment.target.id, "chr22", msg=msg)
        if bedN >= 4:
            self.assertEqual(alignment.query.id, "mRNA1", msg=msg)
        else:
            self.assertIsNone(alignment.query.id, msg=msg)
        if bedN == 12:
            self.assertTrue(
                np.array_equal(
                    alignment.coordinates,
                    # fmt: off
                    np.array([[1000, 1567, 4512, 5000],
                              [   0,  567,  567, 1055]]),
                    # fmt: on
                ),
                msg=msg,
            )
        else:
            self.assertTrue(
                np.array_equal(
                    alignment.coordinates,
                    np.array([[1000, 5000], [0, 4000]]),
                ),
                msg=msg,
            )
        if bedN >= 7:
            self.assertEqual(alignment.thickStart, 1200, msg=msg)
        if bedN >= 8:
            self.assertEqual(alignment.thickEnd, 4900, msg=msg)
        if bedN >= 9:
            self.assertEqual(alignment.itemRgb, "255,0,0", msg=msg)
        alignment = next(alignments)
        if bedN >= 5:
            self.assertEqual(alignment.score, 900, msg=msg)
        self.assertEqual(alignment.shape, (2, 4000), msg=msg)
        self.assertLess(
            alignment.coordinates[0, 0], alignment.coordinates[0, -1], msg=msg
        )
        if bedN >= 6:
            self.assertGreater(
                alignment.coordinates[1, 0],
                alignment.coordinates[1, -1],
                msg=msg,
            )
        else:
            self.assertLess(
                alignment.coordinates[1, 0],
                alignment.coordinates[1, -1],
                msg=msg,
            )
        self.assertEqual(len(alignment), 2, msg=msg)
        self.assertIs(alignment.sequences[0], alignment.target, msg=msg)
        self.assertIs(alignment.sequences[1], alignment.query, msg=msg)
        self.assertEqual(alignment.target.id, "chr22", msg=msg)
        if bedN >= 4:
            self.assertEqual(alignment.query.id, "mRNA2", msg=msg)
        else:
            self.assertIsNone(alignment.query.id, msg=msg)
        if bedN == 12:
            self.assertTrue(
                np.array_equal(
                    alignment.coordinates,
                    # fmt: off
                    np.array([[2000, 2433, 5601, 6000],
                              [ 832,  399,  399,    0]])
                    # fmt: on
                ),
                msg=msg,
            )
        elif bedN >= 6:
            self.assertTrue(
                np.array_equal(
                    alignment.coordinates,
                    np.array([[2000, 6000], [4000, 0]]),
                ),
                msg=msg,
            )
        else:
            self.assertTrue(
                np.array_equal(
                    alignment.coordinates,
                    np.array([[2000, 6000], [0, 4000]]),
                ),
                msg=msg,
            )
        if bedN >= 7:
            self.assertEqual(alignment.thickStart, 2300, msg=msg)
        if bedN >= 8:
            self.assertEqual(alignment.thickEnd, 5960, msg=msg)
        if bedN >= 9:
            self.assertEqual(alignment.itemRgb, "0,255,0", msg=msg)
        with self.assertRaises(StopIteration) as cm:
            next(alignments)
            self.fail("More than two alignments reported")

    def test_reading(self):
        """Test parsing alignments in file formats BED3 through BED12."""
        for bedN in (3, 4, 5, 6, 7, 8, 9, 12):
            filename = "bed%d.bb" % bedN
            path = support.DATA / "Blat" / filename
            msg = "bed%d" % bedN
            with Align.parse(path, "bigbed") as alignments:
                self.check_autosql(alignments.declaration, bedN, msg)
                self.check_alignments(alignments, bedN, msg)

    def test_writing(self):
        """Test Writing alignments in file formats BED3 through BED12."""
        for bedN in (3, 4, 5, 6, 7, 8, 9, 12):
            filename = "bed%d.bb" % bedN
            path = support.DATA / "Blat" / filename
            with tempfile.TemporaryFile() as output:
                with Align.parse(path, "bigbed") as alignments:
                    Align.write(alignments, output, "bigbed", bedN=bedN)
                output.flush()
                output.seek(0)
                alignments = Align.parse(output, "bigbed")
                msg = "bed%d" % bedN
                self.check_autosql(alignments.declaration, bedN, msg)
                self.check_alignments(alignments, bedN, msg)


class TestAlign_extended_bed(unittest.TestCase):
    # The bigBed files bigbed_extended.littleendian.bb and
    # bigbed_extended.bigendian.bb are BED9+2 files, with nine predefined BED
    # fields and 2 extra (custom) fields. It was created by running
    #
    # bedToBigBed -as=bedExample2.as -type=bed9+2 -extraIndex=name,geneSymbol bedExample2.bed hg18.chrom.sizes bigbed_extended.bb
    #
    # where bedExample2.bed contains 10 lines selected from the example BED9+2
    # file bedExample2.bed downloaded from UCSC
    # (https://genome.ucsc.edu/goldenPath/help/examples/bedExample2.bed)
    # and bedExample2.as the associated AutoSQL file, also downloaded from UCSC
    # (https://genome.ucsc.edu/goldenPath/help/examples/bedExample2.as)
    # declaring the nine predefined BED fields and the two extra fields
    #
    # The bigBed file bigbed_extended.littleendian.bb was generated on a
    # little-endian machine; the bigBed file bigbed_extended.bigendian.bb was
    # generated on a big-endian machine.

    def test_reading(self):
        """Test parsing bigbed_extended.bb."""
        path = support.DATA / "Blat" / "bigbed_extended.littleendian.bb"
        with Align.parse(path, "bigbed") as alignments:
            self.assertEqual(alignments.byteorder, "<")
            self.check_alignments(alignments)
        path = support.DATA / "Blat" / "bigbed_extended.bigendian.bb"
        with Align.parse(path, "bigbed") as alignments:
            self.assertEqual(alignments.byteorder, ">")
            self.check_alignments(alignments)

    def check_alignments(self, alignments):
        self.assertEqual(
            str(alignments.declaration),
            """\
table hg18KGchr7
"UCSC Genes for chr7 with color plus GeneSymbol and SwissProtID"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position of feature on chromosome"
   uint    chromEnd;      "End position of feature on chromosome"
   string  name;          "Name of gene"
   uint    score;         "Score"
   char[1] strand;        "+ or - for strand"
   uint    thickStart;    "Coding region start"
   uint    thickEnd;      "Coding region end"
   uint    reserved;      "Green on + strand, Red on - strand"
   string  geneSymbol;    "Gene Symbol"
   string  spID;          "SWISS-PROT protein Accession number"
)
""",
        )
        self.assertEqual(len(alignments.targets), 1)
        self.assertEqual(alignments.targets[0].id, "chr7")
        self.assertEqual(len(alignments.targets[0]), 158821424)
        self.assertEqual(len(alignments), 10)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 60328)
        self.assertEqual(alignment.thickEnd, 60328)
        self.assertEqual(alignment.itemRgb, "255,0,0")
        self.assertEqual(alignment.shape, (2, 1241))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc010krx.1")
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[60328, 61569], [1241, 0]]))
        )
        self.assertEqual(alignment.annotations["geneSymbol"], ".")
        self.assertEqual(alignment.annotations["spID"], "PDGFA")
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (1241 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 1241:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 1241)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 506606)
        self.assertEqual(alignment.thickEnd, 525164)
        self.assertEqual(alignment.itemRgb, "255,0,0")
        self.assertEqual(alignment.shape, (2, 22585))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc003sir.1")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates, np.array([[503422, 526007], [22585, 0]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], "PDGFA")
        self.assertEqual(alignment.annotations["spID"], "P04085")
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (22585 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 22585:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 22585)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 504726)
        self.assertEqual(alignment.thickEnd, 525164)
        self.assertEqual(alignment.itemRgb, "255,0,0")
        self.assertEqual(alignment.shape, (2, 22585))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc003sis.1")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates, np.array([[503422, 526007], [22585, 0]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], "PDGFA")
        self.assertEqual(alignment.annotations["spID"], "P04085-2")
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (22585 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 22585:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 22585)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 507195)
        self.assertEqual(alignment.thickEnd, 518820)
        self.assertEqual(alignment.itemRgb, "255,0,0")
        self.assertEqual(alignment.shape, (2, 12690))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc003sit.1")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates, np.array([[506940, 519630], [12690, 0]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], "PDGFA")
        self.assertEqual(alignment.annotations["spID"], "Q32M96")
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (12690 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 12690:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 12690)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 556592)
        self.assertEqual(alignment.thickEnd, 717668)
        self.assertEqual(alignment.itemRgb, "255,0,0")
        self.assertEqual(alignment.shape, (2, 162747))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc003siu.1")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates, np.array([[555912, 718659], [162747, 0]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], "PRKAR1B")
        self.assertEqual(alignment.annotations["spID"], "Q8N422")
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (162747 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 162747:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 162747)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 556592)
        self.assertEqual(alignment.thickEnd, 717668)
        self.assertEqual(alignment.itemRgb, "255,0,0")
        self.assertEqual(alignment.shape, (2, 163357))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc003siv.1")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates, np.array([[555912, 719269], [163357, 0]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], "PRKAR1B")
        self.assertEqual(alignment.annotations["spID"], "Q8N422")
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (163357 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 163357:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 163357)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 556592)
        self.assertEqual(alignment.thickEnd, 717668)
        self.assertEqual(alignment.itemRgb, "255,0,0")
        self.assertEqual(alignment.shape, (2, 177901))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc003siw.1")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates, np.array([[555912, 733813], [177901, 0]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], "PRKAR1B")
        self.assertEqual(alignment.annotations["spID"], "Q8N422")
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (177901 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 177901:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 177901)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 585418)
        self.assertEqual(alignment.thickEnd, 585418)
        self.assertEqual(alignment.itemRgb, "255,0,0")
        self.assertEqual(alignment.shape, (2, 22329))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc003six.1")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates, np.array([[585418, 607747], [22329, 0]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], ".")
        self.assertEqual(alignment.annotations["spID"], "PRKAR1B")
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (22329 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 22329:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 22329)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 733217)
        self.assertEqual(alignment.thickEnd, 791816)
        self.assertEqual(alignment.itemRgb, "0,255,0")
        self.assertEqual(alignment.shape, (2, 59779))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc003siz.2")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates, np.array([[732863, 792642], [0, 59779]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], ".")
        self.assertEqual(alignment.annotations["spID"], "DKFZp762F1415")
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (59779 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 59779:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 59779)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 0)
        self.assertEqual(alignment.thickStart, 732883)
        self.assertEqual(alignment.thickEnd, 791816)
        self.assertEqual(alignment.itemRgb, "0,255,0")
        self.assertEqual(alignment.shape, (2, 59779))
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr7")
        self.assertEqual(alignment.query.id, "uc010krz.1")
        self.assertTrue(
            np.array_equal(
                alignment.coordinates, np.array([[732863, 792642], [0, 59779]])
            )
        )
        self.assertEqual(alignment.annotations["geneSymbol"], "HEATR2")
        self.assertEqual(alignment.annotations["spID"], "Q86Y56")
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (59779 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 59779:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 59779)

    def test_writing(self):
        """Test writing bigbed_extended.bb."""
        byteorder = sys.byteorder  # "little" or "big"
        path = support.DATA / "Blat" / f"bigbed_extended.{byteorder}endian.bb"
        with open(path, "rb") as stream:
            correct = stream.read()
        with open(support.DATA / "Blat" / "bedExample2.as") as stream:
            autosql_data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(autosql_data)
        with tempfile.TemporaryFile() as output:
            with Align.parse(path, "bigbed") as alignments:
                Align.write(
                    alignments,
                    output,
                    "bigbed",
                    bedN=9,
                    declaration=declaration,
                    extraIndex=["name", "geneSymbol"],
                )
            output.flush()
            output.seek(0)
            alignments = Align.parse(output, "bigbed")
            self.check_alignments(alignments)
        with tempfile.TemporaryFile() as output:
            with Align.parse(path, "bigbed") as alignments:
                targets = alignments.targets
                Align.write(
                    alignments,
                    output,
                    "bigbed",
                    bedN=9,
                    declaration=declaration,
                    targets=targets,
                    extraIndex=["name", "geneSymbol"],
                )
            output.flush()
            output.seek(0)
            alignments = Align.parse(output, "bigbed")
            self.check_alignments(alignments)


class TestAlign_searching(unittest.TestCase):
    path = support.DATA / "Blat" / "bigbedtest.bb"

    # The bigBed file bigbedtest.bb contains the following data:
    # chr1     10     100     name1   1       +
    # chr1     29      39     name2   2       -
    # chr1    200     300     name3   3       +
    # chr2     50      50     name4   6       +
    # chr2    100     110     name5   4       +
    # chr2    200     210     name6   5       +
    # chr2    220     220     name7   6       +
    # chr3      0       0     name8   7       -

    def check_alignments(self, alignments):
        self.assertEqual(
            str(alignments.declaration),
            """\
table bed
"Browser Extensible Data"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position in chromosome"
   uint    chromEnd;      "End position in chromosome"
   string  name;          "Name of item."
   uint    score;         "Score (0-1000)"
   char[1] strand;        "+ or - for strand"
)
""",
        )
        self.assertEqual(len(alignments.targets), 3)
        self.assertEqual(alignments.targets[0].id, "chr1")
        self.assertEqual(len(alignments.targets[0]), 1000)
        self.assertEqual(alignments.targets[1].id, "chr2")
        self.assertEqual(len(alignments.targets[1]), 2000)
        self.assertEqual(alignments.targets[2].id, "chr3")
        self.assertEqual(len(alignments.targets[2]), 1000)
        self.assertEqual(len(alignments), 8)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 1)
        self.assertEqual(alignment.shape, (2, 90))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "name1")
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[10, 100], [0, 90]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (90 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 90:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 90)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 2)
        self.assertEqual(alignment.shape, (2, 10))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "name2")
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[29, 39], [10, 0]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (10 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 10:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 10)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 3)
        self.assertEqual(alignment.shape, (2, 100))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr1")
        self.assertEqual(alignment.query.id, "name3")
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[200, 300], [0, 100]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (100 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 100:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 100)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 6)
        self.assertEqual(alignment.shape, (2, 0))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "name4")
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[50, 50], [0, 0]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (0 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 0:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 0)
        self.assertEqual(counts.identities, 0)
        self.assertEqual(counts.mismatches, 0)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 4)
        self.assertEqual(alignment.shape, (2, 10))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "name5")
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[100, 110], [0, 10]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (10 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 10:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 10)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 5)
        self.assertEqual(alignment.shape, (2, 10))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "name6")
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[200, 210], [0, 10]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (10 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 10:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 10)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 6)
        self.assertEqual(alignment.shape, (2, 0))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr2")
        self.assertEqual(alignment.query.id, "name7")
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[220, 220], [0, 0]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (0 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 0:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 0)
        self.assertEqual(counts.identities, 0)
        self.assertEqual(counts.mismatches, 0)
        alignment = next(alignments)
        self.assertEqual(alignment.score, 7)
        self.assertEqual(alignment.shape, (2, 0))
        self.assertEqual(len(alignment), 2)
        self.assertIs(alignment.sequences[0], alignment.target)
        self.assertIs(alignment.sequences[1], alignment.query)
        self.assertEqual(alignment.target.id, "chr3")
        self.assertEqual(alignment.query.id, "name8")
        self.assertTrue(
            np.array_equal(alignment.coordinates, np.array([[0, 0], [0, 0]]))
        )
        counts = alignment.counts()
        self.assertEqual(
            repr(counts),
            "<AlignmentCounts object (0 aligned letters; 0 identities; 0 mismatches; 0 gaps) at 0x%x>"
            % id(counts),
        )
        self.assertEqual(
            str(counts),
            """\
AlignmentCounts object with
    aligned = 0:
        identities = 0,
        mismatches = 0.
    gaps = 0:
        left_gaps = 0:
            left_insertions = 0:
                open_left_insertions = 0,
                extend_left_insertions = 0;
            left_deletions = 0:
                open_left_deletions = 0,
                extend_left_deletions = 0;
        internal_gaps = 0:
            internal_insertions = 0:
                open_internal_insertions = 0,
                extend_internal_insertions = 0;
            internal_deletions = 0:
                open_internal_deletions = 0,
                extend_internal_deletions = 0;
        right_gaps = 0:
            right_insertions = 0:
                open_right_insertions = 0,
                extend_right_insertions = 0;
            right_deletions = 0:
                open_right_deletions = 0,
                extend_right_deletions = 0.
""",
        )
        self.assertEqual(counts.left_insertions, 0)
        self.assertEqual(counts.left_deletions, 0)
        self.assertEqual(counts.right_insertions, 0)
        self.assertEqual(counts.right_deletions, 0)
        self.assertEqual(counts.internal_insertions, 0)
        self.assertEqual(counts.internal_deletions, 0)
        self.assertEqual(counts.left_gaps, 0)
        self.assertEqual(counts.right_gaps, 0)
        self.assertEqual(counts.internal_gaps, 0)
        self.assertEqual(counts.insertions, 0)
        self.assertEqual(counts.deletions, 0)
        self.assertEqual(counts.gaps, 0)
        self.assertEqual(counts.aligned, 0)
        self.assertEqual(counts.identities, 0)
        self.assertEqual(counts.mismatches, 0)

    def test_reading(self):
        """Test reading bigbedtest.bb."""
        with Align.parse(self.path, "bigbed") as alignments:
            self.check_alignments(alignments)

    def test_writing(self):
        """Test writing bigbedtest.bb."""
        with tempfile.TemporaryFile() as output:
            with Align.parse(self.path, "bigbed") as alignments:
                Align.write(alignments, output, "bigbed", bedN=6)
            output.flush()
            output.seek(0)
            alignments = Align.parse(output, "bigbed")
            self.check_alignments(alignments)

    def test_search_chromosome(self):
        with Align.parse(self.path, "bigbed") as alignments:
            self.assertEqual(
                str(alignments.declaration),
                """\
table bed
"Browser Extensible Data"
(
   string  chrom;         "Reference sequence chromosome or scaffold"
   uint    chromStart;    "Start position in chromosome"
   uint    chromEnd;      "End position in chromosome"
   string  name;          "Name of item."
   uint    score;         "Score (0-1000)"
   char[1] strand;        "+ or - for strand"
)
""",
            )
            selected_alignments = alignments.search("chr2")
            names = [alignment.query.id for alignment in selected_alignments]
        self.assertEqual(names, ["name4", "name5", "name6", "name7"])

    def test_search_region(self):
        with Align.parse(self.path, "bigbed") as alignments:
            selected_alignments = alignments.search("chr2", 105, 1000)
            names = [alignment.query.id for alignment in selected_alignments]
            self.assertEqual(names, ["name5", "name6", "name7"])
            selected_alignments = alignments.search("chr2", 110, 1000)
            names = [alignment.query.id for alignment in selected_alignments]
            self.assertEqual(names, ["name6", "name7"])
            selected_alignments = alignments.search("chr2", 40, 50)
            names = [alignment.query.id for alignment in selected_alignments]
            self.assertEqual(names, ["name4"])
            selected_alignments = alignments.search("chr2", 50, 50)
            names = [alignment.query.id for alignment in selected_alignments]
            self.assertEqual(names, ["name4"])
            selected_alignments = alignments.search("chr2", 50, 200)
            names = [alignment.query.id for alignment in selected_alignments]
            self.assertEqual(names, ["name4", "name5"])
            selected_alignments = alignments.search("chr2", 200, 220)
            names = [alignment.query.id for alignment in selected_alignments]
            self.assertEqual(names, ["name6", "name7"])
            selected_alignments = alignments.search("chr2", 220, 220)
            names = [alignment.query.id for alignment in selected_alignments]
            self.assertEqual(names, ["name7"])

    def test_search_position(self):
        with Align.parse(self.path, "bigbed") as alignments:
            selected_alignments = alignments.search("chr1", 250)
            names = [alignment.query.id for alignment in selected_alignments]
            self.assertEqual(names, ["name3"])

    def test_search_all(self):
        """Search with no chromosome returns every alignment."""
        # This used to raise UnboundLocalError, also for bigMaf and bigPsl
        # files, which inherit search from bigbed.AlignmentIterator.
        with Align.parse(self.path, "bigbed") as alignments:
            names = [alignment.query.id for alignment in alignments.search()]
        self.assertEqual(names, [f"name{i}" for i in range(1, 9)])
        path = support.DATA / "MAF" / "ucsc_test.bb"
        with Align.parse(path, "bigmaf") as alignments:
            expected = [format(alignment, "maf") for alignment in alignments]
            found = [format(alignment, "maf") for alignment in alignments.search()]
        self.assertEqual(len(found), 3)
        self.assertEqual(found, expected)

    def test_three_iterators(self):
        """Create three iterators and use them concurrently."""
        with Align.parse(self.path, "bigbed") as alignments1:
            alignments2 = alignments1.search("chr2")
            alignments3 = alignments1.search("chr2", 110, 1000)
            alignment1 = next(alignments1)
            self.assertEqual(alignment1.query.id, "name1")
            alignment1 = next(alignments1)
            self.assertEqual(alignment1.query.id, "name2")
            alignment2 = next(alignments2)
            self.assertEqual(alignment2.query.id, "name4")
            alignment2 = next(alignments2)
            self.assertEqual(alignment2.query.id, "name5")
            alignment2 = next(alignments2)
            self.assertEqual(alignment2.query.id, "name6")
            alignment3 = next(alignments3)
            self.assertEqual(alignment3.query.id, "name6")
            alignment3 = next(alignments3)
            self.assertEqual(alignment3.query.id, "name7")
            alignment1 = next(alignments1)
            self.assertEqual(alignment1.query.id, "name3")
            alignment1 = next(alignments1)
            self.assertEqual(alignment1.query.id, "name4")
            alignment1 = next(alignments1)
            self.assertEqual(alignment1.query.id, "name5")
            alignment2 = next(alignments2)
            self.assertEqual(alignment2.query.id, "name7")
            self.assertRaises(StopIteration, next, alignments2)
            alignment1 = next(alignments1)
            self.assertEqual(alignment1.query.id, "name6")
            alignment1 = next(alignments1)
            self.assertEqual(alignment1.query.id, "name7")
            self.assertRaises(StopIteration, next, alignments3)
            alignment1 = next(alignments1)
            self.assertEqual(alignment1.query.id, "name8")
            self.assertRaises(StopIteration, next, alignments1)


class BinaryTestBaseClass(unittest.TestCase):
    def assertBinaryEqual(self, file1, file2):
        blocksize = 1024
        n = 0
        while True:
            data1 = file1.read(blocksize)
            data2 = file2.read(blocksize)
            if data1 == b"" and data2 == b"":
                return
            n1 = len(data1)
            if data1 == data2:
                n += n1
                continue
            n2 = len(data2)
            if n1 < n2:
                return self.fail(f"unequal file sizes: {n1} bytes vs >= {n2} bytes")
            if n1 > n2:
                return self.fail(f"unequal file sizes: >= {n1} bytes vs {n2} bytes")
            for i, (c1, c2) in enumerate(zip(data1, data2)):
                if c1 != c2:
                    return self.fail(f"bytes at position {n + i} differ: {c1} vs {c2}")


@unittest.skipUnless(big is True, "big file; use --big to run")
class TestAlign_big(BinaryTestBaseClass):
    # BED files were downloaded from the UCSC table browser:
    #
    # ucsc.bed contains the GENCODE V43 Basic gene annotations for human genome
    # assembly hg38.
    #
    # anoGam3.bed contains the AUGUSTUS gene annotations for genome assembly
    # anoGam3 of Anopheles gambiae (African malaria mosquito).
    #
    # ailMel1.bed contains the NCBI RefSeq All gene annotations for genome
    # assembly ailMel1 of Ailuropoda melanoleuca (giant panda).
    #
    # bisBis1.bed contains the AUGUSTUS gene annotations for genome assembly
    # bisBis1 of Bison bison bison (American bison).
    #
    # bigBed files were generated using bedToBigBed v. 2.9 found in
    # jksrc.v445.zip of the 'kent' source tree provided by UCSC
    # (http://hgdownload.cse.ucsc.edu/admin/).

    def test_a_compressed(self):
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -k1,1 -k2,2n ucsc.clean.bed -o ucsc.clean.bed
        # bedToBigBed -as=Blat/bed12.as ucsc.clean.bed Align/hg38.chrom.sizes ucsc.bb
        with open(support.DATA / "Blat" / "bed12.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_b_uncompressed(self):
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -k1,1 -k2,2n ucsc.clean.bed -o ucsc.clean.bed
        # bedToBigBed -as=Blat/bed12.as -unc ucsc.clean.bed Align/hg38.chrom.sizes ucsc.unc.bb
        with open(support.DATA / "Blat" / "bed12.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.unc.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                declaration=declaration,
                compress=False,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_c_bed3(self):
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -k1,1 -k2,2n ucsc.clean.bed -o ucsc.clean.bed
        # cut -f 1-3 ucsc.clean.bed > ucsc.bed3.bed
        # bedToBigBed -as=Blat/bed3.as -type=bed3 ucsc.bed3.bed Align/hg38.chrom.sizes ucsc.bed3.bb
        with open(support.DATA / "Blat" / "bed3.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed3.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=3,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_d_bed4(self):
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -k1,1 -k2,2n ucsc.clean.bed -o ucsc.clean.bed
        # cut -f 1-4 ucsc.clean.bed > ucsc.bed4.bed
        # bedToBigBed -as=Blat/bed4.as -type=bed4 ucsc.bed4.bed Align/hg38.chrom.sizes ucsc.bed4.bb
        with open(support.DATA / "Blat" / "bed4.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed4.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=4,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_e_bed5(self):
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -k1,1 -k2,2n ucsc.clean.bed -o ucsc.clean.bed
        # cut -f 1-5 ucsc.clean.bed > ucsc.bed5.bed
        # bedToBigBed -as=Blat/bed5.as -type=bed5 ucsc.bed5.bed Align/hg38.chrom.sizes ucsc.bed5.bb
        with open(support.DATA / "Blat" / "bed5.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed5.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=5,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_f_bed6(self):
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -k1,1 -k2,2n ucsc.clean.bed -o ucsc.clean.bed
        # cut -f 1-6 ucsc.clean.bed > ucsc.bed6.bed
        # bedToBigBed -as=Blat/bed6.as -type=bed6 ucsc.bed6.bed Align/hg38.chrom.sizes ucsc.bed6.bb
        with open(support.DATA / "Blat" / "bed6.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed6.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=6,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_g_bed7(self):
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -k1,1 -k2,2n ucsc.clean.bed -o ucsc.clean.bed
        # cut -f 1-7 ucsc.clean.bed > ucsc.bed7.bed
        # bedToBigBed -as=Blat/bed7.as -type=bed7 ucsc.bed7.bed Align/hg38.chrom.sizes ucsc.bed7.bb
        with open(support.DATA / "Blat" / "bed7.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed7.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=7,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_h_bed8(self):
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -k1,1 -k2,2n ucsc.clean.bed -o ucsc.clean.bed
        # cut -f 1-8 ucsc.clean.bed > ucsc.bed8.bed
        # bedToBigBed -as=Blat/bed8.as -type=bed8 ucsc.bed8.bed Align/hg38.chrom.sizes ucsc.bed8.bb
        with open(support.DATA / "Blat" / "bed8.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed8.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=8,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_i_bed9(self):
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -k1,1 -k2,2n ucsc.clean.bed -o ucsc.clean.bed
        # cut -f 1-9 ucsc.clean.bed > ucsc.bed9.bed
        # bedToBigBed -as=Blat/bed9.as -type=bed9 ucsc.bed9.bed Align/hg38.chrom.sizes ucsc.bed9.bb
        with open(support.DATA / "Blat" / "bed9.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        bigBedFileName = "ucsc.bed9.bb"
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                bedN=9,
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_j_extraindex(self):
        # If the names in the BED file are not unique, the binary contents of
        # the generated bigBed file will depend on the exact implementation of
        # the quicksort algorithm in the Standard Library of C (used by
        # bedToBigBed) and the quicksort algorithm in numpy (used by Biopython).
        # To prevent spurious errors when comparing bigBed files created by
        # bedToBigBed and by Biopython, we first remove lines with duplicated
        # names from the BED file.
        # grep -E -v 'fix|alt' Blat/ucsc.bed > ucsc.clean.bed
        # sort -u -k4,4 ucsc.clean.bed | sort -k1,1 -k2,2n > ucsc.unique.bed
        # bedToBigBed -as=Blat/bed12.as -extraIndex=name ucsc.unique.bed Align/hg38.chrom.sizes ucsc.indexed.bb
        with open(support.DATA / "Blat" / "bed12.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ucsc.indexed.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                declaration=declaration,
                extraIndex=["name"],
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_k_anogam(self):
        # sort -k1,1 -k2,2n Blat/anoGam3.bed -o anoGam3.bed
        # bedToBigBed -as=Blat/bed12.as anoGam3.bed Blat/anoGam3.chrom.sizes anoGam3.bb
        with open(support.DATA / "Blat" / "bed12.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "anoGam3.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_l_ailmel(self):
        # sort -k1,1 -k2,2n Blat/ailMel1.bed -o ailMel1.bed
        # bedToBigBed -as=Blat/bed12.as ailMel1.bed Blat/ailMel1.chrom.sizes ailMel1.bb
        with open(support.DATA / "Blat" / "bed12.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "ailMel1.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)

    def test_m_bisbis(self):
        # sort -k1,1 -k2,2n bisBis1.bed -o bisBis1.bed
        # bedToBigBed -as=Blat/bed12.as bisBis1.bed Blat/bisBis1.chrom.sizes bisBis1.bb
        with open(support.DATA / "Blat" / "bed12.as") as stream:
            data = stream.read()
        declaration = bigbed.AutoSQLTable.from_string(data)
        bigBedFileName = "bisBis1.bb"
        alignments = Align.parse(bigBedFileName, "bigbed")
        with tempfile.TemporaryFile() as output, open(bigBedFileName, "rb") as stream:
            Align.write(
                alignments,
                output,
                "bigbed",
                declaration=declaration,
                compress=True,
            )
            output.flush()
            output.seek(0)
            self.assertBinaryEqual(output, stream)


class TestDeclarations(unittest.TestCase):
    def test_declarations(self):
        for length in (3, 4, 5, 6, 7, 8, 9, 12):
            filename = "bed%d.as" % length
            path = support.DATA / "Blat" / filename
            with open(path) as stream:
                data = stream.read()
            declaration = bigbed.AutoSQLTable.from_string(data)
            self.assertEqual(declaration.name, "bed", msg=filename)
            self.assertEqual(
                declaration.comment, "Browser Extensible Data", msg=filename
            )
            self.assertEqual(len(declaration), length, msg=filename)
            field = declaration[0]
            self.assertEqual(field.as_type, "string", msg=filename)
            self.assertEqual(field.name, "chrom", msg=filename)
            self.assertEqual(
                field.comment, "Reference sequence chromosome or scaffold", msg=filename
            )
            field = declaration[1]
            self.assertEqual(field.as_type, "uint", msg=filename)
            self.assertEqual(field.name, "chromStart", msg=filename)
            self.assertEqual(
                field.comment, "Start position in chromosome", msg=filename
            )
            field = declaration[2]
            self.assertEqual(field.as_type, "uint", msg=filename)
            self.assertEqual(field.name, "chromEnd", msg=filename)
            self.assertEqual(field.comment, "End position in chromosome", msg=filename)
            if length == 3:
                return
            field = declaration[3]
            self.assertEqual(field.as_type, "string", msg=filename)
            self.assertEqual(field.name, "name", msg=filename)
            self.assertEqual(field.comment, "Name of item.", msg=filename)
            if length == 4:
                return
            field = declaration[4]
            self.assertEqual(field.as_type, "uint", msg=filename)
            self.assertEqual(field.name, "score", msg=filename)
            self.assertEqual(field.comment, "Score (0-1000)", msg=filename)
            if length == 5:
                return
            field = declaration[5]
            self.assertEqual(field.as_type, "char[1]", msg=filename)
            self.assertEqual(field.name, "strand", msg=filename)
            self.assertEqual(field.comment, "+ or - for strand", msg=filename)
            if length == 6:
                return
            field = declaration[6]
            self.assertEqual(field.as_type, "uint", msg=filename)
            self.assertEqual(field.name, "thickStart", msg=filename)
            self.assertEqual(
                field.comment,
                "Start of where display should be thick (start codon)",
                msg=filename,
            )
            if length == 7:
                return
            field = declaration[7]
            self.assertEqual(field.as_type, "uint", msg=filename)
            self.assertEqual(field.name, "thickEnd", msg=filename)
            self.assertEqual(
                field.comment,
                "End of where display should be thick (stop codon)",
                msg=filename,
            )
            if length == 8:
                return
            field = declaration[8]
            self.assertEqual(field.as_type, "uint", msg=filename)
            self.assertEqual(field.name, "reserved", msg=filename)
            self.assertEqual(
                field.comment, "Used as itemRgb as of 2004-11-22", msg=filename
            )
            if length == 9:
                return
            field = declaration[9]
            self.assertEqual(field.as_type, "int", msg=filename)
            self.assertEqual(field.name, "blockCount", msg=filename)
            self.assertEqual(field.comment, "Number of blocks", msg=filename)
            field = declaration[10]
            self.assertEqual(field.as_type, "int[blockCount]", msg=filename)
            self.assertEqual(field.name, "blockSizes", msg=filename)
            self.assertEqual(
                field.comment, "Comma separated list of block sizes", msg=filename
            )
            field = declaration[11]
            self.assertEqual(field.as_type, "int[blockCount]", msg=filename)
            self.assertEqual(field.name, "chromStarts", msg=filename)
            self.assertEqual(
                field.comment, "Start positions relative to chromStart", msg=filename
            )


class TestAlign_undeclared_fields(unittest.TestCase):
    """Items with more fields than the file header declares."""

    def test_bed4_declared_as_bed3(self):
        with open(support.DATA / "Blat" / "bed4.bb", "rb") as stream:
            data = bytearray(stream.read())
        # fieldCount and definedFieldCount are at byte 32 of the file header.
        # Declare 3 fields, all standard; the items still carry a name.
        struct.pack_into("<HH", data, 32, 3, 3)
        with self.assertRaisesRegex(
            ValueError,
            "^The bigBed item at 1000-5000 has 4 fields, but the header declares 3$",
        ):
            list(Align.parse(BytesIO(data), "bigbed"))

    def test_bed6_declared_with_one_custom_field(self):
        with open(support.DATA / "Blat" / "bed6.bb", "rb") as stream:
            data = bytearray(stream.read())
        # Declare 4 fields, the last one custom; the items carry 6. The two
        # extra values used to be dropped silently.
        struct.pack_into("<HH", data, 32, 4, 3)
        with self.assertRaisesRegex(
            ValueError,
            "^The bigBed item at 1000-5000 has 6 fields, but the header declares 4$",
        ):
            list(Align.parse(BytesIO(data), "bigbed"))


class TestAlign_no_alignments(unittest.TestCase):
    """Writing and reading a bigBed file without alignments."""

    def test_round_trip(self):
        # Writing used to raise AttributeError. As bedToBigBed does, the file
        # has no zoom levels and its R tree index is a header without nodes.
        for kwargs in ({}, {"compress": False}, {"bedN": 4, "extraIndex": ["name"]}):
            with self.subTest(**kwargs):
                alignments = Alignments([])
                alignments.targets = [SeqRecord(Seq(None, length=1000), id="chr1")]
                stream = BytesIO()
                Align.write(alignments, stream, "bigbed", **kwargs)
                data = stream.getvalue()
                # The zoom level count is at byte 6 of the file header, and
                # the offset of the R tree index at byte 24. The index header
                # of 48 bytes is followed by the signature closing the file
                # when there are no extra indices.
                self.assertEqual(struct.unpack_from("<H", data, 6), (0,))
                (index_offset,) = struct.unpack_from("<Q", data, 24)
                if "extraIndex" not in kwargs:
                    self.assertEqual(len(data), index_offset + 48 + 4)
                stream.seek(0)
                with Align.parse(stream, "bigbed") as alignments:
                    self.assertEqual(len(alignments), 0)
                    self.assertEqual(alignments.targets, [])
                    self.assertEqual(list(alignments), [])
                    self.assertEqual(list(alignments.search()), [])

    def test_empty_chromosome_index(self):
        # The chromosome index of a file with no alignments is its 32-byte
        # header alone. Reading it used to read a root node after the header;
        # in a bigBed file the item count 0 that starts the data section made
        # that node look like an empty one, but other bytes there failed.
        dtype = [("name", "S1"), ("id", "=i4"), ("size", "=i4")]
        for after in (b"", b"\xff" * 16):
            with self.subTest(after=after):
                stream = BytesIO()
                bigbed._BPlusTreeFormatter().write(np.zeros(0, dtype), 0, stream)
                self.assertEqual(stream.tell(), 32)
                stream.write(after)
                stream.seek(0)
                self.assertEqual(bigbed._BPlusTreeFormatter().read(stream), [])

    def test_counts_disagree(self):
        # The item count of the file starts the data section (its offset is
        # at byte 16 of the file header); the item count of the R tree index
        # is at byte 8 of its header (its offset is at byte 24). A file of 2
        # items whose index claimed none used to read as empty.
        with open(support.DATA / "Blat" / "bed3.bb", "rb") as stream:
            data = bytes(stream.read())
        (data_offset,) = struct.unpack_from("<Q", data, 16)
        (index_offset,) = struct.unpack_from("<Q", data, 24)
        no_blocks = bytearray(data)
        struct.pack_into("<Q", no_blocks, index_offset + 8, 0)
        no_items = bytearray(data)
        struct.pack_into("<Q", no_items, data_offset, 0)
        for bad, message in (
            (no_blocks, "declares 2 items, but its R tree index has no data blocks"),
            (no_items, "declares 0 items, but its R tree index has some data blocks"),
        ):
            with self.subTest(message=message):
                with self.assertRaisesRegex(ValueError, f"^The file {message}$"):
                    Align.parse(BytesIO(bytes(bad)), "bigbed")


class TestAlign_no_zoom_levels(unittest.TestCase):
    """Searching a bigBed file without zoom levels."""

    def test_search(self):
        # The number of items per data block is read from the zoom level
        # indices only, so search() used to raise AttributeError when there
        # were none.
        with open(support.DATA / "Blat" / "bed3.bb", "rb") as stream:
            data = bytearray(stream.read())
        # The zoom level count is at byte 6 of the file header, and the table
        # of 10 zoom levels of 24 bytes each follows the 64-byte header.
        struct.pack_into("<H", data, 6, 0)
        data[64 : 64 + 240] = bytes(240)
        with Align.parse(BytesIO(data), "bigbed") as alignments:
            self.assertEqual(len(list(alignments)), 2)
            self.assertEqual(len(list(alignments.search("chr22"))), 2)


class TestAlign_search_blocks_read(unittest.TestCase):
    """Searching reads only the data blocks that overlap the query."""

    def test_search(self):
        # After reading a leaf, search() moved on to the next sibling leaf
        # without checking its bounds, so it used to read every later data
        # block under the same parent: 20 blocks for the first query.
        lines = "".join(f"chr1\t{i * 100}\t{i * 100 + 50}\tf{i}\n" for i in range(20))
        alignments = Alignments(list(Align.parse(StringIO(lines), "bed")))
        chromosome = SeqRecord(Seq(None, length=10000), id="chr1")
        for alignment in alignments:
            alignment.sequences[0] = chromosome
        alignments.targets = [chromosome]
        with tempfile.TemporaryFile() as stream:
            Align.write(alignments, stream, "bigbed", bedN=4, itemsPerSlot=1)
            for start, end, names in (
                (0, 150, ["f0", "f1"]),
                (950, 1050, ["f10"]),
                (1900, 2000, ["f19"]),
            ):
                with self.subTest(start=start, end=end):
                    stream.seek(0)
                    alignments = Align.parse(stream, "bigbed")
                    with mock.patch.object(
                        alignments, "_read_block", wraps=alignments._read_block
                    ) as read_block:
                        found = alignments.search("chr1", start, end)
                        found = [alignment.query.id for alignment in found]
                    self.assertEqual(found, names)
                    self.assertLessEqual(read_block.call_count, 2)


class TestAlign_header_offsets(unittest.TestCase):
    """File header offsets that point past the end of the file."""

    def test_past_end_of_file(self):
        # The offsets of the chromosome index, the data count, the R tree
        # index and the first zoom level index are at bytes 8, 16, 24 and
        # 80 of the file. Each used to raise struct.error when it pointed
        # too near the end of the file, and OverflowError when it was too
        # large to seek to.
        with open(support.DATA / "Blat" / "bed12.bb", "rb") as stream:
            data = stream.read()
        size = len(data)
        fields = (
            (8, "chromosome index", 32),
            (16, "data count", 8),
            (24, "R tree index", 48),
            (80, "zoom level index", 48),
        )
        for position, name, length in fields:
            for offset in (size, size - length + 1, 2**63 + 5):
                with self.subTest(name=name, offset=offset):
                    bad = bytearray(data)
                    struct.pack_into("<Q", bad, position, offset)
                    with self.assertRaises(ValueError) as cm:
                        Align.parse(BytesIO(bad), "bigbed")
                    self.assertEqual(
                        str(cm.exception),
                        f"The {name} at offset {offset}, of {length} bytes, "
                        f"extends past the end of the file ({size} bytes)",
                    )


class SeekGuard(BytesIO):
    """A binary stream that fails the test if seeked on and on."""

    def __init__(self, data):
        """Hold the data, and count the seeks made in it."""
        super().__init__(data)
        self.seeks = 0

    def seek(self, *args):
        # Reading a valid bigBed file seeks fewer than 40 times.
        self.seeks += 1
        if self.seeks > 1000:
            raise AssertionError("seek called 1000 times")
        return super().seek(*args)


class TestAlign_malformed_chromosome_index(unittest.TestCase):
    """Chromosome indices whose nodes do not form a tree."""

    def make_stream(self, nodes):
        """Return bed3.bb with its chromosome index replaced, and node offsets.

        The new index is appended to the file. Each node is given as a list
        of the indices of its child nodes in nodes, the first being the root.
        A node without children is written as a leaf holding no chromosomes.
        """
        with open(support.DATA / "Blat" / "bed3.bb", "rb") as stream:
            data = bytearray(stream.read())
        # The offset of the chromosome index is at byte 8 of the file header.
        (index_offset,) = struct.unpack_from("<Q", data, 8)
        # The chromosome index header is 32 bytes; its root node follows it.
        index_header = data[index_offset : index_offset + 32]
        (keySize,) = struct.unpack_from("<I", index_header, 8)
        offsets = []
        offset = len(data) + len(index_header)
        for children in nodes:
            offsets.append(offset)
            offset += 4 + len(children) * (keySize + 8)
        struct.pack_into("<Q", data, 8, len(data))
        data += index_header
        for children in nodes:
            if children:
                data += struct.pack("<?xH", False, len(children))
                for child in children:
                    data += struct.pack(f"<{keySize}sQ", b"chr22", offsets[child])
            else:
                data += struct.pack("<?xH", True, 0)
        return SeekGuard(bytes(data)), offsets

    def check(self, nodes, node):
        stream, offsets = self.make_stream(nodes)
        with self.assertRaises(ValueError) as cm:
            Align.parse(stream, "bigbed")
        self.assertEqual(
            str(cm.exception),
            "Chromosome index is not a tree: the node at offset "
            f"{offsets[node]} is reached more than once",
        )

    def test_unchanged(self):
        # The same index as in bed3.bb, but appended, reads as before.
        with open(support.DATA / "Blat" / "bed3.bb", "rb") as stream:
            data = bytearray(stream.read())
        (index_offset,) = struct.unpack_from("<Q", data, 8)
        struct.pack_into("<Q", data, 8, len(data))
        data += data[index_offset : index_offset + 32 + 4 + 5 + 8]
        with Align.parse(SeekGuard(bytes(data)), "bigbed") as alignments:
            self.assertEqual(len(alignments.targets), 1)
            self.assertEqual(alignments.targets[0].id, "chr22")
            self.assertEqual(len(alignments.targets[0]), 51304566)
            self.assertEqual(len(list(alignments)), 2)

    def test_root_child_of_itself(self):
        # This used to loop for ever, using more memory on each pass.
        self.check([[0]], 0)

    def test_cycle(self):
        self.check([[1], [2], [1]], 1)

    def test_shared_child(self):
        # Reaching each node twice from its parent used to double the work
        # at each level, taking 2**30 steps here.
        nodes = [[i + 1, i + 1] for i in range(30)] + [[]]
        self.check(nodes, 30)

    def test_overlapping_nodes(self):
        # A chain of 100 non-leaf nodes of 100 children each, each node
        # starting inside the one before it, just after that node's first
        # child, which names it. The nodes lie within 3 kB of the file, but
        # each used to be read in full, 130 kB in all.
        with open(support.DATA / "Blat" / "bed3.bb", "rb") as stream:
            data = bytearray(stream.read())
        (index_offset,) = struct.unpack_from("<Q", data, 8)
        index_header = data[index_offset : index_offset + 32]
        (keySize,) = struct.unpack_from("<I", index_header, 8)
        struct.pack_into("<Q", data, 8, len(data))
        data += index_header
        start = len(data)
        stride = 4 + keySize + 8
        for i in range(100):
            data += struct.pack("<?xH", False, 100)
            data += struct.pack(f"<{keySize}sQ", b"chr22", start + stride * (i + 1))
        data += bytes(100 * (keySize + 8))
        with self.assertRaises(ValueError) as cm:
            Align.parse(SeekGuard(bytes(data)), "bigbed")
        self.assertEqual(
            str(cm.exception),
            "Chromosome index is not a tree: its nodes overlap, as together "
            f"they are larger than the file ({len(data)} bytes)",
        )

    def test_past_end_of_file(self):
        # The root's one child is moved to the end of the file, or to the
        # largest offset the index can hold, or the leaf, the last node in
        # the file, claims one entry more than it holds. Reading a node or
        # its entries past the end used to raise struct.error, and seeking
        # to the largest offset OverflowError.
        stream, offsets = self.make_stream([[1], []])
        data = bytearray(stream.getvalue())
        (keySize,) = struct.unpack_from("<I", data, offsets[0] - 32 + 8)
        size = len(data)
        cases = []
        for offset in (size, 2**64 - 1):
            moved = bytearray(data)
            struct.pack_into("<Q", moved, offsets[0] + 4 + keySize, offset)
            cases.append((moved, offset))
        grown = bytearray(data)
        struct.pack_into("<H", grown, offsets[1] + 2, 1)
        cases.append((grown, offsets[1]))
        for bad, offset in cases:
            with self.subTest(offset=offset):
                with self.assertRaises(ValueError) as cm:
                    Align.parse(SeekGuard(bytes(bad)), "bigbed")
                self.assertEqual(
                    str(cm.exception),
                    f"The node at offset {offset} in the chromosome index "
                    f"extends past the end of the file ({size} bytes)",
                )

    def test_wide_node(self):
        # A synthetic index (accepted by the parser; UCSC writes no empty
        # leaves) whose root has 65535 children, the most a node can have:
        # 65534 empty leaves, then the leaf of bed3.bb. Each child used
        # to be taken with list.pop(0), which shifts the children left after
        # it, about 2 billion moves in all, so time was quadratic in the
        # width. Count those moves with a profile hook instead of timing.
        width = 65535
        with open(support.DATA / "Blat" / "bed3.bb", "rb") as stream:
            data = stream.read()
        (index_offset,) = struct.unpack_from("<Q", data, 8)
        leaf = data[index_offset + 32 : index_offset + 32 + 4 + 5 + 8]
        stream, offsets = self.make_stream([list(range(1, width + 1))] + [[]] * width)
        # The last node, an empty leaf, ends the file; swap in bed3.bb's leaf.
        # Each node is a seek, so use a plain stream, not a SeekGuard. The
        # block size in the index header (byte 4 of its 32 bytes) is the
        # most children a node may have; bed3.bb's is 1, so raise it.
        data = bytearray(stream.getvalue()[: offsets[-1]] + leaf)
        struct.pack_into("<I", data, offsets[0] - 32 + 4, width)
        stream = BytesIO(bytes(data))
        moves = 0
        popping = None

        def count_moves(frame, event, arg):
            nonlocal moves, popping
            if getattr(arg, "__name__", None) != "pop":
                return
            items = getattr(arg, "__self__", None)
            if not isinstance(items, list):
                return
            if event == "c_call":
                popping = items[-1] if items else None
            elif event == "c_return" and items and items[-1] == popping:
                # Not taken from the end, so the items after it moved left.
                moves += len(items)

        profile = sys.getprofile()
        sys.setprofile(count_moves)
        try:
            alignments = Align.parse(stream, "bigbed")
        finally:
            sys.setprofile(profile)
        with alignments:
            self.assertEqual(len(alignments.targets), 1)
            self.assertEqual(alignments.targets[0].id, "chr22")
            self.assertEqual(len(list(alignments)), 2)
        self.assertLess(moves, width)


class TestAlign_malformed_data_index(unittest.TestCase):
    """R tree indices of the data whose nodes do not form a tree."""

    not_a_tree = (
        "R tree index is not a tree: the node at offset {} is reached more than once"
    )
    no_children = "The non-leaf node at offset {} in the R tree index has no children"

    def make_stream(self, nodes):
        """Return bed3.bb with its R tree index replaced, and node offsets.

        The new index is appended to the file. Each node is given as a list
        of the indices of its child nodes in nodes, the first being the root,
        or as None for a leaf holding the one data block of bed3.bb.
        """
        with open(support.DATA / "Blat" / "bed3.bb", "rb") as stream:
            data = bytearray(stream.read())
        # The offset of the R tree index is at byte 24 of the file header.
        (index_offset,) = struct.unpack_from("<Q", data, 24)
        # The R tree index header is 48 bytes; its root node follows it. In
        # bed3.bb the root is a leaf holding one data block in 32 bytes.
        index_header = data[index_offset : index_offset + 48]
        leaf = data[index_offset + 48 : index_offset + 48 + 4 + 32]
        # The item count in the header is the number of data blocks, but at
        # least 1, as an index with no items has no nodes to read.
        struct.pack_into("<Q", index_header, 8, max(nodes.count(None), 1))
        offsets = []
        offset = len(data) + len(index_header)
        for children in nodes:
            offsets.append(offset)
            if children is None:
                offset += len(leaf)
            else:
                offset += 4 + len(children) * 24
        struct.pack_into("<Q", data, 24, len(data))
        data += index_header
        for children in nodes:
            if children is None:
                data += leaf
            else:
                data += struct.pack("<?xH", False, len(children))
                for child in children:
                    # The chromosome range of the data block, then the offset.
                    data += leaf[4:20] + struct.pack("<Q", offsets[child])
        return SeekGuard(bytes(data)), offsets

    def check(self, nodes, node, message):
        stream, offsets = self.make_stream(nodes)
        with self.assertRaises(ValueError) as cm:
            Align.parse(stream, "bigbed")
        self.assertEqual(str(cm.exception), message.format(offsets[node]))

    def test_tree(self):
        # A tree of three levels with three leaves, each holding the data
        # block of bed3.bb, so its two alignments are read three times.
        stream, _ = self.make_stream([[1, 2], [3, 4], [5], None, None, None])
        with Align.parse(stream, "bigbed") as alignments:
            self.assertEqual(len(alignments.targets), 1)
            self.assertEqual(alignments.targets[0].id, "chr22")
            self.assertEqual(len(list(alignments)), 6)

    def test_root_child_of_itself(self):
        # This used to loop for ever, using more memory on each pass.
        self.check([[0]], 0, self.not_a_tree)

    def test_cycle(self):
        self.check([[1], [2], [1]], 1, self.not_a_tree)

    def test_shared_child(self):
        # Reaching each node twice from its parent used to double the work
        # at each level, taking 2**30 steps here.
        nodes = [[i + 1, i + 1] for i in range(30)] + [None]
        self.check(nodes, 30, self.not_a_tree)

    def test_no_children(self):
        # These used to raise IndexError.
        self.check([[]], 0, self.no_children)
        self.check([[1], []], 1, self.no_children)

    def test_empty_leaf(self):
        # A leaf with no items may only be the root of a tree with no data
        # blocks. Both of these used to raise IndexError.
        for nodes in ([None], [[1], None]):
            stream, offsets = self.make_stream(nodes)
            data = bytearray(stream.getvalue())
            # Set the count of the leaf (byte 2 of the node) to 0, and for a
            # lone root also the item count in the index header (48 bytes
            # before the root) and the item count of the file, which starts
            # the data section (its offset is at byte 16 of the file header).
            struct.pack_into("<H", data, offsets[-1] + 2, 0)
            if len(nodes) == 1:
                struct.pack_into("<Q", data, offsets[0] - 48 + 8, 0)
                (data_offset,) = struct.unpack_from("<Q", data, 16)
                struct.pack_into("<Q", data, data_offset, 0)
            stream = SeekGuard(bytes(data))
            with self.subTest(nodes=nodes):
                if len(nodes) == 1:
                    with Align.parse(stream, "bigbed") as alignments:
                        self.assertEqual(len(list(alignments)), 0)
                        self.assertEqual(len(list(alignments.search("chr22"))), 0)
                else:
                    with self.assertRaises(ValueError) as cm:
                        Align.parse(stream, "bigbed")
                    self.assertEqual(
                        str(cm.exception),
                        f"The leaf node at offset {offsets[1]} in the R tree "
                        "index has no items",
                    )

    def test_overlapping_nodes(self):
        # A root naming 100 leaves that start 4 bytes apart, so that they
        # overlap, each holding 100 data blocks. The leaves lie within 3.6 kB
        # of the file, but each used to be read in full, giving 10000 data
        # blocks, 320 kB in all; the item count in the header agrees.
        with open(support.DATA / "Blat" / "bed3.bb", "rb") as stream:
            data = bytearray(stream.read())
        (index_offset,) = struct.unpack_from("<Q", data, 24)
        index_header = data[index_offset : index_offset + 48]
        struct.pack_into("<Q", index_header, 8, 100 * 100)
        struct.pack_into("<Q", data, 24, len(data))
        data += index_header
        leaves = len(data) + 4 + 100 * 24
        data += struct.pack("<?xH", False, 100)
        for i in range(100):
            data += struct.pack("<IIIIQ", 0, 0, 0, 0, leaves + 4 * i)
        # Every 4 bytes from here is the header of a leaf of 100 entries of
        # 32 bytes each, and the bytes are long enough for the last leaf.
        data += struct.pack("<?xH", True, 100) * (100 + 8 * 100)
        with self.assertRaises(ValueError) as cm:
            Align.parse(SeekGuard(bytes(data)), "bigbed")
        self.assertEqual(
            str(cm.exception),
            "R tree index is not a tree: its nodes overlap, as together "
            f"they are larger than the file ({len(data)} bytes)",
        )

    def test_past_end_of_file(self):
        # The root's one child is moved to the end of the file, or to the
        # largest offset the index can hold, or the leaf, the last node in
        # the file, claims one entry more than it holds. Reading a node or
        # its entries past the end used to raise struct.error, and seeking
        # to the largest offset OverflowError.
        stream, offsets = self.make_stream([[1], None])
        data = bytearray(stream.getvalue())
        size = len(data)
        cases = []
        for offset in (size, 2**64 - 1):
            moved = bytearray(data)
            struct.pack_into("<Q", moved, offsets[0] + 4 + 16, offset)
            cases.append((moved, offset))
        grown = bytearray(data)
        struct.pack_into("<H", grown, offsets[1] + 2, 2)
        cases.append((grown, offsets[1]))
        for bad, offset in cases:
            with self.subTest(offset=offset):
                with self.assertRaises(ValueError) as cm:
                    Align.parse(SeekGuard(bytes(bad)), "bigbed")
                self.assertEqual(
                    str(cm.exception),
                    f"The node at offset {offset} in the R tree index "
                    f"extends past the end of the file ({size} bytes)",
                )

    def test_wide_node(self):
        # A valid tree whose root has 65535 leaves, the most a node can
        # have, each holding the data block of bed3.bb. Moving on to the
        # next sibling used to search the parent's children from the start,
        # so opening, iterating and searching each took time quadratic in
        # the width, about 5 minutes each. All three now take a few seconds.
        # Run them in a child process, which a timeout can stop.
        width = 65535
        stream, _ = self.make_stream([list(range(1, width + 1))] + [None] * width)
        code = """\
import sys
from Bio import Align
with Align.parse(sys.argv[1], "bigbed") as alignments:
    print(sum(1 for alignment in alignments))
    print(sum(1 for alignment in alignments.search("chr22")))
"""
        env = os.environ.copy()
        env["PYTHONPATH"] = os.pathsep.join(
            filter(
                None,
                [os.path.dirname(os.path.dirname(Bio.__file__)), env.get("PYTHONPATH")],
            )
        )
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "wide_node.bb")
            with open(path, "wb") as output:
                output.write(stream.getvalue())
            try:
                result = subprocess.run(
                    [sys.executable, "-c", code, path],
                    timeout=60,
                    capture_output=True,
                    text=True,
                    env=env,
                )
            except subprocess.TimeoutExpired:
                self.fail("Reading a root node of 65535 leaves took over 60 seconds")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, f"{2 * width}\n{2 * width}\n")


class TestAlign_malformed_data_block(unittest.TestCase):
    """Compressed data blocks larger than uncompressBufSize, or incomplete."""

    # The one data block of bed3.bb, uncompressed. Each of its two items is
    # the chromosome number, start and end, then the other fields, of which
    # bed3 has none, ended by a NUL byte.
    items = (
        struct.pack("<III", 0, 1000, 5000)
        + b"\0"
        + struct.pack("<III", 0, 2000, 6000)
        + b"\0"
    )

    too_large = (
        "The data block at offset {} decompresses to more than "
        "uncompressBufSize ({}) bytes"
    )
    incomplete = "The compressed data block at offset {} is incomplete"
    invalid = "The compressed data block at offset {} is invalid ({})"
    past_end = (
        "The data block at offset {}, of {} bytes, extends past the end of "
        "the file ({} bytes)"
    )

    def make_stream(self, block, uncompressBufSize=None, offset=None, size=None):
        """Return bed3.bb with its data block replaced, and the block offset.

        The new block, given compressed, is appended to the file, and the one
        leaf of the R tree index points to it instead. If uncompressBufSize is
        given, it replaces the value of 16384 in the file header. If offset or
        size is given, the leaf gives it in place of the offset or size of the
        new block.
        """
        with open(support.DATA / "Blat" / "bed3.bb", "rb") as stream:
            data = bytearray(stream.read())
        if uncompressBufSize is not None:
            # uncompressBufSize is at byte 52 of the file header.
            struct.pack_into("<I", data, 52, uncompressBufSize)
        # The offset of the R tree index is at byte 24 of the file header.
        (index_offset,) = struct.unpack_from("<Q", data, 24)
        # The R tree index header is 48 bytes; its root node follows it. In
        # bed3.bb the root is a leaf holding one data block, and the offset
        # and size of the block are the last 16 of its 4 + 32 bytes.
        if offset is None:
            offset = len(data)
        if size is None:
            size = len(block)
        struct.pack_into("<QQ", data, index_offset + 48 + 20, offset, size)
        data += block
        return BytesIO(bytes(data)), offset

    def check(self, block, uncompressBufSize, message):
        stream, offset = self.make_stream(block, uncompressBufSize)
        if uncompressBufSize is None:
            uncompressBufSize = 16384
        message = message.format(offset, uncompressBufSize)
        self.check_message(stream, message)

    def check_message(self, stream, message):
        # Iterating and searching read the data block alike.
        with Align.parse(stream, "bigbed") as alignments:
            with self.assertRaises(ValueError) as cm:
                next(alignments)
            self.assertEqual(str(cm.exception), message)
            with self.assertRaises(ValueError) as cm:
                next(alignments.search("chr22"))
            self.assertEqual(str(cm.exception), message)

    def test_fits(self):
        # A block of exactly uncompressBufSize bytes reads as before.
        stream, _ = self.make_stream(zlib.compress(self.items), len(self.items))
        with Align.parse(stream, "bigbed") as alignments:
            coordinates = [alignment.coordinates.tolist() for alignment in alignments]
            self.assertEqual(
                coordinates, [[[1000, 5000], [0, 4000]], [[2000, 6000], [0, 4000]]]
            )
            coordinates = [
                alignment.coordinates.tolist()
                for alignment in alignments.search("chr22", 4500, 5500)
            ]
            self.assertEqual(
                coordinates, [[[1000, 5000], [0, 4000]], [[2000, 6000], [0, 4000]]]
            )

    def test_chromId_not_in_index(self):
        # bed3.bb has one chromosome, chr22, so its only chromId is 0. Looking
        # up chromId 7 in the chromosome list used to raise IndexError.
        items = struct.pack("<III", 7, 1000, 5000) + self.items[12:]
        stream, _ = self.make_stream(zlib.compress(items))
        message = (
            "The item at 1000-5000 has chromId 7, which is not in the "
            "chromosome index (number of chromosomes: 1)"
        )
        with Align.parse(stream, "bigbed") as alignments:
            with self.assertRaises(ValueError) as cm:
                next(alignments)
            self.assertEqual(str(cm.exception), message)
            # Searching all chromosomes reads every item, as iterating does.
            with self.assertRaises(ValueError) as cm:
                next(alignments.search())
            self.assertEqual(str(cm.exception), message)
            # Searching one chromosome skips items on others, as before.
            coordinates = [
                alignment.coordinates.tolist()
                for alignment in alignments.search("chr22")
            ]
            self.assertEqual(coordinates, [[[2000, 6000], [0, 4000]]])

    def test_one_byte_too_large(self):
        block = zlib.compress(self.items)
        self.check(block, len(self.items) - 1, self.too_large)

    def test_decompression_bomb(self):
        # 1.3 MB of NUL bytes, which zlib compresses to 1.3 kB, are 100000
        # empty items. Blocks used to be decompressed whatever their size, so
        # a file of 274 kB used over 500 MB of memory on the first next().
        block = zlib.compress(bytes(13 * 100000))
        self.check(block, None, self.too_large)

    def test_decompression_bomb_large_header(self):
        # The header's uncompressBufSize comes from the file, so a crafted
        # file could raise it to 4 GB. The reader now stops at 256 MB anyway.
        self.assertEqual(bigbed._MAX_UNCOMPRESS_BUF_SIZE, 256 * 1024 * 1024)
        # Lower that cap to 1 kB here, rather than decompress 256 MB. The
        # block is the two items followed by NUL bytes, 1 kB + 1 bytes in
        # all; it used to decompress fully and the items read as normal.
        limit = 1024
        block = zlib.compress(self.items + bytes(limit + 1 - len(self.items)))
        stream, offset = self.make_stream(block, 0xFFFFFFFF)
        # The message used to give the limit as the header's uncompressBufSize,
        # so a valid but very large file looked corrupt.
        message = (
            f"The data block at offset {offset} decompresses to more than "
            f"{limit} bytes, the largest block Biopython reads (the file "
            "header gives uncompressBufSize = 4294967295)"
        )
        with mock.patch.object(bigbed, "_MAX_UNCOMPRESS_BUF_SIZE", limit):
            self.check_message(stream, message)

    def test_incomplete(self):
        # These used to raise zlib.error.
        block = zlib.compress(self.items)
        self.check(block[:-1], None, self.incomplete)
        self.check(b"", None, self.incomplete)

    def test_truncated_item(self):
        # The block ends with 5 of the 12 bytes of the chromosome number,
        # start and end of its second item. This used to raise struct.error.
        stream, offset = self.make_stream(zlib.compress(self.items[:18]))
        message = f"The data block at offset {offset} ends with a truncated item"
        with Align.parse(stream, "bigbed") as alignments:
            with self.assertRaises(ValueError) as cm:
                list(alignments)
            self.assertEqual(str(cm.exception), message)
            with self.assertRaises(ValueError) as cm:
                list(alignments.search("chr22"))
            self.assertEqual(str(cm.exception), message)
        # With one item per slot, a search reads only the first item of each
        # block; here the block holds 5 of the 12 bytes of that item, then
        # all 12 bytes but not the NUL byte that ends it. The last byte of
        # its end (5000) is 0, which was read as that NUL byte, so the
        # 12-byte block used to give a bed3 alignment.
        for length in (5, 12):
            with self.subTest(length=length):
                block = zlib.compress(self.items[:length])
                stream, offset = self.make_stream(block)
                message = (
                    f"The data block at offset {offset} ends with a truncated item"
                )
                with Align.parse(stream, "bigbed") as alignments:
                    alignments.itemsPerSlot = 1
                    with self.assertRaises(ValueError) as cm:
                        list(alignments.search("chr22"))
                    self.assertEqual(str(cm.exception), message)

    def test_item_without_nul(self):
        # The block ends with all 12 bytes of the chromosome number, start and
        # end of its second item, but not the NUL byte that ends the item.
        # This used to raise ValueError("subsection not found").
        stream, offset = self.make_stream(zlib.compress(self.items[:-1]))
        message = f"The data block at offset {offset} ends with a truncated item"
        with Align.parse(stream, "bigbed") as alignments:
            with self.assertRaises(ValueError) as cm:
                list(alignments)
            self.assertEqual(str(cm.exception), message)
            with self.assertRaises(ValueError) as cm:
                list(alignments.search("chr22"))
            self.assertEqual(str(cm.exception), message)

    def test_invalid(self):
        # These used to raise zlib.error. The last four bytes of the
        # compressed data are the Adler-32 checksum of the block.
        block = bytearray(zlib.compress(self.items))
        block[-1] ^= 1
        stream, offset = self.make_stream(bytes(block))
        error = "Error -3 while decompressing data: incorrect data check"
        self.check_message(stream, self.invalid.format(offset, error))
        stream, offset = self.make_stream(b"not zlib" + zlib.compress(self.items))
        error = "Error -3 while decompressing data: incorrect header check"
        self.check_message(stream, self.invalid.format(offset, error))

    def test_past_end_of_file(self):
        block = zlib.compress(self.items)
        # The file ends with the block, which used to read as it holds all
        # of the compressed data, though the index gives one more byte.
        stream, offset = self.make_stream(block, size=len(block) + 1)
        file_size = offset + len(block)
        message = self.past_end.format(offset, len(block) + 1, file_size)
        self.check_message(stream, message)
        # This used to raise OverflowError, or MemoryError for a file on disk.
        stream, offset = self.make_stream(block, size=2**64 - 1)
        message = self.past_end.format(offset, 2**64 - 1, file_size)
        self.check_message(stream, message)
        stream, _ = self.make_stream(block, offset=file_size + 1)
        message = self.past_end.format(file_size + 1, len(block), file_size)
        self.check_message(stream, message)

    def test_short_read(self):
        # The file loses the end of its block after it was opened. The
        # compressed data are complete, and the short read used to pass.
        block = zlib.compress(self.items) + bytes(8)
        stream, offset = self.make_stream(block)
        with Align.parse(stream, "bigbed") as alignments:
            stream.truncate(offset + len(block) - 8)
            message = (
                f"The data block at offset {offset} is incomplete: read "
                f"{len(block) - 8} of {len(block)} bytes"
            )
            with self.assertRaises(ValueError) as cm:
                next(alignments)
            self.assertEqual(str(cm.exception), message)
            with self.assertRaises(ValueError) as cm:
                next(alignments.search("chr22"))
            self.assertEqual(str(cm.exception), message)

    def test_large_block(self):
        # A valid block of 3000000 empty items, 39 MB in all, that fits in
        # uncompressBufSize, in a file of 51 kB. Each item read used to copy
        # the rest of the block, so reading the first 100000 items copied
        # about 4 TB, and 2000 of them took 57 seconds. They now take about
        # a second. Read them in a child process, which a timeout can stop.
        data = bytes(13 * 3000000)
        stream, _ = self.make_stream(zlib.compress(data), len(data))
        code = """\
import itertools
import sys
from Bio import Align
count = 0
with Align.parse(sys.argv[1], "bigbed") as alignments:
    for alignment in itertools.islice(alignments, 100000):
        assert alignment.coordinates.tolist() == [[0, 0], [0, 0]]
        count += 1
print(count)
"""
        env = os.environ.copy()
        env["PYTHONPATH"] = os.pathsep.join(
            filter(
                None,
                [os.path.dirname(os.path.dirname(Bio.__file__)), env.get("PYTHONPATH")],
            )
        )
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "large_block.bb")
            with open(path, "wb") as output:
                output.write(stream.getvalue())
            try:
                result = subprocess.run(
                    [sys.executable, "-c", code, path],
                    timeout=60,
                    capture_output=True,
                    text=True,
                    env=env,
                )
            except subprocess.TimeoutExpired:
                self.fail("Reading 100000 items of the block took over 60 seconds")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "100000\n")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
