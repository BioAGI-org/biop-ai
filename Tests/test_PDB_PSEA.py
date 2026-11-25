# Copyright 2022 by Valentin Vareskic (valentin.vareskic@gmail.com).
# All rights reserved. This code is part of the Biopython distribution
# and governed by its license.  Please see the LICENSE file that should
# have been included as part of this package.

"""Tests for PDB PSEA."""

import io
import os
import sys
import tempfile
import unittest
from subprocess import getoutput

import support

from Bio import MissingExternalDependencyError
from Bio.PDB import PDBParser
from Bio.PDB.PSEA import PSEA
from Bio.PDB.PSEA import psea
from Bio.PDB.PSEA import psea2HEC
from Bio.PDB.PSEA import run_psea

os.environ["LANG"] = "C"


cmd_output = getoutput("psea -h")
if not cmd_output.startswith("o---"):
    raise MissingExternalDependencyError(
        "Download and install psea from "
        "ftp://ftp.lmcp.jussieu.fr/pub/sincris/software/protein/p-sea/. "
        "Make sure that psea is on path"
    )


def run_in_temp_dir(test):
    """Run the test in a new temporary directory, removed afterwards.

    run_psea writes its output to the working directory, so this keeps the
    tests away from any .sea files already in the directory they start in.
    """
    tmpdir = tempfile.TemporaryDirectory()
    # Cleanups run last first: change back, then remove the directory.
    test.addCleanup(tmpdir.cleanup)
    test.addCleanup(os.chdir, os.getcwd())
    os.chdir(tmpdir.name)


class TestPDBPSEA(unittest.TestCase):
    def setUp(self):
        run_in_temp_dir(self)

    def test_run_psea_verbose(self):
        captured_ouput = io.StringIO()
        sys.stdout = captured_ouput
        psae_run = run_psea(support.DATA / "PDB" / "1A8O.pdb", verbose=True)
        sys.stdout = sys.__stdout__
        self.assertEqual(os.path.basename(psae_run), "1A8O.sea")
        self.assertTrue(captured_ouput.getvalue())

    def test_run_psea_quiet(self):
        captured_ouput = io.StringIO()
        sys.stdout = captured_ouput
        psae_run = run_psea(support.DATA / "PDB" / "1A8O.pdb", verbose=False)
        sys.stdout = sys.__stdout__
        self.assertEqual(os.path.basename(psae_run), "1A8O.sea")
        self.assertFalse(captured_ouput.getvalue())

    def test_psea(self):
        psae_run = psea(support.DATA / "PDB" / "2BEG.pdb")
        self.assertEqual(psae_run, "ccccbbbbbbbccccbbbbbbbbbbc")

    def test_psea_2HEC(self):
        seq = psea(support.DATA / "PDB" / "2BEG.pdb")
        psae_run = psea2HEC(seq)
        self.assertEqual(
            psae_run,
            [
                "C",
                "C",
                "C",
                "C",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "C",
                "C",
                "C",
                "C",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "C",
            ],
        )

    def test_run_psea_tempdir(self):
        psae_run = run_psea(support.DATA / "PDB" / "1A8O.pdb")
        self.assertTrue(os.path.exists(psae_run))
        self.assertTrue(psae_run.endswith(".sea"))

    def test_keeps_existing_sea_files(self):
        """Running a test leaves the .sea files in the directory alone."""
        with open("unrelated.sea", "w") as handle:
            handle.write("keep\n")
        with open("1A8O.sea", "w") as handle:
            handle.write("mine\n")
        result = unittest.TestResult()
        TestPDBPSEA("test_run_psea_quiet").run(result)
        self.assertTrue(result.wasSuccessful(), result.failures + result.errors)
        with open("unrelated.sea") as handle:
            self.assertEqual(handle.read(), "keep\n")
        with open("1A8O.sea") as handle:
            self.assertEqual(handle.read(), "mine\n")


class TestPSEA(unittest.TestCase):
    def setUp(self):
        run_in_temp_dir(self)

    def test_get_seq(self):
        p = PDBParser()
        s = p.get_structure("X", support.DATA / "PDB" / "2BEG.pdb")
        psea_class = PSEA(s[0], support.DATA / "PDB" / "2BEG.pdb")
        self.assertEqual(
            psea_class.get_seq(),
            [
                "C",
                "C",
                "C",
                "C",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "C",
                "C",
                "C",
                "C",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "E",
                "C",
            ],
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
