# Copyright 2022 by Valentin Vareskic (valentin.vareskic@gmail.com).
# All rights reserved. This code is part of the Biopython distribution
# and governed by its license.  Please see the LICENSE file that should
# have been included as part of this package.

"""Tests for PDB PSEA."""

import io
import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import support

from Bio.PDB import PDBParser
from Bio.PDB.PSEA import PSEA
from Bio.PDB.PSEA import psea
from Bio.PDB.PSEA import psea2HEC
from Bio.PDB.PSEA import run_psea

os.environ["LANG"] = "C"


# Only the tests that run the psea program need it, so they skip without it.
have_psea = subprocess.getoutput("psea -h").startswith("o---")


def require_psea(test):
    """Skip the test unless the psea program is on the path."""
    if not have_psea:
        test.skipTest(
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
        require_psea(self)
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


def fake_run(cmd, cwd, **kwargs):
    """Stand in for the psea program: write <name>.sea where it runs.

    The output always holds what psea gives for 2BEG.pdb.
    """
    base = os.path.splitext(os.path.basename(cmd[1]))[0]
    with open(os.path.join(cwd, base + ".sea"), "w") as handle:
        handle.write(">p-sea\nccccbbbbbbbccccbbbbbbbbbbc\n\n")
    return subprocess.CompletedProcess(cmd, 0, "", "")


class TestRunInTempDir(unittest.TestCase):
    """run_in_temp_dir, with a stand-in for psea, so it needs no psea."""

    def setUp(self):
        run_in_temp_dir(self)

    def test_keeps_existing_sea_files(self):
        """Running a test leaves the .sea files in the directory alone."""
        with open("unrelated.sea", "w") as handle:
            handle.write("keep\n")
        for name in ("1A8O.sea", "2BEG.sea"):
            with open(name, "w") as handle:
                handle.write("mine\n")
        result = unittest.TestResult()
        # Run a real test of each class that runs psea, even without psea.
        with (
            mock.patch("Bio.PDB.PSEA.subprocess.run", fake_run),
            mock.patch.object(sys.modules[__name__], "have_psea", True),
        ):
            TestPDBPSEA("test_run_psea_tempdir").run(result)
            TestPSEA("test_get_seq").run(result)
        self.assertTrue(result.wasSuccessful(), result.failures + result.errors)
        self.assertEqual(result.testsRun, 2)
        self.assertFalse(result.skipped)
        with open("unrelated.sea") as handle:
            self.assertEqual(handle.read(), "keep\n")
        for name in ("1A8O.sea", "2BEG.sea"):
            with open(name) as handle:
                self.assertEqual(handle.read(), "mine\n", name)


class TestPsea2HEC(unittest.TestCase):
    """psea2HEC, which needs no psea program."""

    def test_valid_symbols(self):
        self.assertEqual(psea2HEC("abc"), ["H", "E", "C"])

    def test_unknown_symbol(self):
        """A symbol other than a, b or c is an error, wherever it is."""
        for pseq in ("Xab", "aXc", "ab c"):
            with self.subTest(pseq=pseq):
                with self.assertRaisesRegex(ValueError, "symbol"):
                    psea2HEC(pseq)


class TestPSEA(unittest.TestCase):
    def setUp(self):
        require_psea(self)
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
