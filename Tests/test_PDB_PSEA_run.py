# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for Bio.PDB.PSEA.run_psea with a stand-in for the psea program.

These run without psea installed, unlike test_PDB_PSEA.py.
"""

import errno
import ntpath
import os
import subprocess
import sys
import tempfile
import types
import unittest
from unittest import mock

try:
    import numpy as np  # noqa: F401
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB."
    ) from None

from Bio.PDB.PSEA import _join_cwd
from Bio.PDB.PSEA import run_psea

# Like psea: read the input file and write <name>.sea in the working
# directory. The output is a copy of the input, to show which file was read.
FAKE_PSEA = """\
import os, sys
with open(sys.argv[1]) as handle:
    data = handle.read()
base = os.path.splitext(os.path.basename(sys.argv[1]))[0]
with open(base + ".sea", "w") as handle:
    handle.write(data)
"""

real_run = subprocess.run

# The working directory of this process each time fake_run was called
run_cwds = []


def fake_run(cmd, **kwargs):
    run_cwds.append(os.getcwd())
    return real_run([sys.executable, "-c", FAKE_PSEA, *cmd[1:]], **kwargs)


class Interrupted(BaseException):
    """Stand-in for KeyboardInterrupt that pytest will not act on."""


def missing_run(cmd, **kwargs):
    raise FileNotFoundError(2, "No such file or directory", cmd[0])


class TestRunPSEA(unittest.TestCase):
    def setUp(self):
        run_cwds.clear()
        self.olddir = os.getcwd()
        self.tmpdir = tempfile.TemporaryDirectory()
        os.chdir(self.tmpdir.name)
        # As the system reports it, which on Windows can be a short (8.3)
        # name that os.path.realpath would expand.
        self.cwd = os.getcwd()
        with open("tiny.pdb", "w") as handle:
            handle.write("END\n")

    def tearDown(self):
        os.chdir(self.olddir)
        self.tmpdir.cleanup()

    def read(self, path):
        with open(path) as handle:
            return handle.read()

    def test_relative_path_and_working_directory(self):
        """Accept a relative path, and never change the working directory."""
        with mock.patch("Bio.PDB.PSEA.subprocess.run", fake_run):
            output = run_psea("tiny.pdb")
        self.assertEqual(run_cwds, [self.cwd])
        self.assertEqual(output, os.path.join(self.cwd, "tiny.sea"))
        self.assertEqual(self.read("tiny.sea"), "END\n")
        # psea not installed: the error must not leave us in a deleted directory
        with mock.patch("Bio.PDB.PSEA.subprocess.run", missing_run):
            self.assertRaises(FileNotFoundError, run_psea, "tiny.pdb")
        self.assertEqual(os.getcwd(), self.cwd)

    def test_cross_filesystem_move(self):
        """Move the output into place when renaming fails with EXDEV."""
        real_rename, real_replace = os.rename, os.replace
        moves = []

        def cross_device_once(real):
            # Fail the first rename, as across filesystems, whichever
            # function makes it.
            def move(src, dst, *args, **kwargs):
                moves.append(dst)
                if len(moves) == 1:
                    raise OSError(errno.EXDEV, "Invalid cross-device link", src)
                return real(src, dst, *args, **kwargs)

            return move

        with open("tiny.sea", "w") as handle:
            handle.write("old output\n")
        with (
            mock.patch("Bio.PDB.PSEA.subprocess.run", fake_run),
            mock.patch("os.rename", cross_device_once(real_rename)),
            mock.patch("os.replace", cross_device_once(real_replace)),
        ):
            output = run_psea("tiny.pdb")
        self.assertEqual(moves[0], output)
        self.assertEqual(output, os.path.join(self.cwd, "tiny.sea"))
        self.assertEqual(self.read("tiny.sea"), "END\n")
        # No temporary copy left behind
        self.assertEqual(sorted(os.listdir()), ["tiny.pdb", "tiny.sea"])

    def test_output_name_is_a_directory(self):
        """Fail, rather than write the output inside a directory of its name."""
        os.mkdir("tiny.sea")
        with mock.patch("Bio.PDB.PSEA.subprocess.run", fake_run):
            self.assertRaises(OSError, run_psea, "tiny.pdb")
        self.assertEqual(os.listdir("tiny.sea"), [])

    def test_working_directory_changes_during_run(self):
        """Write the output in the directory of the call, not a later one."""
        os.mkdir("other")
        with open(os.path.join("other", "tiny.sea"), "w") as handle:
            handle.write("unrelated\n")

        def chdir_run(cmd, **kwargs):
            os.chdir("other")
            return fake_run(cmd, **kwargs)

        with mock.patch("Bio.PDB.PSEA.subprocess.run", chdir_run):
            output = run_psea("tiny.pdb")
        self.assertEqual(output, os.path.join(self.cwd, "tiny.sea"))
        self.assertEqual(self.read(output), "END\n")
        self.assertEqual(self.read(os.path.join("..", "tiny.sea")), "END\n")
        self.assertEqual(self.read("tiny.sea"), "unrelated\n")

    def test_interrupt_after_cross_filesystem_move(self):
        """Raise the interruption, even if it comes after the output is moved."""
        real_replace = os.replace
        calls = []

        def replace(src, dst, *args, **kwargs):
            calls.append(dst)
            if len(calls) == 1:
                raise OSError(errno.EXDEV, "Invalid cross-device link", src)
            real_replace(src, dst, *args, **kwargs)
            raise Interrupted

        with (
            mock.patch("Bio.PDB.PSEA.subprocess.run", fake_run),
            mock.patch("os.replace", replace),
        ):
            self.assertRaises(Interrupted, run_psea, "tiny.pdb")
        self.assertEqual(len(calls), 2)
        self.assertEqual(self.read("tiny.sea"), "END\n")
        self.assertEqual(sorted(os.listdir()), ["tiny.pdb", "tiny.sea"])

    def test_windows_path_on_another_drive(self):
        """Make a Windows path absolute even when it names another drive."""
        windows = types.SimpleNamespace(path=ntpath)
        with mock.patch("Bio.PDB.PSEA.os", windows):
            self.assertEqual(
                _join_cwd("C:\\caller", "tiny.pdb"), "C:\\caller\\tiny.pdb"
            )
            self.assertEqual(_join_cwd("C:\\caller", "D:\\tiny.pdb"), "D:\\tiny.pdb")
            # Joining keeps "link\\.." for the system to resolve.
            self.assertEqual(
                _join_cwd("C:\\caller", "link\\..\\tiny.pdb"),
                "C:\\caller\\link\\..\\tiny.pdb",
            )
            # "D:tiny.pdb" is relative to the current directory of drive D:,
            # which only os.path.abspath (GetFullPathNameW) knows.
            with mock.patch.object(
                ntpath, "abspath", return_value="D:\\work\\tiny.pdb"
            ) as abspath:
                path = _join_cwd("C:\\caller", "D:tiny.pdb")
            abspath.assert_called_once_with("D:tiny.pdb")
            self.assertEqual(path, "D:\\work\\tiny.pdb")

    @unittest.skipIf(os.name == "nt", "Windows resolves '..' before following links")
    def test_parent_of_symbolic_link(self):
        """Read the file psea would open for a path with link/.. in it."""
        os.makedirs(os.path.join("real", "sub"))
        with open(os.path.join("real", "tiny.pdb"), "w") as handle:
            handle.write("REAL\n")
        os.mkdir("dir")
        with open(os.path.join("dir", "tiny.pdb"), "w") as handle:
            handle.write("WRONG\n")
        try:
            os.symlink(
                os.path.join(self.cwd, "real", "sub"), os.path.join("dir", "link")
            )
        except (OSError, NotImplementedError) as err:
            self.skipTest(f"cannot make a symbolic link: {err}")
        with mock.patch("Bio.PDB.PSEA.subprocess.run", fake_run):
            output = run_psea(os.path.join("dir", "link", "..", "tiny.pdb"))
        self.assertEqual(self.read(output), "REAL\n")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
