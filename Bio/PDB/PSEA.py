# Copyright (C) 2006, Thomas Hamelryck (thamelry@binf.ku.dk)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Wrappers for PSEA, a program for secondary structure assignment.

See this citation for P-SEA, PMID: 9183534

Labesse G, Colloc'h N, Pothier J, Mornon J-P:  P-SEA: a new efficient
assignment of secondary structure from C_alpha.
Comput Appl Biosci 1997 , 13:291-295

ftp://ftp.lmcp.jussieu.fr/pub/sincris/software/protein/p-sea/
"""

import contextlib
import errno
import os
import shutil
import subprocess
import tempfile

from Bio.PDB.Polypeptide import is_aa


def _move(src, dst):
    """Move file src to exactly dst, replacing it, also across filesystems."""
    try:
        os.replace(src, dst)
    except OSError as err:
        if err.errno != errno.EXDEV:
            raise
        # Copy next to dst first, so that dst is replaced in one step.
        fd, tmp = tempfile.mkstemp(dir=os.path.dirname(dst))
        os.close(fd)
        try:
            shutil.copy2(src, tmp)
            os.replace(tmp, dst)
        except BaseException:
            # tmp is gone if the interruption came after it was renamed.
            with contextlib.suppress(FileNotFoundError):
                os.remove(tmp)
            raise


def _join_cwd(cwd, fname):
    """Return fname as an absolute path, taking a relative one from cwd."""
    # Not os.path.abspath, which would resolve "link/.." without following
    # the link, unlike psea opening the file.
    path = os.path.join(cwd, fname)
    if not os.path.isabs(path):
        # A Windows path on another drive, such as "D:tiny.pdb", is still
        # relative, to the current directory of that drive.
        path = os.path.abspath(path)
    return path


def run_psea(fname, verbose=False):
    """Run PSEA and return output filename."""
    # Once, so that the input and output paths stay in the directory of the
    # call even if the working directory changes while psea runs.
    cwd = os.getcwd()
    fname = _join_cwd(cwd, fname)
    last = os.path.basename(fname)
    base = os.path.splitext(last)[0]
    cmd = ["psea", fname]

    with tempfile.TemporaryDirectory() as tmpdir:
        p = subprocess.run(cmd, cwd=tmpdir, capture_output=True, text=True)

        if verbose:
            print(p.stdout)

        output = os.path.join(tmpdir, base + ".sea")

        if not p.stderr.strip() and os.path.exists(output):
            # move output to the current directory
            final_path = os.path.join(cwd, base + ".sea")
            _move(output, final_path)
            return final_path
        else:
            raise RuntimeError(f"Error running p-sea: {p.stderr}")


def psea(pname):
    """Parse PSEA output file."""
    fname = run_psea(pname)
    start = 0
    ss = ""
    with open(fname) as fp:
        for line in fp:
            if line[0:6] == ">p-sea":
                start = 1
                continue
            if not start:
                continue
            if line[0] == "\n":
                break
            ss = ss + line[0:-1]
    return ss


def psea2HEC(pseq):
    """Translate PSEA secondary structure string into HEC."""
    seq = []
    for ss in pseq:
        if ss == "a":
            n = "H"
        elif ss == "b":
            n = "E"
        elif ss == "c":
            n = "C"
        seq.append(n)
    return seq


def annotate(m, ss_seq):
    """Apply secondary structure information to residues in model."""
    c = m.get_list()[0]
    all = c.get_list()
    residues = []
    # Now remove HOH etc.
    for res in all:
        if is_aa(res):
            residues.append(res)
    L = len(residues)
    if not L == len(ss_seq):
        raise ValueError("Length mismatch %i %i" % (L, len(ss_seq)))
    for i in range(L):
        residues[i].xtra["SS_PSEA"] = ss_seq[i]
    # subprocess.call(["rm", fname])


class PSEA:
    """Define PSEA class.

    PSEA object is a wrapper to PSEA program for secondary structure assignment.
    """

    def __init__(self, model, filename):
        """Initialize the class."""
        ss_seq = psea(filename)
        ss_seq = psea2HEC(ss_seq)
        annotate(model, ss_seq)
        self.ss_seq = ss_seq

    def get_seq(self):
        """Return secondary structure string."""
        return self.ss_seq
