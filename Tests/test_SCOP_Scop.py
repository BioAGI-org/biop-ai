# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# Modifications Copyright 2010 Jeffrey Finkelstein. All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Unit test for Scop."""

import os
import subprocess
import sys
import unittest
from contextlib import redirect_stdout
from io import StringIO
from unittest import mock

import support

import Bio
from Bio.SCOP import cmp_sccs
from Bio.SCOP import Node
from Bio.SCOP import parse_domain
from Bio.SCOP import Scop


class ScopTests(unittest.TestCase):
    def _compare_cla_lines(self, cla_line_1, cla_line_2):
        """Compare the two specified Cla lines for equality.

        The order of the key-value pairs in the sixth field of the lines does
        not matter. For more information, see
        http://scop.mrc-lmb.cam.ac.uk/scop/release-notes.html.
        """
        fields1 = cla_line_1.rstrip().split("\t")
        fields2 = cla_line_2.rstrip().split("\t")
        # compare the first five fields in a Cla line, which should be exactly
        # the same
        self.assertEqual(fields1[:5], fields2[:5])
        # compare the hierarchy key-value pairs, which are unordered
        self.assertCountEqual(fields1[5].split(","), fields2[5].split(","))

    def testParse(self):
        with open(support.DATA / "SCOP" / "dir.cla.scop.txt_test") as f:
            cla = f.read()
        with open(support.DATA / "SCOP" / "dir.des.scop.txt_test") as f:
            des = f.read()
        with open(support.DATA / "SCOP" / "dir.hie.scop.txt_test") as f:
            hie = f.read()

        scop = Scop(StringIO(cla), StringIO(des), StringIO(hie))

        cla_out = StringIO()
        scop.write_cla(cla_out)
        lines = zip(cla.rstrip().split("\n"), cla_out.getvalue().rstrip().split("\n"))
        for expected_line, line in lines:
            self._compare_cla_lines(expected_line, line)

        des_out = StringIO()
        scop.write_des(des_out)
        self.assertEqual(des_out.getvalue(), des)

        hie_out = StringIO()
        scop.write_hie(hie_out)
        self.assertEqual(hie_out.getvalue(), hie)

        domain = scop.getDomainBySid("d1hbia_")
        self.assertEqual(domain.sunid, 14996)

        domains = scop.getDomains()
        self.assertEqual(len(domains), 14)
        self.assertEqual(domains[4].sunid, 14988)

        dom = scop.getNodeBySunid(-111)
        self.assertIsNone(dom)
        dom = scop.getDomainBySid("no such domain")
        self.assertIsNone(dom)

    def testHieSunidMissingFromDes(self):
        with open(support.DATA / "SCOP" / "dir.cla.scop.txt_test") as f:
            cla = f.read()
        with open(support.DATA / "SCOP" / "dir.des.scop.txt_test") as f:
            des = f.read()
        with open(support.DATA / "SCOP" / "dir.hie.scop.txt_test") as f:
            hie = f.read() + "99999\t-\t-\n"

        stdout = StringIO()
        with redirect_stdout(stdout):
            with self.assertRaisesRegex(ValueError, "99999"):
                Scop(StringIO(cla), StringIO(des), StringIO(hie))
        self.assertEqual(stdout.getvalue(), "")

    def testClaSunidMissingFromDes(self):
        with open(support.DATA / "SCOP" / "dir.cla.scop.txt_test") as f:
            cla = f.read() + (
                "d9zzza_\t9zzz\tA:\ta.1.1.1\t999999\tcl=46456,cf=46457,"
                "sf=46458,fa=46459,dm=46460,sp=46461,px=999999\n"
            )
        with open(support.DATA / "SCOP" / "dir.des.scop.txt_test") as f:
            des = f.read()
        with open(support.DATA / "SCOP" / "dir.hie.scop.txt_test") as f:
            hie = f.read()

        with self.assertRaisesRegex(ValueError, "sunid 999999 in CLA file"):
            Scop(StringIO(cla), StringIO(des), StringIO(hie))

    def testClaDisagreesWithDes(self):
        with open(support.DATA / "SCOP" / "dir.des.scop.txt_test") as f:
            des = f.read()
        with open(support.DATA / "SCOP" / "dir.hie.scop.txt_test") as f:
            hie = f.read()
        with open(support.DATA / "SCOP" / "dir.cla.scop.txt_test") as f:
            cla = f.read()
        for old, new, message in [
            ("A:\ta.1.1.2\t14984", "A:\ta.1.1.3\t14984", "sccs a.1.1.3"),
            ("d3sdha_\t3sdh", "d3sdhz_\t3sdh", "sid d3sdhz_"),
            # 46465 is an 'sp' node in the DES file, not a 'px' domain
            ("A:\ta.1.1.2\t14984", "A:\ta.1.1.2\t46465", "'sp' node"),
        ]:
            with self.subTest(message=message):
                self.assertIn(old, cla)
                broken = StringIO(cla.replace(old, new, 1))
                with self.assertRaisesRegex(ValueError, message):
                    Scop(broken, StringIO(des), StringIO(hie))

    def testSccsOrder(self):
        self.assertEqual(cmp_sccs("a.1.1.1", "a.1.1.1"), 0)
        self.assertEqual(cmp_sccs("a.1.1.2", "a.1.1.1"), 1)
        self.assertEqual(cmp_sccs("a.1.1.2", "a.1.1.11"), -1)
        self.assertEqual(cmp_sccs("a.1.2.2", "a.1.1.11"), 1)
        self.assertEqual(cmp_sccs("a.1.2.2", "a.5.1.11"), -1)
        self.assertEqual(cmp_sccs("b.1.2.2", "a.5.1.11"), 1)
        self.assertEqual(cmp_sccs("b.1.2.2", "b.1.2"), 1)

    def testParseDomain(self):
        s = ">d1tpt_1 a.46.2.1 (1-70) Thymidine phosphorylase {Escherichia coli}"
        dom = parse_domain(s)

        self.assertEqual(dom.sid, "d1tpt_1")
        self.assertEqual(dom.sccs, "a.46.2.1")
        self.assertEqual(dom.residues.pdbid, "1tpt")
        self.assertEqual(dom.description, "Thymidine phosphorylase {Escherichia coli}")

        s2 = "d1tpt_1 a.46.2.1 (1tpt 1-70) Thymidine phosphorylase {E. coli}"
        self.assertEqual(s2, str(parse_domain(s2)))

        # Genetic domains (See Astral release notes)
        s3 = "g1cph.1 g.1.1.1 (1cph B:,A:) Insulin {Cow (Bos taurus)}"
        self.assertEqual(s3, str(parse_domain(s3)))

        s4 = "e1cph.1a g.1.1.1 (1cph A:) Insulin {Cow (Bos taurus)}"
        self.assertEqual(s4, str(parse_domain(s4)))

        # Raw Astral header
        s5 = ">e1cph.1a g.1.1.1 (A:) Insulin {Cow (Bos taurus)}"
        self.assertEqual(s4, str(parse_domain(s5)))

        self.assertRaises(ValueError, parse_domain, "Totally wrong")

    def test_parse_domain_whitespace(self):
        """parse_domain fails fast on a header with a long whitespace run."""
        # With no sccs, two or more whitespace characters precede the
        # residues, as Domain.__str__ writes them; one is not enough.
        s = "d1tpt_1  (1tpt 1-70) Thymidine phosphorylase {E. coli}"
        dom = parse_domain(s)
        self.assertEqual(dom.sccs, "")
        self.assertEqual(s, str(dom))
        with self.assertRaisesRegex(ValueError, r"^Domain: d1tpt_1 \(1-70\)"):
            parse_domain("d1tpt_1 (1-70) Thymidine phosphorylase")

        # The old pattern tried every way of splitting a whitespace run
        # around an empty sccs before it failed, which took time quadratic
        # in the length of the run. It still raised ValueError in the end,
        # so time it in a child process: these two headers take the old
        # pattern about six minutes each and the new one milliseconds, and
        # the timeout fails the test on the old pattern.
        code = """\
from Bio.SCOP import parse_domain
for term in [
    "d1tpt_1" + " " * 500_000 + "x",
    "d1tpt_1" + " " * 200_000 + "(" + "1" * 200_000,
]:
    try:
        parse_domain(term)
    except ValueError:
        print("ValueError")
"""
        env = os.environ.copy()
        env["PYTHONPATH"] = os.pathsep.join(
            filter(
                None,
                [os.path.dirname(os.path.dirname(Bio.__file__)), env.get("PYTHONPATH")],
            )
        )
        result = subprocess.run(
            [sys.executable, "-c", code],
            timeout=60,
            capture_output=True,
            text=True,
            env=env,
        )
        self.assertEqual(result.stdout, "ValueError\nValueError\n", result.stderr)

        # Shorter headers of the same two kinds, which the old pattern
        # failed on in several seconds each, give the usual error message.
        for term in [
            "d1tpt_1" + " " * 100_000 + "x",
            "d1tpt_1" + " " * 20_000 + "(" + "1" * 20_000,
        ]:
            with self.subTest(length=len(term)):
                with self.assertRaisesRegex(ValueError, "^Domain: d1tpt_1   "):
                    parse_domain(term)

    def testConstructFromDirectory(self):
        scop = Scop(dir_path=support.DATA / "SCOP", version="test")
        self.assertIsInstance(scop, Scop)

        domain = scop.getDomainBySid("d1hbia_")
        self.assertEqual(domain.sunid, 14996)

    def testGetAscendent(self):
        scop = Scop(dir_path=support.DATA / "SCOP", version="test")
        domain = scop.getDomainBySid("d1hbia_")

        # get the fold
        fold = domain.getAscendent("cf")
        self.assertEqual(fold.sunid, 46457)

        # get the superfamily
        sf = domain.getAscendent("superfamily")
        self.assertEqual(sf.sunid, 46458)

        # px has no px ascendent
        px = domain.getAscendent("px")
        self.assertIsNone(px)

        # an sf has no px ascendent
        px2 = sf.getAscendent("px")
        self.assertIsNone(px2)

    def test_get_descendents(self):
        """Test getDescendents method."""
        scop = Scop(dir_path=support.DATA / "SCOP", version="test")
        fold = scop.getNodeBySunid(46457)

        # get px descendents
        domains = fold.getDescendents("px")
        self.assertEqual(len(domains), 14)
        for d in domains:
            self.assertEqual(d.type, "px")

        sfs = fold.getDescendents("superfamily")
        self.assertEqual(len(sfs), 1)
        for d in sfs:
            self.assertEqual(d.type, "sf")

        # cl has no cl descendent
        cl = fold.getDescendents("cl")
        self.assertEqual(cl, [])

    def test_get_descendents_not_a_tree(self):
        """getDescendents raises ValueError if the HIE file gives no tree."""
        with open(support.DATA / "SCOP" / "dir.cla.scop.txt_test") as f:
            cla = f.read()
        with open(support.DATA / "SCOP" / "dir.des.scop.txt_test") as f:
            des = f.read()
        with open(support.DATA / "SCOP" / "dir.hie.scop.txt_test") as f:
            hie = f.read()
        domains = "14984,14985,14986,14987,14988,14989,14990,14991,14992"

        # Without the check the walk never ends; stop it if so.
        get_children = Node.getChildren
        calls = 0

        def bounded_get_children(node):
            nonlocal calls
            calls += 1
            if calls > 1000:
                raise AssertionError("getDescendents did not terminate")
            return get_children(node)

        for old, new, sunid in [
            # sf 46458 lists itself as a child
            ("46458\t46457\t46463\n", "46458\t46457\t46458,46463\n", 46458),
            # fa 46463 lists its child dm 46464 twice
            ("46463\t46458\t46464\n", "46463\t46458\t46464,46464\n", 46464),
            # sp 46465 lists its ancestor cf 46457 as a child
            (f"\t{domains}", f"\t46457,{domains}", 46457),
        ]:
            with self.subTest(new=new):
                self.assertIn(old, hie)
                broken = StringIO(hie.replace(old, new, 1))
                scop = Scop(StringIO(cla), StringIO(des), broken)
                calls = 0
                with mock.patch.object(Node, "getChildren", bounded_get_children):
                    with self.assertRaisesRegex(
                        ValueError, f"not a tree.* reaches node {sunid} twice"
                    ):
                        scop.getRoot().getDescendents("px")

    def test_get_ascendent_not_a_tree(self):
        """getAscendent and toClaRecord raise ValueError on a parent cycle."""
        with open(support.DATA / "SCOP" / "dir.cla.scop.txt_test") as f:
            cla = f.read()
        with open(support.DATA / "SCOP" / "dir.des.scop.txt_test") as f:
            des = f.read()
        with open(support.DATA / "SCOP" / "dir.hie.scop.txt_test") as f:
            hie = f.read()

        # Without the check the walk never ends; stop it if so.
        get_parent = Node.getParent
        calls = 0

        def bounded_get_parent(node):
            nonlocal calls
            calls += 1
            if calls > 1000:
                raise AssertionError("walk up the hierarchy did not terminate")
            return get_parent(node)

        for old, new, node_type, sunid in [
            # sp 46465 gives itself as its parent
            ("46465\t46464\t", "46465\t46465\t", "cf", 46465),
            # dm 46464 gives its child sp 46465 as its parent
            ("46464\t46463\t", "46464\t46465\t", "cf", 46465),
            # cl 46456 gives its child cf 46457 as its parent, so a walk up
            # to the root (as toClaRecord makes) never reaches it
            ("46456\t0\t", "46456\t46457\t", "ro", 46457),
        ]:
            with self.subTest(new=new):
                self.assertIn(old, hie)
                broken = StringIO(hie.replace(old, new, 1))
                scop = Scop(StringIO(cla), StringIO(des), broken)
                domain = scop.getDomainBySid("d1hbia_")
                message = (
                    f"not a tree: walking up from node 14996 reaches node {sunid} twice"
                )
                with mock.patch.object(Node, "getParent", bounded_get_parent):
                    calls = 0
                    with self.assertRaisesRegex(ValueError, message):
                        domain.getAscendent(node_type)
                    calls = 0
                    with self.assertRaisesRegex(ValueError, message):
                        domain.toClaRecord()
                    calls = 0
                    with self.assertRaisesRegex(ValueError, "not a tree"):
                        scop.write_cla(StringIO())


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
