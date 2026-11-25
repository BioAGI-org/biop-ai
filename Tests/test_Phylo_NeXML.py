# Copyright (C) 2013 by Ben Morris (ben@bendmorris.com)
# based on code by Eric Talevich (eric.talevich@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Unit tests for the NeXML and NeXMLIO modules."""

import os
import tempfile
import copy
import unittest
from io import BytesIO

import support

from Bio import BiopythonParserWarning
from Bio import Phylo

# Example NeXML files
nexml_files = (
    "characters.xml",
    "edgelabels.xml",
    "meta_taxa.xml",
    "meta_types.xml",
    "nexml.xml",
    "phenoscape.xml",
    "sets.xml",
    "taxa.xml",
    "timetree.xml",
    "tolweb.xml",
    "treebase-record.xml",
    "trees-uris.xml",
    "trees.xml",
)

tree_counts = {
    "taxa.xml": 0,
    "timetree.xml": 38,
    "phenoscape.xml": 0,
    "nexml.xml": 0,
    "meta_types.xml": 0,
    "meta_taxa.xml": 0,
    "trees.xml": 2,
    "characters.xml": 0,
}


class ParseTests(unittest.TestCase):
    """Tests for proper parsing of example NeXML files."""

    def test_parse(self):
        """Extract and count phylogenetic trees using Phylo.parse."""
        for filename in nexml_files:
            count = tree_counts.get(filename, 1)
            path = support.DATA / "NeXML" / filename
            msg = f"Failed parser test for {path}"
            trees = list(Phylo.parse(path, "nexml"))
            self.assertEqual(len(trees), count, msg=msg)

    def test_entity_declaration(self):
        """Check an XML entity declaration raises ValueError."""
        data = b"""\
<?xml version="1.0"?>
<!DOCTYPE nex:nexml [<!ENTITY a "aaaaaaaaaa">]>
<nex:nexml xmlns:nex="http://www.nexml.org/2009" version="&a;&a;&a;"/>
"""
        with self.assertRaisesRegex(ValueError, "^XML entity declarations are not"):
            list(Phylo.parse(BytesIO(data), "nexml"))

    def test_meta_named_like_clade_attribute(self):
        """Check a meta property cannot replace a Clade attribute."""
        data = b"""\
<?xml version="1.0"?>
<nex:nexml xmlns:nex="http://www.nexml.org/2009" xmlns="http://www.nexml.org/2009">
 <trees id="t" otus="o"><tree id="tr">
  <node id="n1" root="true"/>
  <node id="n2" otu="A">
   <meta property="is_terminal">x</meta><meta property="clades">x</meta>
   <meta property="root">x</meta><meta property="branch_length">x</meta>
   <meta property="self">x</meta><meta property="__deepcopy__">x</meta>
   <meta property="dc:creator">me</meta>
  </node>
  <node id="n3" otu="B"/>
  <edge id="e1" source="n1" target="n2"/>
  <edge id="e2" source="n1" target="n3"/>
 </tree></trees>
</nex:nexml>
"""
        with self.assertWarns(BiopythonParserWarning):
            tree = Phylo.read(BytesIO(data), "nexml")
        self.assertEqual(sorted(c.name for c in tree.get_terminals()), ["A", "B"])
        self.assertEqual(tree.total_branch_length(), 3.0)
        (clade,) = tree.find_clades("A")
        self.assertEqual(getattr(clade, "dc:creator"), "me")
        # A name such as __deepcopy__ would be called by copy.deepcopy
        self.assertEqual(copy.deepcopy(tree).count_terminals(), 2)


class WriterTests(unittest.TestCase):
    """NeXML writer tests."""

    def check(self, path):
        """Parse, rewrite and retest an example phylogeny file."""
        msg = f"Failed NeXMLIO writer test for {path}"

        # Using Phylo.NeXMLIO directly, binary mode handles
        with open(path, "rb") as stream:
            t1 = next(Phylo.NeXMLIO.Parser(stream).parse())
        stream = BytesIO()
        Phylo.NeXMLIO.write([t1], stream)
        stream.seek(0)
        t2 = next(Phylo.NeXMLIO.Parser(stream).parse())
        self.compare(t1, t2, msg)

        # Using Phylo.parse/write, using filenames
        msg = f"Failed Phylo API writer test for {path}"
        t1 = next(Phylo.parse(path, "nexml"))
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = os.path.join(tmpdir, "tree.xml")
            Phylo.write([t1], tmp, "nexml")
            with open(tmp, "rb") as handle:
                t2 = next(Phylo.parse(handle, "nexml"))
        self.compare(t1, t2, msg)

    def compare(self, t1, t2, msg=None):
        """Compare two trees."""
        for prop_name in ("name", "branch_length", "confidence"):
            p1 = sorted(
                getattr(n, prop_name)
                for n in t1.get_terminals()
                if getattr(n, prop_name)
            )
            p2 = sorted(
                getattr(n, prop_name)
                for n in t2.get_terminals()
                if getattr(n, prop_name)
            )
            self.assertEqual(p1, p2, msg=msg)

    def test_write(self):
        """Test for serialization of objects to NeXML format."""
        for filename in nexml_files:
            count = tree_counts.get(filename, 1)
            if count > 0:
                path = support.DATA / "NeXML" / filename
                self.check(path)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
