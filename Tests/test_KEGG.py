# Copyright 2001 by Tarjei Mikkelsen.  All rights reserved.
# Revisions copyright 2007 by Michiel de Hoon. All rights reserved.
# Revisions copyright 2017 by Peter Cock. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests the basic functionality of the KEGG parsers."""

import unittest
from io import StringIO

import support

from Bio.KEGG import Compound
from Bio.KEGG import Enzyme
from Bio.KEGG import Map
from Bio.Pathway import System


class EnzymeTests(unittest.TestCase):
    """Tests for Bio.KEGG.Enzyme."""

    def test_sample(self):
        """enzyme.sample tests."""
        with open(support.DATA / "KEGG" / "enzyme.sample") as handle:
            records = list(Enzyme.parse(handle))
        self.assertEqual(len(records), 8)
        self.assertEqual(records[0].entry, "1.1.1.1")
        self.assertEqual(
            records[0].name,
            [
                "alcohol dehydrogenase",
                "aldehyde reductase",
                "ADH",
                "alcohol dehydrogenase (NAD)",
                "aliphatic alcohol dehydrogenase",
                "ethanol dehydrogenase",
                "NAD-dependent alcohol dehydrogenase",
                "NAD-specific aromatic alcohol dehydrogenase",
                "NADH-alcohol dehydrogenase",
                "NADH-aldehyde dehydrogenase",
                "primary alcohol dehydrogenase",
                "yeast alcohol dehydrogenase",
            ],
        )
        self.assertEqual(
            records[0].pathway,
            [
                ("PATH", "ec00010", "Glycolysis / Gluconeogenesis"),
                ("PATH", "ec00071", "Fatty acid degradation"),
                ("PATH", "ec00260", "Glycine, serine and threonine metabolism"),
                ("PATH", "ec00350", "Tyrosine metabolism"),
                ("PATH", "ec00592", "alpha-Linolenic acid metabolism"),
                ("PATH", "ec00625", "Chloroalkane and chloroalkene degradation"),
                ("PATH", "ec00626", "Naphthalene degradation"),
                ("PATH", "ec00830", "Retinol metabolism"),
                ("PATH", "ec00980", "Metabolism of xenobiotics by cytochrome P450"),
                ("PATH", "ec00982", "Drug metabolism - cytochrome P450"),
                ("PATH", "ec01100", "Metabolic pathways"),
                ("PATH", "ec01110", "Biosynthesis of secondary metabolites"),
                ("PATH", "ec01120", "Microbial metabolism in diverse environments"),
                ("PATH", "ec01130", "Biosynthesis of antibiotics"),
            ],
        )
        self.assertEqual(
            records[0].dblinks,
            [
                ("ExplorEnz - The Enzyme Database", ["1.1.1.1"]),
                ("IUBMB Enzyme Nomenclature", ["1.1.1.1"]),
                ("ExPASy - ENZYME nomenclature database", ["1.1.1.1"]),
                ("UM-BBD (Biocatalysis/Biodegradation Database)", ["1.1.1.1"]),
                ("BRENDA, the Enzyme Database", ["1.1.1.1"]),
                ("CAS", ["9031-72-5"]),
            ],
        )
        self.assertEqual(records[-1].entry, "2.7.2.1")
        self.assertEqual(
            str(records[-1]).replace(" ", "").split("\n")[:10],
            [
                "ENTRYEC2.7.2.1",
                "NAMEacetatekinase",
                "acetokinase",
                "AckA",
                "AK",
                "acetickinase",
                "acetatekinase(phosphorylating)",
                "CLASSTransferases;",
                "Transferringphosphorus-containinggroups;",
                "Phosphotransferaseswithacarboxygroupasacceptor",
            ],
        )

    def test_irregular(self):
        """enzyme.irregular tests."""
        with open(support.DATA / "KEGG" / "enzyme.irregular") as handle:
            records = list(Enzyme.parse(handle))
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].entry, "1.14.18.1")
        self.assertEqual(records[-1].entry, "3.4.21.50")

    def test_new(self):
        """enzyme.new tests."""
        with open(support.DATA / "KEGG" / "enzyme.new") as handle:
            records = list(Enzyme.parse(handle))
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].entry, "6.2.1.25")

    def test_4letter(self):
        """enzyme.4letter tests."""
        with open(support.DATA / "KEGG" / "enzyme.4letter") as handle:
            records = list(Enzyme.parse(handle))
            self.assertEqual(len(records), 1)
        self.assertEqual(records[0].entry, "5.4.2.2")
        self.assertEqual(records[0].genes[0], ("HSA", ["5236", "55276"]))
        self.assertEqual(records[0].genes[8], ("CSAB", ["103224690", "103246223"]))

    def test_exceptions(self):
        """enzyme.exceptions tests."""
        with open(support.DATA / "KEGG" / "enzyme.sample") as handle:
            with self.assertRaises(ValueError) as context:
                list(Enzyme.read(handle))
            self.assertIn(
                "More than one record found in handle", str(context.exception)
            )
            records = Enzyme.parse(handle)
            for i in range(6):
                next(records)
            self.assertRaises(StopIteration, next, records)


class CompoundTests(unittest.TestCase):
    """Bio.KEGG.Compound tests."""

    def test_sample(self):
        """compound.sample tests."""
        with open(support.DATA / "KEGG" / "compound.sample") as handle:
            records = list(Compound.parse(handle))
        self.assertEqual(len(records), 8)
        self.assertEqual(records[1].entry, "C00017")
        self.assertEqual(records[1].mass, "")
        self.assertEqual(records[1].formula, "C2H4NO2R(C2H2NOR)n")
        self.assertEqual(records[1].name, ["Protein"])
        self.assertEqual(
            records[1].pathway, [("PATH", "map00450", "Selenocompound metabolism")]
        )
        self.assertEqual(len(records[1].enzyme), 21)
        self.assertEqual(records[1].enzyme[0], ("2.3.2.6"))
        self.assertEqual(records[1].structures, [])
        self.assertEqual(records[1].dblinks[0], ("PubChem", ["3319"]))
        self.assertEqual(
            str(records[-1]).replace(" ", "").split("\n")[:11],
            [
                "ENTRYC01386",
                "NAMENH2Mec",
                "7-Amino-4-methylcoumarin",
                "FORMULAC10H9NO2",
                "MASS175.0633",
                "DBLINKSCAS:26093-31-2",
                "PubChem:4580",
                "ChEBI:51771",
                "ChEMBL:CHEMBL270672",
                "KNApSAcK:C00048593",
                "PDB-CCD:MCM",
            ],
        )

    def test_pathway_round_trip(self):
        """str() of a record keeps the PATHWAY map ids and names."""
        with open(support.DATA / "KEGG" / "compound.sample") as handle:
            record = list(Compound.parse(handle))[2]
        self.assertEqual(record.entry, "C00099")
        self.assertIn("PATHWAY     map00240  Pyrimidine metabolism\n", str(record))
        reparsed = next(Compound.parse(StringIO(str(record) + "\n")))
        self.assertEqual(reparsed.pathway, record.pathway)
        self.assertEqual(len(reparsed.pathway), 8)

    def test_pathway_round_trip_odd_names(self):
        """PATHWAY names with a double space, or long enough to wrap, read back."""
        long_name = (
            "Biosynthesis of various plant secondary metabolites"
            " and other very long pathway names here"
        )
        record = Compound.Record()
        record.entry = "C00001"
        record.pathway = [
            ("PATH", "map00240", "Foo  bar"),
            ("PATH", "map00999", long_name),
            ("PATH", "map00010", "Glycolysis"),
        ]
        self.assertIn("\n" + " " * 28 + "other very long", str(record))
        reparsed = next(Compound.parse(StringIO(str(record) + "\n")))
        self.assertEqual(reparsed.pathway, record.pathway)

    def test_pathway_round_trip_double_space_at_wrap(self):
        """A double space where a PATHWAY name wraps reads back as two spaces."""
        record = Compound.Record()
        record.entry = "C00001"
        record.pathway = [
            # the wrap leaves one space at the end of the first line
            ("PATH", "map00240", "foo " * 13 + "word  suffix"),
            # the wrap leaves one space at the start of the continuation
            ("PATH", "map00241", "foo " * 13 + "words  suffix"),
        ]
        text = str(record)
        self.assertIn(" word \n" + " " * 28 + "suffix\n", text)
        self.assertIn(" words\n" + " " * 29 + "suffix\n", text)
        reparsed = next(Compound.parse(StringIO(text + "\n")))
        self.assertEqual(reparsed.pathway, record.pathway)

    def test_structures(self):
        """STRUCTURES lines, including a continuation line, are parsed."""
        text = (
            "ENTRY       C00001                      Compound\n"
            "NAME        H2O\n"
            "STRUCTURES  PDB: 1ABC  2DEF  3GHI\n"
            "                 4JKL\n"
            "DBLINKS     PubChem: 3303\n"
            "///\n"
        )
        record = next(Compound.parse(StringIO(text)))
        self.assertEqual(record.structures, [("PDB", ["1ABC", "2DEF", "3GHI", "4JKL"])])
        self.assertEqual(record.dblinks, [("PubChem", ["3303"])])
        reparsed = next(Compound.parse(StringIO(str(record) + "\n")))
        self.assertEqual(reparsed.structures, record.structures)

    def test_structures_without_database(self):
        """An empty STRUCTURES line is skipped; ids with no database raise."""
        text = "ENTRY       C00001\nSTRUCTURES  \n            PDB: 1ABC\n///\n"
        record = next(Compound.parse(StringIO(text)))
        self.assertEqual(record.structures, [("PDB", ["1ABC"])])
        text = "ENTRY       C00001\nSTRUCTURES  1ABC\n///\n"
        with self.assertRaisesRegex(ValueError, "^STRUCTURES ids before any database"):
            next(Compound.parse(StringIO(text)))

    def test_irregular(self):
        """compound.irregular tests."""
        with open(support.DATA / "KEGG" / "compound.irregular") as handle:
            records = list(Compound.parse(handle))
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].entry, "C01454")

    def test_mass(self):
        """record.mass tests."""
        with open(support.DATA / "KEGG" / "compound.sample") as handle:
            records = list(Compound.parse(handle))
        self.assertEqual(records[0].entry, "C00023")
        self.assertEqual(records[0].mass, "55.9349")
        self.assertEqual(records[2].entry, "C00099")
        self.assertEqual(records[2].mass, "89.0477")
        self.assertEqual(records[3].entry, "C00294")
        self.assertEqual(records[3].mass, "268.0808")


class MapTests(unittest.TestCase):
    """Bio.KEGG.Map tests."""

    def test_map00950(self):
        """map00950.rea tests."""
        system = System()
        with open(support.DATA / "KEGG" / "map00950.rea") as handle:
            for reaction in Map.parse(handle):
                system.add_reaction(reaction)
        rxs = system.reactions()
        self.assertEqual(len(rxs), 56)
        # sort the reaction output by the string names, so that the
        # output will be consistent between python versions
        rxs.sort(key=lambda x: str(x))  # noqa: E731
        self.assertEqual(
            str(rxs[0]),
            "(R)-N-Methylcoclaurine + (S)-Coclaurine + NADPH + O2 "
            "<=> 2'-Norberbamunine + 2 H2O + NADP",
        )
        self.assertEqual(str(rxs[-1]), "Tyramine <=> Dopamine")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
