# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for Affy module."""

import io
import os
import struct
import tempfile
import unittest

try:
    import numpy.testing
    from numpy import array
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.Affy.CelFile"
    ) from None

try:
    import tracemalloc
except ImportError:
    # PyPy has no _tracemalloc.
    tracemalloc = None

import support

from Bio.Affy import CelFile


class AffyTest(unittest.TestCase):
    def setUp(self):
        self.affy3 = support.DATA / "Affy" / "affy_v3_example.CEL"
        self.affy4 = support.DATA / "Affy" / "affy_v4_example.CEL"
        self.affy4Bad = support.DATA / "Affy" / "affy_v4_bad_example.CEL"
        with open(self.affy4Bad, "wb") as f:
            self.writeExampleV4(f, bad=True)

    def tearDown(self):
        os.remove(self.affy4Bad)

    # tests the Affymetrix v3 parser
    def testAffy3(self):
        with open(self.affy3) as f:
            record = CelFile.read(f)
            self.assertGreater(len(record.DatHeader), 0)
            self.assertEqual(record.intensities.shape, (5, 5))
            self.assertEqual(record.intensities.shape, record.stdevs.shape)
            self.assertEqual(record.intensities.shape, record.npix.shape)
            self.assertEqual(record.ncols, 5)
            self.assertEqual(record.nrows, 5)
            self.assertEqual(record.version, 3)
            self.assertEqual(record.GridCornerUL, (206, 129))
            self.assertEqual(record.GridCornerUR, (3570, 107))
            self.assertEqual(record.GridCornerLR, (3597, 3470))
            self.assertEqual(record.GridCornerLL, (234, 3492))
            self.assertEqual(record.DatHeader["filename"], "1g_A9AF")
            self.assertEqual(record.DatHeader["CLS"], 3684)
            self.assertEqual(record.DatHeader["RWS"], 3684)
            self.assertEqual(record.DatHeader["XIN"], 1)
            self.assertEqual(record.DatHeader["YIN"], 1)
            self.assertEqual(record.DatHeader["VE"], 30)
            self.assertAlmostEqual(record.DatHeader["laser-power"], 2.0)
            self.assertEqual(record.DatHeader["scan-date"], "08/23/07")
            self.assertEqual(record.DatHeader["scan-time"], "11:23:24")
            self.assertEqual(record.DatHeader["scanner-id"], "50205880")
            self.assertEqual(record.DatHeader["scanner-type"], "M10")
            self.assertEqual(record.DatHeader["array-type"], "Tgondii_SNP1.1sq")
            self.assertEqual(record.DatHeader["filter-wavelength"], 570)
            self.assertAlmostEqual(record.DatHeader["arc-radius"], 25356.509766)
            self.assertAlmostEqual(record.DatHeader["laser-spotsize"], 3.5)
            self.assertAlmostEqual(record.DatHeader["pixel-size"], 1.56)
            self.assertEqual(record.DatHeader["image-orientation"], 6)
            self.assertEqual(record.Algorithm, "Percentile")
            self.assertEqual(len(record.AlgorithmParameters), 16)
            self.assertEqual(record.AlgorithmParameters["Percentile"], 75)
            self.assertEqual(record.AlgorithmParameters["CellMargin"], 2)
            self.assertAlmostEqual(record.AlgorithmParameters["OutlierHigh"], 1.500)
            self.assertAlmostEqual(record.AlgorithmParameters["OutlierLow"], 1.004)
            self.assertEqual(record.AlgorithmParameters["AlgVersion"], "6.0")
            self.assertEqual(
                record.AlgorithmParameters["FixedCellSize"], True
            )  # noqa: A502
            self.assertEqual(record.AlgorithmParameters["FullFeatureWidth"], 7)
            self.assertEqual(record.AlgorithmParameters["FullFeatureHeight"], 7)
            self.assertEqual(
                record.AlgorithmParameters["IgnoreOutliersInShiftRows"], False
            )  # noqa: A502
            self.assertEqual(
                record.AlgorithmParameters["FeatureExtraction"], True
            )  # noqa: A502
            self.assertEqual(record.AlgorithmParameters["PoolWidthExtenstion"], 2)
            self.assertEqual(record.AlgorithmParameters["PoolHeightExtension"], 2)
            self.assertEqual(
                record.AlgorithmParameters["UseSubgrids"], False
            )  # noqa: A502
            self.assertEqual(
                record.AlgorithmParameters["RandomizePixels"], False
            )  # noqa: A502
            self.assertEqual(record.AlgorithmParameters["ErrorBasis"], "StdvMean")
            self.assertAlmostEqual(record.AlgorithmParameters["StdMult"], 1.0)
            self.assertEqual(record.NumberCells, 25)

            global message
            try:
                numpy.testing.assert_allclose(
                    record.intensities,
                    [
                        [234.0, 170.0, 22177.0, 164.0, 22104.0],
                        [188.0, 188.0, 21871.0, 168.0, 21883.0],
                        [188.0, 193.0, 21455.0, 198.0, 21300.0],
                        [188.0, 182.0, 21438.0, 188.0, 20945.0],
                        [193.0, 20370.0, 174.0, 20605.0, 168.0],
                    ],
                )
                message = None
            except AssertionError as err:
                message = str(err)
            if message is not None:
                self.fail(message)
            try:
                numpy.testing.assert_allclose(
                    record.stdevs,
                    [
                        [24.0, 34.5, 2669.0, 19.7, 3661.2],
                        [29.8, 29.8, 2795.9, 67.9, 2792.4],
                        [29.8, 88.7, 2976.5, 62.0, 2914.5],
                        [29.8, 76.2, 2759.5, 49.2, 2762.0],
                        [38.8, 2611.8, 26.6, 2810.7, 24.1],
                    ],
                )
                message = None
            except AssertionError as err:
                message = str(err)
            if message is not None:
                self.fail(message)
            try:
                numpy.testing.assert_array_equal(
                    record.npix,
                    [
                        [25, 25, 25, 25, 25],
                        [25, 25, 25, 25, 25],
                        [25, 25, 25, 25, 25],
                        [25, 25, 25, 25, 25],
                        [25, 25, 25, 25, 25],
                    ],
                )
                message = None
            except AssertionError as err:
                message = str(err)
            if message is not None:
                self.fail(message)
            self.assertEqual(record.nmask, 3)
            try:
                numpy.testing.assert_array_equal(
                    record.mask,
                    [
                        [False, False, False, False, False],
                        [False, False, False, True, True],
                        [False, False, False, False, True],
                        [False, False, False, False, False],
                        [False, False, False, False, False],
                    ],
                )
                message = None
            except AssertionError as err:
                message = str(err)
            if message is not None:
                self.fail(message)
            self.assertEqual(record.noutliers, 3)
            try:
                numpy.testing.assert_array_equal(
                    record.outliers,
                    [
                        [False, False, False, False, False],
                        [False, True, True, False, False],
                        [False, False, False, False, False],
                        [False, True, False, False, False],
                        [False, False, False, False, False],
                    ],
                )
                message = None
            except AssertionError as err:
                message = str(err)
            if message is not None:
                self.fail(message)
            self.assertEqual(record.nmodified, 3)
            try:
                numpy.testing.assert_allclose(
                    record.modified,
                    [
                        [0.0, 0.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 189.0, 220.0],
                        [0.0, 0.0, 0.0, 21775.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0, 0.0],
                    ],
                )
                message = None
            except AssertionError as err:
                message = str(err)
            if message is not None:
                self.fail(message)

    def testAffy4(self):
        with open(self.affy4, "rb") as f:
            record = CelFile.read(f)
        self.assertEqual(record.intensities.shape, (5, 5))
        self.assertEqual(record.intensities.shape, record.stdevs.shape)
        self.assertEqual(record.intensities.shape, record.npix.shape)
        self.assertEqual(record.ncols, 5)
        self.assertEqual(record.nrows, 5)
        global message
        try:
            numpy.testing.assert_allclose(
                record.intensities,
                [
                    [0.0, 1.0, 2.0, 3.0, 4.0],
                    [5.0, 6.0, 7.0, 8.0, 9.0],
                    [10.0, 11.0, 12.0, 13.0, 14.0],
                    [15.0, 16.0, 17.0, 18.0, 19.0],
                    [20.0, 21.0, 22.0, 23.0, 24.0],
                ],
            )
            message = None
        except AssertionError as err:
            message = str(err)
        if message is not None:
            self.fail(message)
        try:
            numpy.testing.assert_allclose(
                record.stdevs,
                [
                    [0.0, -1.0, -2.0, -3.0, -4.0],
                    [-5.0, -6.0, -7.0, -8.0, -9.0],
                    [-10.0, -11.0, -12.0, -13.0, -14.0],
                    [-15.0, -16.0, -17.0, -18.0, -19.0],
                    [-20.0, -21.0, -22.0, -23.0, -24.0],
                ],
            )
            message = None
        except AssertionError as err:
            message = str(err)
        if message is not None:
            self.fail(message)
        try:
            numpy.testing.assert_allclose(
                record.npix,
                [
                    [9, 9, 9, 9, 9],
                    [9, 9, 9, 9, 9],
                    [9, 9, 9, 9, 9],
                    [9, 9, 9, 9, 9],
                    [9, 9, 9, 9, 9],
                ],
            )
            message = None
        except AssertionError as err:
            message = str(err)
        if message is not None:
            self.fail(message)
        self.assertEqual(len(record.AlgorithmParameters), 329)
        self.assertEqual(len(record.GridCornerUL), 7)
        self.assertEqual(record.AlgorithmParameters[-3:], "169")

    def testAffyBadHeader(self):
        with self.assertRaises(CelFile.ParserError):
            with open(self.affy4Bad, "rb") as f:
                record = CelFile.read(f)

    def testAffyDatHeaderWithoutBrackets(self):
        with open(self.affy3) as f:
            data = f.read()
        self.assertIn("DatHeader=[11..65533]", data)
        data = data.replace("DatHeader=[11..65533]", "DatHeader=11..65533", 1)
        with self.assertRaisesRegex(CelFile.ParserError, r"Expected \[min..max\]"):
            CelFile.read(io.StringIO(data))

    def testAffyDatHeaderWrongSubfieldCount(self):
        # The DatHeader of a version 3 file must hold exactly 12 subfields
        # separated by \x14; any other count used to give an IndexError or a
        # misleading error from int() on a shifted field, or pass silently.
        with open(self.affy3) as f:
            data = f.read()
        cases = {
            # Seven subfields: the five numeric ones at the end removed
            7: data.replace(" 570 \x14 25356.509766 \x14", "", 1).replace(
                " \x14 3.500000 \x14 1.5600 \x14 6", "", 1
            ),
            # Eleven subfields: the last one cut off
            11: data.replace(" \x14 6\n", "\n", 1),
            # Thirteen subfields: an extra one in the middle
            13: data.replace(" \x14 1.5600 \x14", " \x14 99 \x14 1.5600 \x14", 1),
        }
        for count, bad in cases.items():
            with self.subTest(count=count):
                self.assertNotEqual(bad, data)
                with self.assertRaisesRegex(
                    CelFile.ParserError,
                    rf"Expected 12 subfields in DatHeader, found {count}$",
                ):
                    CelFile.read(io.StringIO(bad))

    def testAffyDatHeaderScannerIdWithoutType(self):
        # A DatHeader whose first subfield holds a scanner ID but no scanner
        # type, as in GEO sample GSM2810894 (biopython/biopython#4144).
        with open(self.affy3) as f:
            data = f.read()
        self.assertIn(" 50205880  M10   \x14", data)
        data = data.replace(" 50205880  M10   \x14", " 50205880   \x14", 1)
        record = CelFile.read(io.StringIO(data))
        self.assertEqual(record.DatHeader["scanner-id"], "50205880")
        self.assertNotIn("scanner-type", record.DatHeader)
        self.assertEqual(record.DatHeader["array-type"], "Tgondii_SNP1.1sq")
        self.assertEqual(record.DatHeader["filter-wavelength"], 570)
        self.assertEqual(record.DatHeader["image-orientation"], 6)

    def testAffy3BadSizes(self):
        # The version 3 parser used to size its arrays from Rows and Cols in
        # the header alone. A short file declaring a big array came back as
        # that many cells of zeros, and memory grew with the declared size,
        # not with the file. Cells missing from the INTENSITY section were
        # silently left as zero.
        with open(self.affy3) as f:
            data = f.read()
        dims = "Cols=5\nRows=5\n"
        last_cell = "  4\t  4\t168.0\t24.1\t 25\n"
        self.assertIn(dims, data)
        self.assertIn(last_cell, data)
        start = data.index("[INTENSITY]\n")
        intensity = data[start : data.index("[MASKS]")]
        header_only = data[: start + len("[INTENSITY]\n")]
        cases = [
            (
                data.replace(dims, "Cols=2000\nRows=2000\n"),
                "has 2000 rows and 2000 columns but values for 25 cells",
            ),
            (
                data.replace(last_cell, ""),
                "has 5 rows and 5 columns but values for 24 cells",
            ),
            (header_only, "has 5 rows and 5 columns but values for 0 cells"),
            # Ends before the INTENSITY section; used to give no intensities
            (data[:start], "has no INTENSITY section"),
            (data[: data.index("[HEADER]")], "has no INTENSITY section"),
            (
                data.replace(last_cell, "  5\t  4\t168.0\t24.1\t 25\n"),
                "has cell X=5 Y=4 outside its 5 columns and 5 rows",
            ),
            (
                data.replace(last_cell, " -1\t  4\t168.0\t24.1\t 25\n"),
                "has cell X=-1 Y=4 outside its 5 columns and 5 rows",
            ),
            (data.replace(dims, "Cols=5\n"), "gives no Rows or Cols before INTENSITY"),
            (
                data.replace(dims, "Cols=5\nRows=-5\n"),
                "has an invalid size of -5 rows and 5 columns",
            ),
            (
                data.replace(dims, f"Cols={2**40}\nRows={2**40}\n"),
                f"has an invalid size of {2**40} rows and {2**40} columns",
            ),
            # A second header could change the size after it had been checked
            (
                data.replace("[MASKS]", "[HEADER]\nCols=2000\nRows=2000\n[MASKS]"),
                "has a HEADER section after INTENSITY",
            ),
            (data.replace(intensity, ""), "has a MASKS section before INTENSITY"),
        ]
        for bad, message in cases:
            with self.subTest(message=message):
                self.assertNotEqual(bad, data)
                with self.assertRaisesRegex(
                    CelFile.ParserError, f"^CEL version 3 file {message}$"
                ):
                    CelFile.read(io.StringIO(bad))

    def testAffy3MarkedCellOutside(self):
        # The cells listed in the MASKS, OUTLIERS and MODIFIED sections were
        # used as numpy indices unchecked: X=5 raised a bare IndexError, and
        # X=-1 silently marked a cell counted from the end of the row.
        with open(self.affy3) as f:
            data = f.read()
        for section, line in [
            ("MASKS", "4 2\n"),
            ("OUTLIERS", "1 3\n"),
            ("MODIFIED", "3 2 21775.0\n"),
        ]:
            self.assertEqual(data.count(line), 1)
            for x in ("5", "-1"):
                bad = data.replace(line, x + line[1:])
                with self.subTest(section=section, x=x):
                    with self.assertRaisesRegex(
                        CelFile.ParserError,
                        f"^CEL version 3 file has cell X={x} Y={line[2]} "
                        "outside its 5 columns and 5 rows$",
                    ):
                        CelFile.read(io.StringIO(bad))
            for y in ("5", "-1"):
                bad = data.replace(line, line[:2] + y + line[3:])
                with self.subTest(section=section, y=y):
                    with self.assertRaisesRegex(
                        CelFile.ParserError,
                        f"^CEL version 3 file has cell X={line[0]} Y={y} "
                        "outside its 5 columns and 5 rows$",
                    ):
                        CelFile.read(io.StringIO(bad))

    @unittest.skipIf(tracemalloc is None, "tracemalloc is not available")
    def testAffy3SizeBeyondCellsMemory(self):
        # The 3 kB example file declaring 2000 by 2000 cells made the old
        # parser allocate 96 MB of arrays; the fixed one stops at the end of
        # the INTENSITY section, having stored only the 29 cell lines.
        with open(self.affy3) as f:
            bad = f.read().replace("Cols=5\nRows=5\n", "Cols=2000\nRows=2000\n")
        was_tracing = tracemalloc.is_tracing()
        if not was_tracing:
            tracemalloc.start()
        try:
            tracemalloc.reset_peak()
            baseline = tracemalloc.get_traced_memory()[0]
            with self.assertRaisesRegex(CelFile.ParserError, "values for 25 cells$"):
                CelFile.read(io.StringIO(bad))
            peak = tracemalloc.get_traced_memory()[1] - baseline
        finally:
            if not was_tracing:
                tracemalloc.stop()
        self.assertLess(peak, 2**20)

    def readV4Example(self, **changes):
        """Return the version 4 example file with some of its sizes changed.

        cellNo and headerLen replace the values in the binary pre-header.
        Rows and Cols replace the values in the text header, and the header
        length in the pre-header changes to match unless headerLen is given.
        """
        with open(self.affy4, "rb") as f:
            data = f.read()
        cellNo, headerLen = struct.unpack_from("<ii", data, 16)
        for key in ("Rows", "Cols"):
            if key in changes:
                old = b"\n%s=5\n" % key.encode()
                new = b"\n%s=%i\n" % (key.encode(), changes[key])
                self.assertIn(old, data)
                data = data.replace(old, new, 1)
                headerLen += len(new) - len(old)
        cellNo = changes.get("cellNo", cellNo)
        headerLen = changes.get("headerLen", headerLen)
        return data[:16] + struct.pack("<ii", cellNo, headerLen) + data[24:]

    def testAffy4BadSizes(self):
        # A version 4 file gives its number of cells and its header length
        # in the binary pre-header, and its rows and columns in the text
        # header. The parser used to trust them: 2**31 - 1 cells in a 1 kB
        # file made it allocate three arrays of 16 GiB each before finding
        # the cells missing, then fail with struct.error. A row count of -1
        # was accepted silently, and a header length of -1 read the rest of
        # the file as header.
        with open(self.affy4, "rb") as f:
            data = f.read()
        dims = "has 5 rows and 5 columns but"
        short = "ends early: expected"
        cases = [
            (self.readV4Example(cellNo=2**31 - 1), f"{dims} 2147483647 cells"),
            (self.readV4Example(cellNo=-1), f"{dims} -1 cells"),
            (
                self.readV4Example(Rows=-5, Cols=-5),
                "has -5 rows and -5 columns but 25 cells",
            ),
            (self.readV4Example(Rows=-1), "has -1 rows and 5 columns but 25 cells"),
            (self.readV4Example(Cols=6), "has 5 rows and 6 columns but 25 cells"),
            (self.readV4Example(headerLen=-1), "has a negative header length"),
            # Rows and columns that agree with the cell count, but no cells
            (
                self.readV4Example(Rows=46340, Cols=46340, cellNo=46340**2),
                f"{short} {10 * 46340**2} bytes of cell data, found 250",
            ),
            (
                self.readV4Example(headerLen=2**31 - 1),
                f"{short} 2147483647 bytes of header, found 1055",
            ),
            (data[:-5], f"{short} 250 bytes of cell data, found 245"),
        ]
        for bad, message in cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex(
                    CelFile.ParserError, f"^CEL version 4 file {message}$"
                ):
                    CelFile.read(io.BytesIO(bad), version=4)

    @unittest.skipIf(tracemalloc is None, "tracemalloc is not available")
    def testAffy4SizesBeyondEndOfFileMemory(self):
        # From a real file, as opposed to io.BytesIO, a single read(n) sets
        # aside n bytes before it reads any. A 1 kB file declaring 20 GiB of
        # cells or a 2 GiB header must not make the parser ask for that much.
        # The old code asked for 68 GiB for the cells (three arrays and the
        # read) and 2 GiB for the header; the fixed code peaks at about
        # 1 MiB, the most it reads at once.
        cases = {
            "cell data": self.readV4Example(Rows=46340, Cols=46340, cellNo=46340**2),
            "header": self.readV4Example(headerLen=2**31 - 1),
        }
        for what, bad in cases.items():
            with self.subTest(what=what), tempfile.TemporaryFile() as handle:
                handle.write(bad)
                handle.seek(0)
                was_tracing = tracemalloc.is_tracing()
                if not was_tracing:
                    tracemalloc.start()
                try:
                    tracemalloc.reset_peak()
                    baseline = tracemalloc.get_traced_memory()[0]
                    with self.assertRaisesRegex(
                        CelFile.ParserError, f"ends early: expected .* of {what},"
                    ):
                        CelFile.read(handle, version=4)
                    peak = tracemalloc.get_traced_memory()[1] - baseline
                finally:
                    if not was_tracing:
                        tracemalloc.stop()
                self.assertLess(peak, 2**24)

    def testAffyWrongModeReadV3(self):
        with self.assertRaises(ValueError):
            with open(self.affy3, "rb") as f:
                record = CelFile.read(f, version=3)

    def testAffyWrongModeReadV4(self):
        with self.assertRaises(ValueError):
            with open(self.affy4) as f:
                record = CelFile.read(f, version=4)

    # Writes a small example Affymetrix V4 CEL File
    def writeExampleV4(self, f, bad=False):
        preHeaders = {
            "cellNo": 25,
            "columns": 5,
            "headerLen": 752,
            "magic": 64,
            "rows": 5,
            "version": 4,
        }
        goodH = {"Axis-invertX": b"0"}
        badH = {"Axis-invertX": b"1"}

        headers = {
            "Algorithm": b"Percentile",
            "AlgorithmParameters": b"Percentile:75;CellMargin:4;Outlie"
            b"rHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize"
            b":TRUE;FullFeatureWidth:7;FullFeatureHeight:7;IgnoreOutlie"
            b"rsInShiftRows:FALSE;FeatureExtraction:TRUE;PoolWidthExten"
            b"stion:1;PoolHeightExtension:1;UseSubgrids:FALSE;Randomize"
            b"Pixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000;NumDATS"
            b"ubgrids:169",
            "AxisInvertY": b"0",
            "Cols": b"5",
            "DatHeader": b"[0..65534]  20_10N:CLS=19420RWS=19420XIN=0"
            b" YIN=0  VE=30        2.0 05/25/05 23:19:07 50102310  M10 "
            b"  \x14  \x14 HuEx-1_0-st-v2.1sq \x14  \x14  \x14  \x14  "
            b"\x14570 \x14 25540.671875 \x14 3.500000 \x14 0.7000 \x14"
            b" 3",
            "GridCornerLL": b"518 18668",
            "GridCornerLR": b"18800 18825",
            "GridCornerUL": b"659 469",
            "GridCornerUR": b"18942 623",
            "OffsetX": b"0",
            "OffsetY": b"0",
            "Rows": b"5",
            "TotalX": b"2560",
            "TotalY": b"2560",
            "swapXY": b"0",
        }
        if not bad:
            headers.update(goodH)
        else:
            headers.update(badH)
        prePadding = b"this text doesn't matter and is ignored\x04"
        preHeadersOrder = ["magic", "version", "columns", "rows", "cellNo", "headerLen"]
        headersEncoded = struct.pack(
            "<" + "i" * len(preHeadersOrder),
            *(preHeaders[header] for header in preHeadersOrder),
        )

        f.write(headersEncoded)
        for header in headers:
            try:
                f.write(
                    bytes(header, encoding="utf-8") + b"=" + headers[header] + b"\n"
                )
            except TypeError:
                f.write(header + b"=" + headers[header] + b"\n")
        f.write(prePadding)
        f.write(b"\x00" * 15)
        for i in range(25):
            f.write(struct.pack("< f f h", i, -i, 9))  # intensity, sdev, pixel


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
