# Copyright 2004 by Harry Zuzan. All rights reserved.
# Copyright 2016 by Adam Kurkiewicz. All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Reading information from Affymetrix CEL files version 3 and 4."""

import struct
import sys
from array import array

try:
    import numpy as np
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.Affy.CelFile"
    ) from None


class ParserError(ValueError):
    """Affymetrix parser error."""

    def __init__(self, *args):
        """Initialise class."""
        super().__init__(*args)


class Record:
    """Stores the information in a cel file.

    Example usage:

    >>> from Bio.Affy import CelFile
    >>> with open("Affy/affy_v3_example.CEL") as handle:
    ...     c = CelFile.read(handle)
    ...
    >>> print(c.ncols, c.nrows)
    5 5
    >>> print(c.intensities)
    [[  234.   170. 22177.   164. 22104.]
     [  188.   188. 21871.   168. 21883.]
     [  188.   193. 21455.   198. 21300.]
     [  188.   182. 21438.   188. 20945.]
     [  193. 20370.   174. 20605.   168.]]
    >>> print(c.stdevs)
    [[  24.    34.5 2669.    19.7 3661.2]
     [  29.8   29.8 2795.9   67.9 2792.4]
     [  29.8   88.7 2976.5   62.  2914.5]
     [  29.8   76.2 2759.5   49.2 2762. ]
     [  38.8 2611.8   26.6 2810.7   24.1]]
    >>> print(c.npix)
    [[25 25 25 25 25]
     [25 25 25 25 25]
     [25 25 25 25 25]
     [25 25 25 25 25]
     [25 25 25 25 25]]

    """

    def __init__(self):
        """Initialize the class."""
        self.version = None
        self.GridCornerUL = None
        self.GridCornerUR = None
        self.GridCornerLR = None
        self.GridCornerLL = None
        self.DatHeader = None
        self.Algorithm = None
        self.AlgorithmParameters = None
        self.NumberCells = None
        self.intensities = None
        self.stdevs = None
        self.npix = None
        self.nrows = None
        self.ncols = None
        self.nmask = None
        self.mask = None
        self.noutliers = None
        self.outliers = None
        self.modified = None


def read(handle, version=None):
    """Read Affymetrix CEL file and return Record object.

    CEL files format versions 3 and 4 are supported. The format is
    detected from the file contents. If you pass the version argument (3 or
    4), a ValueError is raised if the file is in the other format.

    The Record object returned by this function stores the intensities from
    the CEL file in record.intensities.
    Currently, record.mask and record.outliers are not set in when parsing
    version 4 CEL files.

    Example Usage:

    >>> from Bio.Affy import CelFile
    >>> with open("Affy/affy_v3_example.CEL") as handle:
    ...     record = CelFile.read(handle)
    ...
    >>> record.version == 3
    True
    >>> print("%i by %i array" % record.intensities.shape)
    5 by 5 array

    >>> with open("Affy/affy_v4_example.CEL", "rb") as handle:
    ...     record = CelFile.read(handle, version=4)
    ...
    >>> record.version == 4
    True
    >>> print("%i by %i array" % record.intensities.shape)
    5 by 5 array

    """
    if version not in (None, 3, 4):
        raise ValueError(f"version should be None, 3 or 4, not {version!r}")
    try:
        data = handle.read(0)
    except AttributeError:
        raise ValueError("handle should be a file handle") from None
    data = handle.read(4)
    if not data:
        raise ValueError("Empty file.")
    if data == b"[CEL":
        raise ValueError("CEL file in version 3 format should be opened in text mode")
    if data == "[CEL":
        # Version 3 format. Continue to read the header here before passing
        # control to _read_v3 to avoid having to seek to the beginning of
        # the file. Use readline rather than next, as some file-like objects
        # such as tempfile.NamedTemporaryFile are not iterators (issue #4926).
        data += handle.readline()
        if data.strip() != "[CEL]":
            raise ValueError("Failed to parse Affy Version 3 CEL file.")
        line = handle.readline()
        if "=" not in line:
            raise ValueError("Failed to parse Affy Version 3 CEL file.")
        keyword, value = line.split("=", 1)
        if keyword != "Version":
            raise ValueError("Failed to parse Affy Version 3 CEL file.")
        if int(value) != 3:
            raise ValueError("Incorrect version number in Affy Version 3 CEL file.")
        if version == 4:
            raise ValueError("expected CEL version 4, but the file is version 3")
        return _read_v3(handle)
    try:
        magicNumber = struct.unpack("<i", data)
    except TypeError:
        raise ValueError(
            "CEL file in version 4 format should be opened in binary mode"
        ) from None
    except struct.error:
        raise ValueError(
            "Failed to read magic number from Affy Version 4 CEL file"
        ) from None
    if magicNumber != (64,):
        raise ValueError("Incorrect magic number in Affy Version 4 CEL file")
    if version == 3:
        raise ValueError("expected CEL version 3, but the file is version 4")
    return _read_v4(handle)


def _read_v4_bytes(f, size, what):
    # Read exactly size bytes, at most 1 MiB at a time. On a file handle a
    # single f.read(size) sets aside size bytes before reading any, so a
    # small file declaring a huge size would make it ask for gigabytes.
    chunks = []
    remaining = size
    while remaining > 0:
        chunk = f.read(min(remaining, 1 << 20))
        if not chunk:
            raise ParserError(
                f"CEL version 4 file ends early: expected {size} bytes of "
                f"{what}, found {size - remaining}"
            )
        chunks.append(chunk)
        remaining -= len(chunk)
    return b"".join(chunks)


def _read_v4(f):
    # We follow the documentation here:
    # http://www.affymetrix.com/estore/support/developer/powertools/changelog/gcos-agcc/cel.html.affx
    record = Record()
    preHeaders = ["version", "columns", "rows", "cellNo", "headerLen"]
    preHeadersMap = {}
    headersMap = {}

    # Load pre-headers. The magic number was already parsed in the read
    # function calling _read_v4.
    preHeadersMap["magic"] = 64
    try:
        for name in preHeaders:
            preHeadersMap[name] = struct.unpack("<i", f.read(4))[0]
    except struct.error:
        raise ParserError("Failed to parse CEL version 4 file") from None

    if preHeadersMap["headerLen"] < 0:
        raise ParserError("CEL version 4 file has a negative header length")
    char = _read_v4_bytes(f, preHeadersMap["headerLen"], "header")
    header = char.decode("ascii", "ignore")
    for line in header.split("\n"):
        if "=" in line:
            headline = line.split("=")
            headersMap[headline[0]] = "=".join(headline[1:])

    record.version = preHeadersMap["version"]
    if record.version != 4:
        raise ParserError("Incorrect version number in CEL version 4 file")

    record.GridCornerUL = headersMap["GridCornerUL"]
    record.GridCornerUR = headersMap["GridCornerUR"]
    record.GridCornerLR = headersMap["GridCornerLR"]
    record.GridCornerLL = headersMap["GridCornerLL"]
    record.DatHeader = headersMap["DatHeader"]
    record.Algorithm = headersMap["Algorithm"]
    record.AlgorithmParameters = headersMap["AlgorithmParameters"]
    record.NumberCells = preHeadersMap["cellNo"]
    # record.intensities are set below
    # record.stdevs are set below
    # record.npix are set below
    record.nrows = int(headersMap["Rows"])
    record.ncols = int(headersMap["Cols"])
    if not (
        record.nrows >= 0
        and record.ncols >= 0
        and record.nrows * record.ncols == record.NumberCells
    ):
        raise ParserError(
            f"CEL version 4 file has {record.nrows} rows and {record.ncols} "
            f"columns but {record.NumberCells} cells"
        )

    # These cannot be reliably set in v4, because of discrepancies between real
    # data and the documented format.
    record.nmask = None
    record.mask = None
    record.noutliers = None
    record.outliers = None
    record.modified = None

    # Real data never seems to have anything but zeros here, but we don't want
    # to take chances. Raising an error is better than returning unreliable
    # data.
    def raiseBadHeader(field, expected):
        actual = int(headersMap[field])
        message = f"The header {field} is expected to be 0, not {actual}"
        if actual != expected:
            raise ParserError(message)

    raiseBadHeader("Axis-invertX", 0)

    raiseBadHeader("AxisInvertY", 0)

    raiseBadHeader("OffsetX", 0)

    raiseBadHeader("OffsetY", 0)

    # This is unfortunately undocumented, but it turns out that real data has
    # the record.AlgorithmParameters repeated in the data section, until an
    # EOF, i.e. b"\x04". Give up at the end of the file or after 10,000 bytes.
    for _ in range(10**4):
        char = f.read(1)
        if char in (b"\x04", b""):
            break
    if char != b"\x04":
        raise ParserError(
            "CEL version 4 file has no b'\\x04' ending the short, "
            "undocumented binary blob after the header"
        )

    # After that there are precisely 15 bytes padded. Again, undocumented.
    padding = f.read(15)

    # That's how we pull out the values (triplets of the form float, float,
    # signed short).
    structa = struct.Struct("< f f h")

    # There are 10 bytes in our struct.
    structSize = 10

    # Read the cells before allocating arrays for them, so that a file
    # declaring more cells than it holds fails before asking for memory.
    b = _read_v4_bytes(f, structSize * record.NumberCells, "cell data")

    # We initialize the most important: intensities, stdevs and npixs.
    record.intensities = np.empty(record.NumberCells, dtype=float)
    record.stdevs = np.empty(record.NumberCells, dtype=float)
    record.npix = np.empty(record.NumberCells, dtype=int)

    for i in range(record.NumberCells):
        binaryFragment = b[i * structSize : (i + 1) * structSize]
        intensity, stdevs, npix = structa.unpack(binaryFragment)
        record.intensities[i] = intensity
        record.stdevs[i] = stdevs
        record.npix[i] = npix

    # reshape without copying.
    record.intensities = record.intensities.reshape((record.nrows, record.ncols))
    record.stdevs = record.stdevs.reshape((record.nrows, record.ncols))
    record.npix = record.npix.reshape((record.nrows, record.ncols))

    return record


def _read_v3_cells(record, cells):
    # Make the intensity arrays from the cells stored while reading the
    # INTENSITY section. They are made only once the cells are known to cover
    # the whole array, so that a short file declaring large Rows and Cols
    # cannot make the parser allocate memory for cells it does not hold.
    index, intensities, stdevs, npix = cells
    shape = (record.nrows, record.ncols)
    # Usually each cell is listed once, in order. As every flat index lies in
    # [0, nrows * ncols), that many strictly increasing ones must be exactly
    # 0, 1, 2, ...; then use the stored values as they are, without copying.
    flat = np.asarray(index)
    if len(flat) == record.nrows * record.ncols and (flat[1:] > flat[:-1]).all():
        record.intensities = np.asarray(intensities, float).reshape(shape)
        record.stdevs = np.asarray(stdevs, float).reshape(shape)
        record.npix = np.asarray(npix, int).reshape(shape)
        return
    # A cell may be listed more than once; as before, the last line wins.
    # np.unique gives the first line of each distinct cell, in cell order;
    # on the lines reversed, that is the last.
    last = np.unique(flat[::-1], return_index=True)[1]
    if len(last) < record.nrows * record.ncols:
        raise ParserError(
            f"CEL version 3 file has {record.nrows} rows and {record.ncols} "
            f"columns but values for {len(last)} cells"
        )
    last = len(index) - 1 - last
    record.intensities = np.asarray(intensities, float)[last].reshape(shape)
    record.stdevs = np.asarray(stdevs, float)[last].reshape(shape)
    record.npix = np.asarray(npix, int)[last].reshape(shape)


def _read_v3_cell(record, words):
    # Return the row and column of the cell on a line of the INTENSITY, MASKS,
    # OUTLIERS or MODIFIED section. The file gives the column (X) first.
    # Check both, as numpy would count a negative index from the end.
    y = int(words[0])
    x = int(words[1])
    if not (0 <= x < record.nrows and 0 <= y < record.ncols):
        raise ParserError(
            f"CEL version 3 file has cell X={y} Y={x} outside "
            f"its {record.ncols} columns and {record.nrows} rows"
        )
    return x, y


def _read_v3_shape(record, section):
    # MASKS, OUTLIERS and MODIFIED mark cells of the intensity array, whose
    # size has been checked against the cells in the file.
    if record.intensities is None:
        raise ParserError(
            f"CEL version 3 file has a {section} section before INTENSITY"
        )
    return record.intensities.shape


def _read_v3(handle):
    # Needs error handling.
    # Needs to know the chip design.
    record = Record()
    # The version number was already obtained when the read function calling
    # _read_v3 parsed the CEL section.
    record.version = 3
    section = ""
    cells = None
    for line in handle:
        line = line.rstrip("\r\n")
        if not line:
            continue
        if section == "INTENSITY" and line.startswith("["):
            _read_v3_cells(record, cells)
        # Set current section
        if line.startswith("[HEADER]"):
            if record.intensities is not None:
                # Rows and Cols must not change once checked against the cells
                raise ParserError(
                    "CEL version 3 file has a HEADER section after INTENSITY"
                )
            section = "HEADER"
        elif line.startswith("[INTENSITY]"):
            section = "INTENSITY"
            if record.nrows is None or record.ncols is None:
                raise ParserError(
                    "CEL version 3 file gives no Rows or Cols before INTENSITY"
                )
            # No file holds more than sys.maxsize cells, and the limit keeps
            # the flat index of each cell within a 64-bit integer.
            if not (
                record.nrows >= 0
                and record.ncols >= 0
                and record.nrows * record.ncols <= sys.maxsize
            ):
                raise ParserError(
                    f"CEL version 3 file has an invalid size of {record.nrows} "
                    f"rows and {record.ncols} columns"
                )
            # The flat index, intensity, stdev and pixel count of each line,
            # until the section ends; see _read_v3_cells.
            cells = (array("q"), array("d"), array("d"), array("q"))
        elif line.startswith("[MASKS]"):
            section = "MASKS"
            record.mask = np.zeros(_read_v3_shape(record, section), bool)
        elif line.startswith("[OUTLIERS]"):
            section = "OUTLIERS"
            record.outliers = np.zeros(_read_v3_shape(record, section), bool)
        elif line.startswith("[MODIFIED]"):
            section = "MODIFIED"
            record.modified = np.zeros(_read_v3_shape(record, section))
        elif line.startswith("["):
            raise ParserError("Unknown section found in version 3 CEL file")
        else:  # read the data in a section
            if section == "HEADER":
                # Set record.ncols and record.nrows, remaining data goes into
                # record.header dict
                key, value = line.split("=", 1)
                if key == "Cols":
                    record.ncols = int(value)
                elif key == "Rows":
                    record.nrows = int(value)
                elif key == "GridCornerUL":
                    x, y = value.split()
                    record.GridCornerUL = (int(x), int(y))
                elif key == "GridCornerUR":
                    x, y = value.split()
                    record.GridCornerUR = (int(x), int(y))
                elif key == "GridCornerLR":
                    x, y = value.split()
                    record.GridCornerLR = (int(x), int(y))
                elif key == "GridCornerLL":
                    x, y = value.split()
                    record.GridCornerLL = (int(x), int(y))
                elif key == "DatHeader":
                    # not sure if all parameters here are interpreted correctly
                    record.DatHeader = {}
                    i = value.find(":")
                    if i >= 0:
                        min_max_pixel_intensity, filename = value[:i].split()
                        record.DatHeader["filename"] = filename
                        if not (
                            min_max_pixel_intensity.startswith("[")
                            and min_max_pixel_intensity.endswith("]")
                        ):
                            raise ParserError(
                                "Expected [min..max] pixel intensity in DatHeader, "
                                f"found {min_max_pixel_intensity!r}"
                            )
                        (
                            min_pixel_intensity,
                            max_pixel_intensity,
                        ) = min_max_pixel_intensity[1:-1].split("..")
                        record.DatHeader["min-pixel_intensity"] = int(
                            min_pixel_intensity
                        )
                        record.DatHeader["max-pixel_intensity"] = int(
                            max_pixel_intensity
                        )
                        value = value[i + 1 :]
                        index = 0
                        field = value[index : index + 9]
                        if field[:4] != "CLS=" or field[8] != " ":
                            raise ValueError(
                                "Field does not start with 'CLS=' or have a blank space at position 8"
                            )
                        record.DatHeader["CLS"] = int(field[4:8])
                        index += 9
                        field = value[index : index + 9]
                        if field[:4] != "RWS=" or field[8] != " ":
                            raise ValueError(
                                "Field does not start with 'RWS=' or have a blank space at position 8"
                            )
                        record.DatHeader["RWS"] = int(field[4:8])
                        index += 9
                        field = value[index : index + 7]
                        if field[:4] != "XIN=" or field[6] != " ":
                            raise ValueError(
                                "Field does not start with 'XIN=' or have a blank space at position 6"
                            )
                        record.DatHeader["XIN"] = int(field[4:6])
                        index += 7
                        field = value[index : index + 7]
                        if field[:4] != "YIN=" or field[6] != " ":
                            raise ValueError(
                                "Field does not start with 'YIN=' or have a blank space at poition 6"
                            )
                        record.DatHeader["YIN"] = int(field[4:6])
                        index += 7
                        field = value[index : index + 6]
                        if field[:3] != "VE=" or field[5] != " ":
                            raise ValueError(
                                "Field does not start with 'VE=' or have a blank space at position 5"
                            )
                        record.DatHeader["VE"] = int(field[3:5])
                        index += 6
                        field = value[index : index + 7]
                        if field[6] != " ":
                            raise ValueError(
                                "Field value for position 6 isn't a blank space"
                            )
                        temperature = field[:6].strip()
                        if temperature:
                            record.DatHeader["temperature"] = int(temperature)
                        else:
                            record.DatHeader["temperature"] = None
                        index += 7
                        field = value[index : index + 4]
                        if not field.endswith(" "):
                            raise ValueError("Field doesn't end with a blank space")
                        record.DatHeader["laser-power"] = float(field)
                        index += 4
                        field = value[index : index + 18]
                        if field[8] != " ":
                            raise ValueError(
                                "Field value for position 8 isn't a blank space"
                            )
                        record.DatHeader["scan-date"] = field[:8]
                        if field[17] != " ":
                            raise ValueError(
                                "Field value for position 17 isn't a blank space"
                            )
                        record.DatHeader["scan-date"] = field[:8]
                        record.DatHeader["scan-time"] = field[9:17]
                        index += 18
                        value = value[index:]
                    subfields = value.split("\x14")
                    if len(subfields) != 12:
                        raise ParserError(
                            "Expected 12 subfields in DatHeader, "
                            f"found {len(subfields)}"
                        )
                    subfield = subfields[0]
                    try:
                        scanner_id, scanner_type = subfield.split()
                    except ValueError:
                        scanner_id = subfield.strip()
                    else:
                        record.DatHeader["scanner-type"] = scanner_type
                    record.DatHeader["scanner-id"] = scanner_id
                    record.DatHeader["array-type"] = subfields[2].strip()
                    field = subfields[7].strip()
                    if field:
                        record.DatHeader["filter-wavelength"] = int(field)
                    field = subfields[8].strip()
                    if field:
                        record.DatHeader["arc-radius"] = float(field)
                    field = subfields[9].strip()
                    if field:
                        record.DatHeader["laser-spotsize"] = float(field)
                    field = subfields[10].strip()
                    if field:
                        record.DatHeader["pixel-size"] = float(field)
                    field = subfields[11].strip()
                    if field:
                        record.DatHeader["image-orientation"] = int(field)
                elif key == "Algorithm":
                    record.Algorithm = value
                elif key == "AlgorithmParameters":
                    parameters = value.split(";")
                    values = {}
                    for parameter in parameters:
                        key, value = parameter.split(":", 1)
                        if key in (
                            "Percentile",
                            "CellMargin",
                            "FullFeatureWidth",
                            "FullFeatureHeight",
                            "PoolWidthExtenstion",
                            "PoolHeightExtension",
                            "NumPixelsToUse",
                            "ExtendPoolWidth",
                            "ExtendPoolHeight",
                            "OutlierRatioLowPercentile",
                            "OutlierRatioHighPercentile",
                            "HalfCellRowsDivisor",
                            "HalfCellRowsRemainder",
                            "HighCutoff",
                            "LowCutoff",
                            "featureRows",
                            "featureColumns",
                        ):
                            values[key] = int(value)
                        elif key in (
                            "OutlierHigh",
                            "OutlierLow",
                            "StdMult",
                            "PercentileSpread",
                            "PairCutoff",
                            "featureWidth",
                            "featureHeight",
                        ):
                            values[key] = float(value)
                        elif key in (
                            "FixedCellSize",
                            "IgnoreOutliersInShiftRows",
                            "FeatureExtraction",
                            "UseSubgrids",
                            "RandomizePixels",
                            "ImageCalibration",
                            "IgnoreShiftRowOutliers",
                        ):
                            if value == "TRUE":
                                value = True
                            elif value == "FALSE":
                                value = False
                            else:
                                raise ValueError("Unexpected boolean value")
                            values[key] = value
                        elif key in (
                            "AlgVersion",
                            "ErrorBasis",
                            "CellIntensityCalculationType",
                        ):
                            values[key] = value
                        else:
                            raise ValueError("Unexpected tag in AlgorithmParameters")
                    record.AlgorithmParameters = values
            elif section == "INTENSITY":
                if line.startswith("NumberCells="):
                    key, value = line.split("=", 1)
                    record.NumberCells = int(value)
                elif line.startswith("CellHeader="):
                    key, value = line.split("=", 1)
                    if value.split() != ["X", "Y", "MEAN", "STDV", "NPIXELS"]:
                        raise ParserError(
                            "Unexpected CellHeader in INTENSITY "
                            "section CEL version 3 file"
                        )
                else:
                    words = line.split()
                    x, y = _read_v3_cell(record, words)
                    cells[0].append(x * record.ncols + y)
                    cells[1].append(float(words[2]))
                    cells[2].append(float(words[3]))
                    cells[3].append(int(words[4]))
            elif section == "MASKS":
                if line.startswith("NumberCells="):
                    key, value = line.split("=", 1)
                    record.nmask = int(value)
                elif line.startswith("CellHeader="):
                    key, value = line.split("=", 1)
                    if value.split() != ["X", "Y"]:
                        raise ParserError(
                            "Unexpected CellHeader in MASKS "
                            "section in CEL version 3 file"
                        )
                else:
                    words = line.split()
                    x, y = _read_v3_cell(record, words)
                    record.mask[x, y] = True
            elif section == "OUTLIERS":
                if line.startswith("NumberCells="):
                    key, value = line.split("=", 1)
                    record.noutliers = int(value)
                elif line.startswith("CellHeader="):
                    key, value = line.split("=", 1)
                    if value.split() != ["X", "Y"]:
                        raise ParserError(
                            "Unexpected CellHeader in OUTLIERS "
                            "section in CEL version 3 file"
                        )
                else:
                    words = line.split()
                    x, y = _read_v3_cell(record, words)
                    record.outliers[x, y] = True
            elif section == "MODIFIED":
                if line.startswith("NumberCells="):
                    key, value = line.split("=", 1)
                    record.nmodified = int(value)
                elif line.startswith("CellHeader="):
                    key, value = line.split("=", 1)
                    if value.split() != ["X", "Y", "ORIGMEAN"]:
                        raise ParserError(
                            "Unexpected CellHeader in MODIFIED "
                            "section in CEL version 3 file"
                        )
                else:
                    words = line.split()
                    x, y = _read_v3_cell(record, words)
                    record.modified[x, y] = float(words[2])
    if section == "INTENSITY":
        _read_v3_cells(record, cells)
    if record.intensities is None:
        raise ParserError("CEL version 3 file has no INTENSITY section")
    return record


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
