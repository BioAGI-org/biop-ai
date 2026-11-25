# Adapted from Bio.AlignIO.FastaIO copyright 2008-2011 by Peter Cock.
# Copyright 2012 by Wibowo Arindrarto.
# All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
r"""Bio.SearchIO support for Bill Pearson's FASTA tools.

This module adds support for parsing FASTA outputs. FASTA is a suite of
programs that finds regions of local or global similarity between protein
or nucleotide sequences, either by searching databases or identifying
local duplications.

Bio.SearchIO.FastaIO was tested on the following FASTA flavors and versions:

    - flavors: fasta, ssearch, tfastx
    - versions: 35, 36

Other flavors and/or versions may introduce some bugs. Please file a bug report
if you see such problems to Biopython's bug tracker.

More information on FASTA are available through these links:

    - Website: http://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml
    - User guide: http://fasta.bioch.virginia.edu/fasta_www2/fasta_guide.pdf


Supported Formats
=================

Bio.SearchIO.FastaIO supports parsing and indexing FASTA outputs triggered by
the -m 10 flag. Other formats that mimic other programs (e.g. the BLAST tabular
format using the -m 8 flag) may be parseable but using SearchIO's other parsers
(in this case, using the 'blast-tab' parser).


fasta-m10
=========

Note that in FASTA -m 10 outputs, HSPs from different strands are considered to
be from different hits. They are listed as two separate entries in the hit
table. FastaIO recognizes this and will group HSPs with the same hit ID into a
single Hit object, regardless of strand.

FASTA also sometimes output extra sequences adjacent to the HSP match. These
extra sequences are discarded by FastaIO. Only regions containing the actual
sequence match are extracted.

The following object attributes are provided:

+-----------------+-------------------------+----------------------------------+
| Object          | Attribute               | Value                            |
+=================+=========================+==================================+
| QueryResult     | description             | query sequence description       |
|                 +-------------------------+----------------------------------+
|                 | id                      | query sequence ID                |
|                 +-------------------------+----------------------------------+
|                 | program                 | FASTA flavor                     |
|                 +-------------------------+----------------------------------+
|                 | seq_len                 | full length of query sequence    |
|                 +-------------------------+----------------------------------+
|                 | target                  | target search database           |
|                 +-------------------------+----------------------------------+
|                 | version                 | FASTA version                    |
+-----------------+-------------------------+----------------------------------+
| Hit             | seq_len                 | full length of the hit sequence  |
+-----------------+-------------------------+----------------------------------+
| HSP             | bitscore                | \*_bits line                     |
|                 +-------------------------+----------------------------------+
|                 | evalue                  | \*_expect line                   |
|                 +-------------------------+----------------------------------+
|                 | ident_pct               | \*_ident line                    |
|                 +-------------------------+----------------------------------+
|                 | init1_score             | \*_init1 line                    |
|                 +-------------------------+----------------------------------+
|                 | initn_score             | \*_initn line                    |
|                 +-------------------------+----------------------------------+
|                 | opt_score               | \*_opt line, \*_s-w opt line     |
|                 +-------------------------+----------------------------------+
|                 | pos_pct                 | \*_sim line                      |
|                 +-------------------------+----------------------------------+
|                 | sw_score                | \*_score line                    |
|                 +-------------------------+----------------------------------+
|                 | z_score                 | \*_z-score line                  |
+-----------------+-------------------------+----------------------------------+
| HSPFragment     | aln_annotation          | al_cons block, if present        |
| (also via HSP)  +-------------------------+----------------------------------+
|                 | hit                     | hit sequence                     |
|                 +-------------------------+----------------------------------+
|                 | hit_end                 | hit sequence end coordinate      |
|                 +-------------------------+----------------------------------+
|                 | hit_start               | hit sequence start coordinate    |
|                 +-------------------------+----------------------------------+
|                 | hit_strand              | hit sequence strand              |
|                 +-------------------------+----------------------------------+
|                 | query                   | query sequence                   |
|                 +-------------------------+----------------------------------+
|                 | query_end               | query sequence end coordinate    |
|                 +-------------------------+----------------------------------+
|                 | query_start             | query sequence start coordinate  |
|                 +-------------------------+----------------------------------+
|                 | query_strand            | query sequence strand            |
+-----------------+-------------------------+----------------------------------+

"""

import re

from Bio.SearchIO._index import SearchIndexer
from Bio.SearchIO._model import Hit
from Bio.SearchIO._model import HSP
from Bio.SearchIO._model import HSPFragment
from Bio.SearchIO._model import QueryResult

__all__ = ("FastaM10Parser", "FastaM10Indexer")


# precompile regex patterns
# regex for program name
_RE_FLAVS = re.compile(r"t?fast[afmsxy]|pr[sf][sx]|lalign|[gs]?[glso]search")
# regex for sequence ID and length ~ deals with both \n and \r\n
# Use it with .match() at the first ">>>" of the line, and take the description
# as group(2).rstrip(" "). No two parts of it can match the same characters;
# the old pattern's could, so it took minutes to reject a line without the
# "- <length> aa" ending. The ID is one character and then non-whitespace, all
# that the old "(.+?)\s+" could match, and the lookahead stops "\s+" from
# giving back whitespace.
_PTR_ID_DESC_SEQLEN = r">>>(.\S*)\s+(?!\s)(.*?)- (\d+) (?:aa|nt)\s*$"
_RE_ID_DESC_SEQLEN = re.compile(_PTR_ID_DESC_SEQLEN)
_RE_ID_DESC_SEQLEN_IDX = re.compile(_PTR_ID_DESC_SEQLEN.encode())
# regex for qresult, hit, or hsp attribute value
_RE_ATTR = re.compile(r"^; [a-z]+(_[ \w-]+):\s+(.*)$")
# regex for capturing excess start and end sequences in alignments
_RE_START_EXC = re.compile(r"^-*")
# the lookbehind makes re.search try each run of "-" once, not from every "-"
_RE_END_EXC = re.compile(r"(?<!-)-*$")

# attribute name mappings
_HSP_ATTR_MAP = {
    "_initn": ("initn_score", int),
    "_init1": ("init1_score", int),
    "_opt": ("opt_score", int),
    "_s-w opt": ("opt_score", int),
    "_z-score": ("z_score", float),
    "_bits": ("bitscore", float),
    "_expect": ("evalue", float),
    "_score": ("sw_score", int),
    "_ident": ("ident_pct", float),
    "_sim": ("pos_pct", float),
}

# state flags
_STATE_NONE = 0
_STATE_QUERY_BLOCK = 1
_STATE_HIT_BLOCK = 2
_STATE_CONS_BLOCK = 3


def _ends_query(line):
    """Return True if a line ends the alignments of a query (PRIVATE).

    That is an end marker, '>>><<<' or '>>>///', or the next query line,
    such as '  2>>>ID DESCRIPTION - LENGTH aa'. A '>>' hit line is neither,
    even if its description holds '>>>'.
    """
    return line.strip() in (">>><<<", ">>>///") or (
        not line.startswith(">>") and ">>>" in line
    )


def _set_qresult_hits(qresult, hit_rows=()):
    """Append Hits without alignments into QueryResults (PRIVATE)."""
    for hit_row in hit_rows:
        hit_id, remainder = hit_row.split(" ", 1)
        # TODO: parse hit and hsp properties properly; by dealing with:
        #   - any character in the description (brackets, spaces, etc.)
        #   - possible [f] or [r] presence (for frame info)
        #   - possible presence of E2() column
        #   - possible incomplete hit_id due to column length limit
        # The current method only looks at the Hit ID, none of the things above
        if hit_id not in qresult:
            frag = HSPFragment(hit_id, qresult.id)
            hsp = HSP([frag])
            hit = Hit([hsp])
            qresult.append(hit)

    return qresult


def _set_hsp_seqs(hsp, parsed, program):
    """Set HSPs sequences (PRIVATE).

    :param hsp: HSP whose properties will be set
    :type hsp: HSP
    :param parsed: parsed values of the HSP attributes
    :type parsed: dictionary {string: object}
    :param program: program name
    :type program: string

    """
    # a truncated or malformed HSP block may lack some of these
    required = {"seq", "_start", "_stop", "_type"}
    if "tfast" not in program:
        required.add("_display_start")
    for seq_type in ("query", "hit"):
        missing = sorted(required - parsed[seq_type].keys())
        if missing:
            raise ValueError(
                f"Incomplete HSP alignment block for hit {hsp.hit_id!r}:"
                f" {seq_type} block is missing {', '.join(missing)}"
            )
    # the parser collects the aligned sequence lines in lists
    for seq_type in ("hit", "query"):
        parsed[seq_type]["seq"] = "".join(parsed[seq_type]["seq"])
    if "similarity" in hsp.aln_annotation:
        hsp.aln_annotation["similarity"] = "".join(hsp.aln_annotation["similarity"])
    # get aligned sequences and check if they have equal lengths
    start = 0
    for seq_type in ("hit", "query"):
        if "tfast" not in program:
            pseq = parsed[seq_type]
            # adjust start and end coordinates based on the amount of
            # filler characters
            start, stop = _get_aln_slice_coords(pseq)
            start_adj = len(re.search(_RE_START_EXC, pseq["seq"]).group(0))
            stop_adj = len(re.search(_RE_END_EXC, pseq["seq"]).group(0))
            start = start + start_adj
            stop = stop + start_adj - stop_adj
            parsed[seq_type]["seq"] = pseq["seq"][start:stop]
    if len(parsed["query"]["seq"]) != len(parsed["hit"]["seq"]):
        raise ValueError(
            "Length mismatch: %r %r"
            % (len(parsed["query"]["seq"]), len(parsed["hit"]["seq"]))
        )
    if "similarity" in hsp.aln_annotation:
        # only using 'start' since FASTA seems to have trimmed the 'excess'
        # end part
        hsp.aln_annotation["similarity"] = hsp.aln_annotation["similarity"][start:]
        # hit or query works equally well here
        if len(hsp.aln_annotation["similarity"]) != len(parsed["hit"]["seq"]):
            raise ValueError(
                "Length mismatch between the al_cons line and the aligned"
                f" sequences: {len(hsp.aln_annotation['similarity'])}"
                f" vs {len(parsed['hit']['seq'])}"
            )

    # query and hit sequence types must be the same
    if parsed["query"]["_type"] != parsed["hit"]["_type"]:
        raise ValueError(
            f"Query sequence type {parsed['query']['_type']!r} does not match"
            f" hit sequence type {parsed['hit']['_type']!r}"
        )
    type_val = parsed["query"]["_type"]  # hit works fine too
    molecule_type = "DNA" if type_val == "D" else "protein"
    setattr(hsp.fragment, "molecule_type", molecule_type)

    for seq_type in ("hit", "query"):
        # get and set start and end coordinates
        start = int(parsed[seq_type]["_start"])
        end = int(parsed[seq_type]["_stop"])

        setattr(hsp.fragment, seq_type + "_start", min(start, end) - 1)
        setattr(hsp.fragment, seq_type + "_end", max(start, end))
        # set seq and molecule type
        setattr(hsp.fragment, seq_type, parsed[seq_type]["seq"])

        if molecule_type != "protein":
            # get strand from coordinate; start <= end is plus
            # start > end is minus
            if start <= end:
                setattr(hsp.fragment, seq_type + "_strand", 1)
            else:
                setattr(hsp.fragment, seq_type + "_strand", -1)
        else:
            setattr(hsp.fragment, seq_type + "_strand", 0)


def _get_aln_slice_coords(parsed_hsp):
    """Get HSPs sequences (PRIVATE).

    To get the actual pairwise alignment sequences, we must first
    translate the un-gapped sequence based coordinates into positions
    in the gapped sequence (which may have a flanking region shown
    using leading - characters).  To date, I have never seen any
    trailing flanking region shown in the m10 file, but the
    following code should also cope with that.

    Note that this code seems to work fine even when the "sq_offset"
    entries are present as a result of using the -X command line option.
    """
    seq = parsed_hsp["seq"]
    seq_stripped = seq.strip("-")
    disp_start = int(parsed_hsp["_display_start"])
    start = int(parsed_hsp["_start"])
    stop = int(parsed_hsp["_stop"])

    if start <= stop:
        start = start - disp_start
        stop = stop - disp_start + 1
    else:
        start = disp_start - start
        stop = disp_start - stop + 1
    stop += seq_stripped.count("-")
    if not (0 <= start < stop <= len(seq_stripped)):
        raise ValueError(
            "Problem with sequence start/stop,\n%s[%i:%i]\n%s"
            % (seq, start, stop, parsed_hsp)
        )
    return start, stop


class FastaM10Parser:
    """Parser for Bill Pearson's FASTA suite's -m 10 output."""

    def __init__(self, handle, __parse_hit_table=False):
        """Initialize the class."""
        self.handle = handle
        self._preamble = self._parse_preamble()

    def __iter__(self):
        """Iterate over FastaM10Parser object yields query results."""
        for qresult in self._parse_qresult():
            # re-set desc, for hsp query description
            qresult.description = qresult.description
            yield qresult

    def _parse_preamble(self):
        """Parse the Fasta preamble for Fasta flavor and version (PRIVATE)."""
        preamble = {}
        line = self.handle.readline()
        if not line:
            # an empty file holds no query results
            self.line = line
            return preamble
        while True:
            # this should be the line just before the first qresult
            if line.startswith("Query"):
                break
            # try to match for version line
            elif line.startswith(" version"):
                fields = line.split(" ")
                if len(fields) < 3 or not fields[2].strip():
                    raise ValueError(
                        f"Expected a version number on the FASTA version line:\n{line!r}"
                    )
                preamble["version"] = fields[2]
            else:
                # try to match for flavor line
                flav_match = re.match(_RE_FLAVS, line.lower())
                if flav_match:
                    preamble["program"] = flav_match.group(0)
            line = self.handle.readline()
            if not line:
                raise ValueError(
                    "File ended before the 'Query' line of the FASTA preamble"
                )
        self.line = line

        return preamble

    def __parse_hit_table(self):
        """Parse hit table rows."""
        # parse hit table until we see an empty line
        hit_rows = []
        while True:
            line = self.handle.readline()
            if (not line) or line.strip():
                break
            hit_rows.append("")
        self.line = line
        return hit_rows

    def _parse_qresult(self):
        """Parse query result (PRIVATE)."""
        # initial qresult value
        qresult = None
        query_ended = False
        hit_rows = []
        # state values
        state_QRES_NEW = 1
        state_QRES_HITTAB = 3
        state_QRES_CONTENT = 5
        state_QRES_END = 7

        line = self.line
        # an empty file holds no query results
        if not line:
            return

        while True:
            # one line before the hit table
            if line.startswith("The best scores are:"):
                qres_state = state_QRES_HITTAB
            # the end of a query or the file altogether
            elif line.strip() == ">>>///" or not line:
                qres_state = state_QRES_END
            # the beginning of a new query
            elif not line.startswith(">>>") and ">>>" in line:
                qres_state = state_QRES_NEW
            # the beginning of the query info and its hits + hsps
            elif line.startswith(">>>") and line.strip() != ">>><<<":
                qres_state = state_QRES_CONTENT
            # default qres mark
            else:
                qres_state = None

            if qres_state is not None:
                if qres_state == state_QRES_HITTAB:
                    # parse hit table if flag is set
                    hit_rows = self.__parse_hit_table()
                    line = self.handle.readline()

                elif qres_state == state_QRES_END:
                    if qresult is None:
                        raise ValueError(
                            "File ended before the first query line"
                            " '>>>ID DESCRIPTION - LENGTH aa'"
                        )
                    # FASTA 34 and 35 have no '>>>///', but every complete
                    # file has '>>><<<' after the alignments of its last query
                    if not line and not query_ended:
                        raise ValueError(
                            "File ended before the '>>><<<' end marker"
                            f" of query {qresult.id!r}"
                        )
                    yield _set_qresult_hits(qresult, hit_rows)
                    break

                elif qres_state == state_QRES_NEW:
                    # if qresult is filled, yield it first
                    if qresult is not None:
                        yield _set_qresult_hits(qresult, hit_rows)
                    regx = _RE_ID_DESC_SEQLEN.match(line, line.find(">>>"))
                    if regx is None:
                        raise ValueError(
                            "Expected a query line '>>>ID DESCRIPTION - LENGTH aa'"
                            f" (or nt), found:\n{line!r}"
                        )
                    query_id = regx.group(1)
                    seq_len = regx.group(3)
                    desc = regx.group(2).rstrip(" ")
                    qresult = QueryResult(id=query_id)
                    query_ended = False
                    qresult.seq_len = int(seq_len)
                    # get target from the next line
                    line = self.handle.readline()
                    fields = [x for x in line.split(" ") if x]
                    if len(fields) < 2 or not fields[1].strip():
                        raise ValueError(
                            "Expected the target library line after the query line"
                            f" of query {query_id!r}, found:\n{line!r}"
                        )
                    qresult.target = fields[1].strip()
                    if desc is not None:
                        qresult.description = desc
                    # set values from preamble
                    for key, value in self._preamble.items():
                        setattr(qresult, key, value)
                    line = self.handle.readline()

                elif qres_state == state_QRES_CONTENT:
                    if qresult is None:
                        raise ValueError(
                            "Expected a query line '>>>ID DESCRIPTION - LENGTH aa'"
                            f" before the first hit block, found:\n{line!r}"
                        )
                    if not line[3:].startswith(qresult.id):
                        raise ValueError(
                            f"Expected '>>>{qresult.id}' for the current query,"
                            f" found:\n{line!r}"
                        )
                    for hit, strand in self._parse_hit(query_id):
                        # HACK: re-set desc, for hsp hit and query description
                        hit.description = hit.description
                        hit.query_description = qresult.description
                        # if hit is not in qresult, append it
                        if hit.id not in qresult:
                            qresult.append(hit)
                        # otherwise, it might be the same hit with a different strand
                        else:
                            # make sure strand is different and then append hsp to
                            # existing hit
                            for hsp in hit.hsps:
                                assert strand != hsp.query_strand  # noqa: S101
                                qresult[hit.id].append(hsp)
                    line = self.line

            else:
                if line.strip() == ">>><<<":
                    query_ended = True
                line = self.handle.readline()

        self.line = line

    def _parse_hit(self, query_id):
        """Parse hit on query identifier (PRIVATE)."""
        while True:
            line = self.handle.readline()
            if not line:
                raise ValueError(
                    f"File ended before the first '>>' hit line of query {query_id!r}"
                )
            if _ends_query(line):
                raise ValueError(
                    "Expected a '>>' hit line in the alignments of query"
                    f" {query_id!r}, found:\n{line!r}"
                )
            if line.startswith(">>"):
                break

        state = _STATE_NONE
        strand = None
        hsp_list = []
        hsp = None
        parsed_hsp = None
        hit_desc = None
        seq_len = None
        while True:
            # yield hit and create a new one if we're still in the same query
            if line.startswith(">>"):
                # try yielding,  if we have hsps
                if hsp_list:
                    _set_hsp_seqs(hsp, parsed_hsp, self._preamble["program"])
                    hit = Hit(hsp_list)
                    hit.description = hit_desc
                    hit.seq_len = seq_len
                    yield hit, strand
                    hsp_list = []
                # try to get the hit id and desc, and handle cases without descs
                hit_descr = line[2:].strip()
                if not hit_descr:
                    raise ValueError(f"Hit line with no identifier: {line!r}")
                try:
                    hit_id, hit_desc = hit_descr.split(" ", 1)
                except ValueError:
                    hit_id = hit_descr
                    hit_desc = ""
                # create the HSP object for Hit
                frag = HSPFragment(hit_id, query_id)
                hsp = HSP([frag])
                hsp_list.append(hsp)
                # set or reset the state to none
                state = _STATE_NONE
                parsed_hsp = {"query": {}, "hit": {}}
            # create and append a new HSP if line starts with '>--'
            elif line.startswith(">--"):
                # set seq attributes of previous hsp
                _set_hsp_seqs(hsp, parsed_hsp, self._preamble["program"])
                # and create a new one
                frag = HSPFragment(hit_id, query_id)
                hsp = HSP([frag])
                hsp_list.append(hsp)
                # set the state ~ none yet
                state = _STATE_NONE
                parsed_hsp = {"query": {}, "hit": {}}
            # this is either query or hit data in the HSP, depending on the state
            elif line.startswith(">"):
                if state == _STATE_NONE:
                    # make sure it's the correct query
                    if not query_id.startswith(line[1:].split(" ")[0]):
                        raise ValueError(f"{query_id!r} vs {line!r}")
                    state = _STATE_QUERY_BLOCK
                    parsed_hsp["query"]["seq"] = []
                elif state == _STATE_QUERY_BLOCK:
                    # make sure it's the correct hit
                    if not hit_id.startswith(line[1:].split(" ")[0]):
                        raise ValueError(
                            f"Expected the aligned sequence of hit {hit_id!r},"
                            f" found:\n{line!r}"
                        )
                    state = _STATE_HIT_BLOCK
                    parsed_hsp["hit"]["seq"] = []
            # check for conservation block
            elif line.startswith("; al_cons"):
                state = _STATE_CONS_BLOCK
                hsp.fragment.aln_annotation["similarity"] = []
            elif line.startswith(";"):
                # Fasta outputs do not make a clear distinction between Hit
                # and HSPs, so we check the attribute names to determine
                # whether it belongs to a Hit or HSP
                regx = re.search(_RE_ATTR, line.strip())
                if regx is None:
                    raise ValueError("Unexpected line: %r" % line)
                name = regx.group(1)
                value = regx.group(2)

                # for values before the '>...' query block
                if state == _STATE_NONE:
                    if name in _HSP_ATTR_MAP:
                        attr_name, caster = _HSP_ATTR_MAP[name]
                        if caster is not str:
                            value = caster(value)
                        if name in ["_ident", "_sim"]:
                            value *= 100
                        setattr(hsp, attr_name, value)
                # otherwise, pool the values for processing later
                elif state == _STATE_QUERY_BLOCK:
                    parsed_hsp["query"][name] = value
                elif state == _STATE_HIT_BLOCK:
                    if name == "_len":
                        seq_len = int(value)
                    else:
                        parsed_hsp["hit"][name] = value
                # for values in the hit block
                else:
                    raise ValueError("Unexpected line: %r" % line)
            # otherwise, it must be lines containing the sequences
            else:
                if ">" in line:
                    raise ValueError(
                        f"Expected an aligned sequence line, found:\n{line!r}"
                    )
                # if we're in hit, parse into hsp.hit
                if state == _STATE_HIT_BLOCK:
                    parsed_hsp["hit"]["seq"].append(line.strip())
                elif state == _STATE_QUERY_BLOCK:
                    parsed_hsp["query"]["seq"].append(line.strip())
                elif state == _STATE_CONS_BLOCK:
                    hsp.fragment.aln_annotation["similarity"].append(line.strip("\r\n"))
                # we should not get here!
                else:
                    raise ValueError("Unexpected line: %r" % line)
            # yield the last hit if the next line starts a new query or
            # ends the search
            self.line = line = self.handle.readline()
            if not line:
                raise ValueError(
                    f"File ended inside the alignments of query {query_id!r}"
                )
            if _ends_query(line):
                # process HSP alignment and coordinates
                _set_hsp_seqs(hsp, parsed_hsp, self._preamble["program"])
                hit = Hit(hsp_list)
                hit.description = hit_desc
                hit.seq_len = seq_len
                yield hit, strand
                break


class FastaM10Indexer(SearchIndexer):
    """Indexer class for Bill Pearson's FASTA suite's -m 10 output."""

    _parser = FastaM10Parser

    def __init__(self, filename):
        """Initialize the class."""
        SearchIndexer.__init__(self, filename)

    def __iter__(self):
        """Iterate over FastaM10Indexer; yields query results' keys, start offsets, offset lengths."""
        handle = self._handle
        handle.seek(0)
        start_offset = handle.tell()
        qresult_key = None
        query_mark = b">>>"

        line = handle.readline()
        if not line:
            # an empty file holds no query results
            return
        while True:
            end_offset = handle.tell()

            if not line.startswith(b">>") and query_mark in line:
                regx = _RE_ID_DESC_SEQLEN_IDX.match(line, line.find(query_mark))
                if regx is None:
                    raise ValueError(
                        "Expected a query line '>>>ID DESCRIPTION - LENGTH aa'"
                        f" (or nt), found:\n{line!r}"
                    )
                qresult_key = regx.group(1).decode()
                start_offset = end_offset - len(line)
            # yield whenever we encounter a new query or at the end of the file
            if qresult_key is not None:
                if not line:
                    yield qresult_key, start_offset, end_offset - start_offset
                    break
                line = handle.readline()
                if not line.startswith(b">>") and query_mark in line:
                    yield qresult_key, start_offset, end_offset - start_offset
                    start_offset = end_offset
            else:
                line = handle.readline()
                if not line:
                    raise ValueError(
                        "File ended before the first query line of the FASTA output"
                    )

    def get_raw(self, offset):
        """Return the raw record from the file as a bytes string."""
        handle = self._handle
        lines = []
        query_mark = b">>>"

        # read header first
        handle.seek(0)
        line = handle.readline()
        while True:
            lines.append(line)
            line = handle.readline()
            if not line:
                raise ValueError(
                    "File ended before the first query line after the FASTA preamble"
                )
            if not line.startswith(b">>") and query_mark in line:
                break

        # and read the qresult raw string
        handle.seek(offset)
        line = handle.readline()
        while True:
            # preserve whitespace, don't use read_forward
            if not line:
                break
            lines.append(line)

            line = handle.readline()
            # break when we've reached qresult end
            if not line.startswith(b">>") and query_mark in line:
                # FASTA 34 and 35 write '>>><<<' only after the last query, so
                # add one here. Add none at end of file: a last query cut
                # before its '>>><<<' must fail here as it does in parse().
                lines.append(b">>><<<\n")
                break

        return b"".join(lines)


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
