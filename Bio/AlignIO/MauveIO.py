# Copyright 2015-2015 by Eric Rasche.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.AlignIO support for "xmfa" output from Mauve/ProgressiveMauve.

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).

For example, consider a progressiveMauve alignment file containing the following::

    #FormatVersion Mauve1
    #Sequence1File	a.fa
    #Sequence1Entry	1
    #Sequence1Format	FastA
    #Sequence2File	b.fa
    #Sequence2Entry	2
    #Sequence2Format	FastA
    #Sequence3File	c.fa
    #Sequence3Entry	3
    #Sequence3Format	FastA
    #BackboneFile	three.xmfa.bbcols
    > 1:0-0 + a.fa
    --------------------------------------------------------------------------------
    --------------------------------------------------------------------------------
    --------------------------------------------------------------------------------
    > 2:5417-5968 + b.fa
    TTTAAACATCCCTCGGCCCGTCGCCCTTTTATAATAGCAGTACGTGAGAGGAGCGCCCTAAGCTTTGGGAAATTCAAGC-
    --------------------------------------------------------------------------------
    CTGGAACGTACTTGCTGGTTTCGCTACTATTTCAAACAAGTTAGAGGCCGTTACCTCGGGCGAACGTATAAACCATTCTG
    > 3:9476-10076 - c.fa
    TTTAAACACCTTTTTGGATG--GCCCAGTTCGTTCAGTTGTG-GGGAGGAGATCGCCCCAAACGTATGGTGAGTCGGGCG
    TTTCCTATAGCTATAGGACCAATCCACTTACCATACGCCCGGCGTCGCCCAGTCCGGTTCGGTACCCTCCATGACCCACG
    ---------------------------------------------------------AAATGAGGGCCCAGGGTATGCTT
    =
    > 2:5969-6015 + b.fa
    -----------------------
    GGGCGAACGTATAAACCATTCTG
    > 3:9429-9476 - c.fa
    TTCGGTACCCTCCATGACCCACG
    AAATGAGGGCCCAGGGTATGCTT

This is a multiple sequence alignment with multiple aligned sections, so you
would probably load this using the Bio.AlignIO.parse() function:

    >>> from Bio import AlignIO
    >>> align = AlignIO.parse("Mauve/simple_short.xmfa", "mauve")
    >>> alignments = list(align)
    >>> for aln in alignments:
    ...     print(aln)
    ...
    Alignment with 3 rows and 240 columns
    --------------------------------------------...--- a.fa
    TTTAAACATCCCTCGGCCCGTCGCCCTTTTATAATAGCAGTACG...CTG b.fa/5416-5968
    TTTAAACACCTTTTTGGATG--GCCCAGTTCGTTCAGTTGTG-G...CTT c.fa/9475-10076
    Alignment with 2 rows and 46 columns
    -----------------------GGGCGAACGTATAAACCATTCTG b.fa/5968-6015
    TTCGGTACCCTCCATGACCCACGAAATGAGGGCCCAGGGTATGCTT c.fa/9428-9476

Additional information is extracted from the XMFA file and available through
the annotation attribute of each record::

    >>> for record in alignments[0]:
    ...     print(record.id, len(record))
    ...     print("  start: %d, end: %d, strand: %d" %(
    ...         record.annotations['start'], record.annotations['end'],
    ...         record.annotations['strand']))
    ...
    a.fa 240
      start: 0, end: 0, strand: 1
    b.fa/5416-5968 240
      start: 5416, end: 5968, strand: 1
    c.fa/9475-10076 240
      start: 9475, end: 10076, strand: -1

"""

import re

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .Interfaces import AlignmentIterator
from .Interfaces import SequentialAlignmentWriter

XMFA_HEADER_REGEX = re.compile(
    r"> (?P<id>\d+):(?P<start>\d+)-(?P<end>\d+) (?P<strand>[+-]) (?P<name>.*)"
)
XMFA_HEADER_REGEX_BIOPYTHON = re.compile(
    r"> (?P<id>\d+):(?P<start>\d+)-(?P<end>\d+) (?P<strand>[+-]) (?P<name>[^#]*) # (?P<realname>.*)"
)
ID_LINE_FMT = "> {seq_name}:{start}-{end} {strand} {filename} # {ugly_hack}"


def _identifier_split(identifier):
    """Return (name, start, end) string tuple from an identifier (PRIVATE)."""
    id, loc, strand = identifier.split(":")
    start, end = map(int, loc.split("-"))
    start -= 1
    return id, start, end, strand


def _seq_number(record, record_idx):
    """Return the XMFA sequence number of a record, as an int (PRIVATE).

    This is the record name if it is an integer, else its row number (from 1).
    """
    try:
        return int(record.name)
    except ValueError:
        return record_idx + 1


def _snapshot(alignment):
    """Return a copy of what MauveWriter writes from an alignment (PRIVATE).

    Only the id, name, sequence and start/end/strand annotations are kept,
    so an unused annotation that cannot be copied does no harm. Anything
    that is not an alignment is returned as it is, for write_alignment to
    reject.
    """
    if not isinstance(alignment, MultipleSeqAlignment):
        return alignment
    records = []
    for record in alignment:
        seq = record.seq
        if isinstance(seq, MutableSeq):
            seq = Seq(seq)
        annotations = {
            key: record.annotations[key]
            for key in ("start", "end", "strand")
            if key in record.annotations
        }
        records.append(
            SeqRecord(seq, id=record.id, name=record.name, annotations=annotations)
        )
    return MultipleSeqAlignment(records)


class MauveWriter(SequentialAlignmentWriter):
    """Mauve/XMFA alignment writer.

    The first block written lists every sequence in the file, in sequence
    number order, as the Mauve GUI expects. A sequence missing from the first
    alignment gets a ``N:0-0`` placeholder there, which is read back as an
    all-gap row (as with files written by the Mauve GUI).
    """

    def __init__(self, *args, **kwargs):
        """Initialize the class."""
        super().__init__(*args, **kwargs)
        self._wrote_header = False
        self._numbers = None

    def write_file(self, alignments):
        """Write all the alignments, listing every sequence in the first block.

        The sequence numbers of all the alignments are needed before the
        first block is written, so the alignments are read into a list first.
        What is written from each is copied as it is read, since an iterator
        may reuse one alignment object.
        """
        alignments = [_snapshot(alignment) for alignment in alignments]
        self._numbers = {
            _seq_number(record, idx)
            for alignment in alignments
            if isinstance(alignment, MultipleSeqAlignment)
            for idx, record in enumerate(alignment)
        }
        return super().write_file(alignments)

    def write_alignment(self, alignment):
        """Use this to write (another) single alignment to an open file.

        Note that sequences and their annotation are recorded
        together (rather than having a block of annotation followed
        by a block of aligned sequences).
        """
        count = len(alignment)

        self._length_of_sequences = alignment.get_alignment_length()

        # NOTE - For now, the alignment object does not hold any per column
        # or per alignment annotation - only per sequence.

        if count == 0:
            raise ValueError("Must have at least one sequence")
        if self._length_of_sequences == 0:
            raise ValueError("Non-empty sequences are required")

        rows = [
            (_seq_number(record, idx), record) for idx, record in enumerate(alignment)
        ]
        first = not self._wrote_header
        if first:
            self._wrote_header = True
            if self._numbers is None:
                # write_alignment was called directly, not via write_file
                self._numbers = {number for number, record in rows}
            self.handle.write("#FormatVersion Mauve1\n")
            # There are some more headers, but we ignore those for now.
            # Sequence1File	unknown.fa
            # Sequence1Entry	1
            # Sequence1Format	FastA
            for number in sorted(self._numbers):
                self.handle.write(f"#Sequence{number}Entry\t{number}\n")
            # The first LCB we write out is special, and must list ALL
            # sequences, for the Mauve GUI
            # http://darlinglab.org/mauve/user-guide/files.html#non-standard-xmfa-formatting-used-by-the-mauve-gui
            present = {number for number, record in rows}
            rows += [
                (number, None) for number in self._numbers if number not in present
            ]
            rows.sort(key=lambda row: row[0])

        for number, record in rows:
            if record is None:
                id_line = ID_LINE_FMT.format(
                    seq_name=number,
                    start=0,
                    end=0,
                    strand="+",
                    filename=f"{number}.fa",
                    ugly_hack=number,
                )
                self.handle.write(id_line + "\n\n")
            else:
                self._write_record(record, str(number), first)
        self.handle.write("=\n")

    def _write_record(self, record, seq_name, first):
        """Write a single SeqRecord to the file (PRIVATE)."""
        if self._length_of_sequences != len(record.seq):
            raise ValueError("Sequences must all be the same length")

        if (
            "start" in record.annotations
            and "end" in record.annotations
            and "strand" in record.annotations
        ):
            id_line = ID_LINE_FMT.format(
                seq_name=seq_name,
                start=record.annotations["start"] + 1,
                end=record.annotations["end"],
                strand=("+" if record.annotations["strand"] == 1 else "-"),
                filename=record.name + ".fa",
                ugly_hack=record.id,
            )
            lacking_annotations = False
        else:
            id_line = ID_LINE_FMT.format(
                seq_name=seq_name,
                start=0,
                end=0,
                strand="+",
                filename=record.name + ".fa",
                ugly_hack=record.id,
            )
            lacking_annotations = True

        # If the sequence is an empty one, skip writing it out
        if (":0-0 " in id_line or ":1-0 " in id_line) and not lacking_annotations:
            # Except in the first LCB, which lists every sequence
            if first:
                id_line = ID_LINE_FMT.format(
                    seq_name=seq_name,
                    start=0,
                    end=0,
                    strand="+",
                    filename=record.name + ".fa",
                    ugly_hack=record.id,
                )
                id_line = id_line.replace("\n", " ").replace("\r", " ")
                self.handle.write(id_line + "\n\n")
            # Alignments lacking a start/stop/strand were generated by
            # Biopython on load, and shouldn't exist according to XMFA
        else:
            # In other blocks, we only write sequences if they exist in a given
            # alignment.
            id_line = id_line.replace("\n", " ").replace("\r", " ")
            self.handle.write(id_line + "\n")
            for i in range(0, len(record.seq), 80):
                self.handle.write(f"{record.seq[i : i + 80]}\n")


class MauveIterator(AlignmentIterator):
    """Mauve xmfa alignment iterator."""

    def __init__(self, handle, seq_count=None):
        """Create a MauveIterator; the id cache belongs to this iterator."""
        super().__init__(handle, seq_count)
        self._ids: list[str] = []  # for caching IDs between __next__ calls

    def __next__(self):
        """Parse the next alignment from the handle."""
        handle = self.handle
        line = handle.readline()

        if not line:
            raise StopIteration

        # Strip out header comments
        while line and line.strip().startswith("#"):
            line = handle.readline()

        seqs: dict[str, str] = {}
        seq_regions = {}

        latest_id = None
        while True:
            if not line:
                break  # end of file
            line = line.strip()

            if line.startswith("="):
                # There may be more data, but we've reached the end of this
                # alignment
                break
            elif line.startswith(">"):
                m = XMFA_HEADER_REGEX_BIOPYTHON.match(line)
                if not m:
                    m = XMFA_HEADER_REGEX.match(line)
                    if not m:
                        raise ValueError("Malformed header line: %s" % line)

                parsed_id = m.group("id")
                parsed_data = {}
                for key in ("start", "end", "id", "strand", "name", "realname"):
                    try:
                        value = m.group(key)
                        if key == "start":
                            value = int(value)
                            # Convert to zero based counting
                            if value > 0:
                                value -= 1

                        if key == "end":
                            value = int(value)
                        parsed_data[key] = value
                    except IndexError:
                        # This will occur if we're asking for a group that
                        # doesn't exist. It's fine.
                        pass
                seq_regions[parsed_id] = parsed_data

                if parsed_id not in self._ids:
                    self._ids.append(parsed_id)

                seqs.setdefault(parsed_id, "")
                latest_id = parsed_id
            else:
                if latest_id is None:
                    raise ValueError("Saw sequence before definition line")
                seqs[latest_id] += line
            line = handle.readline()

        assert len(seqs) <= len(self._ids)  # noqa: S101

        self.ids = self._ids
        self.sequences = seqs

        if self._ids and seqs:
            alignment_length = max(map(len, list(seqs.values())))
            records = []
            for id in self._ids:
                # Sometimes we don't see a particular sequence in the
                # alignment, so we skip that record since it isn't present in
                # that LCB/alignment
                if id not in seq_regions:
                    continue

                seq = seqs[id] or "-" * alignment_length
                if alignment_length != len(seq):
                    raise ValueError(
                        "Sequences have different lengths, or repeated identifier"
                    )

                if seq_regions[id]["start"] != 0 or seq_regions[id]["end"] != 0:
                    suffix = "/{start}-{end}".format(**seq_regions[id])
                    if "realname" in seq_regions[id]:
                        corrected_id = seq_regions[id]["realname"]
                    else:
                        corrected_id = seq_regions[id]["name"]
                    if corrected_id.count(suffix) == 0:
                        corrected_id += suffix
                else:
                    if "realname" in seq_regions[id]:
                        corrected_id = seq_regions[id]["realname"]
                    else:
                        corrected_id = seq_regions[id]["name"]

                record = SeqRecord(Seq(seq), id=corrected_id, name=id)

                record.annotations["start"] = seq_regions[id]["start"]
                record.annotations["end"] = seq_regions[id]["end"]
                record.annotations["strand"] = (
                    1 if seq_regions[id]["strand"] == "+" else -1
                )

                records.append(record)
            if (
                self.records_per_alignment is not None
                and self.records_per_alignment != len(records)
            ):
                raise ValueError(
                    "Found %i records in this alignment, told to expect %i"
                    % (len(records), self.records_per_alignment)
                )
            return MultipleSeqAlignment(records)
        else:
            raise StopIteration
