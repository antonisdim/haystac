#!/usr/bin/env python
#
# Based on 'FilterUniqueBAM' by
#   Martin Kircher
#   Martin.Kircher@eva.mpg.de
#
"""paleomix rmdup_collapsed [options] < sorted.bam > out.bam
The rmdup_collapsed filters a BAM file for PCR duplicates unpaired reads under
the assumption that any unpaired read have been generated by the merging of
overlapping paired-end reads, and thereby represent the complete template
sequence. PCR duplicates are therefore detected based on both the 5' and 3'
alignment coordinate.
Paired reads (0x1), unmapped reads (0x4), secondary alignments (0x100),
reads that failed QC (0x200), and chimeric alignments (0x800), as identified
using the BAM record flags, are not filtered, but simply written to the output.
By default, filtered reads are flagged using the "duplicate" flag (0x400), and
written to the output. Use the --remove-duplicates command-line option to
instead remove these records from the output.
"""
import collections
import random
import sys

from argparse import ArgumentParser

import pysam


_FILTERED_FLAGS = 0x1  # PE reads
_FILTERED_FLAGS |= 0x4  # Unmapped
_FILTERED_FLAGS |= 0x100  # Secondary alignment
_FILTERED_FLAGS |= 0x200  # Failed QC
_FILTERED_FLAGS |= 0x800  # Chimeric alignment

_CIGAR_SOFTCLIP = 4
_CIGAR_HARDCLIP = 5


def read_quality(read):
    qualities = read.query_alignment_qualities
    if qualities is None:
        # Generate value in range (-1; 0]
        return -random.random()

    return sum(qualities)


def copy_number(read):
    # has_tag is faster than try/except, since most reads lack the tag.
    if read.has_tag("XP"):
        return read.get_tag("XP")

    return 0


def mark_duplicate_reads(reads):
    """Identifies the best read from a set of PCR duplicates, and marks all
    other reads as duplicates. The best read is selected by quality, among the
    reads sharing the most common CIGAR string.
    """
    by_cigar = collections.defaultdict(list)
    for read in reads:
        key = tuple(read.cigartuples)
        by_cigar[key].append(read)

    # Select the most common CIGAR strings, favoring simple CIGARs
    best_count, best_cigar_len = max(
        (len(values), -len(cigar)) for cigar, values in by_cigar.items()
    )
    best_cigar_len = -best_cigar_len

    best_read = None
    best_quality = -1
    copies = len(reads)

    for cigar, candidates in by_cigar.items():
        if len(cigar) == best_cigar_len and len(candidates) == best_count:
            for read in candidates:
                copies += copy_number(read)
                quality = read_quality(read)

                if quality > best_quality:
                    best_read = read
                    best_quality = quality
        else:
            copies += sum(copy_number(read) for read in candidates)

    best_read.set_tag("XP", copies, "i")
    for read in reads:
        read.is_duplicate = read is not best_read


def write_read(args, out, read_and_alignment, duplicates_by_alignment):
    read, alignment = read_and_alignment
    if alignment is not None:
        duplicates = duplicates_by_alignment.pop(alignment)

        if len(duplicates) > 1:
            # Select the best read and mark the others as duplicates.
            mark_duplicate_reads(duplicates)
        else:
            duplicates[0].is_duplicate = False
            duplicates[0].set_tag("XP", 1, "i")

    if not (args.remove_duplicates and read.is_duplicate):
        out.write(read)


def can_write_read(read_and_alignment, current_position):
    """Returns true if the first read in the cache can safely be written. This
    will be the case if the read was not the first in a set of reads with the
    same alignment, or if the current position has gone beyond the last base
    covered in that alignment.
    """
    _, alignment = read_and_alignment
    if alignment is None:
        return True

    current_ref_id, current_ref_start = current_position
    alignment_ref_id, _, _, alignment_ref_end = alignment

    return alignment_ref_id != current_ref_id or alignment_ref_end < current_ref_start


def clipped_bases_at_front(cigartuples):
    """Returns number of bases soft or hard clipped at start of the CIGAR."""
    total = 0
    for (operation, length) in cigartuples:
        if operation != _CIGAR_SOFTCLIP and operation != _CIGAR_HARDCLIP:
            break

        total += length

    return total


def unclipped_alignment_coordinates(read):
    """Returns tuple describing the alignment, with external coordinates
    modified to account for clipped bases, assuming an ungapped alignment to
    the reference. This is equivalent to the behavior of Picard MarkDuplicates.
    """
    cigartuples = read.cigartuples
    start = read.reference_start - clipped_bases_at_front(cigartuples)
    end = read.reference_end + clipped_bases_at_front(reversed(cigartuples))

    return (read.reference_id, read.is_reverse, start, end)


def process_aligned_read(cache, duplicates_by_alignment, read):
    """Processes an aligned read, either pairing it with an existing read, or
    creating a new alignment block to track copies of this copies.
    """
    alignment = unclipped_alignment_coordinates(read)

    try:
        duplicates_by_alignment[alignment].append(read)
        cache.append((read, None))
    except KeyError:
        # No previous reads with matching alignment; this read will
        # serve to track any other reads with the same alignment.
        duplicates_by_alignment[alignment] = [read]
        cache.append((read, alignment))


def is_trailing_unmapped_read(read):
    return read.is_unmapped and read.reference_id == -1 and read.reference_start == -1


def process(args, infile, outfile):
    cache = collections.deque()
    duplicates_by_alignment = {}
    last_position = (0, 0)
    read_num = 1

    for read_num, read in enumerate(infile, start=read_num):
        current_position = (read.reference_id, read.reference_start)
        if last_position > current_position:
            # Check also catches trailing unmapped reads mapped to (-1, -1).
            if not is_trailing_unmapped_read(read):
                sys.stderr.write(
                    "ERROR: Input file is not sorted by "
                    "coordinates at read %i. Aborting!\n" % (read_num,)
                )
                return 1

            cache.append((read, None))
            break
        elif read.flag & _FILTERED_FLAGS:
            cache.append((read, None))
        else:
            process_aligned_read(cache, duplicates_by_alignment, read)

        last_position = current_position
        while cache and can_write_read(cache[0], current_position):
            write_read(args, outfile, cache.popleft(), duplicates_by_alignment)

    while cache:
        write_read(args, outfile, cache.popleft(), duplicates_by_alignment)

    assert not duplicates_by_alignment, duplicates_by_alignment

    for read_num, read in enumerate(infile, start=read_num + 1):
        if not is_trailing_unmapped_read(read):
            sys.stderr.write(
                "ERROR: Input file is not sorted by "
                "coordinates at read %i. Aborting!\n" % (read_num,)
            )
            return 1

        outfile.write(read)

    return 0


def parse_args(argv):
    parser = ArgumentParser(usage=__doc__)
    parser.add_argument(
        "input",
        default="-",
        nargs="?",
        help="BAM file; if not set, input is read from STDIN.",
    )
    parser.add_argument(
        "--remove-duplicates",
        help="Remove duplicates from output; by default "
        "duplicates are only flagged (flag = 0x400).",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--seed",
        default=None,
        type=int,
        help="Seed used for randomly selecting representative "
        "reads when no reads have quality scores assigned"
        "[default: initialized using system time].",
    )

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)

    # Initialize seed used when selecting among reads without quality scores
    random.seed(args.seed)

    if args.input == "-" and sys.stdin.isatty():
        sys.stderr.write("STDIN is a terminal, terminating!\n")
        return 1
    elif sys.stdout.isatty():
        sys.stderr.write("STDOUT is a terminal, terminating!\n")
        return 1

    with pysam.AlignmentFile(args.input, "rb") as infile:
        with pysam.AlignmentFile("-", "wb", template=infile) as outfile:
            return process(args, infile, outfile)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))