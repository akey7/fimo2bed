import sys
import csv
import tracemalloc
from operator import attrgetter
from argparse import ArgumentParser, BooleanOptionalAction

from interval import Interval


def fimo_to_bed(file_in, file_out, log_out, sort, set_name, shift=False, center=0):
    """
    When reading the source fimo.tsv, it filters out rows with a leading "#"
    character.

    Parameters
    ----------
    file_in
        A file handle to the input stream for example, sys.stdin.

    file_out
        A file handle to the output stream for example, sys.stdout.

    log_out
        A file handle to the loggin stream for example, sys.stderr

    sort
        True if the output should be sorted. If sorted, the output will have
        its serial numbers reassigned to match the sort order.

    set_name
        The "set name" which is appended to each sequence name on the
        output.

    shift
        If True, shifts the fragment to the motif

    center
        If non-zero, shifts the fragment to its center +/- the width
        specified by this parameter.

    Returns
    -------
    None
    """
    log_out.write("action\tinterval\treason\n")

    unique_intervals = {}

    reader = csv.DictReader(
        filter(lambda row: not row[0].strip().startswith("#"), file_in), delimiter="\t"
    )

    serial = 1
    for row_in in reader:
        intrvl = Interval(
            row_in["sequence_name"], row_in["score"], row_in["strand"], set_name, serial
        )
        serial += 1

        if shift:
            start_shift = int(row_in["start"])
            end_shift = int(row_in["stop"])
            intrvl.shift(start_shift=start_shift, end_shift=end_shift)

        if center != 0:
            intrvl.center(center)

        if intrvl in unique_intervals and unique_intervals[intrvl] > intrvl:
            log_out.write(f"skip\t{intrvl.sequence_name}\tscore {intrvl.score} less than existing {unique_intervals[intrvl].score}\n")
        elif intrvl in unique_intervals and unique_intervals[intrvl] < intrvl:
            log_out.write(f"replace\t{intrvl.sequence_name}\tscore {intrvl.score} greater than existing {unique_intervals[intrvl].score}\n")
            unique_intervals[intrvl] = intrvl
        else:
            log_out.write(f"append\t{intrvl.sequence_name}\tnew fragment\n")
            unique_intervals[intrvl] = intrvl

    if sort:
        final_intervals = sorted(unique_intervals.values(), key=attrgetter('chromosome_sort_key', 'start', 'end'))
        sort_serial = 1
        for sort_interval in final_intervals:
            sort_interval.serial = sort_serial
            sort_serial += 1
    else:
        final_intervals = unique_intervals.values()

    for unique_intervals in final_intervals:
        file_out.write(str(unique_intervals) + "\n")


if __name__ == "__main__":
    """
    This code block accepts input from command line arguments, opens 
    STDIN, STDOUT, and STDERR and invokes fimo_to_bed()
    """

    parser = ArgumentParser(
        description="Converts a fimo tsv to a bed file. Optionally shifts each fragment to the motif and/or centers the fragment."
    )
    parser.add_argument("--center", default=0, required=False, type=int)
    parser.add_argument("--shift", action=BooleanOptionalAction, default=False)
    parser.add_argument("--set", default="default", required=False, type=str)
    parser.add_argument("--sort", action=BooleanOptionalAction, default=False)
    args = parser.parse_args()

    tracemalloc.start()
    fimo_to_bed(
        file_in=sys.stdin,
        file_out=sys.stdout,
        log_out=sys.stderr,
        sort=args.sort,
        set_name=args.set,
        shift=args.shift,
        center=args.center,
    )
    sys.stderr.write(f'# Max memory used: {tracemalloc.get_traced_memory()[1]} bytes\n')
    tracemalloc.stop()

    sys.stdin.close()
    sys.stdout.close()
    sys.stderr.close()
