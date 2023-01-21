import sys
import csv
import tracemalloc
from operator import attrgetter
from argparse import ArgumentParser, BooleanOptionalAction


class Interval:
    """
    This class represents a fragment as found in a fimo.tsv file after motif
    calling. After the called motifs are made into fragments, the methods in
    this class will help to elegantly output the data into a bed file.

    This class has magic methods for >, <, >=, <= comparisons. 
    They compare the score of the fragment. Their use is to determine if, given
    a second fragment that points to the same location, the current fragment
    with should be retained of replaced with the new fragment according to
    their scores.
    """

    def __init__(self, fragment_string, score, strand, set_name, serial):
        """
        Parameters
        ----------
        fragment_string: str
            The sequence_name from the time file which has the components to parse
            and fill the attributes of this instance.

        score: str
            Score of the motif match. This score is cast to a float.

        strand: str
            Either "+" or "-"

        set_name: str
            The set name to be appended to each fragment sequence name. The sequence
            name itself is based on the fimo.tsv sequence name and is formatted by
            the sequence name property below.
        
        serial: int
            A serial number to assign to this fragment.
        """
        
        self.set_name = set_name
        self.serial = serial
        self.score = float(score)
        self.strand = strand
        self.chromosome = fragment_string.split(":")[0]
        locations = fragment_string.split(":")[1].split("-")
        self.start = int(locations[0])
        self.end = int(locations[1])

    def shift(self, start_shift=0, end_shift=0):
        """
        Shifts the fragment to the start and end of the motif match.

        Parameters
        ----------
        start_shift: int
            The delta by which to shift the start of the fragment.

        end_shift: int
            The end of the fragment, relative to the start.
        """

        self.start += start_shift
        self.end = self.start + end_shift - start_shift - 1

    def center(self, width):
        """
        Repositions the start and end of the fragment to center over
        the motif. You will probably want to use shift() first.

        Parameters
        ----------
        width: int
            Width on either side of the middle. For example, if you
            want +50 and -50 from the center, this value should be 50.
        """
        midpoint = (self.start + self.end) // 2
        self.start = midpoint - width
        self.end = midpoint + width

    @property
    def sequence_name(self):
        """
        Returns the sequence name to be used in bed file output. This is
        is in the format of the fimo.tsv output. This is used both for bed
        output and hashing the fragments in a dictionary (see __hash__
        below).

        Returns
        -------
        str
        """
        
        return f"{self.chromosome}:{self.start}-{self.end}"

    @property
    def chromosome_sort_key(self):
        """
        Extracts the integer value out of the chromosome string, parses it, and
        returns it. This creates a sort key that can be sorted numericlly
        instead of lexigraphically.

        For example, the following strings will both return 17:
            - "chr17"
            - "chr17_GL000258v2_alt"

        Returns
        -------
        int
            The chromosome number.
        """

        trim_left = self.chromosome[3:]
        trim_right = trim_left.split("_")[0]

        if trim_right == 'X':
            return 100
        elif trim_right == 'Y':
            return 101
        elif trim_right == 'Un':
            return 99
        else:
            return int(trim_right)

    def __hash__(self):
        """
        Returns the hash value of the string of the sequence name to use
        when inserting this fragment into a dictionary and while comparing
        it to other fragments.

        Returns
        -------
        int
        """

        return hash(self.sequence_name)

    def __gt__(self, other):
        """
        Parameters
        ----------
        other: Fragment
            The other fragment to compare the score to.

        Returns
        -------
        bool
            True if the other.score > self.score
        """
        
        return self.score > other.score

    def __lt__(self, other):
        """
        Parameters
        ----------
        other: Fragment
            The other fragment to compare the score to.

        Returns
        -------
        bool
            True if the other.score < self.score
        """
        
        return self.score < other.score

    def __ge__(self, other):
        """
        Parameters
        ----------
        other: Fragment
            The other fragment to compare the score to.

        Returns
        -------
        bool
            True if the other.score >= self.score
        """
        
        return self.score >= other.score

    def __le__(self, other):
        """
        Parameters
        ----------
        other: Fragment
            The other fragment to compare the score to.

        Returns
        -------
        bool
            True if the other.score <= self.score
        """
        
        return self.score <= other.score

    def __eq__(self, other):
        """
        Parameters
        ----------
        other: Fragment
            The other fragment to compare the score to.

        Returns
        -------
        bool
            True if the self.score == other.score
        """

        return self.score != other.score

    def __ne__(self, other):
        """
        Parameters
        ----------
        other: Fragment
            The other fragment to compare the score to.

        Returns
        -------
        bool
            True if the self.score != other.score
        """

        return self.score != other.score

    def __str__(self):
        """
        Formats this fragment to output as a row in a bed file.

        Returns
        -------
        str
        """

        return f"{self.chromosome}\t{self.start}\t{self.end}\t{self.sequence_name}|{self.set_name}_{self.serial}\t{self.score}\t{self.strand}"


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
    
    sys.stdin.close()
    sys.stdout.close()
    sys.stderr.close()
    
    tracemalloc.stop()
