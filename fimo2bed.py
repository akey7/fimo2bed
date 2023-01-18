from argparse import ArgumentParser, BooleanOptionalAction
import sys
import csv


# chr1    91369   91487   chr1:91382-91550|carroll_ctcf_mcf7_v45m`2_GTGGCACCAGGTGGCAGCA_cfdna_1   16.2951 +


class Fragment:
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

    def __hash__(self):
        """
        Returns the hash value of the string of the sequence name to use
        when inserting this fragment into a dictionary.

        Returns
        -------
        int
        """

        return hash(self.sequence_name)

    def __str__(self):
        """
        Formats this fragment to output as a row in a bed file.

        Returns
        -------
        str
        """

        return f"{self.chromosome}\t{self.start}\t{self.end}\t{self.sequence_name}|{self.set_name}_{self.serial}\t{self.score}\t{self.strand}"


def fimo_to_bed(file_in, file_out, log_out, set_name, shift=False, center=0):
    unique_fragments = {}

    reader = csv.DictReader(
        filter(lambda row: not row[0].strip().startswith("#"), file_in), delimiter="\t"
    )

    serial = 1
    for row_in in reader:
        frag = Fragment(
            row_in["sequence_name"], row_in["score"], row_in["strand"], set_name, serial
        )
        serial += 1

        if shift:
            start_shift = int(row_in["start"])
            end_shift = int(row_in["stop"])
            frag.shift(start_shift=start_shift, end_shift=end_shift)

        if center != 0:
            frag.center(center)

        unique_fragments[frag] = frag
        # file_out.write(str(frag) + '\n')

    for unique_frag in unique_fragments.values():
        file_out.write(str(unique_frag) + "\n")

    log_out.write(f"total_unique_fragments={len(unique_fragments)}\n")

    file_in.close()
    file_out.close()
    log_out.close()


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Converts a fimo tsv to a bed file. Optionally shifts each fragment to the motif and/or centers the fragment."
    )
    parser.add_argument("--center", default=0, required=False, type=int)
    parser.add_argument("--shift", action=BooleanOptionalAction, default=False)
    parser.add_argument("--set", default="default", required=False, type=str)
    args = parser.parse_args()

    fimo_to_bed(
        sys.stdin,
        sys.stdout,
        sys.stderr,
        args.set,
        shift=args.shift,
        center=args.center,
    )
