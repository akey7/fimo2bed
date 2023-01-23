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