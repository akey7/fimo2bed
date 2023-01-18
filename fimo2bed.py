from argparse import ArgumentParser, BooleanOptionalAction
import sys
import csv


# chr1    91369   91487   chr1:91382-91550|carroll_ctcf_mcf7_v45m`2_GTGGCACCAGGTGGCAGCA_cfdna_1   16.2951 +


class Fragment:
    def __init__(self, fragment_string, score, strand, set_name, serial):
        self.set_name = set_name
        self.serial = serial
        self.score = float(score)
        self.strand = strand
        self.chromosome = fragment_string.split(':')[0]
        locations = fragment_string.split(':')[1].split('-')
        self.start = int(locations[0])
        self.end = int(locations[1])

    def __str__(self):
        return f"{self.chromosome}\t{self.start}\t{self.end}\t{self.chromosome}:{self.start}-{self.end}|{self.set_name}_{self.serial}\t{self.score}\t{self.strand}"

    def shift(self, start_shift=0, end_shift=0):
        self.start += start_shift
        self.end = self.start + end_shift - start_shift - 1

    def center(self, width):
        midpoint = (self.start + self.end) // 2
        self.start = midpoint - width
        self.end = midpoint + width


def fimo_to_bed(file_in, file_out, log_out, set_name, shift=False, center=0):
    unique_fragments = {}
    
    reader = csv.DictReader(filter(lambda row: not row[0].strip().startswith('#'), file_in), delimiter='\t')

    serial = 1
    for row_in in reader:
        frag = Fragment(row_in['sequence_name'], row_in['score'], row_in['strand'], set_name, serial)
        serial += 1

        if shift:
            start_shift = int(row_in['start'])
            end_shift = int(row_in['stop'])
            frag.shift(start_shift=start_shift, end_shift=end_shift)

        if center != 0:
            frag.center(center)

        # file_out.write(str(frag) + '\n')
    log_out.write(f'total_fragments={serial}\n')

    file_in.close()
    file_out.close()
    log_out.close()


if __name__ == '__main__':
    parser = ArgumentParser(description='Converts a fimo tsv to a bed file. Optionally shifts each fragment to the motif and/or centers the fragment.')
    parser.add_argument('--center', default=0, required=False, type=int)
    parser.add_argument('--shift', action=BooleanOptionalAction, default=False)
    parser.add_argument('--set', default="default", required=False, type=str)
    args = parser.parse_args()

    fimo_to_bed(sys.stdin, sys.stdout, sys.stderr, args.set, shift=args.shift, center=args.center)
