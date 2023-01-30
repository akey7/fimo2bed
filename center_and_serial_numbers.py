import sys
import tracemalloc
from argparse import ArgumentParser


def serial_numbers(file_in, file_out, set, center=50):
    serial = 1

    for l in file_in:
        row = l.strip().split('\t')
        chrom = row[0]
        original_start = int(row[1])
        original_end = int(row[2])
        midpoint = (original_start + original_end) // 2
        start = midpoint - center
        end = midpoint + center
        interval = f'{chrom}:{start}-{end}|{set}_{serial}'
        file_out.write(f'{chrom}\t{start}\t{end}\t{interval}\t0.0\t+\n')
        serial += 1


if __name__ == '__main__':
    """
    This code block accepts input from command line arguments, opens 
    STDIN, STDOUT, and STDERR and invokes serial_numbers()
    """

    parser = ArgumentParser(
        description='Adds serial numbers and centers motifs to a bed file.'
    )
    parser.add_argument('--set', default='default', required=False, type=str)
    parser.add_argument('--center', default=50, required=False, type=int)
    args = parser.parse_args()

    tracemalloc.start()
    serial_numbers(sys.stdin, sys.stdout, args.set, args.center)
    sys.stderr.write(f'# Max memory used: {tracemalloc.get_traced_memory()[1]} bytes\n')
    tracemalloc.stop()

    sys.stdin.close()
    sys.stderr.close()
    sys.stdout.close()
