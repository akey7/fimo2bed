import sys
import tracemalloc
from argparse import ArgumentParser


def serial_numbers(file_in, file_out, set):
    serial = 1

    for l in file_in:
        row = l.strip().split('\t')
        interval = f'{row[0]}:{row[1]}-{row[2]}|{set}_{serial}'
        file_out.write(f'{row[0]}\t{row[1]}\t{row[2]}\t{interval}\t0.0\t+\n')
        serial += 1


if __name__ == '__main__':
    """
    This code block accepts input from command line arguments, opens 
    STDIN, STDOUT, and STDERR and invokes serial_numbers()
    """

    parser = ArgumentParser(
        description='Adds serial numbers to a bed file.'
    )
    parser.add_argument("--set", default="default", required=False, type=str)
    args = parser.parse_args()

    tracemalloc.start()
    serial_numbers(sys.stdin, sys.stdout, args.set)
    sys.stderr.write(f'# Max memory used: {tracemalloc.get_traced_memory()[1]} bytes\n')
    tracemalloc.stop()

    sys.stdin.close()
    sys.stderr.close()
    sys.stdout.close()
