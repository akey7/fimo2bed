import sys
import csv
import tracemalloc
from argparse import ArgumentParser, BooleanOptionalAction

from interval import Interval


def serial_numbers(file_in, file_out, set_name):
    pass


if __name__ == '__main__':
    """
    This code block accepts input from command line arguments, opens 
    STDIN, STDOUT, and STDERR and invokes serial_numbers()
    """

    parser = ArgumentParser(
        description='Adds serial numbers to a bed file.'
    )
    parser.add_argument("--set", default="default", required=False, type=str)

    tracemalloc.start()
    sys.stderr.write(f'# Max memory used: {tracemalloc.get_traced_memory()[1]} bytes\n')
    tracemalloc.stop()

    sys.stdin.close()
    sys.stderr.close()
    sys.stdout.close()
