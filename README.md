# fimo2bed
A script to convert fimo `.tsv` files into `.bed` files, with options to sort, drop duplicates, shift to motif, and center the fragments.

It accepts fimo input on STDIN, writes results to STDOUT, and creates log messages on STDERR.

An example command to invoke it is:

```
cat data/fimo.tsv |  python fimo2bed.py --set dux4rep1 --shift --center 50 > /dev/null 2> data/conversion_log.tsv
```
