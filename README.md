# simple_paired_read_simulator.py

A very simple paired read simulator using stand-alone python2, i.e. no third-party modules.

"Mutates" template sequence(s) from a FASTA file, introducing substitutions only. These mutations could be a simple model of sequencing error and/or variants in the target sequence.

Sequences are selected with the input with uniform probability (i.e., not weighted by length).

R1 and R2 reads are written in FASTQ format.

The mutation rate per read is specified the integer PctId parameter (positional command-line option; see code).

If PctId=99, then there is a 1% mutation rate, etc.

The mutation rate is identical for every read, unlike most simulators which have a per-base mutation probability
which causes the number of mutations per read to vary stochastically.
