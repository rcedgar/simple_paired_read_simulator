#!/usr/bin/python2

# A very simple paired read simulator using stand-alone
# python2, i.e. no third-party modules.

# Author Robert C. Edgar
# email robert@drive5.com

# "Mutates" template sequence(s) from a FASTA file,
# introducing substitions only.
# Sequences are selected with the input with uniform
# probability (i.e., not weighted by length).
# R1 and R2 reads are written in FASTQ format.
# The mutation rate per read is specified the integer PctId
# parameter (positional command-line option, see below).
# If PctId=99, then there is a 1% mutation rate, etc.
# The mutation rate is identical for every read, unlike
# most simulators which have a per-base mutation probability
# which causes the number of mutations per read to vary
# stochastically.

import sys
import random

# Positional parameters specified on the command line.
FastaFileName = sys.argv[1]
ReadLength = int(sys.argv[2])
ReadCount = int(sys.argv[3])
PctId = int(sys.argv[4])
R1FileName = sys.argv[5]
R2FileName = sys.argv[6]

fR1 = open(R1FileName, "w")
fR2 = open(R2FileName, "w")

### Start hard coded parameters ###

# The "construct length" is the number of bases between
# first base in R1 and last base in R2, inclusive.
ConstructMean = 350
ConstructStdDev = 50

# Random number seed, to enable reproducible results.
RandSeed = 1
random.seed(RandSeed)

# Quality score for one base in FASTQ.
# All quals will have this value.
PHRED = "5"
assert len(PHRED) == 1

### End hard coded parameters ###

def RevCompLetter(c):
	c = c.upper()
	if c == 'A':
		return 'T'
	elif c == 'C':
		return 'G'
	elif c == 'G':
		return 'C'
	elif c == 'T':
		return 'A'
	return 'N'

def RevComp(s):
	global Map
	t = ""
	n = len(s)
	for i in range(0, n):
		c = s[n-i-1]
		t += RevCompLetter(c)
	return t

def MutateLetter(c):
	r = random.randint(0, 3)
	NewLetter = "ACGT"[r]
	if NewLetter == c:
		NewLetter = "ACGT"[(r+1)%4]
	return NewLetter

def Mutate(Seq, PctId):
	assert PctId >= 50 and PctId <= 100
	assert len(Seq) == ReadLength
	L = len(Seq)
	MutatedBaseCount = (L*(100 - PctId))/100
	MutVec = []
	for i in range(0, MutatedBaseCount):
		MutVec.append(True)
	for i in range(MutatedBaseCount, L):
		MutVec.append(False)
	random.shuffle(MutVec)
	MutatedSeq = ""
	for i in range(0, L):
		c = Seq[i]
		Mut = MutVec[i]
		if Mut:
			MutatedSeq += MutateLetter(c)
		else:
			MutatedSeq += c
	return MutatedSeq

Seqs = {}
Label = ""
CurrSeqVec = []
for Line in open(FastaFileName):
	Line = Line.strip()
	if len(Line) == 0:
		continue
	if Line[0] == ">":
		if len(CurrSeqVec) > 0:
			Seqs[Label] = "".join(CurrSeqVec)
			CurrSeqVec = []
		Label = Line[1:]
	else:
		if Label == "":
			Die("FASTA file does not start with '>'")
		CurrSeqVec.append(Line)
if len(CurrSeqVec) > 0:
	Seqs[Label] = "".join(CurrSeqVec)

Labels = Seqs.keys()
N = len(Labels)

Quals = PHRED * ReadLength
for ReadIndex in range(0, ReadCount):
	Tries = 0
	while True:
		# Avoid infinite loop
		Tries += 1
		assert Tries < 100

		Label = Labels[random.randrange(0, N)]
		Seq = Seqs[Label]
		L = len(Seq)
		ConstructLength = int(random.gauss(ConstructMean, ConstructStdDev))

		if ConstructLength + 100 > L:
			continue

		Lo = random.randrange(10, L - ConstructLength - 10);
		Hi = Lo + ConstructLength -  1

		R1Seq = Seq[Lo:Lo+ReadLength]
		R2Seq = Seq[Hi-ReadLength:Hi]
		assert len(R1Seq) == ReadLength
		assert len(R2Seq) == ReadLength

		Plus = (random.randint(0, 1000000)%2 == 0)
		if Plus:
			R2Seq = RevComp(R2Seq)
		else:
			R1Seq = RevComp(R1Seq)

		if Plus:
			Strand = "+"
		else:
			Strand = "-"
		ReadLabel = "@Sim%d.%d.%s.%d.%d.%c" % (PctId, ReadIndex+1, Label, Lo, Hi, Strand)

		R1Seq = Mutate(R1Seq, PctId)
		R2Seq = Mutate(R2Seq, PctId)

		assert len(R1Seq) == ReadLength
		assert len(R2Seq) == ReadLength

		print >> fR1, ReadLabel + " 1"
		print >> fR2, ReadLabel + " 2"

		print >> fR1, R1Seq
		print >> fR2, R2Seq

		print >> fR1, "+"
		print >> fR2, "+"

		print >> fR1, Quals
		print >> fR2, Quals
		break
