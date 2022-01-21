#!/bin/env python3
import sys, argparse
import numpy as np
import pandas as pd
from subprocess import Popen, PIPE
import re

parser=argparse.ArgumentParser()
parser.add_argument('--id', help='Identifier for the job; output files will use this as base')
parser.add_argument('--file', help='Input file (the genome as VCF)')
parser.add_argument('--L', help='Fingerprint length', default=120)
parser.add_argument('--C', help='Too close cutoff', default=20)
args=parser.parse_args()

# Parameters:
L = int(args.L)
C = int(args.C)
id = args.id
file = args.file

# Prepare data structures:
def enumerateKeys(alphabet):
	keys = []
	for ref1 in alphabet:
		for var1 in alphabet:
			if var1==ref1: continue
			for ref2 in alphabet:
				for var2 in alphabet:
					if var2==ref2: continue
					keys.append(ref1+var1+ref2+var2)
	return keys

keys = enumerateKeys(['A', 'C', 'G', 'T'])
count = {key: np.zeros(L) for key in keys}
close = {key: np.zeros(C) for key in keys}

prevchrom = ''
prevpos = 0
prevkey = ''
snvPairs = 0
closeCount = 0

# Process VCF:
### must sanitize file before doing the Popen - reject spaces, semicolons, pipes
vcf = Popen(['gunzip -c ' + file + ' | grep -v ^# | cut -f1,2,4,5'], shell=True, stdout=PIPE, bufsize=1)
###
with vcf.stdout:
	for line in iter(vcf.stdout.readline, b''):
		chrom,pos,ref,alt = line.decode("utf-8").split() ## will fail if there aren't enough fields
		if not re.match("^(chr)?\d", chrom): continue
		if not re.match("^[ACGT]$", ref, flags=re.I): continue
		if not re.match("^[ACGT]$", alt, flags=re.I): continue
		key = (ref+alt).upper()
		pos = int(pos)
		if chrom==prevchrom:
			d = pos-prevpos-1
			if d<0: continue
			snvPairs += 1
			#if snvPairs>=10000: break
			pairkey = prevkey+key
			if d<C:
				closeCount += 1
				close[pairkey][d] += 1
			else:
				count[pairkey][d % L] += 1
		prevchrom = chrom
		prevpos = pos
		prevkey = key

cdf = pd.DataFrame.from_dict(count, orient='index')
closedf = pd.DataFrame.from_dict(close, orient='index')

# Output main fingerprint table:
### old format: should add L as first column in each fingerprint row
with open(id+'.out', 'w') as outf:
	outf.write(''.join(['#source\t', file, '\n']))
	outf.write(''.join(['#SNVpairs\t', str(snvPairs), '\n']))
	outf.write(''.join(['#vectorLengths\t', str(L), '\n']))
	outf.write(''.join(['#tooCloseCutoff\t', str(C), '\n']))
	print(cdf.sort_index(axis=0).to_csv(sep='\t', float_format='%.0f', header=False), file=outf)

# Output secondary (short distance) fingerprint table:
with open(id+'.close', 'w') as outf:
	outf.write(''.join(['#source\t', file, '\n']))
	outf.write(''.join(['#SNVpairs\t', str(snvPairs), '\n']))
	outf.write(''.join(['#vectorLengths\t', str(L), '\n']))
	outf.write(''.join(['#tooCloseCutoff\t', str(C), '\n']))
	print(closedf.sort_index(axis=0).to_csv(sep='\t', float_format='%.0f', header=False), file=outf)

# Normalize fingerprint per position (column), then per keypair (row):
for a in (0,1):
	cdf = cdf.sub(cdf.mean(a), axis=1-a) .div(cdf.std(a), axis=1-a)

# Output normalized fingerprint table:
### old format: should add L as first column in each fingerprint row
with open(id+'.outn', 'w') as outf:
	outf.write(''.join(['#source\t', file, '\n']))
	outf.write(''.join(['#SNVpairs\t', str(snvPairs), '\n']))
	outf.write(''.join(['#vectorLengths\t', str(L), '\n']))
	outf.write(''.join(['#tooCloseCutoff\t', str(C), '\n']))
	print(cdf.sort_index(axis=0).to_csv(sep='\t', float_format='%.3f', header=False), file=outf)


