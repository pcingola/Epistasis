#!/usr/bin/env python

import sys
import random

for line in sys.stdin:
	if line.startswith("#"): continue
	line = line.rstrip()

	out = ""
	idx = 0
	for f in line.split('\t'):
		if idx > 0: out += '\t'

		if idx > 8:
			# Randomly select genotype
			r = random.random()

			if r < 0.25 : out += "1/1"
			elif r < 0.5 : out += "0/1"
			else : out += "0/0"

		else:
			out += f

		idx += 1

	print out
	

