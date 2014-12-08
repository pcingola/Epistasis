#!/usr/bin/env python

import sys

def show(vals):
	ret = ''
	if 'MOL_ID' not in vals:	return ret

	for key in ['MOL_ID', 'SYNONYM', 'CHAIN']:
		if key in vals:	ret += vals[key]
		ret += '\t'

	return ret

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

done = False

vals = dict()
out = ''

for l in sys.stdin:
	l = l.rstrip()

	if l.startswith('COMPND'):
		done = True
		field = l[10:].strip()

		# Parse compound sub fields
		if ':' in field:
			(key, val) = field.split(':', 1)
			key = key.strip()
		else:
			val = field

		val = val.strip()
		if val.endswith(';'): val = val[:-1]

		# Add to hash
		if key == 'MOL_ID':
			out += show(vals)	# Show old values
			vals = dict()		# New dictionary

		# Add key or append to previous key (line continuation)
		if key in vals: vals[key] += val
		else: vals[key] = val
	else:
		if done: break;


out += show(vals)
print out
