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
		(key, val) = field.split(':', 2)
		key = key.strip()
		val = val.strip()
		if val.endswith(';'): val = val[:-1]

		# Add to hash
		if key == 'MOL_ID':
			out += show(vals)	# Show old values
			vals = dict()		# New dictionary

		vals[key] = val
	else:
		if done: break;


out += show(vals)
print out
