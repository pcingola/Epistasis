#!/usr/bin/env python

# Read Qhat matrix form file
def readAaMatrix(fileName):
	qmatrix = dict()
	with open(fileName) as qfile:
		header = True

		for line in qfile:
	
			vals = line.rstrip().split('\t')
			if header: 
				aas = vals
				header = False
			else:
				aa = vals[0]

				for idx in range(0, len(vals)-1):
					key = aa + "_" + aas[idx]
					#print "\t", key, '\t', vals[idx+1]
					qmatrix[key] = float( vals[idx+1] ) 

	return qmatrix
		
#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

q = readAaMatrix('Qhat.txt')
q2 = readAaMatrix('Qhat2.txt')

for aas1 in q.keys():
	for aas2 in q.keys():
		if aas1 != aas2 and q[aas1] > 0 and q[aas2] > 0:
			ratio1 = q2[ aas1 + '_' + aas2 ] / ( q[aas1] * q[aas2] )
			ratio2 = q2[ aas2 + '_' + aas1 ] / ( q[aas1] * q[aas2] )
			ratio = min(ratio1, ratio2)
			print ratio, aas1, aas2, ratio1, ratio2, q[aas1], q[aas2], q2[ aas1 + '_' + aas2 ], q2[ aas2 + '_' + aas1 ]
