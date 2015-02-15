#!/bin/sh


# Extract sequence from PDB file
grep "^SEQRES" | tr -s " " | cut -f 5- -d " " \
	| sed "s/ALA/A/g" \
	| sed "s/LEU/L/g" \
	| sed "s/ARG/R/g" \
	| sed "s/LYS/K/g" \
	| sed "s/ASN/N/g" \
	| sed "s/MET/M/g" \
	| sed "s/ASP/D/g" \
	| sed "s/PHE/F/g" \
	| sed "s/CYS/C/g" \
	| sed "s/PRO/P/g" \
	| sed "s/GLN/Q/g" \
	| sed "s/SER/S/g" \
	| sed "s/GLU/E/g" \
	| sed "s/THR/T/g" \
	| sed "s/GLY/G/g" \
	| sed "s/TRP/W/g" \
	| sed "s/HIS/H/g" \
	| sed "s/TYR/Y/g" \
	| sed "s/ILE/I/g" \
	| sed "s/VAL/V/g" \
	| tr -d " "

