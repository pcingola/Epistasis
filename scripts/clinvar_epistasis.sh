#!/bin/sh

# # Annorate clinVar
# cat db/GRCh37/clinVar/clinvar-latest.vcf | snpeff -v hg19 > clinvar-latest.eff.vcf

#-------------------------------------------------------------------------------
# Some stats:
#
# Variants in ClinVar:
#
#	$ grep -v "^#" clinvar-latest.eff.vcf | wc -l 
#	88714
#
# 'Clinically relevant' variants
#		$ cat clinvar-latest.eff.vcf | snpsift filter "(CLNSIG[*] >= 2 & CLNSIG[*] <= 6) " | wc -l
#		25921
#
# 'Clinically relevant' SNPs
#		$ cat clinvar-latest.eff.vcf | snpsift filter "(CLNSIG[*] >= 2 & CLNSIG[*] <= 6) & VC = 'SNV'" | wc -l
#		22392
#
# 'Clinically relevant' SNPs with moderate or high impact
#		$ cat clinvar-latest.eff.vcf | snpsift filter "(CLNSIG[*] >= 2 & CLNSIG[*] <= 6) & VC = 'SNV' & ( EFF[*].IMPACT = 'MODERATE' | EFF[*].IMPACT = 'HIGH')" | wc -l
#		17471
#-------------------------------------------------------------------------------

# Filter:
cat clinvar-latest.eff.vcf | snpsift filter "(CLNSIG[*] >= 2 & CLNSIG[*] <= 6) & VC = 'SNV' & ( EFF[*].IMPACT = 'MODERATE' | EFF[*].IMPACT = 'HIGH')" > clinvar-latest.pass.vcf

# Have 3D structure: 
#		$ grep S3D clinvar-latest.pass.vcf | wc -l
#		7172


