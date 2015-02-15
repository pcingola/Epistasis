#-------------------------------------------------------------------------------
#
# P-values correction for gene set analysis (SnpEff gsa)
#
#
#															Pablo Cingolani 2014
#-------------------------------------------------------------------------------

show <- TRUE		# Used for debugging
savePdf <- FALSE

pvalCoefThreshold = 10^-12	# Adjust only if linearmodel p-value is less than this


#-------------------------------------------------------------------------------
# Perform a linear regression frmo logLikelihood interaction values
# Use log10( number_of_AA ) as covariate (number of calculated AA interactions)
#-------------------------------------------------------------------------------
normLogLik <- function(p, c, title) {
	# Scores (logLikelihood)
	# Note: the conversion 'as.numeric(as.character())' is necesary 
	# because R reads '0.6049666666666666' as a number, 
	# but reads       '0.60496666666666666' as a string (yes, R is amazing!)
	p = as.numeric(as.character(p))
	lc = log10( c )	# Use log10(count), gives better results

	#---
	# Linear model
	#---
	lmfit = lm( p ~ lc )
	sumLmfit =  summary(lmfit)
	pvalCoef = sumLmfit$coefficients[2,4]

	padj = lmfit$residuals		# Residuals (adjusted scores)

	if( show ) {
		print(sumLmfit)
		cat(title, '\tSlope:\t', sumLmfit$coefficients[2], '\tp-value:\t', pvalCoef, 'File:\t', fileName, '\n');

		par(mfrow=c(1,1))
		smoothScatter( lc, p, main="Scores", xlab="Number of scores", ylab="Score", sub=title, ylim = c(-20,50) )
		lines(lowess(lc, p ), col='orange')

		smoothScatter( lc, padj, main="Adjusted Scores", xlab="Number of scores", ylab="Score", sub=title, ylim = c(-20,50) )
		lines(lowess(lc, padj ), col='orange')

		plot( density(p) , main="Scores distribution", xlab="red:Adjusted black:Unadjusted", sub=title)
		lines( density(padj) , col='red')

		plot( density(lc) , main="log10(Number of scores)", xlab="", sub=title)
	}

	#---
	# Create output file
	# Decide whether to use corrected scores or not.
	# Not very significant? Use original scores
	#---
	#if( pvalCoef > pvalCoefThreshold ) { padj = d$score	}
	#dout = data.frame( geneId = d$geneId, score = padj )
	#write.table( dout, file=outFileName, quote=FALSE, row.names = FALSE, col.names = FALSE, sep="\t")

	return (padj)
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

#---
# Read data
#---
if( ! exists('d.clin') ) {
	fileName <- "likelihood.clinvar/top.sorted.clnsig.txt"
	d.clin = read.table(fileName, sep="\t", header=TRUE)
	#names(d.clin) <- c('chr', 'pos', 'ref', 'alt', 'll', 'len', 'clnsig')
	
	fileName <- "likelihood.g1k/top.sorted.txt"
	d.g1k = read.table(fileName, sep="\t", header=FALSE)
	names(d.g1k) <- c('chr', 'pos', 'ref', 'alt', 'll', 'len')

	fileName <- "likelihood.hgmd/top.sorted.txt"
	d.hgmd = read.table(fileName, sep="\t", header=FALSE)
	names(d.hgmd) <- c('chr', 'pos', 'ref', 'alt', 'll', 'len')
}

#---
# Normalize
#---
if( savePdf )	pdf();

# lladj.g1k <- normLogLik( d.g1k$ll, d.g1k$len, '1000 Genomes')
# lladj.clin <- normLogLik( d.clin$ll, d.clin$len, 'ClinVar')
# lladj.hgmd <- normLogLik( d.hgmd$ll, d.hgmd$len, 'HGMD')
# 
# plot( density( lladj.hgmd), col='green', main='Distribution of adjusted log-likelihhods', sub='Black: 1KGenomes, Red: Clinvar, Green: HGMD', xlim = c(-30,50) )
# lines( density( lladj.clin), col='red' )
# lines( density( lladj.g1k) )

#---
# Plot by clinical significance
#---

# Group by clinical significance (CLNSIG from CinVar)
#	Variant Clinical Significance: 
#		0 - Uncertain significance, 
#		1 - not provided, 
#		2 - Benign, 
#		3 - Likely benign, 
#		4 - Likely pathogenic, 
#		5 - Pathogenic, 
#		6 - drug response, 
#		7 - histocompatibility, 
#		255 - other
# So we create 3 groups:
#		Unknown: 0, 1, 255
#		Benign : 2, 3, 6 (something that has a drug response, is good)
#		Pathogenic: 4, 5
#
# Some stats:
#	CLNSIG: 0 		Count: 21 		LL_mean: 14.7234 	LL_median 9.826063 	Len_mean:  4352.1 	Len_meadian: 935 
#	CLNSIG: 1 		Count: 52 		LL_mean: 23.31171 	LL_median 13.61291 	Len_mean:  1435.9 	Len_meadian: 1481 
#	CLNSIG: 2 		Count: 272 		LL_mean: 34.10323 	LL_median 26.27931 	Len_mean:  3925.3 	Len_meadian: 1172.5 
#	CLNSIG: 3 		Count: 258 		LL_mean: 31.56365 	LL_median 25.56079 	Len_mean: 11619.7 	Len_meadian: 4575 
#	CLNSIG: 4 		Count: 562 		LL_mean: 17.53154 	LL_median 12.04729 	Len_mean:  2132.5 	Len_meadian: 838 
#	CLNSIG: 5 		Count: 4206 	LL_mean: 16.90444 	LL_median 11.73378 	Len_mean:  1057.6 	Len_meadian: 688.5 
#	CLNSIG: 6 		Count: 18 		LL_mean: 32.63433 	LL_median 22.10314 	Len_mean:   729.9 	Len_meadian: 492 
#	CLNSIG: 255 	Count: 10 		LL_mean: 20.09807 	LL_median 11.39337 	Len_mean:   4088 	Len_meadian: 767 
#
# 1000 Genomes	:		LL_mean: 23.2	LL_median: 14.8
# HGMD 			:		LL_mean: 19.8	LL_median: 12.1
# ClinVar		:		LL_mean: 18.6	LL_median: 11.9


plot(density( d.clin$ll ), main="Log-Likelihood by Clinical Significance (CLNSIG)", xlab="Log-likelihood", sub="Black: All, Blue: Unknown, Red: Pathogenic, Green: Benign", xlim = c(0, 75) )	# All entries

# CLNSIG: Unknown
keep <- (d.clin$clnsig == 0) | (d.clin$clnsig == 1) | (d.clin$clnsig == 255)
lines(density( d.clin$ll[keep] ), col = 'blue')

keep <- (d.clin$clnsig == 2) | (d.clin$clnsig == 3) | (d.clin$clnsig == 6)
lines(density( d.clin$ll[keep] ), col = 'green')

keep <- (d.clin$clnsig == 4) | (d.clin$clnsig == 4) | (d.clin$clnsig == 6)
lines(density( d.clin$ll[keep] ), col = 'red')

if( savePdf )	dev.off(); 

