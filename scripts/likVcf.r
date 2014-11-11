#-------------------------------------------------------------------------------
#
# P-values correction for gene set analysis (SnpEff gsa)
#
#
#															Pablo Cingolani 2014
#-------------------------------------------------------------------------------

show <- TRUE		# Used for debugging
savePdf <- TRUE

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
		smoothScatter( lc, p, main="Scores", xlab="Number of scores", ylab="Score", sub=title )
		lines(lowess(lc, p ), col='orange')

		smoothScatter( lc, padj, main="Adjusted Scores", xlab="Number of scores", ylab="Score", sub=title )
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
	fileName <- "likelihood.clinvar/top.sorted.txt"
	d.clin = read.table(fileName, sep="\t", header=FALSE)
	names(d.clin) <- c('chr', 'pos', 'ref', 'alt', 'll', 'len')
	
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

lladj.g1k <- normLogLik( d.g1k$ll, d.g1k$len, '1000 Genomes')
lladj.clin <- normLogLik( d.clin$ll, d.clin$len, 'ClinVar')
lladj.hgmd <- normLogLik( d.hgmd$ll, d.hgmd$len, 'HGMD')

plot( density( lladj.hgmd), col='green', main='Distribution of adjusted log-likelihhods', sub='Black: 1KGenomes, Red: Clinvar, Green: HGMD', xlim = c(-30,50) )
lines( density( lladj.clin), col='red' )
lines( density( lladj.g1k) )

if( savePdf )	dev.off(); 
