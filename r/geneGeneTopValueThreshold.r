
savePng <- TRUE

if( savePng ) png(width = 2048, height = 1024)

files <- c('geneGeneTopValue_threshold.1E-10.tsv'	,'geneGeneTopValue_threshold.1E-20.tsv'	,'geneGeneTopValue_threshold.1E-40.tsv'	,'geneGeneTopValue_threshold.1E-60.tsv', 'geneGeneTopValue_threshold.1E-100.tsv'	,'geneGeneTopValue_threshold.1E-30.tsv'	,'geneGeneTopValue_threshold.1E-50.tsv'	,'geneGeneTopValue_threshold.1E-80.tsv')

for( f in files ) {
	fileAlt <- paste('alt/', f, sep='')
	fileNull <- paste('null/', f, sep='')

	ll.alt <- read.table(fileAlt, sep="\t", header=TRUE)
	ll.null <- read.table(fileNull, sep="\t", header=TRUE)

	par(mfrow=c(1,2))

	sub <- paste('LL threashold (Null):', ll.alt$threshold_ll_null[1])
	plot( density( ll.alt$ll_max ), col='red', main="Histogram LL(MSA)", sub=sub )
	lines( density( ll.null$ll_max ), col='green' )

	x <- log(ll.alt$countPass)
	y <- ll.alt$ll_max
	smoothScatter( x, y , col='red', main="LL(MSA) vs log10(count_scores)", sub=sub )
	lines( lowess(x,y), col='red' )

	x <- log(ll.null$countPass)
	y <- ll.null$ll_max
	points( x, y , col='green', cex=0.1 )
	lines( lowess(x,y), col='green' )
}


if( savePng )	dev.off()
