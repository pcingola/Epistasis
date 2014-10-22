
savePlot <- T

if( savePlot ) { pdf(); }

#---
# Analysis of 13K data
#---

if( FALSE ) {
	# Load transposed PCA table
	fileName <- 'pheno.covariates.T2D_13K.transposed.txt';
	d <- read.csv( fileName, sep='\t')

	# Plot all PCA pairs
	#pdf()
	nmin <- 3
	nmax <- 12
	for( n in nmin:(nmax-1) ) {
		plot( d[,n] , d[,n+1], main=paste(names(d)[n] , 'vs', names(d)[n+1]) )
	}
}

#---
# Analysis of 26K data
#---

if( TRUE ) {
	# Load transposed PCA table
	fileName <- 'pheno.covariates.T2D_26K.transposed.txt';
	d <- read.csv( fileName, sep='\t')

	# Plot all PCA pairs
	nmin <- 2
	nmax <- 10
	for( n in nmin:(nmax-1) ) {
		x <- as.numeric( d[,n] )
		y <- as.numeric( d[,n+1] )
		plot( x, y, main=paste(names(d)[n] , 'vs', names(d)[n+1]) )
	}
}

if( savePlot ) { dev.off(); }
