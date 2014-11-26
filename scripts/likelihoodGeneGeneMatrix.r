
savePng <- TRUE

if( savePng ) png(width = 1024, height = 1024)

fileNull <- 'likelihoodGeneGeneMatrix.null.txt'
fileAlt <- 'likelihoodGeneGeneMatrix.alt.txt'

ll.alt <- read.table(fileAlt, sep="\t", header=TRUE)
ll.null <- read.table(fileNull, sep="\t", header=TRUE)

par(mfrow=c(1,1))

for( neigh in 0:9 ) {
	keep <- ll.alt$neigh == neigh
	lla <- as.vector(ll.alt$ll[keep])

	keep <- ll.null$neigh == neigh
	lln <- as.vector(ll.null$ll[keep])

	wt <- wilcox.test(lla,lln,alternative='g')
	#cat('Neighbourhood:', neigh, '\tp-value:', wt$p.value, '\n')
	cat(neigh, '\t', wt$p.value, '\n')

	title <- paste("Best LL(MSA) using", neigh ,"amino acids")
	plot( density( lla ), col='red', main=title, xlab='' )
	lines( density( lln ), col='green' )
}


if( savePng )	dev.off()
