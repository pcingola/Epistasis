
revaa <- function(a) {
	return( paste( rev( substring( a, 1:nchar(a), 1:nchar(a) ) ), collapse="") );
}

tr <- read.table("t.txt", header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
tr <- as.matrix(tr)

trn <- read.table("tnull.txt", header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
trn <- as.matrix(tr)

cn <- colnames(tr)
len <- length(cn)
for( i in 1:(len-1) ) {
	for( j in (i+1):len ) {
		ni <- cn[i]
		nj <- cn[j]
		ri <- revaa( cn[i] )
		rj <- revaa( cn[j] )

		m <- c( tr[ ni, nj ], tr[ni, rj], tr[ri, nj], tr[ri, rj] )
		cat( cn[i] , "\t", cn[j], '\t', ri, '\t', rj, '\t', m, "\n")
	}
}

# heatmap(tr, Rowv=NA, Colv=NA, col = cm.colors(256), scale="none", margins=c(5,10), main="Transitions AA in contact")

