
revaa <- function(a) {
	return( paste( rev( substring( a, 1:nchar(a), 1:nchar(a) ) ), collapse="") );
}

cat('Reading transitions for AA in contact\n')
tr.aa <- read.table("transitions.aa_in_contact.txt", header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
tr.aa <- as.matrix(tr.aa)

cat('Reading transitions null model \n')
tr.bg <- read.table("transitions.aa_in_contact.txt", header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
tr.bg <- as.matrix(tr.bg)

cn <- colnames(tr.aa)
len <- length(cn)
for( i in 1:(len-1) ) {
	for( j in (i+1):len ) {
		ni <- cn[i]
		nj <- cn[j]
		ri <- revaa( cn[i] )
		rj <- revaa( cn[j] )

		m <- c( tr.aa[ ni, nj ], tr.aa[ni, rj], tr.aa[ri, nj], tr.aa[ri, rj] )
		cat( cn[i] , "\t", cn[j], '\t', ri, '\t', rj, '\t', m, "\n")
	}
}

# heatmap(tr, Rowv=NA, Colv=NA, col = cm.colors(256), scale="none", margins=c(5,10), main="Transitions AA in contact")

