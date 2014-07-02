
savePlot <- T

#-------------------------------------------------------------------------------
# Reverse an amino acid string
#-------------------------------------------------------------------------------

revaa <- function(a) {
	return( paste( rev( substring( a, 1:nchar(a), 1:nchar(a) ) ), collapse="") );
}

#-------------------------------------------------------------------------------
# Histogram & density
#-------------------------------------------------------------------------------

histDens <- function( x, title, xlim, q=1.0, breaks = 50 ) {
    # Show only this part of the data
    xmin <- min( xlim )
    xmax <- max( xlim )
    data <- x[ (x >= xmin) & (x <= xmax) ];

    dens <- density(data)

    h <- hist(data, main=title, xlim=xlim, xlab = "data", ylab = "Frequency", freq = T, breaks=breaks);

    # Adjust density height to 'frecuency'
    dens$y <- max(h$counts) * dens$y/max(dens$y)
    lines(dens, col='red', xlim=xlim)

    # Mean & median calculated over the whola data
    abline( v=mean(x), col='blue', lty=2, lwd=2);
    abline( v=median(x), col='green', lty=2, lwd=2);

    legend("topright",c("Mean","Median"),lty=c(1,1),col=c("blue","green"))
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

if( savePlot ) { png( width=800, height=800 ) }

if( ! exists('tr.aa') ) {
	cat('Reading transitions for AA in contact\n')
	tr.aa <- read.table("transitions.aa_in_contact.txt", header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
	tr.aa <- as.matrix(tr.aa)
}

if( ! exists('tr.bg') ) {
	cat('Reading transitions null model \n')
	tr.bg <- read.table("transitions.bg_rand.txt", header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
	tr.bg <- as.matrix(tr.bg)
}

# Transition reversal ratios
if( F ) {
	cn <- colnames(tr.aa)
	len <- length(cn)
	h <- 1
	r.aa <- ( 1:(len*(len-1)/2) ) * 0
	r.bg <- r.aa
	for( i in 1:(len-1) ) {
		for( j in (i+1):len ) {
			ni <- cn[i]
			nj <- cn[j]
			ri <- revaa( cn[i] )
			rj <- revaa( cn[j] )

			m.aa <- c( tr.aa[ ni, nj ], tr.aa[ni, rj], tr.aa[ri, nj], tr.aa[ri, rj] )
			m.bg <- c( tr.bg[ ni, nj ], tr.bg[ni, rj], tr.bg[ri, nj], tr.bg[ri, rj] )
			r.aa[h] <-  tr.aa[ ni, nj ] / tr.aa[ri, rj]
			r.bg[h] <-  tr.bg[ ni, nj ] / tr.bg[ri, rj]
			# cat( ni , "\t", nj, '\t', ri, '\t', rj, '\t', m.aa, '\t', m.bg, "\t", r.aa[h], '\t', r.bg[h], "\n")
			cat( ni , "\t", nj, "\n")

			h <- h+1
		}
	}

	r.aa <- r.aa[ !is.infinite(r.aa) & !is.nan(r.aa) ]
	r.bg <- r.bg[ !is.infinite(r.bg) & !is.nan(r.bg) ]

	par( mfcol=c(2,1) )

	xlim <- c(0,5)
	histDens(r.aa[ r.aa > 0 ], 'Transition reversal ratios AA in contact', xlim )
	histDens(r.bg, 'Transition reversal ratios null distribution', xlim )
}


taa <- scale( tr.aa, center=F, scale=colSums( tr.aa ) )
diag(taa) <- 0
heatmap(taa, Rowv=NA, Colv=NA, col = cm.colors(256), scale="row", margins=c(5,10), main="Transitions AA in contact")

tbg <- scale( tr.bg, center=F, scale=colSums( tr.bg ) )
diag(tbg) <- 0
heatmap(tbg, Rowv=NA, Colv=NA, col = cm.colors(256), scale="row", margins=c(5,10), main="Transitions 'null'")

t <- taa / tbg
diag(t) <- 0
heatmap(t, Rowv=NA, Colv=NA, col = cm.colors(256), scale="row", margins=c(5,10), main="Transitions ratio")

if( savePlot )	{ dev.off() } 
