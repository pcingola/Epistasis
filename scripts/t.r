
library( 'gplots')

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
# Compare heatmaps
#-------------------------------------------------------------------------------
heatComp <- function( tr.aa, tr.bg ) {
	taa <- tr.aa
	diag(taa) <- 0
	taa <- scale( taa, center=F, scale=colSums( taa ) )
	heatmap.2(taa, main = "Transitions: in contact", sub='Note: Column scaled. Main diagonal set to 0', Rowv=F, Colv=F, col = redgreen(100), density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none"); 

	tbg <- tr.bg
	diag(tbg) <- 0
	tbg <- scale( tbg, center=F, scale=colSums( tbg ) )
	heatmap.2(tbg, main = "Transitions: 'null'", sub='Note: Column scaled. Main diagonal set to 0', Rowv=F, Colv=F, col = redgreen(100), density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none"); 

	t <- tr.aa / tr.bg
	t[ is.na(t) ] <- 0
	t <- scale( t, center=F, scale=colSums( t ) )
	heatmap.2(t, main = "Transitions Ratio: 'in contact' / 'null'", sub='Note: Column scaled. Main diagonal set to 0', Rowv=F, Colv=F, col = redgreen(100), density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none"); 
	return(t)
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

if( savePlot ) { png( width=800, height=800 ) }

if( ! exists('tr.aa2') ) {
	cat('Reading transitions for AA in contact\n')
	tr.aa2 <- read.table("transitions.aa_pairs.in_contact.txt", header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
	tr.aa2 <- as.matrix(tr.aa2)
}

if( ! exists('tr.bg2') ) {
	cat('Reading transitions null model \n')
	tr.bg2 <- read.table("transitions.aa_pairs.bg_rand.txt", header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
	tr.bg2 <- as.matrix(tr.bg2)
}

if( ! exists('tr.aa') ) {
	cat('Reading transitions for AA in contact\n')
	tr.aa <- read.table("transitions.aa_single.in_contact.txt", header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
	tr.aa <- as.matrix(tr.aa)
}

if( ! exists('tr.bg') ) {
	cat('Reading transitions null model \n')
	tr.bg <- read.table("transitions.aa_single.bg.txt", header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
	tr.bg <- as.matrix(tr.bg)
}

# Transition reversal ratios
if( F ) {
	cn <- colnames(tr.aa2)
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

			m.aa <- c( tr.aa2[ ni, nj ], tr.aa2[ni, rj], tr.aa2[ri, nj], tr.aa2[ri, rj] )
			m.bg <- c( tr.bg2[ ni, nj ], tr.bg2[ni, rj], tr.bg2[ri, nj], tr.bg2[ri, rj] )
			r.aa[h] <-  tr.aa2[ ni, nj ] / tr.aa2[ri, rj]
			r.bg[h] <-  tr.bg2[ ni, nj ] / tr.bg2[ri, rj]
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


t <- heatComp( tr.aa, tr.bg )
t <- heatComp( tr.aa2, tr.bg2 )

if( savePlot )	{ dev.off() } 

