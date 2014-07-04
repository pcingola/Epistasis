
library( 'gplots')

savePlot <- F

#-------------------------------------------------------------------------------
# Reverse an amino acid string
#-------------------------------------------------------------------------------

revaa <- function(a) {
	return( paste( rev( substring( a, 1:nchar(a), 1:nchar(a) ) ), collapse="") );
}

#-------------------------------------------------------------------------------
# Scale by column or row
#-------------------------------------------------------------------------------

scaleCol <- function(x) {
	return( scale( x, center=F, scale=colSums(x) ) )
}

scaleRow <- function(x) {
	tx <- t(x)
	return( t( scale( tx, center=F, scale=colSums(tx) ) ) )
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
# Chech that names match
#-------------------------------------------------------------------------------
checkNames <- function(a, b) {
	if( any( colnames(a) != colnames(b) ) ) {
		stop('Column names do not match!')
	}

	if( any( rownames(a) != rownames(b) ) ) {
		stop('Row names do not match!')
	}
}

#-------------------------------------------------------------------------------
# Compare heatmaps
#-------------------------------------------------------------------------------
heatComp <- function( tr.aa, tr.bg ) {
	taa <- tr.aa
	taa <- scaleRow( taa )
	diag(taa) <- 0
	heatmap.2(taa, main = "Transitions: in contact", sub='Note: Row scaled. Main diagonal set to 0', Rowv=F, Colv=F, col = redgreen(100), density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none"); 

	tbg <- tr.bg
	tbg <- scaleRow( tbg )
	diag(tbg) <- 0
	heatmap.2(tbg, main = "Transitions: 'null'", sub='Note: Row scaled. Main diagonal set to 0', Rowv=F, Colv=F, col = redgreen(100), density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none"); 

	ta <- tr.aa / sum( as.numeric( tr.aa) )
	tg <- tr.bg / sum( as.numeric( tr.bg) )
	t <- ta / tg
	t[ is.na(t) ] <- 0
	t <- scaleRow(t)
	heatmap.2(t, main = "Transitions Ratio: 'in contact' / 'null'", sub='Note: Row scaled. Main diagonal set to 0', Rowv=F, Colv=F, col = redgreen(100), density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none"); 
	return(t)
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

if( savePlot ) { png( width=800, height=800 ) }

#---
# Load data
#---

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

checkNames(tr.aa, tr.bg)
checkNames(tr.aa2, tr.bg2)

if( T ) {
	par( mfcol=c(2,1) )
	xlim <- c(0,20)
	histDens( log2(as.numeric(tr.aa2)), "Log2[ transition count 'in contact']", xlim )
	histDens( log2(as.numeric(tr.bg2)), "Log2[ transition count 'null']", xlim )
}

#---
# Transition reversal ratios
#---
if( F ) {
	min.count <- 20

	cn <- colnames(tr.aa2)
	len <- length(cn)
	h <- 1
	r.aa <- ( 1:(len*(len-1)/2) ) * 0
	r.bg <- r.aa
	for( i in 1:(len-1) ) {
		ni <- cn[i]
		ri <- revaa( cn[i] )
		cat( ni , "\n")

		for( j in (i+1):len ) {
			nj <- cn[j]
			rj <- revaa( cn[j] )

			#m.aa <- c( tr.aa2[ ni, nj ], tr.aa2[ni, rj], tr.aa2[ri, nj], tr.aa2[ri, rj] )
			#m.bg <- c( tr.bg2[ ni, nj ], tr.bg2[ni, rj], tr.bg2[ri, nj], tr.bg2[ri, rj] )

			if((tr.aa2[ ni, nj ] > min.count) && (tr.aa2[ri, rj] > min.count)) {
				r.aa[h] <- tr.aa2[ ni, nj ] / tr.aa2[ri, rj]
			} else {
				r.aa[h] <- NA
			}

			if((tr.bg2[ ni, nj ] > min.count) && (tr.bg2[ri, rj] > min.count)) {
				r.bg[h] <- tr.bg2[ ni, nj ] / tr.bg2[ri, rj]
			} else {
				r.bg[h] <- NA
			}

			h <- h+1
		}
	}

	# Filter
	r.aa <- r.aa[ !is.infinite(r.aa) & !is.nan(r.aa) & !is.na(r.aa) & (r.aa > 0) ]
	r.bg <- r.bg[ !is.infinite(r.bg) & !is.nan(r.bg) & !is.na(r.bg) & (r.bg > 0)  ]

	# Show histograms
	par( mfcol=c(2,1) )
	xlim <- c(0,5)
	histDens(r.aa, 'Transition reversal ratios AA in contact', xlim)
	histDens(r.bg, 'Transition reversal ratios null distribution', xlim)
	xlim <- c(-4,4)
	histDens(log2(r.aa), 'Log2[ Transition reversal ratios AA in contact ] ', xlim)
	histDens(log2(r.bg), 'Log2[ Transition reversal ratios null distribution ]', xlim)
}

#---
# Transition heatmaps
#---
if( F ) {
	t <- heatComp( tr.aa, tr.bg )
	t <- heatComp( tr.aa2, tr.bg2 )
}

#---
# Transitions ratio histograms
#---
if( F ) {
	ta <- as.numeric( tr.aa )
	tg <- as.numeric( tr.bg )
	t <- as.vector( (ta/sum(ta)) / (tg / sum(tg)) )
	histDens( t, 'Normalized counts ratio: AA in contact / null', c(0,2), breaks=10  )

	ta <- as.numeric( tr.aa2 )
	tg <- as.numeric( tr.bg2 )
	t <- as.vector( (ta/sum(ta)) / (tg / sum(tg)) )
	t <- t[ ! is.na(t) ]
	histDens( t, 'Normalized counts ratio: AA-Pairs in contact / null', c(0,2)  )
}


if( F ) {
	pa <- scaleRow( tr.aa )
	pg <- scaleRow( tr.bg )

	pg2 <- matrix( 0 , nrow=400, ncol=400 )
	colnames( pg2 ) <- colnames( tr.bg2 )
	rownames( pg2 ) <- rownames( tr.bg2 )
	pa2 <- pg2

	cn <- colnames(pg)
	len <- length(cn)
	for( i1 in cn ) {
		cat(i1, '\n')
		for( i2 in cn ) {
			i = paste(i1, '_', i2, sep='')
			for( j1 in cn ) {
				for( j2 in cn ) {
					j = paste(j1, '_', j2, sep='')
					pa2[i,j] <- pa[i1,j1] * pa[i2,j2]
					pg2[i,j] <- pg[i1,j1] * pg[i2,j2]
				}
			}
		}
	}

	ta2 <- scaleRow( tr.aa2 )
	ra <- ta2 / pa2
	lra <- log2( ra )

	tg2 <- scaleRow( tr.bg2 )
	rg <- tg2 / pg2
	lrg <- log2( rg )

	par( mfcol=c(2,1) )
	xlim <-c(-5, 5)
	histDens( lra, "Log2[ P(AB -> XY) / ( P(A -> X) * P(B -> Y) ) ] 'in contact'", xlim )
	histDens( lrg, "Log2[ P(AB -> XY) / ( P(A -> X) * P(B -> Y) ) ] 'null'", xlim )
}


if( savePlot )	{ dev.off() } 

