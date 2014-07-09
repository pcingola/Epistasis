
library('gplots')
library('expm')

savePlot <- T

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

# Keep numbers: remove NA, Inf and NaN
nums <- function(x) {
	k <- !is.na(x) & !is.nan(x) & !is.infinite(x)
	return( x[k] )
}

# Keep numbers: remove NA, Inf and NaN
numsZero <- function(x) {
	k <- is.na(x) | is.nan(x) | is.infinite(x)
	x[k] <- 0
	return( x )
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

	xlab = paste('mean:', mean(x), ',    median:', median(x), ',    stdev:', sd(x))
    h <- hist(data, main=title, xlim=xlim, xlab = xlab, ylab = "Frequency", freq = T, breaks=breaks);

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
heatCompLog <- function( m1, m2, name1='in contact', name2='null', minCount = -1 ) {
	taa <- scaleRow( m1 )
	diag(taa) <- 0
	heatmap.2(taa, main = paste("Transitions:", name1), sub='Note: Row scaled. Main diagonal set to 0', Rowv=F, Colv=F, col = redgreen(100), density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none"); 

	tbg <- scaleRow( m2 )
	diag(tbg) <- 0
	heatmap.2(tbg, main = paste("Transitions:", name2), sub='Note: Row scaled. Main diagonal set to 0', Rowv=F, Colv=F, col = redgreen(100), density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none"); 

	# Ratio
	ta <- scaleRow( m1 )
	tb <- scaleRow( m2 )
	t <- ta / tb
	t <- log2(t) 
	t[ is.na(t) ] <- 0
	t[ is.nan(t) ] <- 0
	t[ is.infinite(t) ] <- 0
	t[ m1 < minCount ] <- 0
	t[ m2 < minCount ] <- 0

	if( minCount > 0 )	{ sub <- paste( 'Note: Row scaled. Min. count: ', minCount ) }
	else 				{ sub <- 'Note: Row scaled.' }

	heatmap.2(t, main = paste("Transitions: Log2[ ", name1, " / ", name2, " ]"), sub=sub, Rowv=F, Colv=F, col = redgreen(100), density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none", na.rm=T); 

	return(t)
}

#-------------------------------------------------------------------------------
# Compare heatmaps
#-------------------------------------------------------------------------------
heatComp <- function( m1, m2, name1='in contact', name2='null', minCount = -1) {
	taa <- scaleRow( m1 )
	diag(taa) <- 0
	heatmap.2(taa, main = paste("Transitions:", name1), sub='Note: Row scaled. Main diagonal set to 0', Rowv=F, Colv=F, col = redgreen(100), density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none"); 

	tbg <- scaleRow( m2 )
	diag(tbg) <- 0
	heatmap.2(tbg, main = paste("Transitions:", name2), sub='Note: Row scaled. Main diagonal set to 0', Rowv=F, Colv=F, col = redgreen(100), density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none"); 

	# Ratio
	t <- m1 / m2

	if( minCount > 0 )	{ sub <- paste( 'Note: Row scaled. Min. count: ', minCount ) }
	else 				{ sub <- 'Note: Row scaled.' }

	heatmap.2(t, main = paste("Transitions: ", name1, " / ", name2), sub=sub, Rowv=F, Colv=F, col = redgreen(100), density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none", na.rm=T); 

	return(t)
}

#-------------------------------------------------------------------------------
# Compare heatmaps
#-------------------------------------------------------------------------------
heatCompLog <- function( m1, m2, name1='in contact', name2='null', minCount = -1) {
	taa <- scaleRow( m1 )
	diag(taa) <- 0
	heatmap.2(taa, main = paste("Transitions:", name1), sub='Note: Row scaled. Main diagonal set to 0', Rowv=F, Colv=F, col = redgreen(100), density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none"); 

	tbg <- scaleRow( m2 )
	diag(tbg) <- 0
	heatmap.2(tbg, main = paste("Transitions:", name2), sub='Note: Row scaled. Main diagonal set to 0', Rowv=F, Colv=F, col = redgreen(100), density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none"); 

	# Log2[ Ratio ]
	ta <- scaleRow( m1 )
	tb <- scaleRow( m2 )
	t <- ta / tb
	t <- log2(t) 
	t[ is.na(t) ] <- 0
	t[ is.nan(t) ] <- 0
	t[ is.infinite(t) ] <- 0
	t[ m1 < minCount ] <- 0
	t[ m2 < minCount ] <- 0

	if( minCount > 0 )	{ sub <- paste( 'Note: Row scaled. Min. count: ', minCount ) }
	else 				{ sub <- 'Note: Row scaled.' }

	heatmap.2(t, main = paste("Transitions: Log2[ ", name1, " / ", name2, " ]"), sub=sub, Rowv=F, Colv=F, col = redgreen(100), density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none", na.rm=T); 

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
	tr.bg2 <- read.table("transitions.aa_pairs.bg_within_prot.txt", header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
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

if( F ) {
	par( mfcol=c(3,1) )
	xlim <- c(0,15)

	ta <- as.numeric(tr.aa2)
	l <- log2(ta)
	histDens( nums( log2(ta) ), "Log2[ transition count 'in contact' ]", xlim )

	ta <- ta[ ta > 20 ]
	histDens( nums( log2(as.numeric(ta)) ), "Log2[ transition count 'in contact' ], only count > 20", xlim )

	xlim <- c(-4,4)
	histDens( scale(nums( log2(as.numeric(tr.bg2)))), "Scale[ Log2( transition count 'null' ) ]", xlim )
}

#---
# Transition reversal ratios
#---
if( F ) {
	min.count <- 20

	cat("Transition reversal ratios\n")
	cn <- colnames(tr.aa2)
	len <- length(cn)
	h <- 1
	r.aa <- ( 1:(len*(len-1)/2) ) * 0
	r.bg <- r.aa
	for( i in 1:(len-1) ) {
		ni <- cn[i]
		ri <- revaa( cn[i] )
		cat("\t", ni , "\n")

		for( j in (i+1):len ) {
			nj <- cn[j]
			rj <- revaa( cn[j] )

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
	t <- heatCompLog( tr.aa, tr.bg )
	t <- heatCompLog( tr.aa2, tr.bg2 )
	t <- heatCompLog( tr.aa2, tr.bg2, minCount = 20 )
}

#---
# Transitions ratio histograms
#---
if( F ) {
	par( mfcol=c(2,1) )

	# Single AA 
	ta <- as.numeric( tr.aa )
	tg <- as.numeric( tr.bg )
	t <- as.vector( (ta/sum(ta)) / (tg / sum(tg)) )
	t <- nums( t )
	histDens( t, 'Normalized counts ratio: SINGLE AA in contact / null', c(0,2), breaks=15  )
	histDens( t[ ta > 20 & tg > 20 ], 'Normalized counts ratio: SINGLE AA in contact / null (counts > 20)', c(0,2), breaks=15  )

	# AA-Pairs
	ta <- as.numeric( tr.aa2 )
	tg <- as.numeric( tr.bg2 )
	t <- as.vector( (ta/sum(ta)) / (tg / sum(tg)) )
	histDens( t, 'Normalized counts ratio: AA-PAIRS in contact / null', c(0,2)  )
	histDens( t[ ta > 20 & tg > 20 ], 'Normalized counts ratio: AA-PAIRS in contact / null (counts > 20)', c(0,2)  )
}

#---
# Transition of AA-pairs compared to Single-AA probabilities
#---
if( F ) {
	cat("Transition AA-Pairs background probability from single AA transition:\n")
	pa <- scaleRow( tr.aa )
	pg <- scaleRow( tr.bg )

	pg2 <- matrix( 0 , nrow=400, ncol=400 )
	colnames( pg2 ) <- colnames( tr.bg2 )
	rownames( pg2 ) <- rownames( tr.bg2 )
	pa2 <- pg2

	cn <- colnames(pg)
	len <- length(cn)
	for( i1 in cn ) {
		cat('\t', i1, '\n')
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
	t <- heatCompLog( ta2, pa2, 'P[AB => XY] in contact', 'P[A => X] P[B => Y] in contact')
	ra <- ta2 / pa2
	lra <- log2( ra )

	tg2 <- scaleRow( tr.bg2 )
	t <- heatCompLog( tg2, pg2, 'P[AB => XY] null', 'P[A => X] P[B => Y] null')
	rg <- tg2 / pg2
	lrg <- log2( rg )

	par( mfcol=c(2,1) )
	xlim <-c(-5, 5)
	histDens( nums(lra), "Log2[ P(AB -> XY) / ( P(A -> X) * P(B -> Y) ) ] 'in contact'", xlim )
	histDens( nums(lrg), "Log2[ P(AB -> XY) / ( P(A -> X) * P(B -> Y) ) ] 'null'", xlim )
}

#---
# Compare Qhat calculation methods
#---
if( F ) {
	files <- c('Q_HAT_METHOD_0.txt', 'Q_PRIME_HAT_METHOD_0.txt', 'Q_HAT_METHOD_1.txt', 'Q_PRIME_HAT_METHOD_1.txt')
	for( file in files ) {
		Qhat <- read.table(file, header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
		Qhat <- as.matrix(Qhat)
		heatmap.2(Qhat, main = "Qhat", sub=file, Rowv=F, Colv=F, col = redgreen(100), density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none"); 
	}

	Qhat.1 <- read.table('Q_PRIME_HAT_METHOD_0.txt', header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
	Qhat.1 <- as.matrix(Qhat.1)

	Qhat.2 <- read.table('Q_PRIME_HAT_METHOD_1.txt', header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
	Qhat.2 <- as.matrix(Qhat.2)

	heatComp(Qhat.1, Qhat.2, 'Qhat.1', 'Qhat.2')
	heatComp(Qhat.2, Qhat.1, 'Qhat.2', 'Qhat.1')
}

#---
# Compare to Qhat to PAM1
#---
if( T ) {
	# Load Qhat
	Qhat <- read.table('Q_HAT_METHOD_1.txt', header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
	Qhat <- as.matrix(Qhat)

	# Load AA frequencies
	aa.freq <- read.table('aa.frequencies.txt', header = F, row.names = 1, sep="\t", na.strings = 'null')
	aa.freq <- as.vector(aa.freq)
	aa.freq <- aa.freq / sum(aa.freq)

	# Create a function to calculate the probability of amino acid mutation, based on Qhat and time
	Pt <- function(t) { expm(t * Qhat); }
	Mut <- function(t) { sum( diag( Pt(t) ) * aa.freq ); }
	Mut.PAM1 <- function(t)	{ Mut(t) - 0.99; }

	# Solve for PAM (i.e. 't' such that there is 1% probability of mutation
	res <- uniroot( Mut.PAM1, c(0,1) )
	t0 <- res$root
	cat('PAM 1:\tt0 =', t0, '\tMutation(t0): ', Mut(t0), '\n')

	# Load PAM1 matrix
	pam <- read.table('pam1.txt', header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
	pam <- as.matrix(pam)

	# Re-organize PAM1 matrix to have the same order than Qhat
	p <- pam * 0
	colnames(p) <- colnames(Qhat)
	rownames(p) <- rownames(Qhat)
	for( rn  in rownames(pam) ) {
		for( cn  in colnames(pam) ) {
			p[ rn, cn ] <- pam[ rn, cn ]
		}
	}
	pam <- p

	Pt0 <- round( 10000 * Pt(t0) )

	d <- max(dim(pam))
	err <- pam * 0
	for( i in 1:d ) {
		for( j in 1:d ) {
			err[i,j] <- abs(pam[i,j] - Pt0[i,j]) / max(pam[i,j], Pt0[i,j])
		}
	}
	err[ is.nan(err) ] <- 0

	heatmap.2(err, main = "PAM1", sub="", Rowv=F, Colv=F, col = redgreen(100), density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none", na.rm=T); 
}

#---
# Qhat2
#---
if( F ) {
	Qhat2 <- read.table('Q_HAT2.txt', header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
	Qhat2<- as.matrix(Qhat2)
	q2 <- Qhat2
	diag(q2) <- 0
	heatmap.2(q2, main = "Qhat2", sub="Diagonal set to zero", Rowv=F, Colv=F, col = redgreen(100), density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none", na.rm=T); 
}

if( savePlot )	{ dev.off() } 

