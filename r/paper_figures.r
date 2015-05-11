
savePng <- T
fig2    <- F
fig3    <- F
fig4    <- F
fig5    <- F
figS1   <- F
figS2   <- F
figS3   <- T

#-------------------------------------------------------------------------------
# Compare two distributions
#-------------------------------------------------------------------------------
compareHist <- function(data1, data2, main, xlim, xlab, ylab, breaks) {
	dens1 <- density(data1)
	h1 <- hist(data1, main=main, xlim=xlim, xlab=xlab, ylab=ylab, freq = F, breaks=breaks, col=rgb(1,0,0,1/4) , add=F );
	lines(dens1, col='red', xlim=xlim)

	dens2 <- density(data2)
	h2 <- hist(data2, xlim=xlim, freq = F, breaks=breaks, col=rgb(0,0,1,1/4), add=T );
	lines(dens2, col='blue', xlim=xlim)
}

#-------------------------------------------------------------------------------
# Cummulative probability ratio
#-------------------------------------------------------------------------------
cummProbRatio <- function(x, ll.int, ll.non) {
	p.int <- sum( ll.int >= x ) / length( ll.int ) 
	p.non <- sum( ll.non >= x ) / length( ll.non ) 
	return(p.int / p.non)
}

#-------------------------------------------------------------------------------
# Figure 2: Compare LL(MSA) of AA 'in contact' vs 'not in contact' using
#           'within protein' data from PDB
#
# Note: Actually the 'not in contact' data is taken as random AA (so there 
#       might be some 'in contact' within the set)
#
# Data calculated using Epistasis.jar (run.bds)
# Result files:
#		./data/likelihood.*.values.txt.gz
#-------------------------------------------------------------------------------
figure2 <- function(lc, la, lalcOdds, x) {

	par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis

	xlim <- c(-40,20)
	title <- "" # "Lig-likelihood (MSA) of amino acids 'in contact' vs 'not in contact'"
	xlab <- 'Log-likelihood'
	breaks <- 200

	cat('Distributions summary:\n')
	cat('\tInteracting\t', summary(lc), '\n')
	cat('\tNon-interacting\t', summary(la), '\n')

	# Histograms
	compareHist(lc, la, title, xlim, xlab, "Frequency", breaks)
	legend("topleft", inset=.05, c('In contact', 'Not in contact', 'Probability ratio'), fill=c( rgb(1,0,0,1/4), rgb(0,0,1,1/4), rgb(0,0,0,1)), horiz=F)

	# Ratio of cummulative probabilities
	keep <- is.finite(lalcOdds)
	lalcOdds[!keep] <- NA
	lor <- log(lalcOdds) 
	par(new = TRUE)
	plot(x, lor, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "")
	axis(side=4, at = pretty(range(lor)))
	mtext("Log odds ratio", side=4, line=3)
	lines( supsmu(x, lor), col='gray', lty=2 )
}

#-------------------------------------------------------------------------------
# Figure 3: Compare LL(MSA) of AA 'in contact' vs 'not in contact' using
#           co-crystalized data from PDB
#
# Data calculated using Epistasis.jar (run.bds)
# Result files:
#		./data/interactions/pdb/likelihood.pdb_compound.neigh_1.*.0.txt.gz
#-------------------------------------------------------------------------------
figure3 <- function(ll.int, ll.non, or, orx) {
	par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis

	#---
	# Show density distributions
	#---
	xlim <- c(-40,20)
	title <- '' # "LL(MSA) of amino acids 'interacting' vs 'non-interacint'"
	xlab <- 'Log-likelihood'
	breaks <- 100

	cat('Distributions summary:\n')
	cat('\tInteracting\t', summary(ll.int), '\n')
	cat('\tNon-interacting\t', summary(ll.non), '\n')

	compareHist(ll.int, ll.non, title, xlim, xlab, "Frequency", breaks)
	legend("topleft", inset=.05, c('Interacting', 'Not in contact', 'Log[Probability ratio]'), fill=c( rgb(1,0,0,1/4), rgb(0,0,1,1/4), rgb(0,0,0,1)), horiz=F)

	#---
	# Figure 3.B: Plot odds ratio (cummulative probability)
	#---
	keep <- is.finite(or)
	or[!keep] <- NA
	lor <- log(or) 

	# Remove inacurrate numbers (too few points to calculate stats)
	maxShow <- 10.7
	lor[ orX > maxShow ] <- NA

	# Ratio of cummulative probabilities
	par(new = TRUE)
	plot( orX, lor, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "")
	axis(side=4, at = pretty(range(lor[keep])))
	mtext("Log odds ratio", side=4, line=3)
	lines( supsmu(orx[keep], lor[keep]), col='gray', lty=2 )
}

#-------------------------------------------------------------------------------
# Figure 4: Histogram of ratios R(ab,cd) in comparisson between Q and Q2
#
# Result files:
#		./data/q_q2_compare.txt
#-------------------------------------------------------------------------------
figure4 <- function() {
    cat('Loading Q_Q2 compare data\n')
    qcmp <- read.table('q_q2_compare.ratio.txt', sep='\t', header=F)

	# Ratio
	r <- as.vector(unlist(qcmp))

	# Log ratio
	l <- log(r)
    cat('Raw:\t\tMean:', mean(r), ' StdDev:', sd(r), ' Median:', median(r), '\n')
    cat('Log:\t\tMean:', mean(l), ' StdDev:', sd(l), ' Median:', median(l), '\n')

	title <- 'Histogram of Log[R(ab,cd)] ratios'
    plot( density(l), main=title, xlab='Log[ R(ab,cd) ], Green: median, Blue: mean')
	abline( v=median(l), col='green', lty=2, lwd=2);
	abline( v=mean(l), col='blue', lty=2, lwd=2);
}

#-------------------------------------------------------------------------------
# Figure 5: Power analysis figures
#-------------------------------------------------------------------------------
figure5 <- function() {

	# Big plot size
	#plotSize <- 1 * 1024
	#if( savePlot ) { png( width=plotSize, height=plotSize ) }

	if( ! exists('pow') ) {
		#pow <- read.table('logisticRegressionPowerGT_parseResults_CoEvolution.txt', sep="\t", header=TRUE)
		pow <- read.table('logisticRegressionPowerGT_parseResults.txt', sep="\t", header=TRUE)
	}

	afs <- c(0.01, 0.05, 0.10)

	af1s <- unique(pow$af1)
	af1s <- afs

	af2s <- unique(pow$af2)
	#af2s <- afs

	betas3 <- unique(pow$beta3)
	samplesK.all <- 2 * unique(pow$n) / 1000

	par( mfrow=c(3,3) )

	for( af1 in af1s ) {
		for( af2 in af2s ) {
			first <- T
			col <- 1
			pch <- 15
			cols <- c()
			pchs <- c()
			b3s <- c()

			cat('af1:', af1, '\taf2:', af2, '\n')
			for( beta3 in betas3 ) {
				cat('\tbeta3:', beta3, '\n')

				keep <- (pow$af1 == af1) & (pow$af2 == af2) & (pow$beta3 == beta3) 
				pk <- pow[keep,]

				samplesK <- 2 * pk$n / 1000
				perc <- pk$perc

				# Remove last item if it is 100% power
				while((length(perc) > 2) && (perc[ length(perc) ] == 100) && (perc[ length(perc) - 1 ] == 100)) {
					n1 <- length(perc) - 1
					perc <- perc[1:n1]
					samplesK <- samplesK[1:n1]
				}

				# Anything to plot?
				if( length(perc) > 2 ) {
					# Add point at zero
					samplesK <- c(0, samplesK)
					perc <- c(0, perc)
					nonzero <- (perc > 0)
				
					if( first ) {
						title <- paste('Power    AF_1:', af1, '    AF_2:', af2)
						title <- paste('AF1:', af1, '    AF2:', af2)

						xlab <- 'Sample size [in thousands]'
						xlab <- ''

						ylab <- 'Power %'
						ylab <- ''

						plot( samplesK, perc, ylim=c(0,100), xlab=xlab, ylab=ylab, main=title, pch=pch, col=col, type='l')
						points( samplesK[nonzero], perc[nonzero], pch=pch, col=col)
						first <- F
					} else {
						points( samplesK[nonzero], perc[nonzero], pch=pch, col=col, cex=0.7)
						lines( samplesK, perc, pch=pch, col=col)
					}
				
					cols <- c(cols, col)
					pchs <- c(pchs, pch)
					b3s <- c(b3s, beta3)
				} else {
					cat('\t\tSkipped: Not enough points\n')
				}

				col <- col + 1
				if( col > 8 ) {
					col <- 1
					pch <- pch + 1
				}
			}

			# Show logend
			legend('bottomright', legend=paste('',b3s), pch=pchs, col=cols )
		}
	}
}

#-------------------------------------------------------------------------------
# Figure S1: Comparisson between PAM(1) and P[Qhat, t=1]
#-------------------------------------------------------------------------------
figureS1 <- function() {
	# Load Qhat
	Qhat <- read.table('Qhat.txt', header = TRUE, row.names = 1, sep="\t")
	Qhat <- as.matrix(Qhat)

	# Load AA frequencies
	aa.freq <- read.table('aa.frequencies.txt', header = F, row.names = 1, sep="\t")
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
	pam <- read.table('pam1.txt', header = TRUE, row.names = 1, sep="\t")
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
		s <- sum(pam[i,]) - pam[i, i]
		cat('Rowsum:[', i ,']', s, '\n')
		for( j in 1:d ) {
			err[i,j] <- abs(pam[i,j] - Pt0[i,j]) / s
		}
	}
	err[ is.nan(err) ] <- 0

	mypalette <- redgreen(100)
	mypalette <- colorRampPalette(c("red", "white", "green"))(n = 100)

	heatmap.2(Qhat, main = "Qhat", sub="Normalized by row", Rowv=F, Colv=F, col = mypalette, density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "row", na.rm=T); 

	#heatmap.2(err, main = "PAM1 vs P[Qhat, t=1]", sub="Error normalized by row", Rowv=F, Colv=F, col = mypalette, density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "none", na.rm=T); 
	heatmap.2(err, main = "", sub="Error normalized by row", Rowv=F, Colv=F, col = mypalette, density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = F, symbreaks = T, scale = "none", na.rm=T); 
}

#-------------------------------------------------------------------------------
# Figure S2: Show Q2's structure
#-------------------------------------------------------------------------------
figureS2 <- function() {
	Qhat2 <- read.table('Qhat2.txt', header = TRUE, row.names = 1, sep="\t")
	Qhat2<- as.matrix(Qhat2)
	q2 <- Qhat2
	diag(q2) <- 0
	mypalette <- colorRampPalette(c("red", "white", "black"))(n = 100)
	heatmap.2(q2, main = "", sub="Diagonal set to zero", Rowv=F, Colv=F, col = mypalette, density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = F, symbreaks = T, scale = "none", na.rm=T); 
}

#-------------------------------------------------------------------------------
# Figure S3: Gene-Gene interaction prediction. Distribution of $avg_3[LL(MSA)]$ 
#            for interacting and non-interacting genes showing that they can be 
#            differentiated.
#-------------------------------------------------------------------------------
figureS3 <- function(ll.alt, ll.null) {
	par(mfrow=c(1,1))

	for( wsize in 2:2 ) {
		keep <- ll.alt$wsize == wsize
		lla <- as.vector(ll.alt$ll[keep])

		keep <- ll.null$wsize == wsize
		lln <- as.vector(ll.null$ll[keep])

		wt <- wilcox.test(lla,lln,alternative='g')
		cat(wsize, '\t', wt$p.value, '\n')

		title <- paste("Best LL(MSA) using", wsize ,"amino acids")
		title <- ""

		xlim <- c(5,20)
		plot( density( lla ), col='red', main=title, xlab='', xlim=xlim )
		lines( density( lln ), col='green' )
	}
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

pngSize <- 1 * 1024
if( savePng )	png(width=pngSize, height=pngSize)

#---
# Figure 2
#---
if( fig2 ) {
	# Load data 
	if( !exists('lc') ) {
		cat('Reading file: likelihood.contact.values.txt\n')
		lc <- read.csv("likelihood.contact.values.txt", header=F)

		cat('Reading file: likelihood.null.values.txt\n')
		la <- read.csv("likelihood.null.values.txt", header=F)

		# Convert to numbers
    	lc <- as.numeric( lc[,1] )
    	la <- as.numeric( la[,1] )

		cat('Calculating probability ratio\n')
		xlim <- c(-40,20)
		lalcOddsX <- seq( min(xlim), max(xlim), 0.05)
		cr <- function(x) { cummProbRatio(x, lc, la) }
		lalcOdds <- sapply(lalcOddsX, cr)
	}

	figure2(lc, la, lalcOdds, lalcOddsX)
}

#---
# Figure 3
#---
if( fig3 ) {
	if( !exists('d.int') ) {
		windowSize <- 1							# Originally we calculated for different "windowSizes", but now we stick to windowSize=1
		cat('Window size:', windowSize, '\n')

		# Load data 'in contact'
		fileNeigh <- paste('likelihood.pdb_compound.neigh_', windowSize, '.3.0.txt', sep='')
		cat('Reading file:', fileNeigh, '\n')
		d.int <- read.table(fileNeigh, header=F, sep='\t')

		# Load data 'not in contact'
		fileNeigh <- paste('likelihood.pdb_compound.neigh_', windowSize, '.-30.0.txt', sep='')
		cat('Reading file:', fileNeigh, '\n')
		d.non <- read.table(fileNeigh, header=F, sep='\t')

		# Extract log-likelyhood
		ll.int <- d.int[,5]
		ll.non <- d.non[,5]

		cat('Calculating probability ratio\n')
		xlim <- c(-40,20)
		orX <- seq( min(xlim), max(xlim), 0.05)
		cr <- function(x) { cummProbRatio(x, ll.int, ll.non) }
		or <- sapply(orX, cr)
	}

	figure3(ll.int, ll.non, or, orX)
}

#---
# Other figures
#---
if( fig4 ) { figure4() }
if( fig5 ) { figure5() }
if( figS1 ) { figureS1() }

if( figS3 ) { 
	# Load data
	if( !exists('ll.genegene.alt') ) {
    	fileNull <- 'likelihoodGeneGeneMatrix.null.txt'
    	fileAlt <- 'likelihoodGeneGeneMatrix.alt.txt'
    	ll.genegene.alt <- read.table(fileAlt, sep="\t", header=TRUE)
    	ll.genegene.null <- read.table(fileNull, sep="\t", header=TRUE)
	}

	figureS3( ll.genegene.alt, ll.genegene.null) 
}

#---
# High definition plots
#---
if( savePng )	{
	# Reopen PNG plot using higher definition
	pngSize <- 5 * 1024
	dev.off()
	png(width=pngSize, height=pngSize)
}

if( figS2 ) { figureS2() }

# Close graphics device
if( savePng )	dev.off()

