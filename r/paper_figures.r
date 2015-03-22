
savePng <- T

fig2 <- T
fig3 <- T
fig4 <- T

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
# Main
#-------------------------------------------------------------------------------

pngSize <- 1024
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
# Figure 4
#---
if( fig4 ) {
	figure4()
}

# Close graphics device
if( savePng )	dev.off()

