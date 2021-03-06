
library('Matrix')
library('gplots')
library('expm')

savePng <- F

# Main figures
fig2    <- F
fig2b   <- F
fig3    <- F
fig4    <- F
fig5    <- F

# Suppleamentary
figS1   <- F
figS2   <- F
figS3   <- F
figS4   <- F

# Thesis final submission
figFS1  <- F
figFS2  <- T
figFS3  <- T

# Physical distance threshold to be considered as 'in-contact'
dist.th <- 8.0

# LL(MSA) threshold to be considered as 'in-contact'
llmsa.th <- 10.0

if( savePng) {
	heatmap.keysize <- 0.4
} else {
	heatmap.keysize <- 1.0
}

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
# True positive rate
#-------------------------------------------------------------------------------
truePositiveRate <- function(x, ll.int, ll.non) {
	return( sum( ll.int >= x ) / length( ll.int ) ) 
}

#-------------------------------------------------------------------------------
# True positive rate
#-------------------------------------------------------------------------------
falsePositiveRate <- function(x, ll.int, ll.non) {
	return( sum( ll.non >= x ) / length( ll.non ) )
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
# Figure 2b: ROC based on histogram from Fig2
#-------------------------------------------------------------------------------
figure2b <- function(fpr, tpr) {
	xlim <- c(0,1)
	plot(fpr, tpr, type = "l", col='red', xlab = "False positive rate", ylab = "True positive rate", xlim=xlim, ylim=xlim)
	abline(0,1, col='grey', lty=5)
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

						xlab <- ''
						xlab <- 'Sample size [in thousands]'

						ylab <- ''
						ylab <- 'Power %'

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

	# Calculate relative error
#	err <- pam * 0
#	for( i in 1:d ) {
#		s <- sum(pam[i,]) - pam[i, i]
#		cat('Rowsum:[', i ,']', s, '\n')
#		for( j in 1:d ) {
#			err[i,j] <- abs(pam[i,j] - Pt0[i,j]) / s
#		}
#	}

	err <- pam * 0
	# Mathieu wants log scale ratio (not relative error)
	for( i in 1:d ) {
		cat('Rowsum:[', i ,']\n')
		for( j in 1:d ) {
			err[i,j] <- log2(Pt0[i,j] / pam[i,j])
		}
		#err[i,i] <- 0
	}
	err[ is.nan(err) ] <- 0
	err[ is.infinite(err) ] <- 0

	mypalette <- redgreen(100)
	mypalette <- colorRampPalette(c("red", "white", "green"))(n = 100)

	heatmap.2(Qhat, main = "Qhat", sub="Normalized by row", Rowv=F, Colv=F, col = mypalette, density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = T, symbreaks = T, scale = "row", na.rm=T, keysize = heatmap.keysize); 

	#heatmap.2(err, main = "", sub="Ratio normalized by row", Rowv=F, Colv=F, col = mypalette, density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = F, symbreaks = T, scale = "none", na.rm=T, keysize = heatmap.keysize); 
	hmcols <- colorRampPalette(c("red","white","blue"))(256)
	heatmap.2(err, main = "", sub="Ratio normalized by row", Rowv=F, Colv=F, col = hmcols, density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = F, symbreaks = T, scale = "none", na.rm=T, keysize = heatmap.keysize); 
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
	heatmap.2(q2, main = "", sub="Diagonal set to zero", Rowv=F, Colv=F, col = mypalette, density.info = "none", trace = "none", dendrogram = "none", symm = F, symkey = F, symbreaks = T, scale = "none", na.rm=T, keysize = heatmap.keysize); 
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
		plot( density( lla ), col='red', main=title, xlab='log(Lc)', xlim=xlim )
		lines( density( lln ), col='green' )
	}
}

#-------------------------------------------------------------------------------
# Figure S4: Clinical significance categories
#
# Group by clinical significance (CLNSIG from CinVar)
#	Variant Clinical Significance: 
#		0 - Uncertain significance, 
#		1 - not provided, 
#		2 - Benign, 
#		3 - Likely benign, 
#		4 - Likely pathogenic, 
#		5 - Pathogenic, 
#		6 - drug response, 
#		7 - histocompatibility, 
#		255 - other
# So we create 3 groups:
#		Unknown: 0, 1, 255
#		Benign : 2, 3, 6 (something that has a drug response, is good)
#		Pathogenic: 4, 5
#
# Some stats:
#	CLNSIG: 0 		Count: 21 		LL_mean: 14.7234 	LL_median 9.826063 	Len_mean:  4352.1 	Len_meadian: 935 
#	CLNSIG: 1 		Count: 52 		LL_mean: 23.31171 	LL_median 13.61291 	Len_mean:  1435.9 	Len_meadian: 1481 
#	CLNSIG: 2 		Count: 272 		LL_mean: 34.10323 	LL_median 26.27931 	Len_mean:  3925.3 	Len_meadian: 1172.5 
#	CLNSIG: 3 		Count: 258 		LL_mean: 31.56365 	LL_median 25.56079 	Len_mean: 11619.7 	Len_meadian: 4575 
#	CLNSIG: 4 		Count: 562 		LL_mean: 17.53154 	LL_median 12.04729 	Len_mean:  2132.5 	Len_meadian: 838 
#	CLNSIG: 5 		Count: 4206 	LL_mean: 16.90444 	LL_median 11.73378 	Len_mean:  1057.6 	Len_meadian: 688.5 
#	CLNSIG: 6 		Count: 18 		LL_mean: 32.63433 	LL_median 22.10314 	Len_mean:   729.9 	Len_meadian: 492 
#	CLNSIG: 255 	Count: 10 		LL_mean: 20.09807 	LL_median 11.39337 	Len_mean:   4088 	Len_meadian: 767 
#
# 1000 Genomes	:		LL_mean: 23.2	LL_median: 14.8
# HGMD 			:		LL_mean: 19.8	LL_median: 12.1
# ClinVar		:		LL_mean: 18.6	LL_median: 11.9
#-------------------------------------------------------------------------------
figureS4 <- function(ll.alt, ll.null) {
	# Read data
	if( ! exists('d.clin') ) {
		fileName <- "likelihood.clinvar/top.sorted.clnsig.txt"
		d.clin = read.table(fileName, sep="\t", header=TRUE)
		#names(d.clin) <- c('chr', 'pos', 'ref', 'alt', 'll', 'len', 'clnsig')
	}

	#plot(density( d.clin$ll ), main="Log-Likelihood by Clinical Significance (CLNSIG)", xlab="Log-likelihood", sub="Black: All, Blue: Unknown, Red: Pathogenic, Green: Benign", xlim = c(0, 75) )	# All entries
	plot(density( d.clin$ll ), main='', xlab="log(Lc)", sub="", xlim = c(0, 40) )	# All entries

	# CLNSIG: Unknown
	keep <- (d.clin$clnsig == 0) | (d.clin$clnsig == 1) | (d.clin$clnsig == 255)
	lines(density( d.clin$ll[keep] ), col = 'blue')

	keep <- (d.clin$clnsig == 2) | (d.clin$clnsig == 3) | (d.clin$clnsig == 6)
	lines(density( d.clin$ll[keep] ), col = 'green')

	keep <- (d.clin$clnsig == 4) | (d.clin$clnsig == 4) | (d.clin$clnsig == 6)
	lines(density( d.clin$ll[keep] ), col = 'red')

	legend("topright", inset=.05, c('All', 'Unknown', 'Pathogenic', 'Benign' )
								, fill=c( rgb(0,0,0,1), rgb(0,0,1,1), rgb(1,0,0,1), rgb(0,1,0,1) )
								, horiz=F)
}

#-------------------------------------------------------------------------------
# Figure FS1: Thesis' final submission
#-------------------------------------------------------------------------------

figureFS1 <- function(fs) {
	# Add 'in-contact' criteria
	ic <- (fs$dist <= 3.0)

	# LL(MSA)
	minx <- -40
	maxx <- 25
	llmsa.ic  <- fs$llmsa[ic & (fs$llmsa >= minx) & (fs$llmsa <= maxx)]
	llmsa.nic <- fs$llmsa[!ic & (fs$llmsa >= minx) & (fs$llmsa <= maxx)]
	compareHist(llmsa.ic, llmsa.nic, "LL(MSA)", c(-40,30), "LL(MSA)", "Density", 100)
	llmsa.test <- wilcox.test(llmsa.ic, llmsa.nic)
	cat('LL(MSA):', llmsa.test$p.value, '\n')

	# McBASC. 
	# Note: Negative values (-1) indicate that the calculation cannot be perfromed, so we filter them out
	corr.ic  <- fs$corr[ic & (fs$corr >= 0) ]
	corr.nic <- fs$corr[!ic & (fs$corr >= 0) ]
	compareHist(corr.ic, corr.nic, "", c(0,1), "Correlation (McBASC)", "Density", 100)
	corr.test <- wilcox.test(corr.ic, corr.nic)
	cat('McBASC:', corr.test$p.value, '\n')

	# McBASC Fodor et. al. implementation
	# Note: Negative values (-1) indicate that the calculation cannot be perfromed, so we filter them out
	corr.ic  <- fs$corr.fodor[ic & (fs$corr.fodor >= 0) ]
	corr.nic <- fs$corr.fodor[!ic & (fs$corr.fodor >= 0) ]
	compareHist(corr.ic, corr.nic, "", c(0,1), "Correlation (McBASC)", "Density", 100)
	corr.test <- wilcox.test(corr.ic, corr.nic)
	cat('McBASC (Fodor):', corr.test$p.value, '\n')

	# MI 
	mi.ic  <- fs$mi[ic & (fs$mi > 0)  & (fs$mi < 1)]
	mi.nic <- fs$mi[!ic & (fs$mi > 0) & (fs$mi < 1) ]
	compareHist(mi.ic, mi.nic, "", c(0,1), "Mutual information", "Density", 100)
	mi.test <- wilcox.test(mi.ic, mi.nic)
	cat('MI:', mi.test$p.value, '\n')

	# Variation of information
	vi.ic  <- fs$vi[ic & (fs$vi >= 0) ]
	vi.nic <- fs$vi[!ic & (fs$vi >= 0) ]
	compareHist(vi.ic, vi.nic, "VI", c(0,4.5), "VI", "Density", 100)
	vi.test <- wilcox.test(vi.ic, vi.nic)
	cat('VI:', vi.test$p.value, '\n')

	# Variation of information
	vi.ic  <- log2( fs$vi[ic & (fs$vi > 0) ] )
	vi.nic <- log2( fs$vi[!ic & (fs$vi > 0) ] )
	compareHist(vi.ic, vi.nic, "VI", c(-5,2.2), "VI", "Density", 100)
	vi.test <- wilcox.test(vi.ic, vi.nic)
	cat('VI log:', vi.test$p.value, '\n')

	# McBASC Fodor et. al. implementation
	# Note: Negative values (-1) indicate that the calculation cannot be perfromed, so we filter them out
	corr.ic  <- log2( fs$corr.fodor[ic & (fs$corr.fodor >= 0) ] )
	corr.nic <- log2( fs$corr.fodor[!ic & (fs$corr.fodor >= 0) ] )
	compareHist(corr.ic, corr.nic, "Log McBASC (Fodor)", c( -15,1), "Log McBASC (Fodor)", "Density", 100)
	corr.test <- wilcox.test(corr.ic, corr.nic)
	cat('Log McBASC (Fodor):', corr.test$p.value, '\n')
}

#-------------------------------------------------------------------------------
# Figure FS2: Thesis' final submission
#-------------------------------------------------------------------------------

figureFS2 <- function(pdb.id, d, ic.paper) {
	# Add 'in-contact' criteria
	ic <- (d$dist <= dist.th)

	#---
	# Create images from data
	#---
	m.size <- max( d$aa1.pos, d$aa2.pos ) 
	md <- matrix(NA, m.size, m.size)
	mll <- matrix(0, m.size, m.size)
	mll.bin <- matrix(0, m.size, m.size)
	mll.top <- matrix(0, m.size, m.size)
	mic.paper <- matrix(0, m.size, m.size)

	#---
	# Create images from data
	#---
	len <- dim(d)[1]
	max.dist <- max(d$dist)
	#cat("Matrix size:", m.size, "\tmax.dist:", max.dist, "\tdata.len:", len, "\n")
	for( i in seq(1, len) ) {
		dd <- d[i,]
		p1 <- dd$aa1.pos
		p2 <- dd$aa2.pos
		dist <- dd$dist

		md[p1,p2] <- dist
		md[p2,p1] <- dist

		mll[p1,p2] <- dd$llmsa
		mll[p2,p1] <- dd$llmsa
	}
	mic <- md <= dist.th

	# Top N LL(MSA) values
	N <- 50
	topN <- sort(d$llmsa, decreasing=T)[N]
	#cat('Top', N,' value', topN, '\tcount:', sum(mll >= topN), '\n')
	mll.top[ mll >= topN ] <- 1
	
	# Binary 'mll' matrix
	mll.bin[ mll >= llmsa.th ] <- 1

	#---
	# Create images from paper's data
	#---
	len <- dim(ic.paper)[1]
	for( i in seq(1, len) ) {
		dd <- ic.paper[i,]
		p1 <- dd$pos1
		p2 <- dd$pos2
		mic.paper[p1,p2] <- 1
		mic.paper[p2,p1] <- 1
	}

	#---
	# Show images
	#---
	md <- md / max.dist
	md <- 1 - md	# Number closer to 1.0 indicates 'in contact' (i.e. closer ditance)
	md[ is.na(md) ] <- 0	# Fill missing data

	# Scale mll to [0, 1]
	mll <- (mll - min(mll)) / (max(mll) - min(mll))	

	col <- grey(seq(1, 0, length = 256))
	image(md, axes = FALSE, col = col, main=paste(pdb.id, 'Distance'))
	image(mic, axes = FALSE, col = col, main=paste(pdb.id, 'In contact (', dist.th, ')'))
	image(mll, axes = FALSE, col = col, main=paste(pdb.id, 'LL(MSA)'))
	image(mll.bin, axes = FALSE, col = col, main=paste(pdb.id, 'LL(MSA) >=', llmsa.th))
	image(mll.top, axes = FALSE, col = col, main=paste(pdb.id, 'Top', N,'[ LL(MSA) ]'))
	image(mic.paper, axes = FALSE, col = col, main=paste(pdb.id, 'In contact (Paper)'))

	return( list(md=md, mic=mic, mll=mll, mll.bin=mll.bin, mll.top=mll.top, mic.paper=mic.paper) )
}

#-------------------------------------------------------------------------------
# Compare results to paper
#-------------------------------------------------------------------------------
compareToPaper <- function(pid, res) {
	cat('\n')
	tit <- paste(pid, dist.th, llmsa.th, sep='\t')
	compareTpFp( paste(tit, 'Top LL(MSA)', sep='\t'), res$mll.top, res$mic)
	compareTpFp( paste(tit, 'LL(MSA) >=', llmsa.th, sep='\t'), res$mll.bin, res$mic)
	compareTpFp( paste(tit, 'paper', sep='\t'), res$mic.paper, res$mic)
}

compareTpFp <- function(title, m, mic) {
	count <- sum( m == 1 )
	tp <- (!is.na(mic)) & (!is.na(m)) & (m == 1) & mic 
	#cat('\tMatches', title, ':', sum(tp), '\trate:', (sum(tp)/count), '\n')

	fp <- (!is.na(mic)) & (!is.na(m)) & (m == 1) & (!mic) 
	#cat('\tFalse positives', title, ':', sum(fp), '\trate:', (sum(fp)/count), '\n')

	cat(title, sum(tp), (sum(tp)/count), sum(fp), (sum(fp)/count), '\n', sep='\t')
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

pngSize <- 1 * 1024
if( savePng )	png(width=pngSize, height=pngSize)

#---
# Figure 2
#---
if( fig2 || fig2b ) {
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

	if( !exists('tpr') ) {
		cat('Calculating true positive rate\n')
		tp.rate <- function(x) { truePositiveRate(x, lc, la); }
		tpr <- sapply(lalcOddsX, tp.rate)
	}

	if( !exists('fpr') ) {
		cat('Calculating false positive rate\n')
		fp.rate <- function(x) { falsePositiveRate(x, lc, la); }
		fpr <- sapply(lalcOddsX, fp.rate)
	}

	if( fig2 )	{ figure2(lc, la, lalcOdds, lalcOddsX) }
	if( fig2b )	{ figure2b(fpr, tpr) }
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

if( figS4 ) { figureS4() }

#---
# Final submission figures
#---
if( figFS1 ) {
	if( !exists('fs') ) {
		fileName <- "fs_R.txt"
		cat('Loading', fileName, '\n')
		fs = read.table(fileName, sep="\t", header=TRUE)
		cat('Done\n')
	}

	figureFS1(fs)
}

if( figFS2 ) {
	if( !exists('pdb.1a17') ) {
		fileName <- "1a17_R.txt"
		cat('Loading', fileName, '\n')
		pdb.1a17 = read.table(fileName, sep="\t", header=TRUE)

		fileName <- "1a17_paper.txt"
		cat('Loading', fileName, '\n')
		pdb.1a17.paper = read.table(fileName, sep="\t", header=TRUE)

		cat('Done\n')
	}

	pid <- '1A17'
	res.1a17 <- figureFS2(pid, pdb.1a17, pdb.1a17.paper)

	# Compare to paper
	compareToPaper(pid, res.1a17)
}

if( figFS3 ) {
	if( !exists('pdb.1ubi') ) {
		fileName <- "1ubi_R.txt"
		cat('Loading', fileName, '\n')
		pdb.1ubi = read.table(fileName, sep="\t", header=TRUE)

		fileName <- "1ubi_paper.txt"
		cat('Loading', fileName, '\n')
		pdb.1ubi.paper = read.table(fileName, sep="\t", header=TRUE)

		cat('Done\n')
	}

	pid <- "1UBI"
	res.1ubi <- figureFS2(pid, pdb.1ubi, pdb.1ubi.paper)

	# Compare to paper
	compareToPaper(pid, res.1ubi)
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


for( dist.th in c(3, 5, 8) ) {
	for( llmsa.th in c(0, 3, 5, 8, 10, 12, 15) ) {
		pid <- '1A17'
		res.1a17 <- figureFS2(pid, pdb.1a17, pdb.1a17.paper)
		compareToPaper(pid, res.1a17)

		pid <- "1UBI"
		res.1ubi <- figureFS2(pid, pdb.1ubi, pdb.1ubi.paper)
		compareToPaper(pid, res.1ubi)
	}
}
