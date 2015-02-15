
savePng <- F

#-------------------------------------------------------------------------------
# Compare two distributions
#-------------------------------------------------------------------------------
compareHist <- function(data1, data2, main, xlim, xlab, ylab, breaks) {
	dens1 <- density(data1)
	h1 <- hist(data1, main=main, xlim=xlim, xlab=xlab, ylab=ylab, freq = F, breaks=breaks, col=rgb(1,0,0,1/4) , add=F );
	lines(dens1, col='red', xlim=xlim)

	dens2 <- density(data2)
	h2 <- hist(data2, xlim=xlim, xlab=xlab, ylab=ylab, freq = F, breaks=breaks, col=rgb(0,0,1,1/4), add=T );
	lines(dens2, col='blue', xlim=xlim)
}

#-------------------------------------------------------------------------------
# Cummulative probability ratio
#-------------------------------------------------------------------------------
cummProbRatio <- function(x) {
	p.int <- sum( ll.int >= x ) / length( ll.int ) 
	p.non <- sum( ll.non >= x ) / length( ll.non ) 
	return(p.int / p.non)
}

#-------------------------------------------------------------------------------
# Probability ratio
#-------------------------------------------------------------------------------
probRatio <- function(x) {
	xmin <- x - 0.5
	xmax <- x + 0.5

	keep <- (xmin <= ll.int) & (ll.int <= xmax)
	p.int <- sum( keep ) / length( ll.int ) 

	keep <- (xmin <= ll.non) & (ll.non <= xmax)
	p.non <- sum( keep ) / length( ll.non ) 

	return(p.int / p.non)
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

pngSize <- 1024
if( savePng )	png(width=pngSize, height=pngSize)

#---
# Figure 3.A: Compare LL(MSA) for proteins 'in contact' versus 'not in contact' (from PDB co-crystalized structures)
#---

# Data calculated using Epistasis.jar 
# Result files:
#		./data/interactions/pdb/likelihood.pdb_compound.neigh_1.*.0.txt.gz

windowSize <- 1							# Originally we calculated for different "windowSizes", but now we stick to windowSize=1
cat('Window size:', windowSize, '\n')

# Load data 'in contact'
fileNeigh <- paste('likelihood.pdb_compound.neigh_', windowSize, '.3.0.txt', sep='')
cat('Reding file:', fileNeigh, '\n')
d.int <- read.table(fileNeigh, header=F, sep='\t')

# Load data 'not in contact'
fileNeigh <- paste('likelihood.pdb_compound.neigh_', windowSize, '.-30.0.txt', sep='')
cat('Reding file:', fileNeigh, '\n')
d.non <- read.table(fileNeigh, header=F, sep='\t')

# Extract log-likelyhood
ll.int <- d.int[,5]
ll.non <- d.non[,5]

#---
# Show density distributions
#---
xlim <- c(-20,15)
xlim <- c(-40,20)
title <- 'LL(MSA) of interacting proteins (PDB)'
xlab <- 'LL(MSA)'
breaks <- 100

cat('Distributions summary:\n')
cat('\tInteracting\t', summary(ll.int), '\n')
cat('\tNon-interacting\t', summary(ll.non), '\n')

compareHist(ll.int, ll.non, title, xlim, xlab, "Frequency", breaks)
legend("topleft", inset=.05, c('Interacting', 'Non-interacting'), fill=c('red','green'), horiz=F)

#---
# Figure 3.B: Plot odds ratio (cummulative probability)
#---
#x <- seq( -10, 10, 0.05)
x <- seq( min(xlim), max(xlim), 0.05)
y <- sapply(x, cummProbRatio)
ly <- log(y) 

# Remove inacurrate numbers (too few points to calculate stats)
maxShow <- 10.7
y[ x > maxShow ] <- NA
ly[ x > maxShow ] <- NA

plot( x, ly, main='Log-Odds ratio (Cumm Prob)', sub='log{ P[ LL(MSA|M1) > X ] / P[ LL(MSA|M0) > X ] }', cex=0.5, xlab='LL(MSA)', ylab='log( Odds )', ylim=c(0, 3))
lines( supsmu(x, ly), col='red' )

# Close graphics device
if( savePng )	dev.off()
