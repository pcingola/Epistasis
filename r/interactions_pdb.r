
savePng <- F

cummProbRatio <- function(x) {
	p.int <- sum( ll.int >= x ) / length( ll.int ) 
	p.non <- sum( ll.non >= x ) / length( ll.non ) 
	return(p.int / p.non)
}

probRatio <- function(x) {
	xmin <- x - 0.5
	xmax <- x + 0.5

	keep <- (xmin <= ll.int) & (ll.int <= xmax)
	p.int <- sum( keep ) / length( ll.int ) 

	keep <- (xmin <= ll.non) & (ll.non <= xmax)
	p.non <- sum( keep ) / length( ll.non ) 

	return(p.int / p.non)
}

if( savePng )	png(width=1024, height=1024)

for( neigh in 1:1 ) {
	cat('\nNeighbours:', neigh, '\n')
	d.int <- read.table(paste('likelihood.pdb_compound.neigh_', neigh, '.3.0.txt', sep=''), header=F, sep='\t')
	d.non <- read.table(paste('likelihood.pdb_compound.neigh_', neigh, '.-30.0.txt', sep=''), header=F, sep='\t')

	ll.int <- d.int[,5]
	ll.non <- d.non[,5]
	xlim <- c(-20,15)
	title <- 'LL(MSA) of interacting proteins (PDB)'
	xlab <- 'LL(MSA)'
	sub <- paste('Neighbours:', 0) # neigh * 2 + 1)

#	# Use MI
#	ll.int <- d.int[,15]
#	ll.int <- ll.int[ ll.int > 0 ]
#	ll.non <- d.non[,15]
#	ll.non <- ll.non[ ll.non > 0 ]
#	xlim <- c(0, 1)
#	title <- 'Mutual information of interacting proteins (PDB)'
#	xlab <- 'MI'
#	sub <- ''

#	# Use variation of information
#	ll.int <- d.int[,16]
#	ll.int <- ll.int[ ll.int > 0 ]
#	ll.non <- d.non[,16]
#	ll.non <- ll.non[ ll.non > 0 ]
#	xlim <- c(0, 3.1)

	plot( density( ll.int ), main=title, sub=sub, xlim=xlim, xlab=xlab, col='red' )
	cat('Interacting\t', summary(ll.int), '\n')

	lines( density( ll.non ), col='green' )
	cat('Non-interacting\t', summary(ll.non), '\n')

	legend("topleft", inset=.05, c('Interacting', 'Non-interacting'), fill=c('red','green'), horiz=F)
}


# cols <- 1:dim(d)[2]
# plot( density( d[,1] ), main='LL(MSA) of interacting proteins (PDB)', xlim= c(-20,15), xlab='LL(MSA)' )
# for( col in cols ) {
# 	name <- names(d)[col]
# 	cat( name, '\t', summary(d[,col]), '\n' )
# 	lines( density( d[,col] ), col=col )
# }
# 
# legend("topleft", inset=.05, names(d), fill=cols, horiz=F)
# 
# if( savePng )	dev.off()

#---
# Plot odds ratio ('density' probability)
#---
x <- seq(-10,10,0.05)
y <- sapply(x, probRatio)
ly <- log(y) 
plot( x, ly, main='Log-Odds ratio', sub='log{ P[ LL(MSA|M1) > X ] / P[ LL(MSA|M0) > X ] }', cex=0.5, xlab='LL(MSA)', ylab='log( Odds )')
lines( supsmu(x, ly), col='red' )

#---
# Plot odds ratio (cummulative probability)
#---
x <- seq(-10,10,0.05)
y <- sapply(x, cummProbRatio)
ly <- log(y) 
plot( x, ly, main='Log-Odds ratio (Cumm Prob)', sub='log{ P[ LL(MSA|M1) > X ] / P[ LL(MSA|M0) > X ] }', cex=0.5, xlab='LL(MSA)', ylab='log( Odds )', ylim=c(0, 3))
lines( supsmu(x, ly), col='red' )

lly <- log(log(y))
llymodel <- lm( lly ~ x )


# Close graphics device
if( savePng )	dev.off()
