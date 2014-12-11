
savePng <- T

if( savePng )	png(width=1024, height=1024)
for( neigh in 1:4 ) {
	cat('\nNeighbours:', neigh, '\n')
	d.int <- read.table(paste('likelihood.pdb_compound.neigh_', neigh, '.3.0.txt', sep=''), header=F, sep='\t')
	d.non <- read.table(paste('likelihood.pdb_compound.neigh_', neigh, '.-30.0.txt', sep=''), header=F, sep='\t')

	ll.int <- d.int[,5]
	ll.non <- d.non[,5]

	sub <- paste('Neighbours:', neigh * 2 + 1)
	plot( density( ll.int ), main='LL(MSA) of interacting proteins (PDB)', sub=sub, xlim= c(-20,15), xlab='LL(MSA)', col='red' )
	cat('Interacting\t', summary(ll.int), '\n')

	lines( density( ll.non ), col='green' )
	cat('Non-interacting\t', summary(ll.non), '\n')

	legend("topleft", inset=.05, c('Interacting', 'Non-interacting'), fill=c('red','green'), horiz=F)
}
if( savePng )	dev.off()


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
