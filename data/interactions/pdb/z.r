
savePng <- T
if( savePng )	png(width=1024, height=1024)

d <- read.table( 'll.txt', header=T, sep='\t')

cols <- 1:dim(d)[2]
plot( density( d[,1] ), main='LL(MSA) of interacting proteins (PDB)', xlim= c(-20,15), xlab='LL(MSA)' )
for( col in cols ) {
	name <- names(d)[col]
	cat( name, '\t', summary(d[,col]), '\n' )
	lines( density( d[,col] ), col=col )
}

legend("topleft", inset=.05, names(d), fill=cols, horiz=F)

if( savePng )	dev.off()
