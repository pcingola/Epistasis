
savePng <- TRUE

if( savePng ) png(width = 1024, height = 1024)

fileNull <- 'likelihood_aa5.null.avgLL.txt'
fileAlt <- 'likelihood_aa5.alt.avgLL.txt'

ll.alt <- read.table(fileAlt, sep="\t", header=FALSE)
ll.null <- read.table(fileNull, sep="\t", header=FALSE)

ll.alt <- as.vector(ll.alt[,1])
ll.null <- as.vector(ll.null[,1])

plot( density( ll.alt ), col='red', main="LL(MSA) using 5 amino acids", xlab='' )
lines( density( ll.null ), col='green' )


if( savePng )	dev.off()
