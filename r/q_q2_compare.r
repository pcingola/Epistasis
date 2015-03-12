
#-------------------------------------------------------------------------------
# Show density
#-------------------------------------------------------------------------------
showDens <- function(y, title) {
	l <- log(y)
    cat('Raw:\t\tMean:', mean(y), ' StdDev:', sd(y), ' Median:', median(y), '\n')
    cat('Log:\t\tMean:', mean(l), ' StdDev:', sd(l), ' Median:', median(l), '\n')

    plot( density(l), main=title, xlab='Log[ R(ab,cd) ], Green: median, Blue: mean')
	abline( v=median(l), col='green', lty=2, lwd=2);
	abline( v=mean(l), col='blue', lty=2, lwd=2);
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

savePlot <- T

if( savePlot ) { png( width=1024, height=1024 ) }

#---
# Load data
#---
if( ! exists('qcmp') ) {
    cat('Loading Q_Q2 compare data\n')
    qcmp <- read.table('q_q2_compare.ratio.txt', sep='\t', header=F)
}

r <- as.vector(unlist(qcmp))
showDens(r, 'Histogram of Log[R(ab,cd)] ratios')


if( savePlot ) { dev.off(); }
