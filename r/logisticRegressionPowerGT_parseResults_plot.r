
savePlot <- T

if( ! exists('pow') ) {
	pow <- read.table('logisticRegressionPowerGT_parseResults.txt', sep="\t", header=TRUE)
}

af1s <- unique(pow$af1)
af2s <- unique(pow$af2)
betas3 <- unique(pow$beta3)
samplesK.all <- 2 * unique(pow$n) / 1000

if( savePlot ) { png( width=1024, height=1024 ) }

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
					plot( samplesK, perc, ylim=c(0,100), xlab='Sample size [in thousands]', ylab='Power %', main=title, pch=pch, col=col, type='l')
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

if( savePlot )  { dev.off() }

