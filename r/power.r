
savePlot = T

# p-value threshold using Bonferoni correction over 1M variants
pvalue.th <- 0.05 / ( 10^6 )^2

# Load data
p <- read.table('power.txt', header=T, sep='\t')

# Filter those having more than 800 lines (80% of 1000 tests)
keep <- p$line > 800
pp <- p[keep,]

# Does the logistic test have 'enough power'
pp$ok <- pp$p80 < pvalue.th
pp$af12 <- pp$af1 * pp$af2

#---
# Plots
#---

if( savePlot ) { png( width=1024, height=1024 ) }

for( dots in c(TRUE, FALSE) ) {
	# Plot without axis so we can add the propper axis later (xaxt='n')
	x <- (pp$n)
	y <- log10(pp$af12)
	plot( x, y, cex=0.01
			, xlab='Number of cases'
			, ylab='AF1 * AF2'
			, main='Power (80% detection)'
			, sub=paste('Note: Curves show exponential fit from simulated values. p-value threshold:'
			, pvalue.th)
			, xaxt = "n"
			, yaxt = "n"
			)

	# Add custom axis
	axis(1, sort(unique(x)))
	axis(2, sort(unique(y)), labels=sort(unique(pp$af12)))

	# Plot for each 'beta3' level
	b3s <- sort(unique(pp$beta3))
	l2 <- length(b3s)/2
	for( i in 1:length(b3s) ) {
		b3 <- b3s[i]
		k <- pp$ok & (pp$beta3 == b3)
		kn <- !pp$ok & (pp$beta3 == b3)

		# Point type
		pch <- 14 + i 

		# Some offset
		dx <- 0 * (i-l2)
		dy <- 0.03 * (i-l2)

		if( dots ) {
			#---
			# Plot dost on each beta level
			#---
			points( dx + x[k] , dy + y[k] , col='green', pch=pch )
			points( dx + x[kn], dy + y[kn], col='red',   pch=pch )
		} else {
			#---
			# Interpolate an exponential curve
			#---

			# Find minimum 'power' for each beta3 anf each 'n'
			xx <- c()
			yy <- c()
			for( nn in sort(unique(pp$n)) ) {
				af <- min( pp[k & (pp$n == nn),'af12'] )
				cat(b3, nn, af, '\n')
				if( is.finite(af) ) {
					xx <- c(xx, nn)
					yy <- c(yy, af)
				}
			}
        
			# Any results? Plot them using a exponential interpolation
			if( length(xx) > 0 ) {
				# Try to fit an exponential fit
				lyy <- -log10(yy)
				fit <- lm( log(lyy) ~ xx )
        
				# Create a smooth curve using exponential fit
				xx.fit <- seq(min(xx), max(xx), 1000)
				lyy.fit <- exp(fit$coefficients[1] + fit$coefficients[2] * xx.fit)
        
				# Plot lines
				lines( xx.fit, -lyy.fit, pch=pch, lty=i, col=i)
        
				# Select some points form the previous lines and show 'dots'
				k <- ( 1:length(xx.fit) ) %% 50 == 0
				points( xx.fit[k], -lyy.fit[k] , col=i, pch=pch )
			}
		}
	}

	legend('bottomleft', legend=paste('beta3 =',b3s), pch=14 + (1:length(b3s)), col=1:length(b3s) )
}

if( savePlot )	{ dev.off() } 

