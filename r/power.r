
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


# Plot
x <- log10(pp$n)
y <- log10(pp$af12)
plot( x, y, cex=0.01, xlab='log10(cases)', ylab='log10(AF1 * AF2)', main='Power (80% detection)', sub=paste('p-value threshold:', pvalue.th))

b3s <- unique(pp$beta3)
l2 <- length(b3s)/2
for( i in 1:length(b3s) ) {
	b3 <- b3s[i]
	k <- pp$ok & (pp$beta3 == b3)
	kn <- !pp$ok & (pp$beta3 == b3)

	pch <- 14 + i 
	dx <- 0.01 * (i-l2)
	points( dx + x[k] , y[k] , col='green', pch=pch )
	points( dx + x[kn], y[kn], col='red',   pch=pch )
}

legend('bottomleft', legend=paste('beta3 =',b3s), pch=14 + (1:length(b3s)), col='green')
