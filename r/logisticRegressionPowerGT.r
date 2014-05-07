
library(epicalc)

n <- 1000					# Samples
m <- 100					# Number of iterations to get an aerage (200)

beta <- c(2, 0, 0, 2)		# Parameters
N <- length(beta) - 2		# Dimensions

af <- 0.1					# Allele frequency
debug <- TRUE

#-------------------------------------------------------------------------------
# Random genotype
#-------------------------------------------------------------------------------
randGt <- function(n, af) {
	r <- runif(n)

	gt <- rep(0, n)
	gt[ r < af ] <- 1
	gt[ r < af^2 ] <- 2

	return(gt)
}

#-------------------------------------------------------------------------------
# Perform Logistic regression test
#-------------------------------------------------------------------------------

testLr <- function(n, beta, af) {
	# Create input
	var1 <- randGt(n, af)
	var2 <- randGt(n, af)
	var12 <- abs( var1 - var2 )
	ones <- rep(1, n)
	X <- cbind( ones, var1, var2, var12 )	# As a matrix

	# Logit(pi)
	logitp = X %*% beta
	p = 1 / ( 1 + exp(-logitp) )

	# P(yi | X)
	r <- runif(n)
	y <- rep(0, n)
	y[ r <= p ] = 1

	# Logistic regression
	lr0 <- glm( y ~ var1 + var2 + var12 , family=binomial) 
	lr1 <- glm( y ~ var1 + var2         , family=binomial) 

	# Likelyhood ratio test
	lrt <- lrtest(lr0, lr1)
	pvalueLr <- lrt$p.value   # p-value from likelihood ration test
	lrSum <- summary( lr0 )
	pvalueWald <- lrSum$coefficients[4,4]
	if( debug )	cat('\t\tp-value (LR):', pvalueLr, '\tp-value (Wald):', pvalueWald, '\n')

	return( c(pvalueLr, pvalueWald) )
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

# for( n in seq(1000, 10000, 100) ) {
for( n in seq(1000, 1000, 100) ) {
	pvals <- matrix( rep(0, 2*m), ncol=2 )
	for( i in 1:m ) {
		pvals[i,] <- testLr(n, beta, af)
	}	
	cat('Iterations:', m, '\tSize:', n, '\tMean p-value (LR):', mean(pvals[,1]), '\tp-value (Wald):', mean(pvals[,2]), '\n')
}

