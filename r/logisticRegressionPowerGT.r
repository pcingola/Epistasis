
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
# Logistic regression test
#-------------------------------------------------------------------------------
testLr <- function(perc, n, af1, af2, beta) {
	# Create random genotypes
	m <- n / prev	# We need to generate a lot of samples
	var1 <- randGt(m, af1)
	var2 <- randGt(m, af2)
	var12 <- var1 * var2
	ones <- rep(1, m)

	# Compute logistic model
	X <- cbind( ones, var1, var2, var12 )
	logitp = X %*% beta
	p = 1 / ( 1 + exp(-logitp) )

	# Assign disease estimates
	r <- runif(m)
	y <- rep(0, m)
	y[ r <= p ] = 1

	# Select controls
	mm <- 1:m
	cases.idx <- (mm[y==1])[1:n]
	ctrls.idx <- (mm[y==0])[1:n]

	# Create a matrix of cases and controls
	CC <- rbind(X[ctrls.idx,] , X[cases.idx,])
	Y <- c( y[ctrls.idx], y[cases.idx] )

	# Logistic regression
	lr0 <- glm( Y ~ CC[,2] + CC[,3]          , family=binomial)		# Null model
	lr1 <- glm( Y ~ CC[,2] + CC[,3] + CC[,4] , family=binomial)		# Alt model

	# Likelyhood ratio test
	lrt <- lrtest(lr0, lr1)
	pvalueLr <- lrt$p.value   # p-value from likelihood ratio test

	af12 <- sum( var12 > 0 ) / length(var12)
	cat(perc * 100, '%\tn:', n, '\taf1:', af1, '\taf2:', af2, '\taf12:', af12, '\tbeta3:', beta3,'\tp-value (LR):', pvalueLr, '\tbeta.null:', lr0$coefficients, '\tbeta.alt:', lr1$coefficients, '\tprevalecense:', (sum(y)/length(y)), '\n')
}

#  #-------------------------------------------------------------------------------
#  # Perform Logistic regression test
#  #-------------------------------------------------------------------------------
#  
#  testLr <- function(n, beta, af) {
#  	# Create input
#  	var1 <- randGt(n, af)
#  	var2 <- randGt(n, af)
#  	var12 <- abs( var1 - var2 )
#  	ones <- rep(1, n)
#  	X <- cbind( ones, var1, var2, var12 )	# As a matrix
#  
#  	# Logit(pi)
#  	logitp = X %*% beta
#  	p = 1 / ( 1 + exp(-logitp) )
#  
#  	# P(yi | X)
#  	r <- runif(n)
#  	y <- rep(0, n)
#  	y[ r <= p ] = 1
#  
#  	# Logistic regression
#  	lr0 <- glm( y ~ var1 + var2 + var12 , family=binomial) 
#  	lr1 <- glm( y ~ var1 + var2         , family=binomial) 
#  
#  	# Likelyhood ratio test
#  	lrt <- lrtest(lr0, lr1)
#  	pvalueLr <- lrt$p.value   # p-value from likelihood ratio test
#  	lrSum <- summary( lr0 )
#  	pvalueWald <- lrSum$coefficients[4,4]
#  	if( debug )	cat('\t\tp-value (LR):', pvalueLr, '\tp-value (Wald):', pvalueWald, '\n')
#  
#  	return( c(pvalueLr, pvalueWald) )
#  }

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

if( F ) {
	# for( n in seq(1000, 10000, 100) ) {
	for( n in seq(1000, 1000, 100) ) {
		pvals <- matrix( rep(0, 2*m), ncol=2 )
		for( i in 1:m ) {
			pvals[i,] <- testLr(n, beta, af)
		}	
		cat('Iterations:', m, '\tSize:', n, '\tMean p-value (LR):', mean(pvals[,1]), '\tp-value (Wald):', mean(pvals[,2]), '\n')
	}
}

#---
# Initial parameters
#---
iter <-	100		# Number of tests
n <- 1000		# Number of cases and controls

af1 = 0.1		# Allele frequency for variant 1
af2 = 0.1		# Allele frequency for variant 2

prev = 0.08		# Disease prevalecense 8% (type II diabetes)
beta0 = log(prev/(1-prev))

# Other parameters
beta1 = 0
beta2 = 0
beta3 = 2		# Interaction parameter

for( beta3 in c(0.1, 0.5, 1, 2, 5) ) {
	for( n in c(1000000, 100000, 50000, 25000, 10000, 1000) ) {
		for( af1 in seq(0.01,0.1,0.01) ) {
			for( af2 in seq(0.01,0.1,0.01) ) {
				beta <- c(beta0, beta1, beta2, beta3)

				for( i in 1:iter ) {
					testLr(i/iter, n, af1, af2, beta)
				}
			}
		}
	}
}

