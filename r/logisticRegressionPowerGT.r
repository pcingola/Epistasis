
library(epicalc)

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

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

#---
# Parse command line arguments
#---

if( !exists('cmdLineArgs') ) cmdLineArgs <- commandArgs(trailingOnly = TRUE);

# Stop if there are no command line argument
if( length(cmdLineArgs) < 1 ) { fatalError('No command line arguments!\n'); }

iter  <- as.integer(cmdLineArgs[1]);	# Number of tests
n     <- as.integer(cmdLineArgs[2]);	# Number of cases and controls
af1   <- as.numeric(cmdLineArgs[3]);	# Allele frequency for variant 1
af2   <- as.numeric(cmdLineArgs[4]);	# Allele frequency for variant 2
beta3 <- as.numeric(cmdLineArgs[5]);	# Interaction term
prev  <- as.numeric(cmdLineArgs[6]);	# Disease prevalecense: e.g. 8% for type II diabetes

cat('Parameters\n')
cat('\titer  :', iter, '\n')
cat('\tn     :', n, '\n')
cat('\taf1   :', af1, '\n')
cat('\taf2   :', af2, '\n')
cat('\tbeta3 :', beta3, '\n')
cat('\tprev  :', prev, '\n')
cat('\tprev  :', (1-af1), '\n')

#---
# Initial parameters
#---
beta0 = log(prev/(1-prev))
beta1 = 0		# Other parameters are assumed to be zero
beta2 = 0
beta <- c(beta0, beta1, beta2, beta3)

#---
# Test logistic regression 'iter' times
#---
for( i in 1:iter ) {
	testLr(i/iter, n, af1, af2, beta)
}

