
#-------------------------------------------------------------------------------
#
# Logistic regression and likelihood ratio test
# Used for debugging data from LogisticeRegression class (Java code)
#
#
# Sample output:
#	Loading samples logReg_test_IRWLS_01.txt 
#	Beta (Real ): -0.75 -3 0.5 	LL:  -383.3362 
#	Beta (R    ): -0.7075655 -2.715519 0.4777281 	LL:  -382.0379 
#	Beta (IRWLS): 0 0 0 	LL:  -693.1472 
#	Beta (IRWLS): -0.3819791 -1.359302 0.2420161 	LL:  -429.0905 
#	Beta (IRWLS): -0.5784837 -2.167474 0.3780719 	LL:  -388.0521 
#	Beta (IRWLS): -0.6831148 -2.612633 0.4582468 	LL:  -382.2269 
#	Beta (IRWLS): -0.7066294 -2.711639 0.4769707 	LL:  -382.0381 
#	Beta (IRWLS): -0.7075641 -2.715513 0.4777269 	LL:  -382.0379 
#	Beta (IRWLS): -0.7075655 -2.715519 0.4777281 	LL:  -382.0379 
#	Beta (IRWLS): -0.7075655 -2.715519 0.4777281 	LL:  -382.0379 
#	Beta (IRWLS): -0.7075655 -2.715519 0.4777281 	LL:  -382.0379 
#	Beta (IRWLS): -0.7075655 -2.715519 0.4777281 	LL:  -382.0379 
#	Beta (IRWLS): -0.7075655 -2.715519 0.4777281 	LL:  -382.0379 
#	
#																Pablo Cingolani
#-------------------------------------------------------------------------------

library(epicalc)

calcGlm <- FALSE
calcGlm <- TRUE

calcIRWLS <- FALSE
calcIRWLS <- TRUE

createSamples <- FALSE
fileName <- "logReg_test_IRWLS_01.txt"

#-------------------------------------------------------------------------------
# Sigmoid finction
#-------------------------------------------------------------------------------

s <- function(h) { 
	1 / (1 + exp(-h)); 
}

sX <- function(beta, X) { 
	eta <- X %*% beta
	mu <- 1 / (1 + exp(- eta));

	n <- length(mu)
	r <- runif(n)
	out <- rep(0, n)
	out[ r < mu ] <- 1
	return (data.frame( out, mu, eta ))
}

log.lik <- function(y, p) {
	sum( log( p^y * (1 - p)^(1 - y) ) );
}

log.lik.X <- function(y, beta, X) {
	s.lik <- sX(beta, X)
	p <- s.lik$mu
	return( log.lik(y, p) )
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

#---
# Real parameters (used to generate data)
#---
beta.0 <- c( -0.75, -3.0, 0.5 )	# Real model

#---
# Create samples using model
#---
N <- 1000	# Numner of samples

# Dimensions
D.0 <- length(beta.0)
D <- length(beta.0) - 1

#---
# Initialize data
#---
if( createSamples ) {
	# Create samples
	X <- matrix( rnorm( N * D ), nrow=N, ncol=D)
	colnames(X) <- c('x1', 'x2')
	X.0 <- cbind( rep(1, N), X)
	colnames(X.0) <- c('x0', 'x1', 'x2')
	s.0 <- sX(beta.0, X.0)
	y <- as.numeric(s.0$out)

	# Save data file
	file <- "logReg_test.txt"
	cat('Saving samples', file, '\n')
	dsave <- data.frame( y, X.0 )
	write.table(dsave, file=file, quote=FALSE, sep="\t", row.names=FALSE)

} else {
	# Load samples
	cat('Loading samples', fileName, '\n')
	d <- read.csv(file=fileName, sep="\t")
	y <- as.numeric( d$y )

	cols <- dim(d)[2]
	X.0 <- as.matrix( d[,2:cols] )
	X <- X.0[,2:(cols-1)]
}

# Calculate output
s.0 <- sX(beta.0, X.0)
cat('Beta (Real ):', beta.0, '\tLL: ', log.lik.X(y, beta.0, X.0) , '\n')

#---
# Calculate logistic regression model (using R's GLM)
#---
if( calcGlm ) {
	# Full model, takes into account genotypes and PCs
	d <- data.frame( out=y, x1=X[,1], x2=X[,2] )
	lr  <- glm( out ~ x1 + x2 , family=binomial, data=d)
	beta.r <- as.numeric( lr$coefficients )
	cat('Beta (R    ):', beta.r, '\tLL: ', log.lik.X(y, beta.r, X.0) ,'\n')
}

#---
# Calculate using IRWLS algorithm
#---
if( calcIRWLS ) {
	# Initialize
	beta.n <- beta.0 * 0
	cat('Beta (IRWLS):', beta.n, '\tLL: ', log.lik.X(y, beta.n, X.0) ,'\n')

	for( i in 1:10 ) {
		# Update parameters
		s.n <- sX(beta.n, X.0)
		eta <- s.n$eta
		mu <- s.n$mu

		nu <- mu * (1 - mu)
		zeta <- eta + (y - mu) / nu
		w <- nu

		cat('\teta  :', eta, '\n')
		cat('\tmu   :', mu, '\n')
		cat('\tw    :', w, '\n')
		cat('\tzeta :', zeta, '\n')

		# Fit weighted model
		dd <- data.frame( X, zeta )
		slm <- lm( zeta ~ x1 + x2 , data=dd, weights=w)

		# Update beta
		beta.n <- as.vector( slm$coefficients )
		cat('Beta (IRWLS):', beta.n, '\tLL: ', log.lik.X(y, beta.n, X.0) ,'\n')
	}
}

