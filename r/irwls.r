
#-------------------------------------------------------------------------------
#
# Logistic regression and likelihood ratio test
# Used for debugging data from LogisticeRegression class (Java code)
#
#																Pablo Cingolani
#-------------------------------------------------------------------------------

library(epicalc)

calcGlm <- FALSE
calcGlm <- TRUE

calcIRWLS <- FALSE
calcIRWLS <- TRUE

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

N <- 1000	# Numner of samples
beta.0 <- c( -0.75, -3.0, 0.5 )	# Real model

# Dimensions
D.0 <- length(beta.0)
D <- length(beta.0) - 1

# Initialize data
X <- matrix( rnorm( N * D ), nrow=N, ncol=D)
colnames(X) <- c('x1', 'x2')
X.0 <- cbind( rep(1, N), X)
colnames(X.0) <- c('x0', 'x1', 'x2')

# Calculate output
s.0 <- sX(beta.0, X.0)
y <- as.numeric(s.0$out)
cat('Beta (Real ):', beta.0, '\tLL: ', log.lik.X(y, beta.0, X.0) , '\n')

#---
# Calculate logistic regression models
#---
if( calcGlm ) {
	# Full model, takes into account genotypes and PCs
	d <- data.frame( out=s.0$out, x1=X[,1], x2=X[,2] )
	lr  <- glm( out ~ x1 + x2 , family=binomial, data=d)
	beta.r <- as.numeric( lr$coefficients )
	cat('Beta (R    ):', beta.r, '\tLL: ', log.lik.X(y, beta.r, X.0) ,'\n')
}

#---
# IRWLS
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

		# Fit weighted model
		dd <- data.frame( X, zeta )
		slm <- lm( zeta ~ x1 + x2 , data=dd, weights=w)

		# Update beta
		beta.n <- as.vector( slm$coefficients )
		cat('Beta (IRWLS):', beta.n, '\tLL: ', log.lik.X(y, beta.n, X.0) ,'\n')
	}
}

