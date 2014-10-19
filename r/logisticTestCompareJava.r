
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

#-------------------------------------------------------------------------------
# Sigmoid finction
#-------------------------------------------------------------------------------
s <- function(h) { 
	1 / (1 + exp(-h)); 
}

log.lik <- function(y, p) {
	sum( log( p^y * (1 - p)^(1 - y) ) );
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

name <- 'lr_test.1_10001.alt'
name <- 'lr_test.1_10002.alt'

# lr_test.1_10002.alt:
#	(Intercept)   in0   in1        in2        in3     in4    in5      in6      in7         in8         in9         in10       in11     in12  
#	-2.657e+01    5.3   -5.8e-09   -4.3e-09   -2e-10  -9e-09 -5.5e-09 -5.1e-09 -2.887e-09  6.914e-09   -1.7e-09    1.2e-09    1.2e-11  1.4e-10  
#
#                 [5.3, 0,         0,         0,      0,     0,       0,       0,          0,          0,          0,         0,       0,       -0.26 ]
#

# Load data from TXT file. See Regression.toStringSample() method
if( ! exists('d') ) {
	dataFile <- paste( name, '.txt', sep='')
	modelFile <- paste( name, '.model.txt', sep='')
	d <- read.csv(dataFile , header = TRUE , sep="\t")

	beta.j <- read.csv(modelFile, header = FALSE, sep="\t")
	beta.j <- as.numeric(beta.j)
}

#---
# Calculate logistic regression models
#---
if( calcGlm ) {
	# Full model, takes into account genotypes and PCs
	lr.alt  <- glm( out ~ in0 + in1 + in2 + in3 + in4 + in5 + in6 + in7 + in8 + in9 + in10 + in11 + in12 , family=binomial, data=d)

	# Reduced model: only covariates, no genotypes
	lr.null <- glm( out ~       in1 + in2 + in3 + in4 + in5 + in6 + in7 + in8 + in9 + in10 + in11 + in12 , family=binomial, data=d)

	# Likelyhood ratio test
	lik.ratio.test <- lrtest(lr.alt, lr.null)
	pvalue <- lik.ratio.test$p.value   # p-value from likelihood ration test

	#---
	# Alt model (Java and R)
	#---
	beta.r <- as.numeric( lr.alt$coefficients )
	beta.r[ is.na(beta.r) ] <- 0	# Remove 'NA'
	beta.r <- c( beta.r[-1] , beta.r[1] )	# Re-order elements: same order as beta.j (intercept is the last element instead of the first one)

	# Transform input data to matrix and calculate logistic regression
	x <- as.matrix( d[,1:14] )
	y <- as.numeric(d$out)

	h.j <- x %*% beta.j
	p.j <- s( h.j )
	ll.alt.j <- log.lik( y, p.j )
	cat('Log likelihood Alt  (Java):', ll.alt.j, '\n')

	h.r <- x %*% beta.r
	p.r <- s( h.r )
	ll.alt.r <- log.lik( y, p.r )
	cat('Log likelihood Alt  (R   ):', ll.alt.r, '\n')

	if( ll.alt.j > ll.alt.r ) cat('Java likelihood is better :', (ll.alt.j - ll.alt.r), '\n')
}

#---
# IRWLS
#---

# Initialize
X <- as.matrix( d[,1:13] )		# Data to fit
Xzero <- as.matrix( d[,1:14] )	# Data to fi: Include intercepts (last columns of '1')

N <- dim(X)[2] + 1
beta <- rep(0, N)

y <- as.numeric(d$out)

for( i in 1:10 ) {
	# Update parameters
	eta <- Xzero %*% beta
	mu <- 1 / (1 + exp( - eta))
	nu <- mu * (1 - mu)
	zeta <- eta + (y - mu) / nu
	w <- nu

	# Fit weighted model
	dd <- data.frame( X, zeta )
	slm <- lm( zeta ~ in0 + in1 + in2 + in3 + in4 + in5 + in6 + in7 + in8 + in9 + in10 + in11 + in12 , data=dd, weights=w)

	# Update beta
	beta <- as.vector( c( slm$coefficients[2:14], slm$coefficients[1] ) )
	cat('Beta:', beta, '\n')
}

