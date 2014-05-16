
library('Matrix')

#-------------------------------------------------------------------------------
# Single digit random numbers
#-------------------------------------------------------------------------------
rrand <- function( n ) { floor( runif(n) * 100 ) / 10; }

#-------------------------------------------------------------------------------
# Create a random symmetric matrix 
#-------------------------------------------------------------------------------
randMatrix <- function( n ) {
	m <-  matrix( rrand(n^2), n,n)
	m <- ( m + t(m) ) / 2
	e <- eigen( m ) 
	V <- e$vectors
	lambda <- e$values
	lambda <- lambda - max(lambda)
	D <- diag( lambda )
	m <- V %*% D %*% t(V)
	return(m)
}

#-------------------------------------------------------------------------------
# Random sparse matrix of 'minCount' non-zero elements (symmetric)
#-------------------------------------------------------------------------------
randMatrixSparse <- function( n, minCount ) {
	# Sparse random matris
	count <- 0
	Q <- randMatrix(n)
	Z <- matrix( runif(n^2), n,n)
	q <- quantile( Z, minCount / n^2 )
	zero <- (Z > q) 
	Q[ zero ] <- 0
	diag( Q ) <- 0

	return( Q + t(Q) )
}

#-------------------------------------------------------------------------------
# Approximate exp(A+B)
#-------------------------------------------------------------------------------
aproxExpAB <- function (eigA, B, m) {
	n <- 2^m

	# Compute X <- expm(A)
	V <- eigA$vectors
	lambda <- eigA$values
	X <- V %*% diag(exp(lambda/n)) %*% t(V)

	# Compute Y <- expm(B)
	Y <- expm(B/n)

	# Compute (X*Y)^(2^m)
	XY <- X %*% Y
	for( i in 1:m ) {
		XY <- XY %*% XY
	}

	return( XY )
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

# Create block matrix
if( ! exists('Q') ) { 
	n <- 3
	n2 <- 3 * n

	Q <- randMatrix(n)
	Qh <- as.matrix( bdiag(Q, Q, Q) )

	# Sparse random matris
	count <- 0
	while( count <= 2 ) {
		Qr <- matrix( rrand( n2^2 ), n2, n2)
		zero <- (Qr < 9.5) | ( abs(Q2) > 0 )
		count <- sum( !zero )
		cat('Count non-zero:', count, '\n')
	}

	Qr <- matrix( rrand( (n2)^2 ), n2, n2 )
	Qr[zero] <- 0
	Qr <- ( t(Qr) + Qr ) / 2 
}

# l <- 1/2
# Q2 <- (1-l) * Qh + l * Qr
# 
# # Diagonalize
# e <- eigen(Q2)
# lambda <- e$values
# V <- e$vectors
# D <- diag( lambda )
# cat('Error: ', (sum(V %*% D %*% t(V) - Q2)), '\n')

n <- 400
A <- randMatrix(n)
eigA <- eigen(A)
B <- randMatrixSparse(n, 5)
eab <- expm(A+B)

for( m in 5:20 ) {
	t <- system.time( apeab <- aproxExpAB(eigA, B, m) )
	cat('Number of iterations:', m, '\terror:', sum(abs(eab-apeab)), '\t', t, '\n')
}

