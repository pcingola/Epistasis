
#-------------------------------------------------------------------------------
#
# Felsenstein algorithm example:
#
# Reference: "Computational Molecular Evolution" Z. Yang, Section 4.2.2, page 104
#
#
#															Pablo Cingolani 2014
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Matrix square root
#-------------------------------------------------------------------------------

matrixSqrt <- function( A ) { 
	e <- eigen(A)
	V <- e$vectors
	return ( V %*% diag( sqrt(e$values) ) %*% t(V) )
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

t1 <- t2 <- t3 <- t4 <- t5 <- 0.2
t6 <- t7 <- t8 <- 0.1

#---
# Create matrices
#---
v1 <- 0.906563
v2 <- 0.045855
v3 <- 0.023791
Pt1 <- matrix( c(v1,v2,v3,v3,v2,v1,v3,v3,v3,v3,v1,v2,v3,v3,v2,v1) , nrow=4, ncol=4 )

v1 <- 0.825092
v2 <- 0.084274
v3 <- 0.045317
Pt2 <- matrix( c(v1,v2,v3,v3,v2,v1,v3,v3,v3,v3,v1,v2,v3,v3,v2,v1) , nrow=4, ncol=4 )

P6 <- P7 <- P8 <- Pt1
P1 <- P2 <- P3 <- P4 <- P5 <- Pt2

#---
# Leaf nodes
#---
L1 <- c(1,0,0,0)	# T
L2 <- c(0,1,0,0)	# C
L3 <- c(0,0,1,0)	# A
L4 <- c(0,1,0,0)	# C
L5 <- c(0,1,0,0)	# C

#---
# Calculate likelihoods
#---

L8 <- (P4 %*% L4) * (P5 %*% L5)
L7 <- (P1 %*% L1) * (P2 %*% L2)
L6 <- (P7 %*% L7) * (P3 %*% L3)
L0 <- (P6 %*% L6) * (P8 %*% L8)

#---
# Calculate likelihood using one equation
#---

LL0 <- (P6 %*% ((P7 %*% ((P1 %*% L1) * (P2 %*% L2))) * (P3 %*% L3))) * (P8 %*% ((P4 %*% L4) * (P5 %*% L5)))

#---
# Expand equation Z0
#---

Qt1 <- matrixSqrt(Pt1)
Qt2 <- matrixSqrt(Pt2)

# Sanity check
check1 <- norm( Pt1 - (Qt1 %*% Qt1) ) 
check2 <- norm( Pt2 - (Qt2 %*% Qt2) )
if( (check1 + check2) > 0.001 )	stop("Matrix square root did not work!")

Q6 <- Q7 <- Q8 <- Qt1
Q1 <- Q2 <- Q3 <- Q4 <- Q5 <- Qt2

#Z0 <- (P6 %*% ((P7 %*% ((P1 %*% L1) * (P2 %*% L2))) * (P3 %*% L3))) * (P8 %*% ((P4 %*% L4) * (P5 %*% L5)))
#Z0 <- (P6 %*% ((P7 %*% ((P1 %*% L1) * (P2 %*% L2))) * (P3 %*% L3))) * ((Q8 %*% Q8) %*% ((P4 %*% L4) * (P5 %*% L5)))
#Z0 <- (P6 %*% ((P7 %*% ((P1 %*% L1) * (P2 %*% L2))) * (P3 %*% L3))) * ((Q8 %*% P4 %*% L4) * (Q8 %*% P5 %*% L5))
Z6 <- (P7 %*% ((P1 %*% L1) * (P2 %*% L2))) * (P3 %*% L3)


if( norm(L6 - Z6) > 0.001 ) {
	cat('ERROR!!!!\nDifference: \n')
	print( L6 - Z6 )
} else {
	cat("OK:", norm(L6 - Z6), '\n')
}
