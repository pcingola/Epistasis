
#-------------------------------------------------------------------------------
#
# Felsenstein algorithm example:
#
# Reference: "Computational Molecular Evolution" Z. Yang, Section 4.2.2, page 104
#
#
#															Pablo Cingolani 2014
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

# Now in one equation
