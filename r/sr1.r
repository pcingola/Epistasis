
#-------------------------------------------------------------------------------
#
# SR1 algorithm
#
#																Pablo Cingolani
#-------------------------------------------------------------------------------

library(epicalc)

#-------------------------------------------------------------------------------
# Sigmoid finction
#-------------------------------------------------------------------------------

s <- function(h) { 
	1 / (1 + exp(-h)); 
}

log.lik <- function(y, p) {
	sum( log( p^y * (1 - p)^(1 - y) ) );
}

sr1 <- function(sample.X, sample.y, xk) {
	xk 
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

name <- 'lr_test.1_10002.alt'
name <- 'lr_test.1_10001.alt'

# Load data from TXT file. See Regression.toStringSample() method
dataFile <- paste( name, '.txt', sep='')
modelFile <- paste( name, '.model.txt', sep='')
d <- read.csv(dataFile , header = TRUE , sep="\t")

beta.j <- read.csv(modelFile, header = FALSE, sep="\t")
beta.j <- as.numeric(beta.j)

#---
# Calculate logistic regression models
#---

# Transform input data to matrix and calculate logistic regression

sample.X <- as.matrix( d[,1:14] )
sample.y <- as.numeric(d$out)

size <- dim( x )[2]
beta <- rep(0, size);

sr1( sample.X, sample.y, beta )
