
#-------------------------------------------------------------------------------
#
# Logistic regression and likelihood ratio test
# Used for debugging data from LogisticeRegression class (Java code)
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

