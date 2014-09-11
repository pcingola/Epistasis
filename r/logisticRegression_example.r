#------------------------------------------------------------------------------
#
# Simulate genotypes and phenotypes, calculate a logistic regression model
#
# Note: This is for two "interacting" genomic loci (epistasis)
#
#															Pablo Cingolani
#------------------------------------------------------------------------------

library(epicalc)

debug <- FALSE				# Debug mode?

n <- 1000					# Samples
m <- 200					# Number of iterations to get an average

beta123 <- c(1, 1, -10)		# Initiali parameters
N <- 2						# Number of parameters


# Model: We avoid having a different result each time we invoke casesControls()
risk <- runif(n)					# Random 'risk' numbers for each person (sample).

var1 <- runif(n)					# THESE SHOULD BE GENOTYPES!
var2 <- runif(n)
var12mult <- var1 * var2			# Used in the model
var12 <- abs( var1 - var2 )			# Real interaction
var12 <- var1 * var2

#------------------------------------------------------------------------------
# Create a distribution of cases and controls
#------------------------------------------------------------------------------

casesControls <- function(beta0, beta123) {
	# Create input
	ones <- rep(1, n)
	X <- cbind( ones, var1, var2, var12 )	# As a matrix

	# Logit(pi)
	beta <- c(beta0, beta123)
	logitp = X %*% beta
	p = 1 / ( 1 + exp(-logitp) )

	# P(yi | X)
	disease <- rep(0, n)
	disease[ risk <= p ] = 1
	numCases = sum(disease)
	cc <- list(var1=var1, var2=var2, var12=var12, X=X, p=p, logitp=logitp, disease=disease, numCases=numCases)
	return(cc)
}

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------

# Create a function that is minimized by finding half cases, half controls
findBeta <- function(x) {
	cc <- casesControls(x, beta123)
	numCases = cc$numCases
	res <- (numCases - n/2)^2
	return(res)
}

# Optimize the function
opt <- optimize( findBeta , c(-30, 10) )
beta0 <- opt$minimum

# Calculate using 'beta' that has equal number of cases and controls
cc <- casesControls(beta0, beta123)
y <- cc$disease
var1 <- cc$var1
var2 <- cc$var2
var12 <- cc$var12

# Show
par(mfrow=c(2,1))
hist( -log10(cc$p) )
hist( cc$p )
cat('Beta0         :', beta0, '\n')
cat('Controls      :', sum( y == 0), '\n')
cat('Cases         :', sum( y ), '\n')
cat('Percent cases :', 100 * sum(y)/length(y), '%\n')

# Logistic regression
lr0 <- glm( y ~ var1 + var2 + var12 , family=binomial) 
lr1 <- glm( y ~ var1 + var2         , family=binomial) 

# Likelihood ratio test
lrt <- lrtest(lr0, lr1)
pvalueLr <- lrt$p.value   # p-value from likelihood ratio test
lrSum <- summary( lr0 )
pvalueWald <- lrSum$coefficients[4,4]
cat('p-value (LR)  :', pvalueLr, '\n')
cat('p-value (Wald):', pvalueWald, '\n')

# Filters (cases, controls)
k0 <- (y == 0)
k1 <- (y == 1)

# Correlations
c <- cor(var1,var2)
c0 <- cor(var1[k0],var2[k0])
c1 <- cor(var1[k1],var2[k1])
cr <- (c1-c0) / c0
cat('\n')
cat('Correlation            :', c , '\n')
cat('Correlation (controls) :', c0, '\n')
cat('Correlation (cases)    :', c1, '\n')
cat('Correlation ratio      :', cr, '\n')
