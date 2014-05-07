#------------------------------------------------------------------------------
# Example of a Bayesian factor calculation for logistic model
#------------------------------------------------------------------------------

N <- 300

#------------------------------------------------------------------------------
# Step size in a vector
#------------------------------------------------------------------------------

step <- function(x) { (max(x) - min(x)) / (length(x)-1); }

#------------------------------------------------------------------------------
# Likelyhood function
#------------------------------------------------------------------------------
progLogit <- function(g, mu, beta) {
	h <- mu + beta * g
	p <- 1 / (1 + exp(-h))
}

#------------------------------------------------------------------------------
# Likelyhood function
#------------------------------------------------------------------------------
lik <- function(y, g, mu, beta) {
	p <- progLogit(g, mu, beta)
	lik <- prod( p^y * ((1-p)^(1-y)) )
}

loglik <- function(y, g, mu, beta) {
	p <- progLogit(g, mu, beta)
	lik <- sum( y * log(p) + (1-y) * log(1-p) )
}

#------------------------------------------------------------------------------
# Priors
#------------------------------------------------------------------------------

pi.mu <- function(x)	{ dnorm(x.mu, 0, 0.05); }
pi.beta <- function(x)	{ dnorm(x, beta, 0.5); }

#------------------------------------------------------------------------------
# Intergrate likelyhood function over a prior
#------------------------------------------------------------------------------
intlik <- function(y, g, x.mu, beta) {
	s <- sapply( x.mu, function(x) { lik(y, g, x, beta); } )
	p.mu <- pi.mu( x.mu )
	st <- step(x.mu)
	return( st * sum(s * p.mu) )
}

#------------------------------------------------------------------------------
# Intergrate likelyhood function over a prior
#------------------------------------------------------------------------------
intintlik <- function(y, g, x.mu, x.beta) {
	s <- sapply( x.beta, function(x) { intlik(y, g, x.mu, x); } )
	p.beta <- pi.beta( x.beta )
	st <- step(x.beta)
	return( st * sum(s * p.beta) )
}

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------

#---
# Create data based on model
#---

r <- runif(N)
g <- 2 * runif(N) - 1 	# Random 'genotypes'. In real data this should be {0,1,2}, but we use uniform random [-1, 1]

# Model:  logit(p) = mu + beta * g
beta <- 1.0
mu <- -mean( beta * g)

# Create output 'y' based on model
p <- progLogit(g, mu, beta)
y <- 1:N * 0
y[ r <= p ] <- 1

#---
# Bayesian factor
#---

par(mfrow=c(2,1))

# Prior mu
x.mu <- (-100:100)/500
p.mu <- pi.mu( x.mu )
plot( x.mu, p.mu, type='l', main='Prior( mu )')

# Prior for beta
x.beta <- (-100:100)/50 + beta
p.beta <- pi.beta( x.beta ) 
plot( x.beta, p.beta, type='l', main='Prior( beta )')

# Calculate likelyhood
lik.M0 <- intlik(y, g, x.mu, 0); 
lik.M1 <- intintlik(y, g, x.mu, x.beta); 

# Bayes factor
bf <- lik.M1 / lik.M0
cat('Likelyhood M0   :', lik.M0, '\n')
cat('Likelyhood M1   :', lik.M1, '\n')
cat('Bayesian factor :', bf, '\n')
cat('log10(BF)       :', log10(bf), '\n')

#---
# Laplace approximation
#---

log.h <- loglik(y, g, mu, beta) + log( pi.mu(mu) ) + log( pi.beta(beta) ) + 2 / 2 * log(2 * pi) - 1 / 2 log()

