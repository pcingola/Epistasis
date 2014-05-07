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

#------------------------------------------------------------------------------
# Parametric model: Logistic regression
#------------------------------------------------------------------------------

N <- 2									# Number of parameters in logistic regression
n <- 1000								# Number of samples
af <- 0.2								# Allele frequency
OR <- 1.3 								# Odds ratio	
lor <- log(OR)							# Log of odds ratio
interaction <- 2 						# Interaction of two loci (negative means protective, positive means risk)

gwas.snps <- 100000						# We assume this number of coding, non-synonimous variants in a GWAS of 'n' samples
gwas.p <- 0.05 / (gwas.snps^2)			# p-value required for GWAS level detection
gwas.q <- -log10( gwas.p )				# q-value (i.e. -log(p-value) ) required for GWAS level detection
single.q <- -log10( 0.05 )				# q-value (i.e. -log(p-value) ) required for single test detection

# Logistic regression parameters: beta = ( beta_0, beta_1, beta_2, beta_3 )
# Note: beta_0 is adjusted to produce an equal number of cases and controls (see optimization)
beta123 <- c( lor, lor, interaction * lor)

#------------------------------------------------------------------------------
# Create a distribution of cases and controls
#------------------------------------------------------------------------------

casesControls <- function(beta0, gt) {
	# Create input
	ones <- rep(1, n)
	X <- cbind( ones, gt$var1, gt$var2, gt$var12 )	# As a matrix

	# Logit(pi)
	beta <- c(beta0, gt$beta123)
	logitp = X %*% beta
	p = 1 / ( 1 + exp(-logitp) )

	# P(yi | X)
	disease <- rep(0, n)
	disease[ gt$risk <= p ] = 1
	numCases = sum(disease)
	
	cc <- list(p=p, logitp=logitp, disease=disease, numCases=numCases)
	return( cc )
}

#------------------------------------------------------------------------------
# Create a function that is minimized by finding half cases, half controls
#------------------------------------------------------------------------------

findBeta <- function(beta0, gt) {
	cc <- casesControls(beta0, gt)
	numCases = cc$numCases
	res <- (numCases - n/2)^2
	return(res)
}

#------------------------------------------------------------------------------
# Calculate the p-value
#------------------------------------------------------------------------------

logRegPval <- function(cc, gt) {
	var1 <- gt$var1
	var2 <- gt$var2
	var12 <- gt$var12


	y <- cc$disease

	# Logistic regression
	lr0 <- glm( y ~ var1 + var2 + var12 , family=binomial) 
	lr1 <- glm( y ~ var1 + var2         , family=binomial) 

	# Likelihood ratio test, p-value
	lrt <- lrTestWarn(lr0, lr1)
	pvalueLr <- lrt$p.value   

	# Wald test p-value
	lrSum <- summary( lr0 )
	if( min(dim(lrSum$coefficients)) == 4 ) { 
		pvalueWald <- lrSum$coefficients[4,4]
	} else {
		pvalueWald <- 1
	}

	# Create a list of results
	ret <- list( lr0=lr0, lr0lr0=lr0, lrSum=lrSum, logRegTest=lrt, pvalueLr=pvalueLr, pvalueWald=pvalueWald )
	return( ret )
}

#------------------------------------------------------------------------------
# A logistic regression test
#------------------------------------------------------------------------------

logRegTest <- function(n, af, show=FALSE) {
	# Model: We avoid having a different result each time we invoke casesControls()
	risk <- runif(n)					# Random 'risk' numbers for each person (sample).
	var1 <- randgt(n, af)				# Genotypes for SNP1
	var2 <- randgt(n, af)				# Genotypes for SNP2
	var12 <- var1 * var2				# Interacting genotypes
	gt <- list(beta123=beta123, risk=risk, var1=var1, var2=var2, var12=var12)

	# Optimize the function to find half cases, half controls
	fbeta <- function(x) { return( findBeta(x, gt) ); }
	opt <- optimize( fbeta , c(-30, 10) )
	beta0 <- opt$minimum

	# Calculate using 'beta' that has equal number of cases and controls
	cc <- casesControls(beta0, gt)
	y <- cc$disease

	# Logistic regression and p-values
	lr <- logRegPval(cc, gt)

	# Show
	if( show ) {
		par(mfrow=c(2,1))
		hist( -log10(cc$p) )
		hist( cc$p )
		cat('AF               :', af, '\n')
		cat('Count(var12)     :', sum( gt$var12 > 0 ), '\n')
		cat('Beta0            :', beta0, '\n')
		cat('Controls         :', sum( y == 0), '\n')
		cat('Cases            :', sum( y ), '\n')
		cat('Percent cases    :', 100 * sum(y)/length(y), '%\n')
		cat('p-value (LR)     :', lr$pvalueLr, '\n')
		cat('p-value (Wald)   :', lr$pvalueWald, '\n')
	}

	res <- list( gt=gt, cc=cc, lr=lr )
	return(res)
}

#------------------------------------------------------------------------------
# Redefine lrtest function from 'epicalc' library
#
# We re-define this function in order to change error condition "Likelihood 
# gets worse with more variables. Test not executed". In this implementation 
# is just a warning
#
#------------------------------------------------------------------------------

lrTestWarn <- function (model1, model2) {
    if (any(class(model1) != class(model2))) {
        stop("Two models have different classes")
    }

	if (suppressWarnings((all(class(model1) == c("glm", "lm")) & all(class(model2) == c("glm", "lm"))) | (any(class(model1) == "negbin") & any(class(model2) == "negbin")))) {
		if (sum(model1$df.null) != sum(model2$df.null)) stop("Number of observation not equal!!")

		df1 <- attributes(logLik(model1))$df
		df2 <- attributes(logLik(model2))$df
		lrt <- 2 * (as.numeric(logLik(model2) - logLik(model1)))
		diff.df <- df2 - df1
		
		if (lrt < 0) {
			lrt <- -lrt
			diff.df <- -diff.df
		}

		if (lrt * diff.df < 0) {
			cat("WARNING: Likelihood gets worse with more variables. Test not executed\n")
			lrt <- 1
		}
	} else {
		stop("Unimplemented model! Try using 'lrtest' from epicalc package")
	}

	output <- list(model1 = model1$call, model2 = model2$call, model.class = class(model1), Chisquared = lrt, df = diff.df, p.value = pchisq(lrt, diff.df, lower.tail = FALSE))
	class(output) <- "lrtest"
	output
}

#------------------------------------------------------------------------------
# Odds / oddsRatio
#------------------------------------------------------------------------------

odds <- function(p) { return ( p / (1-p) ); }

oddsRatio <- function(p1,p2) { return ( odds(p1) / odds(p2) ); }


#------------------------------------------------------------------------------
# Random Genotype
#------------------------------------------------------------------------------

randgt <- function(n, af) {
	r <- runif(n)
	gt <- rep(0, n)
	gt[ r <= af ] <- 1
	gt[ r < (af^2) ] <- 2
	return(gt)
}

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------

savePdf <- TRUE

if( FALSE ) {
	lrt <- logRegTest(n, af, TRUE)
	P <- data.frame( y=lrt$cc$disease, p=lrt$cc$p, var1=lrt$gt$var1, var2=lrt$gt$var2, var12=lrt$gt$var12 )
} else if( TRUE ) {

	K <- 1000
	ns <- c(0.1, 0.5, 1, 2.5, 5, 7.5, 10, 15, 20, 25, 30) * K
	#ns <- c(100)

	if( savePdf ) {
		# Save plots to PDF file
		outFile <- paste('log_p.OR_', OR ,'.pdf', sep='')
		pdf(outFile)
	}
	
	for( n in ns ) {

		afs <- 1:200/1000
		iter <- 10
		qs <- rep(0, length(afs) )
		i <- 1
		for( af in afs ) {

			slpv <- 0
			for( j in 1:iter ) {
				lrt <- logRegTest(n, af, FALSE)
				pval <- lrt$lr$pvalueLr
				lpv <- -log10(pval)
				slpv <- slpv + lpv
			}

			avgq <- slpv/iter
			cat('n:', n, '\taf:', af, '\tavg( log(p_value) ): ', avgq ,'\n')
			qs[i] <- avgq
			i <- i + 1
		}

		title <- paste('log10( p_value ) vs AF\nSample size: ', n )
		ors <- format( exp(beta123) , digits=2)
		subtitle <- paste('OR: ', ors[1], '    Interaction OR:', ors[3])
		plot( afs, qs, main=title, sub=subtitle, xlab='Allele frequency', ylab='log10(p_value)' )

		abline(h=single.q, col='red')
		abline(h=gwas.q, col='green')
	}

	if( savePdf )	dev.off()

} else if( FALSE ) {
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
}
