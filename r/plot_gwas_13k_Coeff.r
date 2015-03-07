
savePlot <- T

#-------------------------------------------------------------------------------
# Show density
#-------------------------------------------------------------------------------
showDens <- function(y, i, type) {
    title <- paste(type, 'distribution: theta', i)
    subtitle <- paste('Mean:', mean(y), " StdDev:", sd(y))

    cat(title, subtitle, '\n')
    plot( density(y), main=title, sub=subtitle)
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

if( savePlot ) { png( width=1024, height=1024 ) }

#---
# Load data
#---
if( ! exists('alt') ) {
    cat('Loading ALT data\n')
    alt <- read.table('gwas.logBf_non_zero.theta.alt.txt', sep='\t', header=F)
}

if( ! exists('null') ) {
    cat('Loading NULL data\n')
    null <- read.table('gwas.logBf_non_zero.theta.null.txt', sep='\t', header=F)
}

# Show all densities
par( mfrow=c(2,3) )
for( i in 1:(dim(alt)[2]) ) {
    showDens( alt[,i], i, 'ALT')
}

par( mfrow=c(2,3) )
for( i in 1:(dim(null)[2]) ) {
    showDens( null[,i], i, 'NULL')
}

# Show parameter correlations
par( mfrow=c(2,2) )
smoothScatter(alt[,1], alt[,2], main="Scatter plot theta1 vs theta2 (ALT)", xlab='theta1', ylab='theta2')
smoothScatter(null[,1], null[,2], main="Scatter plot theta1 vs theta2 (NULL)", xlab='theta1', ylab='theta2')

smoothScatter(alt[,1], alt[,3], main="Scatter plot theta1 vs theta3 (ALT)", xlab='theta1', ylab='theta3')
smoothScatter(alt[,2], alt[,3], main="Scatter plot theta2 vs theta3 (ALT)", xlab='theta2', ylab='theta3')

# Show correlation between ALT and NULL parameters
# Note: We skip paremter number 3 in ALT model (it doesn't exist in NULL model)
for( i in 1:2 ) {
	cat('Theta[', i,'] null vs alt\n')
	smoothScatter(null[,i], alt[,i], main=paste("Scatter plot theta_null vs theta_alt",i), xlab='theta_null', ylab='theta_alt')
	abline(0, 1, col='green')
}
for( i in 3:(dim(null)[2]) ) {
	cat('Theta[', i,'] null vs alt\n')
	smoothScatter(null[,i], alt[,i+1], main=paste("Scatter plot theta_null vs theta_alt",i), xlab='theta_null', ylab='theta_alt')
	abline(0, 1, col='green')
}

if( savePlot )  { dev.off() }
