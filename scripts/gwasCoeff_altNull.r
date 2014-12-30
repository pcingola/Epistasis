
savePng <- T

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
# Show density for a given BF threshold
#-------------------------------------------------------------------------------
showDensBf <- function(y, i, type, bf, bfMin) {
	y <- y[ bf >= bfMin ]
	title <- paste(type, 'distribution: theta', i, 'BF >=', bfMin)
	subtitle <- paste('Mean:', mean(y), " StdDev:", sd(y))

	cat(title, subtitle, '\n')
	plot( density(y), main=title, sub=subtitle)
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

#---
# Load data
#---
if( ! exists('altbf') ) {
	cat('Loading ALT data\n')
	altbf <- read.table('alt.txt', sep='\t', header=F)
}

if( ! exists('nullbf') ) {
	cat('Loading NULL data\n')
	nullbf <- read.table('null.txt', sep='\t', header=F)
}

#---
# Show density plots
#---
alt <- altbf[,2:17]
null <- nullbf[,2:16]
bf <- altbf[,1]

if( savePng )	png( width = 1024 * 1.5, height = 1024 )

# par( mfrow=c(2,3) )
# for( i in 1:(dim(alt)[2]) ) {
# 	showDens( alt[,i], i, 'ALT')
# }
# 
# par( mfrow=c(2,3) )
# for( i in 1:(dim(null)[2]) ) {
# 	showDens( null[,i], i, 'NULL')
# }

#---
# Show Scatter plots
#---
par( mfrow=c(2,3) )
smoothScatter(alt[,1], alt[,2], main="Scatter plot theta1 vs theta2 (ALT)", xlab='theta1', ylab='theta2')
#smoothScatter(null[,1], null[,2], main="Scatter plot theta1 vs theta2 (NULL)", xlab='theta1', ylab='theta2')

smoothScatter(alt[,1], alt[,3], main="Scatter plot theta1 vs theta3 (ALT)", xlab='theta1', ylab='theta3')
smoothScatter(alt[,2], alt[,3], main="Scatter plot theta2 vs theta3 (ALT)", xlab='theta2', ylab='theta3')

for( bfMin in c(2,4,6) ) {
	keep <- ( bf >= bfMin )
	smoothScatter(alt[keep,1], alt[keep,3], main=paste("Scatter plot theta1 vs theta3 (ALT), BF >=", bfMin) , xlab='theta1', ylab='theta3')
}

#---
# Coefficients' distribution for several BF threasholds
#---

# for( i in 1:3 ) {
# 	par( mfrow=c(2,3) )
# 	for( bfMin in 1:6 ) {
# 		showDensBf( alt[,i], i, 'ALT', bf, bfMin)
# 	}
# }

if( savePng )	dev.off()
