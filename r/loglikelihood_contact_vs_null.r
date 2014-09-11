
# How to create the files
#  ( echo llratio ; cut -f 3 likelihood.all.txt ) | grep -v "^$" > likelihood.all.values.txt 

savePlot <- T

if( savePlot ) { png( width=800, height=800 ) }

if( ! exists('lc') ) {
	lc <- read.csv("likelihood.contact.values.txt")
	la <- read.csv("likelihood.null.values.txt")

	lc <- as.numeric( lc[,1] )
	la <- as.numeric( la[,1] )

}

xmin <- -40
xmax <- 25
xlim <- c(xmin, xmax)
breaks <- seq(xmin, xmax, 1)
overlap <- F

# In contact
data <- lc[ (xmin <= lc) & (lc < xmax) ]
dens <- density(data)
hc <- hist(data, main="Log likelyhood ratio histogram", xlim=xlim, xlab = "data", ylab = "Frequency", freq = F, breaks=breaks, col=rgb(1,0,0,1/4) , add=F );
lines(dens, col='red', xlim=xlim)

# All
data <- la[ (xmin <= la) & (la < xmax) ]
dens <- density(data)
ha <- hist(data, xlim=xlim, xlab = "data", ylab = "Frequency", freq = F, breaks=breaks, col=rgb(0,0,1,1/4), add=T );
lines(dens, col='blue', xlim=xlim)

legend("topleft",c("In contact","All"),lty=c(1,1),col=c("red","blue"))

dens.ratio <- hc$density / ha$density
plot( hc$mids, dens.ratio, , main='Log likelyhood density ratio', xlab='Log likelihood', ylab='Density Ratio')
lines( hc$mids, dens.ratio )

if( savePlot )  { dev.off() }

