
savePlot <- T

if( savePlot ) { png( width=800, height=800 ) }

lc <- read.csv("likelihood.contact.values.txt")
la <- read.csv("likelihood.all.values.txt")

lc <- as.numeric( lc[,1] )
la <- as.numeric( la[,1] )

xlim <- c(-400, 50)
breaks <- 50

# In contact
data <- lc
dens <- density(data)
h <- hist(data, main="Log likelyhood ratio histogram", xlim=xlim, xlab = "data", ylab = "Frequency", freq = T, breaks=breaks, col=rgb(1,0,0,1/4) );
dens$y <- max(h$counts) * dens$y/max(dens$y)
lines(dens, col='red', xlim=xlim)

# All
data <- la
dens <- density(data)
h <- hist(data, xlim=xlim, xlab = "data", ylab = "Frequency", freq = T, breaks=breaks, add = T, col=rgb(0,0,1,1/4) );
dens$y <- max(h$counts) * dens$y/max(dens$y)
lines(dens, col='blue', xlim=xlim)

legend("topleft",c("In contact","All"),lty=c(1,1),col=c("red","blue"))

if( savePlot )  { dev.off() }

