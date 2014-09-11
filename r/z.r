
savePlot <- T

#-------------------------------------------------------------------------------
# Histogram & density
#-------------------------------------------------------------------------------

histDens <- function( x, title, xlim, q=1.0, breaks = 50 ) {
    # Show only this part of the data
    xmin <- min( xlim )
    xmax <- max( xlim )
    data <- x[ (x >= xmin) & (x <= xmax) ];

    dens <- density(data)

    h <- hist(data, main=title, xlim=xlim, xlab = "data", ylab = "Frequency", freq = T, breaks=breaks);

    # Adjust density height to 'frecuency'
    dens$y <- max(h$counts) * dens$y/max(dens$y)
    lines(dens, col='red', xlim=xlim)

    # Mean & median calculated over the whola data
    abline( v=mean(x), col='blue', lty=2, lwd=2);
    abline( v=median(x), col='green', lty=2, lwd=2);

    legend("topright",c("Mean","Median"),lty=c(1,1),col=c("blue","green"))

}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

# Background distribution MI
if( ! exists('bg.mi') ) {
	cat("Loading background MI\n")
	bg.mi <- read.csv("bg.mi.details.txt")
}

if( ! exists('bg.mi3') ) {
	cat("Loading background MI 3\n")
	bg.mi3 <- read.csv("bg.mi.3.details.txt")
}

# Background distribution Var Inf
if( ! exists('bg.varinf') ) {
	cat("Loading background VarInf\n")
	bg.varinf <- read.csv("bg.varInf.details.txt")
}

if( ! exists('bg.varinf3') ) {
	cat("Loading background VarInf 3\n")
	bg.varinf3 <- read.csv("bg.varInf.3.details.txt")
}

# Values for AA in contact
if( ! exists('aacont') ) {
	cat("Loading contact MI & VarInf\n")
	aacont <- read.csv("aa.contact.stats.txt", sep="\t")
}

if( ! exists('aacont.mi3') ) {
	cat("Loading contact MI 3\n")
	aacont.mi3 <- read.csv("aa.contact.mi.3.vals.txt", sep="\t")
}

if( ! exists('aacont.vi3') ) {
	cat("Loading contact VarInf 3\n")
	aacont.vi3 <- read.csv("aa.contact.varInf.3.vals.txt", sep="\t")
}

if( ! exists('aacount') ) {
	cat("Loading AA count matrix\n")
	aacount <- read.csv("aa.count.matrix.txt", header = TRUE, sep="\t")
}

if( savePlot ) { png( width=800, height=800 ) }

#---
# Densities
#---

par( mfcol=c(2,1) )

cat("MI\n")
bm <- bg.mi[[1]]
bm <- bm[ ! is.na(bm) ]
aami <- aacont$mi
histDens( bm, title="Mutual Information: 'null' distribution (non-zero)", xlim=c(0,2) )
histDens( aami[ aami > 0 ], title="Mutual Information: AA 'in contact' (non-zero)", xlim=c(0,2) )

cat("MI 3\n")
bm <- bg.mi3[[1]]
bm <- bm[ ! is.na(bm) ]
aami <- aacont.mi3[[1]]
histDens( bm, title="Mutual Information: 'null' distribution, 3 AAs (non-zero)", xlim=c(0,2) )
histDens( aami[ aami > 0 ], title="Mutual Information: AA 'in contact', 3 AAs (non-zero)", xlim=c(0,2) )

cat("VarInf\n")
bv <- bg.varinf[[1]]
bv <- bv[ ! is.na(bv) ]
aavi <- aacont$varinf
histDens( bv, title="Variation of Information: 'null' distribution (non-zero)", xlim=c(0,4) )
histDens( aavi[ aavi > 0 ] , title="Variation of Information: AA 'in contact' (non-zero)", xlim=c(0,4) )

cat("VarInf 3\n")
bv <- bg.varinf3[[1]]
bv <- bv[ ! is.na(bv) ]
aavi <- aacont.vi3[[1]]
histDens( bv, title="Variation of Information: 'null' distribution, 3 AAs (non-zero)", xlim=c(0,4) )
histDens( aavi[ aavi > 0 ] , title="Variation of Information: AA 'in contact', 3 AAs (non-zero)", xlim=c(0,4) )

cat("H(X)\n")
h <- c(aacont$hx, aacont$hy)
histDens( h, title="H(X)", xlim=c(0,2) )
histDens( h[ h > 0 ] , title="H(X) non-zero", xlim=c(0,2) )

#---
# Heatmaps
#---

par( mfcol=c(1,1) )

aac <- as.matrix(aacount)
aac[ aac == 0 ] <- NA
heatmap(aac, Rowv=NA, Colv=NA, col = cm.colors(256), scale="row", margins=c(5,10), main='AA count (in contact)')

aami <- read.csv("aa.mi.matrix.txt", header = TRUE, sep="\t")
aam <- as.matrix(aami)
aam[ aam == 0 ] <- NA
heatmap(aam, Rowv=NA, Colv=NA, col = cm.colors(256), scale="row", margins=c(5,10), main='AA MI (in contact)')

aavi <- read.csv("aa.vi.matrix.txt", header = TRUE, sep="\t")
aav <- as.matrix(aavi)
aav[ aav == 0 ] <- NA
heatmap(aav, Rowv=NA, Colv=NA, col = cm.colors(256), scale="row", margins=c(5,10), main='AA VarInf (in contact)')

#---
# Heatmaps by entropy
#---

par( mfcol=c(1,1) )

step <- 0.1
s <- seq(0, 2-step, by=step)
n <- length(s)
count <- matrix(0, n, n)
mi <- count
vi <- count
i <- 1
for( si in s ) {
	j <- 1 
	ki <- (si <= aacont$hx) & (aacont$hx < (si+step))
	for( sj in s ) {
		kj <- (sj <= aacont$hy) & (aacont$hy < (sj+step))
		k <- ki & kj

		count[i,j] <- sum(k)
		
		if( count[i,j] > 0 ) {
			mi[i,j] <- sum( aacont$mi[k] ) / count[i,j]
			vi[i,j] <- sum( aacont$varinf[k] ) / count[i,j]
		} else {
			mi[i,j] <- 0
			vi[i,j] <- 0
		}

		cat(si, "\t", sj, "\t", sum(ki), "\t", sum(kj), "\t", count[i,j], "\t", mi[i,j], "\t", vi[i,j] , "\n")
		j <- j+1
	}
	i <- i+1
}


colnames(count) <- s
rownames(count) <- s
heatmap(log10(count), Rowv=NA, Colv=NA, col = cm.colors(256), scale="none", margins=c(5,10), main="log(Count) by H(X) and H(Y)")

colnames(mi) <- s
rownames(mi) <- s
heatmap(mi, Rowv=NA, Colv=NA, col = cm.colors(256), scale="none", margins=c(5,10), main="MI by H(X) and H(Y)")

colnames(vi) <- s
rownames(vi) <- s
heatmap(vi, Rowv=NA, Colv=NA, col = cm.colors(256), scale="none", margins=c(5,10), main="VarInf by H(X) and H(Y)")

if( savePlot )	{ dev.off() } 
