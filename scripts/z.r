
library(ggplot2)

# Background distribution MI
if( ! exists('bgmi') ) {
	bgmi <- read.csv("bg.mi.details.txt")
}

# Background distribution Var Inf
if( ! exists('bgvi') ) {
	bgvi <- read.csv("bg.varInf.details.txt")
}

# Values for AA in contact
if( ! exists('aacont') ) {
	aacont <- read.csv("aa.contact.mi.values.txt", sep="\t")
}

#---
# Densities
#---

#png(width = 1000, height = 1200)
#par( mfcol=c(2,1) )

#bm <- bgmi[[1]]
#bm <- bm[ ! is.na(bm) ]
#plot( density(bm), main="Mutual Information: null distribution", xlim=c(0,2) )
#plot( density( aacont$mi), col='green', main="Mutual Information: AA in contact", xlim=c(0,2) )
#
#bm <- bgmi[[1]]
#bm <- bm[ ! is.na(bm) ]
#plot( density(bm), main="Mutual Information: null distribution", xlim=c(0,2) )
#plot( density( aacont$mi[ aacont$mi > 0 ]), col='red', main="Mutual Information: AA in contact (non-zero)", xlim=c(0,2) )
#
#bv <- bgvi[[1]]
#bv <- bv[ ! is.na(bv) ]
#plot( density(bv), main="Variation of Information: null distribution", xlim=c(0,4) )
#plot( density( aacont$varinf ), col='green', main="Variation of Information: AA in contact", xlim=c(0,4) )
#
#bv <- bgvi[[1]]
#bv <- bv[ ! is.na(bv) ]
#plot( density(bv), main="Variation of Information: null distribution", xlim=c(0,4) )
#plot( density( aacont$varinf[ aacont$varinf > 0 ] ), col='red', main="Variation of Information: AA in contact (non-zero)", xlim=c(0,4) )
#
#dev.off()

#---
# Heatmaps
#---

png( width=800, height=800 )
par( mfcol=c(1,1) )

aacount <- read.csv("aa.count.matrix.txt", header = TRUE, sep="\t")
aac <- as.matrix(aacount)
aac[ aac == 0 ] <- NA
heatmap(aac, Rowv=NA, Colv=NA, col = cm.colors(256), scale="row", margins=c(5,10))

aami <- read.csv("aa.mi.matrix.txt", header = TRUE, sep="\t")
aam <- as.matrix(aami)
aam[ aam == 0 ] <- NA
heatmap(aam, Rowv=NA, Colv=NA, col = cm.colors(256), scale="row", margins=c(5,10))

aavi <- read.csv("aa.vi.matrix.txt", header = TRUE, sep="\t")
aav <- as.matrix(aavi)
aav[ aav == 0 ] <- NA
heatmap(aav, Rowv=NA, Colv=NA, col = cm.colors(256), scale="row", margins=c(5,10))

dev.off()
