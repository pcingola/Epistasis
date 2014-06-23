
#library(ggplot2)

savePlot = TRUE

# Background distribution MI
if( ! exists('bg.mi3') ) {
	cat("Loading background MI\n")
	#bg.mi <- read.csv("bg.mi.details.txt")
	bg.mi3 <- read.csv("bg.mi.3.details.txt")
}

# Background distribution Var Inf
if( ! exists('bg.varinf3') ) {
	cat("Loading background VarInf\n")
	#bg.varinf <- read.csv("bg.varInf.details.txt")
	bg.varinf3 <- read.csv("bg.varInf.3.details.txt")
}

# Values for AA in contact
if( ! exists('aacont.mi3') ) {
	cat("Loading contact MI & VarInf\n")
	aacont.mi3 <- read.csv("aa.contact.mi.3.vals.txt", sep="\t")
	aacont.vi3 <- read.csv("aa.contact.varInf.3.vals.txt", sep="\t")
}

#---
# Densities
#---

if( savePlot ) { png( width=800, height=800 ) }
par( mfcol=c(2,1) )

cat("OK 1\n")
bm <- bg.mi3[[1]]
bm <- bm[ ! is.na(bm) ]
aami <- aacont.mi3[[1]]
plot( density(bm), main="Mutual Information: null distribution", xlim=c(0,2) )
plot( density( aami ), col='green', main="Mutual Information: AA in contact", xlim=c(0,2) )

cat("OK 2\n")
plot( density(bm), main="Mutual Information: null distribution", xlim=c(0,2) )
plot( density( aami[ aami > 0 ] ), col='red', main="Mutual Information: AA in contact (non-zero)", xlim=c(0,2) )

cat("OK 3\n")
bv <- bg.varinf3[[1]]
bv <- bv[ ! is.na(bv) ]
aavi <- aacont.vi3[[1]]
plot( density(bv), main="Variation of Information: null distribution", xlim=c(0,4) )
plot( density( aavi ), col='green', main="Variation of Information: AA in contact", xlim=c(0,4) )

cat("OK 4\n")
plot( density(bv), main="Variation of Information: null distribution", xlim=c(0,4) )
plot( density( aavi[ aavi > 0 ] ), col='red', main="Variation of Information: AA in contact (non-zero)", xlim=c(0,4) )

#---
# Heatmaps
#---

# par( mfcol=c(1,1) )
# 
# aacount <- read.csv("aa.count.matrix.txt", header = TRUE, sep="\t")
# aac <- as.matrix(aacount)
# aac[ aac == 0 ] <- NA
# heatmap(aac, Rowv=NA, Colv=NA, col = cm.colors(256), scale="row", margins=c(5,10))
# 
# aami <- read.csv("aa.mi.matrix.txt", header = TRUE, sep="\t")
# aam <- as.matrix(aami)
# aam[ aam == 0 ] <- NA
# heatmap(aam, Rowv=NA, Colv=NA, col = cm.colors(256), scale="row", margins=c(5,10))
# 
# aavi <- read.csv("aa.vi.matrix.txt", header = TRUE, sep="\t")
# aav <- as.matrix(aavi)
# aav[ aav == 0 ] <- NA
# heatmap(aav, Rowv=NA, Colv=NA, col = cm.colors(256), scale="row", margins=c(5,10))

if( savePlot )	{ dev.off() } 
