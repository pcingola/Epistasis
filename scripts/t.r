
tr <- read.table("t.txt", header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
tr <- as.matrix(tr)

trn <- read.table("tnull.txt", header = TRUE, row.names = 1, sep="\t", na.strings = 'null')
trn <- as.matrix(tr)

# heatmap(tr, Rowv=NA, Colv=NA, col = cm.colors(256), scale="none", margins=c(5,10), main="Transitions AA in contact")
