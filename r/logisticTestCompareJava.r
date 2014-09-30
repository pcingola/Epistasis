
#-------------------------------------------------------------------------------
#
# Logistic regression and likelihood ratio test
# Used for debugging data from LogisticeRegression class (Java code)
#
#																Pablo Cingolani
#-------------------------------------------------------------------------------

library(epicalc)

# Load data from TXT file. See Regression.toStringSample() method
d <- read.csv("lr_test.alt.txt", header = TRUE, sep="\t")

#---
# Calculate logistic regression models
#---

# Full model, takes into account genotypes and PCs
lr.alt  <- glm( out ~ in0 + in1 + in2 + in3 + in4 + in5 + in6 + in7 + in8 + in9 + in10 + in11 + in12 , family=binomial, data=d)

# Reduced model: only covariates, no genotypes
lr.null <- glm( out ~       in1 + in2 + in3 + in4 + in5 + in6 + in7 + in8 + in9 + in10 + in11 + in12 , family=binomial, data=d)

# Likelyhood ratio test
lik.ratio.test <- lrtest(lr.alt, lr.null)
pvalue <- lik.ratio.test$p.value   # p-value from likelihood ration test

cat('Iteration', count, '\n')

