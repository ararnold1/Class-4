
source("https://bioconductor.org/biocLite.R")

biocLite("affy")

biocLite("SpikeInSubset")

library(SpikeInSubset)

data(spikein95)

image(spikein95)

ids <- geneNames(spikein95)

ids[1:10]

mas5.eset <- mas5(spikein95)

mas5.e <- log2(exprs(mas5.eset))

boxplot(spikein95)

x11()

boxplot(mas5.e, col = 2:5)

density1 <- density(mas5.e[, 1])

plot(density1, main = "MAS5 expression measure distributions")

density2 <- density(mas5.e[, 2])

lines(density2, col = "red")

density3 <- density(mas5.e[, 3])

lines(density3, col = "blue")

# Making MA plots
# M:difference in average log intensities
#A: Average log intensities
#split pictures exp vs control

d <- rowMeans(mas5.e[,1:3]) - rowMeans(mas5.e[,4:6])
a <- rowMeans(mas5.e)

#plotting the data
plot(a,d, ylim = c(-5,5), main="MAS 5.0 MA plot", xlab = "A", ylab = "M", pch = ".")
abline(h=c(-1, 1))

# name() is a function 
# ylab and xlab used to label axis
# pch is "plot as"
# fold difference is above or below 1.5 between exp and control sample considerened up or down regulated

abline(h=c(-1.5, 1.5), col="red")

#Rhobust multi array chip normalization (RMA)- between chip and within chip 

rma.eset <- rma(spikein95)
rma.e <- exprs(rma.eset)

#compare the different normalization techniques by using head(rma.e) and head(mas5.e)

x11()
boxplot (mas5.e, col=2:5, main="mas 5.0 Norm")
x11()
boxplot (rma.e, col=2:5, main="RMA Norm")

#Only except genes that are upregulated or downregulated according to 3 normalization methods else false discovery

#finding specific genes or probe sets in MA plot

spikedn <- colnames(pData(spikein95)) 
spikedIndex <- match(spikedn, featureNames (mas5.eset))
points(a[spikedIndex], d[spikedIndex], pch=19, col="red")


