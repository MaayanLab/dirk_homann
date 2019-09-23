celpath = "/Users/verena_vdh/Desktop/Verena Van Der Heide/CEL files"
list = list.files(celpath,full.names=TRUE)
data = read.celfiles(list)
library("pd.mta.1.0")

sampleNames=rowNames
sample=index
ph$sample = ph$index

# microarray pictures

for (i in 1:8)
{
  name = paste("image",i,".jpg",sep="")
  jpeg(name)
  image(data[,i],main=ph@data$index[i])
  dev.off()
}

#weight-based pseudo chip image not available for R v3.4.4

for (i in 1:8)
{
  name = paste("pseudoimage",i,".jpg",sep="")
  jpeg(name)
  image(Pset,which=i,main=ph@data$index[i])
  dev.off()
}

#histograms

color=c('green','green','green','green','red','red','red','red')
hist(data[,1:8],lwd=2,which='pm',col=color,ylab='Density',xlab='Log2 intensities',main='Histogram of raw data')

#histogram ggplot all arrays

pmexp = pm(data)

rowNames = vector()
logs = vector()
for (i in 1:8)
{
  rowNames = c(rowNames,rep(ph@data[i,1],dim(pmexp)[1]))
  logs = c(logs,log2(pmexp[,i]))
}

logData = data.frame(logInt=logs,index=rowNames)

dataHist2 = ggplot(logData, aes(logInt, colour = indexName)) 
dataHist2 + geom_density()

#ggplot boxplot not normalized

pmexp = pm(data)

rowNames = vector()
logs = vector()
for (i in 1:8)
{
  rowNames = c(rowNames,rep(ph@data[i,1],dim(pmexp)[1]))
  logs = c(logs,log2(pmexp[,i]))
}

logData = data.frame(logInt=logs,index=rowNames)

dataBox = ggplot(logData,aes(sampleName,logInt))
dataBox + geom_boxplot()

#RMA normalization

data.rma = rma(data)
data.matrix = exprs(data.rma)

#PCA

color=c('green','green','green','green','red','red','red','red')
data.PCA = prcomp(t(data.matrix),scale.=TRUE)
#not used: plot(data.PC$x[1:2],col=color)
# sqrt of eigenvalues
data.PCA$sdev

#loadings
head(data.PCA$rotation)

#PCs (aka scores)
head(data.PCA$x)

ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores))) +
  geom_hline(yintercept = 0, colour = "gray") +
  geom_vline(xintercept = 0, colour = "gray") +
  geom_text(colour = "blue", alpha = 0.8, size = 4) +
  ggtitle("PCA plot microarray")


