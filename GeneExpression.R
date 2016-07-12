library(DESeq2)
library(WGCNA)

##Read in raw counts
counts <- read.delim("RawCounts.txt",header=T,row.names=1)

##Read in Conditions
conds <- read.csv("Individual_Metadata.csv",header=T)
condsOrder <- conds[match(colnames(counts),conds$Individual),]
rownames(condsOrder) <- condsOrder$Individual
condsOrder <- data.frame(apply(condsOrder,2,factor))

##Normalize counts in DESeq2
counts1 <- subset(subCounts,rowMeans(subCounts)>1) # remove low coverage contigs
dds <- DESeqDataSetFromMatrix(counts1,condsOrder,design=~Pool)
dds <- DESeq(dds)
scaledCounts <- counts(dds,normalized=T)
colnames(scaledCounts) <- colnames(counts1)

#########################
###      WGCNA     ######
#########################

options(stringsAsFactors=FALSE)
enableWGCNAThreads()
sub <- scaledCounts
datExpr <- as.data.frame(t(sub))

##Check data for missing values
gsg=goodSamplesGenes(datExpr,verbose=3)
gsg$allOK

##Cluster samples for obvious outliers
sampleTree <- hclust(dist(datExpr),method="average")
plot(sampleTree)

##Trait data that we want
rownames(condsOrder) <- condsOrder[,1]
subConds <- condsOrder[,c("Origin","HVgrowth","MVgrowth","HVsurvival","MVsurvival")]
traitColors <- numbers2colors(subConds,signed=F)
plotDendroAndColors(sampleTree,traitColors,groupLabels=names(subConds))

##Soft-thresholding
powers=c(c(1:10),seq(12,20,by=2))
sft=pickSoftThreshold(datExpr,powerVector=powers,verbose=5)
par(mfrow=c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
	xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
	type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	labels=powers,col="red")	
abline(h=0.95,col="red")
plot(sft$fitIndices[,1],sft$fitIndices[,5],type="n",
 	xlab="Soft Threshold (power)",ylab="Mean Connectivity",
 	main="Mean Connectivity")
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,col="red")

###Network construction
##Note: I did this on a computer with 64GB RAM. For smaller computers, change maxBlockSize
net <- blockwiseModules(datExpr, power=4,
		TOMType="unsigned",minModuleSize=30,
		reassignThreshold=0,mergeCutHeight=0.3,
		numericLabels=T, pamRespectsDendro=F,
		verbose=3, maxBlockSize=24500)
table(bwnet$colors)
mergedColors=labels2colors(net$colors)
		
##Save Module Information
moduleLabels=net$colors
moduleColors=labels2colors(net$colors)
MEs=net$MEs
geneTree=net$dendrograms[[1]]

##module-trait associations
MEs0=moduleEigengenes(datExpr,moduleLabels)$eigengenes
MEs=orderMEs(MEs0)

###Plot and examine correlations with fitness traits - can do this for all 4 fitness traits
trait = condsOrder$MVgrowth
par(mfrow=c(5,5),mar=c(1,1,1,1))
for (i in 1:ncol(MEs)) {
	plot(MEs[,i]~trait,pch=19,cex=2)
	print(i)
	print(summary(lm(MEs[,i]~trait)))
}

##Save the data
write.table(MEs,file="Cluster_Eigenvalues.txt",quote=F,row.names=T,col.names=T,sep="\t")
clusterID <- data.frame(cbind(scaledCounts,cluster=net$colors))
write.table(clusterID,file="Contig_clusters.txt",quote=F,row.names=T,col.names=T,sep="\t")

