##import transplant data
trans <- read.delim("SuppTab1.txt",header=T)

###Permute branches and calculate standard deviation across individuals within each pool of origin
n=1000
sdpermHV <- c()
sdpermMV <- c()
for (i in 1:n) {
surv <- sample(trans$Survival,length(trans$Survival),replace=F)
agg <- aggregate(surv,list(Col=trans$Colony,Or=trans$Origin),mean,na.rm=T)
sdpermHV[i] <- sd(agg$x[which(agg$Or=="HV")])
sdpermMV[i] <- sd(agg$x[which(agg$Or=="MV")])
#print(i)
}

##Calculate real standard deviation within each pool
realSurv <- aggregate(trans$Survival,list(Col=trans$Colony,Or=trans$Origin),mean,na.rm=T)
realsd=sd(ag$x)
realsdHV <- sd(ag$x[which(ag$Or=="HV")])
realsdMV <- sd(ag$x[which(ag$Or=="MV")])

##Calculate permutational p-values
pHV <- length(which(sdpermHV>realsdHV))/n
pHV
pMV <- length(which(sdpermMV>realsdMV))/n
pMV
