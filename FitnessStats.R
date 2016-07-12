library(lme4)
library(lmerTest)
library(WGCNA)

##Read in data
data <- read.delim("SuppTab1.txt",header=T)


##############
##	GLMM	##
##############

###Survival
glmerControl=glmerControl(optimizer="bobyqa")
surv.full <- glmer(Survival~Origin*Location+(1|Crate)+(1|Colony),data=data,
	family="binomial",control=glmerControl) ##Full model
surv.M1 <- glmer(Survival~Origin+(1|Colony)+(1|Crate),data=data,
	family="binomial",control=glmerControl) ##Best model
AIC(surv.full,surv.M1)
summary(surv.M1)

###Growth
growth.full <- lmer(BW~Origin*Location+(1|Crate)+(1|Colony),data=data) ##Full model
growth.M1 <- lmer(BW~Location+(1|Crate)+(1|Colony),data=data)  ##Best model
AIC(growth.full,growth.M1)
summary(growth.M1)

##########################
##	Crate Normalization	##
##########################
rownames(data) <- paste(data$Colony,data$Crate,sep="")

goodData <- data[!is.na(data$Survival),]

##Survival
HVcrates <- lm(scale(goodData$Survival[goodData$Location=="HV"],center=T,scale=F)~goodData$Crate[goodData$Location=="HV"])$residuals# + mean(goodData$BW[goodData$Location=="HV"],na.rm=T) 
MVcrates <- lm(scale(goodData$Survival[goodData$Location=="MV"],center=T,scale=F)~goodData$Crate[goodData$Location=="MV"])$residuals# + mean(goodData$BW[goodData$Location=="MV"],na.rm=T)
surv.crates <- c(HVcrates,MVcrates)
survival <- cbind(goodData,surv.crates)

##Growth 
survivors <- data[!is.na(data$BW),]
HVcrates <- lm(scale(survivors$BW[survivors$Location=="HV"],center=T,scale=F)~survivors$Crate[survivors$Location=="HV"])$residuals# + mean(survivors$BW[survivors$Location=="HV"]) 
MVcrates <- lm(scale(survivors$BW[survivors$Location=="MV"],center=T,scale=F)~survivors$Crate[survivors$Location=="MV"])$residuals# + mean(survivors$BW[survivors$Location=="MV"])
growth.crates <- c(HVcrates,MVcrates)
growth <- cbind(survivors,growth.crates)

##Colony means
frame <- cbind(data,normSurv=survival$surv.crates[match(rownames(data),rownames(survival))],
		normGrowth=growth$growth.crates[match(rownames(data),rownames(growth))])
means <- aggregate(frame[,8:9],list(Colony=frame$Colony,Origin=frame$Origin,Location=frame$Location),mean,na.rm=T)
meansframe <- cbind(meansHV=means[1:21,],meansMV=means[22:42,])

##MV Growth vs. HV survival
summary(lm((meansframe[,10]~meansframe[,4])))
plot(meansframe[,10]~meansframe[,4])


