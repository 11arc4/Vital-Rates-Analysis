#Un-correlated Life Stage Simulation Analysis for the Appendix
#this is updated to use the most accurate survival estimates from MARK

setwd("~/Masters Thesis Project/TRES Data Analysis/Vital Rates Models/RMark Preliminary Survival Analysis")

library(dplyr)
library(popbio)
library(ggplot2)

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

#Read in the data, including all the different types of data from RMark
PopData<- read.csv(file= "file:///C:/Users/11arc/Documents/Masters Thesis Project/TRES Data Analysis/Matrix Pop Estimates/Yearly Vital Rates.csv", na.strings=c(""), as.is=T)
#Need to get rid of all the columns that are specific to SY or ASY
PopData <- PopData[, -c(2, 5, 9, 10, 12, 13)]

FemaleSurvival <- read.csv("file:///C:/Users/11arc/Documents/Masters Thesis Project/TRES Data Analysis/Vital Rates Models/RMark Preliminary Survival Analysis/Model Averaged Female Yearly Apparent Survival_updated.csv", as.is=T)
#These are the best estimates from our model averaged female survival analysis
NestlingSurvival <- read.csv("file:///C:/Users/11arc/Documents/Masters Thesis Project/TRES Data Analysis/Vital Rates Models/RMark Preliminary Survival Analysis/Model Averaged Yearly Apparent Survival for all records to get recruitment.csv", as.is=T)
#These are the best estimates from our model averaged anlysis of ALL birds (male
#and female, over the entire course of the study). This is good for getting a
#better estimate of nestling recruitment
FemaleSurvival$SYReturn[which(FemaleSurvival$seSYReturn > 0.2)]<- NA
FemaleSurvival$ASYReturn[which(FemaleSurvival$seASYReturn > 0.2)]<- NA


PopData$SYReturn<- NA
PopData$ASYReturn<- NA
PopData$Recruitment <- NA
PopData$SYReturn[1:42] <- FemaleSurvival$SYReturn
PopData$ASYReturn[1:42] <- FemaleSurvival$ASYReturn
PopData$Recruitment[1:42] <-NestlingSurvival$Recruitment



recruitmentParameters <- estBetaParams(mu=mean(PopData$Recruitment, na.rm=T), var=var(PopData$Recruitment, na.rm=T))
SYReturnParameters<- estBetaParams(mu=mean(PopData$SYReturn, na.rm=T ), var=var(PopData$SYReturn, na.rm=T))
ASYReturnParameters<- estBetaParams(mu=mean(PopData$ASYReturn, na.rm=T), var=var(PopData$ASYReturn, na.rm=T))
fledgeParameters <- estBetaParams(mu=mean(PopData$fledgeRate), var=var(PopData$fledgeRate))
hatchrateParameters <- estBetaParams(mu=mean(PopData$hatchRate), var=var(PopData$hatchRate))
averageNestParametersSY<- estBetaParams(mu=mean(PopData$averageNestsSY-1, na.rm=T), var=var(PopData$averageNestsSY-1, na.rm=T))
averageNestParametersASY<- estBetaParams(mu=mean(PopData$averageNestsASY-1), var=var(PopData$averageNestsASY-1))
betaClutchSizeParametersSY <- estBetaParams(mu=(mean(PopData$clutchSizeSY, na.rm=T)-min(PopData$clutchSizeSY, na.rm=T))/(max(PopData$clutchSizeSY, na.rm=T)-min(PopData$clutchSizeSY, na.rm=T)), 
                                            var=(var(PopData$clutchSizeSY, na.rm=T)* (1/ (max(PopData$clutchSizeSY, na.rm=T)-min(PopData$clutchSizeSY, na.rm=T)))^2))
betaClutchSizeParametersASY <- estBetaParams(mu=(mean(PopData$clutchSizeASY, na.rm=T)-min(PopData$clutchSizeASY, na.rm=T))/(max(PopData$clutchSizeASY, na.rm=T)-min(PopData$clutchSizeASY, na.rm=T)), 
                                             var=(var(PopData$clutchSizeASY, na.rm=T)* (1/(max(PopData$clutchSizeASY, na.rm=T)-min(PopData$clutchSizeASY, na.rm=T)))^2))




#Vital Rates analysis drawing repeatedly from the known distributions Currently
#we are not seperating out SY and ASY clutch sizes etc, even though ultimately
#we may want to.

#This is where we need to start looping! and collecting data for the vital rates analysis!
vrdat <- as.data.frame(matrix(nrow=10000, ncol=28)) 
#Doing 10,000 because that's what Taylor et al 2012 did)
VitalRates <- c("SYClutch",
                "SYNests", 
                "ASYClutch", 
                "ASYNests", 
                "Hatch", 
                "Fledge", 
                "Recruit", 
                "SYReturn", 
                "ASYReturn")
colnames<- rep(NA, 3*length(VitalRates)+1)
for(i in 1:length(VitalRates)){
  colnames[0+i]<- VitalRates[i]
  colnames[length(VitalRates)+i]<- paste("Sens", VitalRates[i], sep="")
  colnames[2*length(VitalRates)+i]<- paste("Elas", VitalRates[i], sep="")
}


colnames[28]<- "lambda"
colnames[1:7]<-c("averageNestsSY",  "averageNestsASY", "clutchSizeSY"  ,  "clutchSizeASY"   ,"hatchrate"    ,   "fledgerate"      ,"recruitrate" )
colnames(vrdat)<- colnames

MatrixElements <- c("SYFertility", "ASYFertility", "NestlingSurvival", "SYSurvival", "ASYSurvival")


for (i in 1:nrow(vrdat)){
  
  vrdat$averageNestsSY[i] <- 1+ rbeta (n=1, shape1= averageNestParametersSY$alpha , shape2=averageNestParametersSY$beta) #oesn't work because neither should be negative! What's wrong with this?
  vrdat$averageNestsASY[i] <- 1+ rbeta (n=1, shape1= averageNestParametersASY$alpha , shape2=averageNestParametersASY$beta) #oesn't work because neither should be negative! What's wrong with this?
  
  
  
  betavalSY <- rbeta(1, shape1= betaClutchSizeParametersSY$alpha , shape2=betaClutchSizeParametersSY$beta)
  
  vrdat$clutchSizeSY[i]<- 0.5*betavalSY*(max(PopData$clutchSizeSY, na.rm=T)-min(PopData$clutchSizeSY, na.rm=T))+min(PopData$clutchSizeSY, na.rm=T) #as per morris and doak pg 283
  
  betavalASY <- rbeta(1, shape1= betaClutchSizeParametersASY$alpha , shape2=betaClutchSizeParametersASY$beta)
  
  vrdat$clutchSizeASY[i]<- 0.5* betavalASY*(max(PopData$clutchSizeASY, na.rm=T)-min(PopData$clutchSizeASY, na.rm=T))+min(PopData$clutchSizeASY, na.rm=T) #as per morris and doak pg 283
  
  
  #Pull out a a random ratch rate, fledge rate, renest status and clutch size from the vectors I made above. 
  vrdat$hatchrate[i] <- rbeta (n=1, shape1= hatchrateParameters$alpha , shape2=hatchrateParameters$beta)
  
  vrdat$fledgerate[i] <- rbeta (n=1, shape1= fledgeParameters$alpha , shape2=fledgeParameters$beta)
  #recruitment is randomly generated from the data that I brought to RMark
  vrdat$recruitrate[i] <-  rbeta (n=1, shape1= recruitmentParameters$alpha , shape2=recruitmentParameters$beta)
  
  
  
  vrdat$SYReturn[i] <- rbeta (n=1, shape1= SYReturnParameters$alpha , shape2=SYReturnParameters$beta)
  vrdat$ASYReturn[i] <- rbeta (n=1, shape1= ASYReturnParameters$alpha , shape2=ASYReturnParameters$beta)
  
  
  
  vrdat$layrateSY[i]<- vrdat$clutchSizeSY[i]*vrdat$averageNestsASY[i]
  vrdat$layrateASY[i]<- vrdat$clutchSizeASY[i]*vrdat$averageNestsASY[i]
  
  vrdat$eggRecruitRate[i] <- vrdat$hatchrate[i] * vrdat$fledgerate[i] * vrdat$recruitrate[i]
  
  
  #Put everything into a matrix, calculate sensitivity and elasticity
  tres.vr <- list(SYClutch=vrdat$clutchSizeSY[i], 
                  SYNests=vrdat$averageNestsSY[i], 
                  ASYClutch=vrdat$clutchSizeASY[i], 
                  ASYNests=vrdat$averageNestsASY[i], 
                  Hatch=vrdat$hatchrate[i], 
                  Fledge=vrdat$fledgerate[i], 
                  Recruit=vrdat$recruitrate[i], 
                  SYReturn=vrdat$SYReturn[i], 
                  ASYReturn=vrdat$ASYReturn[i])
  
  tres.e1 <- expression(
    0, SYClutch*SYNests, ASYClutch*ASYNests, 
    Hatch*Fledge*Recruit, 0, 0, 
    0, SYReturn, ASYReturn
  )
  matrixEl<-sapply(tres.e1, eval, tres.vr, NULL)
  
  A <- matrix(matrixEl, nrow=sqrt(length(tres.e1)), byrow=TRUE)
  
  x<- vitalsens(tres.e1, tres.vr)
  
  vrdat$lambda[i] <- lambda(A)
  vrdat[i, 10:18]<- x$sensitivity
  vrdat[i,19:27]<- x$elasticity
  
  vrdat$VRMaxSens[i] <- VitalRates[order(vrdat[i,10:18] , decreasing=T)][1]
  vrdat$VR2MaxSens[i] <- VitalRates[order(vrdat[i,10:18] , decreasing=T)][2]
  vrdat$VR3MaxSens[i] <- VitalRates[order(vrdat[i,10:18] , decreasing=T)][3]
  vrdat$MatrixMaxElas[i] <- MatrixElements[which(vrdat[i,c(19, 21, 23, 26, 27)] == max((vrdat[i,c(19, 21, 23, 26, 27)])))]
  
}


vrdat$VRMaxSens <- as.factor(vrdat$VRMaxSens)
vrdat$VR2MaxSens <- as.factor(vrdat$VR2MaxSens)
vrdat$VR3MaxSens <- as.factor(vrdat$VR3MaxSens)
vrdat$MatrixMaxElas <- as.factor(vrdat$MatrixMaxElas)






#Now let's make a plot with all this data as proportional bar graphs
#That is going to involve restucturing the data slightly. 

vrdat2 <- reshape2::melt(vrdat[,32:34], measure.vars=colnames(vrdat)[32:34], variable.name="Rank", value.name="VitalRate")
vrdat3 <- vrdat2 %>% group_by(Rank) %>% summarise(ASYReturnTot= length(which(VitalRate=="ASYReturn")), 
                                                  SYReturnTot=length(which(VitalRate=="SYReturn")),
                                                  RecruitTot= length(which(VitalRate=="Recruit")),
                                                  FledgeTot = length(which(VitalRate=="Fledge")),
                                                  HatchTot = length(which(VitalRate=="Hatch")))
vrdat4 <- reshape2::melt(vrdat3, id="Rank", variable.name="VitalRate", value.name="Total")

vrdat4$VitalRate <- factor(vrdat4$VitalRate, levels=c("RecruitTot", "ASYReturnTot", "SYReturnTot", "FledgeTot", "HatchTot"))

SensRankPlot <- ggplot(vrdat4, aes(x=Rank, y=Total, fill=VitalRate))+
  geom_bar(position=c("fill"), stat="identity")+
  scale_x_discrete("Sensitivity ranking by vital rate", 
                   labels=c("VRMaxSens"="1st", "VR2MaxSens"="2nd", "VR3MaxSens"="3rd"))+
  scale_fill_grey(start=0.3, end=0.8, labels=c("ASYReturnTot"= "Older return", 
                                               "SYReturnTot"= "1-year-old return",
                                               "RecruitTot"= "Recruitment", 
                                               "FledgeTot"="Fledge rate", 
                                               "HatchTot"="Hatch rate"), 
                  name="Vital Rate")+
  ylab("Proportion")+
  ggthemes::theme_few(base_size = 20)

png(filename = "~/Masters Thesis Project/Vital Rates Paper/Figures/Sensitivity Rankings Bar Plot_uncorrelatedUPDATED.png", 
    width=500, height=500)
SensRankPlot
dev.off()


####Lets do the same thing now for elasticity but for each matrix element since
####vital rates within a particular matrix element have equivalent elasticity
vrdat5 <- as.data.frame(matrix(nrow=5, ncol=2, NA))
colnames(vrdat5)<- c("MatrixElement", "Total")
vrdat5$MatrixElement <- MatrixElements
vrdat5$Total[1] <-  length(which(vrdat$MatrixMaxElas=="SYFertility"))
vrdat5$Total[2] <- length(which(vrdat$MatrixMaxElas=="ASYFertility"))
vrdat5$Total[3] <- length(which(vrdat$MatrixMaxElas=="NestlingSurvival"))
vrdat5$Total[4] <- length(which(vrdat$MatrixMaxElas=="SYSurvival"))
vrdat5$Total[5] <- length(which(vrdat$MatrixMaxElas=="ASYSurvival")) 


#make the plot
ElasRankPlot <-  ggplot(vrdat5 %>% filter (Total>0), aes(x=factor(1), y=Total, fill=MatrixElement))+
  geom_bar(position=c("fill"), stat="identity")+
  scale_fill_grey(start=0.3, end=0.8, labels=c("ASYSurvival"= "Older survival", 
                                               "NestlingSurvival"= "Nestling survival"), 
                  name="Matrix Element")+
  scale_x_discrete("Elasticity ranking 1st \n by vital rate", 
                   labels=c("1"=""))+
  
  ylab("Proportion")+
  
  ggthemes::theme_few(base_size = 20)

png(filename = "~/Masters Thesis Project/Vital Rates Paper/Figures/Elasticity Rankings Bar Plot_uncorrelatedUPDATED.png", 
    width=300, height=500)
ElasRankPlot
dev.off()




############Now lets do a deterministic sensitivity analysis on this data

#Calculate sensitivity according to Taylor et al. 2012 using the linear regression
modASYReturn <- lm(lambda~ ASYReturn, data=vrdat)
sumASY <- summary(modASYReturn) #Adjusted R^2 = 0.003741 so quite unimportant

modSYReturn <- lm(lambda~ SYReturn, data=vrdat)
sumSY <- summary(modSYReturn) # Adjusted R^2 = 0.002058 so even less important

modClutchSizeSY <- lm(lambda ~ clutchSizeSY, data=vrdat)
sumClutchSY<- summary(modClutchSizeSY) #Adjusted R^2 = 0.2067 so it's not that!

modClutchSizeASY <- lm(lambda ~ clutchSizeASY, data=vrdat)
sumClutchASY<- summary(modClutchSizeASY) #Adjusted R^2 = 0.1071 so it's not that!


modAverageNestsSY<- lm(lambda ~ averageNestsSY, data=vrdat) 
sumAverageNestsSY <- summary(modAverageNestsSY) #Adjusted R^2 =0.01887

modAverageNestsASY<- lm(lambda ~ averageNestsASY, data=vrdat) 
sumAverageNestsASY <- summary(modAverageNestsASY) #Adjusted R^2 =0.1284

modEggtoAdult <- lm(lambda~eggRecruitRate, data=vrdat)
sumEggtoAdult<- summary(modEggtoAdult) #Adjusted R^2 = 0.9091 so explains almost all the variation!
#What within this explains variation? 

modHatch <- lm(lambda ~ hatchrate, data=vrdat)
sumHatch<- summary(modHatch) # Adjussted R^2 =0.0231

modFledge <- lm (lambda ~fledgerate, data=vrdat)
sumFledge <- summary(modFledge) #Adjusted R^2 =0.092 so another large portion of the variation, 

modRecruit <- lm(lambda~recruitrate, data=vrdat)
sumRecruit <- summary(modRecruit) #R^2= 0.7669



#Calculate Elasticity according to Taylor et al. 2012 using the linear regression of the log log transformation
modASYReturn_elasticity <- lm(log(lambda)~ log(ASYReturn), data=vrdat)
sumASY_e<- summary(modASYReturn_elasticity) #R^2=0.0001127

modSYReturn_elasticity <- lm(log(lambda) ~ log(SYReturn), data=vrdat)
sumSY_e <- summary(modSYReturn_elasticity) #R^2 = 0.0002805

modClutchSizeSY_elasticity <- lm(log(lambda) ~ log(clutchSizeSY), data=vrdat)
sumClutchSY_e<- summary(modClutchSizeSY_elasticity) # R^2 = 0.001324 so it's not that!

modClutchSizeASY_elasticity <- lm(log(lambda) ~ log(clutchSizeASY), data=vrdat)
sumClutchASY_e<- summary(modClutchSizeASY_elasticity) # R^2 = 0.01159 so it's not that!

modAverageNestsSY_elasticity<- lm(log(lambda) ~ log(averageNestsSY), data=vrdat) 
sumAverageNestsSY_e<- summary(modAverageNestsSY_elasticity) #Adjusted R^2 =0.0002755

modAverageNestsASY_elasticity<- lm(log(lambda) ~ log(averageNestsASY), data=vrdat) 
sumAverageNestsASY_e<- summary(modAverageNestsASY_elasticity) #Adjusted R^2 =0.01694, better but not a whole lot of the variation

modEggtoAdult_elasticity <- lm(log(lambda)~log(eggRecruitRate), data=vrdat)
sumEggtoAdult_e <- summary(modEggtoAdult_elasticity) #R^2 = 0.5758

modhatch_elasticity <- lm(log(lambda)~log(hatchrate), data=vrdat)
sumHatch_e <- summary(modhatch_elasticity) #R^2 =0.01803 so low elasticity

modFledge_elasticity <- lm (log(lambda) ~log(fledgerate), data=vrdat)
sumFledge_e <- summary(modFledge_elasticity) #Adjusted R^2 =0.08591 so another large portion of the variation, 

modRecruit_elasticity <- lm(log(lambda)~log(recruitrate), data=vrdat)
sumRecruit_e<- summary(modRecruit_elasticity) #r^2=0.5271

#Put all of this analysis into a seperate dataframe
SensitivityAnalysis<- as.data.frame(matrix(nrow=9, ncol=5))
colnames(SensitivityAnalysis)<- c("VitalRate", "Sensitivity", "Elasticity", "R2_S", "R2_E")

SensitivityAnalysis$VitalRate<- c("Recruitment", 
                                  "SY Return",
                                  "ASY Return", 
                                  "SY Nests", 
                                  "ASY Nests",
                                  "SY Clutch Size", 
                                  "ASY Clutch Size", 
                                  "Hatch Rate", 
                                  "Fledge Rate"
)

SensitivityAnalysis$Sensitivity[1]<- modRecruit$coefficients[2]
SensitivityAnalysis$R2_S[1] <- sumRecruit$adj.r.squared
SensitivityAnalysis$Sensitivity[2]<- modSYReturn$coefficients[2]
SensitivityAnalysis$R2_S[2] <- sumSY$adj.r.squared

SensitivityAnalysis$Sensitivity[3]<- modASYReturn$coefficients[2]
SensitivityAnalysis$R2_S[3] <- sumASY$adj.r.squared

SensitivityAnalysis$Sensitivity[4]<- modAverageNestsSY$coefficients[2]
SensitivityAnalysis$R2_S[4] <- sumAverageNestsSY$adj.r.squared

SensitivityAnalysis$Sensitivity[5]<- modAverageNestsASY$coefficients[2]
SensitivityAnalysis$R2_S[5] <- sumAverageNestsASY$adj.r.squared

SensitivityAnalysis$Sensitivity[6]<- modClutchSizeSY$coefficients[2]
SensitivityAnalysis$R2_S[6] <- sumClutchSY$adj.r.squared

SensitivityAnalysis$Sensitivity[7]<- modClutchSizeASY$coefficients[2]
SensitivityAnalysis$R2_S[7] <- sumClutchASY$adj.r.squared

SensitivityAnalysis$Sensitivity[8]<- modHatch$coefficients[2]
SensitivityAnalysis$R2_S[8] <- sumHatch$adj.r.squared

SensitivityAnalysis$Sensitivity[9]<- modFledge$coefficients[2]
SensitivityAnalysis$R2_S[9] <- sumFledge$adj.r.squared



SensitivityAnalysis$Elasticity[1]<- modRecruit_elasticity$coefficients[2]
SensitivityAnalysis$R2_E[1] <- sumRecruit_e$adj.r.squared
SensitivityAnalysis$Elasticity[2]<- modSYReturn_elasticity$coefficients[2]
SensitivityAnalysis$R2_E[2] <- sumSY_e$adj.r.squared

SensitivityAnalysis$Elasticity[3]<- modASYReturn_elasticity$coefficients[2]
SensitivityAnalysis$R2_E[3] <- sumASY_e$adj.r.squared

SensitivityAnalysis$Elasticity[4]<- modAverageNestsSY_elasticity$coefficients[2]
SensitivityAnalysis$R2_E[4] <- sumAverageNestsSY_e$adj.r.squared

SensitivityAnalysis$Elasticity[5]<- modAverageNestsASY_elasticity$coefficients[2]
SensitivityAnalysis$R2_E[5] <- sumAverageNestsASY_e$adj.r.squared

SensitivityAnalysis$Elasticity[6]<- modClutchSizeSY_elasticity$coefficients[2]
SensitivityAnalysis$R2_E[6] <- sumClutchSY_e$adj.r.squared

SensitivityAnalysis$Elasticity[7]<- modClutchSizeASY_elasticity$coefficients[2]
SensitivityAnalysis$R2_E[7] <- sumClutchASY_e$adj.r.squared

SensitivityAnalysis$Elasticity[8]<- modhatch_elasticity$coefficients[2]
SensitivityAnalysis$R2_E[8] <- sumHatch_e$adj.r.squared

SensitivityAnalysis$Elasticity[9]<- modFledge_elasticity$coefficients[2]
SensitivityAnalysis$R2_E[9] <- sumFledge_e$adj.r.squared

SensitivityAnalysis$Label <- c(7,8,9,1,2,3,4,5,6)
SensitivityAnalysis$Location[1:3]<- "Over wintering"

SensitivityAnalysis$Location[4:9]<- "Breeding"



SensitivityPlot <- ggplot(SensitivityAnalysis, aes(x=Sensitivity, y=R2_S))+
  xlab("Sensitivity")+
  ylab(expression(paste("Adjusted R"^ 2)))+
  geom_point(size=3)+ 
  ggrepel::geom_label_repel(data=SensitivityAnalysis %>% filter(Sensitivity>0.5), aes(label=VitalRate), size=6, min.segment.length = 0, box.padding = 0.5)+
  ggrepel::geom_label_repel(data=SensitivityAnalysis %>% filter(Sensitivity<0.5), aes(label=VitalRate), size=6, min.segment.length = 0, box.padding = 1)+
  geom_vline(xintercept = 0, linetype="dashed")+
  theme_classic(base_size = 20)+
  theme(axis.title.y=element_text(angle=0, vjust=0.5))
ElasticityPlot <- ggplot(SensitivityAnalysis, aes(x=Elasticity, y=R2_E))+
  xlab("Elasticity")+
  ylab(expression(paste("Adjusted R"^ 2)))+
  geom_point(size=3)+ 
  ggrepel::geom_label_repel(data=SensitivityAnalysis %>% filter(R2_E>0.3) ,aes(label=VitalRate), size=6, min.segment.length = 0,  box.padding = .5)+
  ggrepel::geom_label_repel(data=SensitivityAnalysis %>% filter(R2_E<0.3 & Elasticity<0), aes(label=VitalRate), size=6, min.segment.length = 0,  box.padding = .8)+
  ggrepel::geom_label_repel(data=SensitivityAnalysis %>% filter(R2_E<0.3 & Elasticity>0.2), aes(label=VitalRate), size=6, min.segment.length = 0,  box.padding = 1.5)+
  ggrepel::geom_label_repel(data=SensitivityAnalysis %>% filter(R2_E<0.3 & Elasticity>0 & Elasticity<0.2), aes(label=VitalRate), size=6, min.segment.length = 0,  box.padding = .8)+
  geom_vline(xintercept = 0, linetype="dashed")+
  theme_classic(base_size = 20)+
  theme(axis.title.y=element_text(angle=0, vjust=0.5))


#Yes this is still quite different than when we correlate the vital rates more
#appropriately. We will relegate this to the supplemental



SensitivityPlot2 <- ggplot(SensitivityAnalysis, aes(x=Sensitivity, y=R2_S))+
  geom_vline(xintercept = 0, linetype="dashed")+
  xlab("Sensitivity")+
  ylab(expression(paste("Adjusted R"^ 2)))+
  ggrepel::geom_label_repel(data=SensitivityAnalysis %>% filter(Sensitivity>0.5), aes(label=VitalRate), size=8,  box.padding = 0.5, segment.colour = "black", show.legend = F)+
  ggrepel::geom_label_repel(data=SensitivityAnalysis %>% filter(Sensitivity<0.5), aes(label=VitalRate), size=8, box.padding=2, force=5, segment.colour = "black", show.legend = F)+
  geom_point(size=5)+ 
  theme_classic(base_size = 20)+
  theme(axis.title.y=element_text(angle=0, vjust=0.5))+
  xlim(-.5, 2.5)
SensitivityPlot2


ElasticityPlot2 <- ggplot(SensitivityAnalysis, aes(x=Elasticity, y=R2_E))+
  xlab("Elasticity")+
  ylab(expression(paste("Adjusted R"^ 2)))+
  geom_vline(xintercept = 0, linetype="dashed")+
  ggrepel::geom_label_repel(data=SensitivityAnalysis , aes(label=VitalRate), size=8, min.segment.length = 0,  box.padding = 1.5, segment.color="black", show.legend = F)+
  geom_point(size=5)+ 
  theme_classic(base_size = 20)+
  theme(axis.title.y=element_text(angle=0, vjust=0.5))+
  xlim(-1, 0.8)

ElasticityPlot2

png(filename = "~/Masters Thesis Project/Vital Rates Paper/Figures/Sensitivity Analysis_uncorrelated2.png", width=1000, height=1000)
cowplot::plot_grid(SensitivityPlot2, ElasticityPlot2, align="v",nrow=2, ncol=1, labels=c("A", "B"))
dev.off()
