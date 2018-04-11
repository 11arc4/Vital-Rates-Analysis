#Correlated Life Stage Simulation Analysis
#this is updated to use the most accurate survival estimates from MARK

setwd("~/Masters Thesis Project/TRES Data Analysis/Vital Rates Models/RMark Preliminary Survival Analysis")

library(dplyr)
library(popbio)
library(ggplot2)
#Function to estimate the parameters for a beta distribution
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}
#this is just the original eigen function from base R that has been modified so 
#it doesn't order the eigen values and vectors by size of eigen value (this is
#what matlab does and therefore what I need to do in order to match the Morris
#and Doak matlab code properly)
eigenUnOrdered <- function (x, symmetric, only.values = FALSE, EISPACK = FALSE) 
{
  x <- unname(as.matrix(x))
  n <- nrow(x)
  if (!n) 
    stop("0 x 0 matrix")
  if (n != ncol(x)) 
    stop("non-square matrix in 'eigen'")
  n <- as.integer(n)
  if (is.na(n)) 
    stop("invalid nrow(x)")
  complex.x <- is.complex(x)
  if (!all(is.finite(x))) 
    stop("infinite or missing values in 'x'")
  if (missing(symmetric)) 
    symmetric <- isSymmetric.matrix(x)
  if (symmetric) {
    z <- if (!complex.x) 
      .Internal(La_rs(x, only.values))
    else .Internal(La_rs_cmplx(x, only.values))
    ord <- rev(seq_along(z$values))
  }
  else {
    z <- if (!complex.x) 
      .Internal(La_rg(x, only.values))
    else .Internal(La_rg_cmplx(x, only.values))
    ord <- sort.list(Mod(z$values), decreasing = TRUE)
  }
  return(list(values = z$values, vectors = if (!only.values) z$vectors))
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


write.csv(PopData, row.names=F, na="", file="file:///C:/Users/11arc/Documents/Masters Thesis Project/Vital Rates Paper/Yearly Vital Rates Estimates.csv")


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



#Follow directions in Morris and Doak 2002 Quantitative Conservation Biology pg 284 

#1. get the correlation matrix of all the vital rates (not inluding the year!)
c <- cor(PopData[,-c(1)], method = "spearman", use="na.or.complete")


#2. make a matrix (W) who's columns are all the possible right eigenvectors of C
eigC <- eigenUnOrdered(c)
W<- eigC$vectors

#3. make a matrix (D) who's has the eigen values along the diagonal and all else 0
#check to make sure that all eigen values are >=0 first, if less than they are likely artificial and should be changed to 0
#All of my eigen values are >0 so life is good and we can carry on
sqrt(abs(eigC$vectors))
D<- matrix(data=0, nrow=9, ncol=9)
for (i in 1:9){
  D[i,i]<- eigC$values[i]
  
}


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


for (i in 1:nrow(vrdat)){
  
  
  #3. Generate a vector (m) of uncorrelation standard normal variables. 
  m <- rnorm(n=9, mean=0, sd=1)
  
  #4. Multiply m by C 1/2 (ie W * D 1/2 * W') to get a set of correlated standard normal variables (y)
  c12<- W%*%sqrt(D)%*% t(W)
  
  y <- m%*%c12
  
  #5 None of our actual distributions were standard normal though so now we have
  #to find the F stat for the standard normal value given, and then find the value
  #for the equivalent F stat in the appropriate distribution.
  colnames(c)
  
  #averageNests goes first
  FaverageNestsSY <- pnorm(y[1], mean=0, sd=1)
  vrdat$averageNestsSY[i] <- 1+qbeta(FaverageNestsSY, shape1= averageNestParametersSY$alpha , shape2=averageNestParametersSY$beta, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  
  FaverageNestsASY <- pnorm(y[2], mean=0, sd=1)
  vrdat$averageNestsASY[i] <- 1+qbeta(FaverageNestsASY, shape1= averageNestParametersASY$alpha , shape2=averageNestParametersASY$beta, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  
  #clutchSizes 
  FClutchSizeSY <- pnorm(y[3], mean=0, sd=1)
  betavalSY <- qbeta(FClutchSizeSY, shape1= betaClutchSizeParametersSY$alpha , shape2=betaClutchSizeParametersSY$beta, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  
  vrdat$clutchSizeSY[i]<- 0.5*betavalSY*(max(PopData$clutchSizeSY, na.rm=T)-min(PopData$clutchSizeSY, na.rm=T))+min(PopData$clutchSizeSY, na.rm=T) #as per morris and doak pg 283
  
  FClutchSizeASY <- pnorm(y[4], mean=0, sd=1)
  betavalASY <- qbeta(FClutchSizeASY, shape1= betaClutchSizeParametersASY$alpha , shape2=betaClutchSizeParametersASY$beta, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  
  vrdat$clutchSizeASY[i]<- 0.5* betavalASY*(max(PopData$clutchSizeASY, na.rm=T)-min(PopData$clutchSizeASY, na.rm=T))+min(PopData$clutchSizeASY, na.rm=T) #as per morris and doak pg 283
  
  #hatchRate 
  FHatchRate <- pnorm(y[5], mean=0, sd=1)
  vrdat$hatchrate[i] <- qbeta(FHatchRate, shape1= hatchrateParameters$alpha , shape2=hatchrateParameters$beta, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  
  
  #fledgeRate 
  FFledgeRate <- pnorm(y[6], mean=0, sd=1)
  vrdat$fledgerate[i] <- qbeta(FFledgeRate, shape1= fledgeParameters$alpha , shape2=fledgeParameters$beta, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  
  
  #SYReturn
  FSYReturn <- pnorm(y[7], mean=0, sd=1)
  vrdat$SYReturn[i] <- qbeta(FSYReturn, shape1= SYReturnParameters$alpha , shape2=SYReturnParameters$beta, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  
  #ASYReturn
  FASYReturn <- pnorm(y[8], mean=0, sd=1)
  vrdat$ASYReturn[i] <- qbeta(FASYReturn, shape1= ASYReturnParameters$alpha , shape2=ASYReturnParameters$beta, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  
  #Recruitment
  FRecruitment <- pnorm(y[9], mean=0, sd=1)
  vrdat$recruitrate[i] <- qbeta(FRecruitment, shape1= recruitmentParameters$alpha , shape2=recruitmentParameters$beta, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  
  
  
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
  
  
}





p1 <- ggplot(data=vrdat, aes(SensSYClutch))+
  geom_density()+
  xlab("Sensitivity to SY Clutchsize")+
  ylab("Density")+
  ggthemes::theme_few(15)+
  xlim(0,4)



p2<- ggplot(data=vrdat, aes(SensASYClutch))+
  geom_density()+
  xlab("Sensitivity to ASY Clutchsize")+
  ylab("Density")+
  ggthemes::theme_few(15)+
  xlim(0,4)

p3 <- ggplot(data=vrdat, aes(SensSYNests))+
  geom_density()+
  xlab("Sensitivity to SY renests")+
  ylab("Density")+
  ggthemes::theme_few(15)+
  xlim(0,4)

p4 <- ggplot(data=vrdat, aes(SensASYNests))+
  geom_density()+
  xlab("Sensitivity to ASY renests")+
  ylab("Density")+
  ggthemes::theme_few(15)+
  xlim(0,4)

p5 <- ggplot(data=vrdat, aes(SensHatch))+
  geom_density()+
  xlab("Sensitivity to hatching")+
  ylab("Density")+
  ggthemes::theme_few(15)+
  xlim(0,4)

p6 <- ggplot(data=vrdat, aes(SensFledge))+
  geom_density()+
  xlab("Sensitivity to fledging")+
  ylab("Density")+
  ggthemes::theme_few(15)+
  xlim(0,4)

p7 <- ggplot(data=vrdat, aes(SensRecruit))+
  geom_density()+
  xlab("Sensitivity to recruitment")+
  ylab("Density")+
  ggthemes::theme_few(15)+
  xlim(0,4)

p8 <- ggplot(data=vrdat, aes(SensSYReturn))+
  geom_density()+
  xlab("Sensitivity to SY Return")+
  ylab("Denstity")+
  ggthemes::theme_few(15)+
  xlim(0,4)

p9 <- ggplot(data=vrdat, aes(SensASYReturn))+
  geom_density()+
  xlab("Sensitivity to ASY Return")+
  ylab("Density")+
  ggthemes::theme_few(15)+
  xlim(0,4)

library(cowplot)



png(filename = "~/Masters Thesis Project/Vital Rates Paper/Figures/Sensitivity Distributions.png", width=1000, height=1000)
plot_grid(p1, p2, p3, p4,p5, p6, p7, p8, p9, align="h")

dev.off()




p10 <- ggplot(data=vrdat, aes(ElasSYClutch))+
  geom_density()+
  xlab("Elasticity to SY Clutchsize")+
  ylab("Density")+
  ggthemes::theme_few(15)+
  xlim(0,1)



p11<- ggplot(data=vrdat, aes(ElasASYClutch))+
  geom_density()+
  xlab("Elasticity to ASY Clutchsize")+
  ylab("Density")+
  ggthemes::theme_few(15)+
  xlim(0,1)

p12 <- ggplot(data=vrdat, aes(ElasSYNests))+
  geom_density()+
  xlab("Elasticity to SY renests")+
  ylab("Density")+
  ggthemes::theme_few(15)+
  xlim(0,1)

p13 <- ggplot(data=vrdat, aes(ElasASYNests))+
  geom_density()+
  xlab("Elasticity to ASY renests")+
  ylab("Density")+
  ggthemes::theme_few(15)+
  xlim(0,1)

p14 <- ggplot(data=vrdat, aes(ElasHatch))+
  geom_density()+
  xlab("Elasticity to hatching")+
  ylab("Density")+
  ggthemes::theme_few(15)+
  xlim(0,20)

p15 <- ggplot(data=vrdat, aes(ElasFledge))+
  geom_density()+
  xlab("Elasticity to fledging")+
  ylab("Density")+
  ggthemes::theme_few(15)+
  xlim(0,20)

p16 <- ggplot(data=vrdat, aes(ElasRecruit))+
  geom_density()+
  xlab("Elasticity to recruitment")+
  ylab("Density")+
  ggthemes::theme_few(15)+
  xlim(0,1)

p17 <- ggplot(data=vrdat, aes(ElasSYReturn))+
  geom_density()+
  xlab("Elasticity to SY Return")+
  ylab("Denstity")+
  ggthemes::theme_few(15)+
  xlim(0,1)

p18 <- ggplot(data=vrdat, aes(ElasASYReturn))+
  geom_density()+
  xlab("Elasticity to ASY Return")+
  ylab("Density")+
  ggthemes::theme_few(15)+
  xlim(0,1)





png(filename = "~/Masters Thesis Project/Vital Rates Paper/Figures/Elasticity Distributions.png", width=1000, height=1000)
plot_grid(p10, p11, p12, p13,p14, p15, p16, p17, p18,  align="h")

dev.off()



ggplot(vrdat, aes(x=recruitrate, y=ElasRecruit, color=ASYReturn))+
  geom_point(alpha=0.7)
#if ou have lower adult return rates, recruitment rate isn't very elastic. If adult return is high, then recruitment is elastic


ggplot(vrdat, aes(x=hatchrate, y=ElasHatch, color=ASYReturn))+
  geom_point(alpha=0.7)+
  ylim(0,23)


#Figure out which vital rate has the highest sensitivity and which matrix
#element has the highest elasticity in each of the matrices.
MatrixElements <- c("SYFertility", "ASYFertility", "NestlingSurvival", "SYSurvival", "ASYSurvival")
for (i in 1:nrow(vrdat)){
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

vrdat2 <- reshape2::melt(vrdat[,29:31], measure.vars=colnames(vrdat)[29:31], variable.name="Rank", value.name="VitalRate")
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
  scale_fill_grey(start=0.3, end=0.8, labels=c("ASYReturnTot"= "ASY return", 
                                               "SYReturnTot"= "SY return",
                                               "RecruitTot"= "Recruitment", 
                                               "FledgeTot"="Fledge rate", 
                                               "HatchTot"="Hatch rate"), 
                  name="Vital Rate")+
  ylab("Proportion \n of simulations")+
  ggthemes::theme_few(base_size = 20)

png(filename = "~/Masters Thesis Project/Vital Rates Paper/Figures/Sensitivity Rankings Bar Plot_updated.png", 
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
  scale_fill_grey(start=0.3, end=0.8, labels=c("ASYSurvival"= "ASY return", 
                                               "NestlingSurvival"= "Nestling Survival"), 
                  name="Matrix Element")+
   scale_x_discrete("Elasticity ranking 1st", 
                    labels=c("1"=""))+
  ylab("Proportion \n of simulations")+
  ggthemes::theme_few(base_size = 20)
  #theme(legend.key.height = unit(2,"cm"))

png(filename = "~/Masters Thesis Project/Vital Rates Paper/Figures/Elasticity Rankings Bar Plot_updated.png", 
    width=350, height=500)
ElasRankPlot
dev.off()


#Make a combined plot for my extended abstract
png(filename = "~/Masters Thesis Project/Vital Rates Paper/Figure for conference.png", 
    width=800, height=500)

cowplot::plot_grid(SensRankPlot, ElasRankPlot, align="h", rel_widths = c(5, 3.3), labels=c("a", "b"), label_size = 23)

dev.off()
############Now lets do a deterministic sensitivity analysis on this data

#################IN REALITY THIS IS ACTUALLY NOT CALCULATING SENSITIVITY. WE'VE
#################DONE THAT ABOVE FOR EACH MATRIX. INSTEAD WHAT THIS IS DOING IS CALCULATING A
#################PROPORTION OF VARIATION IN POPULATION GROWTH RATE EXPLAINED BY EACH VITAL
#RATE--ALSO IMPORTANT
#Calculate "sensitivity" according to Taylor et al. 2012
#using the linear regression
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


######WE'VE STOPPED CALCULATING ELASTICITY LIKE THIS AND AM INSTEAD CALCULATING
#IT INDIVIDUALLY FOR EACH MATRIX. THIS IS QUITE DIFFERENT MATHEMATICALLY BUT
#TRUER TO WHAT ELASTICITY ACTUALLY IS
#Calculate Elasticity according to Taylor
#et al. 2012 using the linear regression of the log log transformation
#modASYReturn_elasticity <- lm(log(lambda)~ log(ASYReturn), data=vrdat)
#sumASY_e<- summary(modASYReturn_elasticity) #R^2=0.0001127
#
#modSYReturn_elasticity <- lm(log(lambda) ~ log(SYReturn), data=vrdat) sumSY_e
#<- summary(modSYReturn_elasticity) #R^2 = 0.0002805
#
#modClutchSizeSY_elasticity <- lm(log(lambda) ~ log(clutchSizeSY), data=vrdat)
#sumClutchSY_e<- summary(modClutchSizeSY_elasticity) # R^2 = 0.001324 so it's
#not that!
#
#modClutchSizeASY_elasticity <- lm(log(lambda) ~ log(clutchSizeASY), data=vrdat)
#sumClutchASY_e<- summary(modClutchSizeASY_elasticity) # R^2 = 0.01159 so it's
#not that!
#
#modAverageNestsSY_elasticity<- lm(log(lambda) ~ log(averageNestsSY),
#data=vrdat) sumAverageNestsSY_e<- summary(modAverageNestsSY_elasticity)
##Adjusted R^2 =0.0002755
#
#modAverageNestsASY_elasticity<- lm(log(lambda) ~ log(averageNestsASY),
#data=vrdat) sumAverageNestsASY_e<- summary(modAverageNestsASY_elasticity)
##Adjusted R^2 =0.01694, better but not a whole lot of the variation
#
#modEggtoAdult_elasticity <- lm(log(lambda)~log(eggRecruitRate), data=vrdat)
#sumEggtoAdult_e <- summary(modEggtoAdult_elasticity) #R^2 = 0.5758
#
#modhatch_elasticity <- lm(log(lambda)~log(hatchrate), data=vrdat) sumHatch_e <-
#summary(modhatch_elasticity) #R^2 =0.01803 so low elasticity
#
#modFledge_elasticity <- lm (log(lambda) ~log(fledgerate), data=vrdat)
#sumFledge_e <- summary(modFledge_elasticity) #Adjusted R^2 =0.08591 so another
#large portion of the variation,
#
#modRecruit_elasticity <- lm(log(lambda)~log(recruitrate), data=vrdat)
#sumRecruit_e<- summary(modRecruit_elasticity) #r^2=0.5271
#

#Put all of this analysis into a seperate dataframe 
SensitivityAnalysis<-as.data.frame(matrix(nrow=9, ncol=4))
colnames(SensitivityAnalysis)<-c("VitalRate", "PropVarExplained")

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

SensitivityAnalysis$PropVarExplained[1]<-  sumRecruit$adj.r.squared
SensitivityAnalysis$PropVarExplained[2]<-  sumSY$adj.r.squared

SensitivityAnalysis$PropVarExplained[3]<-  sumASY$adj.r.squared

SensitivityAnalysis$PropVarExplained[4]<-  sumAverageNestsSY$adj.r.squared

SensitivityAnalysis$PropVarExplained[5]<-  sumAverageNestsASY$adj.r.squared

SensitivityAnalysis$PropVarExplained[6]<-  sumClutchSY$adj.r.squared

SensitivityAnalysis$PropVarExplained[7]<-  sumClutchASY$adj.r.squared

SensitivityAnalysis$PropVarExplained[8]<- sumHatch$adj.r.squared

SensitivityAnalysis$PropVarExplained[9]<-  sumFledge$adj.r.squared



# SensitivityAnalysis$Elasticity[1]<- modRecruit_elasticity$coefficients[2]
# SensitivityAnalysis$R2_E[1] <- sumRecruit_e$adj.r.squared
# SensitivityAnalysis$Elasticity[2]<- modSYReturn_elasticity$coefficients[2]
# SensitivityAnalysis$R2_E[2] <- sumSY_e$adj.r.squared
# 
# SensitivityAnalysis$Elasticity[3]<- modASYReturn_elasticity$coefficients[2]
# SensitivityAnalysis$R2_E[3] <- sumASY_e$adj.r.squared
# 
# SensitivityAnalysis$Elasticity[4]<- modAverageNestsSY_elasticity$coefficients[2]
# SensitivityAnalysis$R2_E[4] <- sumAverageNestsSY_e$adj.r.squared
# 
# SensitivityAnalysis$Elasticity[5]<- modAverageNestsASY_elasticity$coefficients[2]
# SensitivityAnalysis$R2_E[5] <- sumAverageNestsASY_e$adj.r.squared
# 
# SensitivityAnalysis$Elasticity[6]<- modClutchSizeSY_elasticity$coefficients[2]
# SensitivityAnalysis$R2_E[6] <- sumClutchSY_e$adj.r.squared
# 
# SensitivityAnalysis$Elasticity[7]<- modClutchSizeASY_elasticity$coefficients[2]
# SensitivityAnalysis$R2_E[7] <- sumClutchASY_e$adj.r.squared
# 
# SensitivityAnalysis$Elasticity[8]<- modhatch_elasticity$coefficients[2]
# SensitivityAnalysis$R2_E[8] <- sumHatch_e$adj.r.squared
# 
# SensitivityAnalysis$Elasticity[9]<- modFledge_elasticity$coefficients[2]
# SensitivityAnalysis$R2_E[9] <- sumFledge_e$adj.r.squared

SensitivityAnalysis$Label <- c(7,8,9,1,2,3,4,5,6)
SensitivityAnalysis$Location[1:3]<- "Over wintering"

SensitivityAnalysis$Location[4:9]<- "Breeding"



# SensitivityPlot <- ggplot(SensitivityAnalysis, aes(x=Sensitivity, y=R2_S))+
#   geom_vline(xintercept = 0, linetype="dashed")+
#   xlab("Sensitivity")+
#   ylab(expression(paste("Adjusted R"^ 2)))+
#   ggrepel::geom_label_repel(data=SensitivityAnalysis %>% filter(Sensitivity>0.5), aes(label=VitalRate, color=Location), size=5,  box.padding = 0.5, segment.colour = "black", show.legend = F)+
#   ggrepel::geom_label_repel(data=SensitivityAnalysis %>% filter(Sensitivity<0.5), aes(label=VitalRate, color=Location), size=5, box.padding=2, force=5, segment.colour = "black", show.legend = F)+
#   geom_point(size=5, aes(color=Location))+ 
#   theme_classic(base_size = 20)+
#   theme(axis.title.y=element_text(angle=0, vjust=0.5))+
#   scale_color_manual(values=c("springgreen4", "steelblue4"))+
#   xlim(-.5, 2.5)
# SensitivityPlot
# 
# 
# ElasticityPlot <- ggplot(SensitivityAnalysis, aes(x=Elasticity, y=R2_E))+
#   xlab("Elasticity")+
#   ylab(expression(paste("Adjusted R"^ 2)))+
#   geom_vline(xintercept = 0, linetype="dashed")+
#   ggrepel::geom_label_repel(data=SensitivityAnalysis , aes(label=VitalRate, color=Location), size=5, min.segment.length = 0,  box.padding = 1.5, segment.color="black", show.legend = F)+
#   geom_point(size=5, aes(color=Location))+ 
#   theme_classic(base_size = 20)+
#   theme(axis.title.y=element_text(angle=0, vjust=0.5))+
#   xlim(-1, 0.8)+
#   scale_color_manual(values=c("springgreen4", "steelblue4"))
#   
# 
# ElasticityPlot
# setwd("~/Masters Thesis Project/Committee Meetings/Dec 15 2017")
# 
# png(filename = "Elasticity Plot.png", 
#      width=1000, height=600)
# ElasticityPlot
# dev.off()
# 
# png(filename = "Sensitivity Plot.png", 
#     width=1000, height=600)
# SensitivityPlot
# dev.off()




#Make the correlation plot for my paper-- take out recruitment!
library(corrplot)
library(RColorBrewer)
M <- cor(PopData[,-c(1)], use="pairwise.complete.obs")
colnames(M) <- c("SY renests", "ASY renests", "SY clutch size", "ASY clutch size", "Hatch rate", "Fledge rate", "SY return", "ASY return", "Recruitment")
rownames(M) <- c("SY renests", "ASY renests", "SY clutch size", "ASY clutch size", "Hatch rate", "Fledge rate", "SY return", "ASY return", "Recruitment")

CorMatrix <- corrplot(M, method = "circle", type="lower", order="alphabet", col=brewer.pal(n = 8, name = "RdYlBu"), tl.col = "black", tl.srt = 45)


png(filename = "~/Masters Thesis Project/Vital Rates Paper/Figures/Correlation Matrix Plot_updated.png", 
     width=500, height=500)
corrplot(M, method = "circle", type="lower", order="alphabet", col=brewer.pal(n = 8, name = "RdBu"), tl.col = "black", tl.srt = 45)

dev.off()

#Sensitivity analysis plot for my paper
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

png(filename = "~/Masters Thesis Project/Vital Rates Paper/Figures/Sensitivity Analysis_correlated2.png", width=1000, height=1000)
cowplot::plot_grid(SensitivityPlot2, ElasticityPlot2, align="v",nrow=2, ncol=1, labels=c("A", "B"))
dev.off()









########## Make some presentation quality graphs
ggplot(vrdat4, aes(x=Rank, y=Total, fill=VitalRate))+
  geom_bar(position=c("fill"), stat="identity")+
  scale_fill_brewer(palette= "Paired",
                   labels=c("RecruitTot"="Juvenile surivival", "ASYReturnTot"= "Older female survival", "SYReturnTot"="1-yr-old female survival", "FledgeTot"="Fledge rate", "HatchTot"="Hatch rate"))+
  scale_x_discrete(labels=c("VRMaxSens"="1st", "VR2MaxSens"="2nd", "VR3MaxSens"="3rd"))+
  
  labs(y="Proportion \n of simulations", fill="", x= "Sensitivity ranking by vital rate")+
  ggthemes::theme_few(base_size = 22)+
  theme(axis.title.y=element_text(angle=0, vjust=0.5)) 
  

ggsave(filename='~/Masters Thesis Project/BGRS symposium presentation/Sensitivity Rankings Plot.jpeg', width=9, height=6, units="in", device="jpeg")



ggplot(vrdat5 %>% filter (Total>0), aes(x=factor(1), y=Total, fill=MatrixElement))+
  geom_bar(position=c("fill"), stat="identity")+
  scale_fill_brewer(labels=c("ASYSurvival"= "Older female survival", 
                                               "NestlingSurvival"= "Nestling Survival"), palette="Accent")+
  scale_x_discrete("Elasticity ranking 1st", 
                   labels=c("1"=""))+
  ylab("Proportion \n of simulations")+
  ggthemes::theme_few(base_size = 20)
#theme(legend.key.height = unit(2,"cm"))