#How do our results change if we include the correlation structure for the KNOWN
#correlates, and don't let recruitment correlate with anything, because it's
#unclear whether recruitment does correlate with anything

setwd("~/Masters Thesis Project/TRES Data Analysis/RMark Preliminary Survival Analysis")
library(RMark)
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
  
  
nestlingMark <- readRDS( "Best Nestling MARK Results_capture rate unfixed.rda")

#The best model doesn't have recruitment vary by year so we'll put the same number in for all years, but use the se to create the distribution
PopData$recruitment<- nestlingMark$results$real$estimate[1]

#The best female model DOES let survival vary by year and by age
FemaleMark <- readRDS("Best Female MARK Results without fixing capture.rda")
summaryF <- summary(FemaleMark)
# I need to pull the number out on the diagonal for the SY birds
for(i in 1:42){
  PopData$SYReturn[i] <- summaryF$reals$Phi$`Group:age1`[[1]][i,i]
}
#remove 1975, 1976 and 2017 from both of these columns-- it's set to 1 for both the SY and ASY
#but so few birds were banded it's a terrible estimate
PopData$SYReturn[c(1, 43)]<- NA
PopData$ASYReturn <- c(summaryF$reals$Phi$`Group:age2`[[1]][1,], NA)
PopData$ASYReturn[c(1,  43)]<- NA

#Calculate parameters necessary later
recruitmentParameters <- estBetaParams(mu=nestlingMark$results$real$estimate[1], var=nestlingMark$results$real$se[1]^2)
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

#1. get the correlation matrix of all the vital rates that show yearly variation (ie not recruitment!!) 

c <- cor(PopData[,c(2:7, 9, 10)], method = "spearman", use="na.or.complete")
#Recruitment can't be correlated with anything here. so we will fill in those NAs with 0

#2. make a matrix (W) who's columns are all the possible right eigenvectors of C
eigC <- eigenUnOrdered(c)
W<- eigC$vectors

#3. make a matrix (D) who's has the eigen values along the diagonal and all else 0
#check to make sure that all eigen values are >=0 first, if less than they are likely artificial and should be changed to 0
#All of my eigen values are >0 so life is good and we can carry on
sqrt(abs(eigC$vectors))
D<- matrix(data=0, nrow=8, ncol=8)
for (i in 1:8){
  D[i,i]<- eigC$values[i]
  
}
#CHange D[1,1] to 0 because it's artificially negative (it's really just SUPER small)
D[1,1]<- 0


#This is where we need to start looping! and collecting data for the vital rates analysis!
vrdat <- as.data.frame(matrix(nrow=10000, ncol=12)) 
#Doing 10,000 because that's what Taylor et al 2012 did)
colnames(vrdat)<- c("hatchrate", "fledgerate", "recruitrate", "eggRecruitRate", "clutchSizeSY", "averageNestsSY","clutchSizeASY", "averageNestsASY", "SYReturn", "ASYReturn", "lambda", "GenTime" )
stages <- c("egg", "SY", "ASY")
A <- matrix(0, nrow=3, ncol=3, dimnames = list(stages, stages))
N0 <- c(0, 16 , 38)

for (i in 1:nrow(vrdat)){
  
  
  #3. Generate a vector (m) of uncorrelation standard normal variables. 
  m <- rnorm(n=8, mean=0, sd=1)
  
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
  
    #Phew I think we've now generated correlated data!! That's so great. 
  
  #Just need to generate the uncorrelated recruitment data!
  #Recruitment
  vrdat$recruitrate[i] <- rbeta(1, shape1=recruitmentParameters$alpha, shape2=recruitmentParameters$beta)

  
  
  vrdat$layrateSY[i]<- vrdat$clutchSizeSY[i]*vrdat$averageNestsASY[i]
  vrdat$layrateASY[i]<- vrdat$clutchSizeASY[i]*vrdat$averageNestsASY[i]
  
  vrdat$eggRecruitRate[i] <- vrdat$hatchrate[i] * vrdat$fledgerate[i] * vrdat$recruitrate[i]
  
  
  #Put those vital rates into the matrix
  #Rows= what stage the birds are in now
  #Columns=what stage the birds are going to
  A[1, 2] <- vrdat$layrateSY[i] 
  A[1, 3] <- vrdat$layrateASY [i]
  A[2, 1] <- vrdat$eggRecruitRate[i]
  A[3, 2] <- vrdat$SYReturn[i]
  A[3, 3] <- vrdat$ASYReturn[i]
  p<- pop.projection(A=A, n=N0, iterations=42)
  
  vrdat$lambda[i] <- p$lambda
  vrdat$GenTime[i] <- generation.time(A)
  
}


cor(vrdat[,1:8], method="spearman", use="complete.obs"
    )
cor(PopData, use="complete.obs", method="spearman")

#Fantastic! The correlation structure is right and the numbers look fine. I'm stoked!
plot(lambda~ASYReturn, data=vrdat)
plot(lambda~SYReturn, data=vrdat)
plot(lambda~layrateSY, data=vrdat)
plot(lambda~layrateASY, data=vrdat)

plot(lambda~eggRecruitRate, data=vrdat) #Ooooo we have a winner!!
#What within the eggRecruitRate is important? \
plot(lambda~recruitrate, data=vrdat) #recruitment looks super important but I'm still not sure I really like those numbers very much so I want to do some more thinking about that. 
plot(lambda~hatchrate, data=vrdat) #hatch rate looks much less important now. 
plot(lambda~fledgerate, data=vrdat)  #fledge rate is important


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




#Let's put all of those points onto a plot for visualizing just like Taylor et al 2012 did

#color things in by type of rate. 
#survival= dark green
#fertility= blue
#nestling success =purple

#Make a dataframe
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
library(ggrepel)

SensitivityPlot <- ggplot(SensitivityAnalysis, aes(x=Sensitivity, y=R2_S))+
  xlab("Sensitivity")+
  ylab(expression(paste("Adjusted R"^ 2)))+
  geom_point(size=3)+ 
  geom_label_repel(data=SensitivityAnalysis %>% filter(Sensitivity>0.5), aes(label=VitalRate), size=6, min.segment.length = 0, box.padding = 0.5)+
  geom_label_repel(data=SensitivityAnalysis %>% filter(Sensitivity<0.5), aes(label=VitalRate), size=6, min.segment.length = 0, box.padding = 1)+
  geom_vline(xintercept = 0, linetype="dashed")+
  theme_classic(base_size = 20)+
  theme(axis.title.y=element_text(angle=0, vjust=0.5))
  



png(filename = "~/TAships/Biol 302_ Populaiton Ecology/302 Guest Lecture/Sensitivity Plot_correlated.png", 
     width=800, height=500)
SensitivityPlot
dev.off()



#Calculate Elasticity according to Taylor et al. 2012 using the linear regression of the log log transformation

ElasticityPlot <- ggplot(SensitivityAnalysis, aes(x=Elasticity, y=R2_E))+
  xlab("Elasticity")+
  ylab(expression(paste("Adjusted R"^ 2)))+
  geom_point(size=3)+ 
  geom_label_repel(data=SensitivityAnalysis %>% filter(R2_E>0.3) ,aes(label=VitalRate), size=6, min.segment.length = 0,  box.padding = .5)+
  geom_label_repel(data=SensitivityAnalysis %>% filter(R2_E<0.3 & Elasticity<0), aes(label=VitalRate), size=6, min.segment.length = 0,  box.padding = .8)+
  geom_label_repel(data=SensitivityAnalysis %>% filter(R2_E<0.3 & Elasticity>0.2), aes(label=VitalRate), size=6, min.segment.length = 0,  box.padding = 1.5)+
  geom_label_repel(data=SensitivityAnalysis %>% filter(R2_E<0.3 & Elasticity>0 & Elasticity<0.2), aes(label=VitalRate), size=6, min.segment.length = 0,  box.padding = .8)+
  geom_vline(xintercept = 0, linetype="dashed")+
  theme_classic(base_size = 20)+
  theme(axis.title.y=element_text(angle=0, vjust=0.5))


png(filename = "~/TAships/Biol 302_ Populaiton Ecology/302 Guest Lecture/Elasticity Plot_correlated.png", 
     width=800, height=500)
ElasticityPlot
dev.off()


#Let's make

SensitivityPlot2 <- ggplot(SensitivityAnalysis, aes(x=Sensitivity, y=R2_S))+
  xlab("Sensitivity")+
  ylab(expression(paste("Adjusted R"^ 2)))+
  geom_point(size=3)+ 
  geom_label_repel(data=SensitivityAnalysis %>% filter (Sensitivity<0.5), aes(label=VitalRate), size=3, min.segment.length = 0, max.iter = 3000, box.padding = 1.2)+
  geom_label_repel(data=SensitivityAnalysis %>% filter (Sensitivity>0.5), aes(label=VitalRate), size=3, min.segment.length = 0, max.iter = 3000, box.padding = 0.5)+
  geom_vline(xintercept = 0, linetype="dashed")+
  theme_classic(base_size = 20)

ElasticityPlot2 <- ggplot(SensitivityAnalysis, aes(x=Elasticity, y=R2_E))+
  xlab("Elasticity")+
  ylab(expression(paste("Adjusted R"^ 2)))+
  geom_point(size=3)+ 
  geom_label_repel(aes(label=VitalRate), size=3, min.segment.length = 0, max.iter = 3000, box.padding = .5)+
  geom_vline(xintercept = 0, linetype="dashed")+
  theme_classic(base_size = 20)

library(cowplot)

CorrSensPlot <- plot_grid(SensitivityPlot2, ElasticityPlot2, align="v",nrow=2, labels=c("A", "B"))


png("~/Masters Thesis Project/Vital Rates Paper/Sensitivity Analysis Plot with correlation.png", height=800, width=500)
CorrSensPlot
dev.off()


#WHat is the mean and sd on lambda?
hist(vrdat$lambda)
median(vrdat$lambda) 
#median population growth rate is 0.801, therefore it must be true that to have a
#stable population, we would need to have immigration supplying 0.2 of the
#lambda.
mean(vrdat$lambda)
x<- vrdat$lambda
se <- sd(x)/sqrt(length(x))



#What is the mean and sd on generation time?

hist(vrdat$GenTime)
mean(vrdat$GenTime)
median(vrdat$GenTime)
sd(vrdat$GenTime)


#Check in to see if rankings of elasticity and sensitivity change when you do the correlations

RankSensitivity <- SensitivityAnalysis$VitalRate[order(SensitivityAnalysis$Sensitivity, decreasing=T)]
#[1] "Recruitment"     "SY Return"       "ASY Return"      "Fledge Rate"     "Hatch Rate"      "ASY Nests"       "ASY Clutch Size" "SY Clutch Size" [9] "SY Nests

RankElasticity <- SensitivityAnalysis$VitalRate[order(SensitivityAnalysis$Elasticity, decreasing=T)]
#[1] "ASY Clutch Size" "ASY Return"      "Fledge Rate"     "SY Return"       "Recruitment"     "Hatch Rate"      "ASY Nests"       "SY Nests" [9] "SY Clutch Size" 



SensitivityAnalysis[order(SensitivityAnalysis$Elasticity, decreasing=T),]


#let's make a correlation plot
library(corrplot)
library(RColorBrewer)
colnames(c)<- c("SY Nests", "ASY Nests", "SY Clutch Size", "ASY Clutch Size", "Hatch Rate", "Fledge Rate", "SY Return", "ASY Return")
rownames(c) <- colnames(c)


png("~/Masters Thesis Project/Vital Rates Paper/Correlation Matrix.png", height=500, width=600)
corrplot(c, method = "circle", type="upper", col=brewer.pal(n = 8, name = "RdYlBu"), tl.col="black")
dev.off()
