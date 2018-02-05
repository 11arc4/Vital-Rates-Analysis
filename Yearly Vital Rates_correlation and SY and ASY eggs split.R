#Let's improve on the yearly vital rates analysis with correlations by seperating things out into SY eggs and ASY eggs


setwd("~/Masters Thesis Project/TRES Data Analysis/RMark Preliminary Survival Analysis")
library(RMark)
library(dplyr)
library(popbio)
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
PopData<- read.csv(file= "file:///C:/Users/Amelia/Documents/Masters Thesis Project/TRES Data Analysis/Matrix Pop Estimates/Yearly Vital Rates.csv", na.strings=c(""), as.is=T)
nestlingMark <- readRDS( "Best MARK Results for Vital Rates Nestlings Analysis with discrete time variable.rda")
summary <- summary(nestlingMark)




 PopData$recruitment<- c(summary$reals$Phi[[1]][[1]][1,], NA)

#there are a number of years that weren't parameterized properly-- you can tell because it's set to 1
#remove thos THIS IS ONE OF THE PROBLEMS IN OUR DATA THAT I WOULD PREFER TO HAVE DEALT WITH
PopData$recruitment[which(PopData$recruitment>0.98)]<- NA

FemaleMark <- readRDS("Best MARK Results for Vital Rates Adults Analysis.rda")
summaryF <- summary(FemaleMark)
# I need to pull the number out on the diagonal for the SY birds
for(i in 1:42){
  PopData$SYReturn[i] <- summaryF$reals$Phi$`Group:age1`[[1]][i,i]
}
#remove 1975 from both of these columns-- it's set to 1 for both the SY and ASY
#but so few birds were banded it's a terrible estimate. Also take out 2017 because we don't know whether 2017 returned or not
PopData$SYReturn[1]<- NA
PopData$SYReturn[nrow(PopData)]<- NA

PopData$ASYReturn <- c(summaryF$reals$Phi$`Group:age2`[[1]][1,], NA)
PopData$ASYReturn[1]<- NA
PopData$ASYReturn[nrow(PopData)]<- NA


#Calculate parameters necessary later
recruitmentParameters <- estBetaParams(mu=mean(PopData$recruitment, na.rm=T), var=var(PopData$recruitment, na.rm=T))
SYReturnParameters<- estBetaParams(mu=mean(PopData$SYReturn, na.rm=T ), var=var(PopData$SYReturn, na.rm=T))
ASYReturnParameters<- estBetaParams(mu=mean(PopData$ASYReturn, na.rm=T), var=var(PopData$ASYReturn, na.rm=T))
fledgeParametersSY <- estBetaParams(mu=mean(PopData$fledgeRateSY, na.rm=T), var=var(PopData$fledgeRateSY, na.rm=T))
fledgeParametersASY <- estBetaParams(mu=mean(PopData$fledgeRateASY), var=var(PopData$fledgeRateASY))
hatchrateParametersSY <- estBetaParams(mu=mean(PopData$hatchRateSY, na.rm=T), var=var(PopData$hatchRateSY, na.rm=T))
hatchrateParametersASY <- estBetaParams(mu=mean(PopData$hatchRateASY), var=var(PopData$hatchRateASY))
averageNestParametersSY<- estBetaParams(mu=mean(PopData$averageNestsSY-1, na.rm=T), var=var(PopData$averageNestsSY-1, na.rm=T))
averageNestParametersASY<- estBetaParams(mu=mean(PopData$averageNestsASY-1), var=var(PopData$averageNestsASY-1))

#clutch size works using a normal distribution but only if we use clutch size averaged across the ages and that's not useful!
shapiro.test(PopData$clutchSizeSY)
shapiro.test(PopData$clutchSizeASY)
#I will instead use a stretched beta distribution
#max clutch ever seen was 12, min was 1, so we'll go with those as the max and mmins for the stretch beta


betaClutchSizeParametersSY <- estBetaParams(mu=(mean(PopData$clutchSizeSY, na.rm=T)-min(PopData$clutchSizeSY, na.rm=T))/(max(PopData$clutchSizeSY, na.rm=T)-min(PopData$clutchSizeSY, na.rm=T)), 
                                            var=(var(PopData$clutchSizeSY, na.rm=T)* (1/ (max(PopData$clutchSizeSY, na.rm=T)-min(PopData$clutchSizeSY, na.rm=T)))^2))

betaClutchSizeParametersASY <- estBetaParams(mu=(mean(PopData$clutchSizeASY, na.rm=T)-min(PopData$clutchSizeASY, na.rm=T))/(max(PopData$clutchSizeASY, na.rm=T)-min(PopData$clutchSizeASY, na.rm=T)), 
                                             var=(var(PopData$clutchSizeASY, na.rm=T)* (1/(max(PopData$clutchSizeASY, na.rm=T)-min(PopData$clutchSizeASY, na.rm=T)))^2))


#Follow directions in Morris and Doak 2002 Quantitative Conservation Biology pg 284 

#1. get the correlation matrix of all the vital rates (not inluding the year!)
c <- cor(PopData[,c(3, 4, 6, 7, 9, 10, 12, 13, 14, 15, 16)], method = "spearman", use="na.or.complete")

#2. make a matrix (W) who's columns are all the possible right eigenvectors of C
eigC <- eigenUnOrdered(c)
W<- eigC$vectors

#3. make a matrix (D) who's has the eigen values along the diagonal and all else 0
#check to make sure that all eigen values are >=0 first, if less than they are likely artificial and should be changed to 0
#All of my eigen values are >0 so life is good and we can carry on
sqrt(abs(eigC$vectors))
D<- matrix(data=0, nrow=11, ncol=11)
for (i in 1:11){
  D[i,i]<- eigC$values[i]
  
}


#This is where we need to start looping! and collecting data for the vital rates analysis!
vrdat <- as.data.frame(matrix(nrow=10000, ncol=16)) 
#Doing 10,000 because that's what Taylor et al 2012 did)
colnames(vrdat)<- c("hatchrateSY",
                    "hatchrateASY",
                    "fledgerateSY", 
                    "fledgerateASY",
                    "recruitrate", 
                    "eggRecruitRateSY",
                    "eggRecruitRateASY",
                    "clutchSizeSY",
                    "clutchSizeASY",
                    "averageNestsSY" ,
                    "averageNestsASY",
                    "SYReturn", 
                    "ASYReturn", 
                    "layrateSY", 
                    "layrateASY", 
                    "lambda" )

stages <- c("eggSY","eggASY", "SY", "ASY")
A <- matrix(0, nrow=4, ncol=4, dimnames = list(stages, stages))
N0 <- c(0,0, 16 , 38)

for (i in 1:nrow(vrdat)){
  
  
  #3. Generate a vector (m) of uncorrelation standard normal variables. 
  m <- rnorm(n=11, mean=0, sd=1)
  
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
  
  vrdat$clutchSizeSY[i]<- 0.5 * betavalSY*(max(PopData$clutchSizeSY, na.rm=T)-min(PopData$clutchSizeSY, na.rm=T))+min(PopData$clutchSizeSY, na.rm=T) #as per morris and doak pg 283
  
  FClutchSizeASY <- pnorm(y[4], mean=0, sd=1)
  betavalASY <- qbeta(FClutchSizeASY, shape1= betaClutchSizeParametersASY$alpha , shape2=betaClutchSizeParametersASY$beta, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  
  vrdat$clutchSizeASY[i]<- 0.5* betavalASY*(max(PopData$clutchSizeASY, na.rm=T)-min(PopData$clutchSizeASY, na.rm=T))+min(PopData$clutchSizeASY, na.rm=T) #as per morris and doak pg 283
  
   #as per morris and doak pg 283
  
  #hatchRate 
  FHatchRateSY <- pnorm(y[5], mean=0, sd=1)
  vrdat$hatchrateSY[i] <- qbeta(FHatchRateSY, shape1= hatchrateParametersSY$alpha , shape2=hatchrateParametersSY$beta, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  FHatchRateASY <- pnorm(y[6], mean=0, sd=1)
  vrdat$hatchrateASY[i] <- qbeta(FHatchRateASY, shape1= hatchrateParametersASY$alpha , shape2=hatchrateParametersASY$beta, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  
  
  #fledgeRate 
  FFledgeRateSY <- pnorm(y[7], mean=0, sd=1)
  vrdat$fledgerateSY[i] <- qbeta(FFledgeRateSY, shape1= fledgeParametersSY$alpha , shape2=fledgeParametersSY$beta, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  FFledgeRateASY <- pnorm(y[8], mean=0, sd=1)
  vrdat$fledgerateASY[i] <- qbeta(FFledgeRateASY, shape1= fledgeParametersASY$alpha , shape2=fledgeParametersASY$beta, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  
  #SYReturn
  FSYReturn <- pnorm(y[10], mean=0, sd=1)
  vrdat$SYReturn[i] <- qbeta(FSYReturn, shape1= SYReturnParameters$alpha , shape2=SYReturnParameters$beta, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  
  #ASYReturn
  FASYReturn <- pnorm(y[11], mean=0, sd=1)
  vrdat$ASYReturn[i] <- qbeta(FASYReturn, shape1= ASYReturnParameters$alpha , shape2=ASYReturnParameters$beta, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  
  #Recruitment
  FRecruitment <- pnorm(y[9], mean=0, sd=1)
  vrdat$recruitrate[i] <- qbeta(FRecruitment, shape1= recruitmentParameters$alpha , shape2=recruitmentParameters$beta, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  #Phew I think we've now generated correlated data!! That's so great. 
  
  
  vrdat$layrateSY[i]<- vrdat$clutchSizeSY[i]*vrdat$averageNestsSY[i]
  vrdat$layrateASY[i]<- vrdat$clutchSizeASY[i]*vrdat$averageNestsASY[i]
  
  vrdat$eggRecruitRateSY[i] <- vrdat$hatchrateSY[i] * vrdat$fledgerateSY[i] * vrdat$recruitrate[i]
  
  vrdat$eggRecruitRateASY[i] <- vrdat$hatchrateASY[i] * vrdat$fledgerateASY[i] * vrdat$recruitrate[i]
  
  #Put those vital rates into the matrix
  #Rows= what stage the birds are in now
  #Columns=what stage the birds are going to
  A[1, 3] <- vrdat$layrateSY[i] 
  A[1, 4] <- vrdat$layrateASY [i]
  A[3, 1] <- vrdat$eggRecruitRateSY[i]
  A[3, 2] <- vrdat$eggRecruitRateASY[i]
  A[4, 3] <- vrdat$SYReturn[i]
  A[4, 4] <- vrdat$ASYReturn[i]
  p<- pop.projection(A=A, n=N0, iterations=43)
  
  vrdat$lambda[i] <- p$lambda
  
}


cor(vrdat[,1:8], method="spearman", use="complete.obs"
)
cor(PopData, use="complete.obs", method="spearman")

#Fantastic! The correlation structure is right and the numbers look fine. I'm stoked!
plot(lambda~ASYReturn, data=vrdat)
plot(lambda~SYReturn, data=vrdat)
plot(lambda~layrateSY, data=vrdat)
plot(lambda~layrateASY, data=vrdat)

plot(lambda~eggRecruitRateSY, data=vrdat) 
plot(lambda~eggRecruitRateASY, data=vrdat)

#What within the eggRecruitRate is important? \
plot(lambda~recruitrate, data=vrdat) 
plot(lambda~hatchrateSY, data=vrdat) 
plot(lambda~fledgerateSY, data=vrdat) 
plot(lambda~hatchrateASY, data=vrdat)  
plot(lambda~fledgerateASY, data=vrdat)  


############Now lets do a deterministic sensitivity analysis on this data

#Calculate sensitivity according to Taylor et al. 2012 using the linear regression
modASYReturn <- lm(lambda~ ASYReturn, data=vrdat)
sumASY <- summary(modASYReturn) #Adjusted R^2 = 0.005868 so quite unimportant

modSYReturn <- lm(lambda~ SYReturn, data=vrdat)
sumSY <- summary(modSYReturn) # Adjusted R^2 = 0.004281 so even less important

modClutchSizeSY <- lm(lambda ~ clutchSizeSY, data=vrdat)
sumClutchSY<- summary(modClutchSizeSY) #Adjusted R^2 = 0.007057

modClutchSizeASY <- lm(lambda ~ clutchSizeASY, data=vrdat)
sumClutchASY<- summary(modClutchSizeASY) #Adjusted R^2 = 0.02036 


modAverageNestsSY<- lm(lambda ~ averageNestsSY, data=vrdat) 
sumAverageNestsSY <- summary(modAverageNestsSY) #Adjusted R^2 =-0.0009415

modAverageNestsASY<- lm(lambda ~ averageNestsASY, data=vrdat) 
sumAverageNestsASY <- summary(modAverageNestsASY) #Adjusted R^2 =0.00283

modEggtoAdultSY <- lm(lambda~eggRecruitRateSY, data=vrdat)
sumEggtoAdultSY<- summary(modEggtoAdultSY) #Adjusted R^2 = 0.8306 so explains almost all the variation!
#What within this explains variation? 

modHatchSY <- lm(lambda ~ hatchrateSY, data=vrdat)
sumHatchSY<- summary(modHatchSY) # Adjussted R^2 =0.03876

modFledgeSY <- lm (lambda ~fledgerateSY, data=vrdat)
sumFledgeSY <- summary(modFledgeSY) #Adjusted R^2 =0.1319 so another large portion of the variation, 

modRecruit <- lm(lambda~recruitrate, data=vrdat)
sumRecruit <- summary(modRecruit) #R^2= 0.5667


modEggtoAdultSY <- lm(lambda~eggRecruitRateSY, data=vrdat)
sumEggtoAdultSY<- summary(modEggtoAdultSY) #Adjusted R^2 = 0.8306 so explains almost all the variation!
#What within this explains variation? 

modHatchASY <- lm(lambda ~ hatchrateASY, data=vrdat)
sumHatchASY<- summary(modHatchASY) # Adjussted R^2 =0.000592

modFledgeASY <- lm (lambda ~fledgerateASY, data=vrdat)
sumFledgeASY <- summary(modFledgeASY) #Adjusted R^2 =0.08612 so another large portion of the variation, 

#Let's put all of those points onto a plot for visualizing just like Taylor et al 2012 did
ggplot()+
  xlab("Sensitivity")+
  ylab(expression(paste("R"^ 2)))+
  geom_point(aes(x=modRecruit$coefficients[2], y=sumRecruit$r.squared))+ #Add Recruitment
  annotate("text",x=modRecruit$coefficients[2]-0.08, y=sumRecruit$r.squared, label="Recruitment")+
  geom_point(aes(x=modASYReturn$coefficients[2], y=sumASY$r.squared))+ #add ASY return
  annotate("text", x=modASYReturn$coefficients[2]-0.07, y=sumASY$r.squared, label="ASY Return")+
  geom_point(aes(x=modSYReturn$coefficients[2], y=sumSY$r.squared))+ #add SY return
  annotate("text",x=modSYReturn$coefficients[2]-0.07, y=sumSY$r.squared, label="SY Return" )+
  geom_point(aes(x=modAverageNestsSY$coefficients[2], y=sumAverageNestsSY$r.squared)) + #add average nests
  annotate("text",x=modAverageNestsSY$coefficients[2]+0.1, y=sumAverageNestsSY$r.squared, label="SY Nest Attempts" )+
  geom_point(aes(x=modAverageNestsASY$coefficients[2], y=sumAverageNestsASY$r.squared)) + #add average nests
  annotate("text",x=modAverageNestsASY$coefficients[2]-0.11, y=sumAverageNestsASY$r.squared, label="ASY Nest Attempts" )+
  geom_point(aes(x=modClutchSizeSY$coefficients[2], y=sumClutchSY$r.squared))+ #add clutch size
  annotate("text", x=modClutchSizeSY$coefficients[2]-0.08, y=sumClutchSY$r.squared, label="SY Clutch Size")+
  geom_point(aes(x=modClutchSizeASY$coefficients[2], y=sumClutchASY$r.squared))+ #add clutch size
  annotate("text", x=modClutchSizeASY$coefficients[2]-0.09, y=sumClutchASY$r.squared, label="ASY Clutch Size")+
  geom_point(aes(x=modHatchSY$coefficients[2], y=sumHatchSY$r.squared))+ #add Hatch size
  annotate("text", x=modHatchSY$coefficients[2]-0.085, y=sumHatchSY$r.squared, label="SY Hatch Rate")+
  geom_point(aes(x=modHatchASY$coefficients[2], y=sumHatchASY$r.squared))+ #add Hatch size
  annotate("text", x=modHatchASY$coefficients[2]+0.09, y=sumHatchASY$r.squared, label="ASY Hatch Rate")+
  geom_point(aes(x=modFledgeSY$coefficients[2], y=sumFledgeSY$r.squared))+ #add fledge rate
  annotate("text",x=modFledgeSY$coefficients[2]-0.085, y=sumFledgeSY$r.squared, label="SY Fledge Rate" )+
  
  geom_point(aes(x=modFledgeASY$coefficients[2], y=sumFledgeASY$r.squared))+ #add fledge rate
  annotate("text",x=modFledgeASY$coefficients[2]-0.09, y=sumFledgeASY$r.squared, label="ASYFledge Rate" )+
  theme_classic()+
  scale_y_log10(breaks=c(0.01, 0.05, 0.1, 0.2,0.3, 0.4,0.5, 0.7))+ #I need a much better log scale
  theme(axis.title.y = element_text(angle=0, vjust=0.5) )
  
#Calculate Elasticity according to Taylor et al. 2012 using the linear regression of the log log transformation
modASYReturn_elasticity <- lm(log(lambda)~ log(ASYReturn), data=vrdat)
sumASY_e<- summary(modASYReturn_elasticity) #R^2=0.01992,

modSYReturn_elasticity <- lm(log(lambda) ~ log(SYReturn), data=vrdat)
sumSY_e <- summary(modSYReturn_elasticity) #R^2 = 0.02011

modClutchSizeSY_elasticity <- lm(log(lambda) ~ log(clutchSizeSY), data=vrdat)
sumClutchSY_e<- summary(modClutchSizeSY_elasticity) # R^2 = 0.001324 so it's not that!

modClutchSizeASY_elasticity <- lm(log(lambda) ~ log(clutchSizeASY), data=vrdat)
sumClutchASY_e<- summary(modClutchSizeASY_elasticity) # R^2 = 0.01159 so it's not that!

modAverageNestsSY_elasticity<- lm(log(lambda) ~ log(averageNestsSY), data=vrdat) 
sumAverageNestsSY_e<- summary(modAverageNestsSY_elasticity) #Adjusted R^2 =0.0002755

modAverageNestsASY_elasticity<- lm(log(lambda) ~ log(averageNestsASY), data=vrdat) 
sumAverageNestsASY_e<- summary(modAverageNestsASY_elasticity) #Adjusted R^2 =0.01694, better but not a whole lot of the variation


modEggtoAdultSY_elasticity <- lm(log(lambda)~log(eggRecruitRateSY), data=vrdat)
sumEggtoAdultSY_e <- summary(modEggtoAdultSY_elasticity) #R^2 = 0.5328

modEggtoAdultASY_elasticity <- lm(log(lambda)~log(eggRecruitRateASY), data=vrdat)
sumEggtoAdultASY_e <- summary(modEggtoAdultASY_elasticity) #R^2 = 0.4564

modhatchSY_elasticity <- lm(log(lambda)~log(hatchrateSY), data=vrdat)
sumHatchSY_e <- summary(modhatchSY_elasticity) #R^2 =0.05163 so low elasticity

modhatchASY_elasticity <- lm(log(lambda)~log(hatchrateASY), data=vrdat)
sumHatchASY_e <- summary(modhatchASY_elasticity) #R^2 =0.00109 so low elasticity

modFledgeSY_elasticity <- lm (log(lambda) ~log(fledgerateSY), data=vrdat)
sumFledgeSY_e <- summary(modFledgeSY_elasticity) #Adjusted R^2 =0.1621 so another large portion of the variation, 

modFledgeASY_elasticity <- lm (log(lambda) ~log(fledgerateASY), data=vrdat)
sumFledgeASY_e <- summary(modFledgeASY_elasticity) #Adjusted R^2 =0.05954 so another large portion of the variation, 

modRecruit_elasticity <- lm(log(lambda)~log(recruitrate), data=vrdat)
sumRecruit_e<- summary(modRecruit_elasticity) #r^2=0.4087



ggplot()+
  xlab("Elasticity")+
  ylab(expression(paste("R"^ 2)))+
  geom_point(aes(x=modRecruit_elasticity$coefficients[2], y=sumRecruit_e$r.squared))+ #Add Recruitment
  annotate("text",x=modRecruit_elasticity$coefficients[2], y=sumRecruit_e$r.squared-.01, label="Recruitment")+
  geom_point(aes(x=modASYReturn_elasticity$coefficients[2], y=sumASY_e$r.squared))+ #add ASY return
  annotate("text", x=modASYReturn_elasticity$coefficients[2], y=sumASY_e$r.squared+0.01, label="ASY Return")+
  geom_point(aes(x=modSYReturn_elasticity$coefficients[2], y=sumSY_e$r.squared))+ #add SY return
  annotate("text",x=modSYReturn_elasticity$coefficients[2], y=sumSY_e$r.squared+0.01, label="SY Return" )+
  geom_point(aes(x=modAverageNestsSY_elasticity$coefficients[2], y=sumAverageNestsSY_e$r.squared)) + #add average nests
  annotate("text",x=modAverageNestsSY_elasticity$coefficients[2]+0.03, y=sumAverageNestsSY_e$r.squared-.01, label="SY Nest Attempts" )+
  geom_point(aes(x=modAverageNestsASY_elasticity$coefficients[2], y=sumAverageNestsASY_e$r.squared)) + #add average nests
  annotate("text",x=modAverageNestsASY_elasticity$coefficients[2], y=sumAverageNestsASY_e$r.squared+.01, label="ASY Nest Attempts" )+
  geom_point(aes(x=modClutchSizeSY_elasticity$coefficients[2], y=sumClutchSY_e$r.squared))+ #add clutch size
  annotate("text", x=modClutchSizeSY_elasticity$coefficients[2], y=sumClutchSY_e$r.squared+.01, label="SY Clutch Size")+
  geom_point(aes(x=modClutchSizeASY_elasticity$coefficients[2], y=sumClutchASY_e$r.squared))+ #add clutch size
  annotate("text", x=modClutchSizeASY_elasticity$coefficients[2], y=sumClutchASY_e$r.squared-.01, label="ASY Clutch Size")+
  geom_point(aes(x=modhatchSY_elasticity$coefficients[2], y=sumHatchSY_e$r.squared))+ #add Hatch size
  annotate("text", x=modhatchSY_elasticity$coefficients[2], y=sumHatchSY_e$r.squared+.01, label="SY Hatch Rate")+
  
  geom_point(aes(x=modhatchASY_elasticity$coefficients[2], y=sumHatchASY_e$r.squared))+ #add Hatch size
  annotate("text", x=modhatchASY_elasticity$coefficients[2]-0.03, y=sumHatchASY_e$r.squared+.01, label="ASY Hatch Rate")+
  geom_point(aes(x=modFledgeSY_elasticity$coefficients[2], y=sumFledgeSY_e$r.squared))+ #add fledge rate
  annotate("text",x=modFledgeSY_elasticity$coefficients[2], y=sumFledgeSY_e$r.squared+.01, label="SY Fledge Rate" )+
  
  geom_point(aes(x=modFledgeASY_elasticity$coefficients[2], y=sumFledgeASY_e$r.squared))+ #add fledge rate
  annotate("text",x=modFledgeASY_elasticity$coefficients[2], y=sumFledgeASY_e$r.squared+.01, label="ASY Fledge Rate" )+
  theme_classic()+
  theme(axis.title.y = element_text(angle=0, vjust=0.5) )+
  scale_x_continuous(limits=c(-0.22, 0.7))


