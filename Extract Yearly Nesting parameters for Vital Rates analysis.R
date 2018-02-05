#I need to get the population level statistics for each year of our data. 
library(dplyr)
library(popbio)
library(popdemo)

#To do that I need to pull out all the data for each of the nests

firsteggdate<- rep(NA, length(as.list(globalData$nests)))
clutch<- rep(NA, length(as.list(globalData$nests)))
hatch <- rep(NA, length(as.list(globalData$nests)))
fledge <- rep(NA, length(as.list(globalData$nests)))
FAge <- rep(NA, length(as.list(globalData$nests)))
FBand <- rep(NA, length(as.list(globalData$nests)))
year <- rep(NA, length(as.list(globalData$nests)))
parameters <- data.frame(firsteggdate, clutch, hatch, fledge, FBand, FAge, year)

i=0
for(nest in as.list(globalData$nests)){
  i=i+1
  parameters$clutch[i] <- nest$clutchSize
  parameters$firsteggdate[i]<- nest$firstEggDate
  parameters$hatch[i] <- nest$hatchSize
  parameters$fledge[i] <- nest$fledgeSize
  parameters$year[i]<- nest$year
  parameters$renestStatus[i] <- nest$renestStatus
  
  parameters$experiment[i] <- nest$experiment
  
  if (!is.na(nest$femaleID$m_key)){
    #message("Known female", nest$femaleID$m_key)
    
    bird <- get(nest$femaleID$m_key, globalData$birds)
    for (year in bird$yearsSeen$as.list()){
      if(nest$year==year$year){
        parameters$FAge[i]<- year$age
        parameters$FBand[i]<- bird$bandID
      }
    }
  }
}
parameters$FAge[which(parameters$FAge=="HY")]<- NA  
parameters$FAge2 <- rep(NA, nrow(parameters))
parameters$FAge2[which(parameters$FAge!="SY" & !is.na(parameters$FAge))]<- "ASY"
parameters$FAge2[which(parameters$FAge=="SY")]<- "SY"

#Also need to pull out how many nesting attempts are average for each year
nestsinyear <- rep(NA, length(as.list(globalData$birds)))
year <- rep(NA, length(as.list(globalData$birds)))
sex <- rep(NA, length(as.list(globalData$birds)))
age <- rep(NA, length(as.list(globalData$birds)))
birdband <- rep(NA, length(as.list(globalData$birds)))
reprod<- data.frame(nestsinyear, year, sex, age, birdband)
a=0
for(bird in as.list(globalData$birds)){
  #if a bird was seen as a nestling BUT was seen in multiple years or the bird wasn't a nestling
  if((!is.na(bird$hatchnest$m_key) & bird$yearsSeen$length>1)| is.na(bird$hatchnest$m_key)){
    for(year in bird$yearsSeen$as.list()){
      if(year$nest$length>0){
        a=a+1
        if(is.na(year$hatchNest$m_key)){
          reprod$nestsinyear[a]<- year$nest$length
          reprod$age[a] <- year$age
          reprod$sex[a] <- bird$sex
          reprod$birdband[a] <- bird$bandID
          reprod$year[a]<- year$year
        }
        
      }
      
    }
  }
}
reprod <- reprod[1:a,]


#Fill in a population level data sheet for each year
year<- seq(1975, 2017, 1)
PopData <- data.frame(year)

a <- 0
for (y in 1975:2017){
  a<- a+1
  #fill in the average number of renests (for females) known in that year
  PopData$averageNests[a] <- 
    sum(reprod$nestsinyear[ which(reprod$year==y & reprod$sex=="F")]) / length(which(reprod$year==y & reprod$sex=="F")) 
  PopData$averageNestsSY[a] <- 
    sum(reprod$nestsinyear[ which(reprod$year==y & reprod$sex=="F" & reprod$age=="SY")]) / length(which(reprod$year==y & reprod$sex=="F" & reprod$age=="SY")) 
  PopData$averageNestsASY[a] <- 
    sum(reprod$nestsinyear[ which(reprod$year==y & reprod$sex=="F"& reprod$age!="SY" & reprod$age!="AHY")]) / length(which(reprod$year==y & reprod$sex=="F" & reprod$age!="SY" & reprod$age!="AHY")) 
  
  
  #Fill in average clutch size that year
  clutchnests <- parameters %>% filter (!is.na(clutch) & year==y)
  clutchnestsSY <- parameters %>% filter (!is.na(clutch) & year==y & FAge=="SY")
  clutchnestsASY <- parameters %>% filter (!is.na(clutch) & year==y & FAge!="SY" & FAge!="AHY")
  
  PopData$clutchSize[a] <-  sum(clutchnests$clutch) / nrow(clutchnests)
  PopData$clutchSizeSY[a] <-  sum(clutchnestsSY$clutch) / nrow(clutchnestsSY)
  PopData$clutchSizeASY[a] <-  sum(clutchnestsASY$clutch) / nrow(clutchnestsASY)
  
  
  #Fill in average hatchrate
  hatchnests <- parameters %>% filter (clutch>0 & !is.na(hatch) & year==y)
  hatchnestsSY <- parameters %>% filter (clutch>0 & !is.na(hatch) & year==y & FAge=="SY")
  hatchnestsASY <- parameters %>% filter (clutch>0 & !is.na(hatch) & year==y & FAge!="SY" & FAge!="AHY")
  
  
  PopData$hatchRate[a] <- sum(hatchnests$hatch) / sum (hatchnests$clutch)
  PopData$hatchRateSY[a] <- sum(hatchnestsSY$hatch) / sum (hatchnestsSY$clutch)
  PopData$hatchRateASY[a] <- sum(hatchnestsASY$hatch) / sum (hatchnestsASY$clutch)
  
  #fill in average fledge size
  fledgenests <- parameters %>% filter (!is.na(fledge) & hatch>0 & year==y)
  fledgenestsSY <- parameters %>% filter (!is.na(fledge) & hatch>0 & year==y & FAge=="SY")
  fledgenestsASY <- parameters %>% filter (!is.na(fledge) & hatch>0 & year==y & FAge!="SY" & FAge != "AHY")
  
  
  PopData$fledgeRate[a] <- sum(fledgenests$fledge) / sum(fledgenests$hatch)
  PopData$fledgeRateSY[a] <- sum(fledgenestsSY$fledge) / sum(fledgenestsSY$hatch)
  PopData$fledgeRateASY[a] <- sum(fledgenestsASY$fledge) / sum(fledgenestsASY$hatch)
  
}


write.csv(PopData, file= "file:///C:/Users/Amelia/Documents/Masters Thesis Project/TRES Data Analysis/Matrix Pop Estimates/Yearly Vital Rates.csv", na="", row.names = F)

#lets just see if we have enough data about the differences between SY and ASY
#female reproduction to get estimates of clutch size, hatch rate, or fledge rate
#based on adult age.

t.test(paired=T,PopData$clutchSizeSY, PopData$clutchSizeASY) 
#there are significant differencecs so they should be treated seperately
t.test(paired=T,PopData$hatchRateSY, PopData$hatchRateASY) #also should be seperate
t.test(paired=T,PopData$fledgeRateSY, PopData$fledgeRateASY) #just boarderline but should be different

t.test(paired=T,PopData$averageNestsSY, PopData$averageNestsASY) #different