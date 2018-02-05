#Because records are only used FORWARD not backward, let's try to use all the
#adults to figure out what survival for nestlings is.
library(RMark)
library(ggplot2)
setwd("~/Masters Thesis Project/TRES Data Analysis/Vital Rates Models/RMark Preliminary Survival Analysis")
datMark <- readRDS("All birds MARK Data.rda")



nestlingprocess <- process.data(datMark,model="CJS",begin.time=1975, groups= c("age"), initial.ages =c(0,1,2))
#initial age for all the hatchlings is 0 because here I am only looking at those birds we first saw as nestlings

RELEASEresults <- release.gof(nestlingprocess)

#Hmm that looks bad still. I wonder if including some sort of transient movement would help? 
chat <- RELEASEresults$Chi.square[3]/RELEASEresults$df[3]
#This is below 3, which was the goal. It's still not perfect, but probably good enough to be going with

nestling.ddl<- make.design.data(nestlingprocess, parameters=list(Phi=list(age.bins=c(0,0.8, 1.8, 44)),
                                                                 p=list(age.bins=c(0,0.8,1.8,44)))) 
#again we are binning the ages into HY, SY and ASY (AHY will be assigned 1 just like SY) ie age 0, 1, and 2-13


p.time <- list(formula= ~time)
p.dot <- list(formula= ~1)
p.timeplusage <- list(formula= ~time+age) #what if we don't do a good job catching recruits (ie maybe they are subordinant and aren't allowed boxes)?
p.age <- list(formula= ~age)


Phi.dot <- list(formula= ~1) #if survival is constant in all cases
Phi.Age <- list(formula=~age) #If older birds (based on having seen them previously, not their known age) have lower survival
Phi.time <- list(formula=~time) #if survival rates change over the years
Phi.Age.time <- list(formula=~time*age) #if both age (based on sight history) have an effecct
Phi.Age.plus.time <- list(formula=~time+age)


cml <- create.model.list("CJS")
ts.cjs.results <- mark.wrapper(cml, data=nestlingprocess, ddl=nestling.ddl, output=F, adjust=F)



ts.cjs.results[[6]]$pims$Phi[[1]]$pim

AvResults <- model.average(ts.cjs.results)

AllResults <- as.data.frame(matrix(nrow=42, ncol=7))
colnames(AllResults)= c("Year", "Recruitment", "seRecruitment", "SYReturn", "seSYReturn", "ASYReturn", "seASYReturn")
AllResults$Year<- seq(1975,2016, 1)



RecruitmentPIMS <- rep(NA, 42)
for(i in 1:42){
  #PIMS for Recruitment
  RecruitmentPIMS[i]<- ts.cjs.results[[12]]$pims$Phi[[1]]$pim[i,i]
}
AllResults[,2:3]<- AvResults[RecruitmentPIMS,2:3]
#PIMS for SY  Return--only 41 years where there are SYs returning bause everone starts as a nestling in 1975
SYReturnPIMs <- rep(NA, 41)
for(i in 1:41){
  #PIMS for Recruitment
  SYReturnPIMs[i]<- ts.cjs.results[[6]]$pims$Phi[[1]]$pim[i,i+1]
}
AllResults[2:42,4:5]<- AvResults[SYReturnPIMs,2:3]

#PIMS for ASY Return
#PIMS for SY  Return--only 40 years where there are ASYs returning bause everone starts as a nestling in 1975
ASYReturnPIMs <- seq(3,42,1)
AllResults[3:42,6:7]<- AvResults[ASYReturnPIMs,2:3]


ggplot(AllResults, aes(x=Year))+
  #geom_point(aes(y=SYReturn))+
  #geom_segment(aes(xend=Year, y=SYReturn-seSYReturn, yend=SYReturn+seSYReturn))+
  #geom_point(aes(y=ASYReturn), color="blue")+
  #geom_segment(aes(xend=Year, y=ASYReturn-seASYReturn, yend=ASYReturn+seASYReturn), color="blue")+
  geom_point(aes(y=Recruitment), color="black")+
  geom_segment(aes(xend=Year, y=Recruitment-seRecruitment, yend=Recruitment+seRecruitment), color="black")+
  geom_smooth(aes(y=Recruitment), method="lm")+
  ggthemes::theme_few()
  
  
#This looks good and really quite reasonable! I am going to go with this data
#for the vital rates analysis.



write.csv(AllResults, file="Model Averaged Yearly Apparent Survival for all records to get recruitment.csv", row.names = F, na="")

