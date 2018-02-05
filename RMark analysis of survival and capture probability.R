#RMark analysis of birds
setwd("~/Masters Thesis Project/TRES Data Analysis/Vital Rates Models/RMark Preliminary Survival Analysis")
#Make sure you set the working directory to something seperarate otherwise MARK
#will make 1000 files that you'll want to get rid of later.

datMark_F <- readRDS("Adult Female MARK Data.rda") 
#this data is only for females. All males have been removed as have observations
#of nestlings (if the bird was seen as a nestling the nestling observation is 
#removed and it shows up first at age 1). I did this because I don't know which 
#nestlings are female and which are not if they didn't recruit so I have to 
#consider nestling recruitment independent of sex. We could maye make the
#assumption that 50% of fledgelings are female because 50% of eggs are female,
#but there could easily be unequal survival since that occurs in other species.I
#thought fewer assumptions the better, so I will do a seperate analysis for the nestlings

library(RMark)
#I've created multiple groups based on sex and at what age we first saw the
#bird. Also set the initial year to 1975.

# We're going to want a CJS #model because our data is recapture data, and we
# can't tell whether you've #actually died. (Although ultimately I may want to
# incorporate information about #the birds that were killed for experiments and
# therefore will want to use a #different type of model for that) 
tsprocess <-process.data(datMark_F,model="CJS",
                         begin.time=1975, 
                         groups= ("age"),
                         initial.ages =c(1, 2))
  
  #age vector is the initial ages, "initial.ages" assigns
# a dummy variable to #count up from for each of those. So I could use AHY ASY
# and SY and assign c(1, #2, 1) instead (it's done alphabetically) That
# condensed the data into that form
#where you have the different numbers of times you've seen a particular capture
#probability.
ts.ddl<- make.design.data(tsprocess, parameters=list(Phi=list(age.bins=c(1, 1.5, 43)),
                                                     p=list(age.bins=c(1,1.5, 43)))) 

#Settin age bins to look at birds from 1-2 (SY) and then over 2 (ASY)

#Adding in yearly box Occupancy as a covariant
boxOcc <- read.csv("~/Masters Thesis Project/Tree Swallow Data/Amelia TRES data 1975-2016/Improved and Cleaned Data/Box Occupancy using nestdata renest function.csv")
for (i in 1:nrow(ts.ddl$Phi)){
  ts.ddl$p$BoxOccupancy[i] <- boxOcc$BoxOccTotal[which(boxOcc$Year==ts.ddl$Phi$time[i])]
}



#Now there are a number of different models that I would really like to run

#I will want to include Time (continuous), Age (discrete but based on the number
#of times we've captured the bird),Sex as those factors are
#what might be influencing both capture success and survival. Also make AgeAtFirstSight ()
#I'm going to simplify this even further so that you are either an adult or a neslting at first sighting. 

#I KNOW that we caught different numbers of birds each year so time should go in
#all the capture probabilities. We might also be more likely to catch a ASY bird
#than a SY birds (that's what the nestling records show) but since SY is the
#first time we'd catch these birds we can't look at that. Capture might also
#depend on box occupancy! If box occupancy is high, birds might not have as high
#a capture record. The effect of box occupancy should be stable through time
#though so there's no need to test for an interaction

p.time <- list(formula= ~time)
p.dot <- list(formula= ~1)
#p.boxOcc <- list(formula = ~ BoxOccupancy)
#p.timeboxOcc <- list(formula= ~ time + BoxOccupancy)


#Survival probabilities might also vary by time but now they don't HAVE to, and 
#that time function is continuous (Time not time). These continuous time
#survival probabilities havve been commented out for now, but will more than
#likely be important later on

#I will use discrete time (time) in my RMark alysis for our vital rates anlaysis
#though since that will be a better estimator of variation in survival over the
#years

Phi.dot <- list(formula= ~1) #if survival is constant in all cases
Phi.Age <- list(formula=~age) #If older birds (based on having seen them previously, not their known age) have lower survival
Phi.time <- list(formula=~time) #if survival rates change over the years
Phi.Age.time <- list(formula=~time*age) #if both age (based on sight history) have an effecct
Phi.Age.plus.time <- list(formula=~time+age)




#Let's build all our possible options for models, and then compare them all!
cml <- create.model.list("CJS")
ts.cjs.results <- mark.wrapper(cml, data=tsprocess, ddl=ts.ddl, output=F, adjust=F)

#Check goodness of fit on our most saturated model to make sure the model is OK

RELEASEresults <- release.gof(tsprocess)
#Test 2 isn't violated-- this means that all marked individuals are equally
#likely to be captured. 

#Test 3 is violated quite badly-- this means that all
#marked individuals are NOT eqully likely to survive. That's quite an issue
#acctually

chat <- RELEASEresults$Chi.square[3]/RELEASEresults$df[3]
#Chat isn't 1 so we should adjust it


ts.cjs.results.adjusted <- adjust.chat(chat, ts.cjs.results)
#Minimum adequate model is still phi= time + age and p=time but p=time+ box occupancy is equivalent. 
#I wonder if we shouldn't not adjust chat because out chat is less than 3 and therefore not a problem?



saveRDS(ts.cjs.results.adjusted[[4]], "Best Female MARK Results without fixing capture.rda")


ts.cjs.results.adjusted[[6]]$pims$Phi[[2]]$pim

AvResults_female <- model.average(ts.cjs.results)

ASYPIMs <- ts.cjs.results.adjusted[[6]]$pims$Phi[[1]]$pim[1,]

#along the diagonal for SY
SYPIMs <- rep(NA, 42)
for(i in 1:42){
  SYPIMs[i] <- ts.cjs.results.adjusted[[6]]$pims$Phi[[1]]$pim[i,i]
}
FemaleResults <- as.data.frame(matrix(nrow=42, ncol=5))
colnames(FemaleResults)= c("Year",  "SYReturn", "seSYReturn", "ASYReturn", "seASYReturn")
FemaleResults$Year<- seq(1975,2016, 1)

FemaleResults[,2:3]<- AvResults_female[SYPIMs,2:3]
FemaleResults[,4:5]<- AvResults_female[ASYPIMs,2:3]


ggplot(FemaleResults, aes(x=Year, y=SYReturn))+
  geom_point(aes(y=SYReturn), color="brown")+
  geom_segment(aes(x=Year, xend=Year, y=SYReturn+seSYReturn, yend=SYReturn-seSYReturn), color="brown")+
  geom_point(aes(y=ASYReturn), color="blue")+
  geom_segment(aes(x=Year, xend=Year, y=ASYReturn+seASYReturn, yend=ASYReturn-seASYReturn), color="blue")+
  xlab("Year")+
  ylab("Apparent Survival")+
  ylim(0,1)+
  geom_smooth( formula=lm(y~x, weights=seSYReturn))

write.csv(FemaleResults, file="Model Averaged Female Yearly Apparent Survival.csv", row.names = F, na="")
  
png(filename = "~/Masters Thesis Project/Committee Meetings/Dec 15 2017/Apparent Survival with Discrete year.png",
    width =600, height = 300)

ggplot(FemaleResults %>% filter(Year>1989), aes(x=Year, y=SYReturn))+
  geom_point(aes(y=SYReturn), color="chocolate4")+
  geom_segment(aes(x=Year, xend=Year, y=SYReturn+seSYReturn, yend=SYReturn-seSYReturn), color="chocolate4")+
  geom_point(aes(y=ASYReturn), color="deepskyblue4")+
  geom_segment(aes(x=Year, xend=Year, y=ASYReturn+seASYReturn, yend=ASYReturn-seASYReturn), color="deepskyblue4")+
  xlab("Year")+
  ylab("Apparent \n Survival")+
  ylim(0,1)+
  ggthemes::theme_few(base_size = 20)+
  theme(axis.title.y=element_text(angle=0, vjust=0.5))
dev.off()

################## Now lets do an analysis of the nestlings. I will again bin
#the ages. The nestling data has been created by taking all the fledgelings, and
#following them thorugh the years. Birds that were seen first as adults are not
#included because we can (almost) guarentee that they weren't nestlings at QUBS
#nestboxes since all nestlings are banded
datMark_nestling<- readRDS("Nestling MARK Data.rda")

nestlingprocess <- process.data(datMark_nestling,model="CJS",begin.time=1975, groups= c("age"), initial.ages =c(0))
#initial age for all the hatchlings is 0 because here I am only looking at those birds we first saw as nestlings
nestling.ddl<- make.design.data(nestlingprocess, parameters=list(Phi=list(age.bins=c(0,0.8, 1.8, 43)),
                                                     p=list(age.bins=c(0,0.8,1.8,43)))) 
#again we are binning the ages into HY, SY and ASY (AHY will be assigned 1 just like SY) ie age 0, 1, and 2-13
boxOcc <- read.csv("~/Masters Thesis Project/Tree Swallow Data/Amelia TRES data 1975-2016/Improved and Cleaned Data/Box Occupancy using nestdata renest function.csv")

for (i in 1:nrow(nestling.ddl$p)){
  nestling.ddl$p$BoxOccupancy[i] <- boxOcc$BoxOccTotal[which(boxOcc$Year==nestling.ddl$Phi$time[i])]
}

#Now we need to make the capture functions (p). Capture will also depend on sex 
#but I'm a bit unsure how to take that into account since I don't know the sex 
#for so many of the birds..... it's also super biased because I only know the
#sex for birds that survived..... If birds are likely to be kicked off the grid
#when they are SY because they're subordinant, then we might expect age to
#affect capture function as well
p.time <- list(formula= ~time)
p.dot <- list(formula= ~1)
p.timeplusage <- list(formula= ~time+age) #what if we don't do a good job catching recruits (ie maybe they are subordinant and aren't allowed boxes)?
p.age <- list(formula= ~age)
#p.boxOcc <- list(formula = ~ BoxOccupancy)
#p.timeboxOcc <- list(formula= ~ time + BoxOccupancy)
#p.timeageboxOcc <- list(formula= ~ time + age + BoxOccupancy)
#p.timeagebyboxOcc <- list(formula= ~ time + age * BoxOccupancy)



#Again, time could be continuous or discrete but for the vital rates analysis we
#need discete so all continuous times are commented out
Phi.dot <- list(formula= ~1) #if survival is constant in all cases
Phi.age <- list(formula=~age) #If older birds (based on having seen them previously, not their known age) have lower survival
Phi.time <- list(formula=~time)
Phi.age.time <- list(formula=~time*age)
Phi.age.plus.time <- list(formula=~time+age)

cml2 <- create.model.list("CJS")
nestling.cjs.results <- mark.wrapper(cml2, data=nestlingprocess, ddl=nestling.ddl, output=F, adjust=F)


#Let's just check the model assumptions

RELEASEnestling <- release.gof(nestlingprocess)
#WOW the nestlings are just terrible. Much much worse than the adults
#I wonder if perhaps this is due to permenent emmigration out the of system.
#Actually Gary White says this is much more likely to be due to the sparsity of the data and the number of capture occasions that we have-- this method gets driven up substantially due to that

nestlingChat <- RELEASEnestling$Chi.square[3]/RELEASEnestling$df[3] #That's not 1 at all!!! Severely over dispersed. 


#Let's try to correct the models to take that into account

nestling.cjs.results.adjusted <- adjust.chat(nestlingChat, nestling.cjs.results)


#Now there is only one top model BUT it doesn't show any variation between years,
#just within years. Therefore we can't figure out what is going on in terms of
#correlation structure #also no changes so we can just carry on as though we didn't do this bullshit in the paper
saveRDS(nestling.cjs.results.adjusted[[1]], "Best Nestling MARK Results_capture rate unfixed.rda")



AvResults_nestling<- model.average(nestling.cjs.results)
#use the global model to help us pull out numbers
nestling.cjs.results.adjusted[[12]]$pims$Phi

NestlingResults <- as.data.frame(matrix(nrow=42, ncol=7))
colnames(NestlingResults)= c("Year", "Recruitment", "seRecruitment", "SYReturn", "seSYReturn", "ASYReturn", "seASYReturn")
NestlingResults$Year<- seq(1975,2016, 1)

#PIMS for Recruitment
#have 42 year where there's recruitment
RecruitmentPIMS <- rep(NA, 42)
for(i in 1:42){
  #PIMS for Recruitment
  RecruitmentPIMS[i]<- nestling.cjs.results[[12]]$pims$Phi[[1]]$pim[i,i]
}
NestlingResults[,2:3]<- AvResults_nestling[RecruitmentPIMS,2:3]
#PIMS for SY  Return--only 41 years where there are SYs returning bause everone starts as a nestling in 1975
SYReturnPIMs <- rep(NA, 41)
for(i in 1:41){
  #PIMS for Recruitment
  SYReturnPIMs[i]<- nestling.cjs.results[[12]]$pims$Phi[[1]]$pim[i,i+1]
}
NestlingResults[2:42,4:5]<- AvResults_nestling[SYReturnPIMs,2:3]

#PIMS for ASY Return
#PIMS for SY  Return--only 40 years where there are ASYs returning bause everone starts as a nestling in 1975
ASYReturnPIMs <- seq(3,42,1)
NestlingResults[3:42,6:7]<- AvResults_nestling[ASYReturnPIMs,2:3]

#Here are my yearly results from the nestlings!!

write.csv(NestlingResults, file="Model Averaged Nestling Yearly Apparent Survival.csv", row.names = F, na="")
