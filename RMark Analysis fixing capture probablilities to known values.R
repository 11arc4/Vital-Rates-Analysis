
setwd("~/Masters Thesis Project/TRES Data Analysis/RMark Preliminary Survival Analysis")

###################################### Let's do the nestling analysis, using the
#best capture estimates from the female analysis, and setting capture of
#nestlings at 100% because that's a fairly reasonable assumption.

datMark_nestling<- readRDS("Nestling MARK Data.rda")

nestlingprocess <- process.data(datMark_nestling,model="CJS",begin.time=1975, groups= c("age"), initial.ages =c(0))
#initial age for all the hatchlings is 0 because here I am only looking at those birds we first saw as nestlings
nestling.ddl<- make.design.data(nestlingprocess, parameters=list(Phi=list(age.bins=c(0,1, 2, 43)),
                                                                 p=list(age.bins=c(0,1,2,43)))) 
#again we are binning the ages into HY, SY and ASY (AHY will be assigned 1 just like SY) ie age 0, 1, and 2-13

nestling.ddl$p$fix<- NA

#Fix capture history of adults as the % of birds known for each year. We'll
#assume SY and ASY birds are the same, since we don't know for males.
BirdsSeen <- read.csv("file:///C:/Users/11arc/Documents/Masters Thesis Project/TRES Data Analysis/Extracted Statistics/CaptureHistory.csv",
                      na.strings=c(""), as.is=T)

for (year in 1975:2017){
  nestling.ddl$p$fix[nestling.ddl$p$t==year]<- BirdsSeen$CaptureRate[which(BirdsSeen$year==year)]
}
#This set all the neslting capture histories as well. I'd like to overwrite that and 
#Fix capture history of nestlings at 100%
nestling.ddl$p$fix[nestling.ddl$p$age=="[0,1]"]<- 1



#Again, time could be continuous or discrete but for the vital rates analysis we
#need discete so all continuous times are commented out
phi.dot <- list(formula= ~1) #if survival is constant in all cases
#phi.Time <- list(formula=~Time) #if survival rates change over the years
phi.age <- list(formula=~age) #If older birds (based on having seen them previously, not their known age) have lower survival
#phi.age.Time <- list(formula=~Time*age) #if both age (based on sight history) have an effecct
#phi.age.plus.Time <- list(formula=~Time+age)
#For the vital rates analysis we will want time to not be continuous
phi.time <- list(formula=~time)
phi.age.time <- list(formula=~time*age)
phi.age.plus.time <- list(formula=~time+age)



n1_fixed<- mark(nestlingprocess, nestling.ddl, model.parameters = list(Phi=phi.dot), output = F, adjust=T)
#YESSSSS this is setting all the fixed values for capture. I'm totally impressed. 
n2_fixed<- mark(nestlingprocess, nestling.ddl, model.parameters = list(Phi=phi.age), output = F, adjust=T)
n3_fixed<- mark(nestlingprocess, nestling.ddl, model.parameters = list(Phi=phi.time), output = F, adjust=T)
n4_fixed<- mark(nestlingprocess, nestling.ddl, model.parameters = list(Phi=phi.age.time), output = F, adjust=T)
n5_fixed <- mark(nestlingprocess, nestling.ddl, model.parameters = list(Phi=phi.age.plus.time), output = F, adjust=T)


nestling.cjs.results <- collect.models()
#based on this there is a clear winner and that's the time + age model! Hallelujah. 
bestnestlingmod <- n5
bestnestlingmod$reals$Phi

saveRDS(n5, "Best Nestling MARK Results_fixed.rda")

#how do these results compare to the results I got when I let capture rate vary? 


#These results are quite different from the previous results I got. I'd like to
#compare the capture rates I assumed, to those that the best female model chose
#when I let it choose it's own.

 unfixedFMod <- readRDS("Best Female MARK Results without fixing capture.rda")
 summaryCaptureF <- summary(unfixedFMod)
FCapture <- c(NA, summaryCaptureF$reals$p$`Group:age1`$pim[1,]) #can't estimate 1975 capture rate so we drop that

plot(FCapture~BirdsSeen$FCaptureRate)
abline(a=0, b=1) 
#Appears that we underestimate our capture rates (because we have more renests
#than we think?) Maybe it would work better if I used the capture rates we
#estimated for the adults to parameterize the nestlings


##################################
#Now lets do the females, similarly fixing capture rates!

datMark_F <- readRDS("Adult Female MARK Data.rda") 
#this data is only for females. All males have been removed as have observations
#of nestlings (if the bird was seen as a nestling the nestling observation is 
#removed and it shows up first at age 1). I did this because I don't know which 
#nestlings are female and which are not if they didn't recruit so I have to 
#consider nestling recruitment independent of sex. We could maye make the
#assumption that 50% of fledgelings are female because 50% of eggs are female,
#but there could easily be unequal survival since that occurs in other species.I
#thought fewer assumptions the better, so I will do a seperate analysis for the nestlings
tsprocess <-process.data(datMark_F,model="CJS",
                         begin.time=1975, 
                         groups= ("age"),
                         initial.ages =c(1, 2))
ts.ddl<- make.design.data(tsprocess, parameters=list(Phi=list(age.bins=c(1, 2, 43)),
                                                     p=list(age.bins=c(1,2, 43)))) 

ts.ddl$p$fix<- NA

for (year in 1975:2017){
 ts.ddl$p$fix[ts.ddl$p$t==year]<- BirdsSeen$FCaptureRate[which(BirdsSeen$year==year)]
}


ts1_fixed<- mark(tsprocess, ts.ddl, model.parameters = list(Phi=phi.dot), output = F, adjust=T)
ts2_fixed<- mark(tsprocess, ts.ddl, model.parameters = list(Phi=phi.age), output = F, adjust=T)
ts3_fixed<- mark(tsprocess, ts.ddl, model.parameters = list(Phi=phi.age.time), output = F, adjust=T)
ts4_fixed<- mark(tsprocess, ts.ddl, model.parameters = list(Phi=phi.age.plus.time), output = F, adjust=T)
ts5_fixed<- mark(tsprocess, ts.ddl, model.parameters = list(Phi=phi.time), output = F, adjust=T)

tsF.cjs.results <- collect.models()
#Nice job! This is fantastic. Both the nesltings and the adult survival models
#pop out the same (time+age being important). This makes a whole lot of sense
#actually, and is the most biologically meaningful. 

bestFmod_fixed <- ts4
summaryF_fixed <- summary(bestFmod_fixed)

saveRDS(ts4, "Best Female MARK Results capture fixed.rda")

bestFmod_unfixed <- readRDS("Best Female MARK Results without fixing capture.rda")
summaryF_unfixed <- summary(bestFmod_unfixed)
#How do these results compare to the results indicated earlier when I let the model determine capture rates? 

plot(bestFmod_fixed$results$real$estimate~bestFmod_unfixed$results$real$estimate)
abline(a=0, b=1) 
#OK so this isn't actually too bad. THey're not super tight to the line, but we tend to do get similar answers either way. 




