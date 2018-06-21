#Make the raw data plots of vital rates against calculated populaiton growth rate. 
library(tidyverse)
setwd("~/Masters Thesis Project/Vital Rates Paper/Figures")

vdat <- read.csv("file:///C:/Users/11arc/Documents/Masters Thesis Project/Vital Rates Paper/Yearly Vital Rates Estimates.csv") %>% arrange(year)
boxocc <- read.csv("file:///C:/Users/11arc/Documents/Masters Thesis Project/Tree Swallow Data/Amelia TRES data 1975-2016/Improved and Cleaned Data/Box Occupancy using database.csv")%>% arrange(Year)
vdat$BoxOcc <- boxocc$BoxOccTotal
vdat$calcLambda <- NA


#They're all arranged in descending order 
for(i in 1:42){
  vdat$calcLambda[i] <- vdat$BoxOcc[i+1]/vdat$BoxOcc[i]
}


PanelA <- ggplot(vdat, aes(y=averageNestsSY, x=calcLambda ))+
  geom_point()+
  labs(x="", y="Nests per SY female")+
  theme_classic(base_size = 18)

PanelD <- ggplot(vdat, aes(y=averageNestsASY, x=calcLambda ))+
  geom_point()+
  labs(x="", y="Nests per ASY female")+
  theme_classic(base_size = 18)

PanelB <- ggplot(vdat, aes(y=clutchSizeSY, x=calcLambda ))+
  geom_point()+
  labs(x="", y="SY clutch size")+
  theme_classic(base_size = 18)

PanelE <- ggplot(vdat, aes(y=clutchSizeASY, x=calcLambda ))+
  geom_point()+
  labs(x="", y="ASY clutch size")+
  theme_classic(base_size = 18)

PanelG <- ggplot(vdat, aes(y=hatchRate, x=calcLambda ))+
  geom_point()+
  labs(x="Population growth rate", y="Hatch rate")+
  theme_classic(base_size = 18)

PanelH <- ggplot(vdat, aes(y=fledgeRate, x=calcLambda ))+
  geom_point()+
  labs(x="Population growth rate", y="Fledge rate")+
  theme_classic(base_size = 18)

PanelC <- ggplot(vdat, aes(y=SYReturn, x=calcLambda ))+
  geom_point()+
  labs(x="", y="SY female return")+
  theme_classic(base_size = 18)

PanelF <- ggplot(vdat, aes(y=ASYReturn, x=calcLambda ))+
  geom_point()+
  labs(x="", y="ASY female return")+
  theme_classic(base_size = 18)

PanelI <- ggplot(vdat, aes(y=Recruitment, x=calcLambda ))+
  geom_point()+
  labs(x="Population growth rate", y="Recruitment")+
  theme_classic(base_size = 18)

cowplot::plot_grid(PanelA, PanelB, PanelC, PanelD, PanelE, PanelF, PanelG, PanelH, PanelI, 
                   align="v",
                   nrow=3, 
                   ncol=3, 
                   labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I"))
ggsave(filename="VitalRates against calculated lambda.png", width=10, height=10, units="in", device="png")






#########3Make them as historgrams because those damn graphs didn't show anything
PanelA <- ggplot(vdat, aes(averageNestsSY ))+
  geom_histogram(bins=30)+
  geom_vline(xintercept = mean(vdat$averageNestsSY, na.rm=T), linetype="dashed")+
  labs( x="Nests per SY female", y="Count")+
  scale_y_continuous(breaks=c(0,6,12,18,24))+
  
  theme_classic(base_size = 18)+
  xlim(0.95,2)

PanelD <- ggplot(vdat, aes(averageNestsASY ))+
  geom_histogram()+
  geom_vline(xintercept = mean(vdat$averageNestsASY, na.rm=T), linetype="dashed")+
  labs( x="Nests per ASY female", y="Count")+
  theme_classic(base_size = 18)+
  scale_y_continuous(breaks=c(0,3,6,9,12))+
  xlim(0.95,2)

PanelB <- ggplot(vdat, aes(clutchSizeSY))+
  geom_histogram()+
  geom_vline(xintercept = mean(vdat$clutchSizeSY, na.rm=T), linetype="dashed")+
  labs(y="", x="SY clutch size")+
  theme_classic(base_size = 18)+
  scale_y_continuous(breaks=c(0,2,4,6, 8,10))+
  xlim(4,7.2)

PanelE <- ggplot(vdat, aes(clutchSizeASY))+
  geom_histogram()+
  labs(y="", x="ASY clutch size")+
  geom_vline(xintercept = mean(vdat$clutchSizeASY, na.rm=T), linetype="dashed")+
  theme_classic(base_size = 18)+
  scale_y_continuous(breaks=c(0,2,4,6,8))+
  
  xlim(4,7.2)

PanelG <- ggplot(vdat, aes(hatchRate ))+
  geom_histogram(bins=25)+
  geom_vline(xintercept = mean(vdat$hatchRate, na.rm=T), linetype="dashed")+
  labs(y="Count", x="Hatch rate")+
  theme_classic(base_size = 18)

PanelH <- ggplot(vdat, aes(fledgeRate))+
  geom_histogram(bins=25)+
  geom_vline(xintercept = mean(vdat$fledgeRate, na.rm=T), linetype="dashed")+
  labs(y="", x="Fledge rate")+
  theme_classic(base_size = 18)

PanelC <- ggplot(vdat, aes(SYReturn))+
  geom_histogram(bins=25)+
  geom_vline(xintercept = mean(vdat$SYReturn, na.rm=T), linetype="dashed")+
  labs(x="SY female return", y="")+
  theme_classic(base_size = 18)+
  scale_x_continuous(breaks=c(0,0.3,0.3,0.9))+
  xlim(0,1)

PanelF <- ggplot(vdat, aes(ASYReturn))+
  geom_histogram(bins=25)+
  geom_vline(xintercept = mean(vdat$ASYReturn, na.rm=T), linetype="dashed")+
  labs(x="ASY female return", y="")+
  theme_classic(base_size = 18)+ 
  xlim(0,1)

PanelI <- ggplot(vdat, aes(Recruitment ))+
  geom_histogram(bins=25)+
  geom_vline(xintercept = mean(vdat$Recruitment, na.rm=T), linetype="dashed")+
  labs(x="Recruitment", y="")+
  theme_classic(base_size = 18)

cowplot::plot_grid(PanelA, PanelB, PanelC, PanelD, PanelE, PanelF, PanelG, PanelH, PanelI, 
                   align="v",
                   nrow=3, 
                   ncol=3, 
                   labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I"))
ggsave(filename="Histograms of vital rates.png", width=11, height=10, units="in", device="png")
ggsave(filename="Fig 2 Histograms.eps", width=11, height=10, units="in", device="eps")
