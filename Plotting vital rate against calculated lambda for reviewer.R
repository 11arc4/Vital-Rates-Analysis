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
