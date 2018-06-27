## ----message=FALSE-------------------------------------------------------
# Direct vignette code to use the growthTools package. When ready for release, use the library command; during development, calling devtools::load_all() is suggested:
# http://stackoverflow.com/questions/35727645/devtools-build-vignette-cant-find-functions

#library(priceTools)
devtools::load_all()

library(ggplot2)
library(dplyr)
library(tidyr)

## ------------------------------------------------------------------------
# Construct example data set:
sdat<-data.frame(trt=c(rep('A',10),rep('B',10),rep('C',10)),dtime=rep(seq(1,10),3),ln.fluor=c(c(1,1.1,0.9,1,2,3,4,5,5.2,4.7),c(1.1,0.9,1,2,3,4,4.1,4.2,3.7,4)+0.3,c(3.5,3.4,3.6,3.5,3.2,2.2,1.2,0.5,0.4,0.1)))
sdat$id<-sdat$trt

## ------------------------------------------------------------------------
ggplot(sdat,aes(x=dtime,y=ln.fluor))+
  geom_point(aes(colour=trt))+theme_bw()+
  scale_x_continuous('Time')+
  scale_y_continuous('ln(Abundance)')

## ------------------------------------------------------------------------

# all methods
gdat <- sdat %>% group_by(trt) %>% do(grs=get.growth.rate(.$dtime,.$ln.fluor,.$id,plot.best.Q=T,fpath=NA)) %>% summarise(trt,mu=grs$best.slope,best.model=grs$best.model)
gdat



