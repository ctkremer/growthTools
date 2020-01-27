## ----echo=FALSE---------------------------------------------------------------
library(knitr)
opts_chunk$set(tidy.opts=list(width.cutoff=60),tidy=TRUE)

## ----message=FALSE------------------------------------------------------------
# To direct vignette code to use the growthTools package during development, calling devtools::load_all() is suggested. When ready for release, use the library command. 
# http://stackoverflow.com/questions/35727645/devtools-build-vignette-cant-find-functions
# when switching, also remember to remove %\VignetteDepends{devtools} from header
# see https://community.rstudio.com/t/using-packages-that-are-not-imported-as-part-of-examples-in-a-vignette/2317

devtools::load_all()
#library(growthTools)

library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

## -----------------------------------------------------------------------------
# Construct example data set:
sdat<-data.frame(trt=c(rep('A',10),rep('B',10),rep('C',10),rep('D',10)),
                 dtime=rep(seq(1,10),4),
                 ln.fluor=c(c(1,1.1,0.9,1,2,3,4,5,5.2,4.7),
                            c(1.1,0.9,1,2,3,4,4.1,4.2,3.7,4)+0.3,
                            c(3.5,3.4,3.6,3.5,3.2,2.2,1.2,0.5,0.4,0.1),
                            c(5.5,4.5,3.5,2.5,1.5,0,0.2,-0.1,0,-0.1)))

## ----width=6,heigh=4.5--------------------------------------------------------
ggplot(sdat,aes(x=dtime,y=ln.fluor))+
  geom_point(aes(colour=trt))+theme_bw()+
  scale_x_continuous(name = 'Time',breaks = seq(0,10,2))+
  scale_y_continuous(name = 'ln(Abundance)')

## ----error=FALSE,width=6,heigh=4.5--------------------------------------------
# subset the data, focusing on population A:
sdat2<-sdat[sdat$trt=='A',]

# calculate growth rate using all available methods:
res<-get.growth.rate(sdat2$dtime,sdat2$ln.fluor,plot.best.Q = T,id = 'Population A')

## -----------------------------------------------------------------------------
res$ictab

## -----------------------------------------------------------------------------
res$best.model
res$best.slope

## -----------------------------------------------------------------------------
res$best.model.rsqr

## -----------------------------------------------------------------------------
res$best.se

## -----------------------------------------------------------------------------
res$best.model.slope.n
res$best.model.slope.r2

## -----------------------------------------------------------------------------
summary(res$models$gr.lag)

## -----------------------------------------------------------------------------
res$slopes
res$ses
res$slope.rs
res$slope.ns

## -----------------------------------------------------------------------------
res$best.se
res$ses

## ----error=FALSE--------------------------------------------------------------
gdat <- sdat %>% group_by(trt) %>% 
        do(grs=get.growth.rate(x=.$dtime,y=.$ln.fluor,
                               id=.$trt,plot.best.Q=T,fpath=NA))  

## -----------------------------------------------------------------------------
gdat %>% summarise(trt,mu=grs$best.slope,best.model=grs$best.model,
                   best.se=grs$best.se)

## -----------------------------------------------------------------------------
gdat %>% summarise(trt,mu=grs$best.slope,best.model=grs$best.model,
                   best.se=grs$best.se,best.R2=grs$best.model.rsqr,
                   nobs.exp=grs$best.model.slope.n)

## -----------------------------------------------------------------------------
# Only use the linear method:
gdat <- sdat %>% group_by(trt) %>% 
        do(grs=get.growth.rate(x=.$dtime,y=.$ln.fluor,id=.$trt,plot.best.Q=F,
                               methods=c('linear'))) %>% 
        summarise(trt,mu=grs$best.slope,best.model=grs$best.model)
gdat

## ----error=FALSE--------------------------------------------------------------
# Only use the linear method:
gdat <- sdat %>% group_by(trt) %>% 
        do(grs=get.growth.rate(x=.$dtime,y=.$ln.fluor,id=.$trt,plot.best.Q=F,
                               methods=c('lag','sat','flr'))) %>% 
        summarise(trt,mu=grs$best.slope,best.model=grs$best.model)
gdat

## ----error=FALSE--------------------------------------------------------------
gdat <- sdat %>% group_by(trt) %>% 
        do(grs=get.growth.rate(x=.$dtime,y=.$ln.fluor,id=.$trt,plot.best.Q=F,
                               model.selection=c('AIC'))) %>% 
        summarise(trt,mu=grs$best.slope,best.model=grs$best.model)
gdat

## -----------------------------------------------------------------------------
data("example_TPC_data")
head(example_TPC_data)

