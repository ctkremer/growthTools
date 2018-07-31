## ----message=FALSE-------------------------------------------------------
# To direct vignette code to use the growthTools package during development, calling devtools::load_all() is suggested. When ready for release, use the library command. 
# http://stackoverflow.com/questions/35727645/devtools-build-vignette-cant-find-functions

devtools::load_all()
#library(growthTools)

library(ggplot2)
library(dplyr)
library(tidyr)
library(mleTools)

## ------------------------------------------------------------------------
# Construct example data set:
sdat<-data.frame(trt=c(rep('A',10),rep('B',10),rep('C',10)),dtime=rep(seq(1,10),3),ln.fluor=c(c(1,1.1,0.9,1,2,3,4,5,5.2,4.7),c(1.1,0.9,1,2,3,4,4.1,4.2,3.7,4)+0.3,c(3.5,3.4,3.6,3.5,3.2,2.2,1.2,0.5,0.4,0.1)))

## ------------------------------------------------------------------------
ggplot(sdat,aes(x=dtime,y=ln.fluor))+
  geom_point(aes(colour=trt))+theme_bw()+
  scale_x_continuous('Time')+
  scale_y_continuous('ln(Abundance)')

## ------------------------------------------------------------------------
# subset the data, focusing on population A:
sdat2<-sdat[sdat$trt=='A',]

# calculate growth rate using all available methods:
res<-get.growth.rate(sdat2$dtime,sdat2$ln.fluor,plot.best.Q = T,id = 'Population A')

## ------------------------------------------------------------------------
res$best.model
res$best.slope

## ------------------------------------------------------------------------
summary(res$models$gr.lag)

## ------------------------------------------------------------------------
res$slopes

## ------------------------------------------------------------------------
res$bes.se
res$ses

## ------------------------------------------------------------------------
gdat <- sdat %>% group_by(trt) %>% do(grs=get.growth.rate(x=.$dtime,y=.$ln.fluor,id=.$trt,plot.best.Q=T,fpath=NA))  

## ------------------------------------------------------------------------
gdat %>% summarise(trt,mu=grs$best.slope,best.model=grs$best.model,best.se=grs$best.se)

## ------------------------------------------------------------------------
# Only use the linear method:
gdat <- sdat %>% group_by(trt) %>% do(grs=get.growth.rate(x=.$dtime,y=.$ln.fluor,id=.$trt,plot.best.Q=T,methods=c('linear'))) %>% summarise(trt,mu=grs$best.slope,best.model=grs$best.model)
gdat

## ------------------------------------------------------------------------
head(example_TPC_data)

## ------------------------------------------------------------------------
# Single data set:
sp1 <- example_TPC_data %>% filter(isolate.id=='CH30_4_RI_03' & dilution==1)

# obtain Norberg curve parameters, using a grid search algorithm to consider a range of possible initial parameter values, and plot the results
nbcurve.traits<-get.nbcurve.tpc(sp1$temperature,sp1$mu,method='grid.mle2',plotQ=T,conf.bandQ = T,fpath=NA)
data.frame(nbcurve.traits)

## ----echo=F--------------------------------------------------------------
# First, let's look at an example with multiple dilutions but the same strain:
sp1b <- example_TPC_data %>% filter(isolate.id=='CH30_4_RI_03')

# apply get.nbcurve to the entire data set, grouping by isolate and dilution
res <- sp1b %>% group_by(isolate.id,dilution) %>% do(tpcs=get.nbcurve.tpc(.$temperature,.$mu,method='grid.mle2',plotQ=T,conf.bandQ=T,fpath=NA,id=.$dilution))


## ------------------------------------------------------------------------
# or saving resulting plots:
#fpath<-'/Users/colin/Research/Software/growthTools/user/'

# provide an explicit fpath to invoke plot saving; when `id` is also provided, this column will be used to produce the plot's title, as well as included in the file name.
#res <- sp1b %>% group_by(isolate.id,dilution) %>% do(tpcs=get.nbcurve.tpc(.$temperature,.$mu,method='grid.mle2',plotQ=T,conf.bandQ=T,fpath=fpath,id=.$dilution))

## ------------------------------------------------------------------------
nb.res <- sp1b %>% group_by(isolate.id,dilution) %>% do(tpcs=get.nbcurve.tpc(.$temperature,.$mu,method='grid.mle2',plotQ=T,conf.bandQ=F,fpath=NA,id=.$dilution))

sp1b %>% group_by(isolate.id,dilution) %>% do(tpcs=get.nbcurve.tpc(.$temperature,.$mu,method='grid.mle2',plotQ=F,conf.bandQ=F,fpath=NA,id=.$dilution))

## ------------------------------------------------------------------------
# process results
nb.res2<-nb.res %>% summarise(isolate.id,dilution,topt=tpcs$o,tmin=tpcs$tmin,tmax=tpcs$tmax,rsqr=tpcs$rsqr,a=exp(tpcs$a),b=tpcs$b,w=tpcs$w)

nb.res2

## ----warning=F-----------------------------------------------------------
de.res <- sp1b %>% group_by(isolate.id,dilution) %>% do(tpcs=get.decurve.tpc(.$temperature,.$mu,method='grid.mle2',plotQ=T,conf.bandQ=T,fpath=NA,id=.$dilution))

de.res2 <- de.res %>% summarise(isolate.id,dilution,topt=tpcs$topt,tmin=tpcs$tmin,tmax=tpcs$tmax,rsqr=tpcs$rsqr,b1=tpcs$b1,b2=tpcs$b2,d0=tpcs$d0,d2=tpcs$d2)

de.res2

## ------------------------------------------------------------------------
library(reshape2)

nb.res3<-melt(nb.res2,id.vars=c('isolate.id','dilution'))
de.res3<-melt(de.res2,id.vars=c('isolate.id','dilution'))
res2<-rbind(data.frame(type='nb',nb.res3),data.frame(type='de',de.res3))

ggplot(res2[res2$variable=='rsqr',],aes(x=type,y=value))+
  geom_boxplot(aes(fill=type))+
  scale_y_continuous(limits=c(0,1))

## ----warning=F-----------------------------------------------------------
table(example_TPC_data[,c('isolate.id','dilution')])

# create informative ID column for each combo of unique strain and dilution period:
example_TPC_data$id<-paste(example_TPC_data$isolate.id,example_TPC_data$dilution)

res2 <- example_TPC_data %>% group_by(isolate.id,dilution) %>% do(tpcs=get.nbcurve.tpc(.$temperature,.$mu,method='grid.mle2',plotQ=F,conf.bandQ=F,fpath=NA,id=.$id))

clean.res <- res2 %>% summarise(isolate.id,dilution,topt=tpcs$o,tmin=tpcs$tmin,tmax=tpcs$tmax,rsqr=tpcs$rsqr,a=exp(tpcs$a),b=tpcs$b,w=tpcs$w)
data.frame(clean.res)

