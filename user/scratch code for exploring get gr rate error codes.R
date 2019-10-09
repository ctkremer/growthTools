
library(dplyr)
library(growthTools)

fdat<-read.csv('/Users/colin/Research/Software/growthTools/data/example_time_series.csv')
#fdat<-fdat %>% mutate(isolate.id.media=paste(isolate.id,media,sep='-'))
head(fdat)

uids<-unique(fdat$id)
length(uids)

# New approach, using growthTools and improvied piecewise functions for fits
# - to run w/out making graphs, specify: plot.best.Q=F
grs2 <- fdat %>% filter(id %in% uids[5:10]) %>% group_by(experiment.ID,isolate.id,genus,species,source,plate.temperature,media,temperature,dilution,replicate) %>% do(grs=get.growth.rate(x=.$dtime,y=.$log.fluor,id=.$id,plot.best.Q=T,fpath=NA,model.selection = 'AICc'))

grs <- fdat %>% group_by(experiment.ID,isolate.id,genus,species,source,plate.temperature,media,temperature,dilution,replicate) %>% do(grs=get.growth.rate(x=.$dtime,y=.$log.fluor,id=.$id,plot.best.Q=T,fpath=NA,model.selection = 'AICc'))

# getting situations where plots fail... maybe due to new try-catch errors?

length(uids)

# also fails for uids[2]

#devtools::load_all()

# focal data set to tinker with:
tmp<-fdat %>% filter(id %in% uids[2])

x<-tmp$dtime
y<-tmp$log.fluor

m1<-get.growth.rate(x=x,y=y,id='test1',plot.best.Q=T,fpath=NA,model.selection = 'AICc',methods = c('linear'))
m1<-get.growth.rate(x=x,y=y,id='test1',plot.best.Q=T,fpath=NA,model.selection = 'AICc',methods = c('linear','lag'))  # throws singular gradient matrix at initial parameter estimates
m1<-get.growth.rate(x=x,y=y,id='test1',plot.best.Q=T,fpath=NA,model.selection = 'AICc',methods = c('linear','sat'))  # throws singular gradient matrix at initial parameter estimates
m1<-get.growth.rate(x=x,y=y,id='test1',plot.best.Q=T,fpath=NA,model.selection = 'AICc',methods = c('linear','flr'))  # throws sqrt(b * (4 * s + b * x^2)) : NaNs produced
m1<-get.growth.rate(x=x,y=y,id='test1',plot.best.Q=T,fpath=NA,model.selection = 'AICc',methods = c('linear','lagsat'))


flr(x,a=a.guess,b=min(c(-0.1,round(min(slopes),5))),B2=mean(x),s=1E-10)

8.663835e-01 1.490116e-08 5.710689e+00

flr(x,a=0.866,b=-1.490116e-08,B2=mean(x),s=1E-10)

fit.flr<-nlsLM(y ~ flr(x,a,b,B2,s=1E-10),
               start=c(a=a.guess,b=min(c(-0.1,round(min(slopes),5))),B2=mean(x)),
               data = data,
               upper = c(a=Inf,b=-0.00000001,B2=Inf),
               control = nls.control(maxiter=1000, warnOnly=TRUE))



get.gr.lag(x,y,plotQ=T,fpath=NA,id='')

curve(lag(x,a=0.5,b=0.4,B1=2,s=1E-10),min(x),max(x),col='blue')
curve(lag(x,a=0.5,b=0.4,B1=5,s=1E-10),min(x),max(x),col='red',add=T)

fit.lag<-nlsLM(y ~ lag(x,a,b,B1,s=1E-10),
               #start = c(B1=mean(x)-(mean(x)-min(x))/2, a=min(y), b=1),data = data,
               #start = c(B1=1.96, a=0.469, b=1),data = data,
               start = c(B1=1.96, a=0.469, b=0.36),data = data,
               lower = c(B1=-Inf,a=-Inf,b=0.0001),
               control = nls.control(maxiter=1000, warnOnly=TRUE))
plot(fit.lag)
abline(a = 0.95,b=-0.36)
#try(,silent=T)
#lag(x,B1=1.96, a=0.469, b=1,s=1E-10)

summary(fit.lag)

?get.growth.rate

plot(ys~xs)

methods = c("linear", "lag", "sat", "flr", "lagsat")
