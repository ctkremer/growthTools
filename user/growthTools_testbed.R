
library(devtools)
library(bbmle)
library(emdbook)
library(ggplot2)

devtools::load_all()

joe<-read.csv("/Users/colin/Research/Software/growthTools/example_TPC_data.csv")
head(joe)

example_TPC_data<-joe

devtools::use_data(example_TPC_data)
?use_data()

# Single TPC data set:
devtools::load_all()
sp1<-joe[joe$isolate.id=='CH30_4_RI_03' & joe$dilution==1,]
plot(mu~temperature,data=sp1)

get.nbcurve.tpc(sp1$temperature,sp1$mu,method='grid.mle2',m,plotQ=T,fpath=NA)

r1<-get.decurve.tpc(sp1$temperature,sp1$mu,method='grid.mle2',plotQ=T,fpath=NA)

temp<-sp1$temperature
mu<-sp1$mu
plot(mu~temp,ylim=c(-0.2,3))
curve(decurve(x,topt=20.9,b1 = 0.208,b2=0.123,d0=-0.08,d2=0.153),0,30,col='blue',add=T)

curve(decurve(x,topt=20.72,b1 = 0.203,b2=0.0848,d0=1.4556e-8,d2=0.18),0,30,col='blue',add=T)


#temp<-sp1$temperature
#mu<-sp1$mu

sp1<-joe[joe$isolate.id=='CH30_4_RI_03' & joe$dilution==2,]
plot(mu~temperature,data=sp1)

get.nbcurve.tpc(sp1$temperature,sp1$mu,method='grid.mle2',plotQ=T,fpath=NA)

get.decurve.tpc(sp1$temperature,sp1$mu,method='grid.mle2',plotQ=T,fpath=NA)

#temp<-sp1$temperature
#mu<-sp1$mu

# Now apply to multiple curves:


sp1b<-joe[joe$isolate.id=='CH30_4_RI_03',]

# apply get.nbcurve to data set, grouping by isolate and dilution
devtools::load_all()
res <- sp1b %>% group_by(isolate.id,dilution) %>% do(tpcs=get.nbcurve.tpc(.$temperature,.$mu,method='grid.mle2',plotQ=T,conf.bandQ=T,fpath=NA,id=.$dilution))

# or without confidence bands
res <- sp1b %>% group_by(isolate.id,dilution) %>% do(tpcs=get.nbcurve.tpc(.$temperature,.$mu,method='grid.mle2',plotQ=T,conf.bandQ=F,fpath=NA,id=.$dilution))

# or saving resulting plots:
fpath<-'/Users/colin/Research/Software/growthTools/user/'
res <- sp1b %>% group_by(isolate.id,dilution) %>% do(tpcs=get.nbcurve.tpc(.$temperature,.$mu,method='grid.mle2',plotQ=T,conf.bandQ=T,fpath=fpath,id=.$dilution))

# process results
res %>% summarise(isolate.id,dilution,topt=tpcs$o,tmin=tpcs$tmin,tmax=tpcs$tmax,rsqr=tpcs$rsqr,a=exp(tpcs$a),b=tpcs$b,w=tpcs$w)

res$tpcs

### Or with de model instead of norberg:

res <- sp1b %>% group_by(isolate.id,dilution) %>% do(tpcs=get.decurve.tpc(.$temperature,.$mu,method='grid.mle2',plotQ=T,conf.bandQ=T,fpath=NA,id=.$dilution))

res %>% summarise(isolate.id,dilution,topt=tpcs$topt,tmin=tpcs$tmin,tmax=tpcs$tmax,rsqr=tpcs$rsqr,b1=tpcs$b1,b2=tpcs$b2,d0=tpcs$d0,d2=tpcs$d2)





#### More development:

# apply get.nbcurve to data set, grouping by isolate and dilution
devtools::load_all()

res <- sp1b %>% group_by(isolate.id,dilution) %>% do(tpcs=get.nbcurve.tpc(.$temperature,.$mu,method='grid.mle2',plotQ=T,conf.bandQ=T,fpath=NA,id=.$dilution))

res2 <- res %>% summarise(isolate.id,dilution,topt=tpcs$o,tmin=tpcs$tmin,tmax=tpcs$tmax,rsqr=tpcs$rsqr,a=exp(tpcs$a),b=tpcs$b,w=tpcs$w)

# what would it take to plot all of these, with confidence bands?

# for a given row of coef's, make a row of predictions...

pred.nbcurve<-function(cf,xlim=c(0,35),Sigma,conf.bandQ=F){
  xs<-seq(xlim[1],xlim[2],length.out = 101)
  ys<-nbcurve2(xs,cf$o,cf$w,cf$a,cf$b)
  res<-data.frame(temp=xs,mu=ys)
  
  if(conf.bandQ){
    st<-paste("nbcurve2(c(",paste(xs,collapse=','),"),o,w,a,b)",sep='')
    dvs0<-suppressWarnings(deltavar2(fun=parse(text=st),meanval=cf,Sigma=Sigma))
    res$ml.ci<-1.96*sqrt(dvs0)
  }
  return(res)
}

head(res2)

#pds <- res %>% group_by(isolate.id,dilution) %>% do(preds=pred.nbcurve(cf=c(.$tpcs$o,)))

#  tpcs=get.nbcurve.tpc(.$temperature,.$mu,method='grid.mle2',plotQ=T,conf.bandQ=T,fpath=NA,id=.$dilution))

# Need to read a good tutorial on R data structures, slots, @ and $. It would be nice to store the model's coefficients as a single united list, at the same level as the vcov matrix is stored. but only if I can then access that information using the summarise command.


##############################################

# Construct example data set:
sdat<-data.frame(trt=c(rep('A',10),rep('B',10),rep('C',10),rep('D',10)),dtime=rep(seq(1,10),4),ln.fluor=c(c(1,1.1,0.9,1,2,3,4,5,5.2,4.7),c(1.1,0.9,1,2,3,4,4.1,4.2,3.7,4)+0.3,c(3.5,3.4,3.6,3.5,3.2,2.2,1.2,0.5,0.4,0.1),c(5.5,4.5,3.5,2.5,1.5,0,0.2,-0.1,0,-0.1)))


# subset the data, focusing on population D:
sdat2<-sdat[sdat$trt=='D',]

plot(ln.fluor~dtime,data=sdat2)

get.gr.flr(x=sdat2$dtime,y=sdat2$ln.fluor)

devtools::load_all()

# subset the data, focusing on population A:
sdat2<-sdat[sdat$trt=='A',]

plot(ln.fluor~dtime,data=sdat2)

# calculate growth rate using all available methods:
res<-get.growth.rate(sdat2$dtime,sdat2$ln.fluor,plot.best.Q = T,id = 'Population A')

devtools::load_all()
res<-get.growth.rate(sdat2$dtime,sdat2$ln.fluor,plot.best.Q = T,id = 'Population A',model.selection = 'AICc',method='lag')

res<-get.growth.rate(sdat2$dtime,sdat2$ln.fluor,plot.best.Q = T,id = 'Population A',model.selection = 'AICc')



gdat <- sdat %>% group_by(trt) %>% do(grs=get.growth.rate(x=.$dtime,y=.$ln.fluor,id=.$trt,plot.best.Q=T,fpath=NA)) 
gdat %>% summarise(trt,mu=grs$best.slope,best.model=grs$best.model,best.se=grs$best.se,best.n=grs$best.model.slope.n,best.r2=grs$best.model.slope.r2)


res$best.se
res$ses

x<-sdat2$dtime
y<-2*sdat2$ln.fluor

y<-0.5*sdat2$ln.fluor

# fix se's
library(minpack.lm)
lag2<-function(x,a,b,B1,s=1E-10){
  sqfunc(B1-x,b,s)-(b/2)*(B1-x)+a
}


# with exp(log(b))
fit.lag<-nlsLM(y ~ lag(x,a,logb,B1,s=1E-10),
               start = c(B1=mean(x)-(mean(x)-min(x))/2, a=min(y), logb=0),data = data.frame(x,y),
               control = nls.control(maxiter=1000, warnOnly=TRUE))

# with raw b value
fit.lag2<-nlsLM(y ~ lag2(x,a,b,B1,s=1E-10),
               start = c(B1=mean(x)-(mean(x)-min(x))/2, a=min(y), b=1),data = data.frame(x,y),
               lower = c(B1=-Inf,a=-Inf,b=0),
               control = nls.control(maxiter=1000, warnOnly=TRUE))

coef(fit.lag2)['b']
exp(coef(fit.lag)['logb'])

sqrt(diag(vcov(fit.lag)))['logb']
sqrt(diag(vcov(fit.lag2)))['b']

# How to obtain confidence intervals for nls fits?
class(res$models[[2]])
confint(res$models[[1]])
confint(res$models[[2]],param='logb')
res$models[[2]]

install.packages("nlstools")
library(nlstools)

nlstools::confint2(res$models[[2]],parm='logb',method='profile')
?confint2
