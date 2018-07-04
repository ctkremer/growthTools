

library(bbmle)
library(emdbook)

devtools::load_all()

joe<-read.csv("/Users/colin/Research/Software/growthTools/example_TPC_data.csv")
head(joe)

example_TPC_data<-joe

devtools::use_data(example_TPC_data)

?use_data()

# Single TPC data set:
sp1<-joe[joe$isolate.id=='CH30_4_RI_03' & joe$dilution==1,]
plot(mu~temperature,data=sp1)

get.nbcurve.tpc(sp1$temperature,sp1$mu,method='grid.mle2',plotQ=T,fpath=NA)

#temp<-sp1$temperature
#mu<-sp1$mu

sp1<-joe[joe$isolate.id=='CH30_4_RI_03' & joe$dilution==2,]
plot(mu~temperature,data=sp1)

get.nbcurve.tpc(sp1$temperature,sp1$mu,method='grid.mle2',plotQ=T,fpath=NA)


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

