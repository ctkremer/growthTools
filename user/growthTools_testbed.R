

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


head(sp1b)


########## OLD:


####
i<-1
uis<-unique(gr.ests2$isolate.id)
cfs<-NULL
pred.mat<-c()
cfs<-data.frame(matrix(NA,nrow=length(uis),ncol=12))
names(cfs)<-c('isolate.id','w','topt','a','b','s','rsqr','tmin','tmax','topt.ci','tmin.ci','tmax.ci')
for(i in 1:length(uis)){
  print(i/length(uis))
  joe<-ungroup(gr.ests2[gr.ests2$isolate.id==uis[i],])
  
  # grid search alternative:
  grids<-list(w=seq(15,40,5),o=seq(15,25,5))
  start<-list(w=NA,o=NA,a=-1.11,b=0.05,s=10)
  
  fit2<-grid.mle2(minuslogl=mu~dnorm(mean=nbcurve(temperature,o,w,a,b),sd=s),grids=grids,start=start,data=joe)
  summary(fit2$res.best)
  fit.grid<-fit2$res.best
  cfg<-coef(fit.grid)
  # problem: we lose 'predict' functionality with the grid.mle2 approach. 
  # could re-write grid.mle2, or do a second quick fit?
  
  guesses<-as.list(cfg)
  fit<-mle2(mu~dnorm(mean=nbcurve(temperature,o,w,a,b),sd=s),
            start=guesses,data=joe)
  #summary(fit)
  #confint.FI(fit)
  
  cf<-as.list(coef(fit))
  
  # additional responses/traits:
  tmin<-get.tlim(cf$o,cf$w,cf$a,cf$b,type='tmin')
  tmax<-get.tlim(cf$o,cf$w,cf$a,cf$b,type='tmax')
  
  # calculate R2
  expc<-predict(fit)
  obs<-joe$mu
  rsqr<-1-sum((obs-expc)^2)/sum((obs-mean(obs))^2)
  
  # visualize results  
  xs<-seq(min(joe$temperature),max(joe$temperature),0.1)
  new.data<-data.frame(temperature=xs)
  new.data$mu<-predict(fit,newdata=new.data)
  new.data$isolate.id<-uis[i]
  new.data$sps<-joe$sps[1]
  
  # Follow the delta method (see Bolker book, pg. 255):
  # not 100% convinced that this is correct, given parameter rescaling tricks employed above... think about this.
  dvs<-suppressWarnings(deltavar(fun=nbcurve(xs,o,w,a,b),meanval=cf,Sigma=vcov(fit)))
  new.data$ml.ci<-1.96*sqrt(dvs)
  
  focal.dvs<-suppressWarnings(deltavar(fun=nbcurve(c(cf$o,tmin,tmax),o,w,a,b),meanval=cf,Sigma=vcov(fit)))
  focal.ci<-1.96*sqrt(focal.dvs)
  
  # plot it out
  ggplot(joe,aes(x=temperature,y=mu))+
    geom_point()+
    geom_ribbon(data=new.data,aes(ymin=mu-ml.ci,ymax=mu+ml.ci),alpha=0.2)+
    geom_line(data=new.data)+
    geom_hline(yintercept = 0)+
    theme_bw()
  
  # save output:
  vec<-c(as.character(uis[i]),cf,rsqr,tmin,tmax,focal.ci)
  cfs[i,]<-vec
  pred.mat<-rbind(pred.mat,new.data)
}

cfs

ggplot(gr.ests2,aes(x=temperature,y=mu))+
  geom_point()+
  geom_ribbon(data=pred.mat,aes(ymin=mu-ml.ci,ymax=mu+ml.ci),alpha=0.2)+
  geom_line(data=pred.mat)+
  geom_hline(yintercept = 0)+
  facet_wrap(~isolate.id)+
  coord_cartesian(ylim=c(-0.5,0.8))+
  theme_bw()

ggplot(gr.ests2,aes(x=temperature,y=mu))+
  geom_ribbon(data=pred.mat,aes(ymin=mu-ml.ci,ymax=mu+ml.ci,fill=sps,group=isolate.id),alpha=0.2)+
  geom_line(data=pred.mat,aes(colour=sps,group=isolate.id))+
  geom_hline(yintercept = 0)+
  coord_cartesian(ylim=c(-0.5,0.8))+
  theme_bw()

ggplot(cfs,aes(x=isolate.id,y=topt))+
  geom_point()+
  geom_errorbar(aes(ymin=topt-topt.ci,ymax=topt+topt.ci),width=0)

ggplot(cfs,aes(x=isolate.id,y=tmax))+
  geom_point()+
  geom_errorbar(aes(ymin=tmax-tmax.ci,ymax=tmax+tmax.ci),width=0.05)

cfs3<-merge(cfs,unique(gr.ests2[,c('isolate.id','sps')]))
m1<-melt(cfs3[,-1],id.vars = c('sps'))

ggplot(m1[m1$variable %in% c('topt','tmax'),],aes(x=sps,y=value))+
  geom_boxplot(aes(fill=sps))+
  facet_wrap(~variable)+
  theme_bw()

