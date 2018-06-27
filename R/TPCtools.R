

##########################################

# Conduct actual Norberg fits:

##########################################

library(emdbook)
source("/Users/colin/Teaching/ELME/ELME 2017/labs/Lab4 - strategies for tough fits/grid.mle2 053114.R")

# Jump in with growth rate fitting:
head(gr.ests2)

# Preliminary plotting:
ggplot(data.frame(gr.ests2[gr.ests2$mu<0.8,]),aes(x=temperature,y=mu))+
  geom_point(aes(shape=replicate))+
  geom_hline(yintercept = 0)+
  stat_smooth(method='gam',formula=y~s(x,k=4))+
  facet_wrap(~isolate.id,scales='free_y')+
  theme_bw()


###########

confint.FI<-function(model){
  cfs<-coef(model)
  ses<-sqrt(diag(vcov(model)))	# standard errors
  lw<-cfs-1.96*ses
  up<-cfs+1.96*ses
  res<-cbind(lw,up)
  dimnames(res)<-list(names(cfs),c("2.5 %","97.5 %"))	
  res
}

# Norberg
nbcurve<-function(x,opt,w,a,b){
  #  res<-exp(b*x+a)*(1-((x-exp(opt))/(10*w/2))^2)
  res<-exp(b*x+a)*(1-((x-opt)/(w/2))^2)
  res
}

# re-parameterized so that opt is explicitly the growth optima:
nbcurve<-function(x,opt,w,a,b){
  num <- -2-2*b*opt+2*b*x+sqrt(4+(b*w)^2)
  res<-exp(a+b*x)*(1-(num/(b*w))^2)
  res
}

get.topt<-function(topt,w,a,b){
  (1/(2*b))*(-2+2*b*topt+sqrt(4+b^2*w^2))
}

get.tlim<-function(topt,w,a,b,type='tmin'){
  if(!(type %in% c('tmin','tmax'))) print('Error in get.tlim; unknown type requested')
  
  if(type=='tmax'){
    res <- (2+2*b*topt+b*w-sqrt(4+b^2*w^2))/(2*b)
  }
  if(type=='tmin'){
    res <- (2+2*b*topt-b*w-sqrt(4+b^2*w^2))/(2*b)
  }
  return(res)
}

# May need to implement grid searching here, which is really going to slow things down, alas.

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




