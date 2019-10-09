
library(dplyr)
library(mleTools)
library(growthTools)
library(bbmle)

tmp<-read.csv('./user/example_data_RI.csv')
head(tmp)

# design custom distribution, mixing binomial and normal distributions
dnormbinom<-function(x,mean,sd,p,log=FALSE){
  if(log){
    res<-ifelse(x<=0,log(p),log((1-p)*dnorm(x,mean,sd)))
  }else{
    res<-ifelse(x<=0,p,(1-p)*dnorm(x,mean,sd))
  }
  res
}

expit<-function(x){exp(x)/(1+exp(x))}


#### CH4_15_RI_03-ESAW ####

tmp1<-tmp %>% filter(isolate.id.media == 'CH4_15_RI_03-ESAW' & outlier.QC==0)
tmp1$growQ<-ifelse(tmp1$mu>0,1,0)
head(tmp1)

# convergence issues due to complete separability. something to watch for.
m1<-glm(growQ~temperature,family = binomial,data=tmp1)
summary(m1)

plot(growQ~temperature,data=tmp1)
curve(expit(296.44 + -10.96*x),0,35,col='red',add=T,n = 1000)

# full data set
ml1<-mle2(mu~dnorm(mean=nbcurve(temperature,opt,w,a,b),sd=exp(s)),start = list(opt=16,w=32,a=log(0.46),b=0.05,s=0),data=tmp1)
summary(ml1)

# full data set + dnormbinom censored
ml2<-mle2(mu~dnormbinom(mean=nbcurve(temperature,opt,w,a,b),sd=exp(s),p=expit(p0+p1*tmp1$temperature)),start = list(opt=16,w=32,a=log(0.46),b=0.05,s=0,p0=0,p1=0),data=tmp1)
summary(ml2)

# subset (Temp < 26 C)
ml3<-mle2(mu~dnorm(mean=nbcurve(temperature,opt,w,a,b),sd=exp(s)),start = list(opt=16,w=32,a=log(0.46),b=0.05,s=0),data=tmp1[tmp1$temperature<26,])
summary(ml3)

# subset (Temp < 30 C)
ml4<-mle2(mu~dnorm(mean=nbcurve(temperature,opt,w,a,b),sd=exp(s)),start = list(opt=16,w=32,a=log(0.46),b=0.05,s=0),data=tmp1[tmp1$temperature<30,])
summary(ml4)

#AICtab(ml1,ml2)

cfs1<-coef(ml1)
cfs2<-coef(ml2)
cfs3<-coef(ml3)
cfs4<-coef(ml4)

plot(mu~temperature,data=tmp1,main='CH4_15_RI_03-ESAW')
curve(nbcurve(x,topt=cfs1['opt'],w = cfs1['w'],a=cfs1['a'],b=cfs1['b']),0,35,col='blue',add=T)
curve(nbcurve(x,topt=cfs2['opt'],w = cfs2['w'],a=cfs2['a'],b=cfs2['b']),0,35,col='red',add=T)
curve(nbcurve(x,topt=cfs3['opt'],w = cfs3['w'],a=cfs3['a'],b=cfs3['b']),0,35,col='green',lty=3,add=T)
curve(nbcurve(x,topt=cfs4['opt'],w = cfs4['w'],a=cfs4['a'],b=cfs4['b']),0,35,col='purple',add=T)
legend(x = 3,y = -0.05,legend = c('all data','censored model','< 26 C','< 30 C'),col=c('blue','red','green','purple'),lty=c(1,1,3,1))
abline(0,0)

# visualize censored model:

summary(ml2)

plot(mu~temperature,data=tmp1,main='CH4_15_RI_03-ESAW',xlab='Temperature',ylab='Growth rate (1/day)')
abline(0,0,col='gray')
curve(nbcurve(x,topt=cfs2['opt'],w = cfs2['w'],a=cfs2['a'],b=cfs2['b']),0,35,col='red',add=T)
curve(expit(cfs2['p0']+cfs2['p1']*x),0,35,lty=2,col='red',add=T)
legend(x = 3,y = -0.15,legend = c('Norberg TPC','P(dying)'),col=c('red','red'),lty=c(1,2))



######  CH8_15_RI_03-ESAW  ########

tmp1<-tmp %>% filter(isolate.id.media == 'CH8_15_RI_03-ESAW' & outlier.QC==0)
tmp1$growQ<-ifelse(tmp1$mu>0,1,0)
head(tmp1)

# convergence issues due to complete separability. something to watch for.
m1<-glm(growQ~temperature,family = binomial,data=tmp1)
summary(m1)

plot(growQ~temperature,data=tmp1)
curve(expit(coef(m1)[1] + coef(m1)[2]*x),0,35,col='red',add=T,n = 1000)


## Fit models:

curve(nbcurve2(x,opt=16,w = cfs1['w'],a=cfs1['a'],b=cfs1['b']),0,35,col='blue',add=T)

ml1<-mle2(mu~dnorm(mean=nbcurve2(temperature,opt,w,a,b),sd=exp(s)),start = list(opt=12,w=24.5,a=log(0.56),b=0.05,s=0),data=tmp1)
summary(ml1)

ml2<-mle2(mu~dnormbinom(mean=nbcurve2(temperature,opt,w,a,b),sd=exp(s),p=expit(p0+p1*tmp1$temperature)),start = list(opt=16,w=32,a=log(0.46),b=0.05,s=0,p0=0,p1=0),control=list(maxit=1000),data=tmp1)
summary(ml2)

ml3<-mle2(mu~dnorm(mean=nbcurve2(temperature,opt,w,a,b),sd=exp(s)),start = list(opt=14,w=25,a=-1.5,b=0.08,s=0),data=tmp1[tmp1$temperature<25,])
summary(ml3)

ml4<-mle2(mu~dnorm(mean=nbcurve2(temperature,opt,w,a,b),sd=exp(s)),start = list(opt=16,w=32,a=log(0.46),b=0.05,s=0),data=tmp1[tmp1$temperature<26,])
summary(ml4)

#AICtab(ml1,ml2)

cfs1<-coef(ml1)
cfs2<-coef(ml2)
cfs3<-coef(ml3)
cfs4<-coef(ml4)

plot(mu~temperature,data=tmp1,main='CH8_15_RI_03-ESAW')
curve(nbcurve2(x,opt=cfs1['opt'],w = cfs1['w'],a=cfs1['a'],b=cfs1['b']),0,35,col='blue',add=T)
curve(nbcurve2(x,opt=cfs2['opt'],w = cfs2['w'],a=cfs2['a'],b=cfs2['b']),0,35,col='red',add=T)
curve(nbcurve2(x,opt=cfs3['opt'],w = cfs3['w'],a=cfs3['a'],b=cfs3['b']),0,35,col='green',lty=3,add=T)
curve(nbcurve2(x,opt=cfs4['opt'],w = cfs4['w'],a=cfs4['a'],b=cfs4['b']),0,35,col='purple',add=T)
legend(x = 21,y = 0.68,legend = c('all data','new model','<25','<26'),col=c('blue','red','green','purple'),lty=c(1,1,3,1))
abline(0,0)




######  TH7_4_RI_03-TM  ########

tmp1<-tmp %>% filter(isolate.id.media == 'TH7_4_RI_03-TM' & outlier.QC==0)
tmp1$growQ<-ifelse(tmp1$mu>0,1,0)
head(tmp1)

# convergence issues due to complete separability. something to watch for.
m1<-glm(growQ~temperature,family = binomial,data=tmp1)
summary(m1)

plot(growQ~temperature,data=tmp1)
curve(expit(coef(m1)[1] + coef(m1)[2]*x),0,35,col='red',add=T,n = 1000)


## Fit models:

curve(nbcurve2(x,opt=16,w = cfs1['w'],a=cfs1['a'],b=cfs1['b']),0,35,col='blue',add=T)

ml1<-mle2(mu~dnorm(mean=nbcurve2(temperature,opt,w,a,b),sd=exp(s)),start = list(opt=12,w=24.5,a=log(0.56),b=0.05,s=0),data=tmp1)
summary(ml1)

ml2<-mle2(mu~dnormbinom(mean=nbcurve2(temperature,opt,w,a,b),sd=exp(s),p=expit(p0+p1*tmp1$temperature)),start = list(opt=16,w=32,a=log(0.46),b=0.05,s=0,p0=0,p1=0),control=list(maxit=1000),data=tmp1)
summary(ml2)

ml3<-mle2(mu~dnorm(mean=nbcurve2(temperature,opt,w,a,b),sd=exp(s)),start = list(opt=16,w=32,a=log(0.46),b=0.05,s=0),data=tmp1[tmp1$temperature<26,])
summary(ml3)

ml4<-mle2(mu~dnorm(mean=nbcurve2(temperature,opt,w,a,b),sd=exp(s)),start = list(opt=16,w=32,a=log(0.46),b=0.05,s=0),data=tmp1[tmp1$temperature<30,])
summary(ml4)

#AICtab(ml1,ml2)

cfs1<-coef(ml1)
cfs2<-coef(ml2)
cfs3<-coef(ml3)
cfs4<-coef(ml4)

plot(mu~temperature,data=tmp1,main='TH7_4_RI_03-TM')
curve(nbcurve2(x,opt=cfs1['opt'],w = cfs1['w'],a=cfs1['a'],b=cfs1['b']),0,35,col='blue',add=T)
curve(nbcurve2(x,opt=cfs2['opt'],w = cfs2['w'],a=cfs2['a'],b=cfs2['b']),0,35,col='red',add=T)
curve(nbcurve2(x,opt=cfs3['opt'],w = cfs3['w'],a=cfs3['a'],b=cfs3['b']),0,35,col='green',lty=3,add=T)
curve(nbcurve2(x,opt=cfs4['opt'],w = cfs4['w'],a=cfs4['a'],b=cfs4['b']),0,35,col='purple',add=T)
legend(x = 0,y = 1.1,legend = c('all data','new model','<26','<30'),col=c('blue','red','green','purple'),lty=c(1,1,3,1))
abline(0,0)



######  CH30_4_RI_03-ESAW  ########

tmp1<-tmp %>% filter(isolate.id.media == 'CH30_4_RI_03-ESAW' & outlier.QC==0) %>% filter(!(temperature==0 & mu<=0))
tmp1$growQ<-ifelse(tmp1$mu>0,1,0)
head(tmp1)

# convergence issues due to complete separability. something to watch for.
m1<-glm(growQ~temperature,family = binomial,data=tmp1)
summary(m1)

plot(growQ~temperature,data=tmp1)
curve(expit(coef(m1)[1] + coef(m1)[2]*x),0,35,col='red',add=T,n = 1000)


## Fit models:

curve(nbcurve2(x,opt=16,w = cfs1['w'],a=cfs1['a'],b=cfs1['b']),0,35,col='blue',add=T)

ml1<-mle2(mu~dnorm(mean=nbcurve2(temperature,opt,w,a,b),sd=exp(s)),start = list(opt=12,w=24.5,a=log(0.56),b=0.05,s=0),data=tmp1)
summary(ml1)

ml2<-mle2(mu~dnormbinom(mean=nbcurve2(temperature,opt,w,a,b),sd=exp(s),p=expit(p0+p1*tmp1$temperature)),start = list(opt=16,w=32,a=log(0.46),b=0.05,s=0,p0=0,p1=0),control=list(maxit=1000),data=tmp1)
summary(ml2)

ml3<-mle2(mu~dnorm(mean=nbcurve2(temperature,opt,w,a,b),sd=exp(s)),start = list(opt=16,w=32,a=log(0.46),b=0.05,s=0),data=tmp1[tmp1$temperature<26,])
summary(ml3)

ml4<-mle2(mu~dnorm(mean=nbcurve2(temperature,opt,w,a,b),sd=exp(s)),start = list(opt=16,w=32,a=log(0.46),b=0.05,s=0),data=tmp1[tmp1$temperature<30,])
summary(ml4)

#AICtab(ml1,ml2)

cfs1<-coef(ml1)
cfs2<-coef(ml2)
cfs3<-coef(ml3)
cfs4<-coef(ml4)

plot(mu~temperature,data=tmp1,main='CH30_4_RI_03-ESAW')
curve(nbcurve2(x,opt=cfs1['opt'],w = cfs1['w'],a=cfs1['a'],b=cfs1['b']),0,35,col='blue',add=T)
curve(nbcurve2(x,opt=cfs2['opt'],w = cfs2['w'],a=cfs2['a'],b=cfs2['b']),0,35,col='red',add=T)
curve(nbcurve2(x,opt=cfs3['opt'],w = cfs3['w'],a=cfs3['a'],b=cfs3['b']),0,35,col='green',lty=3,add=T)
curve(nbcurve2(x,opt=cfs4['opt'],w = cfs4['w'],a=cfs4['a'],b=cfs4['b']),0,35,col='purple',add=T)
legend(x = 0,y = -0.08,legend = c('all data','new model','<26','<30'),col=c('blue','red','green','purple'),lty=c(1,1,3,1))
abline(0,0)



##### Two-sided logistic #####

curve(expit(x+10),-20,20)
curve(expit(-1*(x-5)),-10,10,add=T,col='red')

curve(expit(10+x)*expit(5-x),-20,20,col='red',add=T)



######  CH20_15_RI_03-TM  ########

tmp1<-tmp %>% filter(isolate.id.media == 'CH20_15_RI_03-TM' & outlier.QC==0)
tmp1$growQ<-ifelse(tmp1$mu>0,1,0)
head(tmp1)

plot(mu~temperature,data=tmp1,main='CH20_15_RI_03-TM')
abline(0,0)

# convergence issues due to complete separability. something to watch for.
m1<-glm(growQ~temperature+I(temperature*temperature),family = binomial,data=tmp1)
summary(m1)

plot(growQ~temperature,data=tmp1)
curve(expit(coef(m1)[1] + coef(m1)[2]*x + coef(m1)[3]*x*x),0,35,col='red',add=T,n = 1000)


## Fit models:

plot(mu~temperature,data=tmp1,main='CH20_15_RI_03-TM')
curve(nbcurve2(x,opt=19,w = 39,a=log(0.14),b=0.1),0,35,col='blue',add=T)

ml1<-mle2(mu~dnorm(mean=nbcurve2(temperature,opt,w,a,b),sd=exp(s)),start = list(opt=19,w=39,a=log(0.14),b=0.1,s=0),control=list(maxit=1000),data=tmp1)
summary(ml1)

ml2<-mle2(mu~dnormbinom(mean=nbcurve2(temperature,opt,w,a,b),sd=exp(s),p=expit(p0+p1*tmp1$temperature)),start = list(opt=16,w=32,a=log(0.46),b=0.05,s=0,p0=0,p1=0),control=list(maxit=1000),data=tmp1)
summary(ml2)

ml3<-mle2(mu~dnorm(mean=nbcurve2(temperature,opt,w,a,b),sd=exp(s)),start = list(opt=19,w=39,a=log(0.14),b=0.1,s=0),data=tmp1[tmp1$temperature<26,])
summary(ml3)

ml4<-mle2(mu~dnorm(mean=nbcurve2(temperature,opt,w,a,b),sd=exp(s)),start = list(opt=19,w=39,a=log(0.14),b=0.1,s=0),data=tmp1[tmp1$temperature<30,])
summary(ml4)

ml5<-mle2(mu~dnormbinom(mean=nbcurve2(temperature,opt,w,a,b),sd=exp(s),p=expit(p0+p1*tmp1$temperature+p2*tmp1$temperature*tmp1$temperature)),start = list(opt=16,w=32,a=log(0.46),b=0.05,s=0,p0=-65.9,p1=25.2,p2=-0.9),control=list(maxit=1000),data=tmp1)
summary(ml5)

plot(growQ~temperature,data=tmp1)
curve(expit(1*x-3),0,35,add=T)
curve(expit(-1*(x-26)),0,35,add=T,col='red')
curve(expit())

curve(1-expit(coef(m1)[1] + coef(m1)[2]*x + coef(m1)[3]*x*x),0,35,col='red',add=T,n = 1000)


curve(expit(cfs5['p0']+cfs5['p1']*x)*expit(-1*(cfs5['p0b']+cfs5['p1b']*x)),0,35,col='red',ylim=c(0,1))

#AICtab(ml1,ml2)

cfs1<-coef(ml1)
cfs2<-coef(ml2)
cfs3<-coef(ml3)
cfs4<-coef(ml4)
cfs5<-coef(ml5)

plot(mu~temperature,data=tmp1,main='CH30_4_RI_03-ESAW')
curve(nbcurve2(x,opt=cfs1['opt'],w = cfs1['w'],a=cfs1['a'],b=cfs1['b']),0,35,col='blue',add=T)
curve(nbcurve2(x,opt=cfs2['opt'],w = cfs2['w'],a=cfs2['a'],b=cfs2['b']),0,35,col='red',add=T)
curve(nbcurve2(x,opt=cfs3['opt'],w = cfs3['w'],a=cfs3['a'],b=cfs3['b']),0,35,col='green',lty=3,add=T)
curve(nbcurve2(x,opt=cfs4['opt'],w = cfs4['w'],a=cfs4['a'],b=cfs4['b']),0,35,col='purple',add=T)
curve(nbcurve2(x,opt=cfs5['opt'],w = cfs5['w'],a=cfs5['a'],b=cfs5['b']),0,35,col='orange',add=T)

legend(x = 0,y = -0.25,legend = c('all data','new model','<26','<30'),col=c('blue','red','green','purple'),lty=c(1,1,3,1))
abline(0,0)


