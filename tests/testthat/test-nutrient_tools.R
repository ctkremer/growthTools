context("Testing Monod curve regression tools")
library(growthTools)

test_that("Monod fitting in detail", {
  monod.tmp<-data.frame(nutrients=c(0.1,0.5,1,2,5,10),mu=c(0.02,0.1,0.2,0.4,0.5,0.52))
  #plot(mu~nutrients,data=monod.tmp)
  
  k.guess <- max(monod.tmp$nutrients)[1]/3
  umax.guess <- max(c(max(monod.tmp$mu)[1],0.01))
  
  fit<-bbmle::mle2(mu~dnorm(mean=monod_curve(nutrients,umax,k,log.k=TRUE,log.umax=TRUE),sd=exp(s)),
                   start=list(umax=log(umax.guess),k=log(k.guess),s=log(2)),data=monod.tmp)
  
  expect_equal(coef(fit)[[1]],-0.4300399,tolerance=1E-7)
})

test_that("Monod type 2 fitting in detail", {
  monod.tmp<-data.frame(nutrients=c(0.1,0.5,1,2,5,10),mu=c(0.02,0.1,0.2,0.4,0.5,0.52))
  #plot(mu~nutrients,data=monod.tmp)
  
  k.guess <- max(monod.tmp$nutrients)[1]/3
  z.guess <- max(c(0.000001,abs(min(monod.tmp$mu)[1])))
  umax.guess <- max(c(max(monod.tmp$mu)[1],0.01))+z.guess
  
  fit<-bbmle::mle2(mu~dnorm(mean=monod_curve_type2(nutrients,umax,k,z,log.pars=TRUE),sd=exp(s)),
                   start=list(umax=log(umax.guess),k=log(k.guess),z=log(z.guess),s=log(2)),
                   data=monod.tmp)
  
  expect_equal(coef(fit)[[1]],-0.3890017,tolerance=1E-7)
})


test_that("Monod type 2 fitting in detail", {
  monod.tmp<-data.frame(nutrients=c(0.1,0.5,1,2,5,10),mu=c(0.02,0.1,0.2,0.4,0.5,0.52)-0.1)
  #plot(mu~nutrients,data=monod.tmp)
  
  k.guess <- max(monod.tmp$nutrients)[1]/3
  z.guess <- max(c(0.000001,abs(min(monod.tmp$mu)[1])))
  umax.guess <- max(c(max(monod.tmp$mu)[1],0.01))+z.guess
  
  fit<-bbmle::mle2(mu~dnorm(mean=monod_curve_type2(nutrients,umax,k,z,log.pars=TRUE),sd=exp(s)),
                   start=list(umax=log(umax.guess),k=log(k.guess),z=log(z.guess),s=log(2)),
                   data=monod.tmp)
  
  # reframe result on non-logged scale
  model_cf<-as.list(coef(fit))
  tmp.cfs<-model_cf
  tmp.cfs$umax<-exp(tmp.cfs$umax)
  tmp.cfs$k<-exp(tmp.cfs$k)
  tmp.cfs$z<-exp(tmp.cfs$z)
  fit0<-bbmle::mle2(mu~dnorm(mean=monod_curve_type2(nutrients,umax,k,z,log.pars=FALSE),sd=exp(s)),
                    start=tmp.cfs,control=list(maxit=0),
                    data=monod.tmp)
  bbmle::vcov(fit)
  bbmle::vcov(fit0)
  bbmle::summary(fit0)
  ?bbmle::vcov
  
  expect_equal(coef(fit)[[1]],-0.3890017,tolerance=1E-7)
})

test_that("Monod fitting with get.monod, intercept fixed at 0", {
  dat<-data.frame(nutrients=c(0.1,0.5,1,2,5,10),mu=c(0.02,0.1,0.2,0.4,0.5,0.52))
  #plot(mu~nutrients,data=dat)
  
  fit.monod<-get.monod(nutrients=dat$nutrients,mu=dat$mu,method='mle2')
  cis<-fit.monod$display_ci
  #fit.monod$display_cf
  #plot(fit.monod)
  
  tmp.cis<-matrix(c(0.551041,1.021674,-3.892985,
                    0.7499253,2.5946229,-2.7613773),nrow=3,ncol=2)
  attr(tmp.cis,"dimnames")<-list(c("umax","k","s"),c("2.5 %","97.5 %"))
  
  expect_equal(fit.monod$display_cf$umax,0.6504832,tolerance=1E-7)
  expect_equal(fit.monod$display_cf$k,1.808149,tolerance=1E-7)
  expect_equal(fit.monod$display_cf$s,-3.327181,tolerance=1E-7)
  
  expect_equal(fit.monod$rsqr,0.9658253,tolerance=1E-7)
  expect_equal(as.numeric(fit.monod$logLik),11.44946,tolerance=1E-5)
  expect_equal(fit.monod$aic,-16.89892,tolerance=1E-7)
  expect_equal(cis,tmp.cis,tolerance=1E-7)
})

test_that("Monod fitting with get.monod2, estimating z, true intercept ~=0", {
  dat<-data.frame(nutrients=c(0.1,0.5,1,2,5,10),mu=c(0.02,0.1,0.2,0.4,0.5,0.52))
  #plot(mu~nutrients,data=dat)
  
  fit.monod2<-get.monod2(nutrients=dat$nutrients,mu=dat$mu,method='mle2')
  cis<-fit.monod2$display_ci
  #plot(fit.monod2)
  
  tmp.cis<-matrix(c(0.5510423,1.0216824,NaN,-3.8929874,
                    0.7499246,2.5945958,NaN,-2.7613863),nrow=4,ncol=2)
  attr(tmp.cis,"dimnames")<-list(c("umax","k","z","s"),c("2.5 %","97.5 %"))
  
  expect_equal(fit.monod2$umax,0.6504835,tolerance=1E-7)
  expect_equal(fit.monod2$k,1.808139,tolerance=1E-7)
  expect_equal(exp(fit.monod2$cf$z),1e-06,tolerance=1E-7)
  expect_equal(fit.monod2$cf$s,-3.327187,tolerance=1E-7)
  
  expect_equal(fit.monod2$rsqr,0.9658256,tolerance=1E-7)
  expect_equal(as.numeric(fit.monod2$logLik),11.44949,tolerance=1E-5)
  expect_equal(fit.monod2$aic,-14.89898,tolerance=1E-7)
  expect_equal(cis,tmp.cis,tolerance=1E-7)
})

test_that("Monod fitting with get.monod, estimating z, true intercept below 0", {
  dat<-data.frame(nutrients=c(0.1,0.5,1,2,5,10),mu=c(0.02,0.1,0.2,0.4,0.5,0.52)-0.1)
  #plot(mu~nutrients,data=dat)
  
  fit.monod2<-get.monod2(nutrients=dat$nutrients,mu=dat$mu,method='mle2')
  cis<-fit.monod2$display_ci
  #plot(fit.monod2)
  #fit.monod2
  
  tmp.cis<-matrix(c(0.5510423,1.0216824,NaN,-3.8929874,
                    0.7499246,2.5945958,NaN,-2.7613863),nrow=4,ncol=2)
  attr(tmp.cis,"dimnames")<-list(c("umax","k","z","s"),c("2.5 %","97.5 %"))
  
  expect_equal(fit.monod$umax,0.6504835,tolerance=1E-7)
  expect_equal(fit.monod$k,1.808139,tolerance=1E-7)
  expect_equal(exp(fit.monod$cf$z),1e-06,tolerance=1E-7)
  expect_equal(fit.monod$cf$s,-3.327187,tolerance=1E-7)
  
  expect_equal(fit.monod$rsqr,0.9658256,tolerance=1E-7)
  expect_equal(as.numeric(fit.monod$logLik),11.44949,tolerance=1E-5)
  expect_equal(fit.monod$aic,-14.89898,tolerance=1E-7)
  expect_equal(cis,tmp.cis,tolerance=1E-7)
})
