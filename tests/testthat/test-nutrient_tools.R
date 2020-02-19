context("Testing Monod curve regression tools")
library(growthTools)

test_that("Monod fitting in detail", {
  monod.tmp<-data.frame(nutrients=c(0.1,0.5,1,2,5,10),mu=c(0.02,0.1,0.2,0.4,0.5,0.52))
  #plot(mu~nutrients,data=dat)
  
  umax.guess <- max(monod.tmp$mu)[1]
  k.guess <- max(monod.tmp$nutrients)[1]/3
  z.guess <- 0
  
  fit<-bbmle::mle2(mu~dnorm(mean=monod_curve(nutrients,umax,k,z),sd=exp(s)),
                   start=list(umax=umax.guess,k=k.guess,z=z.guess,s=log(2)),
                   fixed=list(z=0),data=monod.tmp)
  
  expect_equal(coef(fit)[[1]],0.650484,tolerance=1E-7)
})

test_that("Monod fitting with get.monod, fixed intercept", {
  dat<-data.frame(nutrients=c(0.1,0.5,1,2,5,10),mu=c(0.02,0.1,0.2,0.4,0.5,0.52))
  #plot(mu~nutrients,data=dat)
  
  fit.monod<-get.monod(nutrients=dat$nutrients,mu=dat$mu,method='mle2',fix_intercept = T)
  cis<-fit.monod$cf_ciFI
  #fit.monod$cf
  #plot(fit.monod)
  
  tmp.cis<-matrix(c(0.5510411,1.0216724,-3.8929826,
                    0.7499269,2.5946336,-2.7613688),nrow=3,ncol=2)
  attr(tmp.cis,"dimnames")<-list(c("umax","k","s"),c("2.5 %","97.5 %"))
  
  #fit.monod.B<-get.monod(nutrients=dat$nutrients,mu=dat$mu,method='mle2',fix_intercept = F)
  #plot(fit.monod.B)
  
  expect_equal(fit.monod$cf$umax,0.650484,tolerance=1E-7)
  expect_equal(fit.monod$cf$k,1.808153,tolerance=1E-7)
  expect_equal(fit.monod$cf$z,0,tolerance=1E-7)
  expect_equal(fit.monod$cf$s,-3.327176,tolerance=1E-7)
  
  expect_equal(fit.monod$rsqr,0.9658253,tolerance=1E-7)
  expect_equal(as.numeric(fit.monod$logLik),11.44946,tolerance=1E-5)
  expect_equal(fit.monod$aic,-16.89892,tolerance=1E-7)
  expect_equal(cis,tmp.cis,tolerance=1E-7)
})

test_that("Monod fitting with get.monod, estimated intercept", {
  dat<-data.frame(nutrients=c(0.1,0.5,1,2,5,10),mu=c(0.02,0.1,0.2,0.4,0.5,0.52))
  #plot(mu~nutrients,data=dat)
  
  fit.monod<-get.monod(nutrients=dat$nutrients,mu=dat$mu,method='mle2',fix_intercept = F)
  cis<-fit.monod$ciFI
  #plot(fit.monod)
  
  tmp.cis<-matrix(c(0.5863497,0.6995717,-0.1187275,-4.0116577,
                    0.76911820,2.15870318,0.02728491,-2.88005030),nrow=4,ncol=2)
  attr(tmp.cis,"dimnames")<-list(c("umax","k","z","s"),c("2.5 %","97.5 %"))
  
  expect_equal(fit.monod$cf$umax,0.6777339,tolerance=1E-7)
  expect_equal(fit.monod$cf$k,1.429137,tolerance=1E-7)
  expect_equal(fit.monod$cf$z,-0.04572128,tolerance=1E-7)
  expect_equal(fit.monod$cf$s,-3.445854,tolerance=1E-7)
  
  expect_equal(fit.monod$rsqr,0.9730458,tolerance=1E-7)
  expect_equal(as.numeric(fit.monod$logLik),12.1615,tolerance=1E-5)
  expect_equal(fit.monod$aic,-16.323,tolerance=1E-7)
  expect_equal(cis,tmp.cis,tolerance=1E-7)
})
