context("Testing Thermal Performance Curve regression tools")
library(growthTools)

test_that("internals of get.nbcurve.tpc work",{
  temperature<-c(0, 0, 0, 0, 10, 10, 10, 10, 16, 16, 16, 16, 2, 2, 2, 2, 20,
                20, 20, 20, 25, 25, 25, 25, 29, 29, 29, 29, 4, 4, 4, 4)
  mu<-c(0.29570190, 0.31715749, 0.23086634, 0.26492142, 0.34044182, 0.31923753,
        0.30635082, 0.31520954, 0.54979405, 0.50926419, 0.46722819, 0.71566211,
        0.12797674, 0.18134072, 0.15040029, 0.20437558, 0.73161025, 0.72734410,
        0.60578314, 0.65573113, 0.54062615, 0.34272669, 0.71122612, 0.23908063,
        -0.07640984, -0.07830364, -0.08734117, -0.07812054, 0.24811323,
        0.33837116, 0.24791146, 0.27200425)
  
  tpc.tmp<-stats::na.omit(data.frame(temperature,mu))
  ntemps<-length(unique(tpc.tmp$temperature))
  
  # set up search of a grid of parameter guesses
  grids<-list(topt=seq(15,35,5),w=seq(10,40,5),a=seq(-0.5,-3,-0.5),b=c(-0.05,0,0.05))
  start<-list(topt=NA,w=NA,a=NA,b=NA,s=log(2))
  
  fit0<-suppressWarnings(mleTools::grid.mle2(minuslogl=mu~dnorm(mean=nbcurve(temperature,topt,w,a,b),sd=exp(s)),
                                   grids=grids,start=start,data=tpc.tmp))
  
  cfg<-bbmle::coef(fit0$res.best) # this seemed to be throwing problems b/c of an issue with accessing mle2...?
  target.cfg<-c(20.3754421,104.1459933,-1.4929778,0.1109893,-2.3665936)
  attr(target.cfg,"names")<-c("topt","w","a","b","s")
  
  expect_equal(cfg,target.cfg,tolerance=1E-6)
})

test_that("get.nbcurve.tpc works as expected",{
  
  tmp<-data.frame(temperature=c(0, 0, 0, 0, 10, 10, 10, 10, 16, 16, 16, 16, 2, 2, 2, 2, 20,
                                20, 20, 20, 25, 25, 25, 25, 29, 29, 29, 29, 4, 4, 4, 4), 
                  mu=c(0.29570190, 0.31715749, 0.23086634, 0.26492142, 0.34044182, 0.31923753,
                       0.30635082, 0.31520954, 0.54979405, 0.50926419, 0.46722819, 0.71566211,
                       0.12797674, 0.18134072, 0.15040029, 0.20437558, 0.73161025, 0.72734410,
                       0.60578314, 0.65573113, 0.54062615, 0.34272669, 0.71122612, 0.23908063,
                       -0.07640984, -0.07830364, -0.08734117, -0.07812054, 0.24811323,
                       0.33837116, 0.24791146, 0.27200425))
  
  nbts<-get.nbcurve.tpc(tmp$temperature,tmp$mu,method='grid.mle2',plotQ=FALSE,
                                  conf.bandQ = FALSE,fpath=NA)
  cis<-nbts$ciF
  
  # target confidence interval values
  tmp.cis<-matrix(c(19.51491008,104.14532361,-1.71386004,0.09852027,-2.61171380,
                    21.2359741,104.1466630,-1.2720955,0.1234584,-2.1214734),nrow=5,ncol=2)
  attr(tmp.cis,"dimnames")<-list(c("topt","w","a","b","s"),c("2.5 %","97.5 %"))
  
  expect_equal(nbts$topt,20.3754421,tolerance=1E-7)
  expect_equal(nbts$w,104.1459933,tolerance=1E-7)
  expect_equal(nbts$a,-1.4929778,tolerance=1E-7)
  expect_equal(nbts$b,0.1109893,tolerance=1E-7)
  expect_equal(nbts$s,-2.3665936,tolerance=1E-7)
  expect_equal(nbts$rsqr,0.8405597,tolerance=1E-7)
  expect_equal(nbts$tmin,-75.5343893,tolerance=1E-7)
  expect_equal(nbts$tmax,28.6116040,tolerance=1E-7)
  expect_equal(nbts$umax,0.6282159,tolerance=1E-7)
  expect_equal(nbts$aic,-50.68104,tolerance=1E-7)
  expect_equal(cis,tmp.cis,tolerance=1E-7)
})


test_that("get.decurve.tpc works as expected",{
  
  tmp<-data.frame(temperature=c(0, 0, 0, 0, 10, 10, 10, 10, 16, 16, 16, 16, 2, 2, 2, 2, 20,
                                20, 20, 20, 25, 25, 25, 25, 29, 29, 29, 29, 4, 4, 4, 4), 
                  mu=c(0.29570190, 0.31715749, 0.23086634, 0.26492142, 0.34044182, 0.31923753,
                       0.30635082, 0.31520954, 0.54979405, 0.50926419, 0.46722819, 0.71566211,
                       0.12797674, 0.18134072, 0.15040029, 0.20437558, 0.73161025, 0.72734410,
                       0.60578314, 0.65573113, 0.54062615, 0.34272669, 0.71122612, 0.23908063,
                       -0.07640984, -0.07830364, -0.08734117, -0.07812054, 0.24811323,
                       0.33837116, 0.24791146, 0.27200425))
  
  dets<-suppressWarnings(get.decurve.tpc(tmp$temperature,tmp$mu,method='grid.mle2',plotQ=FALSE,
                        conf.bandQ = FALSE,fpath=NA))
  cis<-dets$ciF
  
  # target confidence interval values
  tmp.cis<-matrix(c(19.30674384,NaN,NaN,NaN,NaN,0.07115396,
                    22.1291309,NaN,NaN,NaN,NaN,0.1153003),nrow=6,ncol=2)
  attr(tmp.cis,"dimnames")<-list(c("topt","b1","b2","d0","d2","s"),c("2.5 %","97.5 %"))
  
  expect_equal(dets$topt,20.71794,tolerance=1E-5)
  expect_equal(dets$b1,0.2036075,tolerance=1E-7)
  expect_equal(dets$b2,0.08491238,tolerance=1E-7)
  expect_equal(dets$d0,2.73132e-07,tolerance=1E-7)
  expect_equal(dets$d2,0.1800175,tolerance=1E-7)
  expect_equal(dets$s,0.09322711,tolerance=1E-7)
  expect_equal(dets$rsqr,0.8423461,tolerance=1E-7)
  expect_equal(dets$tmin,-159.2436,tolerance=1E-5)
  expect_equal(dets$tmax,28.61903,tolerance=1E-5)
  expect_equal(dets$umax,0.6247376,tolerance=1E-7)
  expect_equal(dets$aic,-49.04162,tolerance=1E-7)
  expect_equal(cis,tmp.cis,tolerance=1E-7)
})


