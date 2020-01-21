context("Testing growth rate estimation tools")
library(growthTools)

test_that("square root smoother function", {
  expect_equal(sqfunc(0,1,0), 0)
  expect_equal(sqfunc(0,1,1), 1)
  expect_equal(sqfunc(1,1,1), 1.118034)
})

test_that("piecewise linear functions behave",{
  expect_equal(lag(0,5,1,4,s=1E-10), 5)
  expect_equal(sat(0,0.9,1,8,s=1E-10), 0.90001)
  expect_equal(lagsat(0,5.1,1,4,8,s=1E-10), 5.1)
  expect_equal(flr(0,10,-1,8,s=1E-10), 9.99999)
  
  expect_equal(lag(10,5,1,4,s=1E-10), 11)
  expect_equal(sat(10,0.9,1,8,s=1E-10), 8.9)
  expect_equal(lagsat(10,5.1,1,4,8,s=1E-10), 9.1)
  expect_equal(flr(10,10,-1,8,s=1E-10), 2)
})

test_that("growth rate estimates are stable",{
  
  # Construct example data set:
  sdat<-data.frame(trt=c(rep('A',10),rep('B',10),rep('C',10),rep('D',10)),
                   dtime=rep(seq(1,10),4),
                   ln.fluor=c(c(1,1.1,0.9,1,2,3,4,5,5.2,4.7),
                              c(1.1,0.9,1,2,3,4,4.1,4.2,3.7,4)+0.3,
                              c(3.5,3.4,3.6,3.5,3.2,2.2,1.2,0.5,0.4,0.1),
                              c(5.5,4.5,3.5,2.5,1.5,0,0.2,-0.1,0,-0.1)))
  sdatA<-sdat[sdat$trt=='A',]
  sdatD<-sdat[sdat$trt=='D',]
  
  # calculate growth rate using all available methods:
  resA<-get.growth.rate(sdatA$dtime,sdatA$ln.fluor,plot.best.Q = F,id = 'Population A')
  resD<-get.growth.rate(sdatD$dtime,sdatD$ln.fluor,plot.best.Q = F,id = 'Population D')
  
  
  # format expectation objects
  target.slopesA<-c(0.5606061,0.6964286,0.5961677,1.0000030)
  attr(target.slopesA,"names")<-c("gr","gr.lag","gr.sat","gr.lagsat")
  target.sesA<-c(0.07157291,0.10093405,0.09890251,0.10033133)
  attr(target.sesA,"names")<-c("gr","gr.lag","gr.sat","gr.lagsat")
  
  target.slopesD<-c(-0.6563636, -1.0714301)
  attr(target.slopesD,"names")<-c("gr","gr.flr")
  target.sesD<-c(0.08827393, 0.03823087)
  attr(target.sesD,"names")<-c("gr","gr.flr")
  
  # expectations to test
  expect_match(resA$best.model, "gr.lagsat")
  expect_equal(resA$slopes, target.slopesA, tolerance=1E-7, scale=1)
  expect_equal(resA$ses, target.sesA, tolerance=1E-7, scale=1)
  expect_equal(resA$best.model.rsqr, 0.9949958, tolerance=1E-7, scale=1)
  
  expect_match(resD$best.model, "gr.flr")
  expect_equal(resD$slopes, target.slopesD, tolerance=1E-7, scale=1)
  expect_equal(resD$ses, target.sesD, tolerance=1E-7, scale=1)
  expect_equal(resD$best.model.rsqr, 0.9955992, tolerance=1E-7, scale=1)
})

