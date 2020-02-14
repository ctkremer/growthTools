context("Testing Monod curve regression tools")
library(growthTools)

test_that("Monod fitting", {
  dat<-data.frame(nutrients=c(0.1,0.5,1,2,5,10),mu=c(0.02,0.1,0.2,0.4,0.5,0.52))
  plot(mu~nutrients,data=dat)
  
  fit.monod<-get.monod(nutrients=dat$nutrients,mu=dat$mu,method='mle2')
  
  #expect_equal(sqfunc(0,1,0), 0)
  #expect_equal(sqfunc(0,1,1), 1)
  #expect_equal(sqfunc(1,1,1), 1.118034)
})
