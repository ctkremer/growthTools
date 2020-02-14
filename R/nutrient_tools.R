


#' Fit Monod curve to growth rate vs. nutrient concentration data
#' 
#' Note: this function currently does not use grid.mle2, as the regressions are usually pretty stable; add this 
#' feature later.
#' 
#' @param nutrients Nutrient concentration
#' @param mu Exponential growth rate
#' @param method Specify which fitting algorithm to use, 'mle2' or 'grid.mle2'
#' @param ... Additional arguments passed to grid.mle2 (e.g., control=list(maxit=2000))
#' 
#' @examples 
#' 
#' @export
#' @import emdbook
#' @import ggplot2
get.monod<-function(nutrients,mu,method='mle2',...){
  monod.tmp<-stats::na.omit(data.frame(nutrients,mu))
  nnutr<-length(unique(monod.tmp$nutrients))
  
  if(nnutr<=3){
    print("Caution in get.monod - focal data set has <=3 unique nutrients, risk of overfitting is high!")
  }
  id<-id[1]
  
  if(method=='grid.mle2'){
    print("Caution: this option not fully beta tested... 2/13/20")
  }
  
  if(method=='mle2'){
    
    umax.guess <- max(monod.tmp$mu)[1]
    k.guess <- max(monod.tmp$nutrients)[1]/3
    fit<-bbmle::mle2(mu~dnorm(mean=umax*nutrients/(nutrients+k),sd=exp(s)),
                     start=list(umax=umax.guess,k=k.guess,s=log(2)),data=monod.tmp)
  }
  
  # pull out parameters
  cf<-as.list(coef(fit))
  vcov.mat<-vcov(fit)
  
  # calculate R2
  rsqr<-get.R2(predict(fit),monod.tmp$mu)
  
  # simple Fisher confidence intervals:
  ciFI<-mleTools::ci.FI(fit)
  
  # save output:
  vec<-as.list(c(type='monod'))
  vec$cf<-cf
  vec$ciFI<-ciFI
  vec$vcov<-vcov.mat
  vec$nobs<-nrow(tpc.tmp)
  vec$ntemps<-ntemps
  vec$rsqr<-rsqr
  vec$logLik<-logLik(fit)
  vec$aic<-stats::AIC(fit)
  vec$data<-monod.tmp
  
  # Finished, return relevant stats.  
  return(vec)
}
