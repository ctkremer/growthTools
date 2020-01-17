#THERMAL RXN NORM FIGURE

library(ggplot2)

#function for predicting y with x and parameters from TPC output
nbcurve<-function(x,opt,w,a,b){
  res<-a*exp(b*x)*(1-((x-opt)/(w/2))^2)
  res
}

# reparameterized norberg function for predicting growth rate as a function of temperature (x) and parameters
nbcurve2<-function(x,opt,w,a,b){
  num <- -2-2*b*opt+2*b*x+sqrt(4+(b*w)^2)
  res<-exp(a)*exp(b*x)*(1-(num/(b*w))^2)
  res
}

#example data
C6 <- data.frame(x=seq(-1,10,0.1), y=nbcurve2(x=seq(-1,10,0.1),a=log(0.123244431),b=0.033536245,w=8.621281231,opt=3.335923056))
D7 <- data.frame(x=seq(-1,10,0.1), y=nbcurve2(x=seq(-1,10,0.1),a=log(0.077125957),b=0.200952094,w=8.612147265,opt=4.768969))
E5 <- data.frame(x=seq(-1,10,0.1), y=nbcurve2(x=seq(-1,10,0.1),a=log(0.51594585),b=-0.352517187,w=10.81390632,opt=0.645300787))
F5 <- data.frame(x=seq(-1,10,0.1), y=nbcurve2(x=seq(-1,10,0.1),a=log(0.064084082),b=0.246734025,w=8.791615871,opt=5.02975128))

mu<-nbcurve2(ts,opt = tmp.nb.fits.clean$topt,w = tmp.nb.fits.clean$w,a = log(tmp.nb.fits.clean$a),b = tmp.nb.fits.clean$b)

# data2 <- data.frame(x=seq(0,50,0.1)/2.8, y=nbcurve(x=seq(0,30,0.1),16,20,0.02,0.12))

ggplot() + 
  geom_line(data=C6, aes(x,y), size=2, color="black") +
  geom_line(data=D7, aes(x,y), size=2, color="red") +
  geom_line(data=E5, aes(x,y), size=2, color="green") +
   geom_line(data=F5, aes(x,y), size=2, color="blue") +
  scale_x_continuous(limits=c(-1,10), expand=c(0,0)) +
  xlab("Temperature") +
  ylab("mu") +
  scale_y_continuous(limits=c(0,0.4), expand=c(0,0)) +
  theme_linedraw(base_size=22)

     

ggsave("rnx_norm_example.png", plot = last_plot(), device = "png", path = NULL,
       scale = 1, width = 25, height = 15, units = "cm", dpi = 300)