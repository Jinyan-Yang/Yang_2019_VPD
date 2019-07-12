# function to get BIC value for each model
get.bic.func <- function(model.vec,data.vec,n.fd){
  d.vec <- model.vec
  m.vec <- data.vec
  
  # resi.vec <- model.vec - data.vec
  # d.sd <- sd(data.vec,na.rm=TRUE)
  # obs.n <- length(data.vec)
  # 
  # w <- rep(1,obs.n)
  # # 0.5 * (sum(log(w)) - n * (log(2 * pi) + 1 - log(n) + log(sum(w * res^2))))
  # df.ll <- 0.5 * (sum(log(w)) - obs.n * (log(2 * pi) + 1 - log(obs.n) + log(sum(w * (resi.vec)^2))))
  # # xlogπ+(n−x)log(1−π)
  # # -2 * ll + log(n) * df.ll
  # bic <- -2*df.ll + log(obs.n) * (n.fd + 1)
  # 
  
  
  res <- model.vec - data.vec
  n <- length(data.vec)  
  w <- rep(1,n) #not applicable
  
  ll<-0.5 * (sum(log(w)) - n * (log(2 * pi) + 1 - log(n) + log(sum(w * res^2))))
  # Cll-logLik(m)==0 #TRUE
  k.original<-n.fd
  df.ll<-k.original+1 
  bic<- -2 * ll + log(n) * df.ll
  
  return(bic)
  # -0.5 *  sum((m.vec - d.vec)^2)/d.sd - n.fd*log(length(m.vec))
}

# function to get BIC value for each model
get.bic.both.func <- function(model.a,model.gs,data.a,data.gs,n.fd){
  
  a.norm <- (model.a - data.a) / sd(data.a,na.rm=TRUE)
  gs.norm <- (model.gs - data.gs) / sd(data.gs,na.rm=TRUE)
  
  # resi.vec <- model.vec - data.vec
  # d.sd <- sd(data.vec,na.rm=TRUE)
  # obs.n <- length(data.vec)
  # 
  # w <- rep(1,obs.n)
  # # 0.5 * (sum(log(w)) - n * (log(2 * pi) + 1 - log(n) + log(sum(w * res^2))))
  # df.ll <- 0.5 * (sum(log(w)) - obs.n * (log(2 * pi) + 1 - log(obs.n) + log(sum(w * (resi.vec)^2))))
  # # xlogπ+(n−x)log(1−π)
  # # -2 * ll + log(n) * df.ll
  # bic <- -2*df.ll + log(obs.n) * (n.fd + 1)
  # 
  
  
  res <- abs(a.norm) + abs(gs.norm)
  n <- length(a.norm)  
  w <- rep(1,n) #not applicable
  
  ll <- 0.5 * (sum(log(w)) - n * (log(2 * pi) + 1 - log(n) + log(sum(w * res^2))))
  # Cll-logLik(m)==0 #TRUE
  k.original<-n.fd
  df.ll<-k.original+1 
  bic<- -2 * ll + log(n) * df.ll
  
  return(bic)
  # -0.5 *  sum((m.vec - d.vec)^2)/d.sd - n.fd*log(length(m.vec))
}

get.r.func <- function(x.vec,y.vec){
  # x.vec = optbb.spots$ALEAF
  # y.vec=spot.amb$Photo
  fit.lm <- lm(y.vec~x.vec)
  return(summary(fit.lm)$r.squared)
}
# get bic for photo
get.bic.both.func(optbb.spots$ALEAF,spot.amb$Photo,optbb.spots$GS,spot.amb$Cond,n.fd=2)
get.bic.both.func(opt.d$ALEAF,spot.amb$Photo,opt.d$GS,spot.amb$Cond,n.fd=3)
get.bic.both.func(leuning.df$ALEAF,spot.amb$Photo,leuning.df$GS,spot.amb$Cond,n.fd=4)
get.bic.both.func(tz.spots$ALEAF,spot.amb$Photo,tz.spots$GS,spot.amb$Cond,n.fd=4)
get.bic.both.func(tz.psi.spots$ALEAF,spot.amb$Photo,tz.psi.spots$GS,spot.amb$Cond,n.fd=5)

get.bic.both.func(leuning.cable.df$ALEAF,spot.amb$Photo,leuning.cable.df$GS,spot.amb$Cond,n.fd=3)
# get.bic.both.func(tz.d.spots$ALEAF,spot.amb$Photo,tz.d.spots$GS,spot.amb$Cond,n.fd=6)
# sum(abs((optbb.spots$ALEAF-spot.amb$Photo)/sd(spot.amb$Photo)) + 
# abs((optbb.spots$GS-spot.amb$Cond)/sd(spot.amb$Cond)))
# log(188)
get.bic.func(optbb.spots$ALEAF,spot.amb$Photo,n.fd=2)
get.r.func(optbb.spots$ALEAF,spot.amb$Photo)
get.bic.func(opt.d$ALEAF,spot.amb$Photo,n.fd=3)
get.r.func(opt.d$ALEAF,spot.amb$Photo)

get.bic.func(leuning.df$ALEAF,spot.amb$Photo,n.fd=3)
get.r.func(leuning.df$ALEAF,spot.amb$Photo)
get.bic.func(leuning.cable.df$ALEAF,spot.amb$Photo,n.fd=3)

get.bic.func(tz.spots$ALEAF,spot.amb$Photo,n.fd=4)
get.r.func(tz.spots$ALEAF,spot.amb$Photo)
get.bic.func(tz.d.spots$ALEAF,spot.amb$Photo,n.fd=5)
get.r.func(tz.d.spots$ALEAF,spot.amb$Photo)
get.bic.func(tz.psi.spots$ALEAF,spot.amb$Photo,n.fd=6)
get.r.func(tz.psi.spots$ALEAF,spot.amb$Photo)
# get bic for gs
get.bic.func(optbb.spots$GS*100,spot.amb$Cond*100,n.fd=2)
get.r.func(optbb.spots$GS*100,spot.amb$Cond*100)
get.bic.func(opt.d$GS*100,spot.amb$Cond*100,n.fd=3)
get.r.func(opt.d$GS*100,spot.amb$Cond*100)

get.bic.func(leuning.df$GS*100,spot.amb$Cond*100,n.fd=3)
get.r.func(leuning.df$GS*100,spot.amb$Cond*100)
get.bic.func(leuning.cable.df$GS*100,spot.amb$Cond*100,n.fd=3)

get.bic.func(tz.spots$GS*100,spot.amb$Cond*100,n.fd=4)
get.r.func(tz.spots$GS*100,spot.amb$Cond*100)
get.bic.func(tz.d.spots$GS*100,spot.amb$Cond*100,n.fd=5)
get.r.func(tz.d.spots$GS*100,spot.amb$Cond*100)
get.bic.func(tz.psi.spots$GS*100,spot.amb$Cond*100,n.fd=6)
get.r.func(tz.psi.spots$GS*100,spot.amb$Cond*100)


# # 
# fit.test <- lm(seq(1,10,0.1)~c(2*seq(1,10,0.1)))
# 
# 
# fit.out <- lm(spot.amb.nov.2013$Cond~opt.d$GS)
# abline(fit.out,col="grey80")
# mylabel <- bquote(italic(R)^2 == .(format(summary(fit.out)$r.squared, digits = 2)))
# 
# BIC(fit.test)
# spot.amb$Photo - tz.psi.spots$ALEAF + res
# y<-spot.amb$Photo
# x<-tz.psi.spots$ALEAF
# m<-lm(y ~ x) 
# BIC(m)
# get.bic.func(y,x,2)
# see <- logLik(m)
# # summary(lm(spot.amb.nov.2013$Photo~tz.d.spots$ALEAF))
# # 
# # get.bic.func(tz.d.spots$ALEAF,spot.amb$Photo,n.fd=5)
# 
# y<-seq(0.001,0.002,0.0001)
# x<-seq(0.001,0.002,0.0001)
# y<-rnorm(100)
# x<-rnorm(100)
# m<-lm(y ~ x) 
# res<-m$residuals
# n<-nrow(m$model)    
# w<-rep(1,n) #not applicable
# 
# ll <- 0.5 * (sum(log(w)) - n * (log(2 * pi) + 1 - log(n) + log(sum(w * res^2))))
# ll-logLik(m)==0 #TRUE
# k.original<-m$rank
# df.ll<-k.original+1 
# bic<- -2 * ll + log(n) * df.ll
# 
# df.ll <- 0.5 * (sum(log(w)) - obs.n * (log(2 * pi) + 1 - log(obs.n) + log(sum(w * (resi.vec)^2))))
# # xlogπ+(n−x)log(1−π)
# # -2 * ll + log(n) * df.ll
# bic <- -2*df.ll + log(obs.n) * (n.fd + 1)
# 
# lm1 <- lm(Fertility ~ Agriculture , data = swiss)
# AIC(lm1)
# stopifnot(all.equal(AIC(lm1),
#                     AIC(logLik(lm1))))
# BIC(lm1)
# lm1$residuals
