# load packages and functions#########
library(plantecophys)
library(rootSolve)

source("r/readSpotAmb.R")
# function to solve for ci given Anet
# get.gm.function(20)
get.gm.function <- function(a.in,...){
  obj.func <- function (gm,...){
    abs(Photosyn(gmeso = gm,...)$ALEAF - a.in)
    
  }
  
  # topt <- uniroot.all(obj.func,interval=c(0.0005,2),...)
  # out <- max(topt,na.rm=TRUE)
  topt <- optim(0.1,obj.func,lower=0.0005,upper=2, method = "L-BFGS-B",...)
  out <- max(topt$par,na.rm=TRUE)
  return(out)
}
# allows as above function to take vector inputs 
get.gm.m.func <- function(a.vec.in,...){
  result.m <- mapply(get.gm.function, 
                     a.vec.in,
                     ..., SIMPLIFY=FALSE)
  do.call(rbind, result.m)
}

# read data
dat <- spot.amb
# get gm for spots according to measured assimilation rate 
gm.vec <- get.gm.m.func(dat$Photo,
                        
                        Ca=dat$CO2R,
                        VPD=dat$VpdL,
                        PPFD=1800,
                        Tleaf=dat$Tleaf,
                        
                        # Jmax=145,Vcmax=90,
                        Jmax=dat$Jmax.aci,Vcmax=dat$Vcmax.aci,
                        alpha = 0.3, theta = 0.4756,
                        
                        EaV = 74189.7435218429, 
                        EdVC = 2e+05, 
                        delsC = 641.989, 
                        
                        EaJ = 39513,
                        EdVJ = 2e+05, 
                        delsJ = 640.2658)

gm.df=data.frame(gm=gm.vec,
                 vpd=dat$VpdL,
                 rh=dat$RH_S)
saveRDS(gm.df,'cache/gm.fit.df.rds')


par(mar=c(5,5,1,1))
plot(log(gm)~log(vpd),data = gm.df[gm.df$gm<2,],pch=16,col="grey",
     xlab=expression(D~(kPa)),ylab=expression(g[m]~(mol~m^-2~s^-1)),
     axes = FALSE,
     ylim=c(-6,0),
     xlim=c(0,2))
fit.gm <- (lm(log(gm)~log(vpd),data = gm.df[gm.df$gm<2,]))
axis(2,at = log(c(0.002,seq(0,0.5,0.1))),labels = c(0,seq(0,0.5,0.1)))
axis(1,at = log(1:6),labels = 1:6)
abline(fit.gm)
summary(fit.gm)$r.squared

