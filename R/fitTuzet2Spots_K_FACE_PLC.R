rm(list=ls())
cat("\014")
s.time <- Sys.time()
library(DEoptim)
source("R/photosynTest.R")
source("r/readSpotAmb.R")

spot.amb$Vcmax.aci
spot.amb.nov.2013 <- (spot.amb)
spot.amb.nov.2013$LWP.pd <- spot.amb.nov.2013$Predawn.wp

tuzets.cost <- function(pars,dat){
  
  #- pull out the parameters from the pars vector
  psif <- pars[1]
  sf <- pars[2]
  Kmax <- pars[3]
  g1 <- pars[4]
  A <- euc.plc[[1]]
  B <- -euc.plc[[2]]
  


  Kactual <- Kmax
  
  out <- PhotosynTuzet(g1=g1,
                       sf=sf, 
                       psif=psif,
                       kl=Kactual,
                       A = A,
                       B = B,
                       Jmax = dat$Jmax.aci,
                       Vcmax = dat$Vcmax.aci,
                       Ca=dat$CO2S,psis=dat$LWP.pd,
                       VPD=dat$VpdL,
                       k.test = TRUE,
                       PPFD=1800,Tleaf=dat$Tleaf)
  
  #
  max.Gs <- sd(dat$Cond)
  max.A <- sd(dat$Photo)
  
  resid.gs <- ((out$GS - dat$Cond)/max.Gs)^2
  resid.A <- ((out$ALEAF - dat$Photo)/max.A)^2
  
  resid.gs[is.na(resid.gs)] <- 10
  resid.A[is.na(resid.A)] <- 10
 
  resid.sum <- sum(resid.gs,resid.A)
  
  return(resid.sum)
}

# pars are psiv, SF, K,g1 
lower <- c(-5,   1,1,2) 
upper <- c(-0.1,25,5,15)
NPmax <- 100
maxiter <- 100

#- set seed for repeatability
set.seed(1234)

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#- model with no change in photosynthetic capacity
spots.fit <- DEoptim(fn=tuzets.cost,lower=lower,upper=upper,
                     dat=spot.amb.nov.2013,
                     DEoptim.control(NP = NPmax,itermax=maxiter,trace=T,parallelType = 1,
                                           packages=list("plantecophys"),
                                           parVar = list("PhotosynTuzet","PhotosynTuzet_f","KPfnc",
                                                         "fsig_tuzet","psil_e","euc.plc")))

spots.fit.best.spots.k <- unname(spots.fit$optim$bestmem)
saveRDS(spots.fit,"Tuzet2spots_g1_K_FACE_PLC.rds")
e.time <- Sys.time() - s.time
print(e.time)




