rm(list=ls())
cat("\014")
s.time <- Sys.time()
library(DEoptim)
library(plantecophys)
source("R/photosynTest.R")

tuzets.cost <- function(pars,dat){
  #- pull out the parameters from the pars vector
  psif <- pars[1]
  sf <- pars[2]
  Kmax <- pars[3]
  psif.v <- pars[4]
  sf.v <- pars [5]
  g1 <- pars [6]
  
  A <- euc.plc[[1]]
  B <- -euc.plc[[2]]
 
  Kactual <- Kmax
 
  out <- PhotosynTuzet(sf=sf, psif=psif, 
                       psif.v = psif.v,
                       sf.v = sf.v,
                       g1=g1,kl=Kactual,
                       Ca=dat$CO2S,
                       psis=dat$LWP.pd,
                       Jmax = dat$Jmax.aci,
                       Vcmax = dat$Vcmax.aci,
                       VPD=dat$VpdL,PPFD=1800,Tleaf=dat$Tleaf,
                       v.psi.test = TRUE,
                       k.test = TRUE,
                       A = A,
                       B = B)
  
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

# PhotosynTuzet <- PhotosynTuzet
# pars are psiv, SF, K, b, and c
lower <- c(-4,  1, 1,  -4,0.5,1) # changed minimum Kmax from 2 to 5
upper <- c(-0.1,25,4,-0.1,25, 15)
NPmax <- 100
maxiter <- 50

#- set seed for repeatability
set.seed(1234)

#------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------
#- model with no change in photosynthetic capacity
spots.fit <- DEoptim(fn=tuzets.cost,lower=lower,upper=upper,
                     dat=spot.amb.nov.2013,
                     DEoptim.control(NP = NPmax,itermax=maxiter,trace=T,parallelType = 1,
                                           packages=list("plantecophys"),
                     parVar = list("PhotosynTuzet","PhotosynTuzet_f","KPfnc","fsig_tuzet","psil_e","euc.plc")))

spots.fit.best <- unname(spots.fit$optim$bestmem)
saveRDS(spots.fit,"Tuzet2spots_psi.rds")

e.time <- Sys.time() - s.time
print(e.time)

