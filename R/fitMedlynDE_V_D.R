# Start of main code###########
s.time <- Sys.time()
library(DEoptim)
library(plantecophys)

spot.amb.nov.2013 <- spot.amb
spot.amb.nov.2013$LWP.pd <- spot.amb.nov.2013$Predawn.wp

# Function to be passed to DE optim
optBB.cost <- function(pars,dat,Vcmax=80,Jmax=120){
  #- pull out the parameters from the pars vector
  g1.m <- pars[1]
  factor.b <- pars[2]
  factor.d <- pars[3]
  
  g1.vec <- g1.m * exp(dat$LWP.pd * factor.b)
  out <- Photosyn(g1 = g1.vec,
                  Vcmax = dat$Vcmax.aci * (1 - factor.d * dat$VpdL),
                  Jmax =dat$Jmax.aci * (1 - factor.d * dat$VpdL),
                  
                  Ca=dat$CO2S,
                  VPD=dat$VpdL,
                  PPFD=1800,
                  Tleaf=dat$Tleaf,

                  alpha = 0.3, theta = 0.4756,
                  
                  EaV = 74189.7435218429, 
                  EdVC = 2e+05, 
                  delsC = 641.989, 
                  
                  EaJ = 39513,
                  EdVJ = 2e+05, 
                  delsJ = 640.2658)
  
  #- Belinda suggested to weight by the standard deviation, rather than the maximum. It doesn't really matter.
  max.Gs <- sd(dat$Cond)
  max.A <- sd(dat$Photo)
  
  resid.gs <- ((out$GS - dat$Cond)/max.Gs)^2
  resid.A <- ((out$ALEAF - dat$Photo)/max.A)^2

  resid.gs[is.na(resid.gs)] <- 10
  resid.A[is.na(resid.A)] <- 10
 
  resid.sum <- sum(resid.A, resid.gs) 
  return(resid.sum)
}

# setting control parameters and limits to values
lower <- c(1,-1,0.1) 
upper <- c(20,5,0.5)
NPmax <- 100
maxiter <- 100

#- set seed for repeatability
set.seed(1234)

#- Call to DEoptim
OptBB.de.fit.vd <- DEoptim(fn=optBB.cost,lower=lower,upper=upper,
                        dat=spot.amb.nov.2013,
                        DEoptim.control(VTR = 1,
                                              NP = NPmax,itermax=maxiter,trace=T,parallelType = 1,
                                              packages=list("plantecophys","rootSolve")))

OptBB.de.fit.k.vd <- unname(OptBB.de.fit.vd$optim$bestmem)

saveRDS(OptBB.de.fit.k.vd,"Medlyn_V_D.rds")

e.time <- Sys.time() - s.time
print(e.time)

