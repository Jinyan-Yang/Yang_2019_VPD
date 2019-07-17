library(plantecophys)
source("r/photosynTest.R")
source("r/readSpotAmb.R")
source('r/load.r')

# 
spot.amb.nov.2013 <- (spots.amb)#[spot.amb$Campaign == "2013Nov",]
# spot.amb.nov.2013$LWP.pd <- mean(gs.wp.sap.lai.df$WP[gs.wp.sap.lai.df$Measurement == "Predawn"])
spot.amb.nov.2013$LWP.pd <- spot.amb.nov.2013$Predawn.wp
euc.plc <- list(SX = 32.389,P50 = 4.055)
# Tuzet PLC####
spots.fit <- unname(readRDS("Tuzet2spots_g1_K_FACE_PLC.rds")$optim$bestmem)
tz.spots <- PhotosynTuzet(Ca =spot.amb.nov.2013$CO2S, 
                          psis = spot.amb.nov.2013$LWP.pd,
                          VPD = spot.amb.nov.2013$VpdL, 
                          Tleaf = spot.amb.nov.2013$Tleaf, 
                          Jmax = spot.amb.nov.2013$Jmax.aci, Vcmax = spot.amb.nov.2013$Vcmax.aci,
                          
                          kl = spots.fit[3],
                          sf = spots.fit[2],
                          psif = spots.fit[1],
                          g1 = spots.fit[4], 
                          
                          PPFD = 1800,
                          k.test = TRUE,
                          A = euc.plc[[1]],
                          B = -euc.plc[[2]]
)

# optBB####
OptBB.de.best <- readRDS("Medlyn.rds")
spot.amb.nov.2013$g1.m <- OptBB.de.best[[1]]*exp(OptBB.de.best[[2]]*spot.amb.nov.2013$LWP.pd)
optbb.spots <- Photosyn(g1 = spot.amb.nov.2013$g1.m,
                        gsmodel = "BBOpti",
                        
                        Ca =spot.amb.nov.2013$CO2R, 
                        VPD = spot.amb.nov.2013$VpdL, 
                        Tleaf = spot.amb.nov.2013$Tleaf, 
                        
                        Jmax = spot.amb.nov.2013$Jmax.aci, Vcmax = spot.amb.nov.2013$Vcmax.aci,
                        alpha = 0.3, theta = 0.4756,
                        PPFD = 1800,
                        
                        EaV = 74189.7435218429, 
                        EdVC = 2e+05, 
                        delsC = 641.989, 
                        
                        EaJ = 39513,
                        EdVJ = 2e+05, 
                        delsJ = 640.2658
)
# Tuzet V-PSI####
spots.psi.fit <- unname(readRDS("Tuzet2spots_psi.rds")$optim$bestmem)

tz.psi.spots <- PhotosynTuzet(Ca =spot.amb.nov.2013$CO2R, 
                              psis = spot.amb.nov.2013$LWP.pd,
                              VPD = spot.amb.nov.2013$VpdL, 
                              Tleaf = spot.amb.nov.2013$Tleaf, 
                              Jmax = spot.amb.nov.2013$Jmax.aci, 
                              Vcmax = spot.amb.nov.2013$Vcmax.aci,
                              
                              kl = spots.psi.fit[3],
                              sf = spots.psi.fit[2],
                              psif = spots.psi.fit[1],
                              sf.v = spots.psi.fit[5],
                              psif.v = spots.psi.fit[4],
                              g1 = spots.psi.fit[6],
                              v.psi.test=TRUE,
                              PPFD = 1800,
                              
                              k.test = TRUE,
                              A=euc.plc[[1]],
                              B = -euc.plc[[2]])


# optBB V-D####
OptBB.VD.de.best <- readRDS("Medlyn_V_D.rds")
spot.amb.nov.2013$g1.m <- OptBB.VD.de.best[[1]]*exp(OptBB.VD.de.best[[2]]*spot.amb.nov.2013$LWP.pd)
opt.d <- Photosyn(g1 = spot.amb.nov.2013$g1.m, 
                  
                  Ca = spot.amb.nov.2013$CO2R, 
                  VPD = spot.amb.nov.2013$VpdL, 
                  Tleaf = spot.amb.nov.2013$Tleaf, 
                  
                  Jmax = spot.amb.nov.2013$Jmax.aci,
                  Vcmax = spot.amb.nov.2013$Vcmax.aci * (1.0 - OptBB.VD.de.best[[3]] * spot.amb.nov.2013$VpdL),
                  
                  gsmodel = "BBOpti",
                  PPFD = 1800,
                  
                  alpha = 0.3, theta = 0.4756,
                  EaV = 74189.7435218429, 
                  EdVC = 2e+05, 
                  delsC = 641.989, 
                  
                  EaJ = 39513,
                  EdVJ = 2e+05, 
                  delsJ = 640.2658)


# Leuning####
# leuning.best.spots.k <- unname(readRDS("Leuning2spots.rds")$optim$bestmem)
leuning.DE <- readRDS("Leuning.rds")
spot.amb.nov.2013$g1.l <-  leuning.DE[[1]]*exp(leuning.DE[[2]]*spot.amb.nov.2013$LWP.pd)
leuning.df <- Photosyn(g1 = spot.amb.nov.2013$g1.l,
                       D0 = leuning.DE[[3]],
                       gsmodel = "BBLeuning",
                       Ca =spot.amb.nov.2013$CO2S,
                       VPD = spot.amb.nov.2013$VpdL,
                       Tleaf = spot.amb.nov.2013$Tleaf,
                       Jmax = spot.amb.nov.2013$Jmax.aci,
                       Vcmax = spot.amb.nov.2013$Vcmax.aci,
                       
                       alpha = 0.3, theta = 0.4756,
                       PPFD = 1800,
                       EaV = 74189.7435218429, 
                       EdVC = 2e+05, 
                       delsC = 641.989, 
                       
                       
                       EaJ = 39513,
                       EdVJ = 2e+05, 
                       delsJ = 640.2658)

leuning.cable.df <- Photosyn(g1 = 9, #from De kauwe 2015 
                             D0 = 1.5,
                             gsmodel = "BBLeuning",
                             Ca =spot.amb.nov.2013$CO2S,
                             VPD = spot.amb.nov.2013$VpdL,
                             Tleaf = spot.amb.nov.2013$Tleaf,
                             Jmax = spot.amb.nov.2013$Jmax.aci,
                             Vcmax = spot.amb.nov.2013$Vcmax.aci,
                             
                             alpha = 0.3, theta = 0.4756,
                             PPFD = 1800,
                             EaV = 74189.7435218429, 
                             EdVC = 2e+05, 
                             delsC = 641.989, 
                             
                             
                             EaJ = 39513,
                             EdVJ = 2e+05, 
                             delsJ = 640.2658)

