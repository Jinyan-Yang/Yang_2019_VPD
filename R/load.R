source("r/readSpotAmb.R")
spot.amb.nov.2013 <- spot.amb
spot.amb.nov.2013$LWP.pd <- spot.amb.nov.2013$Predawn.wp
euc.plc <- list(SX = 32.389,P50 = 4.055)

if(!file.exists('Medlyn.rds')) source('r/fitMedlynDE.R')

if(!file.exists('Leuning.rds')) source('r/fitLeuningDE.R')

if(!file.exists('Medlyn_V_D.rds')) source('r/fitMedlynDE_V_D.R')

if(!file.exists('Tuzet2spots_psi.rds')) source('r/fitTzuet2spots_psi.R')

if(!file.exists('Tuzet2spots_g1_K_FACE_PLC.rds')) source('r/fitTuzet2Spots_K_FACE_PLC.R')