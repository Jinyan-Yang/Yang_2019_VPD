psil_e <- function(ELEAF, kl, psis){
  
  psil <- psis - ELEAF / kl
  
  psil[!is.finite(psil)] <- psis
  
  psil
}

fsig_tuzet <- function(psil, sf=3.2, psif=-1.9){
  (1+exp(sf*psif))/(1+exp(sf*(psif-psil)))
}
# KPfnc(-3,23,-4)
KPfnc <- function(P, SX, PX,psi.cruit=-50,crit.psi.test = FALSE) {
  if(crit.psi.test == TRUE){
    if(P > psi.cruit){
      X <- 50
      V <- (X - 100) * log(1 - X/100)
      p <- (P/PX)^((PX * SX)/V)
      relk <- (1 - X/100)^p
     
    }else{
      0.1+0.0001*P
    }
  }else{
    X <- 50
    V <- (X - 100) * log(1 - X/100)
    p <- (P/PX)^((PX * SX)/V)
    relk <- (1 - X/100)^p
   
  }
   return(1-relk)
}

PhotosynTuzet_f <- function(g1=4,
                            Ca=400,
                            psis=0, 
                            kl=2, 
                            sf=3, 
                            psif=-2,
                            VPD=1,
                            Jmax = 145,
                            Vcmax = 90,
                            v.a = -0.1,
                            v.psi.2=FALSE,
                            v.psi.test=FALSE,
                            k.test =FALSE,
                            A=4,
                            B=-2, 
                            psi.cruit = -3,
                            crit.psi.test = FALSE,
                            v.crit.test =FALSE,
                            psi.v.crit =-50,
                            sf.v = 3, 
                            psif.v = -4,
                            ...){
  
  O <- function(psil, psis, sf, psif, g1, Ca,  
                kl = kl, 
                Jmax = Jmax,
                Vcmax = Vcmax,
                A,
                B,
                v.a,
                v.psi.test,
                v.psi.2,
                k.test,
                psi.cruit,
                crit.psi.test,
                sf.v, psif.v,
                v.crit.test ,
                psi.v.crit ,
                v.min,
                ...){
    # change vcmax and jmax as a function of psi
    Jmax.use = Jmax
    Vcmax.use = Vcmax
    
    if(v.psi.test==TRUE){
      Jmax.use = Jmax #* fsig_tuzet(psil, sf.v, psif.v)
      Vcmax.use = Vcmax* fsig_tuzet(psil, sf.v, psif.v)
    }
    if(v.psi.2==TRUE){
      Jmax.use = Jmax * (1-v.a*psil)
      Vcmax.use = Vcmax* (1-v.a*psil)
    }
    # change k as a function of psi
    if(k.test==TRUE){
      kl.use = kl * KPfnc(psil,A,B)
    }else{
      kl.use = kl
    }
    
    p <- Photosyn(g1=g1, Ca=Ca, gsmodel="BBdefine", BBmult=(g1/Ca)*fsig_tuzet(psil, sf, psif),
                  Jmax = Jmax.use,
                  Vcmax = Vcmax.use,
                  
                  alpha = 0.3, theta = 0.4756,
                  EaV = 74189.7435218429, 
                  EdVC = 2e+05, 
                  delsC = 641.989, 
                  
                  EaJ = 39513,
                  EdVJ = 2e+05, 
                  delsJ = 640.2658,
                  ...)
    psilout <- psil_e(p$ELEAF,kl.use, psis)
    
    psil - psilout  # objective function: psil in = psil out.
  }
  
  topt <- uniroot(O, c(-20,-0.05), psis=psis, kl=kl, sf=sf, psif=psif, g1=g1, Ca=Ca,VPD=VPD, 
                  Jmax = Jmax ,
                  Vcmax = Vcmax,
                  A=A,B=B,
                  v.a=v.a,
                  v.psi.test=v.psi.test,
                  v.psi.2=v.psi.2,
                  k.test = k.test,
                  psi.cruit=psi.cruit,crit.psi.test=crit.psi.test,
                  sf.v = sf.v, psif.v = sf.v,
                  v.crit.test =v.crit.test,
                  psi.v.crit =v.crit.test,
                  ...)
  psil.root <- topt$root
  
    Jmax.use = Jmax
    Vcmax.use = Vcmax

  if(v.psi.test==TRUE){
    Jmax.use = Jmax #* fsig_tuzet(topt$root, sf.v, psif.v)
    Vcmax.use = Vcmax* fsig_tuzet(topt$root,sf.v, psif.v)
  }

  if(v.psi.2==TRUE){
    Jmax.use = Jmax #* (1-v.a*topt$root)
    Vcmax.use = Vcmax* (1-v.a*topt$root)
  }
    if(v.crit.test==TRUE){
    # if(psil.root<=psi.v.crit){
      ov <- function(Vcmax,
                     psil , 
                     psis, sf, psif, g1, Ca,  
                    kl , 
                    Jmax,
                    ...){
        # change vcmax and jmax as a function of psi
        Jmax.use = Jmax
        Vcmax.use = Vcmax
        p <- Photosyn(g1=g1, Ca=Ca, gsmodel="BBdefine", 
                      BBmult=(g1/Ca)*fsig_tuzet(psil, sf, psif),
                      Jmax = Jmax.use,
                      Vcmax = Vcmax.use,
                      
                      alpha = 0.3, theta = 0.4756,
                      EaV = 74189.7435218429, 
                      EdVC = 2e+05, 
                      delsC = 641.989, 
                      
                      EaJ = 39513,
                      EdVJ = 2e+05, 
                      delsJ = 640.2658,
                      ...)
        psilout <- psil_e(p$ELEAF,kl,psis)
        
        psil - psilout  # objective function: psil in = psil out.
      # }
      
      topt.v <- uniroot(ov, c(1,120),psil = psi.v.crit,psis=psis, kl=kl, sf=sf, psif=psif, g1=g1, Ca=Ca,VPD=VPD, 
                      Jmax = Jmax ,
                      ...)
      vc.root <- topt.v$root
      
      psil.root <- psi.v.crit
      Vcmax.use <- vc.root
    }
  }    
    
  p <- Photosyn(g1=g1, Ca=Ca, VPD=VPD,gsmodel="BBdefine",
                BBmult=(g1/Ca)*fsig_tuzet(psil.root, sf, psif),
                Jmax = Jmax.use,
                Vcmax = Vcmax.use,
                
                alpha = 0.3, theta = 0.4756,
                EaV = 74189.7435218429, 
                EdVC = 2e+05, 
                delsC = 641.989, 
                
                EaJ = 39513,
                EdVJ = 2e+05, 
                delsJ = 640.2658,
                ...)
  p2 <- cbind(p, data.frame(PSIL=psil.root,
                            vcmax = Vcmax.use))
  
  return(p2)
}

PhotosynTuzet <- function(g1=8,
                          Ca=400,
                          psis=0, 
                          kl=2, 
                          sf=3, 
                          psif=-2,
                          vpd.test = FALSE,
                          Jmax = 145,
                          Vcmax = 90,
                          VPD=1,
                          v.a = -0.1,
                          v.psi.2=FALSE,
                          v.psi.test=FALSE,
                          k.test = FALSE,
                          A=0.8,
                          B=-4,
                          sf.v = 3, psif.v = -4,
                          psi.cruit=-10,
                          crit.psi.test=FALSE,
                          v.crit.test =FALSE,
                          psi.v.crit =-3.2,
                          cd=0.17,
                          ...){
  if(vpd.test == TRUE){
    Vcmax.use = Vcmax * (1.0 - cd * VPD)
    Jmax.use = Jmax * (1.0 - cd * VPD)
  }else{
    Jmax.use = Jmax
    Vcmax.use = Vcmax 
  }
  m <- mapply(PhotosynTuzet_f, g1=g1, Ca=Ca,VPD =VPD,
              psis=psis, kl=kl, sf=sf, psif=psif,
              Jmax = Jmax.use,
              Vcmax = Vcmax.use,
              v.psi.test=v.psi.test,
              v.a=v.a,
              v.psi.2=v.psi.2,
              k.test = k.test,
              A=A,
              B=B,
              sf.v = sf.v, psif.v = sf.v,
              psi.cruit=psi.cruit,crit.psi.test=crit.psi.test,
              v.crit.test =v.crit.test,
              psi.v.crit = psi.v.crit,
              ..., SIMPLIFY=FALSE)
  do.call(rbind, m)
}

# PhotosynTuzet(VPD=10,
#               crit.psi.test = T,
#               psi.cruit=-3,
#               k.test = T)
#   
# PhotosynTuzet(VPD=5,
#               g1=15,
#               psis=-0.1, 
#               kl=2, 
#               sf=3, 
#               psif=-2,
#               Vcmax = 90,
#               # k.test = T,
#               # A=19.982634 ,
#               # B=-2.334635,
#               v.crit.test = T,psi.v.crit = -2)
