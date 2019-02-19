par_obj <- function(){
  
  rAB=0.05;
  rBC=0.15;
  rCD=0.3;
  C1=2;
  C2=0;
  C3=1;
  C4=1;
  rAC=rAB+rBC-2*C1*rAB*rBC
  rBD=rBC+rCD-2*C2*rBC*rCD
  rAD=rAB+rBC+rCD-2*C1*rAB*rBC-2*C2*rBC*rCD-2*C3*rAB*rCD+4*C4*rAB*rBC*rCD
  list(rAB=rAB,rBC=rBC,rCD=rCD,rAC=rAC,rBD=rBD,rAD=rAD,C1=C1,C2=C2,C3=C3,C4=C4)
}


simu_BC <- function(par_obj,n=400){
  
  rAB=par_obj$rAB;rBC=par_obj$rBC;rCD=par_obj$rCD;
  C1=par_obj$C1;C2=par_obj$C2;C3=par_obj$C3;C4=par_obj$C4;
  
  g111 <- C4*rAB*rBC*rCD
  g110 <- C1*rAB*rBC-C4*rAB*rBC*rCD
  g101 <- C3*rAB*rCD-C4*rAB*rBC*rCD
  g100 <- rAB -C1*rAB*rBC-C3*rAB*rCD-C4*rAB*rBC*rCD
  g011 <- C2*rBC*rCD-C4*rAB*rBC*rCD
  g010 <- rBC-C1*rAB*rBC-C2*rBC*rCD+C4*rAB*rBC*rCD
  g001 <- rCD-C2*rBC*rCD-C3*rAB*rCD+C4*rAB*rBC*rCD
  g000 <- 1-g111-g110-g101-g100-g011-g010-g001
  
  prob <- c(g000,g001,g011,g010,g110,g111,g101,g100,
            g100,g101,g111,g110,g010,g011,g001,g000)/2
  prob.1 <- cumsum(prob)
  Gam.tmp1 <- c("AABB","AABb","AaBB","AaBb")
  Gam.tmp2 <- c("CCDD","CCDd","CcDD","CcDd")
  Gam.tmp <- c()
  for(i in 1:length(Gam.tmp1)){
    Gam.tmp <- c(Gam.tmp,paste(Gam.tmp1[i],Gam.tmp2,sep=""))
  }
  gam <- matrix(NA,n,1)
  for(i in 1:n){
    xx <- runif(1)
    xx.ord <- min(which( xx < prob.1) )
    gam[i] <- c(Gam.tmp)[xx.ord]
  }
  gam
}

simu_F2 <- function(par_obj,n=200){
  
  rAB=par_obj$rAB;rBC=par_obj$rBC;rCD=par_obj$rCD;
  C1=par_obj$C1;C2=par_obj$C2;C3=par_obj$C3;C4=par_obj$C4;
  
  g111 <- C4*rAB*rBC*rCD
  g110 <- C1*rAB*rBC-C4*rAB*rBC*rCD
  g101 <- C3*rAB*rCD-C4*rAB*rBC*rCD
  g100 <- rAB -C1*rAB*rBC-C3*rAB*rCD-C4*rAB*rBC*rCD
  g011 <- C2*rBC*rCD-C4*rAB*rBC*rCD
  g010 <- rBC-C1*rAB*rBC-C2*rBC*rCD+C4*rAB*rBC*rCD
  g001 <- rCD-C2*rBC*rCD-C3*rAB*rCD+C4*rAB*rBC*rCD
  g000 <- 1-g111-g110-g101-g100-g011-g010-g001
  
  prob <- c(g000^2,2*g000*g001,g001^2,2*g000*g011,2*(g000*g010+g001*g011),2*g001*g010,g011^2,2*g011*g010,g010^2,
            2*g000*g110,2*(g000*g111+g001*g110),2*g001*g111,2*(g000*g101+g011*g110),2*(g000*g100+g001*g101+g010*g110+g011*g111),2*(g001*g100+g010*g111),2*g011*g101,2*(g010*g101+g011*g100),2*g010*g100,
            g110^2,2*g110*g111,g111^2,2*g110*g101,2*(g110*g100+g111*g101),2*g111*g100,g101^2,2*g101*g100,g100^2,
            2*g000*g100,2*(g000*g101+g001*g100),2*g001*g101,2*(g000*g111+g011*g100),2*(g000*g110+g001*g111+g010*g100+g011*g101),2*(g001*g110+g010*g101),2*g011*g111,2*(g010*g111+g011*g110),2*g010*g110,
            2*(g000*g010+g110*g100),2*(g000*g011+g001*g010+g110*g101+g111*g100),2*(g001*g011+g111*g101),2*(g000*g001+g011*g010+g110*g111+g101*g100),2*(g000^2+g001^2+g010^2+g011^2+g100^2+g101^2+g110^2+g111^2),2*(g000*g001+g011*g010+g110*g111+g101*g100),2*(g001*g011+g111*g101),2*(g000*g011+g001*g010+g110*g101+g111*g100),2*(g000*g010+g110*g100),
            2*g010*g110,2*(g010*g111+g011*g110),2*g011*g111,2*(g001*g110+g010*g101),2*(g000*g110+g001*g111+g010*g100+g011*g101),2*(g000*g111+g011*g100),2*g001*g101,2*(g000*g101+g001*g100),2*g000*g100,
            g100^2,2*g101*g100,g101^2,2*g111*g100,2*(g110*g100+g111*g101),2*g110*g101,g111^2,2*g110*g111,g110^2,
            2*g010*g100,2*(g010*g101+g011*g100),2*g011*g101,2*(g001*g100+g010*g111),2*(g000*g100+g001*g101+g010*g110+g011*g111),2*(g000*g101+g011*g110),2*g001*g111,2*(g000*g111+g001*g110),2*g000*g110,
            g010^2,2*g011*g010,g011^2,2*g001*g010,2*(g000*g010+g001*g011),2*g000*g011,g001^2,2*g000*g001,g000^2)/4
  prob.1 <- cumsum(prob)
  Gam.tmp1 <- c("AABB","AABb","AAbb","AaBB","AaBb","Aabb","aaBB","aaBb","aabb")
  Gam.tmp2 <- c("CCDD","CCDd","CCdd","CcDD","CcDd","Ccdd","ccDD","ccDd","ccdd")
  Gam.tmp <- c()
  for(i in 1:length(Gam.tmp1)){
    Gam.tmp <- c(Gam.tmp,paste(Gam.tmp1[i],Gam.tmp2,sep=""))
  }
  gam <- matrix(NA,n,1)
  for(i in 1:n){
    xx <- runif(1)
    xx.ord <- min(which( xx < prob.1) )
    gam[i] <- c(Gam.tmp)[xx.ord]
  }
  gam
}