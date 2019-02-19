four.est <- function(dat,Mindex=c(1,2,3,4),lp=c(1)){
  
  
  mrk1 <- dat$geno[,Mindex[1]]
  mrk2 <- dat$geno[,Mindex[2]]
  mrk3 <- dat$geno[,Mindex[3]]
  mrk4 <- dat$geno[,Mindex[4]]
  
  I1.tmp <- I.M(segr.type = dat$segr.type[Mindex[1]])
  I1 <- I1.tmp$I
  m1 <- I1.tmp$m
  
  I2.tmp <- I.M(segr.type = dat$segr.type[Mindex[2]])
  I2 <- I2.tmp$I
  m2 <- I2.tmp$m
  
  I3.tmp <- I.M(segr.type = dat$segr.type[Mindex[3]])
  I3 <- I3.tmp$I
  m3 <- I3.tmp$m
  
  I4.tmp <- I.M(segr.type = dat$segr.type[Mindex[4]])
  I4 <- I4.tmp$I
  m4 <- I4.tmp$m
  
  mrk1 <- dat$geno[,Mindex[1]]
  mrk2 <- dat$geno[,Mindex[2]]
  mrk3 <- dat$geno[,Mindex[3]]
  mrk4 <- dat$geno[,Mindex[4]]
  
  
  n <- numeric(nrow(m1)*nrow(m2)*nrow(m3)*nrow(m4))
  
  k <- 1
  for (p4 in 1:nrow(m4)) for (p1 in 1:nrow(m1)) for (p2 in 1:nrow(m2)) for (p3 in 1:nrow(m3)) {
    n[k] <- length(which(mrk1 == p1 & mrk2 == p2 & mrk3 == p3& mrk4 == p4))
    k <- k + 1
  }
  
  two_n <- 2*sum(n)
  
  rfall <- c()
  for(j in c(lp)){
    rffour <- rf.four(I1,I2,I3,I4,n,two_n,lpindex=j)
    rfall <- rbind(rfall,unlist(rffour))
  }
  H11 <- H(init.g=rep(1/8,8),lpindex=1)
  v <- as.vector(kronecker(kronecker(t(I1),t(I2)),t(I3))%*%H11$HALL%*%I4)
  log_like_null <- sum(n*log(v))
  LOD <- (rfall[,11]-log_like_null)/log(10)
  
  result <- cbind(rfall,LOD=LOD)
  rownames(result) <- NULL
  return(result)
}



rf.four <- function(I1,I2,I3,I4,n,two_n,lpindex=1){
  
  ninit.g <- rep(1/8,8)
  HALL.tmp <- H(init.g=ninit.g,lpindex)
  HALL <- HALL.tmp$HALL
  DALL <- D(H.M=HALL.tmp,lpindex)
  k1 <- kronecker(t(I1),t(I2))
  k2 <- kronecker(k1,t(I3))
  k3 <- k2 %*%HALL%*%I4
  
  nk <- as.vector(k3)
  
  log_Lnew <- sum(n*log(nk),na.rm=T)
  
  MAXITER <- 1000
  THRESH <- log(1E-4)
  diff <- 1
  iter <- 1
  
  while(diff>THRESH && iter<MAXITER) {
    iter <- iter + 1
    log_Lold <- log_Lnew
    
    mid1 <- DALL[,1:4]*HALL
    g000.temp <- as.vector(k2%*%mid1%*%I4)
    
    mid2 <- DALL[,5:8]*HALL
    g001.temp <- as.vector(k2%*%mid2%*%I4)
    
    mid3 <- DALL[,9:12]*HALL
    g010.temp <- as.vector(k2%*%mid3%*%I4)
    
    mid4 <- DALL[,13:16]*HALL
    g011.temp <- as.vector(k2%*%mid4%*%I4)
    
    mid5 <- DALL[,17:20]*HALL
    g100.temp <- as.vector(k2%*%mid5%*%I4)
    
    mid6 <- DALL[,21:24]*HALL
    g101.temp <- as.vector(k2%*%mid6%*%I4)
    
    mid7 <- DALL[,25:28]*HALL
    g110.temp <- as.vector(k2%*%mid7%*%I4)
    
    mid8 <- DALL[,29:32]*HALL
    g111.temp <- as.vector(k2%*%mid8%*%I4)
    
    
    mult <- n/nk
    g000 <- sum(mult*g000.temp)/two_n
    g001 <- sum(mult*g001.temp)/two_n 
    g010 <- sum(mult*g010.temp)/two_n
    g011 <- sum(mult*g011.temp)/two_n
    g100 <- sum(mult*g100.temp)/two_n
    g101 <- sum(mult*g101.temp)/two_n
    g110 <- sum(mult*g110.temp)/two_n
    g111 <- sum(mult*g111.temp)/two_n
    
    ninit.g <- c(g000,g001,g010,g011,g100,g101,g110,g111)
    HALL <- H(init.g=ninit.g,lpindex)$HALL
    nk <- as.vector(k2%*%HALL%*%I4)
    
    log_Lnew <- sum(n*log(nk),na.rm=T)
    diff <-  log_sub(log_Lnew,log_Lold) - log_Lnew #
  }
  rAB  <- g111 + g110 + g101 + g100
  rBC  <- g111 + g110 + g011 + g010
  rCD  <- g111 + g101 + g011 + g001
  
  C4 <- g111/(rAB*rBC*rCD)
  C1 <- (g111+g110)/(rAB*rBC)
  C2 <- (g111+g011)/(rBC*rCD)
  C3 <- (g111+g101)/(rAB*rCD)
  
  rAC=rAB+rBC-2*C1*rAB*rBC
  rBD=rBC+rCD-2*C2*rBC*rCD
  rAD=rAB+rBC+rCD-2*C1*rAB*rBC-2*C2*rBC*rCD-2*C3*rAB*rCD+4*C4*rAB*rBC*rCD
  list(rAB=rAB,rBC=rBC,rCD=rCD,rAC=rAC,rBD=rBD,rAD=rAD,
       C1=C1,C2=C2,C3=C3,C4=C4,log_like_assign=log_Lnew)
}





