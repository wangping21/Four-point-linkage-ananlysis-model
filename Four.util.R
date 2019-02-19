
D <- function(H.M,lpindex){
  
  hap1 <- H.M$hap1
  hap2 <- H.M$hap2
  genotype <- H.M$genotype
  lp1 <- H.M$lp1
  lp2 <- H.M$lp2
  ncom <- dim(H.M$genotype)
  DALL <- list()
  for(ii in c(lpindex)){
    G000.tmp <- matrix(0,nrow=ncom[1],ncol=ncom[2]);G001.tmp <- matrix(0,nrow=ncom[1],ncol=ncom[2]);
    G010.tmp <- matrix(0,nrow=ncom[1],ncol=ncom[2]);G011.tmp <- matrix(0,nrow=ncom[1],ncol=ncom[2]);
    G100.tmp <- matrix(0,nrow=ncom[1],ncol=ncom[2]);G101.tmp <- matrix(0,nrow=ncom[1],ncol=ncom[2]);
    G110.tmp <- matrix(0,nrow=ncom[1],ncol=ncom[2]);G111.tmp <- matrix(0,nrow=ncom[1],ncol=ncom[2]);
    for(jj in 1:ncom[1]){
      for(kk in 1:ncom[2]){
        geno <- strsplit(as.character(genotype[jj,kk]),"")[[1]]
        hap1.tmp <- as.numeric(paste(geno[1],geno[3],geno[5],geno[7],sep=""))
        hap2.tmp <- as.numeric(paste(geno[2],geno[4],geno[6],geno[8],sep=""))
        hap.tmp <- c(lp1[ii,which(hap1.tmp==hap1)],lp2[ii,which(hap2.tmp==hap2)])
        G000.tmp[jj,kk] <- length(which(hap.tmp=="g000"))
        G001.tmp[jj,kk] <- length(which(hap.tmp=="g001"))
        G010.tmp[jj,kk] <- length(which(hap.tmp=="g010"))
        G011.tmp[jj,kk] <- length(which(hap.tmp=="g011"))
        G100.tmp[jj,kk] <- length(which(hap.tmp=="g100"))
        G101.tmp[jj,kk] <- length(which(hap.tmp=="g101"))
        G110.tmp[jj,kk] <- length(which(hap.tmp=="g110"))
        G111.tmp[jj,kk] <- length(which(hap.tmp=="g111"))
      }
    }
    DALL <- cbind(G000.tmp,G001.tmp,G010.tmp,G011.tmp,G100.tmp,G101.tmp,G110.tmp,G111.tmp)
  }
  
  DALL
}


I.M <- function(segr.type){
  
  switch(EXPR = segr.type, A.1 = , A.2 = , A.3 = , A.4 = {
    I <- m <- diag(4)
  }, B1.5 = {
    I <- matrix(c(1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1), 4, 
                3, byrow = TRUE)
    m <- diag(3)
  }, B2.6 = {
    I <- matrix(c(1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1), 4, 
                3, byrow = TRUE)
    m <- diag(3)
  }, B3.7 = {
    I <- matrix(c(1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1), 4, 
                3, byrow = TRUE)
    m <- diag(3)
  }, C.8 = {
    I <- matrix(c(1, 0, 1, 0, 1, 0, 0, 1), 4, 2, byrow = TRUE)
    m <- diag(2)
  }, D1.9 = , D1.10 = , D1.11 = , D1.12 = , D1.13 = {
    I <- matrix(c(1, 0, 1, 0, 0, 1, 0, 1), 4, 2, byrow = TRUE)
    m <- diag(2)
  }, D2.14 = , D2.15 = , D2.16 = , D2.17 = , D2.18 = {
    I <- matrix(c(1, 0, 0, 1, 1, 0, 0, 1), 4, 2, byrow = TRUE)
    m <- diag(2)
  })
  
  list(I=I,m=m)
  
}

















H <- function(init.g = rep(1/8,8),lpindex=1){
  
  g000 <- init.g[1];g001 <- init.g[2];g010 <- init.g[3];g011 <- init.g[4];
  g100 <- init.g[5];g101 <- init.g[6];g110 <- init.g[7];g111 <- init.g[8];
  h1 <- hap()
  hap1lp <- matrix(c("g000","g001","g011","g010","g110","g111","g101","g100",
                     "g001","g000","g010","g011","g111","g110","g100","g101",
                     "g011","g010","g000","g001","g101","g100","g110","g111",
                     "g010","g011","g001","g000","g100","g101","g111","g110",
                     "g110","g111","g101","g100","g000","g001","g011","g010",
                     "g111","g110","g100","g101","g001","g000","g010","g011",
                     "g101","g100","g110","g111","g011","g010","g000","g001",
                     "g100","g101","g111","g110","g010","g011","g001","g000"),ncol=8)
  lp1 <- c()
  lp2 <- c()
  for(i in 1:dim(hap1lp)[1]){
    for(j in 1:dim(hap1lp)[1]){
      lp1 <- rbind(lp1,c(hap1lp[i,],rev(hap1lp[i,])));
      lp2 <- rbind(lp2,c(hap1lp[j,],rev(hap1lp[j,])));
    }
  }
  
  ncom <- dim(h1$genotyp)
  #HALL <- list()
  for(ii in c(lpindex)){
    h.tmp <- matrix(NA,nrow=ncom[1],ncol=ncom[2])
    for(jj in 1:ncom[1]){
      for(kk in 1:ncom[2]){
        geno <- strsplit(as.character(h1$genotyp[jj,kk]),"")[[1]]
        hap1.tmp <- as.numeric(paste(geno[1],geno[3],geno[5],geno[7],sep=""))
        hap2.tmp <- as.numeric(paste(geno[2],geno[4],geno[6],geno[8],sep=""))
        h.tmp[jj,kk] <- eval(parse(text=lp1[ii,which(hap1.tmp==h1$hap1)]))*
          eval(parse(text=lp2[ii,which(hap2.tmp==h1$hap2)]))
      }
    }
  HALL <- h.tmp
  }
  list(HALL=HALL,lp1=lp1,lp2=lp2,hap1=h1$hap1,hap2=h1$hap2,genotype=h1$genotyp)
}

log_sub <- function(x,y){
  if(y > x + 300) return(y)
  if(x > y + 300) return(x)
  if(x <y ){
    tt <- x
    x <- y
    y <- tt   
  }
  return(x + log(1 - exp(y-x)))
}




hap <- function(){
  
  index1 <- c()
  for (i in 1:2){
    for (ii in 1:2){
      for (iii in 1:2){
        for(iiii in 1:2){
          index1 <- c(index1,i*1000+ii*100+iii*10+iiii)
        }
      }
    }
  }
  
  index2 <- c()
  for (i in 3:4){
    for (ii in 3:4){
      for (iii in 3:4){
        for(iiii in 3:4){
          index2 <- c(index2,i*1000+ii*100+iii*10+iiii)
        }
      }
    }
  }
  
  index3 <- c()
  for(j in c(13,14,23,24)){
    for(jj in c(13,14,23,24)){
      for(jjj in c(13,14,23,24)){
        for(jjjj in c(13,14,23,24)){
          index3 <- c(index3,as.numeric(paste(j,jj,jjj,jjjj,sep="")))
        }
      }

    }
  }
  
  nindex3 <- matrix(index3,ncol=4,byrow=T)
  list(hap1=index1,hap2=index2,genotyp=nindex3)
}