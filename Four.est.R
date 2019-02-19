Four.BC <- function(dat){
  n <- length(dat)
  Gam.tmp1 <- c("AABB","AABb","AaBB","AaBb")
  Gam.tmp2 <- c("CCDD","CCDd","CcDD","CcDd")
  Gam.tmp <- c()
  for(i in 1:length(Gam.tmp1)){
    Gam.tmp <- c(Gam.tmp,paste(Gam.tmp1[i],Gam.tmp2,sep=""))
  }
  n.index <- rep(0,length(Gam.tmp))
  for(i in 1:length(Gam.tmp)){
    n.index[i] <- sum((Gam.tmp[i]==dat))
  }
  g000 <- (n.index[1]+n.index[16])/n;
  g001 <- (n.index[2]+n.index[15])/n;
  g010 <- (n.index[4]+n.index[13])/n;
  g011 <- (n.index[3]+n.index[14])/n;
  g100 <- (n.index[8]+n.index[9])/n;
  g101 <- (n.index[7]+n.index[10])/n;
  g110 <- (n.index[5]+n.index[12])/n;
  g111 <- (n.index[6]+n.index[11])/n;
  
  rAB  <- g111+g110+g101+g100
  rBC  <- g111+g110+g011+g010
  rCD  <- g111+g101+g011+g001
  rAC  <- g101+g100+g011+g010
  rBD  <- g110+g010+g101+g001
  rAD  <- g111+g010+g100+g001
  C4 <- g111/(rAB*rBC*rCD)
  C1 <- (g111+g110)/(rAB*rBC)
  C2 <- (g111+g011)/(rBC*rCD)
  C3 <- (g111+g101)/(rAB*rCD)
  
  rAC=rAB+rBC-2*C1*rAB*rBC
  rBD=rBC+rCD-2*C2*rBC*rCD
  rAD=rAB+rBC+rCD-2*C1*rAB*rBC-2*C2*rBC*rCD-2*C3*rAB*rCD+4*C4*rAB*rBC*rCD
  
  list(rAB=rAB,rBC=rBC,rCD=rCD,rAC=rAC,rBD=rBD,rAD=rAD,C4=C4,C1=C1,C2=C2,C3=C3)
}  


Two.BC <- function(dat){
  
  datt <- dat.tran1(dat)
  
  r <- rf.2pts(datt,verbose = F)
  
  r.m <- r$analysis[,1,1]
  
  
  rAB <- r.m[1]
  rAC <- r.m[2]
  rBC <- r.m[3]
  rAD <- r.m[4]
  rBD <- r.m[5]
  rCD <- r.m[6]
  
  list(rAB=rAB,rBC=rBC,rCD=rCD,rAC=rAC,rBD=rBD,rAD=rAD)
}

Three.BC <- function(dat){
  
  datt <- dat.tran1(dat)
  
  r1 <- def.rf.3pts(datt,"M1","M2","M3")
  r1.m <- r1$analysis[,1:3]
  r2 <- def.rf.3pts(datt,"M1","M2","M4")
  r2.m <- r2$analysis[,1:3]
  r3 <- def.rf.3pts(datt,"M1","M3","M4")
  r3.m <- r3$analysis[,1:3]
  r4 <- def.rf.3pts(datt,"M2","M3","M4")
  r4.m <- r4$analysis[,1:3]
  
  rAB <- min(c(r1.m[,1],r2.m[,1]))
  rBC <- min(c(r1.m[,2],r4.m[,1]))
  rCD <- min(c(r3.m[,2],r4.m[,2]))
  rAC <- min(c(r1.m[,3],r3.m[,1]))
  rAD <- min(c(r2.m[,3],r3.m[,3]))
  rBD <- min(c(r2.m[,2],r4.m[,3]))
  
  list(rAB=rAB,rBC=rBC,rCD=rCD,rAC=rAC,rBD=rBD,rAD=rAD)
}

Three.F2 <- function(dat){
  
  datt <- dat.tran1(dat,type="F2")
  
  r1 <- def.rf.3pts(datt,"M1","M2","M3")
  r1.m <- r1$analysis[,1:3]
  r2 <- def.rf.3pts(datt,"M1","M2","M4")
  r2.m <- r2$analysis[,1:3]
  r3 <- def.rf.3pts(datt,"M1","M3","M4")
  r3.m <- r3$analysis[,1:3]
  r4 <- def.rf.3pts(datt,"M2","M3","M4")
  r4.m <- r4$analysis[,1:3]
  
  rAB <- min(c(r1.m[,1],r2.m[,1]))
  rBC <- min(c(r1.m[,2],r4.m[,1]))
  rCD <- min(c(r3.m[,2],r4.m[,2]))
  rAC <- min(c(r1.m[,3],r3.m[,1]))
  rAD <- min(c(r2.m[,3],r3.m[,3]))
  rBD <- min(c(r2.m[,2],r4.m[,3]))
  
  list(rAB=rAB,rBC=rBC,rCD=rCD,rAC=rAC,rBD=rBD,rAD=rAD)
}

Four.F2 <- function(dat,init.g){
  
  Gam.tmp1 <- c("AABB","AABb","AAbb","AaBB","AaBb","Aabb","aaBB","aaBb","aabb")
  Gam.tmp2 <- c("CCDD","CCDd","CCdd","CcDD","CcDd","Ccdd","ccDD","ccDd","ccdd")
  Gam.tmp <- c()
  for(i in 1:length(Gam.tmp1)){
    Gam.tmp <- c(Gam.tmp,paste(Gam.tmp1[i],Gam.tmp2,sep=""))
  }
  n <- length(dat)
  gam.n <- c(rep(0,length(Gam.tmp)))
  for(s in 1:length(Gam.tmp)){
    gam.n[s] <- sum(Gam.tmp[s]==dat)
  }
  n.matrix <- c()
  for(ii in 1:9){
    n.matrix <- rbind(n.matrix,matrix(rep(gam.n[(9*(ii-1)+1):(ii*9)],8),8,byrow=T))
  }
  for(l in 1:1){
    g000 <- init.g[1]
    g001 <- init.g[2]
    g010 <- init.g[3]
    g011 <- init.g[4]
    g100 <- init.g[5]
    g101 <- init.g[6]
    g110 <- init.g[7]
    g111 <- init.g[8]
    probi <- linkage.phase(lpindex=l)
    index.g <- c("g000","g001","g010","g011","g100","g101","g110","g111")
    
    
    old.g000 <- 0;old.g001 <- 0;old.g010 <- 0;old.g011 <- 0;
    old.g100 <- 0;old.g101 <- 0;old.g110 <- 0;old.g111 <- 0;
    loop_k <- 1
    while ((abs(g000 - old.g000) > 1e-05 || abs(g001 - old.g001) > 1e-05 ||
            abs(g010 - old.g010) > 1e-05||
            abs(g011 - old.g011) > 1e-05||
            abs(g100 - old.g100) > 1e-05||
            abs(g101 - old.g101) > 1e-05||
            abs(g110 - old.g110) > 1e-05||
            abs(g111 - old.g111) > 1e-05)&& loop_k < 1000){
      old.g000 <- g000;old.g001 <- g001;old.g010 <- g010;old.g011 <- g011;
      old.g100 <- g100;old.g101 <- g101;old.g110 <- g110;old.g111 <- g111;
      index.m <- matrix(0,72,9)
      for(i in 1:9){
        for(j in 1:9){
          AA <- probi[i,j]
          AA1 <- strsplit(AA,"-")[[1]]
          if(length(AA1)==2){
            for(k in 1:length(AA1)){
              index <- which(AA1[k]==index.g)
              index.m[(i-1)*8+index,j] <- index.m[(i-1)*8+index,j] + length(index)
            }
          }
          if(length(AA1)==4){
            phi <- eval(parse(text=AA1[1]))*eval(parse(text=AA1[2]))/
              (eval(parse(text=AA1[1]))*eval(parse(text=AA1[2]))+eval(parse(text=AA1[3]))*eval(parse(text=AA1[4])))
            vphi <- c(phi,phi,1-phi,1-phi)
            for(k in 1:length(AA1)){
              index <- which(AA1[k]==index.g)
              index.m[(i-1)*8+index,j] <- index.m[(i-1)*8+index,j] + vphi[k]
            }
          }
          if(length(AA1)==8){
            phi1 <- eval(parse(text=AA1[1]))*eval(parse(text=AA1[2]))/
              (eval(parse(text=AA1[1]))*eval(parse(text=AA1[2]))+eval(parse(text=AA1[3]))*eval(parse(text=AA1[4]))+eval(parse(text=AA1[5]))*eval(parse(text=AA1[6]))+eval(parse(text=AA1[7]))*eval(parse(text=AA1[8])))
            phi2 <- eval(parse(text=AA1[3]))*eval(parse(text=AA1[4]))/
              (eval(parse(text=AA1[1]))*eval(parse(text=AA1[2]))+eval(parse(text=AA1[3]))*eval(parse(text=AA1[4]))+eval(parse(text=AA1[5]))*eval(parse(text=AA1[6]))+eval(parse(text=AA1[7]))*eval(parse(text=AA1[8])))
            phi3 <- eval(parse(text=AA1[5]))*eval(parse(text=AA1[6]))/
              (eval(parse(text=AA1[1]))*eval(parse(text=AA1[2]))+eval(parse(text=AA1[3]))*eval(parse(text=AA1[4]))+eval(parse(text=AA1[5]))*eval(parse(text=AA1[6]))+eval(parse(text=AA1[7]))*eval(parse(text=AA1[8])))
            phi4 <- eval(parse(text=AA1[7]))*eval(parse(text=AA1[8]))/
              (eval(parse(text=AA1[1]))*eval(parse(text=AA1[2]))+eval(parse(text=AA1[3]))*eval(parse(text=AA1[4]))+eval(parse(text=AA1[5]))*eval(parse(text=AA1[6]))+eval(parse(text=AA1[7]))*eval(parse(text=AA1[8])))
            vphi <- c(phi1,phi1,phi2,phi2,phi3,phi3,phi4,phi4)
            for(k in 1:length(AA1)){
              index <- which(AA1[k]==index.g)
              index.m[(i-1)*8+index,j] <- index.m[(i-1)*8+index,j] + vphi[k]
            }
          }
          if(length(AA1)==16){
            phi1 <- eval(parse(text=AA1[1]))*eval(parse(text=AA1[2]))/
              (eval(parse(text=AA1[1]))*eval(parse(text=AA1[2]))+eval(parse(text=AA1[3]))*eval(parse(text=AA1[4]))+eval(parse(text=AA1[5]))*eval(parse(text=AA1[6]))+eval(parse(text=AA1[7]))*eval(parse(text=AA1[8]))+
                 eval(parse(text=AA1[9]))*eval(parse(text=AA1[10]))+eval(parse(text=AA1[11]))*eval(parse(text=AA1[12]))+eval(parse(text=AA1[13]))*eval(parse(text=AA1[14]))+eval(parse(text=AA1[15]))*eval(parse(text=AA1[16])))
            phi2 <- eval(parse(text=AA1[3]))*eval(parse(text=AA1[4]))/
              (eval(parse(text=AA1[1]))*eval(parse(text=AA1[2]))+eval(parse(text=AA1[3]))*eval(parse(text=AA1[4]))+eval(parse(text=AA1[5]))*eval(parse(text=AA1[6]))+eval(parse(text=AA1[7]))*eval(parse(text=AA1[8]))+
                 eval(parse(text=AA1[9]))*eval(parse(text=AA1[10]))+eval(parse(text=AA1[11]))*eval(parse(text=AA1[12]))+eval(parse(text=AA1[13]))*eval(parse(text=AA1[14]))+eval(parse(text=AA1[15]))*eval(parse(text=AA1[16])))
            phi3 <- eval(parse(text=AA1[5]))*eval(parse(text=AA1[6]))/
              (eval(parse(text=AA1[1]))*eval(parse(text=AA1[2]))+eval(parse(text=AA1[3]))*eval(parse(text=AA1[4]))+eval(parse(text=AA1[5]))*eval(parse(text=AA1[6]))+eval(parse(text=AA1[7]))*eval(parse(text=AA1[8]))+
                 eval(parse(text=AA1[9]))*eval(parse(text=AA1[10]))+eval(parse(text=AA1[11]))*eval(parse(text=AA1[12]))+eval(parse(text=AA1[13]))*eval(parse(text=AA1[14]))+eval(parse(text=AA1[15]))*eval(parse(text=AA1[16])))
            phi4 <- eval(parse(text=AA1[7]))*eval(parse(text=AA1[8]))/
              (eval(parse(text=AA1[1]))*eval(parse(text=AA1[2]))+eval(parse(text=AA1[3]))*eval(parse(text=AA1[4]))+eval(parse(text=AA1[5]))*eval(parse(text=AA1[6]))+eval(parse(text=AA1[7]))*eval(parse(text=AA1[8]))+
                 eval(parse(text=AA1[9]))*eval(parse(text=AA1[10]))+eval(parse(text=AA1[11]))*eval(parse(text=AA1[12]))+eval(parse(text=AA1[13]))*eval(parse(text=AA1[14]))+eval(parse(text=AA1[15]))*eval(parse(text=AA1[16])))
            phi5 <- eval(parse(text=AA1[9]))*eval(parse(text=AA1[10]))/
              (eval(parse(text=AA1[1]))*eval(parse(text=AA1[2]))+eval(parse(text=AA1[3]))*eval(parse(text=AA1[4]))+eval(parse(text=AA1[5]))*eval(parse(text=AA1[6]))+eval(parse(text=AA1[7]))*eval(parse(text=AA1[8]))+
                 eval(parse(text=AA1[9]))*eval(parse(text=AA1[10]))+eval(parse(text=AA1[11]))*eval(parse(text=AA1[12]))+eval(parse(text=AA1[13]))*eval(parse(text=AA1[14]))+eval(parse(text=AA1[15]))*eval(parse(text=AA1[16])))
            phi6 <- eval(parse(text=AA1[11]))*eval(parse(text=AA1[12]))/
              (eval(parse(text=AA1[1]))*eval(parse(text=AA1[2]))+eval(parse(text=AA1[3]))*eval(parse(text=AA1[4]))+eval(parse(text=AA1[5]))*eval(parse(text=AA1[6]))+eval(parse(text=AA1[7]))*eval(parse(text=AA1[8]))+
                 eval(parse(text=AA1[9]))*eval(parse(text=AA1[10]))+eval(parse(text=AA1[11]))*eval(parse(text=AA1[12]))+eval(parse(text=AA1[13]))*eval(parse(text=AA1[14]))+eval(parse(text=AA1[15]))*eval(parse(text=AA1[16])))
            phi7 <- eval(parse(text=AA1[13]))*eval(parse(text=AA1[14]))/
              (eval(parse(text=AA1[1]))*eval(parse(text=AA1[2]))+eval(parse(text=AA1[3]))*eval(parse(text=AA1[4]))+eval(parse(text=AA1[5]))*eval(parse(text=AA1[6]))+eval(parse(text=AA1[7]))*eval(parse(text=AA1[8]))+
                 eval(parse(text=AA1[9]))*eval(parse(text=AA1[10]))+eval(parse(text=AA1[11]))*eval(parse(text=AA1[12]))+eval(parse(text=AA1[13]))*eval(parse(text=AA1[14]))+eval(parse(text=AA1[15]))*eval(parse(text=AA1[16])))
            phi8 <- eval(parse(text=AA1[15]))*eval(parse(text=AA1[16]))/
              (eval(parse(text=AA1[1]))*eval(parse(text=AA1[2]))+eval(parse(text=AA1[3]))*eval(parse(text=AA1[4]))+eval(parse(text=AA1[5]))*eval(parse(text=AA1[6]))+eval(parse(text=AA1[7]))*eval(parse(text=AA1[8]))+
                 eval(parse(text=AA1[9]))*eval(parse(text=AA1[10]))+eval(parse(text=AA1[11]))*eval(parse(text=AA1[12]))+eval(parse(text=AA1[13]))*eval(parse(text=AA1[14]))+eval(parse(text=AA1[15]))*eval(parse(text=AA1[16])))
            vphi <- c(phi1,phi1,phi2,phi2,phi3,phi3,phi4,phi4,phi5,phi5,phi6,phi6,
                      phi7,phi7,phi8,phi8)
            for(k in 1:length(AA1)){
              index <- which(AA1[k]==index.g)
              index.m[(i-1)*8+index,j] <- index.m[(i-1)*8+index,j] + vphi[k]
            }
          }
        }
      }
      new.matrix <- n.matrix*index.m
      new.matrix1 <- cbind(new.matrix[1:8,],new.matrix[9:16,],new.matrix[17:24,],new.matrix[25:32,],
                           new.matrix[33:40,],new.matrix[41:48,],new.matrix[49:56,],new.matrix[57:64,],
                           new.matrix[65:72,])
      g <- rowSums(new.matrix1)/(2*n)
      g000 <- g[1];g001 <- g[2];g010 <- g[3];g011 <- g[4];
      g100 <- g[5];g101 <- g[6];g110 <- g[7];g111 <- g[8]
      
      cat("Inter:",loop_k,"\n")
      loop_k <- loop_k + 1
    }
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
  list(rAB=rAB,rBC=rBC,rCD=rCD,rAC=rAC,rBD=rBD,rAD=rAD,C4=C4,C1=C1,C2=C2,C3=C3)
}


Two.F2 <- function(dat){
  
  datt <- dat.tran1(dat,type="F2")
  
  r <- rf.2pts(datt,verbose = F)
  
  r.m <- apply(r$analysis[,,1],1,A<- function(x){x[which(x==min(x))][1]})
  rAB <- r.m[1]
  rAC <- r.m[2]
  rBC <- r.m[3]
  rAD <- r.m[4]
  rBD <- r.m[5]
  rCD <- r.m[6]
  
  list(rAB=rAB,rBC=rBC,rCD=rCD,rAC=rAC,rBD=rBD,rAD=rAD)
}

linkage.phase <- function(lpindex=1){
  
  lp <- list()
  prob1 <- matrix(NA,9,9)
  prob1[1,] <- c("g000-g000","g000-g001","g001-g001","g000-g011","g000-g010-g001-g011",
                 "g001-g010","g011-g011","g011-g010","g010-g010")
  prob1[2,] <- c("g000-g110","g000-g111-g001-g110","g001-g111","g000-g101-g011-g110",
                 "g000-g100-g001-g101-g010-g110-g011-g111","g001-g100-g010-g111",
                 "g011-g101","g010-g101-g011-g100","g010-g100")
  prob1[3,] <- c("g110-g110","g110-g111","g111-g111","g110-g101","g110-g100-g111-g101",
                 "g111-g100","g101-g101","g101-g100","g100-g100")
  prob1[4,] <- c("g000-g100","g000-g101-g001-g100","g001-g101","g000-g111-g011-g100",
                 "g000-g110-g001-g111-g010-g100-g011-g101","g001-g110-g010-g101",
                 "g011-g111","g010-g111-g011-g110","g010-g110")
  prob1[5,] <- c("g000-g010-g110-g100","g000-g011-g001-g010-g110-g101-g111-g100",
                 "g001-g011-g111-g101",
                 "g000-g001-g011-g010-g110-g111-g101-g100",
                 "g000-g000-g001-g001-g010-g010-g011-g011-g100-g100-g101-g101-g110-g110-g111-g111",
                 "g000-g001-g011-g010-g110-g111-g101-g100","g001-g011-g111-g101",
                 "g000-g011-g001-g010-g110-g101-g111-g100","g000-g010-g110-g100")
  prob1[6,] <- c("g010-g110","g010-g111-g011-g110","g011-g111","g001-g110-g010-g101",
                 "g000-g110-g001-g111-g010-g100-g011-g101","g000-g111-g011-g100",
                 "g001-g101","g000-g101-g001-g100","g000-g100")
  prob1[7,] <- c("g100-g100","g101-g100","g101-g101","g111-g100","g110-g100-g111-g101",
                 "g110-g101","g111-g111","g110-g111","g110-g110")
  prob1[8,] <- c("g010-g100","g010-g101-g011-g100","g011-g101","g001-g100-g010-g111",
                 "g000-g100-g001-g101-g010-g110-g011-g111","g000-g101-g011-g110",
                 "g001-g111","g000-g111-g001-g110","g000-g110")
  prob1[9,] <- c("g010-g010","g011-g010","g011-g011","g001-g010","g000-g010-g001-g011",
                 "g000-g011","g001-g001","g000-g001","g000-g000")
  lp[[1]] <- prob1
  return(lp[[lpindex]])
}

dat.tran1 <- function(dat,type="BC"){
  
  nmarker <- c()
  for(i in 1:4){
    
    A <- substr(dat,(2*i-1),2*i)
    AA <- names(table(A))
    for(j in 1:length(AA)){
      A[which(A==AA[j])] <- j-1
    }
    m.tmp <- as.numeric(A)
    nmarker <- cbind(nmarker,m.tmp)
  }
  colnames(nmarker) <- c(paste("M",1:4,sep=""))
  if(type=="BC"){
    segr.type=rep("D2.15",4)
    segr.type.num=rep(7,4)
  }
  if(type=="F2"){
    segr.type=rep("B3.7",4)
    segr.type.num=rep(4,4)
  }
  dat1 <- structure(list(n.ind=length(dat),n.mar=4,geno=nmarker+1,segr.type=segr.type
                         ,segr.type.num=segr.type.num),class = "outcross")
  return(dat1)
}

two.F2 <- function(dat){
  
  datt <- dat.tran1(dat,type="F2")
  
  rAB <-  def.rf.3pts(datt,"M2","M3","M4")
  
  
  r <-  rf.2pts(datt)
}