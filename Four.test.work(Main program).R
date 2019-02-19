setwd("/home/dudusdau/ccb/Four-linkage")

library(onemap)
source("Four.sim.R")
source("Four.est.R")
source("Four.est1.R")
source("Four.util.R")

load("data.RData")
map.tmp <- read.table("Populus euphratica.outcross.map")
nlg <- length(table(map.tmp[,1]))
iter <- 0
res1 <- c()
for(i in 1:nlg){
  
  lg <- which(map.tmp[,1]==i)
  nt <- length(lg)
  nnt <- floor(nt/4)
  res <- c()
  for(j in 1:nnt){
    index <- as.character(map.tmp[lg[(4*j-3):(4*j)],2])
    index1 <- c()
    for(k in c(index)){
      index1 <- c(index1,which(k==colnames(dat$geno)))
    }
    datt <- dat
    datt$geno <- dat$geno[,index1];datt$n.mar <- 4;datt$segr.type <- dat$segr.type[index1];
    datt$segr.type.num <- dat$segr.type.num[index1]
    fourpts <- four.est(datt,lp = c(1:64))
    
    select.1 <- fourpts[which(fourpts[,12]==max(fourpts[,12])),]
    newx <- c()
    for(kk in 1:dim(select.1)[1]){
      if(all(select.1[kk,1:6]<0.55)){
        newx <- c(newx,kk)
      }
    }
    if(is.null(newx))
      break;
    
    final <- c(select.1[newx[1],],lg=i,M1=index[1],M2=index[2],M3=index[3],M4=index[4])
    res <- rbind(res,final)
    iter <- iter + 1
    cat("iter =",iter,final,"\n")
    
  }
  
  res1 <- rbind(res1,res)
  
}

