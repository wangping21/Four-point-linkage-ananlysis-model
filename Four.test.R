setwd("/home/dudusdau/ccb/Four-linkage")

library(onemap)
source("Four.sim.R")
source("Four.est.R")
source("Four.est1.R")
source("Four.util.R")

par_st <- par_obj()
####simulation
###################################BC################################
nS <- c()
for(ii in 1:200){
  BC_gam <- simu_BC(par_obj=par_st,n=400)
  S <- Four.BC(dat=BC_gam)
  nS <- rbind(nS,unlist(S))
}

nS1 <- c()
for(ii in 1:200){
  BC_gam <- simu_BC(par_obj=par_st,n=400)
  S <- Two.BC(dat=BC_gam)
  nS1 <- rbind(nS1,unlist(S))
}

nS2 <- c()
for(ii in 1:200){
  BC_gam <- simu_BC(par_obj=par_st,n=400)
  S <- Three.BC(dat=BC_gam)
  nS2 <- rbind(nS2,unlist(S))
}


###################################F2################################

nS <- c()
for(ii in 1:100){
  F2_gam <- simu_F2(par_obj=par_st,n=400)
  S <- Four.F2(dat=F2_gam,init.g=c(0.5345,0.2895,0.1395,-0.0045,0.0155,0.0105,0.0105,0.0045))
  nS <- rbind(nS,unlist(S))
}



nS1 <- c()
for(ii in 1:100){
  F2_gam <- simu_F2(par_obj=par_st,n=400)
  S <- Two.F2(dat=F2_gam)
  nS1 <- rbind(nS1,unlist(S))
}


nS2 <- c()
for(ii in 1:100){
  F2_gam <- simu_F2(par_obj=par_st,n=400)
  S <- Three.F2(dat=F2_gam)
  nS2 <- rbind(nS2,unlist(S))
}


#Work example
##################################################################################
load("data.RData")

datt <- dat
datt$geno <- dat$geno[,c(247,238,242,236)];datt$n.mar <- 4;datt$segr.type <- dat$segr.type[c(247,238,242,236)];
datt$segr.type.num <- dat$segr.type.num[c(247,238,242,236)]
twopts <- rf.2pts(datt,LOD=3,max.rf=0.5) # perform two-point analyses
r.m <- apply(twopts$analysis[,,1],1,A<- function(x){x[which(x==min(x))][1]})
rAB <- r.m[1]
rAC <- r.m[2]
rBC <- r.m[3]
rAD <- r.m[4]
rBD <- r.m[5]
rCD <- r.m[6]


fourpts <- four.est(datt,lp = c(1:64))

select.1 <- fourpts[which(fourpts[,12]==max(fourpts[,12])),]
newx <- c()
for(i in 1:dim(select.1)[1]){
  if(all(select.1[i,1:6]<0.5)){
    newx <- c(newx,i)
  }
}

final <- select.1[newx[1],]





