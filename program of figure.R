==========结果1============================
setwd("C:/data")
outcross <- read.csv("res11.csv")
map <- read.csv("Populus euphratica.outcross.csv",header = F)
lkgroup <- c(which(map[,3]==0))
xmax <- sum(map[lkgroup[2]-1,3],map[lkgroup[3]-1,3],map[lkgroup[4]-1,3],map[lkgroup[5]-1,3],map[lkgroup[6]-1,3],map[lkgroup[7]-1,3],map[lkgroup[8]-1,3],
            map[lkgroup[9]-1,3],map[lkgroup[10]-1,3],map[lkgroup[11]-1,3],map[lkgroup[12]-1,3],map[lkgroup[13]-1,3],map[lkgroup[14]-1,3],map[lkgroup[15]-1,3],
            map[lkgroup[16]-1,3],map[lkgroup[17]-1,3],map[lkgroup[18]-1,3],map[lkgroup[19]-1,3],map[nrow(map),3])
xmap <- c(map[lkgroup[1],3],map[lkgroup[2]-1,3],map[lkgroup[3]-1,3],map[lkgroup[4]-1,3],map[lkgroup[5]-1,3],map[lkgroup[6]-1,3],map[lkgroup[7]-1,3],map[lkgroup[8]-1,3],
          map[lkgroup[9]-1,3],map[lkgroup[10]-1,3],map[lkgroup[11]-1,3],map[lkgroup[12]-1,3],map[lkgroup[13]-1,3],map[lkgroup[14]-1,3],map[lkgroup[15]-1,3],map[lkgroup[16]-1,3],
          map[lkgroup[17]-1,3],map[lkgroup[18]-1,3],map[lkgroup[19]-1,3],map[nrow(map),3])
mapmap <- cumsum(xmap)
#pdf("figure1.pdf")
plot(0,-8,xlim=c(-250,xmax),ylim=c(-100,320),axes=FALSE,xlab="",ylab="",type="l")
segments(0,-8,xmax,-8,lwd = 2)
segments(0,320,xmax,320,lwd = 4)
segments(0,-8,0,320,lwd = 4)
segments(xmax,-8,xmax,-8,lwd = 4)
segments(0,1,xmax,1,lty="dotted",lwd = 2)

segments(-50,0,0,0);text(-150,0,cex=1.7,srt=90,"0")
#segments(-50,40,0,40);text(-150,40,cex=1.5,srt=90,"40")
segments(-50,80,0,80);text(-150,80,cex=1.7,srt=90,"80")
#segments(-50,120,0,120);text(-150,125,cex=1.5,srt=90,"120")
segments(-50,160,0,160);text(-150,160,cex=1.7,srt=90,"160")
#segments(-50,200,0,200);text(-150,200,cex=1.5,srt=90,"200")
segments(-50,240,0,240);text(-150,240,cex=1.7,srt=90,"240")
#segments(-50,280,0,280);text(-150,280,cex=1.5,srt=90,"280")
segments(-50,320,0,320);text(-150,320,cex=1.7,srt=90,"320")
#####-----marker-----#####
for(i in 1:(lkgroup[2]-1)){
  segments(map[i,3],-8,map[i,3],-2,lwd=0.000001)
}

for(i in lkgroup[2]:(lkgroup[3]-1)){
  segments(map[i,3]+mapmap[2],-8,map[i,3]+mapmap[2],-2,lwd=0.000001)
}

for(i in lkgroup[3]:(lkgroup[4]-1)){
  segments(map[i,3]+mapmap[3],-8,map[i,3]+mapmap[3],-2,lwd=0.000001)
}

for(i in lkgroup[4]:(lkgroup[5]-1)){
  segments(map[i,3]+mapmap[4],-8,map[i,3]+mapmap[4],-2,lwd=0.000001)
}

for(i in lkgroup[5]:(lkgroup[6]-1)){
  segments(map[i,3]+mapmap[5],-8,map[i,3]+mapmap[5],-2,lwd=0.000001)
}

for(i in lkgroup[6]:(lkgroup[7]-1)){
  segments(map[i,3]+mapmap[6],-8,map[i,3]+mapmap[6],-2,lwd=0.000001)
}

for(i in lkgroup[7]:(lkgroup[8]-1)){
  segments(map[i,3]+mapmap[7],-8,map[i,3]+mapmap[7],-2,lwd=0.000001)
}


for(i in lkgroup[8]:(lkgroup[9]-1)){
  segments(map[i,3]+mapmap[8],-8,map[i,3]+mapmap[8],-2,lwd=0.000001)
}

for(i in lkgroup[9]:(lkgroup[10]-1)){
  segments(map[i,3]+mapmap[9],-8,map[i,3]+mapmap[9],-2,lwd=0.000001)
}

for(i in lkgroup[10]:(lkgroup[11]-1)){
  segments(map[i,3]+mapmap[10],-8,map[i,3]+mapmap[10],-2,lwd=0.000001)
}

for(i in lkgroup[11]:(lkgroup[12]-1)){
  segments(map[i,3]+mapmap[11],-8,map[i,3]+mapmap[11],-2,lwd=0.000001)
}

for(i in lkgroup[12]:(lkgroup[13]-1)){
  segments(map[i,3]+mapmap[12],-8,map[i,3]+mapmap[12],-2,lwd=0.000001)
}

for(i in lkgroup[13]:(lkgroup[14]-1)){
  segments(map[i,3]+mapmap[13],-8,map[i,3]+mapmap[13],-2,lwd=0.000001)
}
for(i in lkgroup[14]:(lkgroup[15]-1)){
  segments(map[i,3]+mapmap[14],-8,map[i,3]+mapmap[14],-2,lwd=0.000001)
}

for(i in lkgroup[15]:(lkgroup[16]-1)){
  segments(map[i,3]+mapmap[15],-8,map[i,3]+mapmap[15],-2,lwd=0.000001)
}

for(i in lkgroup[16]:(lkgroup[17]-1)){
  segments(map[i,3]+mapmap[16],-8,map[i,3]+mapmap[16],-2,lwd=0.000001)
}

for(i in lkgroup[17]:(lkgroup[18]-1)){
  segments(map[i,3]+mapmap[17],-8,map[i,3]+mapmap[17],-2,lwd=0.000001)
}
for(i in lkgroup[18]:(lkgroup[19]-1)){
  segments(map[i,3]+mapmap[18],-8,map[i,3]+mapmap[18],-2,lwd=0.000001)
}

for(i in lkgroup[19]:(nrow(map))){
  segments(map[i,3]+mapmap[19],-8,map[i,3]+mapmap[19],-2,lwd=0.000001)
}

######-----linkage group-------########
for(i in 1:19){
  segments((mapmap[i]+(mapmap[i+1]-mapmap[i])/2),-8,(mapmap[i]+(mapmap[i+1]-mapmap[i])/2),-11)
}
text((mapmap[1]+(mapmap[2]-mapmap[1])/2),-20,"1",cex = 1.5)
text((mapmap[2]+(mapmap[3]-mapmap[2])/2),-20,"2",cex = 1.5)
text((mapmap[3]+(mapmap[4]-mapmap[3])/2),-20,"3",cex = 1.5)
text((mapmap[4]+(mapmap[5]-mapmap[4])/2),-20,"4",cex = 1.5)
text((mapmap[5]+(mapmap[6]-mapmap[5])/2),-20,"5",cex = 1.5)
text((mapmap[6]+(mapmap[7]-mapmap[6])/2),-20,"6",cex = 1.5)
text((mapmap[7]+(mapmap[8]-mapmap[7])/2),-20,"7",cex = 1.5)
text((mapmap[8]+(mapmap[9]-mapmap[8])/2),-20,"8",cex = 1.5)

text((mapmap[9]+(mapmap[10]-mapmap[9])/2),-20,"9",cex = 1.5)
text((mapmap[10]+(mapmap[11]-mapmap[10])/2),-20,"10",cex = 1.5)
text((mapmap[11]+(mapmap[12]-mapmap[11])/2),-20,"11",cex = 1.5)
text((mapmap[12]+(mapmap[13]-mapmap[12])/2),-20,"12",cex = 1.5)
text((mapmap[13]+(mapmap[14]-mapmap[13])/2),-20,"13",cex = 1.5)
text((mapmap[14]+(mapmap[15]-mapmap[14])/2),-20,"14",cex = 1.5)
text((mapmap[15]+(mapmap[16]-mapmap[15])/2),-20,"15",cex = 1.5)
text((mapmap[16]+(mapmap[17]-mapmap[16])/2),-20,"16",cex = 1.5)
text((mapmap[17]+(mapmap[18]-mapmap[17])/2),-20,"17",cex = 1.5)
text((mapmap[18]+(mapmap[19]-mapmap[18])/2),-20,"18",cex = 1.5)
text((mapmap[19]+(mapmap[20]-mapmap[19])/2),-20,"19",cex = 1.5)
text(2300,-40,cex=2,"chromosome")
text(-300,160,srt=90,cex = 1.7,"Coefficients of coincidence")

outcross_map <- matrix(NA,nrow(outcross),4)
outcross_C <- outcross[,8:11]
for(j in 1:4){	
  for(i in 1:nrow(outcross)){
    if(outcross_C[i,j]>300){
      outcross_C[i,j] <- 300
    }
  }
}

l1 <- map[1:lkgroup[2]-1,3]
l2 <- map[lkgroup[2]:(lkgroup[3]-1),3]+mapmap[2]
l3 <- map[lkgroup[3]:(lkgroup[4]-1),3]+mapmap[3]
l4 <- map[lkgroup[4]:(lkgroup[5]-1),3]+mapmap[4]
l5 <- map[lkgroup[5]:(lkgroup[6]-1),3]+mapmap[5]
l6 <- map[lkgroup[6]:(lkgroup[7]-1),3]+mapmap[6]
l7 <- map[lkgroup[7]:(lkgroup[8]-1),3]+mapmap[7]
l8 <- map[lkgroup[8]:(lkgroup[9]-1),3]+mapmap[8]


l9 <- map[lkgroup[9]:(lkgroup[10]-1),3]+mapmap[9]
l10 <- map[lkgroup[10]:(lkgroup[11]-1),3]+mapmap[10]
l11 <- map[lkgroup[11]:(lkgroup[12]-1),3]+mapmap[11]
l12<- map[lkgroup[12]:(lkgroup[13]-1),3]+mapmap[12]
l13<- map[lkgroup[13]:(lkgroup[14]-1),3]+mapmap[13]
l14<- map[lkgroup[14]:(lkgroup[15]-1),3]+mapmap[14]
l15 <- map[lkgroup[15]:(lkgroup[16]-1),3]+mapmap[15]
l16 <- map[lkgroup[16]:(lkgroup[17]-1),3]+mapmap[16]
l17<- map[lkgroup[17]:(lkgroup[18]-1),3]+mapmap[17]
l18<- map[lkgroup[18]:(lkgroup[19]-1),3]+mapmap[18]
l19 <- map[lkgroup[19]:nrow(map),3]+mapmap[19]

C1 <- c(l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19)
x0 <- c()
y0 <- c()
for(i in 1:nrow(outcross_C)){
  #points(C1[3+4*(i-1)],outcross_C[i,2],pch=20,cex=0.1,col="red")
  x0 <- c(x0,C1[3+4*(i-1)])
  y0 <- c(y0,outcross_C[i,4])
}

#x1 <- c()
#y1 <- c()
#for(i in 1:nrow(outcross_C)){
  #points(C1[2+4*(i-1)],outcross_C[i,1],pch=20,cex=0.1,col="red")
  #x1 <- c(x1,C1[1+4*(i-1)])
  #y1 <- c(y1,outcross_C[i,1])
#}

#x2 <- c()
#y2 <- c()
#for(i in 1:nrow(outcross_C)){
  #points(C1[3+4*(i-1)],outcross_C[i,2],pch=20,cex=0.1,col="red")
  #x2 <- c(x2,C1[3+4*(i-1)])
  #y2 <- c(y2,outcross_C[i,3])
#}

#x <- c()
#for(i in 1:nrow(outcross_C)){
# x <- c(x,c(x1[i],x2[i]))
#}
#y <- c()
#for(i in 1:nrow(outcross_C)){
#y <- c(y,c(y1[i],y2[i]))
#}

lines(x0,y0,col=c("green",alpha=5),lwd=1.5)
#lines(x1,y1,col="blue",lwd=1.5)
#lines(x2,y2,col="red",lwd=1.5)
for(i in 2:20){
  segments(mapmap[i],-8,mapmap[i],320,lwd = 4)
}
#lines(x2,y2,col="yellow",lwd=0.5)
#dev.off()

==============结果2================================
setwd("C:/data")
data <- read.csv("res1.csv",header = T)
map <- read.csv("Populus euphratica.outcross.csv",header = F)
lkgroup <- which(map[,3]==0)
data[data>100] <- 0;
lk1_r <- c(data[which(data[,12]==1),2],data[which(data[,12]==1),3],data[which(data[,12]==1),4])
lk2_r <- c(data[which(data[,12]==2),2],data[which(data[,12]==2),3],data[which(data[,12]==2),4])
lk3_r <- c(data[which(data[,12]==3),2],data[which(data[,12]==3),3],data[which(data[,12]==3),4])
lk4_r <- c(data[which(data[,12]==4),2],data[which(data[,12]==4),3],data[which(data[,12]==4),4])
lk5_r <- c(data[which(data[,12]==5),2],data[which(data[,12]==5),3],data[which(data[,12]==5),4])
lk6_r <- c(data[which(data[,12]==6),2],data[which(data[,12]==6),3],data[which(data[,12]==6),4])
lk7_r <- c(data[which(data[,12]==7),2],data[which(data[,12]==7),3],data[which(data[,12]==7),4])
lk8_r <- c(data[which(data[,12]==8),2],data[which(data[,12]==8),3],data[which(data[,12]==9),4])
lk9_r <- c(data[which(data[,12]==9),2],data[which(data[,12]==9),3],data[which(data[,12]==9),4])
lk10_r <- c(data[which(data[,12]==10),2],data[which(data[,12]==10),3],data[which(data[,12]==10),4])
lk11_r <- c(data[which(data[,12]==11),2],data[which(data[,12]==11),3],data[which(data[,12]==11),4])
lk12_r <- c(data[which(data[,12]==12),2],data[which(data[,12]==12),3],data[which(data[,12]==12),4])
lk13_r <- c(data[which(data[,12]==13),2],data[which(data[,12]==13),3],data[which(data[,12]==13),4])
lk14_r <- c(data[which(data[,12]==14),2],data[which(data[,12]==14),3],data[which(data[,12]==14),4])
lk15_r <- c(data[which(data[,12]==15),2],data[which(data[,12]==15),3],data[which(data[,12]==15),4])
lk16_r <- c(data[which(data[,12]==16),2],data[which(data[,12]==16),3],data[which(data[,12]==16),4])
lk17_r <- c(data[which(data[,12]==17),2],data[which(data[,12]==17),3],data[which(data[,12]==17),4])
lk18_r <- c(data[which(data[,12]==18),2],data[which(data[,12]==18),3],data[which(data[,12]==18),4])
lk19_r <- c(data[which(data[,12]==19),2],data[which(data[,12]==19),3],data[which(data[,12]==19),4])
r.means <- c(mean(lk1_r),mean(lk2_r),mean(lk3_r),mean(lk4_r),mean(lk5_r),mean(lk6_r),mean(lk7_r),mean(lk8_r),
             mean(lk9_r),mean(lk10_r),mean(lk11_r),mean(lk12_r),mean(lk13_r),mean(lk14_r),mean(lk15_r),
             mean(lk16_r),mean(lk17_r),mean(lk18_r),mean(lk19_r))
dis1 <- c(map[lkgroup[2]-1,3],map[lkgroup[3]-1,3],map[lkgroup[4]-1,3],map[lkgroup[5]-1,3],map[lkgroup[6]-1,3],map[lkgroup[7]-1,3],map[lkgroup[8]-1,3],
          map[lkgroup[9]-1,3],map[lkgroup[10]-1,3],map[lkgroup[11]-1,3],map[lkgroup[12]-1,3],map[lkgroup[13]-1,3],map[lkgroup[14]-1,3],map[lkgroup[15]-1,3],
          map[lkgroup[16]-1,3],map[lkgroup[17]-1,3],map[lkgroup[18]-1,3],map[lkgroup[19]-1,3],map[nrow(map),3])
#b <- read.csv("r.means.csv",header = T)
#z <- lm(b[,3]~b[,4]);
#lines(b[,4],fitted(z))
plot(dis1 ,r.means)
#axis(2,seq(0.10,0.17,0.01),seq(0.10,0.17,0.01))
b <- read.csv("r.means.csv",header = T)
x <- as.vector(b[1:9,4]);
y <-  as.vector(b[1:9,3]);
z <- lm(y~x);
#lines(x,fitted(z))
x1 <- as.vector(b[1:8,2]);
y1 <-  as.vector(b[1:8,1]);
z1 <- lm(y1~x1);
#lines(x1,fitted(z1))
x <- c(0:300)
predict(z,newdata=x)
lines(x,0.0005077*x+0.0405563,col="red")
x1 <- c(0:400)
lines(x1,0.0002399*x1+0.0612807,col="blue")
============结果三（matlab）==================================================
data = xlsread('C:\data\res1.csv');
regroupratio = data(:,11:12);
differences = zeros(19,3);
maxvalue = 100;
for i = 1:19
type1 = regroupratio(regroupratio(:,2)==i);
type1 =type1(type1<100 & type1>0.1);
num = size(type1,1);
num = floor(num/3);
before = type1(1:num);
center = type1((num+1):(2*num));
after = type1((2*num+1):(3*num));
before1 = [];
center1 =[];
after1 = [];
for j = 10:10:maxvalue
before1 = [ before1 sum(before>(j-10) & before<j) ] ;
center1 = [ center1 sum(center>(j-10) & center<j) ] ;
after1 = [ after1 sum(after>(j-10) & after<j) ] ;
end
differences(i,1) = sum(abs(before1-center1))/(num);
differences(i,2) = sum(abs(after1-center1))/(num);
differences(i,3) = sum(abs(after1-before1))/(num);

end



