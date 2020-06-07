#specify your folder 
setwd("C:/Programs/Low-dim mat files")

# Libraries -----------------------------------------------
library(factoextra)
library(imputeTS)
library(readr)
library(tsbox)
library(base)
library(timeSeries)
library(stats)
library(rgl)
library(randomForest)
library(caret)
library(EpiModel)
library(R.matlab)

# all functions --------------------------------------------

#normalize 
normalizef <- function(x) {
  return ((x - min(x,na.rm=T)) / (max(x,na.rm=T) - min(x,na.rm=T)))
}

#Elbow (Wilhelm)
findClusters <- function(f) {
  asw <- numeric(10)
  for (k in 2:10)
    asw[[k]] <- pam(f,k)$silinfo$avg.width
  
  k.best <- which.max(asw)
  #  cl <- kmeans(f,k.best)
  return(k.best)
}

#mode function 
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#mode function - 0 if less than 50% agreement
getmode2 <- function(v) {
  uniqv <- unique(v)
  moda<-max(tabulate(match(v, uniqv)))
  if(moda>=5){
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  else{return(0)}
}


#gap_stat 
number_clusters<-function (dissimilarity){
  fclust<-fviz_nbclust(dissimilarity, hcut, method = "gap_stat",nboot=30,
                       hc_method = "complete")
  nclust<-which(diff(fclust$data[,"gap"])<0)[1]
  return(nclust)
}


removeNoise<-function(th,labelhc,nclust){
  newlabelhc<-matrix(0,nrow=length(labelhc),ncol=1)
  for (j in 1:nclust){
    icluster<-as.numeric(labelhc==j)
    #icluster<-smooth(icluster)
    icluster2<-matrix(0,nrow=length(icluster),ncol=1)
     for (k in seq(1,length(icluster),th)){
       icluster2[k:(k+th)]<-getmode2(icluster[k:(k+th)])  
     }
    icluster2<-icluster    
    changeLabel<-cumsum(abs(diff(icluster2)))
    A<-table(changeLabel)>th
    B<-as.numeric(names(A)[A==T])  
    if (length(B)>1){
      C<-match(B,changeLabel)
      #lines(C,rep(.5,length(C)),type="p")
      new<-rep(0,length(icluster))
      
      for (l in 1:(length(C)-1)){
        new[C[l]:C[l+1]]<-getmode2(icluster2[C[l]:C[l+1]])
      }
      newlabelhc[new==1]<-j
    }
    
  }
  return(newlabelhc)}

findSegments<-function(timePoints,newlabelhc){
  th2<-(max(carFrame$timePoint))/60
  Q<- approx(carFrame$timePoint,newlabelhc,n=ceiling((max(carFrame$timePoint))/th2),method="constant")
  lapses<-cbind(Q$x,Q$y)
  
  grides<-abs(diff(Q$y))
  cutTime<-c(0,Q$x[grides!=0],max(carFrame$timePoint))
  
  E<-rle(Q$y)
  po0<-which(E$values==0)
  # timelaps0<-matrix(0,nrow=length(po0),ncol=2)
  # for (j in 1:length(po0)){
  #   if(j==1){
  #     if (po0[1]==1){
  #       timelaps0[1,1]<-1
  #     timelaps0[1,2]<-E$length[1]
  #     }else{
  #       timelaps0[1,1]<-cumsum(E$lengths)[po0[1]]
  #       timelaps0[1,2]<-timelaps0[1,1]+E$lengths[1]-1
  #       
  #     }
  #   }else{
  #     timelaps0[j,1]<-cumsum(E$lengths)[po0[j]]
  #     timelaps0[j,2]<-timelaps0[j,1]+cumsum(E$lengths)[po0[j]]-1
  #   }
  # }
  
  # timelaps0<-cbind(timelaps0,(timelaps0+(rle(Q$y)$lengths[rle(Q$y)$values==0])-1))*10
  segmentTime<-list("lapses"=lapses,"time0"=po0,"cutTime"=cutTime)
  return(segmentTime)
}

#method 2 (Tobias)
knee2Line<-function(fclust){
  k_knee<-fclust
  k_knee = k_knee$data; y_knee = k_knee[,2]; x_knee = 1:length(y_knee); xx=x_knee; yy=y_knee # create xx,yy just for plotting
  
  E = c(Inf)
for(i in 2:(length(x_knee)-1)){
  e = 0
  k1 = (y_knee[i]-y_knee[1])/(x_knee[i]-x_knee[1])  # y(i) = k1x(i) + m1
  m1 = y_knee[1]-k1*x_knee[1]
  k10 = (y_knee[i]-y_knee[length(x_knee)])/(x_knee[i]-x_knee[length(x_knee)])  # y(i) = k10x(i) + m10
  m10 = y_knee[length(x_knee)]-k10*x_knee[length(x_knee)]
  for(j in 2:i){
    e = e + abs(k1*x_knee[j]+m1-y_knee[j])
  }
  for(j in (i+1):length(x_knee)){
    e = e + abs(k10*x_knee[j]+m10-y_knee[j])
  }
  E[i] = e 
}
indices = which.min(E)
k_knee = round(x_knee[indices])
return(k_knee)}


# Opening/saving data ---------------------------------------------------------------
number_Ref<-1
name2<-(4992+(number_Ref-1)*60)
name2<-toString(formatC(name2,width=5,format="d",flag="0"))
name3<-".mat"
Ref_name<-paste(c(name2,name3),collapse="")

allFile<-readMat(Ref_name)
varNames<-rownames(allFile$iState)

#size of file 
nObjects<-as.numeric(ncol(allFile$iState[[1]]))
nVars<-length(varNames)

#matrix [variable x objects], save which ones are completely NA 
missingInfo<-matrix(data=NA,nrow=nVars,ncol=nObjects)
missingInfo<-data.frame(missingInfo,row.names=varNames) 

for(q in 1:nVars){
  iState<-as.matrix(as.data.frame(allFile$iState[[q]]))
  Tn<-nrow(iState)
  missingInfo[q,]<-apply((is.na(iState)),2,sum) #count the numbers of NA in objects per variable
  iVar<-ts(iState,start=1,end=Tn)
  assign(varNames[q],iVar)
}


#track of variables
bestVars<-data.frame(rep(1,nVars),row.names=varNames)
colnames(bestVars)<-"selected"

# hand cleaning  ----------------------------------------------------------

bestVars[,]<-0
bestVars[c(4,5,6,7,9,10,11,13),1]<-1

# segment per Class ---------------------------------------------------------------

objectClass<-as.matrix(as.data.frame(allFile$iAttr[[2]])+1)
objectClass[is.nan(objectClass)==1]<-max(objectClass,na.rm=T)+1
objectString<-as.character(unlist(allFile$iAttr[[4]]))
objectString[1]<-"Undet-Misc"
objectString[11]<-"Unident-Vehicle"
objectString[max(objectClass)]<-"NA"
#cumObjectClass<-as.matrix(as.data.frame(allFile$iAttr[[5]]))
statObjects<-matrix(0,nrow=length(objectString),ncol=1)
rownames(statObjects)<-objectString
for (j in 1:max(objectClass)){
  iClass<-apply(objectClass==j,2,as.numeric)
  statObjects[j,1]<-sum(colSums(iClass,na.rm = TRUE))
}

par(mfrow=c(1,1))
A<-barplot(t(statObjects), main="Objects distribution",horiz=F, 
           xlab="",las=2,cex.names=.6,cex.axis =.7,ylim=c(0,16000))
text(A,statObjects+400,statObjects,cex=.7)


# CAR data frame --------------------------------------------------

CarClass<-apply(objectClass==2,2,as.numeric)
CarClass[is.na(CarClass)==1]<-0
carTimes<-sum(colSums(CarClass))

carFrame<-as.data.frame(matrix(0,nrow=carTimes,ncol=(sum(bestVars)+2)))
varNames2<-c(varNames[bestVars==1],"numberCars","timePoint")
colnames(carFrame)<-varNames2
irow<-1
for (Ti in 1:Tn){
  Ttimes<-sum(CarClass[Ti,]) #how many cars are in time Ti
  j<-1 
  while (j<=Ttimes){
    for (q in 1:sum(bestVars)){
      carFrame[irow,q]<-get(varNames[bestVars==1][q])[Ti,which(CarClass[Ti,]==1)[j]]
    }
    carFrame[irow,(sum(bestVars)+2)]<-Ti
    carFrame[irow,(sum(bestVars)+1)]<-Ttimes
    j<-j+1
    irow<-irow+1
  }
}

timeinfo<-ncol(carFrame)
numberCars<-ncol(carFrame)-1
normDat<-apply(carFrame,2,normalizef)
par(mfrow=c(1,1))
boxplot(normDat[,-c(timeinfo,numberCars)],main="Car Object",las=2,cex.axis=.7)
summary(carFrame)


#data transformation
carFrame2<-carFrame 
carFrame2[,"speed"]<-log(carFrame2[,"speed"])
colnames(carFrame2)[5]<-"logSpeed"
carFrame2[,"acceleration"]<-log(carFrame2[,"acceleration"])
colnames(carFrame2)[8]<-"logAcc"

normDat2<-apply(carFrame2,2,normalizef)
par(mfrow=c(1,1))
boxplot(normDat2[,-c(timeinfo,numberCars)],main="Car Object",las=2,cex.axis=.7)
summary(carFrame)

# 

par(mfrow=c(3,3))
j<-1
for (q in 1:(ncol(carFrame2)-1)){
  hist(carFrame2[,q],main=colnames(carFrame2)[q])
  j<-j+1
}

raw_data<-normDat2[,-timeinfo]

# finding segments (corr -hc) --------------------------------------------------------


#dissimilarity
dist_cor <- 1 - cor(t(raw_data))
distancehc <- as.dist(dist_cor)

#number of clusters (takes time)
fclust<-fviz_nbclust(dist_cor, hcut, method = "wss",nboot=100,
                     hc_method = "complete")
#plot(fclust)
nclust<-knee2Line(fclust)
print(nclust)
hc_cor<-hclust(distancehc)
labelhc<-cutree(hc_cor,nclust)

par(mfrow=c(1,1))
plot(labelhc)
plot(carFrame$timePoint,labelhc)

newlabelhc<-removeNoise(10,labelhc,nclust)

par(mfrow=c(1,1))
plot(labelhc,type="p")
lines(newlabelhc,type="p",col="blue")

plot(carFrame$timePoint,labelhc)
lines(carFrame$timePoint,newlabelhc,type="p",col="blue")

segmentTime<-findSegments(carFrame$timePoint,newlabelhc)

cuts<-segmentTime$cutTime
time0<-segmentTime$time0
plot(segmentTime$lapses[,1],segmentTime$lapses[,2],type="p",col="blue",xlab="Time points",ylab="cluster")

#segmentation image
par(mfrow=c(1,1))
plot(carFrame2$timePoint,raw_data[,1],col="gray",type="p",cex=.5)
for (l in 2:(8)){
  lines(carFrame2$timePoint,raw_data[,l],col="gray",type="p",cex=.5)
}
sep1<-axis(1,at=cuts,tck = 1, lty = 2, col = "blue", labels = NA)


for (j in 1:length(time0)){
  rect(cuts[time0[j]]+5, 0, cuts[time0[j]+1]-5,1,col=rgb(1,0,0,alpha=0.2), border=F)}

cutscor<-cuts
time0cor<-time0
# save(cutscor,time0cor, file = "lapses_cor.RData")

# finding segments (abs corr -hc) --------------------------------------------------------

#dissimilarity
dist_cor <- 1 - abs(cor(t(raw_data)))
distancehc <- as.dist(dist_cor)


#number of clusters (takes time)
fclust<-fviz_nbclust(dist_cor, hcut, method = "wss",nboot=100,
                     hc_method = "complete")
plot(fclust)
nclust<-knee2Line(fclust)
print(nclust)
hc_cor2<-hclust(distancehc)
labelhc<-cutree(hc_cor2,nclust)

par(mfrow=c(1,1))
plot(labelhc)
plot(carFrame$timePoint,labelhc)

newlabelhc<-removeNoise(10,labelhc,nclust)

par(mfrow=c(1,1))
plot(labelhc,type="p")
lines(newlabelhc,type="p",col="blue")
plot(carFrame$timePoint,newlabelhc)

segmentTime<-findSegments(carFrame$timePoint,newlabelhc)

time0<-segmentTime$time0
cuts<-segmentTime$cutTime
#plot(segmentTime$lapses[,1],segmentTime$lapses[,2],type="p",col="blue",xlab="Time points",ylab="cluster")

#segmentation image
par(mfrow=c(1,1))
plot(carFrame2$timePoint,raw_data[,1],col="gray",type="p",cex=.5)
for (l in 2:(8)){
  lines(carFrame2$timePoint,raw_data[,l],col="gray",type="p",cex=.5)
}
sep1<-axis(1,at=cuts,tck = 1, lty = 2, col = "blue", labels = NA)


for (j in 1:length(time0)){
  rect(cuts[time0[j]]+5, 0, cuts[time0[j]+1]-5,1,col=rgb(1,0,0,alpha=0.2), border=F)}


cutscor2<-cuts
time0cor2<-time0
save(cutscor2,time0cor2, file = "lapses_cor2.RData")

# finding segments (corr^2 -hc) --------------------------------------------------------

#dissimilarity
dist_cor <- sqrt(1 - cor(t(raw_data)^2))
distancehc <- as.dist(dist_cor)

#number of clusters (takes time)
fclust<-fviz_nbclust(dist_cor, hcut, method = "wss",nboot=100,
                     hc_method = "complete")
plot(fclust)
nclust<-knee2Line(fclust)
print(nclust)

hc_cor3<-hclust(distancehc)
labelhc<-cutree(hc_cor3,nclust)

par(mfrow=c(1,1))
plot(labelhc)
plot(carFrame$timePoint,labelhc)

newlabelhc<-removeNoise(10,labelhc,nclust)

par(mfrow=c(1,1))
plot(labelhc,type="p")
lines(newlabelhc,type="p",col="blue")
plot(carFrame$timePoint,newlabelhc)

segmentTime<-findSegments(carFrame$timePoint,newlabelhc)

time0<-segmentTime$time0
cuts<-segmentTime$cutTime
#plot(segmentTime$lapses[,1],segmentTime$lapses[,2],type="p",col="blue",xlab="Time points",ylab="cluster")

#segmentation image
par(mfrow=c(1,1))
plot(carFrame2$timePoint,raw_data[,1],col="gray",type="p",cex=.5)
for (l in 2:(8)){
  lines(carFrame2$timePoint,raw_data[,l],col="gray",type="p",cex=.5)
}
sep1<-axis(1,at=cuts,tck = 1, lty = 2, col = "blue", labels = NA)


for (j in 1:length(time0)){
  rect(cuts[time0[j]]+5, 0, cuts[time0[j]+1]-5,1,col=rgb(1,0,0,alpha=0.2), border=F)}


cutscor3<-cuts
time0cor3<-time0
save(cutscor3,time0cor3, file = "lapses_cor3.RData")


# finding segments (dtw-hc) --------------------------------------------------------

distdtw<-dist(raw_data,method="dtw")

#number of clusters (takes time)
# fclust<-fviz_nbclust(raw_data, hcut,dist(raw_data,method="dtw"), method = "silhouette",nboot=100,
#                      hc_method = "complete")
# 
fclust<-fviz_nbclust(raw_data, hcut,dist(raw_data,method="dtw"), method = "gap_stat",
                     hc_method = "complete")
nclust<-which(diff(fclust$data[,"gap"])<0)[1]
plot(fclust)
#nclust<-knee2Line(fclust)



hc_dtw<-hclust(distdtw)
labelhc<-cutree(hc_dtw,nclust)

par(mfrow=c(1,1))
plot(labelhc)
plot(carFrame$timePoint,labelhc)

newlabelhc<-removeNoise(10,labelhc,nclust)

par(mfrow=c(1,1))
plot(labelhc,type="p")
lines(newlabelhc,type="p",col="blue")

segmentTime<-findSegments(carFrame$timePoint,newlabelhc)

cuts<-segmentTime$cutTime
time0<-segmentTime$time0
#plot(segmentTime$lapses[,1],segmentTime$lapses[,2],type="p",col="blue",xlab="Time points",ylab="cluster")

#segmentation image
par(mfrow=c(1,1))
plot(carFrame2$timePoint,raw_data[,1],col="gray",type="p",cex=.5)
for (l in 2:(8)){
  lines(carFrame2$timePoint,raw_data[,l],col="gray",type="p",cex=.5)
}
sep1<-axis(1,at=cuts,tck = 1, lty = 2, col = "blue", labels = NA)

for (j in 1:length(time0)){
  rect(cuts[time0[j]]+5, 0, cuts[time0[j]+1]-5,1,col=rgb(1,0,0,alpha=0.2), border=F)}


cutsdtw<-cuts
time0dtw<-time0
save(cutsdtw,time0dtw, file = "lapses_dtw.RData")

# finding segments (RF 2-hc) --------------------------------------------------------

#dissimilarity
dataR<-raw_data
dataS<-dataR
for (zz in (1:ncol(dataR))) {
  dataS[,zz]<-sample(dataR[,zz],nrow(dataR))}

dataall<-rbind(dataR,dataS)
dataall<-cbind(dataall,c(rep(1,nrow(dataR)),rep(0,nrow(dataR)))) ### synth or original 
colnames(dataall)[ncol(dataall)]<-"type"
# svdData<-svd(dataall)
# plot3d(dataall[,1:3],col=(1+dataall[,ncol(dataall)]), pch=1)
# plot3d(svdData$u[,1:3],col=(1+dataall[,ncol(dataall)]), pch=1)

clust.rf<-randomForest(as.factor(type)~.,data=dataall,proximity=T)

dist_rf<-matrix(1-clust.rf$proximity[1:nrow(dataR),1:nrow(dataR)],nrow(dataR),nrow(dataR))

# fclust<-fviz_nbclust(dist_rf, hcut, method = "silhouette",nboot=80,
#                      hc_method = "complete")
# 
# nbclust<-which(fclust$data[2]==max(fclust$data[2]))  
# plot(fclust)

fclust<-fviz_nbclust(dist_rf, hcut, method = "wss",
                     hc_method = "complete")

nclust<-knee2Line(fclust)

nclust<-5
#  nclust<-which(diff(fclust$data[,"gap"])<0)[1]

#number of clusters (takes time)
#nclust<-number_clusters(dist_rf)

hc_rf<-hclust(as.dist(dist_rf))
labelhc<-cutree(hc_rf,nclust)

par(mfrow=c(1,1))
plot(labelhc)
plot(carFrame2$timePoint,labelhc)

newlabelhc<-removeNoise(10,labelhc,nclust)

par(mfrow=c(1,1))
plot(labelhc,type="p")
lines(newlabelhc,type="p",col="blue")

segmentTime<-findSegments(carFrame2$timePoint,newlabelhc)

cuts<-segmentTime$cutTime
time0<-segmentTime$time0
#plot(segmentTime$lapses[,1],segmentTime$lapses[,2],type="p",col="blue",xlab="Time points",ylab="cluster")

#segmentation image
par(mfrow=c(1,1))
plot(carFrame2$timePoint,raw_data[,1],col="gray",type="p",cex=.5)
for (l in 2:(8)){
  lines(carFrame2$timePoint,raw_data[,l],col="gray",type="p",cex=.5)
}
sep1<-axis(1,at=cuts,tck = 1, lty = 2, col = "blue", labels = NA)


for (j in 1:length(time0)){
  rect(cuts[time0[j]]+5, 0, cuts[time0[j]+1]-5,1,col=rgb(1,0,0,alpha=0.2), border=F)}


cutsrf<-cuts
time0rf<-time0
save(cutsrf,time0rf, file = "lapses_rf.RData")


# finding segments (dtw -hc) --------------------------------------------------------

distdtw2<-dist(raw_data[,-timeinfo],method="dtw")

#finding number of clusters 
fclust<-fviz_nbclust(raw_data[,-timeinfo], hcut,dist(raw_data[,-timeinfo],method="dtw"), method = "gap_stat",
                     hc_method = "complete")
nclust<-which(diff(fclust$data[,"gap"])<0)[1]
# # # # # # # # # # # # #
 

hc<-hclust(distdtw2)
nclust<-3
labelhc<-cutree(hc,nclust)
plot(labelhc)
plot(carFrame$timePoint,labelhc)
# more1<-table(carFrame$timePoint)>1
# more1<-as.numeric(names(more1)[more1==T])
# lines(more1,rep(1.5,length(more1)),type="p",col="blue",pch=8)

th<-10
newlabelhc<-matrix(0,nrow=length(labelhc),ncol=1)
for (j in 1:nclust){
  icluster<-as.numeric(labelhc==j)
  #icluster<-smooth(icluster)
  changeLabel<-cumsum(abs(diff(icluster)))
  A<-table(changeLabel)>10
  B<-as.numeric(names(A)[A==T])  
  if (length(B)>1){
    C<-match(B,changeLabel)
    #lines(C,rep(.5,length(C)),type="p")
    new<-rep(0,length(icluster))
    
    for (l in 1:(length(C)-1)){
      new[C[l]:C[l+1]]<-getmode2(icluster[C[l]:C[l+1]])
    }
    newlabelhc[new==1]<-j
  }
  
}

par(mfrow=c(1,1))
plot(labelhc)
lines(newlabelhc,type="p",col="blue")

Q<- approx(carFrame$timePoint,newlabelhc,n=ceiling(length(icluster)/th),method="constant")
plot(Q$x,Q$y,type="s")


plot(carFrame$timePoint,newlabelhc)

par(mfrow=c(1,1))
plot(carFrame2$timePoint,raw_data[,1],col="gray",type="p",cex=.5)
for (l in 2:(8)){
  lines(carFrame2$timePoint,raw_data[,l],col="gray",type="p",cex=.5)
}

grides<-abs(diff(Q$y))
sep1<-axis(1,at=Q$x[grides==1],tck = 1, lty = 2, col = "blue", labels = NA)


#ts.plot(ts(X),ts(varPatts[,q]),gpars=list(col=c(rep("gray",sum(availObjects[q,])),2),lty=c(rep(2,sum(availObjects[q,])),1)),main=varNames[q])



# finding segments (corr2 -hc) --------------------------------------------------------

dissimilarity <- sqrt(1 - cor(t(carFrame2[,-timeinfo])^2))
distance <- as.dist(dissimilarity)

hc2<-hclust(distance)
par(mfrow=c(1,1))

#finding number of clusters 
fclust<-fviz_nbclust(raw_data[,-timeinfo], hcut,dist(raw_data[,-timeinfo],method="dtw"), method = "gap_stat",
                     hc_method = "complete")
nclust<-which(diff(fclust$data[,"gap"])<0)[1]


labelhc<-cutree(hc2,nclust)
plot(labelhc)
plot(carFrame$timePoint,labelhc)
# more1<-table(carFrame$timePoint)>1
# more1<-as.numeric(names(more1)[more1==T])
# lines(more1,rep(1.5,length(more1)),type="p",col="blue",pch=8)

th<-10
newlabelhc<-matrix(0,nrow=length(labelhc),ncol=1)
for (j in 1:nclust){
  icluster<-as.numeric(labelhc==j)
  #icluster<-smooth(icluster)
  changeLabel<-cumsum(abs(diff(icluster)))
  A<-table(changeLabel)>10
  B<-as.numeric(names(A)[A==T])  
  if (length(B)>1){
    C<-match(B,changeLabel)
    #lines(C,rep(.5,length(C)),type="p")
    new<-rep(0,length(icluster))
    
    for (l in 1:(length(C)-1)){
      new[C[l]:C[l+1]]<-getmode2(icluster[C[l]:C[l+1]])
    }
    newlabelhc[new==1]<-j
  }
  
}

par(mfrow=c(1,1))
lines(newlabelhc,type="p",col="blue")
plot(carFrame2$timePoint,newlabelhc)

Q<- approx(carFrame$timePoint,newlabelhc,n=ceiling(length(icluster)/th),method="constant")
plot(Q$x,Q$y,type="p",col="blue")


plot(carFrame$timePoint,newlabelhc)

par(mfrow=c(1,1))
plot(carFrame2$timePoint,raw_data[,1],col="gray",type="p",cex=.5)
for (l in 2:(8)){
  lines(carFrame2$timePoint,raw_data[,l],col="gray",type="p",cex=.5)
}

grides<-abs(diff(Q$y))
sep1<-axis(1,at=Q$x[grides==1],tck = 1, lty = 2, col = "blue", labels = NA)








# finding segments (RF no -hc) --------------------------------------------------------


dataR<-apply(carFrame2[,-timeinfo],2,normalizef)
dataS<-dataR
for (zz in (1:ncol(dataR))) {
  dataS[,zz]<-sample(dataR[,zz],nrow(dataR))}

dataall<-rbind(dataR,dataS)
dataall<-cbind(dataall,c(rep(1,nrow(dataR)),rep(0,nrow(dataR)))) ### synth or original 
colnames(dataall)[ncol(dataall)]<-"type"
svdData<-svd(dataall)
plot3d(dataall[,1:3],col=(1+dataall[,ncol(dataall)]), pch=1)
plot3d(svdData$u[,1:3],col=(1+dataall[,ncol(dataall)]), pch=1)

clust.rf<-randomForest(as.factor(type)~.,data=dataall,proximity=T)

distrf<-matrix(1-clust.rf$proximity[1:nrow(dataR),1:nrow(dataR)],nrow(dataR),nrow(dataR))

#xx<-cmdscale(1-clust.rf$proximity[1:nrow(dataR),1:nrow(dataR)])
#plot(xx)
hc3<-hclust(as.dist(distrf))
A<-as.dist(distrf)
par(mfrow=c(1,1))

#finding number of clusters 
wss <- function(d) {
  sum(scale(d, scale = FALSE)^2)
}
wrap <- function(i, hc, x) {
  cl <- cutree(hc, i)
  spl <- split(x, cl)
  wss <- sum(sapply(spl, wss))
  wss
}
res <- sapply(seq.int(1, 20), wrap, h = hc3, x = dataR)
plot(res)

labelhc<-cutree(hc3,5)
plot(labelhc)
plot(carFrame$timePoint,labelhc)
# more1<-table(carFrame$timePoint)>1
# more1<-as.numeric(names(more1)[more1==T])
# lines(more1,rep(1.5,length(more1)),type="p",col="blue",pch=8)

th<-10
newlabelhc<-matrix(0,nrow=length(labelhc),ncol=1)
for (j in 1:nclust){
  icluster<-as.numeric(labelhc==j)
  icluster<-smooth(icluster)
  changeLabel<-cumsum(abs(diff(icluster)))
  A<-table(changeLabel)>10
  B<-as.numeric(names(A)[A==T])  
  if (length(B)>1){
    C<-match(B,changeLabel)
    #lines(C,rep(.5,length(C)),type="p")
    new<-rep(0,length(icluster))
    
    for (l in 1:(length(C)-1)){
      new[C[l]:C[l+1]]<-getmode2(icluster[C[l]:C[l+1]])
    }
    newlabelhc[new==1]<-j
  }
  
}

par(mfrow=c(1,1))
plot(labelhc)
lines(newlabelhc,type="p",col="blue")
plot(carFrame2$timePoint,newlabelhc)

Q<- approx(carFrame$timePoint,newlabelhc,n=ceiling(length(icluster)/th),method="constant")
plot(Q$x,Q$y,type="p",col="blue")


plot(carFrame$timePoint,newlabelhc)

par(mfrow=c(1,1))
plot(carFrame2$timePoint,raw_data[,1],col="gray",type="p",cex=.5)
for (l in 2:(8)){
  lines(carFrame2$timePoint,raw_data[,l],col="gray",type="p",cex=.5)
}

grides<-abs(diff(Q$y))
sep1<-axis(1,at=Q$x[grides==1],tck = 1, lty = 2, col = "blue", labels = NA)





# final plot --------------------------------------------------------------


par(mfrow=c(1,1))
plot(carFrame2$timePoint,raw_data[,1],col="gray",type="p",cex=.5,xlim=c(0,598),xlab="Time points",ylab="",main="Time series segmentation for cars",cex.main=1)
for (l in 2:(8)){
  lines(carFrame2$timePoint,raw_data[,l],col="gray",type="p",cex=.5)
}
sep1<-axis(1,at=cutscor-4,tck = 1, lty = 3, col = "blue", lwd=2,labels = NA)
sep2<-axis(1,at=cutscor2,tck =1, lty = 3, col = "green",  lwd=2,labels = NA)
sep3<-axis(1,at=cutsdtw-4,tck = 1, lty = 3, col = "red",  lwd=2,labels = NA)
sep4<-axis(1,at=cutsrf,tck = 1, lty = 3, col = "cyan",  lwd=2,labels = NA)
leg<-c("(1-cor)","abs(1-cor)","DTW","RF")
legend(-2,.2,leg,cex=.7,col=c("blue","green","red","cyan"),lty=rep(3,4),lwd=rep(1.5,4))


# final plot 2 --------------------------------------------------------------


par(mfrow=c(1,1))
plot(carFrame2$timePoint,raw_data[,1],col="gray",type="p",cex=.5,xlim=c(0,598),xlab="Time points",ylab="",main="Time series segmentation for cars",cex.main=1)
for (l in 2:(8)){
  lines(carFrame2$timePoint,raw_data[,l],col="gray",type="p",cex=.5)
}
sep1<-axis(1,at=cutscor-4,tck = 1, lty = 3, col = "blue", lwd=2,labels = NA)
sep2<-axis(1,at=cutscor2,tck =1, lty = 3, col = "green",  lwd=2,labels = NA)
sep3<-axis(1,at=cutsdtw-4,tck = 1, lty = 3, col = "red",  lwd=2,labels = NA)
sep4<-axis(1,at=cutsrf,tck = 1, lty = 3, col = "cyan",  lwd=2,labels = NA)
leg<-c("(1-cor)","abs(1-cor)","DTW","RF")
legend(-2,.2,leg,cex=.7,col=c("blue","green","red","cyan"),lty=rep(3,4),lwd=rep(1.5,4))



par(mfrow=c(2,2))

##correlation 
plot(carFrame2$timePoint,raw_data[,1],col="gray",type="p",cex.main=.8,main="Correlation dissimilarity",xlab="Time points",ylab="",cex=.5)
for (l in 2:(8)){
  lines(carFrame2$timePoint,raw_data[,l],col="gray",type="p",cex=.5)
}
sep1<-axis(1,at=cutscor,tck = 1, lty = 2, col = "blue", labels = NA)


time0<-c(1:(length(cutscor)-1))[-time0cor]
for (j in 1:length(time0)){
  rect(cutscor[time0[j]]+5, 0, cutscor[time0[j]+1]-5,1,col=rgb(0,1,0,alpha=0.2), border=F)}

##correlation 2
plot(carFrame2$timePoint,raw_data[,1],col="gray",type="p",cex.main=.8,main="abs(correlation) dissimilarity",xlab="Time points",ylab="",cex=.5)
for (l in 2:(8)){
  lines(carFrame2$timePoint,raw_data[,l],col="gray",type="p",cex=.5)
}
sep1<-axis(1,at=cutscor2,tck = 1, lty = 2, col = "blue", labels = NA)


time0<-c(1:(length(cutscor2)-1))[-time0cor2]
for (j in 1:length(time0)){
  rect(cutscor2[time0[j]]+5, 0, cutscor2[time0[j]+1]-5,1,col=rgb(0,1,0,alpha=0.2), border=F)}

##DTW 
plot(carFrame2$timePoint,raw_data[,1],col="gray",type="p",cex.main=.8,main="DTW dissimilarity",xlab="Time points",ylab="",cex=.5)
for (l in 2:(8)){
  lines(carFrame2$timePoint,raw_data[,l],col="gray",type="p",cex=.5)
}
sep1<-axis(1,at=cutsdtw,tck = 1, lty = 2, col = "blue", labels = NA)


time0<-c(1:(length(cutscor)-1))[-time0dtw]
for (j in 1:length(time0)){
  rect(cutsdtw[time0[j]]+5, 0, cutsdtw[time0[j]+1]-5,1,col=rgb(0,1,0,alpha=0.2), border=F)}


##DTW 
plot(carFrame2$timePoint,raw_data[,1],col="gray",type="p",cex.main=.8,main="RF dissimilarity",xlab="Time points",ylab="",cex=.5)
for (l in 2:(8)){
  lines(carFrame2$timePoint,raw_data[,l],col="gray",type="p",cex=.5)
}
sep1<-axis(1,at=cutsrf,tck = 1, lty = 2, col = "blue", labels = NA)


time0<-c(1:(length(cutsrf)-1))[-time0rf]
for (j in 1:length(time0)){
  rect(cutsrf[time0[j]]+5, 0, cutsrf[time0[j]+1]-5,1,col=rgb(0,1,0,alpha=0.2), border=F)}