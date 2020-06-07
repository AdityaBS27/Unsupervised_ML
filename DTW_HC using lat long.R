rm(list = ls())
#Data Clustering using

#specify your folder 
setwd("~/Volvo Project/Programs/Low-dim mat files")# Libraries -----------------------------------------------
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
library(dtwclust)
library(TSclust)
library(dtw)

#-----------------General data-------------------------------------------

number_Ref<-1
name2<-(4992+(number_Ref-1)*60)
name2<-toString(formatC(name2,width=5,format="d",flag="0"))
name3<-".mat"
Ref_name<-paste(c(name2,name3),collapse="")
allFile<-readMat(Ref_name)
varNames<-rownames(allFile$iState)

#-----------------------------------------------------------------------------
LongPos<-allFile$iState[[1]]
LatPos<-allFile$iState[[2]]
x<-LongPos
x1<-LatPos

#-----------------------Tobias (Remove nan)-------------------------------------
# >>>> Remove columns with more than 50% NaN:
nan = rep(NA,length(x[1,]))
y = matrix(data = NA, nrow = length(x[,1]), 1) # Data Cleaned from "50% NA" columns
for (i in 1:length(x[1,])){
  nan[i] = sum(!is.finite(x[,i]))  # number of  NaN elements in column i
  if (nan[i]/length(x[,1]) < 0.5){
    y = cbind(y,x[,i])
  }
}
x = y[,-1]

# >>>> Remove all Na:
for(l in 1:length(x[1, ])) {
  if (is.na(x[1,l])) {        # If the first element is Na
    NonNAindex <- which(!is.na(x[, l]))
    firstNonNA <- min(NonNAindex)
    x[1:(firstNonNA - 1), l] = x[firstNonNA, l]
  }
  for (p in 2:length(x[, 1])){
    if (is.na(x[p, l])) {
      x[p, l] = x[p - 1, l]
    }
  }
}

#nan for lat
nan = rep(NA,length(x1[1,]))
y1 = matrix(data = NA, nrow = length(x1[,1]), 1) # Data Cleaned from "50% NA" columns
for (i in 1:length(x1[1,])){
  nan[i] = sum(!is.finite(x[,i]))  # number of  NaN elements in column i
  if (nan[i]/length(x1[,1]) < 0.5){
    y1= cbind(y1,x1[,i])
  }
}
x1= y1[,-1]

# >>>> Remove all Na:
for(l in 1:length(x1[1, ])) {
  if (is.na(x1[1,l])) {        # If the first element is Na
    NonNAindex <- which(!is.na(x1[, l]))
    firstNonNA <- min(NonNAindex)
    x1[1:(firstNonNA - 1), l] = x1[firstNonNA, l]
  }
  for (p in 2:length(x1[, 1])){
    if (is.na(x1[p, l])) {
      x1[p, l] = x1[p - 1, l]
    }
  }
}

#----------------------------------------------------------------------
nLongPos<-x
nLatpos<-x1

#New Data frame
nLong<-as.data.frame(nLongPos)
nLat<-as.data.frame(nLatpos);
names<-seq(1,24,1)

colnames(nLong)<-names
colnames(nLat)<-names

Q1<-matrix(0,nrow=1,ncol = 1)
Q2<-matrix(0,nrow=1,ncol = 1)
Q3<-matrix(0,nrow=1,ncol = 1)
Q4<-matrix(0,nrow=1,ncol = 1) 


for(i in 1:length(nLong[1,]))
{
  count<-i
  if(nLong[1,i] >0 && nLat[1,i]>0){
    Q1<-cbind(Q1[],count)
  }else if(nLong[1,i] >0 && nLat[1,i]<0){
    Q2<-cbind(Q2,count)
  }else if(nLong[1,i] <0 && nLat[1,i]>0){
    Q3<-cbind(Q3,count)
  }else{
    Q4<-cbind(Q4,count)
  }
  
}

#neglected q3 and q3 due less data for 4992

Q1LongPos<-t(as.data.frame(nLong[, Q1[,-1]]))
Q1LatPos<-t(as.data.frame(nLat[, Q1[,-1]]))
Q1Postion<-ts(sqrt(Q1LongPos^2+Q1LatPos^2))
Q1Postion


Q2LongPos<-t(as.data.frame(nLong[, Q2[,-1]]))
Q2LatPos<-t(as.data.frame(nLat[, Q2[,-1]]))
Q2Postion<-ts(sqrt(Q2LongPos^2+Q2LatPos^2))
Q2Postion



#Quantile1 clusterin using sbd and DTW
hc_sbdQ1 <- tsclust(Q1Postion, type = "h", k = 4,preproc = zscore,seed = 555,
                    distance = "sbd", centroid = shape_extraction,control = hierarchical_control(method = "average"))

plot(hc_sbdQ1, xlab="Cluster members")
plot(hc_sbdQ1, type = "sc")
#plot(hc_sbdQ1, type = "series", clus = 2)
#plot(hc_sbdQ1, type = "centroids", clus = 2)
zQ1Postion<- zscore(Q1Postion)
hc_dtwQ1<- tsclust(zQ1Postion, k = 4L,
                   distance = "dtw_basic", centroid = "dba",
                   trace = TRUE, seed = 555,
                   norm = "L2",
                   args = tsclust_args(cent = list(trace = TRUE)))

hc_dtwQ1 <- tsclust(Q1Postion, type = "h", k = 4,
                    preproc = zscore, seed = 555,
                    distance = "dtw", centroid = shape_extraction,
                    control = hierarchical_control(method = "average"))
plot(hc_dtwQ1, type = "sc")

#Quantile2 clusterin using sbd and DTW
hc_sbdQ2 <- tsclust(Q1Postion, type = "h", k = 4,
                    preproc = zscore, seed = 555,
                    distance = "sbd", centroid = shape_extraction,
                    control = hierarchical_control(method = "average"))

plot(hc_sbdQ1, xlab="Cluster menbers")
plot(hc_sbdQ1, type = "sc")


zQ2Postion<- zscore(Q2Postion)
hc_dtwQ2<- tsclust(zQ2Postion, k = 4L,
                   distance = "dtw_basic", centroid = "dba",
                   trace = TRUE, seed = 555,
                   norm = "L2",
                   args = tsclust_args(cent = list(trace = TRUE)))

hc_dtwQ1 <- tsclust(Q1Postion, type = "h", k = 4,
                    preproc = zscore, seed = 555,
                    distance = "dtw", centroid = shape_extraction,
                    control = hierarchical_control(method = "average"))
plot(hc_dtwQ1, type = "sc")





#Clusttering between quantiles dtw
Quantdistance<- proxy::dist(Q1Postion, Q2Postion,method = "dtw_basic")
#4 clusters
hc_dtw <- tsclust(D1, type = "h", k = 4,
                  preproc = zscore, seed = 555,
                  centroid = shape_extraction,
                  control = hierarchical_control(method = "average"))

plot(hc_dtw, xlab="Cluster members")
plot(hc_dtw, type = "sc")
plot(hc_sbd, type = "series", clus = 2)
plot(hc_sbd, type = "centroids", clus = 3)


#-------------------------------------------------------------------


