
# ## settin the workspace  ------------------------------------------------

setwd("D:\\gis certificate\\R-programing\\gully suceptibility\\data")
getwd()

# ##loading the necessary packge  -----------------------------------------

library(sp)
library(rgdal)
library(raster)
library(randomForest)
library(RStoolbox)
library(rgeos)
library(maptools)
library(caret)
library(gdalUtils)
library(RaTlist)
library(gdal_translate)
library(GDAL)
library(MODIStsp)
options(max.print=1000000000)
install.packages("rgdal")
library(MODIStsp)
library(MODIS)
library(covid19.analytics)
library(raster)
install.packages(modis)
# accessing  imageruy data ------------------------------------------------

Band2_PCN<- raster("D:\\gis certificate\\R-programing\\gully suceptibility\\data\\Geo-tiff\\PCN\\T37PCN_20180226T074901_B02_10m1.tif")
Band2_PCN
Band3_PCN<- raster("D:\\gis certificate\\R-programing\\gully suceptibility\\data\\Geo-tiff\\PCN\\T37PCN_20180226T074901_B03_10m1.tif")
Band3_PCN
Band4_PCN<- raster("D:\\gis certificate\\R-programing\\gully suceptibility\\data\\Geo-tiff\\PCN\\T37PCN_20180226T074901_B04_10m1.tif")
Band4_PCN
Band8_PCN<- raster("D:\\gis certificate\\R-programing\\gully suceptibility\\data\\Geo-tiff\\PCN\\T37PCN_20180226T074901_B08_10m1.tif")
Band8_PCN
BAnd11_Pbn<-raster("D:\\gis certificate\\R-programing\\gully suceptibility\\data\\S2A_MSIL2A_20180226T074901_N0206_R135_T37PBN_20180226T110356.SAFE\\GRANULE\\L2A_T37PBN_A014003_20180226T075914\\IMG_DATA\\R20m\\T37PBN_20180226T074901_B11_20m.jp2")
band11_pcn<-raster("D:/gis certificate/R-programing/gully suceptibility/data/S2A_MSIL2A_20180226T074901_N0206_R135_T37PCN_20180226T110356.SAFE/GRANULE/L2A_T37PCN_A014003_20180226T075914/IMG_DATA/R20m/T37PCN_20180226T074901_B11_20m.jp2")
Sentinel_Raster_stack_PCN<-stack(Band2_PCN,Band3_PCN,Band4_PCN,Band8_PCN)
## moscid band 11of PBN and PCN
SWIR<-mosaic(BAnd11_Pbn,band11_pcn,fun=mean)
SWIR
## resamplingupscaling 20m cell size in to 10m 
SWIR<-raster::resample(SWIR,sentinel_mosiced,method="bilinear")
 
SWIR


## ploting true colure composite ( band 2, band3, and  band 4 , using the band order r=3,g=2, and b=1  by consiodering the order of the band at raster stack)
plotRGB(Sentinel_stack_TCC,r=3,g=2,b=1,scale=32767,stretch="lin")
## ploting False colur composite , by replacing Red by NIR,green by red, and blue by green; simply  using band 4, 3, and 2 order
plotRGB(Sentinel_Raster_stack_PCN,r=4,g=3,b=2,scale=32767,stretch="lin")
## wriing the file  in to the disk
writeRaster(Sentinel_Raster_stack_PCN,filename = "Raster_stack_scene_PCN",format="GTiff")

 ## the same procedure   would goes  for scene2-PBN

Band2_PBN<- raster("D:\\gis certificate\\R-programing\\gully suceptibility\\data\\Geo-tiff\\PBN\\T37PBN_20180226T074901_B02_10m1.tif")
Band3_PBN<- raster("D:\\gis certificate\\R-programing\\gully suceptibility\\data\\Geo-tiff\\PBN\\T37PBN_20180226T074901_B03_10m1.tif")
Band3_PBN
Band4_PBN<- raster("D:\\gis certificate\\R-programing\\gully suceptibility\\data\\Geo-tiff\\PBN\\T37PBN_20180226T074901_B04_10m1.tif")
Band4_PBN
Band8_PBN<- raster("D:\\gis certificate\\R-programing\\gully suceptibility\\data\\Geo-tiff\\PBN\\T37PBN_20180226T074901_B08_10m1.tif")
Band8_PBN
sentinel_raster_stack_PBN<-stack(Band2_PBN,Band3_PBN,Band4_PBN,Band8_PBN)

plotRGB(sentinel_raster_stack_PBN,r=3,g=2,b=1,scale=32767,stretch="lin")
plotRGB(sentinel_raster_stack_PBN,r=4,g=3,b=2,scale=65535,stretch="lin")
sentinel_raster_stack_PBN

writeRaster(sentinel_raster_stack_PBN,filename = "Raster_stack_scene_PBN",format="GTiff")
## mosiacing, so that calling the two overlaping images which needs to be mosiaced , and the result is  raster brick 
sentinel_mosiced<-mosaic(Sentinel_Raster_stack_PCN,sentinel_raster_stack_PBN,fun=mean)
sentinel_mosiced

plotRGB(mosiced,r=4,g=3,b=2,scale=32767,stretch="lin")
##write  the mosiaced image to the disc as follows 
writeRaster(sentinel_mosiced,filename = "sentinel2_mosiaced_koga",format="GTiff")

# now to minimize the processing time , it is preferebale to spati --------

E<-extent(275024,330428,1226618,1268688)
## now crop the mosiced Sentinel2-image by E using crop function to retun  geographic subset of an object
Koga_sentinel2_subseted<-crop(sentinel_mosiced,E)

## now let it be to the extenet of   the study area 
SWIR_subsetd <-crop(SWIR,E)## it has to be now  stacked with aother
SWIR_subsetd
## stacking band 11 top other 
Koga_sentinel2_subseted<-stack(Koga_sentinel2_subseted,SWIR_subsetd)

##  names the layer
names(Koga_sentinel2_subseted)<-c("Band2","Band3","Band4","Band8","SWIR_Band")
Koga_sentinel2_subseted
## ploting the substetd image usin the plotRGB fubction 
plotRGB(Koga_sentinel2_subseted,r=4,g=3,b=2,scale=32767,stretch="lin")
## write the subseted image to the disk
raster::writeRaster(Koga_sentinel2_subseted,filename = "sentinel2_koga_finale_substed",format="GTiff",overwrite=TRUE)

# ## loading sample points  -----------------------------------------------

options(max.print = 1000000000)
##step1: generting refernce/smaple data for training the classifier, in this case sample points has been collected from arial photographs)
## refernce  GCP points were collected from 5cm resolution aerial photograph
setwd("D:\\gis certificate\\R-programing\\gully suceptibility\\data")
training<-read.csv("reference sptialpoints.csv")
## randomisation of the trainig  sample 
training1<-training[ ,-1]
training1
Training_randomised<- training1[sample(1:nrow(training1)), ]
Training_randomised



## the cateogorical value should be converted to numeric  first,convert landusetype col to numeric 
Training_randomised$landusetype<-as.numeric(Training_randomised$landusetype)
Training_randomised$landusetype
## for calssification IN RF the response variable has to be converted to factor 
trainingset<-as.factor(Training_randomised$landusetype)
table(trainingset)
length(trainingset)

####drop the ID column

##converting the table in to dataframe 
training.df<-data.frame(Training_randomised)
training.df
##data partion ; make random seed to make the analysis repeatable

## assigning  coordinates  for the df
coordinates(training.df) <- c("lon", "lat")

## check wether preojection is defined or not 
is.projected(training.df)
##  defining cRS
proj4string(training.df) <- CRS("+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs")
## checking whether project is defined ( now it reply True)
is.projected(training.df)
training.df
plot(training.df)
## first all the layers should have the same dimension , using the AOI(shapefile)
##calling the  sentinel 1 image  to be considered as variable 
E1<-extent(285560,313180,1234680,1262490)
sentinel1_image<-raster("D:\\gis certificate\\R-programing\\gully suceptibility\\data\\sentinel1 image _lULC\\S1A_IW_GRDH_1SDV_20181230T031657_subseted_Cal_Spk_TC.tif")

## cropiing and mascking sentinel2&sentinel  to the extent of the images 
Koga_sentinel2_subseted1<-crop(Koga_sentinel2_subseted,extent(Koga_shape))
sentinel2_Masked<-mask(Koga_sentinel2_subseted1,Koga_shape)## mascking to the extenet of the study watershed 
sentinel2_Masked
names(sentinel2_Masked)<-c("Band2","Band3","Band4","Band8","SWIR_Band")
## project raster 
VVVHradar_cband<-projectRaster(sentinel1_image,sentinel2_Masked,res = 10)##  Sentinel1 image to be fused as apredictor 
plot(VVVHradar_cband)
VVVHradar_cband
# ## defining other  covariets /predictors ;exctracted from processed  sentinel2A-Level2A image--------------------------------
NDVI<-overlay(sentinel2_Masked$Band8,sentinel2_Masked$Band4,fun=function(x,y){(x-y)/(x+y)})## NDVI  
plot(NDVI)
NDWI<-overlay(sentinel2_Masked$Band3,sentinel2_Masked$Band8,fun=function(x,y){(x-y)/(x+y)})## NDWI 
plot(NDWI)
NDBI<-overlay(sentinel2_Masked$SWIR_Band,sentinel2_Masked$Band8,fun=function(x,y){(x-y)/(x+y)})## NDBI
plot(NDBI)## using SWIR  of sentinel image 

# ## adding covariet as one raster brick -------------------------------

covs <- addLayer(sentinel2_Masked, NDVI, NDWI,NDBI,VVVHradar_cband)
covs
names(covs)<-c("Band2","Band3","Band4","Band8","SWIR_Band","NDVI","NDWI","NDBI","VV/VH radar_cband")## giving explict name to each covariets/predictors 
## ploting the predictor/covarietes 
plot(covs)


## plot the point 
plot(training.df)
summary(training.df)
str(training.df)
length(training.df)

# ##extract the values of the predictors at the locations of the p --------

extractedvalue1<-as.data.frame(extract(covs,training.df))
extractedvalue1

training.df@data=data.frame(training.df@data,extractedvalue1[match(rownames(training.df@data),rownames(extractedvalue1)),])
training.df@data
library(RStoolbox)
resam
## partitioning train and test data  sets
train.idx <- sample(nrow(covs), 2/3 * nrow(covs))
covs.train <- covs[train.idx, ]
covs.test <- covs[-train.idx, ]

# ## model fiting & prediction  -------------------------------------------
## model selection 
library(caret)
set.seed(123)
finalModel <-randomForest(x = training.df@data[,2:ncol(training.df@data)],
                          y =as.factor(training.df@data$landusetype),proximity=TRUE,
    
                                                ntree = 2000,mtry=3,replace = TRUE, importance = TRUE)
plot(finalModel)
## ploting sample tree 

library(rpart)
cart <- rpart(as.factor(training.df@data$landusetype)~., data=training.df@data[,2:ncol(training.df@data)], method = 'class', minsplit = 9)
-
        


# ## Model Estimation  using the  fitted model  ---------------------------

predict<-predict(covs,model=finalModel,format="GTiff", progress='text')


preiction<-predict(covs, finalModel, "Koga_landuselandcover4",format="GTiff",
                       type="response", index=2,
                       na.rm=TRUE, overwrite=TRUE, progress="text")## landuse land cover predisction based on selected model

plot(pre2)


# ##modelfit --------------------------------------------------------------

rf.pred <- predict(finalModel, training.df@data[,2:ncol(training.df@data)], type="response",na.action=na.pass)
rf.pred
rf.prob <- as.data.frame(predict(finalModel, training.df@data[,2:ncol(training.df@data)], type="prob"))
rf.prob
obs.pred <- data.frame(cbind(Observed=as.numeric(as.character(training.df@data[,"landusetype"])),
                             PRED=as.numeric(as.character(rf.pred)), Prob1=rf.prob[,2],Prob2=rf.prob[,3],Prob3=rf.prob[,4],Prob4=rf.prob[,5],Prob5=rf.prob[,6],
                             Prob0=rf.prob[,1]) )
library(caret)## needed to estaimate confusion matrix 
op <- (obs.pred$Observed == obs.pred$PRED)
op 
tablepred<-table(obs.pred$Observed,obs.pred$PRED)## confusion matrix
tablepred
confusionMatrix(op)
confusionMatrix(tablepred)
( pcc <- (length(op[op == "TRUE"]) / length(op))*100 )
install.packages("verification")
library(verification)
library(rfUtilities)## needed to estimate acuracy 
library(OOBCurve)
accuracy
accuracy(obs.pred$Observed,obs.pred$PRED)## calssification accuracy
acc
obs.pred$Observed
str(obs.pred$PRED)


# labling and ploting RF product ------------------------------------------

rf.prob
## ploting 
plot(predict)
## 
plot(finalModel, main="Bootstrap Error Convergence")
p <- as.matrix(finalModel$importance)   
plot(p)

## calling the rf product-final classification 
Koag_landuselandcover<-raster("koga-landuselandcover.tif")
varImpPlot(finalModel)

plot(Koag_landuselandcover)

##  study area shape 
Koga_shape<-readOGR("d:\\gis certificate\\R-programing\\gully suceptibility\\data\\mosaic\\KOGA_SHAPE.shp")
## crop the result to the extenet of koga watershed
## transformin projection
koga_shape<-proj4string(Koga_shape) <- CRS("+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs")
koga_shape
## extract cell value to the extyenet of koga watershed
library(rgeos)
## Clipping raster using shapefile in R, but keeping the geometry of the shapefile
## first doing crooping using the extenet 

KOgaLUL<-crop(Koag_landuselandcover,extent(Koga_shape))##croping
koga_lulc_mask<-mask(KOgaLUL,Koga_shape)## landus eland cover to the exetent of koga watrehsed 

plot(koga_lulc_mask)


f<-freq(koga_lulc_mask, useNA='no')## frequesncy of cellvalue at each category
f
produt_resolution<-prod(res(koga_lulc_mask))
produt_resolution
landuselandcover_area<-cbind(f, area=f[,2] * produt_resolution)
landuselandcover_area

##the predicted land use land cover rastre layer has to be chnaged in to category for ploting and leveling 
lulc<-as.factor(koga_lulc_mask)
lulc
rat <- levels(lulc)[[1]]
rat
##assigning lnad cover land use class names for each value using the follwing code 
rat[["landcover"]] <- c("Built_up_residential", "cultivated", "Forest", "Grazing","Shrub","Water")
rat[["landcover"]]
rat
levels(lulc) <- rat
library(rasterVis)

levelplot(lulc, col.regions=c("pink","brown","darkgreen","cyan","green","blue"), xlab="", ylab="")


# ## clasiffication using ML classifiers ##  --------------


library(rasclass)
help("rasclass")
library(RStoolbox)
library(ggplot2)
library(lattice)
library(rasterVis)

## Maximum Likelihood classification using RStoolbox pack 
##fitting the model by splitting the data in to 70%(trainig set),and 30% test set

training.df.idx = sort(sample(nrow(training.df), nrow(training.df)*.7))

#creating training data set by selecting the output row values
train<-training.df[training.df.idx,]
test<-training.df[-training.df.idx,]

#creating test data set by not selecting the output row values
test<-data[-data1,]
MLC<-superClass(Koga_sentinel2_subseted,trainData=train,responseCol="landusetype",minDist = 1,
                model = "mlc",valData=test, kfold = 5,format="GTiff",mode = "classification",
                predict = TRUE,predType="raw",filename = "mximum_classifiaction")


raster::writeRaster(MLC,filename="maximmululc",format="GTiff")
MLC<-ratify(MLC)
## validating the result 
getValidation(MLC, from = "testset", metrics = "overall")

## plot the result 
MLC

levelplot(MLC, col.regions=c("pink","brown","darkgreen","cyan","green","blue"), xlab="", ylab="")
ggR(MLC,stretch="lin")

# calling and leveling and ploting maximum likeli hood lulc product -------

max<-raster("mximum_classifiaction.grd")

raster::writeRaster(max,filename="maximmululc",format="GTiff")

MAXLULU.tiff<-raster("maximmululc.tif")

plot(MAXLULU.tiff)

## crop to the extent of the study area 
kogalulcMAx<-crop(MAXLULU.tiff,extent(Koga_shape))## croping 
KogaLulc_max_mask<-mask(kogalulcMAx,Koga_shape)## masking to the exetenet of kogawatershed
## ploting and checking the result of masking
plot(KogaLulc_max_mask)

# ## statistical analysis and colouring  ----------------------------------
MaxLulc<-as.factor(KogaLulc_max_mask)
rat<-levels(MaxLulc)[[1]]
rat

####assigning lnad cover land use class names for each value using the follwing code 
rat[["landcovertype"]]<-c("Built_up_residential", "cultivated", "Forest", "Grazing","Shrub","Water")

rat[["landcovertype"]]
levels(MaxLulc)<-rat## attyaching the cell value to each land cover types 

## now is better to colouring each category as what is done in  RF classifiers using library rasterVis 
library(rasterVis)

levelplot(MaxLulc, col.regions=c("pink","brown","darkgreen","cyan","green","blue"), xlab="", ylab="")
## claculating area representation of each class using maximuym likelihood classifir

f<-freq(MaxLulc, useNA='no')## frequesncy of cellvalue at each category
f
produt_resolution<-prod(res(MaxLulc))
produt_resolution
lulc_area_max<-cbind(f, area=f[,2] * produt_resolution)

lulc_area_max
## adding column percentage 
rat[["percntcoverage"]]<-(lulc_area_max[[area]]/10000 * (1/28152))


## creating confusion matrixs as follows 
rf.pred1 <- predict(MLC$model,train[,2:ncol(train)],type="raw",na.action=na.pass)
obs.pred1 <- data.frame(cbind(Observed=as.numeric(as.character(train@data[,"landusetype"])),
                             PRED=as.numeric(as.character(rf.pred))))
library(caret)
                             
tablepred1<-table(obs.pred1$Observed,obs.pred1$PRED)## confusion matrix
tablepred1          
## then following the confusioon matrix, estimate the user and producer accuracy as follows
# User accuracy
diag(tablepred1) /rowSums(tablepred1)## user accuracy 

diag(tablepred1) /colSums(tablepred1)## producer accuracy 

















lulc2<-extract(KOgaLUL,Koga_shape)
plot(lulc2)
xx<-ranger(as.factor(landusetype)~.,data=training.df@data,keep.inbag = TRUE,mtry=2,num.trees = 501,replace = FALSE,sample.fraction = 0.632,save.memory = TRUE,verbose = TRUE,oob.error = TRUE,min.node.size =3)
xx



## roc and AUC
library(pROC)
pRoc

training.df@data
training.df@data[,]

##  lableing and 

## model construction and selection 
library(randomForest)
library(caret)
library(e1071)
library(PerformanceAnalytics)
library(rfUtilities)
library(ranger)
## using ranger packages ; the fastest application of random forest
## modelfit
RFmodelbets<-randomForest()


train.idx <- sample(nrow(training.df@data), 2/3 * nrow(training.df@data))## subsampling with out replacemnet
length(train.idx) 
Koga.train.idx<-training.df@data[train.idx,]  ## training data

Koga.test.idx<-training.df@data[-train.idx,]  ## randomized test data
Koga.test.idx
set.seed(1234)
Model.fitted<-ranger(as.factor(landusetype)~.,data =Koga.train.idx,mtry=4,oob.error = TRUE,num.trees =100,keep.inbag = TRUE) ## model to be fitt

Model.fitted
predicted.koga<-predict(Model.fitted,data=Koga.test.idx,progress="text")## model prediction 

plot(predicted.koga)
##model selection using rfUtilities 
buffer.distance 

#  parametr tuning using tunrf  function ; using caret package 
myseed<-set.seed(432)
my_tune<-tuneRF(training.df@data,as.factor(training.df@data$landusetype), mtryStart=1, ntreeTry=500, stepFactor=2, improve=0.05,
                plot=TRUE, doBest=FALSE)
## class ballancing using Rf.utilities class balanace 
library(stackoverflow)
library(party)

## the lowr the OOB error the better the classification accuracy; oObB error gets zero at mtry=4
set.seed(1234)
modelse<-rf.modelSel(x=training.df@data[,2:ncol(training.df@data)], 
                     y=as.factor(training.df@data[,"landusetype"]),imp.scale = "mir",r=c(0.1, 0.2, 0.5, 0.7, 0.9),final.model = TRUE,seed = 2345,parsimony =0.5,ntree=1000,mtry=c(2,3,4,5,6,7,8,9,10),replace=TRUE,proximity=TRUE)
modelse
##We can then subset the parameters and run a fit model.



plot(predict)
finalModel$importance
variable<-modelse$parameters[[1]]
class(finalModel)
library(ranger)
## the otherwayround
set.seed(423535)
fitRFmanual <- randomForest (as.factor(landusetype) ~ .,data=training.df@data,
                             mtry =4,ntree=1000,importance=TRUE,replace = TRUE,proximity = TRUE,norm.votes = TRUE)
fitRFmanual

finalmodel3 <- randomForest(as.factor(training.df@data[,"landusetype"]),data=training.df@data,
                  importance=TRUE,ntree = 1000,replace = TRUE,proximity = TRUE,norm.votes = TRUE)



finalmodel3
op
## make prictionn using s3 class'randomForest" as follows 
training
training.df@data

koga_sentinel_stack1<-paste(rownames(finalModel$importance))  


#testing the multi-coliniarity of the features 

cl <- multi.collinear(x=training.df@data[,2:ncol(training.df@data)], p=0.05)
plot(cl)
overlay<-over(training.df@data@data["landusetype"],koga_sentinel_stack)

X<-predict(finalmodel3,training.df@data,filename="koga-landuselandcover2",type="response",
      index=1, na.rm=TRUE,format="GTiff", progress='text')
m<-raster("koga-landuselandcover.tif")
plot(m)
predicted2<-predict(training.df@data,finalModel,type="response",na.rm=TRUE,index=1,format="GTiff",na.action = na.omit,progress="text")
plot(predicted2)
Koga_sentinel2_subseted
newdata1<-as.matrix(Koga_sentinel2_subseted)
newdata1
makelea
plotObsVsPred()



colnames(training.df@data)
colnames(as.matrix(Koga_sentinel2_subseted))
Koga_sentinel2_subseted
colnames(finalModel$importance)
plot(finalModel$err.rate[,1])
## plot OOB CUrve
library(OOBCurve)
library(mlr)
library(ranger)
library(sp)
library(rgdal)
library(raster)
library(rasterVis)
predic

##  final prediction 
set.seed(myseed)
predicted2<-predict(Koga_sentinel2_subseted,finalModel,type="response",na.rm=TRUE)
plot(predicted2)

Koga_sentinel2_subseted
training.df@data[,vars]



#Let's try the build the model with the default values.

  




model<-train(as.factor(landusetype)~.,training.df@data,method="rf",TuneLength=3,trcontrol=trainControl(method = "00b",number = 10,classProbs = TRUE))

model$results
model$bestTune
model$modelInfo
set.seed(myseed)
fitRFmanual <- randomForest (as.factor(landusetype) ~ .,data=training.df@data,
                            mtry =1,ntree=500)
fitRFmanual
fitRFmanual$votes
mySeeds <- sapply(simplify = FALSE, 1:26, function(u) sample(10^4, 3))
cvCtrl = trainControl(method = "repeatedcv", number = 5, repeats = 5,
                      classProbs = TRUE, summaryFunction = sixClassSummary,
                      seeds = mySeeds)
cvctrl=trai

rf.mdl<-train(as.factor(landusetype) ~ .,method="rf",mtry=2,ntree=100,data=training.df@data)
rf.mdl

# data set of model predictions on training data vs. actual observations
results <- data.frame(pred = predict(rf.mdl, training.df@data),
                      obs = as.factor(training.df@data$landusetype))
table(results)


getTree(rf.mdl)
predicted1<-predict(rf.mdl,training.df@data)
predicted1
set.seed(256)
responsedat<-as.factor(training.df@data$landusetype)
model1<-train(x=training.df@data,y=responsedat,method="rf")
model1
randomforestmodel<-randomForest(as.factor(landusetype) ~ .,mtry=2,ntree=1000,data = training.df@data)



confusionMatrix(predicted1,as.factor(training.df@data$landusetype))

print(rf_mtry)




summary(results_tree)

predicted<-predict(defualt_rf,training.df@data)








## partitioning the data for training and testing 
set.seed(123456)
ind<-sample(2,nrow(training.df@data),replace=TRUE,prob=c(0.8,0.2))
train<-training.df@data[ind==1, ]
train
test<-training.df@data[ind==2, ]
test
training.df
## constract random forest model ;## model fit
set.seed(1267)
randomforestmodel2<-randomForest(as.factor(landusetype) ~ .,mtry=3,ntree=50,data = training.df@data,varImp=TRUE)
randomforestmodel2
plot(randomforestmodel2)
table(randomforestmodel2)
library(OOBCurve)

install.packages(mlr)
remove.packages(glue)
devtools::install_github("PhilippPro/OOBCurve")
library(OOBCurve)
## to built tree in random forest use library "party" as follows 
library(caret)
plot.RO
library(party)
library(rpart)
library(data.table)
library(mlr)
library(h2o)
library(randomForest)
library(raster)
predi
## checking randomforeswt algorithms with defualt value 
control<-trainControl(method = "repeatedcv",number = 10,repeats = 3,search = "random")
seed<-7
metric<-"Accuracy"
set.seed(seed)
mtry<-sqrt(ncol(training.df))
tuneGrid<-expand.grid(.mtry=mtry)
rf_defualt<-train(as.factor(landusetype) ~ .,data=training.df@data,method="rf",metric=metric,tuneGrid=tuneGrid,ntree=500)
rf_defualt


set.seed(23465)
mtry_optimization<-tuneRF(training.df@data,as.factor(training.df@data$landusetype),stepFactor = 1.5,improve = 0.5,trace = TRUE,plot = TRUE,missing=TRUE)

mtrytunin<- tuneRF(training.df@data,landusetype,stepFactor = 1.5,improve = 0.5,trace = TRUE,plot = TRUE)

mtrytunin
cart <- rpart(as.factor(landusetype)~., data=training.df@data, method = 'class', minsplit = 3)
##Much cleaner way is to plot the trained classification tree
plot(cart, uniform=TRUE, main="Classification Tree")
text(cart, cex = 0.6)
table(cart)

predicted2<-as.data.frame(predict(Koga_sentinel2_subseted,finalModel,type="response",na.rm=TRUE))
plot(predicted2)

f<-freq(predicted2, useNA='no')## frequesncy of cellvalue at each category
f
produt_resolution<-prod(res(predicted2))
produt_resolution
landuselandcover_area<-cbind(f, area=f[,2] * produt_resolution)
landuselandcover_area
plot(koga_lulc)
str(predicted2)  
##the predicted land use land cover rastre layer has to be chnaged in to category for ploting and leveling 
lulc<-as.factor(koga_lulc)
lulc
rat <- levels(lulc)[[1]]
rat
##assigning lnad cover land use class names for each value using the follwing code 
rat[["landcover"]] <- c("Built_up_residential", "cultivated", "Forest", "Grazing","Shrub","Water")
rat[["landcover"]]
rat
levels(lulc) <- rat
library(rasterVis)

levelplot(lulc, col.regions=c("pink","brown","green","yellow","lightgreen","blue"), xlab="", ylab="")
RFmarkerDetector::tuneNTREE(as.factor(training.df@data$landusetype),2,20,minNTREE = 500,length=1)

x<-readOGR(dsn="mosaic",layer="KOGA_SHAPE")
x

confusionMatrix(predicted2,as.factor(training.df@data$landusetype))

library(cluster)
library(snow)
                   
                   names(getModelInfo())
                   defualt_rf
                   defualt_rf
                   defualt_rf
                   
                   head(modelfit.rf)
                   rf.mdl<-randomForest(as.factor(landusetype) ~ .,training.df@data)
                   rf.mdl$err.rate
                   print(rf.mdl)
                   names(rf.mdl)
                   plot(rf.mdl$forest)         
                   
                   beginCluster()
                   
                   preds_rf <- cluster(Koga_sentinel2_subseted,raster::predict,args= list (model = finalModel ))
                   
                   endCluster()
                   
                   plot(preds_rf)
                   (preds_rf)
                   
                  
               
               
               library(caret)   
                     p1<-predict( modelfit.rf,train)
                     confusionMatrix(p1,as.factor(train$landusetype))
                     getTree(p1,k=1,labelVar = FALSE)
                
                     head(p1)
                     predict(Koga_sentinel2_subseted, finalmodel3, filename="RfClassPred.img", type="response",
                             index=1, na.rm=TRUE, progress="window", overwrite=TRUE)

                     

# PART2: WITH CLUSTER TRAINING DATA ---------------------------------------

                     
                     

# RF_Supervised Machine learning with  cluster training points  -----------

sample3<-readOGR("D:\\gis certificate\\R-programing\\gully suceptibility\\data\\data\\Training_clusters.shp") 
                  
                    sample3
         
##  assign the same projection as  the image              
                    
  sample_final1<- spTransform(fores,CRS(proj4string( Koga_sentinel2_subseted)))              
                    

 sample_final1

 
                                                      

## the cateogorical value should be converted to numeric  first,convert landusetype col to numeric
 sample_final1$landusetyp<-as.numeric(sample_final1@data$landusetyp)
   table( sample_final1$landusetyp)
## for calssification IN RF the response variable has to be converted to factor 
   sample_final1$landusetyp<-as.factor(sample_final1@data$landusetyp)


# ##extract the values of the predictors at the locations of the p --------
extractedvalue11<-as.data.frame(extract(covs, sample_final1,df=TRUE))
extractedvalue11
 
ab1@data=data.frame(ab1@data,extractedvalue9[match(rownames(ab1@data),rownames(extractedvalue9)),])## row and column maching b/n extracted value and dataframe 
ab1@data


   sample22.df<-ab1@data[ -2]
   sample22.df1<-na.omit(sample22.df)
   class(sample22.df)
   sample22.df1
  
table(sample2_polygon$landusetyp)
 
  
  

 
   Impsample<-sample22.df[ -2]
   Impsample
   selectedVar<- Impsample[ -2]
   selectedVar
   # ## model fiting & prediction  -------------------------------------------

   ## model selection 
   library(caret)
   ## accuracy at defualt value of ntree&mtry with all involved variables 
set.seed(123)
RFModel1 <-randomForest(x=  training.df[,2:ncol(training.df)],
                             y = as.factor(   training.df$landusetyp),proximity=TRUE,
                             
                             ntree = 500,mtry=3,replace = TRUE, importance = TRUE)
   
   RFModel1
   
   ## prediction @scenario1
   rf.pred1<-predict(covs,RFModel1,format="GTiff", progress='text',type="response",na.rm=TRUE)
   rf.pred1
  
   ## accuarcy assesesmntye for model1/ modelfit
   rf.pred2 <- predict(RFModel1,  sample22.df[,2:ncol(sample22.df)], type="response")
 
   obs.pred1 <- data.frame(cbind(Observed=as.numeric(as.factor(sample22.df[,"landusetyp"]))),
                                PRED=as.numeric(as.character(rf.pred2)))
   op <- (obs.pred$Observed == obs.pred$PRED)
   
   library(caret)## needed to estaimate confusion matrix 
   
   tablepred1<-table(obs.pred1$Observed,obs.pred1$PRED)## confusion matrix
   tablepred1
   confusionMatrix(op)
   confusionMatrix(tablepred1)
   
   library(verification)
   library(rfUtilities)## needed to estimate acuracy 
   library(OOBCurve)
   accuracy
   accuracy(obs.pred1$Observed,obs.pred1$PRED)## calssification accuracy
   acc
   obs.pred$Observed
   str(obs.pred$PRED)
                          
   ##assigning lnad cover land use class names for each value using the follwing code 
   rat[["landcover1"]] <- c("Built_up_residential", "cultivated", "Forest", "Grazing","Shrub","Water")
   rat[["landcover"]]
   rat
   levels(lulc) <- rat
   library(rasterVis)
   
   levelplot(lulc, col.regions=c("pink","brown","darkgreen","cyan","green","blue"), xlab="", ylab="")
   
   
   
   
         
                                rf.prob <- as.data.frame(predict(rf.fit, sdata@data[,sel.vars], type="prob"))
   ##scenario2,ntree=600,mtry=3,tuned parameter
   
set.seed(123)
RFModel2 <-randomForest(x= training.df[,2:ncol(training.df)],
                          y = as.factor(training.df$landusetyp),proximity=TRUE,
                          
                          ntree = 700,mtry=3,replace = TRUE, importance = TRUE)
   
   
   RFModel2
  
   
   ##prediction at scenario2
   
   predict2<-predict(covs,model=  RFModel2,format="GTiff", filename="Rf_koga_lULC.tif",progress='text')
   plot(predict2)
   
rf.pred<- predict(covs, RFModel2, type="response",norm.votes=TRUE,format="GTiff",
         filename="Rf_lULC.tif3",progress='text',index=2, na.rm=TRUE)



   plot(p2)
   ## Accuracy for secnario2 
   rf.pred<- predict(RFModel2,  sample22.df[,2:ncol(sample22.df)], type="response")
   
   obs.pred2 <- data.frame(cbind(Observed=as.numeric(as.factor(sample22.df[,"landusetyp"]))),
                           PRED=as.numeric(as.character(rf.pred3)))
   
   tablepred2<-table(obs.pred2$Observed,obs.pred2$PRED)## confusion matrix
   tablepred2
   confusionMatrix(tablepred2)
   koga_classification<-raster("D:\\gis certificate\\R-programing\\gully suceptibility\\data\\Rf_koga_lULC.tif")
   ## labling 
   ##the predicted land use land cover rastre layer has to be chnaged in to category for ploting and leveling 
   lulc<-as.factor( koga_classification)
   lulc
   rat <- levels(lulc)[[1]]
   rat
   ##assigning lnad cover land use class names for each value using the follwing code 
   rat[["landcover"]] <- c( "cultivated", "Forest", "Grazing","Shrub","Built_up_residential","Water")
   rat[["landcover"]]
   rat
   levels(lulc) <- rat
   library(rasterVis)
   
   levelplot(lulc, col.regions=c("brown","darkgreen","cyan","lightgreen","pink","blue"), xlab="", ylab="")
   
  
   ## scenario3,ntree=700&mtry=3,with all variable except band 2 and band 3(removing less important variables)
set.seed(123)
RFModel3 <-randomForest(x= selectedVar[,2:ncol(selectedVar)],
                           y = as.factor(selectedVar$landusetyp),proximity=TRUE,
                           
                           ntree = 700,mtry=3,replace = TRUE, importance = TRUE)
   
   RFModel3 
   ## prediction @ scenario3
   
   preiction4<-predict(covs,RFModel3,format="GTiff","koga_finalLULC",
                      type="response", index=2,
                      na.rm=TRUE, overwrite=TRUE, progress="text")
   
   
   plot(preiction4)
   
   
   
   plot(preiction31)
   preiction31
   str(preiction31)
# ## acuracy assessment ---------------------------------------------------

   
  
   layout(matrix(c(1,2),nrow=1),
          width=c(4,1)) 
   par(mar=c(5,4,4,0)) #No margin on the right side
   plot(RFModel1)
   par(mar=c(5,0,4,2)) #No margin on the left side
   plot(c(0,1),type="n", axes=F, xlab="", ylab="")
 
   legend("topleft", legend=unique(sample22.df$landusetyp), col=unique(as.numeric(sample22.df$landusetyp)), pch=19)
   
   legend("top", colnames(RFModel1$err.rate),col=1:4,cex=0.8,fill=1:4)
   ### Tune randomForest for the optimal mtry parameter
   set.seed(123)
   tuneRF(sample22.df1[,2:ncol(sample22.df1)], as.factor(sample22.df1$landusetype), c(1:9), ntree=700, stepFactor=2, improve=0.05,
          trace=TRUE, plot=TRUE, doBest=FALSE)

# ## RF predictioon  ------------------------------------------------------

   set.seed(123)
   predicted2<-predict(Koga_sentinel2_subseted,finalModel,type="response",na.rm=TRUE)
## plotbvaribakle importance 
   varImpPlot(finalModel2)
   ## estimating the frequency of variable used 
   varUsed(finalModel2,by.tree = FALSE,count = TRUE)


# prediction and lable at three scenario ; accuracy assessmnet  -----------

   modelfit --------------------------------------------------------------
     
     rf.pred <- predict(finalModel, training.df@data[,2:ncol(training.df@data)], type="response",na.action=na.pass)
   rf.pred
   rf.prob <- as.data.frame(predict(finalModel, training.df@data[,2:ncol(training.df@data)], type="prob"))
   rf.prob
   obs.pred <- data.frame(cbind(Observed=as.numeric(as.character(training.df@data[,"landusetype"])),
                                PRED=as.numeric(as.character(rf.pred)), Prob1=rf.prob[,2],Prob2=rf.prob[,3],Prob3=rf.prob[,4],Prob4=rf.prob[,5],Prob5=rf.prob[,6],
                                Prob0=rf.prob[,1]))
   class(covs)
   library(caret)## needed to estaimate confusion matrix 
   op <- (obs.pred$Observed == obs.pred$PRED)
   op 
   tablepred<-table(obs.pred$Observed,obs.pred$PRED)## confusion matrix
   tablepred
   confusionMatrix(op)
   confusionMatrix(tablepred)
   ( pcc <- (length(op[op == "TRUE"]) / length(op))*100 )
   install.packages("verification")
   library(verification)
   library(rfUtilities)## needed to estimate acuracy 
   library(OOBCurve)
   accuracy
   accuracy(obs.pred$Observed,obs.pred$PRED)## calssification accuracy
   acc
   obs.pred$Observed
   str(obs.pred$PRED)
   
   
   # labling and ploting RF product ------------------------------------------
   
   rf.prob
   ## ploting 
   plot(predict)
   ## 
   plot(finalModel, main="Bootstrap Error Convergence")
   p <- as.matrix(finalModel$importance)   
   plot(p)
   
   ## calling the rf product-final classification 
   Koag_landuselandcover<-raster("koga-landuselandcover.tif")
   varImpPlot(finalModel)
   
   plot(Koag_landuselandcover)

   
   
   
## Maximum Likelihood classification using RStoolbox pack 
##fitting the model by splitting the data in to 70%(trainig set),and 30% test set
## assiginig the same projection  for thje  validation and test set as follows
sampl_super<- spTransform(fores,CRS(proj4string( Koga_sentinel2_subseted)))
## Partition the sample in to trainingand test data set 
ind<-sample(2,nrow(sampl_super),replace = TRUE,prob = c(0.7,0.3))
train.df<-    sampl_super[ind==1, ]

test.df<-    sampl_super[ind==2, ]
str(train.df)
#creating test data set by not selecting the output row values
test<-data[-data1,]
MLC2<-superClass(Koga_sentinel2_subseted,trainData=train.df,responseCol="landusetyp",minDist = 1,
                   model = "rf",valData=  test.df, kfold = 5,format="GTiff",mode = "classification",
                   predict = TRUE,predType="raw",filename = "mximum_classifiaction")

super
   
##  study area shape 
Koga_shape<-readOGR("d:\\gis certificate\\R-programing\\gully suceptibility\\data\\mosaic\\KOGA_SHAPE.shp")
## crop the result to the extenet of koga watershed
## transformin projection
koga_shape<-proj4string(Koga_shape) <- CRS("+proj=utm +zone=37 +datum=WGS84 +units=m +no_defs")
koga_shape
## extract cell value to the extyenet of koga watershed
library(rgeos)
## Clipping raster using shapefile in R, but keeping the geometry of the shapefile
## first doing crooping using the extenet 

# lableing and ploting  ---------------------------------------------------

# labling 
##the predicted land use land cover rastre layer has to be chnaged in to category for ploting and leveling 
lulc<-as.factor( MLC$map)
lulc
rat <- levels(lulc)[[1]]
rat
##assigning lnad cover land use class names for each value using the follwing code 
rat[["landcover"]] <- c( "cultivated", "Forest", "Grazing","Shrub","Builtup","Water")
rat[["landcover"]]
rat
levels(lulc) <- rat
library(rasterVis)

levelplot(koga_lulc_mask1, col.regions=c("brown","darkgreen","cyan","lightgreen","pink","blue"), xlab="", ylab="")




# ## statistic and acuracy assessment  ------------------------------------


KOgaLUL1<-crop(lulc,extent(Koga_shape))##croping
koga_lulc_mask1<-mask(KOgaLUL1,Koga_shape)## landus eland cover to the exetent of koga watrehsed 

plot(koga_lulc_mask1)
MLC$classMapping

f1<-freq(koga_lulc_mask1, useNA='no')## frequesncy of cellvalue at each category
f1
produt_resolution<-prod(res(koga_lulc_mask1))
produt_resolution
landuselandcover_area<-cbind(f1, area=f1[,2] * produt_resolution)
landuselandcover_area


## accuracy Assessment
library(caret)## needed to estaimate confusion matrix 
library(rfUtilities)## needed to estimate acuracy 

 ## extract validation result from supoer class such as  user and producer accuracyusing get validation
table1<-getValidation(MLC, from = "testset", metrics='confmat')

## then following the confusioon matrix, estimate the user and producer accuracy as follows
# User accuracy
diag(table1) /rowSums(table1)## user accuracy 

diag(table1) /colSums(table1)## producer accuracy 

















   
## creating confusion matrix for   Maximum likelihood classifires 
   creating confusion matrixs as follows 
   rf.pred1 <- predict(MLC$model,train.df,type="raw",na.action=na.pass)
   obs.pred1 <- data.frame(cbind(Observed=as.numeric(as.character(train.df@data[,"landusetyp"])),
                                 PRED=as.numeric(as.character(rf.pred1))))


getValidation(MLC, from = "testset", metrics = "classwise")
   confusionMatrix(xx)
getva
   MLC$modelFit
   
   ##   validation  IN RS toolbox 
XX<-validateMap(MLC$map, test.df, "landusetyp", 
               mode = "classification", classMapping = MLC$classMapping)

 
 plot(MLC$map)
   library(caret)

   tablepred1<-table(obs.pred1$Observed,obs.pred1$PRED)## confusion matrix
   tablepred1          
## then following the confusioon matrix, estimate the user and producer accuracy as follows
# User accuracy
   diag(tablepred1) /rowSums(tablepred1)## user accuracy 
   
   diag(tablepred1) /colSums(tablepred1)## producer accuracy 
     
   shr<-readOGR("D:\\gis certificate\\R-programing\\gully suceptibility\\data\\data\\shrubland .shp")
   cu<-readOGR("D:\\gis certificate\\R-programing\\gully suceptibility\\data\\data\\cultivated land .shp")
   fores<-readOGR("D:\\gis certificate\\R-programing\\gully suceptibility\\data\\data\\Forest.shp")
   gra<-readOGR("D:\\gis certificate\\R-programing\\gully suceptibility\\data\\data\\grazing land .shp")
   wat<-readOGR("D:\\gis certificate\\R-programing\\gully suceptibility\\data\\data\\Water_body .shp")
   samp_best1<-rbind(shr,cu,fores,gra,wat,Bui1)
   Bui1<-readOGR("D:\\gis certificate\\R-programing\\gully suceptibility\\data\\data\\builtup.shp")
   samp_best3<- samp_best1[,-1]
   ## projection transform,ation
   Bui12<-spTransform(Bui1,CRS(proj4string( Koga_sentinel2_subseted)))   
   cu2<-spTransform(cu,CRS(proj4string( Koga_sentinel2_subseted)))
   fores2<-spTransform(fores,CRS(proj4string( Koga_sentinel2_subseted)))
   gra2<-spTransform(gra,CRS(proj4string( Koga_sentinel2_subseted)))
  wat2<-spTransform(wat,CRS(proj4string( Koga_sentinel2_subseted)))
  shr2<-spTransform(shr,CRS(proj4string( Koga_sentinel2_subseted)))
   
  
   
   ##ploting shaps 
   plot(wat, col='blue')
   plot(shr, add=TRUE, col='red')
   plot(gra,add=TRUE,col="lightgreen")
   plot(cu,add=TRUE,col="brown")
   plot(Bui1,add=TRUE,col="yellow")
   plot(fores,add=TRUE,col="darkgreen")
   summary(shr2)
   summary(wat2)
   summary(Bui12)
   summary(fores2)
   summary(cu2)
   summary(gra2)
   ab1 <- rbind(shr2,wat2,cu2, Bui12,fores2,gra2, makeUniqueIDs = TRUE)
   ab1
   
   ab1@data

   
   library(devtools)
   install_git("git://github.com/gsk3/taRifx.geo.git")
   library(taRifx.geo)
   rbind(a,b, fix.duplicated.IDs=TRUE)
   rbind(a,b, fix.duplicated.IDs=TRUE)
   rbind(a,b, fix.duplicated.IDs=TRUE)
   ##
   setwd("C:\\Users\\SDI\\Desktop\\EOforcrop and pasture")
   getwd()
   
   koga_region geroup<-raster()
   
   xx<-dir(path="C:\\Users\\SDI\\Desktop\\EOforcrop and pasture",pattern = ".hdf",all.files = TRUE)
   xx
   Afpar<-gdal_translate(xx[1],dst_dataset="MODIS.tif")
  M<-raster("MODIS.tif")
Afpar<-raster("C:\\Users\\SDI\\Desktop\\EOforcrop and pasture\\MOD15A2H.A2016361.h22v08.006.2017010135407.hdf")
   
install.packages(GDAL)
gdal_setInstallation

Modtiff<-dir(path = "D:\\mODIS",pattern=".tif",all.files = TRUE)
Modtiff
M<-raster::resample(Modtiff[1],250,method="bilinear")
install.packages("hegWin")
MODIStsp()
x<-runMrt(product="MOD15A2H",outproj= "UTM",Jobs="fPAR",zone =37,36,38 ,datum = "WGS84",mosaic = TRUE,begin= "2016169",end  = "2016361",extent = "Ethiopia")

   y<-runMrt(product="MOD15A2H",outproj= "UTM",Jobs="fPAR",zone =37,datum = "WGS84",mosaic = TRUE,begin= "2016169",end  = "2016361",extent = "Ethiopia")
   
   )

library(MODIS)
runGdal(
  job        = 'fPAR',
  product    = 'MOD15A2H', 
  SDSstring  = "1",         
  collection = '006', 
  extent     = 'Ethiopia', 
  begin      = "2015.06.01", 
  end        = "2015.12.28",
  mosaic     = TRUE,
  outDirPath = '.'      # so that we don't have to look for the PROCESSED folder.
)

library(MODIS)
checkDeps()
devtools::install_github("MatMatt/MODIS", ref = "develop")
library(MODIS)

lap = file.path(tempdir(), "MODIS_ARC")
odp = file.path(lap, "PROCESSED")
MODISoptions(lap, odp)
MODIS:::checkDeps()
MODIS:::checkTools()
MODISoptions
install.packages(GDAL)
MODIS:::checkTools('GDAL')
MODISoptions(gdalPath='/Path/to/XXGDAL/bin')

setRepositories() # activate CRAN, R-forge, and Omegahat and then: 
install.packages(c(' ptw '),dependencies=TRUE)
install.packages("rgdal", dependencies = TRUE)
1# activate CRAN, R-forge, and Omegahat and then: 
install.packages(c(' ptw '),dependencies=TRUE)

To install all required and suggested packages run:
  setRepositories()
  # activate CRAN, R-forge, and Omegahat and then: 
install.packages(c(' ptw '),dependencies=TRUE)

sh: gdalinfo: command not found
'MRT_HOME' not set/found! MRT is NOT enabled! See: 'https://lpdaac.usgs.gov/tools/modis_reprojection_tool'

Sys.getenv("PATH")

Sys.setenv("C:\\Program Files\\R\\R-3.5.0\\bin\\x64;C:\\Windows \\system32;C:\\$SEN2COR_BIN\\aux_data\\Sen2Cor-02.05.05-win64;C:\\Program Files\\snap\\bin;C:\\Users\\SDI\\Documents\\sen2cor;C:\\$SEN2COR_BIN\\aux_data\\Sen2Cor-02.05.05-win64;C:\\\\Program Files (x86)\\\\bin\\\\gdalinfo.exe;C:\\Program Files (x86)\\bin\\gdal_translate (2).exe")


##
landusegen<-raster("C:\\Users\\SDI\\Desktop\\data_code _selamawit\\landusegen.tif")
landusegen
landuse_ni_ver<-raster("C:\\Users\\SDI\\Desktop\\data_code _selamawit\\lanuse_clipe1.tif")
##calculate area
f<-freq(landusegen, useNA='no')## frequesncy of cellvalue at each category
f
produt_resolution<-prod(res(landusegen))
produt_resolution
landuselandcover_area<-cbind(f, area=f[,2] * produt_resolution)
landuselandcover_area
##calculate area2
f1<-freq(landuse_ni_ver, useNA='no')## frequesncy of cellvalue at each category
f1
produt_resolution1<-prod(res(landuse_ni_ver))
produt_resolution1
landuselandcover_area1<-cbind(f1, area=f1[,2] * produt_resolution1)
landuselandcover_area1
##############
Lc_vertI_rast1
landuse_nverti<-raster("C:\\Users\\SDI\\Desktop\\data_code _selamawit\\Lc_vertI_rast1.tif")
landuse_nito<-raster("C:\\Users\\SDI\\Desktop\\data_code _selamawit\\LC_Nito_ras1.tif")
##calculate area_verti sols 
f3<-freq(landuse_nverti, useNA='no')## frequesncy of cellvalue at each category
f3
produt_resolution3<-prod(res(landuse_nverti))
produt_resolution3
landuselandcover_area_verti<-cbind(f3, area=f3[,2] * produt_resolution3)
landuselandcover_area_verti
##calculate area_Nito
f4<-freq(landuse_nito, useNA='no')## frequesncy of cellvalue at each category
f4
produt_resolution4<-prod(res(landuse_nito))
produt_resolution4
landuselandcover_area_nito<-cbind(f4, area=f4[,2] * produt_resolution4)
landuselandcover_area_nito

### vulnerability 
sanitation<-raster("E:\\dire_vuln\\san_recl1.tif")
Formof_sett<-raster("E:\\dire_vuln\\form_setl_rec1.tif")
shelter_status<-raster("E:\\dire_vuln\\shelter_rec2.tif")
hh_resid_rec1<-raster("E:\\dire_vuln\\hh_resid_rec1.tif")
electric_rec1<-raster("E:\\dire_vuln\\electric_rec1.tif")
WashandBasicervic<-addLayer(sanitation,Formof_sett,shelter_status,hh_resid_rec1,electric_rec1)<-addLayer(sanitation,Formof_sett,shelter_status,hh_resid_rec1,electric_rec1)
plot(WashandBasicervic)
### coivid report 
data.cov<-covid19.data(case='aggregated')
data.cov
tsa<-covid19.data(case = 'ts-ALL')

tsc<-covid19.data(case = 'ts-confirmed')
tsce<-covid19.data(case = 'ts-confirmed')

tot<-tots.per.location(tsce,geo.loc = 'Ethiopia')
## growth rate 
gr1<-growth.rate(tsa,geo.loc = 'Ethiopia')

## totla 
total<-covid19.data(case='ts-All')
plt<-totals.plt(tot)
totw<-totals.plt(total,c('Ethiopia'))
t<-totals.plt(tsc1,c('Ethiopia'))
l<-live.map(tsc)
## SIR model 
generate.SIR.model(tsc,'Ethiopia',tot.population = 114580493)

####
?(covid19.analytics)
