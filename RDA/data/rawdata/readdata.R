#
setwd("/home/lu/Documents/Github/Multi_NNGP/RDA/data/rawdata")
#MyData <- read.csv(file="GFB_Abm_ba_USA_10192017.csv", header=TRUE, sep=",")
#install.packages("data.table")
library("data.table")
#t <- proc.time()
#data <- fread('GFB_Abm_ba_USA_10192017.csv', header = T)
#proc.time() - t
#summary(data$Ps)
#summary(data$Yr)
#summary(data$B)
rm(list = ls())
library("MODIS")
library("rgdal")
#EarthdataLogin(usr = "****", pwd = "****")

############################
# download data from MODIS #
############################

# MCD64A1 - Combined Level 3 Direct Broadcast Burned Area Monthly Global 500m SIN Grid
search4map(pattern="USA",plot=TRUE)
search4map(database="state",plot=TRUE)
t <- proc.time() 
burnarea <- getHdf(product = "MCD64A1", begin = "2018313", end = "2018333",
                   extent = "california")
proc.time() - t

# MCD12Q1 - MODIS/Terra+Aqua Land Cover Type Yearly L3 Global 500m SIN Grid
t <- proc.time() 
landcover <- getHdf(product = "MCD12Q1", begin = "2016001", end = "2017001",
                   extent = "california")
proc.time() - t
landcover

# MOD13A1 - MODIS/Terra Vegetation Indices 16-Day L3 Global 500m SIN Grid
t <- proc.time() 
VegeInd<- getHdf(product = "MOD13A1", begin = "2018304", end = "2018337",
                    extent = "california")
proc.time() - t
VegeInd

t <- proc.time() 
VegeInd2<- getHdf(product = "MOD13A1", begin = "2016100", end = "2016120",
                 extent = "california")
proc.time() - t
VegeInd2

#GPP
t <- proc.time() 
GPPdata<- getHdf(product = "MOD17A2H", begin = "2016100", end = "2016120",
                  extent = "california")
proc.time() - t
GPPdata

# atomesphere ??
# t <- proc.time() 
# vaperozone<- getHdf(product = "MYD08_M3", begin = "2018304", end = "2018337",
#                  extent = "california")
# proc.time() - t
# vaperozone

t <- proc.time() 
Net_Evapo<- getHdf(product = "MOD16A2", begin = "2016100", end = "2016120",
                 extent = "california")
proc.time() - t
Net_Evapo

# also check 1km by 1km 
# MOD14A2 - MODIS/Terra Thermal Anomalies/Fire 8-Day L3 Global 1km SIN Grid
# MOD44B - MODIS/Terra Vegetation Continuous Fields Yearly L3 Global 500m SIN Grid
# MOD11A2 - MODIS/Terra Land Surface Temperature/Emissivity 8-Day L3 Global 1km SIN Grid
# MOD17A2H - MODIS/Terra Gross Primary Productivity 8-Day L4 Global 500m SIN Grid
# MOD15A2H - MODIS/Terra Leaf Area Index/FPAR 8-Day L4 Global 500m SIN Grid
# MOD21A2 - MODIS/Terra Land Surface Temperature/Emissivity 8-Day L3 Global 1km SIN Grid
# MOD21A1D - MODIS/Terra Daytime Land Surface Temperature/Emissivity Daily L3 Global 1km SIN Grid
# MOD21A1N - MODIS/Terra Nighttime Land Surface Temperature/Emissivity Daily L3 Global 1km SIN Grid
# MOD16A2 - MODIS/Terra Net Evapotranspiration 8-Day L4 Global 500m SIN Grid

##################
# check hdf data #
##################
library(rgdal) 
library(gdalUtils) 
library(raster)
library(corrplot)

burnarea1 <- get_subdatasets('MCD64A1.006/2018.11.01/MCD64A1.A2018305.h08v05.006.2019011073844.hdf')
burnarea1
par(mfrow= c(1,1))
BA1date <- readGDAL(burnarea1[1])
sum(BA1date@data > -2.00 & !is.na(BA1date@data))
summary(BA1date@data[BA1date@data > -2.00  & !is.na(BA1date@data)])
summary(BA1date@data[BA1date@data > 0.00  & !is.na(BA1date@data)])
sum(BA1date@data > 0.00 & !is.na(BA1date@data))
plot(raster(BA1date))

BA1Fdate <- readGDAL(burnarea1[4])
sum(BA1Fdate@data > 0.00 & !is.na(BA1Fdate@data))
table(BA1Fdate@data$band1)
plot(raster(BA1Fdate))
plot(raster(BA1Fdate))

BA1Ldate <- readGDAL(burnarea1[5])
sum(BA1Ldate@data > 0.00 & !is.na(BA1Ldate@data))
duration <- BA1Ldate@data -BA1Fdate@data 
summary(duration)
sum(duration > 0.00 & !is.na(duration))

# good!
landcover1 <- get_subdatasets('MCD12Q1.006/2016.01.01/MCD12Q1.A2016001.h08v05.006.2018055070314.hdf')
landcover1
LC_Type1<- readGDAL(landcover1[1])
table(LC_Type1@data)

LC_Type2<- readGDAL(landcover1[2])
table(LC_Type2@data)

LC_Type3<- readGDAL(landcover1[3])
table(LC_Type3@data)

LC_Type4<- readGDAL(landcover1[4])
table(LC_Type4@data)

LC_Type5<- readGDAL(landcover1[5])
table(LC_Type5@data)

LC_Prop1_A<- readGDAL(landcover1[6])
str(LC_Prop1_A)
summary(LC_Prop1_A@data) # NA's   :1106969 
hist(LC_Prop1_A@data$band1)

LC_Prop2_A<- readGDAL(landcover1[7])
summary(LC_Prop2_A@data) # NA's   :1136228
hist(LC_Prop2_A@data$band1)

LC_Prop3_A<- readGDAL(landcover1[8])
summary(LC_Prop3_A@data) # NA's   :1106969 
hist(LC_Prop3_A@data$band1)

LC_Prop1<- readGDAL(landcover1[9])
summary(LC_Prop1@data) # NA's   :1106969 
table(LC_Prop1@data$band1)

LC_Prop2<- readGDAL(landcover1[10])
summary(LC_Prop2@data)
table(LC_Prop2@data$band1)

LC_Prop3<- readGDAL(landcover1[11])
summary(LC_Prop3@data) 
table(LC_Prop3@data$band1)

LC_QC <- readGDAL(landcover1[12])
summary(LC_QC@data) 
table(LC_QC@data$band1)
#plot(raster(LC_QC))
LC_LW <- readGDAL(landcover1[13])
summary(LC_LW@data) 
table(LC_LW@data$band1)

#~4653031

par(mfrow = c(3, 4))
plot(raster(LC_Type1), main = "LC_Type1")
plot(raster(LC_Type2), main = "LC_Type2")
plot(raster(LC_Type3), main = "LC_Type3")
plot(raster(LC_Type4), main = "LC_Type4")
plot(raster(LC_Type5), main = "LC_Type5")
plot(raster(LC_Prop1_A), main = "LC_Prop1_Assessment")
plot(raster(LC_Prop2_A), main = "LC_Prop2_Assessment")
plot(raster(LC_Prop3_A), main = "LC_Prop3_Assessment")
plot(raster(LC_Prop1), main = "LC_Prop1")
plot(raster(LC_Prop2), main = "LC_Prop2")
plot(raster(LC_Prop3), main = "LC_Prop3")
plot(raster(LC_LW), main = "LW")


#read in LC data
LC_data <- matrix(NA, nrow = 2400*2400, ncol = 15)
LC_data[, 1:2] =  coordinates(LC_Type1)[, 1:2] # read in the coordinates 
for (i in 1:13){
  LC_data[, i + 2] = readGDAL(landcover1[i])@data$band1
}
colnames(LC_data) = c("x", "y", "LC_Type1", "LC_Type2", "LC_Type3", "LC_Type4",
                      "LC_Type5", "LC_Prop1_Assessment", "LC_Prop2_Assessment",
                      "LC_Prop3_Assessment", "LC_Prop1", "LC_Prop2", "LC_Prop3",
                      "QC", "LW")
CORLC = cor(LC_data, use = "pairwise.complete.obs")
library(corrplot)
par(mfrow = c(1,1))
corrplot(CORLC, method="number")



# Vegetable indice data summarize
VegeInd1 <- get_subdatasets('MOD13A1.006/2016.04.06/MOD13A1.A2016097.h08v05.006.2016114041345.hdf')
VegeInd1

VI_NDVI <- readGDAL(VegeInd1[1])
str(VI_NDVI) 
summary(VI_NDVI@data$band1)
hist(VI_NDVI@data$band1)

VI_EVI <- readGDAL(VegeInd1[2])
summary(VI_EVI@data$band1)
hist(VI_EVI@data$band1)

VI_VI_Quality <- readGDAL(VegeInd1[3]) # cannot be response
summary(VI_VI_Quality@data$band1)
hist(VI_VI_Quality@data$band1)

VI_red <- readGDAL(VegeInd1[4])
summary(VI_red@data$band1)
hist(VI_red@data$band1)

VI_NIR <- readGDAL(VegeInd1[5])
summary(VI_NIR@data$band1)
hist(VI_NIR@data$band1)

VI_blue <- readGDAL(VegeInd1[6])
summary(VI_blue@data$band1)
hist(VI_blue@data$band1)

VI_MIR <- readGDAL(VegeInd1[7])
summary(VI_MIR@data$band1)
hist(VI_MIR@data$band1)

VI_VZA<- readGDAL(VegeInd1[8])
summary(VI_VZA@data$band1)
hist(VI_VZA@data$band1)

VI_SZA<- readGDAL(VegeInd1[9])
summary(VI_SZA@data$band1)
hist(VI_SZA@data$band1)

VI_RAA<- readGDAL(VegeInd1[10])
summary(VI_RAA@data$band1)
hist(VI_RAA@data$band1)

VI_CDY<- readGDAL(VegeInd1[11])
summary(VI_CDY@data$band1)
hist(VI_CDY@data$band1)

VI_PR<- readGDAL(VegeInd1[12])
summary(VI_PR@data$band1)
table(VI_PR@data$band1)

# check plots
par(mfrow = c(3, 4))
plot(raster(VI_NDVI), main = "NDVI")
plot(raster(VI_EVI), main = "EVI")
plot(raster(VI_VI_Quality), main = "VI Quality") # clean data with low quality
plot(raster(VI_red), main = "red reflectance")
plot(raster(VI_NIR), main = "NIR reflectance")
plot(raster(VI_blue), main = "blue reflectance")
plot(raster(VI_MIR), main = "MIR reflectance")
plot(raster(VI_VZA), main = "view zenith angle")
plot(raster(VI_SZA), main = "sun zenith angle")
plot(raster(VI_RAA), main = "relative azimuth angle") # categorical data?
plot(raster(VI_CDY), main = "composite day of the year")
plot(raster(VI_PR), main = "pixel reliability") # categorical data


#read in VI data
VI_data <- matrix(NA, nrow = 2400*2400, ncol = 14)
scale_factor <- c(0.0001, 0.0001, 1, 0.0001,0.0001,0.0001, 0.0001, 
                  0.01, 0.01, 0.01, 1, 1, 1)
VI_data[, 1:2] =  coordinates(VI_NDVI)[, 1:2] # read in the coordinates 
for (i in 1:12){
  VI_data[, i + 2] = readGDAL(VegeInd1[i])@data$band1 * scale_factor[i]^2 
  # scale the data 
}
colnames(VI_data) = c("x", "y", "NDVI", "EVI", "VI Quality", "red reflectance",
                      "NIR reflectance", "blue reflectance", "MIR reflectance", 
                      "view zenith angle", "sun zenith angle", 
                      "relative azimuth angle", "composite day of the year",
                      "pixel reliability")

CORVI = cor(VI_data, use = "pairwise.complete.obs")
library(corrplot)
par(mfrow = c(1,1))
corrplot(CORVI, method="circle")

# Check GPP data #
GPP1 <- get_subdatasets('MOD17A2H.006/2016.04.06/MOD17A2H.A2016097.h08v05.006.2016204054521.hdf')
GPP1
GPP2 <- get_subdatasets('MOD17A2H.006/2016.04.14/MOD17A2H.A2016105.h08v05.006.2016204154944.hdf')
GPP2
par(mfrow = c(1, 2))
Gpp_Gpp1 <- readGDAL(GPP1[1])
summary(Gpp_Gpp1@data$band1)
table(Gpp_Gpp1@data$band1[which(Gpp_Gpp1@data$band1 > 3.00)])
#  3.2762  3.2765  3.2766 
#  92337  396397 1256876 
# change the unclassified as NA and the rest as 0.0
length(which(Gpp_Gpp1@data$band1 == 3.2766))
length(which(Gpp_Gpp1@data$band1 == 3.2765))
Gpp_Gpp1@data$band1[which(Gpp_Gpp1@data$band1 == 3.2766)] <- NA
Gpp_Gpp1@data$band1[which(Gpp_Gpp1@data$band1 == 3.2765)] <- 0.0
length(which(Gpp_Gpp1@data$band1 > 3.2761))
Gpp_Gpp1@data$band1[which(Gpp_Gpp1@data$band1 > 3.2761)] <- 0.0
sum(!is.na(Gpp_Gpp1@data))
plot(raster(Gpp_Gpp1))
Gpp_Gpp2 <- readGDAL(GPP2[1])
summary(Gpp_Gpp2@data$band1)
table(Gpp_Gpp2@data$band1[which(Gpp_Gpp2@data$band1 > 3.00)])
length(which(Gpp_Gpp2@data$band1 == 3.2766))
length(which(Gpp_Gpp2@data$band1 == 3.2765))
Gpp_Gpp2@data$band1[which(Gpp_Gpp2@data$band1 == 3.2766)] <- NA
Gpp_Gpp2@data$band1[which(Gpp_Gpp2@data$band1 == 3.2765)] <- 0.0
length(which(Gpp_Gpp2@data$band1 > 3.2761))
Gpp_Gpp2@data$band1[which(Gpp_Gpp2@data$band1 > 3.2761)] <- 0.0
sum(!is.na(Gpp_Gpp2@data))
hist(Gpp_Gpp2@data$band1[Gpp_Gpp2@data$band1 < 3])
sum(Gpp_Gpp2@data < 3.00 & !is.na(Gpp_Gpp2@data))
plot(raster(Gpp_Gpp2))

#Gpp_Gpp12 = Gpp_Gpp1
#Gpp_Gpp22 = Gpp_Gpp2
#Gpp_Gpp12@data$band1[Gpp_Gpp12@data$band1 >= 3] <- NA  #46656
#Gpp_Gpp22@data$band1[Gpp_Gpp22@data$band1 >= 3] <- NA  #46656
#plot(raster(Gpp_Gpp12))
#plot(raster(Gpp_Gpp22))

Gpp_PsnNet1 <- readGDAL(GPP1[2])
summary(Gpp_PsnNet1@data$band1)
table(Gpp_PsnNet1@data$band1[which(Gpp_PsnNet1@data$band1 > 3.00)])
#  3.2762  3.2765  3.2766 
#  92337  396397 1256876 
# change the unclassified as NA and the rest as 0.0
length(which(Gpp_PsnNet1@data$band1 == 3.2766))
length(which(Gpp_PsnNet1@data$band1 == 3.2765))
Gpp_PsnNet1@data$band1[which(Gpp_PsnNet1@data$band1 == 3.2766)] <- NA
Gpp_PsnNet1@data$band1[which(Gpp_PsnNet1@data$band1 == 3.2765)] <- 0.0
length(which(Gpp_PsnNet1@data$band1 > 3.2761))
Gpp_PsnNet1@data$band1[which(Gpp_PsnNet1@data$band1 > 3.2761)] <- 0.0
hist(Gpp_PsnNet1@data$band1)
plot(raster(Gpp_PsnNet1))
Gpp_PsnNet2 <- readGDAL(GPP2[2])
summary(Gpp_PsnNet2@data$band1)
table(Gpp_PsnNet2@data$band1[which(Gpp_PsnNet2@data$band1 > 3.00)])
length(which(Gpp_PsnNet2@data$band1 == 3.2766))
length(which(Gpp_PsnNet2@data$band1 == 3.2765))
Gpp_PsnNet2@data$band1[which(Gpp_PsnNet2@data$band1 == 3.2766)] <- NA
Gpp_PsnNet2@data$band1[which(Gpp_PsnNet2@data$band1 == 3.2765)] <- 0.0
length(which(Gpp_PsnNet2@data$band1 > 3.2761))
Gpp_PsnNet2@data$band1[which(Gpp_PsnNet2@data$band1 > 3.2761)] <- 0.0
hist(Gpp_PsnNet2@data$band1)
plot(raster(Gpp_PsnNet2))

Gpp_PsnQC1 <- readGDAL(GPP1[3])
table(Gpp_PsnQC1@data$band1)
plot(raster(Gpp_PsnQC1))
Gpp_PsnQC2 <- readGDAL(GPP2[3])
table(Gpp_PsnQC2@data$band1)
plot(raster(Gpp_PsnQC2))
summary(Gpp_PsnQC2@data$band1 - Gpp_PsnQC1@data$band1) # different ... just ignore it..




#read in GPP data, average over two 8-days periods 
Gpp_data <- matrix(NA, nrow = 2400*2400, ncol = 2 + 2)
Gpp_data[, 1:2] =  coordinates(Gpp_Gpp1)[, 1:2] # read in the coordinates 
Gpp_data[, 3] = (Gpp_Gpp1@data$band1 + Gpp_Gpp2@data$band1)/2.0
Gpp_data[, 4] = (Gpp_PsnNet1@data$band1 + Gpp_PsnNet2@data$band1)/2.0
colnames(Gpp_data) = c("x", "y", "GPP", "PsnNet")

CORGpp = cor(Gpp_data, use = "pairwise.complete.obs")
par(mfrow = c(1,1))
corrplot(CORGpp, method="ellipse")


# check Net Evapo data
NetEvapo1<- get_subdatasets('MOD16A2.006/2016.04.06/MOD16A2.A2016097.h08v05.006.2017119005111.hdf')
NetEvapo1
NetEvapo2<- get_subdatasets('MOD16A2.006/2016.04.14/MOD16A2.A2016105.h08v05.006.2017119030738.hdf')
NetEvapo2

par(mfrow = c(1, 1))
ET1 <- readGDAL(NetEvapo1[1]) 
ET1@data$band1 = ET1@data$band1 
summary(ET1@data$band1)
table(ET1@data$band1[which(ET1@data$band1 > 3276.05)])
ET1@data$band1[which(ET1@data$band1 > 3276.05)] = NA # Unclassified
summary(ET1@data$band1)
hist(ET1@data$band1)
hist(log(ET1@data$band1))
sum(!is.na(ET1@data$band1)) #3212852
plot(raster(ET1))

ET2 <- readGDAL(NetEvapo2[1]) 
ET2@data$band1 = ET2@data$band1 
summary(ET2@data$band1)
table(ET2@data$band1[which(ET2@data$band1 > 3276.05)])
ET2@data$band1[which(ET2@data$band1 > 3276.05)] = NA # Unclassified
summary(ET2@data$band1)
hist(ET2@data$band1)
hist(log(ET2@data$band1))
sum(!is.na(ET2@data$band1)) #3212852
plot(raster(ET2))

PET1 <- readGDAL(NetEvapo1[3]) 
PET1@data$band1 = PET1@data$band1 
summary(PET1@data$band1)
table(PET1@data$band1[which(PET1@data$band1 > 3276.05)])
PET1@data$band1[which(PET1@data$band1 > 3276.05)] = NA # Unclassified
summary(PET1@data$band1)
hist(PET1@data$band1)
sum(!is.na(PET1@data$band1)) #3212852
plot(raster(PET1))

PET2 <- readGDAL(NetEvapo2[3]) 
PET2@data$band1 = PET2@data$band1 
summary(PET2@data$band1)
table(PET2@data$band1[which(PET2@data$band1 > 3276.05)])
PET2@data$band1[which(PET2@data$band1 > 3276.05)] = NA # Unclassified
summary(PET2@data$band1)
hist(PET2@data$band1)
sum(!is.na(PET2@data$band1)) #3212852
plot(raster(PET2))

LE1 <- readGDAL(NetEvapo1[2])
LE1@data$band1 <- LE1@data$band1 * 0.0001
summary(LE1@data$band1)
LE1@data$band1[which(LE1@data$band1 > 32760.5)] = NA # Unclassified
hist(LE1@data$band1)
sum(!is.na(LE1@data$band1))
summary(LE1@data$band1)
plot(raster(LE1))

LE2 <- readGDAL(NetEvapo2[2])
LE2@data$band1 <- LE2@data$band1 * 0.0001
summary(LE2@data$band1)
LE2@data$band1[which(LE2@data$band1 > 32760.5)] = NA # Unclassified
hist(LE2@data$band1)
sum(!is.na(LE2@data$band1))
summary(LE2@data$band1)
plot(raster(LE2))

PLE1 <- readGDAL(NetEvapo1[4])
PLE1@data$band1 <- PLE1@data$band1 * 0.0001
summary(PLE1@data$band1)
PLE1@data$band1[which(PLE1@data$band1 > 32760.5)] = NA # Unclassified
hist(PLE1@data$band1)
sum(!is.na(PLE1@data$band1))
summary(PLE1@data$band1)
plot(raster(PLE1))

PLE2 <- readGDAL(NetEvapo2[4])
PLE2@data$band1 <- PLE2@data$band1 * 0.0001
summary(PLE2@data$band1)
PLE2@data$band1[which(PLE2@data$band1 > 32760.5)] = NA # Unclassified
hist(PLE2@data$band1)
sum(!is.na(PLE2@data$band1))
summary(PLE2@data$band1)
plot(raster(PLE2))


#read in GPP data, average over two 8-days periods 
LEET_data <- matrix(NA, nrow = 2400*2400, ncol = 2 + 4)
LEET_data[, 1:2] =  coordinates(LE1)[, 1:2] # read in the coordinates 
LEET_data[, 3] = (LE1@data$band1 + LE2@data$band1)/2.0
LEET_data[, 4] = (ET1@data$band1 + ET2@data$band1)/2.0
LEET_data[, 5] = (PLE1@data$band1 + PLE2@data$band1)/2.0
LEET_data[, 6] = (PET1@data$band1 + PET2@data$band1)/2.0
colnames(LEET_data ) = c("x", "y", "LE", "ET", "PLE", "PET")

CORLEET = cor(LEET_data, use = "pairwise.complete.obs")
par(mfrow = c(1,1))
corrplot(CORLEET, method="number")



# check coordinates #
tt1 <- coordinates(LC_Type1); tt2 <- coordinates((VI_NDVI));
tt3 <- coordinates(Gpp_Gpp1); tt4 <- coordinates(ET1)
summary(tt1 - tt2) # same coordinates
summary(tt2 - tt3) # same coordinates
summary(tt3 - tt4) # same coordinates
#4512319
str(VI_NDVI@grid@cellcentre.offset)
# check the area for prediction
VI_NDVI2 = VI_NDVI
Pred_indx = (VI_data[, "x"] > -10400000 & VI_data[, "x"] < -10300000 &
               VI_data[, "y"] > 3800000 & VI_data[, "y"]< 3900000)
# take off a small area 
VI_NDVI2@data$band1[Pred_indx] <- NA  #46656
plot(raster(VI_NDVI2))



# union all data and scale the grids
scaled_coords = cbind(VI_data[, "x"] - VI_NDVI@grid@cellcentre.offset["x"],
                      VI_data[, "y"] - VI_NDVI@grid@cellcentre.offset["y"]) / 1000000
#the distance scale is in 1000km
summary(scaled_coords)
colnames(scaled_coords) = c("scaled_x", "scaled_y")
raw_data = as.data.frame(cbind(VI_data, Gpp_data[, -(1:2)], LC_data[, -(1:2)], 
                               LEET_data[, -(1:2)], scaled_coords))

#save(raw_data, file = "raw_data.RData")
#load("raw_data.RData")
#Raw_corr = cor(raw_data, use = "pairwise.complete.obs")
Raw_corr_resp = 
  cor(raw_data[, c( c("NDVI", "EVI","red reflectance", 
                      "NIR reflectance", "blue reflectance", "MIR reflectance", 
                      "GPP", "PsnNet", 
                      #"view zenith angle", "sun zenith angle", 
                      #"relative azimuth angle"
                      "LE", "ET", "PLE", "PET"
                      ))], use = "pairwise.complete.obs")
par(mfrow = c(1,1))
#corrplot(Raw_corr , method="number")
corrplot(Raw_corr_resp , method="number")
corrplot(Raw_corr_resp , method="circle")
# clean data 
NA_index = which(raw_data[, "LC_Type4"] == 255 | raw_data[, "LC_Type4"] == 0 |
                 is.na(raw_data[, "LC_Type4"]) |
                 is.na(raw_data[, "EVI"]) | 
                 is.na(raw_data[, "MIR reflectance"]) |
                 is.na(raw_data[, "relative azimuth angle"]) | 
                 is.na(raw_data[, "GPP"]) | raw_data[, "LW"] == 1.0 |
                 is.na(raw_data[, "LE"]) | is.na(raw_data[, "ET"]) | 
                   is.na(raw_data[, "PLE"]) | is.na(raw_data[, "PET"]))
data_cleaned = 
  raw_data[-NA_index, c("x", "y", "scaled_x", "scaled_y", "NDVI", "EVI", 
                        "red reflectance", "NIR reflectance", 
                        "blue reflectance", "MIR reflectance", "GPP",
               "PsnNet", "view zenith angle", "sun zenith angle", 
               "relative azimuth angle", "composite day of the year", 
               "LE", "ET", "PLE", "PET",
               "LC_Type4")]
dim(data_cleaned)
table(data_cleaned[, "LC_Type4"])
# "0" Water Bodies (At least 60% of area is covered by permanent water bodies)
# "1" Evergreen Needleleaf Vegetation: The reference group
# "2" Evergreen Broadleaf Vegetation
# "3" Deciduous Needleleaf Vegetation
# "4" Deciduous Broadleaf Vegetation
# "5" Annual Broadleaf Vegetation
# "6" Annual Grass Vegetation
# "7" Non-Vegetated Lands
# "8" Urban and Built-up Lands

#data_cleaned$Water_Bodies = as.numeric(data_cleaned$LC_Type4 == 0)
data_cleaned$Evergreen_Needleleaf_Vegetation = 
  as.numeric(data_cleaned$LC_Type4 == 1)
data_cleaned$Evergreen_Broadleaf_Vegetation = 
  as.numeric(data_cleaned$LC_Type4 == 2)
data_cleaned$Deciduous_Needleleaf_Vegetation = 
  as.numeric(data_cleaned$LC_Type4 == 3)
data_cleaned$Deciduous_Broadleaf_Vegetation = 
  as.numeric(data_cleaned$LC_Type4 == 4)
data_cleaned$Annual_Broadleaf_Vegetation = 
  as.numeric(data_cleaned$LC_Type4 == 5)
data_cleaned$Annual_Grass_Vegetation = 
  as.numeric(data_cleaned$LC_Type4 == 6)
data_cleaned$Non_Vegetated_or_Builtup_Lands = 
  as.numeric(data_cleaned$LC_Type4 == 7 | data_cleaned$LC_Type4 == 8)
#data_cleaned$Urban_and_Builtup_Lands = 
#  as.numeric(data_cleaned$LC_Type4 == 8)

save(data_cleaned, file = "cleaned_data_expanded.RData")
