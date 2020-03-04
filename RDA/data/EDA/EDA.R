setwd("/home/lu/Documents/Github/Multi_NNGP/RDA/data/rawdata")
library("data.table")
rm(list = ls())
library(rgdal) 
library(gdalUtils) 
library(raster)
library(corrplot)
#library(coda)
#library(spBayes)
#library(MBA)
#library(fields)
#library(classInt)
#library(RColorBrewer)
#library(sp)


load("cleaned_data_expanded.RData")
colnames(data_cleaned)

library(geoR)
library(gstat)
N = dim(data_cleaned)[1]
set.seed(123)
subind <- sample.int(N, round(N*0.01))
d.max <- sqrt((max(data_cleaned$scaled_x) - min(data_cleaned$scaled_x))^2 + 
                (max(data_cleaned$scaled_y) - min(data_cleaned$scaled_y))^2)
d.max # around 1572 KM


# check NDVI variogram
hist(data_cleaned[, c("NDVI")])
hist(log(data_cleaned[, c("NDVI")] + 1))
summary(log(data_cleaned[, c("NDVI")] + 1))
outliner.NDVI <- which(data_cleaned[, c("NDVI")] < -0.01)
length(outliner.NDVI)
hist(log(data_cleaned[-outliner.NDVI, c("NDVI")] + 1))
summary(log(data_cleaned[-outliner.NDVI, c("NDVI")] + 1))
v.NDVI<- variog(coords = data_cleaned[subind, c("scaled_x", "scaled_y")], 
                data = log(data_cleaned[subind, c("NDVI")] + 1), 
                uvec = (seq(0, 0.5 * d.max, length = 30))) # 30
par(mfrow=c(1,1))
varioNDVI.fit <- variofit(v.NDVI, cov.model="exponential")
summary(varioNDVI.fit)
plot(v.NDVI, xlab = "Distance (1000km)", cex = 0.1)
lines(varioNDVI.fit, col='red') 
abline(h = c(varioNDVI.fit$nugget, varioNDVI.fit$nugget + varioNDVI.fit$cov.pars[1]))
abline(v = 3 * varioNDVI.fit$cov.pars[2])
variofitphiNDVI <- 1 / varioNDVI.fit$cov.pars[2]; variofitphiNDVI # 3.685804

data.NDVI = data_cleaned[subind, c("scaled_x", "scaled_y", "NDVI", 
                                        "Non_Vegetated_or_Builtup_Lands")]
coordinates(data.NDVI) = ~scaled_x+scaled_y
v.NDVI2 <- variogram( NDVI~scaled_x+scaled_y+Non_Vegetated_or_Builtup_Lands, data.NDVI) # 30
v.NDVI2.fit = fit.variogram(v.NDVI2, vgm("Exp"));v.NDVI2.fit
#model       psill     range
#1   Nug 0.006001477 0.0000000
#2   Exp 0.018345890 0.2219178 # phi = 4.7
plot(v.NDVI2) 
range(data_cleaned[-outliner.NDVI, c("NDVI")])
sqrt(var(data_cleaned[-outliner.NDVI, c("NDVI")]))
#-0.0098  0.9476

v.NDVI3 <- variogram( NDVI~Non_Vegetated_or_Builtup_Lands, data.NDVI) # 30
v.NDVI3.fit = fit.variogram(v.NDVI3, vgm("Exp"));v.NDVI3.fit # 1.54
plot(v.NDVI3)

# check EVI variogram
hist( data_cleaned[, c("EVI")])
summary( data_cleaned[, c("EVI")])
hist(log(data_cleaned[, c("EVI")] + 1))
summary(log(data_cleaned[, c("EVI")] + 1))
v.EVI<- variog(coords = data_cleaned[subind, c("scaled_x", "scaled_y")], 
                data = log(data_cleaned[subind, c("EVI")] + 1), 
                uvec = (seq(0, 0.5 * d.max, length = 30))) # 30
par(mfrow=c(1,1))
varioEVI.fit <- variofit(v.EVI, cov.model="exponential")
summary(varioEVI.fit)
plot(v.EVI, xlab = "Distance (1000km)", cex = 0.1)
lines(varioEVI.fit, col='red') 
abline(h = c(varioEVI.fit$nugget, varioEVI.fit$nugget + varioEVI.fit$cov.pars[1]))
abline(v = 3 * varioEVI.fit$cov.pars[2])
variofitphiEVI <- 1 / varioEVI.fit$cov.pars[2]; variofitphiEVI # 2.5

data.EVI = data_cleaned[subind, c("scaled_x", "scaled_y", "EVI", 
                                   "Non_Vegetated_or_Builtup_Lands")]
coordinates(data.EVI) = ~scaled_x+scaled_y
v.EVI2 <- variogram( EVI~scaled_x+scaled_y+Non_Vegetated_or_Builtup_Lands, data.EVI) # 30
v.EVI2.fit = fit.variogram(v.EVI2, vgm("Exp")); 1 / v.EVI2.fit$range[2]
v.EVI2.fit
plot(v.EVI2) # phi 4.36

v.EVI3 <- variogram( EVI~Non_Vegetated_or_Builtup_Lands, data.EVI) # 30
v.EVI3.fit = fit.variogram(v.EVI3, vgm("Exp")); 1 / v.EVI3.fit$range[2] # 0.88
v.EVI3.fit
plot(v.EVI2) # phi 4.36

# check GPP variogram
hist(data_cleaned[, c("GPP")])
summary(data_cleaned[, c("GPP")])
hist(log(data_cleaned[, c("GPP")]))
summary(log(data_cleaned[, c("GPP")]))
v.GPP <- variog(coords = data_cleaned[subind, c("scaled_x", "scaled_y")], 
                data = log(data_cleaned[subind, c("GPP")]), 
                uvec = (seq(0, 0.5 * d.max, length = 30))) # 30
par(mfrow=c(1,1))
varioGPP.fit <- variofit(v.GPP, cov.model="exponential")
summary(varioGPP.fit)
plot(v.GPP, xlab = "Distance (1000km)", cex = 0.1)
lines(varioGPP.fit, col='red') 
abline(h = c(varioGPP.fit$nugget, varioGPP.fit$nugget + varioGPP.fit$cov.pars[1]))
abline(v = 3 * varioGPP.fit$cov.pars[2])
variofitphiGPP <- 1 / varioGPP.fit$cov.pars[2]; variofitphiGPP # 4.73

data.GPP = data_cleaned[subind, c("scaled_x", "scaled_y", "GPP", 
                                  "Non_Vegetated_or_Builtup_Lands")]
coordinates(data.GPP) = ~scaled_x+scaled_y
v.GPP2 <- variogram( GPP~scaled_x+scaled_y+Non_Vegetated_or_Builtup_Lands, data.GPP) # 30
v.GPP2.fit = fit.variogram(v.GPP2, vgm("Exp")); 1 / v.GPP2.fit$range[2]
plot(v.GPP2) # phi 4.36


# check PsnNet variogram
hist(data_cleaned[, c("PsnNet")])
summary(data_cleaned[, c("PsnNet")])
outliner.PsnNet <- which(data_cleaned[, c("PsnNet")]  < 0.0)
hist(log(data_cleaned[-outliner.PsnNet, c("PsnNet")] + 0.02))
summary(log(data_cleaned[, c("PsnNet")] + 0.02))

v.PsnNet <- variog(coords = data_cleaned[subind, c("scaled_x", "scaled_y")], 
                data = log(data_cleaned[subind, c("PsnNet")] + 0.02), 
                uvec = (seq(0, 0.5 * d.max, length = 30))) # 30
par(mfrow=c(1,1))
varioPsnNet.fit <- variofit(v.PsnNet, cov.model="exponential")
summary(varioPsnNet.fit)
plot(v.PsnNet, xlab = "Distance (1000km)", cex = 0.1)
lines(varioPsnNet.fit, col='red') 
abline(h = c(varioPsnNet.fit$nugget, varioPsnNet.fit$nugget + varioPsnNet.fit$cov.pars[1]))
abline(v = 3 * varioPsnNet.fit$cov.pars[2])
variofitphiPsnNet <- 1 / varioPsnNet.fit$cov.pars[2]; variofitphiPsnNet # 2.65

data.PsnNet  = data_cleaned[subind, c("scaled_x", "scaled_y", "PsnNet", 
                                  "Non_Vegetated_or_Builtup_Lands")]
coordinates(data.PsnNet ) = ~scaled_x+scaled_y
v.PsnNet2 <- variogram( PsnNet~scaled_x+scaled_y+Non_Vegetated_or_Builtup_Lands, data.PsnNet ) # 30
v.PsnNet2.fit = fit.variogram(v.PsnNet2, vgm("Exp")); 1 / v.PsnNet2.fit$range[2]
plot(v.PsnNet2) # phi 4.36



# check red_reflectance variogram
summary(data_cleaned[, c("red reflectance")])
hist(data_cleaned[, c("red reflectance")])
outliner.red = which(data_cleaned[, c("red reflectance")] > 0.40); length(outliner.red)
summary(data_cleaned[-outliner.red, c("red reflectance")])
hist(data_cleaned[-outliner.red, c("red reflectance")])
v.red <- variog(coords = data_cleaned[subind, c("scaled_x", "scaled_y")], 
                   data = data_cleaned[subind, c("red reflectance")], 
                   uvec = (seq(0, 0.5 * d.max, length = 30))) # 30
par(mfrow=c(1,1))
variored.fit <- variofit(v.red, cov.model="exponential")
summary(variored.fit)
plot(v.red, xlab = "Distance (1000km)", cex = 0.1)
lines(variored.fit, col='red') 
abline(h = c(variored.fit$nugget, variored.fit$nugget + variored.fit$cov.pars[1]))
abline(v = 3 * variored.fit$cov.pars[2])
variofitphired <- 1 / variored.fit$cov.pars[2]; variofitphired # 8.67

data.red  = data_cleaned[subind, c("scaled_x", "scaled_y", "red reflectance", 
                                      "Non_Vegetated_or_Builtup_Lands")]
names(data.red)[3] <- "red"
coordinates(data.red ) = ~scaled_x+scaled_y
v.red2 <- variogram( red~scaled_x+scaled_y+Non_Vegetated_or_Builtup_Lands, data.red ) # 30
v.red2.fit = fit.variogram(v.red2, vgm("Exp")); 1 / v.red2.fit$range[2]
plot(v.red2) # phi 7.3

v.red3 <- variogram( red~Non_Vegetated_or_Builtup_Lands, data.red ) # 30
v.red3.fit = fit.variogram(v.red3, vgm("Exp")); 1 / v.red3.fit$range[2]
plot(v.red3) # phi 4.755371


# check blue_reflectance variogram
summary(data_cleaned[, c("blue reflectance")])
hist(data_cleaned[, c("blue reflectance")])
outliner.blue<- which(data_cleaned[, c("blue reflectance")]  > 0.2)
length(outliner.blue)
summary(data_cleaned[-outliner.blue, c("blue reflectance")])
hist(data_cleaned[-outliner.blue, c("blue reflectance")])
v.blue <- variog(coords = data_cleaned[subind, c("scaled_x", "scaled_y")], 
                data = data_cleaned[subind, c("blue reflectance")], 
                uvec = (seq(0, 0.5 * d.max, length = 30))) # 30
par(mfrow=c(1,1))
varioblue.fit <- variofit(v.blue, cov.model="exponential")
summary(varioblue.fit)
plot(v.blue, xlab = "Distance (1000km)", cex = 0.1)
lines(varioblue.fit, col='red') 
abline(h = c(varioblue.fit$nugget, varioblue.fit$nugget + varioblue.fit$cov.pars[1]))
abline(v = 3 * varioblue.fit$cov.pars[2])
variofitphiblue <- 1 / varioblue.fit$cov.pars[2]; variofitphiblue # 9.18


# check LE variogram
summary(data_cleaned[, c("LE")])
hist(data_cleaned[, c("LE")])
summary(log(data_cleaned[, c("LE")]))
hist(log(data_cleaned[, c("LE")]))
v.LE <- variog(coords = data_cleaned[subind, c("scaled_x", "scaled_y")], 
                 data = log(data_cleaned[subind, c("LE")]), 
                 uvec = (seq(0, 0.8 * d.max, length = 30))) # 30
par(mfrow=c(1,1))
varioLE.fit <- variofit(v.LE, cov.model="exponential")
summary(varioLE.fit)
plot(v.LE, xlab = "Distance (1000km)", cex = 0.1)
lines(varioLE.fit, col='red') 
abline(h = c(varioLE.fit$nugget, varioLE.fit$nugget + varioLE.fit$cov.pars[1]))
abline(v = 3 * varioLE.fit$cov.pars[2])
variofitphiLE <- 1 / varioLE.fit$cov.pars[2]; variofitphiLE # 0.28

data.LE  = data_cleaned[subind, c("scaled_x", "scaled_y", "LE", 
                                   "Non_Vegetated_or_Builtup_Lands")]
coordinates(data.LE ) = ~scaled_x+scaled_y
v.LE2 <- variogram( log(LE)~scaled_x+scaled_y+Non_Vegetated_or_Builtup_Lands, data.LE) 
v.LE2.fit = fit.variogram(v.LE2, vgm("Exp")); 1 / v.LE2.fit$range[2]
plot(v.LE2) # phi 5.227206
summary(lm(LE~scaled_x+scaled_y, data = data_cleaned))

# check ET variogram
summary(data_cleaned[, c("ET")])
hist(data_cleaned[, c("ET")])
summary(log(data_cleaned[, c("ET")]))
hist(log(data_cleaned[, c("ET")]))
v.ET <- variog(coords = data_cleaned[subind, c("scaled_x", "scaled_y")], 
               data = log(data_cleaned[subind, c("ET")]), 
               uvec = (seq(0, 0.8 * d.max, length = 30))) # 30
par(mfrow=c(1,1))
varioET.fit <- variofit(v.ET, cov.model="exponential")
summary(varioET.fit)
plot(v.ET, xlab = "Distance (1000km)", cex = 0.1)
lines(varioET.fit, col='red') 
abline(h = c(varioET.fit$nugget, varioET.fit$nugget + varioET.fit$cov.pars[1]))
abline(v = 3 * varioET.fit$cov.pars[2])
variofitphiET <- 1 / varioET.fit$cov.pars[2]; variofitphiET # 0.36

data.ET  = data_cleaned[subind, c("scaled_x", "scaled_y", "ET", 
                                  "Non_Vegetated_or_Builtup_Lands")]
coordinates(data.ET) = ~scaled_x+scaled_y
v.ET2 <- variogram( log(ET)~scaled_x+scaled_y+Non_Vegetated_or_Builtup_Lands, data.ET) # 30
v.ET2.fit = fit.variogram(v.ET2, vgm("Exp")); 1 / v.ET2.fit$range[2]
v.ET2.fit
plot(v.ET2) # phi 5.242543


# check PLE variogram
summary(data_cleaned[-outliner.PLE, c("PLE")] )
hist(data_cleaned[-outliner.PLE, c("PLE")] )
outliner.PLE <- which(log(data_cleaned[, c("PLE")]) < 6.9); length(outliner.PLE)
summary(data_cleaned[-outliner.PLE, c("PLE")]*0.001)
hist(data_cleaned[-outliner.PLE, c("PLE")]*0.001)
v.PLE <- variog(coords = data_cleaned[subind, c("scaled_x", "scaled_y")], 
               data = data_cleaned[subind, c("PLE")]*0.001, 
               uvec = (seq(0, 0.8 * d.max, length = 30))) # 30
par(mfrow=c(1,1))
varioPLE.fit <- variofit(v.PLE, cov.model="exponential")
summary(varioPLE.fit)
plot(v.PLE, xlab = "Distance (1000km)", cex = 0.1)
lines(varioPLE.fit, col='red') 
abline(h = c(varioPLE.fit$nugget, varioPLE.fit$nugget + varioPLE.fit$cov.pars[1]))
abline(v = 3 * varioPLE.fit$cov.pars[2])
variofitphiPLE <- 1 / varioPLE.fit$cov.pars[2]; variofitphiPLE # 1.15

data.PLE  = data_cleaned[subind, c("scaled_x", "scaled_y", "PLE", 
                                  "Non_Vegetated_or_Builtup_Lands")]
coordinates(data.PLE) = ~scaled_x+scaled_y
v.PLE2 <- variogram( PLE~1+scaled_x+scaled_y+Non_Vegetated_or_Builtup_Lands, data.PLE) # 30
v.PLE2.fit = fit.variogram(v.PLE2, vgm("Exp")); 1 / v.PLE2.fit$range[2]
v.PLE2.fit
plot(v.PLE2) # phi 2.391995

# check PET variogram
summary(data_cleaned[, c("PET")]*0.1)
hist(data_cleaned[, c("PET")]*0.1)
outliner.PET <- which(data_cleaned[, c("PET")]*0.1 < 3.5); length(outliner.PET)
v.PET <- variog(coords = data_cleaned[subind, c("scaled_x", "scaled_y")], 
               data = data_cleaned[subind, c("PET")], 
               uvec = (seq(0, 0.8 * d.max, length = 30))) # 30
par(mfrow=c(1,1))
varioPET.fit <- variofit(v.PET, cov.model="exponential")
summary(varioPET.fit)
plot(v.PET, xlab = "Distance (1000km)", cex = 0.1)
lines(varioPET.fit, col='red') 
abline(h = c(varioPET.fit$nugget, varioPET.fit$nugget + varioPET.fit$cov.pars[1]))
abline(v = 3 * varioPET.fit$cov.pars[2])
variofitphiPET <- 1 / varioPET.fit$cov.pars[2]; variofitphiPET # 6.111299e-05

data.PET = data_cleaned[subind, c("scaled_x", "scaled_y", "PET", 
                                        "Non_Vegetated_or_Builtup_Lands")]
coordinates(data.PET) = ~scaled_x+scaled_y
v.PET2 <- variogram( PET~scaled_x+scaled_y+Non_Vegetated_or_Builtup_Lands, data.PET) # 30
v.PET2.fit = fit.variogram(v.PET2, vgm("Exp"));v.PET2.fit
plot(v.PET2)

outliner = union(union(union(union(union(outliner.blue, outliner.NDVI), 
                                   outliner.PET), outliner.PLE), 
                       outliner.PsnNet), outliner.red)
length(outliner)

data_cleaned2 = data_cleaned[-outliner, ]
data_cleaned2$NDVI = log(data_cleaned[-outliner, c("NDVI")] + 1)
data_cleaned2$EVI = log(data_cleaned[-outliner, c("EVI")] + 1)
data_cleaned2$GPP = log(data_cleaned[-outliner, c("GPP")])
data_cleaned2$PsnNet = log(data_cleaned[-outliner, c("PsnNet")] + 0.02)
data_cleaned2$LE = log(data_cleaned[-outliner, c("LE")])
data_cleaned2$ET = log(data_cleaned[-outliner, c("ET")])
data_cleaned2$PLE = (data_cleaned[-outliner, c("PLE")] *0.001)
data_cleaned2$PET = (data_cleaned[-outliner, c("PET")] *0.1)

save(data_cleaned2, file = "cleaned_data2_expanded.RData")
#load("cleaned_data2_expanded.RData")
colnames(data_cleaned2)
corr_resp = 
  cor(data_cleaned2[, c("NDVI", "EVI","GPP", "PsnNet", "red reflectance", 
                        "blue reflectance", 
                        #"view zenith angle", "sun zenith angle", 
                        #"relative azimuth angle"
                        "LE", "ET", "PLE", "PET"
  )], use = "pairwise.complete.obs")
par(mfrow = c(1,1))
#corrplot(Raw_corr , method="number")
corrplot(corr_resp , method="number")
corrplot(corr_resp , method="circle")

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

dim(data_cleaned2)
# The withheld area location index
hold_area_indx = which(data_cleaned2$x > -10400000 & data_cleaned2$x  < -10300000 &
                         data_cleaned2$y > 3800000 & data_cleaned2$y  < 3900000)
length(hold_area_indx) # 36398
# The witheld random location
set.seed(1)
N = dim(data_cleaned2)[1]
hold_ind1 <- union(sample.int(N, floor(0.1 * N)), hold_area_indx)
hold_ind2 <- union(sample.int(N, floor(0.1 * N)), hold_area_indx)
hold_ind3 <- union(sample.int(N, floor(0.1 * N)), hold_area_indx)
hold_ind4 <- union(sample.int(N, floor(0.1 * N)), hold_area_indx)
hold_ind5 <- union(sample.int(N, floor(0.1 * N)), hold_area_indx)
hold_ind6 <- union(sample.int(N, floor(0.1 * N)), hold_area_indx)
hold_ind7 <- union(sample.int(N, floor(0.1 * N)), hold_area_indx)
hold_ind8 <- union(sample.int(N, floor(0.1 * N)), hold_area_indx)
hold_ind9 <- union(sample.int(N, floor(0.1 * N)), hold_area_indx)
hold_ind10 <- union(sample.int(N, floor(0.1 * N)), hold_area_indx)
save(hold_ind1, hold_ind2, hold_ind3, hold_ind4, hold_ind5, hold_ind6, 
     hold_ind7, hold_ind8, hold_ind9, hold_ind10, hold_area_indx, 
     file = "hold_index_expanded.RData")



load("cleaned_data2_expanded.RData")
colnames(data_cleaned2)[1:21] <- 
  c("x", "y", "scaled_x", "scaled_y", "NDVI", "EVI","red_reflectance", 
    "NIR_reflectance", "blue_reflectance", "MIR_reflectance", "GPP",
    "PsnNet", "view_zenith_angle", "sun_zenith_angle", 
    "relative_azimuth_angle", "composite_day_of_the_year", 
    "LE", "ET", "PLE", "PET", "LC_Type4")

tt2 = lm(NDVI~Non_Vegetated_or_Builtup_Lands, data = data_cleaned)
summary(tt2)

tt3 = lm(lGPP~red_reflectance, data = data_cleaned2)
summary(tt3)

tt4 = lm(lET~red_reflectance, data = data_cleaned2)
summary(tt4)


# check PLE variogram
par(mfrow=c(1,1))
N = dim(data_cleaned2)[1]
set.seed(123)
subind <- sample.int(N, round(N*0.01))
data.lGPP  = data_cleaned2[subind, c("scaled_x", "scaled_y", "lNDVI", "NDVI",
                                     "red_reflectance", "NIR_reflectance",
                                   "Non_Vegetated_or_Builtup_Lands", "lGPP")]
coordinates(data.lGPP) = ~scaled_x+scaled_y
v.lGPP <- variogram( lGPP~scaled_x+scaled_y+red_reflectance+NIR_reflectance, 
                     data.lGPP) # 4.00668
v.lGPP.fit = fit.variogram(v.lGPP, vgm("Exp")); 1 / v.lGPP.fit$range[2]
v.lGPP.fit
#model      psill     range
#1   Nug 0.04558131 0.0000000
#2   Exp 0.05972745 0.2221582
plot(v.lGPP) # phi 4.501298

data.lET  = data_cleaned2[subind, c("scaled_x", "scaled_y", "lNDVI", "NDVI",
                                    "red_reflectance", "NIR_reflectance",
                                     "Non_Vegetated_or_Builtup_Lands", "lET",
                                    "ET")]
coordinates(data.lET) = ~scaled_x+scaled_y
v.lET <- variogram( lET~scaled_x+scaled_y+red_reflectance+NIR_reflectance, 
                    data.lET) # 4.318681
v.lET.fit = fit.variogram(v.lET, vgm("Exp"), fit.sills = TRUE, 
                          fit.ranges = TRUE); 1 / v.lET.fit $range[2]
v.lET.fit
#model     psill     range
#1   Nug 0.0365742 0.0000000
#2   Exp 0.1051269 0.2616677
plot(v.lET) # phi 3.821641



load("hold_index_expanded.RData")
U = intersect(intersect(intersect(intersect(intersect(intersect(intersect(
  intersect(intersect(hold_ind1, hold_ind2),  hold_ind3), hold_ind4), 
  hold_ind5), hold_ind6), hold_ind7), hold_ind8), hold_ind9), hold_ind10)
length(U)
# length 36398

# shrink the samplesize for BSLMC model try 500,000#
set.seed(321)
N = dim(data_cleaned2)[1]
small_sample_ind <- sort(sample.int(N, 1020000))
data_cleaned_small <- data_cleaned2[small_sample_ind, ]
save(data_cleaned_small, file = "data_cleaned_small_expanded.RData")
#load("data_cleaned_small_expanded.RData")
hold_area_indx2 = which(data_cleaned_small$x > -10400000 & data_cleaned_small$x  < -10300000 &
                          data_cleaned_small$y > 3800000 & data_cleaned_small$y  < 3900000)
length(hold_area_indx2) # 11900
# The witheld random location
set.seed(1)
N2 = dim(data_cleaned_small)[1]
hold_ind1 <- union(sample.int(N2, floor(0.1 * N2)), hold_area_indx2)
hold_ind2 <- union(sample.int(N2, floor(0.1 * N2)), hold_area_indx2)
hold_ind3 <- union(sample.int(N2, floor(0.1 * N2)), hold_area_indx2)
hold_ind4 <- union(sample.int(N2, floor(0.1 * N2)), hold_area_indx2)
hold_ind5 <- union(sample.int(N2, floor(0.1 * N2)), hold_area_indx2)
hold_ind6 <- union(sample.int(N2, floor(0.1 * N2)), hold_area_indx2)
hold_ind7 <- union(sample.int(N2, floor(0.1 * N2)), hold_area_indx2)
hold_ind8 <- union(sample.int(N2, floor(0.1 * N2)), hold_area_indx2)
hold_ind9 <- union(sample.int(N2, floor(0.1 * N2)), hold_area_indx2)
hold_ind10 <- union(sample.int(N2, floor(0.1 * N2)), hold_area_indx2)
save(hold_ind1, hold_ind2, hold_ind3, hold_ind4, hold_ind5, hold_ind6, 
     hold_ind7, hold_ind8, hold_ind9, hold_ind10, hold_area_indx2, 
     file = "hold_index_small_expanded.RData")
#load("hold_index_small_expanded.RData")

data_cleaned_small_copy = data_cleaned_small
data_cleaned_small_copy[hold_ind1, 1] <- NA
data_cleaned_small_copy[hold_ind2, 2] <- NA
data_cleaned_small_copy[hold_ind3, 3] <- NA
data_cleaned_small_copy[hold_ind4, 4] <- NA
data_cleaned_small_copy[hold_ind5, 5] <- NA
data_cleaned_small_copy[hold_ind6, 6] <- NA
data_cleaned_small_copy[hold_ind7, 7] <- NA
data_cleaned_small_copy[hold_ind8, 8] <- NA
data_cleaned_small_copy[hold_ind9, 9] <- NA
data_cleaned_small_copy[hold_ind10, 10] <- NA

data_response_factor = 
  data_cleaned_small_copy[, c("NDVI", "EVI","GPP", "PsnNet", 
                              "red_reflectance", "blue_reflectance", 
                              "LE", "ET", "PLE", "PET")]

colnames(data_cleaned_small_copy)
corr_resp = 
  cor(data_cleaned_small_copy[, c("NDVI", "EVI","GPP", "PsnNet", 
                                  "red_reflectance", "blue_reflectance", 
                        #"view zenith angle", "sun zenith angle", 
                        #"relative azimuth angle"
                        "LE", "ET", "PLE", "PET"
  )], use = "pairwise.complete.obs")
par(mfrow = c(1,1))
colnames(corr_resp) <- c("NDVI", "EVI","GPP", "PsnNet", 
                         "red refl", "blue refl", 
                         "LE", "ET", "PLE", "PET")
rownames(corr_resp) <- c("NDVI", "EVI","GPP", "PsnNet", 
                         "red reflectance", "blue reflectance", 
                         "LE", "ET", "PLE", "PET")
#corrplot(Raw_corr , method="number")

library(GGally)
library(corrgram)

width <- 720
height <- 720
pointsize <- 16

png(paste("BSLMC_factor_corr_plot_raw.png", sep = ""), 
    width = width, height = height, pointsize = pointsize, family = "Courier")
par(mfrow = c(1, 1))
corrgram(corr_resp, order=FALSE, lower.panel=panel.shade, gap = 0.2,
         upper.panel=panel.pie, text.panel=panel.txt, main=" ",
         col.regions = colorRampPalette(c( "darkseagreen3",
                                           "white", "cadetblue3")))
dev.off()

corrgram(corr_resp, type = "cor",order=FALSE, lower.panel=panel.shade, 
         upper.panel=panel.pie, text.panel=panel.txt, main=" ", gap = 0.2,
         col.regions = colorRampPalette(c( "white",
                                          "pink", "salmon"))) 
corrgram(corr_resp, order=FALSE, lower.panel=panel.shade, gap = 0.2,
         upper.panel=panel.pie, text.panel=panel.txt, main=" ",
         col.regions = colorRampPalette(c( "darkseagreen3",
                                           "white", "cadetblue3"))) 



# check the multivariate regression:
load("cleaned_data2_expanded.RData")
N = dim(data_cleaned2)[1]
set.seed(123)
subind <- sample.int(N, round(N*0.01))
data_test = data_cleaned2[subind, ]
names(data_test)[7] = "red_refl"
coordinates(data_test) = ~scaled_x+scaled_y
g = gstat(NULL, "NDVI_v", NDVI~Non_Vegetated_or_Builtup_Lands+scaled_x+scaled_y, 
          data_test)
g = gstat(g, "red_refl_v", 
          red_refl~Non_Vegetated_or_Builtup_Lands+scaled_x+scaled_y, data_test)
t <- proc.time()
v = variogram(g)
proc.time() - t
plot(variogram(g, cutoff=0.7, width=0.03, map=TRUE), 
     main = "(cross) semivariance maps")

plot(variogram(g, cutoff=0.7, width=0.03, map=TRUE), np=TRUE,
     main = "number of point pairs")

meuse.g <- gstat(id = "NDVI", 
                 formula = NDVI ~ Non_Vegetated_or_Builtup_Lands+scaled_x+scaled_y, 
                 data = data_test)
meuse.g <- gstat(meuse.g, "red_refl", 
                 red_refl~Non_Vegetated_or_Builtup_Lands+scaled_x+scaled_y,
                 data_test)
meuse.g <- gstat(meuse.g, model = vgm(2, "Exp", 6, 0.1), fill.all = TRUE)
x <- variogram(meuse.g, cutoff = 0.7)
meuse.fit = fit.lmc(x, meuse.g)
out = gstat.cv(meuse.fit, nmax = 40, nfold = 5)
summary(out)

X = as.matrix(cbind(rep(1.0, length(subind)), 
          data_cleaned2[subind, c("scaled_x", "scaled_y")],
          data_cleaned2[subind, c("scaled_x", "scaled_y")]^2, 
          data_cleaned2[subind, "scaled_x"] * 
            data_cleaned2[subind, "scaled_y"]))
X = as.matrix(cbind(rep(1.0, length(subind)), 
                    data_cleaned2[subind, c("Non_Vegetated_or_Builtup_Lands")]))
Y = as.matrix(data_cleaned2[subind, c("NDVI", "red reflectance")])
summary(Y); cov(Y)
resid = Y - X %*% solve(crossprod(X))%*%(crossprod(X, Y))
summary(resid); cov(resid)

v.resid1 <- variog(coords = data_cleaned2[subind, c("scaled_x", "scaled_y")], 
                data = resid[, 1], 
                uvec = (seq(0, 1.0, length = 30))) # 30
par(mfrow=c(1,1))
vario.resid1.fit <- variofit(v.resid1, cov.model="exponential")
summary(vario.resid1.fit)
plot(v.resid1, xlab = "Distance (1000km)", cex = 0.1)
lines(vario.resid1.fit, col='red') 
abline(h = c(vario.resid1.fit$nugget, vario.resid1.fit$nugget + 
               vario.resid1.fit$cov.pars[1]))
abline(v = 3 * vario.resid1.fit$cov.pars[2])
variofitphi.resid1 <- 1 / vario.resid1.fit$cov.pars[2]; variofitphi.resid1 

v.resid2 <- variog(coords = data_cleaned2[subind, c("scaled_x", "scaled_y")], 
                   data = resid[, 2], 
                   uvec = (seq(0, 0.8, length = 50))) # 30
par(mfrow=c(1,1))
vario.resid2.fit <- variofit(v.resid2, cov.model="exponential")
summary(vario.resid2.fit)
plot(v.resid2, xlab = "Distance (1000km)", cex = 0.1)
lines(vario.resid2.fit, col='red') 
abline(h = c(vario.resid2.fit$nugget, vario.resid2.fit$nugget + 
               vario.resid2.fit$cov.pars[1]))
abline(v = 3 * vario.resid2.fit$cov.pars[2])
variofitphi.resid2 <- 1 / vario.resid2.fit$cov.pars[2]; variofitphi.resid2 


v.NDVI <- variog(coords = data_cleaned2[subind, c("scaled_x", "scaled_y")], 
                   data = data_cleaned2[subind, "NDVI"], 
                   uvec = (seq(0, 0.8, length = 50))) # 30
par(mfrow=c(1,1))
vario.NDVI.fit <- variofit(v.NDVI, cov.model="exponential")
summary(vario.NDVI.fit)
plot(v.NDVI, xlab = "Distance (1000km)", cex = 0.1)
lines(vario.NDVI.fit, col='red') 
abline(h = c(vario.NDVI.fit$nugget, vario.NDVI.fit$nugget + 
               vario.NDVI.fit$cov.pars[1]))
abline(v = 3 * vario.NDVI.fit$cov.pars[2])
variofitphi.NDVI <- 1 / vario.NDVI.fit$cov.pars[2]; variofitphi.NDVI 

v.red <- variog(coords = data_cleaned2[subind, c("scaled_x", "scaled_y")], 
                 data = data_cleaned2[subind, "red reflectance"], 
                 uvec = (seq(0, 1.0, length = 50))) # 30
par(mfrow=c(1,1))
vario.red.fit <- variofit(v.red, cov.model="exponential")
summary(vario.red.fit)
plot(v.red, xlab = "Distance (1000km)", cex = 0.1)
lines(vario.red.fit, col='red') 
abline(h = c(vario.red.fit$nugget, vario.red.fit$nugget + 
               vario.red.fit$cov.pars[1]))
abline(v = 3 * vario.red.fit$cov.pars[2])
variofitphi.red <- 1 / vario.red.fit$cov.pars[2]; variofitphi.red 



