rm(list=ls(all=TRUE))
##########################################
##                                      ##
##      KALMAN FILTER SIMULATION        ##
##      MP Roeling, miniproject 1       ##
##      Oxford University               ##   
##      April 2016                      ##
##                                      ##
##########################################

library(geosphere)
# http://www.movable-type.co.uk/scripts/latlong.html
# https://en.wikipedia.org/wiki/Great-circle_distance
# https://en.wikipedia.org/wiki/Geographic_coordinate_system

setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/")
alldata = read.table("kalman.txt", sep="|", h=T)
alldata = alldata[,-c(1,ncol(alldata))]

# test coordinates
coor1 = alldata[1,c(5,6)]
coor2 = alldata[20,c(5,6)]

# online tool http://www.tytai.com/gmap/distance/
# Approximate Distance: 34,933.24 meters

# calculated with Haversine formula 
distHaversine(coor1, coor2, r=6378137)
# output = 34972.44 meters

# calculated with Vincenty (ellipsoid) method
distVincentyEllipsoid(coor1, coor2, a=6378137, b=6356752.314140, f=1/298.257222101)

" Ellipsoid used from European Terrestrial Reference System 1989
  The European Terrestrial Reference System 1989 (ETRS89) is an ECEF (Earth-Centered, Earth-Fixed) 
  geodetic Cartesian reference frame, in which the Eurasian Plate as a whole is static. 
  The coordinates and maps in Europe based on ETRS89 are not subject to change due to the continental drift.
  https://en.wikipedia.org/wiki/European_Terrestrial_Reference_System_1989
  a = Equatorial axis of ellipsoid = 6378137.000
  b = Polar axis of ellipsoid = 6356752.314140
  f = Inverse flattening of ellipsoid = 298.257222101 "

####################################################################
# select first flight 
data = alldata[alldata$time<1444740000,]

data$distance = -1
data$delta.time = 0
for (i in 1:(nrow(data)-1)){
  # distance 
  coor1 = data[i,c(5,6)]
  coor2 = data[i+1,c(5,6)]
  data[i+1,]$distance = distVincentyEllipsoid(coor1, coor2, 
                                              a=6378137, 
                                              b=6356752.314140, 
                                              f=1/298.257222101)
  print(i)
}

# many non informative rows
table(data$distance==0)

# altitude is integer, change to double: typeof(data$altitude)
data$altitude_doub = as.double(paste(data$altitude))

# remove non informative data
data_info = data[data$distance != 0,]
data_info = data_info[,c(1,5,6,16,11,14,15)]

data_info$delta.lon = 0
data_info$delta.lat = 0
for (i in 1:(nrow(data_info)-1)){
  # time 
  data_info[i+1,]$delta.time = data_info[i+1,1] - data_info[i,1]
  # longitude lattitude 
  data_info[i+1,]$delta.lon = data_info[i+1,2] - data_info[i,2]
  data_info[i+1,]$delta.lat = data_info[i+1,3] - data_info[i,3]
  print(i)
}

# standardize distance / time
data_info$distime = data_info$distance / data_info$delta.time

# plot with http://www.gpsvisualizer.com/
plot_gps_data = cbind("Aircraft","OHY925",data_info$latitude, data_info$longitude, 1)
write.table(plot_gps_data, "plot_gps_data.txt", quote = F, sep = ",", row.names = F, col.names = F)

# linear model data
linear.model.data = data_info[c(2:500),]
# select raw data
radartest = cbind("Aircraft","OHY925",linear.model.data$latitude, linear.model.data$longitude, 1)  
# change last gps coordinate to radar station
radartest[nrow(radartest),5] <- 2 
write.table(linear.model.data, "linearmodeldata.txt", quote = F, sep = ",", row.names = F, col.names = T)
write.table(radartest, "linearradartest.txt", quote = F, sep = ",", row.names = F, col.names = F)

# non linear model 
nonlinear.model.data = data_info[c(501:nrow(data_info)),]
# select raw data
radartest = cbind("Aircraft","OHY925",nonlinear.model.data$latitude, nonlinear.model.data$longitude, 1)   
# add radar station
radartest = rbind(radartest,c("Radar", "AOCS NM", 52.233904, 5.742352, 2))
write.table(nonlinear.model.data, "nonlinearmodeldata.txt", quote = F, sep = ",", row.names = F, col.names = T)
write.table(radartest, "nonlinearradartest.txt", quote = F, sep = ",", row.names = F, col.names = F)

###############################################################################################
rm(list=ls(all=TRUE))
data = read.table("linearmodeldata.txt" ,sep = ",", h = T)

# set radar station (last data point here)
radar = c(10.24911, 51.31821)
data$radar.distance = 0
data$velocity = 0

for(i in 1:nrow(data)){
  airplane = data[i,c(2,3)]
  # for more info see KF_1.R
  # calculated with Vincenty (ellipsoid) method 
  # distance in meters
  radar_distance = distVincentyEllipsoid(radar, airplane, a=6378137, b=6356752.314140, f=1/298.257222101)
  data[i,]$radar.distance = radar_distance
  if(i > 1){
    # calculate velocity
    distance_step = radar_distance - old_radar_distance
    data[i,]$velocity = distance_step
  }
  old_radar_distance = radar_distance
}

# sometimes the data step is > 1, because multiple timesteps resulted in the same estimate  
# this loop takes the score from timepoint unique timepoint 1 and 2, and depending on the amount of identical 
# steps in between it takes the average

for(i in 2:nrow(data)){
  delta.time = data[i,]$delta.time
  mean.radar.distance = (data[i,]$radar.distance - data[i-1,]$radar.distance) / delta.time
  
  test = data.frame(c(1:delta.time),data[i-1,]$radar.distance)
  test[c(1:delta.time),2] <- test[c(1:delta.time),2] + ((1:delta.time) * mean.radar.distance)
  test[delta.time,2] <- data[i,]$radar.distance
  test = cbind(i, test)
  write.table(test, "output.txt", sep="\t", col.names=F, row.names=F, append=T)
}

# end

########################### NON LINEAR ###########################################
rm(list=ls(all=TRUE))
setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/")
data = read.table("nonlinearmodeldata.txt" ,sep = ",", h = T)

# set radar station
radar = c(5.742352, 52.233904)
data$radar.distance = 0
data$velocity = 0

for(i in 1:nrow(data)){
  airplane = data[i,c(2,3)]
  # for more info see KF_1.R
  # calculated with Vincenty (ellipsoid) method 
  # distance in meters
  radar_distance = distVincentyEllipsoid(radar, airplane, a=6378137, b=6356752.314140, f=1/298.257222101)
  data[i,]$radar.distance = radar_distance
  if(i > 1){
    # calculate velocity
    distance_step = radar_distance - old_radar_distance
    data[i,]$velocity = distance_step
  }
  old_radar_distance = radar_distance
}

# sometimes the data step is > 1, because multiple timesteps resulted in the same estimate  
# this loop takes the score from timepoint unique timepoint 1 and 2, and depending on the amount of identical 
# steps in between it takes the average

for(i in 2:nrow(data)){
  delta.time = data[i,]$delta.time
  mean.radar.distance = (data[i,]$radar.distance - data[i-1,]$radar.distance) / delta.time
  
  test = data.frame(c(1:delta.time),data[i-1,]$radar.distance)
  test[c(1:delta.time),2] <- test[c(1:delta.time),2] + ((1:delta.time) * mean.radar.distance)
  test[delta.time,2] <- data[i,]$radar.distance
  test = cbind(i, test, data[i,]$altitude_doub)
  write.table(test, "nonlinearoutput.txt", sep="\t", col.names=F, row.names=F, append=T)
}

# end
########################### GPS ###########################################
# simulate noise
# for every GPS position, we add up to 100 meters noise in either direction.
# First, we calculate maximum GPS deviations with the vincenti method, so that the distance 
# between GPS coordinates cannot exceed 100 meters. Then, within that area we randomly select a datapoint.
rm(list=ls(all=TRUE))
setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/")
data = read.table("nonlinearmodeldata.txt" ,sep = ",", h = T)
# distance.threshold is the amount of noise in meters.
data$noiselong = 0
data$noiselat = 0 
distance.threshold = 100

pb <- txtProgressBar(min = 0, max = nrow(data), style = 3)

for(i in 1:nrow(data)){
  # Longitude
  # positive deviation 
  GPS.long1 = data[i,]$longitude
  GPS.lat = data[i,]$latitude
  GPS.original = cbind(GPS.long1,GPS.lat)
  distance = 0
  while (distance < distance.threshold){
    GPS.long1 = GPS.long1 + 0.000001
    GPS.simulate = cbind(GPS.long1,GPS.lat)
    distance = distVincentyEllipsoid(GPS.original, GPS.simulate, a=6378137, b=6356752.314140, f=1/298.257222101)
    # print(distance)
  }
  # negative deviation 
  GPS.long2 = data[i,]$longitude
  GPS.lat = data[i,]$latitude
  GPS.original = cbind(GPS.long2,GPS.lat)
  distance = 0
  while (distance < distance.threshold){
    GPS.long2 = GPS.long2 - 0.000001
    GPS.simulate = cbind(GPS.long2,GPS.lat)
    distance = distVincentyEllipsoid(GPS.original, GPS.simulate, a=6378137, b=6356752.314140, f=1/298.257222101)
    #print(distance)
  }
  
  # Lattitude 
  # positive deviation 
  GPS.long = data[i,]$longitude
  GPS.lat1 = data[i,]$latitude
  GPS.original = cbind(GPS.long,GPS.lat1)
  distance = 0
  while (distance < distance.threshold){
    GPS.lat1 = GPS.lat1 + 0.000001
    GPS.simulate = cbind(GPS.long,GPS.lat1)
    distance = distVincentyEllipsoid(GPS.original, GPS.simulate, a=6378137, b=6356752.314140, f=1/298.257222101)
    # print(distance)
  }
                     
  # negative deviation 
  GPS.long = data[i,]$longitude
  GPS.lat2 = data[i,]$latitude
  GPS.original = cbind(GPS.long,GPS.lat2)
  distance = 0
  while (distance < distance.threshold){
    GPS.lat2 = GPS.lat2 - 0.000001
    GPS.simulate = cbind(GPS.long,GPS.lat2)
    distance = distVincentyEllipsoid(GPS.original, GPS.simulate, a=6378137, b=6356752.314140, f=1/298.257222101)
    # print(distance)
  }
  
  # take longitude
  sample.space.long = seq(GPS.long1, GPS.long2, -0.00001)
  sample.space.lat = seq(GPS.lat1, GPS.lat2, -0.00001)
  
  data[i,]$noiselong = sample(sample.space.long, 1) 
  data[i,]$noiselat = sample(sample.space.lat, 1)   
  setTxtProgressBar(pb, i)
}                      
close(pb)
                      
write.table(data, "data_with_noise.txt", col.names = T, row.names = F, quote = F, sep = "\t")                      
       
library(plot3D)
library(rgl)
plot3d(data$noiselong, data$noiselat, data$altitude_doub)                      
                      
                      
