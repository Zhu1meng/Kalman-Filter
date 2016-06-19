##########################################
##                                      ##
##      KALMAN FILTER SIMULATION        ##
##      MP Roeling, miniproject 1       ##
##      Oxford University               ##   
##      April 2016                      ##
##                                      ##
##########################################

setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/")

################ RESIDUALS ###################
# position
residuals = read.table("testresidualsposition.txt")
residuals$scale = scale(residuals$V1)
yrange <- c(min(residuals$scale),max(residuals$scale))
# percentage within 95% CI 
table(residuals$scale > -1.96 & residuals$scale < 1.96)
perc.withinCI = round(table(residuals$scale > -1.96 & residuals$scale < 1.96)[2]/nrow(residuals)*100,2)
# plot
xrange = c(1,900)
yrange = c(-5,5)
plot(xrange, yrange, type = "n", xlab = "Time (sec)", ylab = "Position")
lines(residuals$scale)
abline(h = 1.96, col = "red")
abline(h = -1.96, col = "red")
title("Residuals and 95% Confidence Interval (red lines)")
text(200, -3, paste(perc.withinCI, "% within 95% CI", sep=""))

# velocity
residuals = read.table("testresidualsvelocity.txt")
residuals$scale = scale(residuals$V1)
yrange <- c(min(residuals$scale),max(residuals$scale))
# percentage within 95% CI 
table(residuals$scale > -1.96 & residuals$scale < 1.96)
perc.withinCI = round(table(residuals$scale > -1.96 & residuals$scale < 1.96)[2]/nrow(residuals)*100,2)
# plot
xrange = c(1,900)
yrange = c(-5,5)
plot(xrange, yrange, type = "n", xlab = "Time (sec)", ylab = "Velocity")
lines(residuals$scale)
abline(h = 1.96, col = "red")
abline(h = -1.96, col = "red")
title("Residuals and 95% Confidence Interval (red lines)")
text(200, -3, paste(perc.withinCI, "% within 95% CI", sep=""))

############ INNOVATION FACTORS #############
# position
innovationfactors = read.table("normalized_if_pos.txt")
hist(innovationfactors$V1)
plot(innovationfactors$V1)
mean(innovationfactors$V1)
sd(innovationfactors$V1)
# velocity
innovationfactors = read.table("normalized_if_vel.txt")
hist(innovationfactors$V1)
plot(innovationfactors$V1)
mean(innovationfactors$V1)
sd(innovationfactors$V1)

############# KALMAN GAIN ###################
kalmangain = read.table("testK.txt") 

######## ESIMATED ACCELATION SPEED Ax ##########
testAx = read.table("testAx.txt")
hist(testAx$V1)
mean(testAx$V1)
sd(testAx$V1)

############# PERFORMANCE ###################
rm(list=ls(all=TRUE))
setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/performance/")

Jk_position0 = read.table("Jk_position_a0.txt")
colnames(Jk_position0) = c("JK", "estimated measurement", "estimated state", "noisy real measurement")
Jk_position1 = read.table("Jk_position_a1.txt")
colnames(Jk_position1) = c("JK", "estimated measurement", "estimated state", "noisy real measurement")
Jk_position2 = read.table("Jk_position_a2.txt")
colnames(Jk_position2) = c("JK", "estimated measurement", "estimated state", "noisy real measurement")
Jk_position3 = read.table("Jk_position_a3.txt")
colnames(Jk_position3) = c("JK", "estimated measurement", "estimated state", "noisy real measurement")
Jk_position4 = read.table("Jk_position_a4.txt")
colnames(Jk_position4) = c("JK", "estimated measurement", "estimated state", "noisy real measurement")

# calculate ratios again
Jk_position$ratio1 = abs(Jk_position$`estimated measurement` - Jk_position$`estimated state`)
Jk_position$ratio2 = abs(Jk_position$`estimated measurement` - Jk_position$`noisy real measurement`)
# plot
# Jk_position0$JK = Jk_position0$JK + 0.02

plot(Jk_position0$JK, ylim = c(0.30,0.35),col = "black", main = "Distribution of performance parameter J", xlab = "Time", ylab = "J")
lines(Jk_position1$JK, ylim = c(0.30,0.35),col = "blue")
lines(Jk_position2$JK, ylim = c(0.30,0.35),col = "red")
lines(Jk_position3$JK, ylim = c(0.30,0.35),col = "green")
lines(Jk_position4$JK, ylim = c(0.30,0.35),col = "orange")

##############################################
############### NON LINEAR ###################
##############################################
setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/3dperformance/attack0/")

################ RESIDUALS ###################
# longitude
residuals = read.table("testresidualslongitude.txt")
residuals$scale = scale(residuals$V1)
yrange <- c(min(residuals$scale),max(residuals$scale))
# percentage within 95% CI 
table(residuals$scale > -1.96 & residuals$scale < 1.96)
perc.withinCI = round(table(residuals$scale > -1.96 & residuals$scale < 1.96)[2]/nrow(residuals)*100,2)
# plot
xrange = c(1,900)
yrange = c(-5,5)
plot(xrange, yrange, type = "n", xlab = "Time (sec)", ylab = "Position")
lines(residuals$scale)
abline(h = 1.96, col = "red")
abline(h = -1.96, col = "red")
title("Residuals and 95% Confidence Interval (red lines)")
text(200, -3, paste(perc.withinCI, "% within 95% CI", sep=""))

# latitide
residuals = read.table("testresidualslatitude.txt")
residuals$scale = scale(residuals$V1)
yrange <- c(min(residuals$scale),max(residuals$scale))
# percentage within 95% CI 
table(residuals$scale > -1.96 & residuals$scale < 1.96)
perc.withinCI = round(table(residuals$scale > -1.96 & residuals$scale < 1.96)[2]/nrow(residuals)*100,2)
# plot
xrange = c(1,900)
yrange = c(-5,5)
plot(xrange, yrange, type = "n", xlab = "Time (sec)", ylab = "Velocity")
lines(residuals$scale)
abline(h = 1.96, col = "red")
abline(h = -1.96, col = "red")
title("Residuals and 95% Confidence Interval (red lines)")
text(200, -3, paste(perc.withinCI, "% within 95% CI", sep=""))

# altitide
residuals = read.table("testresidualsaltitude.txt")
residuals$scale = scale(residuals$V1)
yrange <- c(min(residuals$scale),max(residuals$scale))
# percentage within 95% CI 
table(residuals$scale > -1.96 & residuals$scale < 1.96)
perc.withinCI = round(table(residuals$scale > -1.96 & residuals$scale < 1.96)[2]/nrow(residuals)*100,2)
# plot
xrange = c(1,900)
yrange = c(-5,5)
plot(xrange, yrange, type = "n", xlab = "Time (sec)", ylab = "Velocity")
lines(residuals$scale)
abline(h = 1.96, col = "red")
abline(h = -1.96, col = "red")
title("Residuals and 95% Confidence Interval (red lines)")
text(200, -3, paste(perc.withinCI, "% within 95% CI", sep=""))


# non linear model performance longitude
setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/3dperformance/attack0/")
Jk_position0 = read.table("Jk_longitude.txt")
colnames(Jk_position0) = c("JK", "estimated measurement", "estimated state", "noisy real measurement")
setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/3dperformance/attack1/")
Jk_position1 = read.table("Jk_longitude.txt")
colnames(Jk_position1) = c("JK", "estimated measurement", "estimated state", "noisy real measurement")
setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/3dperformance/attack2/")
Jk_position2 = read.table("Jk_longitude.txt")
colnames(Jk_position2) = c("JK", "estimated measurement", "estimated state", "noisy real measurement")
setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/3dperformance/attack3/")
Jk_position3 = read.table("Jk_longitude.txt")
colnames(Jk_position3) = c("JK", "estimated measurement", "estimated state", "noisy real measurement")
setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/3dperformance/attack4/")
Jk_position4 = read.table("Jk_longitude.txt")
colnames(Jk_position4) = c("JK", "estimated measurement", "estimated state", "noisy real measurement")

xrange <- c(0,1050)
yrange <- c(3623.9679,3623.968)
# set up the plot
plot(xrange, yrange, type = "n", xlab = "Time (sec)", ylab = "J", main = "Distribution of performance parameter J") 
lines(Jk_position0[c(2:nrow(Jk_position0)),]$JK, col = "black")
lines(Jk_position1[c(2:nrow(Jk_position1)),]$JK, col = "orange")
lines(Jk_position2[c(2:nrow(Jk_position2)),]$JK, col = "blue")
lines(Jk_position3[c(2:nrow(Jk_position3)),]$JK, col = "red")
lines(Jk_position4[c(2:nrow(Jk_position4)),]$JK, col = "green")

# non linear model performance latitude
setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/3dperformance/attack0/")
Jk_latitude0 = read.table("Jk_latitude.txt")
colnames(Jk_latitude0) = c("JK", "estimated measurement", "estimated state", "noisy real measurement")
setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/3dperformance/attack1/")
Jk_latitude1 = read.table("Jk_latitude.txt")
colnames(Jk_latitude1) = c("JK", "estimated measurement", "estimated state", "noisy real measurement")
setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/3dperformance/attack2/")
Jk_latitude2 = read.table("Jk_latitude.txt")
colnames(Jk_latitude2) = c("JK", "estimated measurement", "estimated state", "noisy real measurement")
setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/3dperformance/attack3/")
Jk_latitude3 = read.table("Jk_latitude.txt")
colnames(Jk_latitude3) = c("JK", "estimated measurement", "estimated state", "noisy real measurement")
setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/3dperformance/attack4/")
Jk_latitude4 = read.table("Jk_latitude.txt")
colnames(Jk_latitude4) = c("JK", "estimated measurement", "estimated state", "noisy real measurement")

xrange <- c(0,1050)
yrange <- c(14492.68,14492.69)
plot(xrange, yrange, type = "n", xlab = "Time (sec)", ylab = "J", main = "Distribution of performance parameter J") 
lines(Jk_latitude0[c(2:nrow(Jk_latitude0)),]$JK, col = "black")
lines(Jk_latitude1[c(2:nrow(Jk_latitude1)),]$JK, col = "orange")
lines(Jk_latitude2[c(2:nrow(Jk_latitude2)),]$JK, col = "blue")
lines(Jk_latitude3[c(2:nrow(Jk_latitude3)),]$JK, col = "red")
lines(Jk_latitude4[c(2:nrow(Jk_latitude4)),]$JK, col = "green")

############# END ###############