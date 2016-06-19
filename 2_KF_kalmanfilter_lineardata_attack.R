rm(list=ls(all=TRUE))
##########################################
##                                      ##
##      KALMAN FILTER SIMULATION        ##
##      MP Roeling, miniproject 1       ##
##      Oxford University               ##   
##      April 2016                      ##
##                                      ##
##########################################

setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/")

##################################################################################
output = read.table("output.txt", sep="\t")
  colnames(output) = c("iteration", "timestep", "radar.distance")
# Rplot_linear_distance_time  
plot(output$radar.distance, xlab = "Time (sec)", ylab = "Distance (meters)", main = "Raw data distance to radar over time")

# calculate velocity
output$velocity = 0
for(j in 2:nrow(output)){
  output[j,]$velocity = output[j,3] - output[j-1,3] 
}
  
startdistance = output[1,]$radar.distance
# QC step remove first row of data
data = output[c(2:nrow(output)),c(3,4)]

#### Kalman filter process ####
  
# H = transformation matrix, R = error in the observation, Y = observation matrix
# 1. Xkp = AXk-1 + BUk + Wk  -> predicted state
# 2. get the starting values for Pkp
# 3. Pkp = APk-1AT + Qk -> process cov matrix
# 4. K = Pkp HT / KPkpHT + R
# 5. Yk = CYkm + Zm
# 6. Xk = Xkp + K[Y-HXkp]
# 7. Pk = (I-KH)Pkp & Xk

nvar = 2
number = nrow(data)
graphdata=matrix(0,(nrow(data)-1),nvar*2+1)
if(exists("i")){rm(i)}
# eigenlijk begin je pas data op te halen bij de 2e regel bij observed state schatting (stap 5) vs predicted
# die eerste regel is al gegeven.
  X.0 = startdistance # initial position
  Vx.0 = mean(data$velocity) # initial velocity
  Vy.0 = 120 ; Y.0 = 3000
  # initial conditions
  # Ax = 0 # acceleration speed up / down 10 m/s, defined in the loop
  Delta.T = 1 # timestep 1 second (unclear if the timesteps are seconds but we assume every step is a second
              # otherwise the dynamical model for velocity would not be accurate)
  Vx=mean(data$velocity)# initial velocity
  Delta.X = 500 # uncertainty in measurement
  # process error in process covariance matrix, is updated at the end of iteration 1.
  Delta.Px = 1000 # 1 km start error
  Delta.Pvx = 300 # 300 m/s start error
  # observation errors
  Delta.X = 1000 # 1000 meter sd
  Delta.Vx = 250 # 250 m/s sd

  A = diag(nvar)
  C = diag(nvar)
  I = diag(nvar)
  H = diag(nvar) 
  var_velocity = var(data$velocity)
  Qr = 500 # error in the process of calculating the process cov matrix, between 0.5 and 1 Ax 
  Zk = 0 # observation errors in mechanism (eg. electronic delays)
  # simulate white noise
  white_noise_gaussian = rnorm(nrow(data), 0, sd(data[,2]))

set.seed(100)
for (i in 1:nrow(data)){ 
# 1. Xkp = AXk-1 + BUk + Wk -> predicted state
  # A allows to update position and velocity
  A[row(A) == (col(A) - 1)]=Delta.T
  if(i==1){
    X = matrix(c(X.0,Vx.0),2,1) # first iteration take start values
  } else {
    X = Xk # second or more iterations take previous estimate
  }
  # B translates acceleration (Ax) into adjustment to position and velocity,
  B = matrix(c(1/2*Delta.T^2,Delta.T),nvar,1)
  # acceleration is not constant (decrease and increase, get estimate from velocity
  # which itself is an estimate from distance over time
  Wk = matrix(c(sample(white_noise_gaussian, 1),
                sample(white_noise_gaussian, 1)),nvar,1)
  if(i < nrow(data)){
    Ax = data[i+1,2]-data[i,2]
    Xkp = A%*%X + B%*%Ax + Wk
  } else {
    Ax = 0
    Xkp = A%*%X + B%*%Ax + Wk
  }
  # I let the loop continue to this point to get the extra prediction Xkp
  if(i==nrow(data)){break} 

# 2. Pkp.i = APk-1AT + Qk -> initial process cov matrix
  if(i==1){ # first iteration take start values
    Pkp.i = matrix(c(Delta.Px^2, Delta.Px*Delta.Pvx,
                   Delta.Px*Delta.Pvx, Delta.Pvx^2),2,2)
    # setting the cross terms to zero
    Pkp.i[row(Pkp.i) == (col(Pkp.i) - 1)]=0 # upper offdiagonal
    Pkp.i[row(Pkp.i) == (col(Pkp.i) + 1)]=0 # lower offdiagonal
  } else {
    Pkp.i = Pk # second or more iterations take previous estimate
    Qr = matrix(c(Delta.T^4/4,Delta.T^3/2,
                  Delta.T^3/2,Delta.T^2),2,2) # formula see page 234 of Kalman Filter Python book
    Qr = Qr * var_velocity
  }
  
# 3. Pkp = APkp.iAT + QR >- Predicted / adjusted process cov matrix
  Pkp = A%*%Pkp.i%*%t(A) + Qr
  # setting the cross terms to zero
  Pkp[row(Pkp) == (col(Pkp) - 1)]=0 # upper offdiagonal
  Pkp[row(Pkp) == (col(Pkp) + 1)]=0 # lower offdiagonal

# 4. K = Pkp HT / KPkpHT + R -> Calculate the Kalman Gain
  R = matrix(c(Delta.X^2,0,0,Delta.Vx^2),2,2)
  # Pkp%*%t(H) <- upper term
  # H%*%Pkp*t(H)+R <- lower term
  # dividing equals multiplication with inverse (solve function)
  K = Pkp%*%t(H)*solve(H%*%Pkp*t(H)+R)
  rho2 = H%*%Pkp*t(H)+R

# 5. Yk = CYkm + Zm -> represents observed state
  Ykm = matrix(c(data[i+1,1],data[i+1,2]),nvar,1)
  Yk = C%*%Ykm+Zk
  
# 6. Xk = Xkp + K[Y-HXkp] -> calculate current state
  Xk = Xkp + K%*%(Yk-H%*%Xkp)
  innovationfactor = Yk-H%*%Xkp
  normalized_innovationfactor = innovationfactor / matrix(c(sqrt(rho2)[1],sqrt(rho2)[4]),2,1)
  
# 7. Pk = (I-KH)Pkp & Xk -> update process cov matrix
  Pk = (I - K*H) %*% Pkp
  
  graphdata[i,1]=i
  graphdata[i,2]=Xkp[1] # predicted position
  graphdata[i,3]=Xkp[2] # predicted velocity
  graphdata[i,4]=Xk[1] # Kalman filter state position
  graphdata[i,5]=Xk[2] # Kalman filter state velocity
  write.table(K[1], "testK.txt", append = T, col.names = F, row.names = F)
  write.table(Pkp[1], "testPkp.txt", append = T, col.names = F, row.names = F)
  write.table(innovationfactor[1], "testresidualsposition.txt", append = T, col.names = F, row.names = F)
  write.table(innovationfactor[2], "testresidualsvelocity.txt", append = T, col.names = F, row.names = F)
  write.table(Ax, "testAx.txt", append = T, col.names = F, row.names = F)
  write.table(normalized_innovationfactor[1], "normalized_if_pos.txt", append = T, col.names = F, row.names = F)
  write.table(normalized_innovationfactor[2], "normalized_if_vel.txt", append = T, col.names = F, row.names = F)
  ow <- options(warn = 2) 
}

## GRAPH ## 
# there will always be an extra prediction, which is unreliable because it 
# does not take the Nth measurement into account at Ykm, but takes the previous N-1 Ykm

data.select = cbind(seq(1,nrow(data),1),data,1)
  colnames(data.select)<-c("Time","Position","Velocity","Group")
predicted = cbind(seq(1,(nrow(graphdata)+2),1),rbind(c(X.0,Vx.0),graphdata[1:nrow(graphdata),c(2,3)],c(Xkp[1],Xkp[2])),2) # predicted
  colnames(predicted)<-c("Time","Position","Velocity","Group")
kalman = cbind(seq(1,(nrow(graphdata)+1),1),rbind(c(X.0,Vx.0),graphdata[1:nrow(graphdata),c(4,5)]),3) # kalman
  colnames(kalman)<-c("Time","Position","Velocity","Group")


graphdata = rbind(data.select, predicted, kalman)
# zoom
graphdata = rbind(data.select[c(1:100),], predicted[c(1:100),], kalman[c(1:100),])

# convert factor to numeric for convenience 
ngroups=nvar+1
# get the range for the x and y axis 
xrange <- c(min(graphdata[,1]),max(graphdata[,1]))
yrange <- c(min(graphdata[,2]),max(graphdata[,2]))
# set up the plot
plot(xrange, yrange, type = "n", xlab = "Time (sec)", ylab = "Position") 
colors = rainbow(ngroups) 
linetype = c(2,3,1) 
group1 = graphdata[graphdata[,ncol(graphdata)] == 1,]
group2 = graphdata[graphdata[,ncol(graphdata)] == 2,]
group3 = graphdata[graphdata[,ncol(graphdata)] == 3,]
lines(group1[,1],group1[,2], type="b", lwd = 1.5, col = colors[1], lty = linetype[1])
lines(group2[,1],group2[,2], type="b", lwd = 1.5, col = colors[2], lty = linetype[2])
lines(group3[,1],group3[,2], type="b", lwd = 1.5, col = colors[3], lty = linetype[3])
# add a title and subtitle 
title("Observed (red), predicted (green), and Kalman (blue)")
# add a legend 
legend(1, 27000, 1:ngroups, cex = 0.8, col=colors, lty = linetype, title = "Groups")

############################################################
# position
  residuals = read.table("testresidualsposition.txt")
  residuals$scale = scale(residuals$V1)
  yrange <- c(min(residuals$scale),max(residuals$scale))
  # percentage within 95% CI 
  table(residuals$scale > -1.96 & residuals$scale < 1.96)
  perc.withinCI = round(table(residuals$scale > -1.96 & residuals$scale < 1.96)[1]/nrow(residuals)*100,2)
  
  plot(xrange, yrange, type = "n", xlab = "Time (sec)", ylab = "Position")
    lines(residuals$scale)
    abline(h = 1.96, col = "red")
    abline(h = -1.96, col = "red")
    title("Residuals and 95% Confidence Interval (red lines)")
    text(200, -3, paste(perc.withinCI, "% within 95% CI", sep=""))
  table(residuals$scale > -1.96 & residuals$scale < 1.96)
  hist(residuals$V1,main = "Distribution Innovation Factors Position", xlab = "Deviation Predicted - Measured Position (meters)")
# velocity
  residuals = read.table("testresidualsvelocity.txt")
  hist(residuals$V1,main = "Distribution Innovation Factors Velocity", xlab = "Deviation Predicted - Measured Velocity (meters/s)")
############################################################

test = read.table("normalized_if_pos.txt")
  
