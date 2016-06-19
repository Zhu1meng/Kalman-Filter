rm(list=ls(all=TRUE))
##########################################
##                                      ##
##      KALMAN FILTER SIMULATION        ##
##      MP Roeling, miniproject 1       ##
##      Oxford University               ##   
##      June 2016                       ##
##                                      ##
##########################################

setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/")
library("numDeriv")
library("rgl")
source("linearizedfunctions.r")

 file.remove("Jk_longitude.txt")
 file.remove("Jk_latitude.txt")
 file.remove("Jk_altitude.txt")
 file.remove("testK.txt")
 file.remove("testPkp.txt")
 file.remove("testresidualslongitude.txt")
 file.remove("testresidualslatitude.txt")
 file.remove("testresidualsaltitude.txt")
 file.remove("normalized_if_lon.txt")
 file.remove("normalized_if_lat.txt")
 file.remove("normalized_if_alt.txt")
 
attack.model = 4

###################################################################################
output = read.table("data_with_noise.txt", sep="\t", h = T)

# Rplot_linear_distance_time  
# plot(output$radar.distance, xlab = "Time (sec)", ylab = "Distance (meters)", main = "Raw data distance to radar over time")

# calculate velocity
# normally you would calculate the velocity on the fly, but here  we already have the data
#output$velocity = 0
#for(j in 2:nrow(output)){
#  output[j,]$velocity = output[j,3] - output[j-1,3] 
#}
  
start.lon = output[1,]$noiselong
start.lat = output[1,]$noiselat
start.alt = output[1,]$altitude_doub
start.long.vel = output[2,]$noiselong - output[1,]$noiselong
start.lat.vel = output[2,]$noiselat - output[1,]$noiselat
start.alt.vel = output[2,]$altitude_doub - output[1,]$altitude_doub
  
# QC step remove first row of data
data = output[c(2:nrow(output)),c(11,12,4)]
nvar = ncol(data)
#### Kalman filter process ####
  
# H = transformation matrix, R = error in the observation, Y = observation matrix
# 1. Xkp = AXk-1 + BUk + Wk  -> predicted state
# 2. get the starting values for Pkp
# 3. Pkp = APk-1AT + Qk -> process cov matrix
# 4. K = Pkp HT / KPkpHT + R
# 5. Yk = CYkm + Zm
# 6. Xk = Xkp + K[Y-HXkp]
# 7. Pk = (I-KH)Pkp & Xk
# A allows to update position and velocity
# A[row(A) == (col(A) - 1)]=Delta.T
# B translates acceleration (Ax) into adjustment to position and velocity,
# B = matrix(c(1/2*Delta.T^2,Delta.T),nvar,1)
# acceleration is not constant (decrease and increase, get estimate from velocity
# which itself is an estimate from distance over time

number = nrow(data)
graphdata=matrix(0,(nrow(data)-1),nvar*2+1)
if(exists("i")){rm(i)}

# initial conditions
# Ax = 0 # acceleration speed up / down 10 m/s, defined in the loop
Delta.T = 1 # timestep 1 second (unclear if the timesteps are seconds but we assume every step is a second
              # otherwise the dynamical model for velocity would not be accurate)

Delta.X = 0.0015 # uncertainty in measurement lon, approximately 100 meters
Delta.Y = 0.0010 # uncertainty in measurement lat, approximately 100 meters
Delta.Z = 100 # uncertainty in measurement alt, approximately 100 meters

# process error in process covariance matrix, is updated at the end of iteration 1.
Delta.Px = 0.015 # 1 km start error lon
Delta.Py = 0.015 # 1 km start error lat
Delta.Pz = 0.015 # 1 km start error alt

Delta.Pvx = 0.0045 # 300 m/s start error lon
Delta.Pvy = 0.0045 # 300 m/s start error lat
Delta.Pvz = 300 # 300 m/s start error alt

# observation errors
Delta.X = 0.015 # 1000 meter sd
Delta.Y = 0.015 # 1000 meter sd
Delta.Z = 1000 # 1000 meter sd
Delta.Vx = 0.00375 # 250 m/s sd
Delta.Vy = 0.00375 # 250 m/s sd
Delta.Vz = 250 # 250 m/s sd
  
A = diag(nvar*2)
  A[row(A) == (col(A) - nvar)] = 1

C = diag(nvar*2)
I = diag(nvar*2)
H = diag(nvar*2) 

var.velocity.lon = var(data[,1])
var.velocity.lat = var(data[,2])
var.velocity.alt = var(data[,3])

Qr = 500 # error in the process of calculating the process cov matrix, between 0.5 and 1 Ax 
Zerror = 0 # observation errors in mechanism (eg. electronic delays)

# simulate white noise
white_noise_gaussian_lon = rnorm(nrow(data), 0, sd(data[,1]))
white_noise_gaussian_lat = rnorm(nrow(data), 0, sd(data[,2]))
white_noise_gaussian_alt = rnorm(nrow(data), 0, sd(data[,3]))

white_noise_gaussian_lon = rnorm(nrow(data), 0, 0.0015)
white_noise_gaussian_lat = rnorm(nrow(data), 0, 0.0015)
white_noise_gaussian_alt = rnorm(nrow(data), 0, 100)

# attack model
attack.model = 0
# 0 = no attack
# 1 = maximum magnitude 
# 2 = wave-based
# 3 = positive deviation
# 4 = negative deviation
# anomaly detection threshold
lambda.max = matrix(c(1.5),6,1) # anomaly detection threshold

set.seed(100)
for (i in 1:nrow(data)){
# for (i in 1:100){
# 1. Xkp = AXk-1 + BUk + Wk -> predicted state
  if(i==1){
    # first iteration take start values
    long = start.lon # initial position long
    lat = start.lat # initial position lat
    alt = start.alt # initial position altitude
    long.vel = start.long.vel # initial velocity lon
    lat.vel = start.lat.vel # initial velocity lat
    alt.vel = start.alt.vel # initial velocity alt
    
    # accelaration
    Acc.lon = start.long.vel
    Acc.lat = start.lat.vel
    Acc.alt = start.alt.vel
  } else {
    # second or more iterations take previous estimate
    long = Xk[1]
    lat = Xk[2]
    alt = Xk[3]
    long.vel = Xk[4]
    lat.vel = Xk[5]
    alt.vel = Xk[6]
  }
  
  # white gaussian noise
  Wk.long = matrix(c(sample(white_noise_gaussian_lon, 1),
                     sample(white_noise_gaussian_lon, 1)),2,1)
  Wk.lat = matrix(c(sample(white_noise_gaussian_lat, 1),
                    sample(white_noise_gaussian_lat, 1)),2,1)
  Wk.alt = matrix(c(sample(white_noise_gaussian_alt, 1),
                    sample(white_noise_gaussian_alt, 1)),2,1)
  Wk = matrix(c(Wk.long, Wk.lat, Wk.alt),6,1)
  #Wk = matrix(c(1,1,1,1,1,1),6,1)
  
  if(i < nrow(data)){
    # acceleration 
    Acc.lon = data[i+1,1]-data[i,1] 
    Acc.lat = data[i+1,2]-data[i,2] 
    Acc.alt = data[i+1,3]-data[i,3] 

    A.X = matrix(c(long + long.vel * Delta.T,
                   lat + lat.vel * Delta.T,
                   alt + alt.vel * Delta.T,
                   long.vel,
                   lat.vel,
                   alt.vel),6,1)
    
    B.Acc = matrix(c(Acc.lon * 0.5 * Delta.T^2,
                     Acc.lat * 0.5 * Delta.T^2,
                     Acc.alt * 0.5 * Delta.T^2,
                     Acc.lon * Delta.T,
                     Acc.lat * Delta.T,
                     Acc.alt * Delta.T),6,1) 
    
    Xkp = nonlinear.prediction(A.X,B.Acc,Wk,Delta.T,nvar)
  
  } else {
    # acceleration 
    Acc.lon = 0
    Acc.lat = 0
    Acc.alt = 0
    
    A.X = matrix(c(long + long.vel * Delta.T,
                   lat + lat.vel * Delta.T,
                   alt + alt.vel * Delta.T,
                   long.vel,
                   lat.vel,
                   alt.vel),6,1)
    
    B.Acc = matrix(c(Acc.lon * 0.5 * Delta.T^2,
                     Acc.lat * 0.5 * Delta.T^2,
                     Acc.alt * 0.5 * Delta.T^2,
                     Acc.lon * Delta.T,
                     Acc.lat * Delta.T,
                     Acc.alt * Delta.T),6,1) 
    
    Xkp = nonlinear.prediction(A.X,B.Acc,Wk,Delta.T,nvar)
  }
  
  # I let the loop continue to this point to get the extra prediction Xkp
  if(i==nrow(data)){break} 

# 2. Pkp.i = APk-1AT + Qk -> initial process cov matrix
  if(i==1){ # first iteration take start values
    Pkp.i=rbind(matrix(c(Delta.Px*Delta.Px,  Delta.Px*Delta.Py, Delta.Px*Delta.Pz, Delta.Px*Delta.Pvx, Delta.Px*Delta.Pvy, Delta.Px*Delta.Pvz),1,6),
                matrix(c(Delta.Py*Delta.Px,  Delta.Py*Delta.Py, Delta.Py*Delta.Pz, Delta.Py*Delta.Pvx, Delta.Py*Delta.Pvy, Delta.Py*Delta.Pvz),1,6),
                matrix(c(Delta.Pz*Delta.Px,  Delta.Pz*Delta.Py, Delta.Pz*Delta.Pz, Delta.Pz*Delta.Pvx, Delta.Pz*Delta.Pvy, Delta.Pz*Delta.Pvz),1,6),
                matrix(c(Delta.Pvx*Delta.Pvx,Delta.Pvx*Delta.Py,Delta.Pvx*Delta.Pz,Delta.Pvx*Delta.Pvx,Delta.Pvx*Delta.Pvy,Delta.Pvx*Delta.Pvz),1,6),
                matrix(c(Delta.Pvy*Delta.Pvy,Delta.Pvy*Delta.Py,Delta.Pvy*Delta.Pz,Delta.Pvy*Delta.Pvx,Delta.Pvy*Delta.Pvy,Delta.Pvy*Delta.Pvz),1,6),
                matrix(c(Delta.Pvz*Delta.Pvz,Delta.Pvz*Delta.Py,Delta.Pvz*Delta.Pz,Delta.Pvz*Delta.Pvx,Delta.Pvz*Delta.Pvy,Delta.Pvz*Delta.Pvz),1,6))

    # setting the cross terms to zero
    Pkp.i[row(Pkp.i) == (col(Pkp.i) - 1)]=0 # upper offdiagonal
    Pkp.i[row(Pkp.i) == (col(Pkp.i) - 2)]=0 # upper offdiagonal
    Pkp.i[row(Pkp.i) == (col(Pkp.i) - 3)]=0 # upper offdiagonal
    Pkp.i[row(Pkp.i) == (col(Pkp.i) - 4)]=0 # upper offdiagonal
    Pkp.i[row(Pkp.i) == (col(Pkp.i) - 5)]=0 # upper offdiagonal
    Pkp.i[row(Pkp.i) == (col(Pkp.i) + 1)]=0 # lower offdiagonal
    Pkp.i[row(Pkp.i) == (col(Pkp.i) + 2)]=0 # lower offdiagonal
    Pkp.i[row(Pkp.i) == (col(Pkp.i) + 3)]=0 # lower offdiagonal
    Pkp.i[row(Pkp.i) == (col(Pkp.i) + 4)]=0 # lower offdiagonal
    Pkp.i[row(Pkp.i) == (col(Pkp.i) + 5)]=0 # lower offdiagonal
    
  } else {
    Pkp.i = Pk # second or more iterations take previous estimate
    Qr.1 = matrix(c(Delta.T^4/4,Delta.T^3/2,
                  Delta.T^3/2,Delta.T^2),2,2) # formula see page 234 of Kalman Filter Python book
    Qr.long = Qr.1 * var.velocity.lon
    Qr.lat = Qr.1 * var.velocity.lat
    Qr.alt = Qr.1 * var.velocity.alt
    #Qr = 0
    Qr = rbind(matrix(c(Qr.long[1,],0,0,0,0),1,6),
               matrix(c(Qr.long[2,],0,0,0,0),1,6),
               matrix(c(0,0,Qr.lat[1,],0,0),1,6),
               matrix(c(0,0,Qr.lat[2,],0,0),1,6),
               matrix(c(0,0,0,0,Qr.alt[1,]),1,6),
               matrix(c(0,0,0,0,Qr.alt[2,]),1,6))
  }
  
# 3. Pkp = APkp.iAT + QR >- Predicted / adjusted process cov matrix
  Pkp = A%*%Pkp.i%*%t(A) + Qr
  # setting the cross terms to zero
  Pkp[row(Pkp) == (col(Pkp) - 1)]=0 # upper offdiagonal
  Pkp[row(Pkp) == (col(Pkp) - 2)]=0 # upper offdiagonal
  Pkp[row(Pkp) == (col(Pkp) - 3)]=0 # upper offdiagonal
  Pkp[row(Pkp) == (col(Pkp) - 4)]=0 # upper offdiagonal
  Pkp[row(Pkp) == (col(Pkp) - 5)]=0 # upper offdiagonal
  Pkp[row(Pkp) == (col(Pkp) + 1)]=0 # lower offdiagonal
  Pkp[row(Pkp) == (col(Pkp) + 2)]=0 # lower offdiagonal
  Pkp[row(Pkp) == (col(Pkp) + 3)]=0 # lower offdiagonal
  Pkp[row(Pkp) == (col(Pkp) + 4)]=0 # lower offdiagonal
  Pkp[row(Pkp) == (col(Pkp) + 5)]=0 # lower offdiagonal
  
# 4. K = Pkp HT / KPkpHT + R -> Calculate the Kalman Gain
  R = rbind(matrix(c(Delta.X^2,0,0,0,0,0),1,6),
            matrix(c(0,Delta.Y^2,0,0,0,0),1,6),
            matrix(c(0,0,Delta.Z^2,0,0,0),1,6),
            matrix(c(0,0,0,Delta.Vx^2,0,0),1,6),
            matrix(c(0,0,0,0,Delta.Vy^2,0),1,6),
            matrix(c(0,0,0,0,0,Delta.Vz^2),1,6))
  # Pkp%*%t(H) <- upper term
  # H%*%Pkp*t(H)+R <- lower term
  # dividing equals multiplication with inverse (solve function)
  K = Pkp%*%t(H)*solve(H%*%Pkp*t(H)+R)

# 5. Yk = CYkm + Zm -> represents observed state
  # Zerror = matrix(c(sample(white_noise_gaussian, 1),
  #                   sample(white_noise_gaussian, 1)),nvar,1)
  Ykm = matrix(c(data[i+1,1],data[i+1,2],data[i+1,3],Acc.lon,Acc.lat,Acc.alt),nvar*2,1)
  Yk = nonlinear.state(Ykm,C,Zerror,Delta.T,nvar)

# 6. Xk = Xkp + K[Y-HXkp] -> calculate current state
  Xk = Xkp + K%*%(Yk-H%*%Xkp)
 
  ############# DATA INJECTION ###############
  ## Yk = matrix with distance and velocity
  innovationfactor = Yk-H%*%Xkp   # formula 15 yang et al.
  rho2 = H%*%Pkp*t(H)+R           # formula 17 yang etal.
  # formula 16 yang et al.
  normalized_innovationfactor = innovationfactor / matrix(c(sqrt(rho2[1]), sqrt(rho2[8]), sqrt(rho2[15]), 
                                                            sqrt(rho2[22]), sqrt(rho2[29]), sqrt(rho2[36])),6,1)
  
  ### BLOCK FOR ATTACK MODEL 3+4 ###
  ## anomaly detection threshold = lambda.max 
  Zk.lowerbound = Xkp + (lambda.max * matrix(c(sqrt(rho2[1]), sqrt(rho2[8]), sqrt(rho2[15]), 
                                               sqrt(rho2[22]), sqrt(rho2[29]), sqrt(rho2[36])),6,1)) # formula 21 Yang et al.
  Zk.upperbound = Xkp - (lambda.max * matrix(c(sqrt(rho2[1]), sqrt(rho2[8]), sqrt(rho2[15]), 
                                               sqrt(rho2[22]), sqrt(rho2[29]), sqrt(rho2[36])),6,1)) # formula 21 Yang et al.
  
  ## substract the measurement from Z value to get attack vector c -> c  = z - y
  attack.vector.c.upperrange = Zk.upperbound - Xkp 
  attack.vector.c.lowerrange = Xkp - Zk.upperbound
  ## Multiply the attack vector c with the Kalman Gain, formula 37 Yang et al.: a = Kc 
  attack.vector.a.lower = K%*%attack.vector.c.lowerrange
  attack.vector.a.upper = K%*%attack.vector.c.upperrange
  ## prediction, formula 38 yang et al.
  Xkp.attack.low = Xkp + attack.vector.a.lower 
  Xkp.attack.up = Xkp + attack.vector.a.upper
  ## estimate, formula formula 40 yang et al
  Yk.attack.low = Yk + attack.vector.c.lowerrange 
  Yk.attack.up = Yk + attack.vector.c.upperrange  
  ### END BLOCK FOR ATTACK MODEL 3+4 ### 
  
  if (attack.model == 0){
    # 6. Xk = Xkp + K[Y-HXkp] -> calculate current state
    Xk = Xkp + K%*%(Yk-H%*%Xkp)
    # create performance index Jk
    # Jk =  abs(estimated measurement vector - true vector of measurements) / 
    #       abs(noise (real) measurement vector - true vector of measurements)
    Jk =  abs(Xkp - Xk)  / abs(Yk - Xk)
    lonperformance = cbind(Jk[1], Xkp[1], Xk[1], Yk[1])
    latperformance = cbind(Jk[2], Xkp[2], Xk[2], Yk[2])
    altperformance = cbind(Jk[3], Xkp[3], Xk[3], Yk[3])
    graphdata[i,1]=i
    graphdata[i,2]=Xkp[1] # predicted long
    graphdata[i,3]=Xkp[2] # predicted lat
    graphdata[i,4]=Xkp[3] # predicted alt
    graphdata[i,5]=Xk[1] # Kalman filter state position
    graphdata[i,6]=Xk[2] # Kalman filter state velocity
    graphdata[i,7]=Xk[3] # Kalman filter state velocity

  # 1 = maximum magnitude
  } else if (attack.model == 1){   
    # determine whether innovation vector is positive or negative 
    # formula 23 Yang et al.
    if(innovationfactor[1]>=0){
      ## anomaly detection threshold = lambda.max 
      Zk = Xkp + (lambda.max * matrix(c(sqrt(rho2[1]), sqrt(rho2[8]), sqrt(rho2[15]), 
                                                   sqrt(rho2[22]), sqrt(rho2[29]), sqrt(rho2[36])),6,1))  # formula 22, 24 Yang et al.
   
      ## substract the measurement from Z value to get attack vector c -> c  = z - y
      attack.vector.c = Zk - Yk
      ## Multiply the attack vector c with the Kalman Gain, formula 37 Yang et al.: a = Kc 
      attack.vector.a = K%*%attack.vector.c
      ## prediction, formula 38 yang et al.
      Xkp.attack = Xkp - attack.vector.a
      ## estimate, formula formula 40 yang et al
      Yk.attack = Yk - attack.vector.c  
      Xk = Xkp + K%*%(Yk-H%*%Xkp) # the normal run, needed for Xkp
      Xkp.a = Xkp.attack
      Yk.a = Yk.attack
      Xk.a = Xkp.a + K%*%(Yk.a-H%*%Xkp.a)
      
      # create performance index Jk
      # Jk =  abs(estimated measurement vector - true vector of measurements) / 
      #       abs(noise (real) measurement vector - true vector of measurements)
      Jk =  abs(Xkp.a - Xk.a)  / abs(Yk.a - Xk.a)
      pperformance = cbind(Jk[1], Xkp.a[1], Xk.a[1], Yk.a[1])
      vperformance = cbind(Jk[2], Xkp.a[2], Xk.a[2], Yk.a[2])
      graphdata[i,1]=i
      graphdata[i,2]=Xkp.a[1] # predicted long
      graphdata[i,3]=Xkp.a[2] # predicted lat
      graphdata[i,4]=Xkp.a[3] # predicted alt
      graphdata[i,5]=Xk.a[1] # Kalman filter state position
      graphdata[i,6]=Xk.a[2] # Kalman filter state velocity
      graphdata[i,7]=Xk.a[3] # Kalman filter state velocity
    } else if(innovationfactor[1]<0){
      ## anomaly detection threshold = lambda.max 
      Zk = Xkp - (lambda.max * matrix(c(sqrt(rho2[1]), sqrt(rho2[8]), sqrt(rho2[15]), 
                                        sqrt(rho2[22]), sqrt(rho2[29]), sqrt(rho2[36])),6,1))  # formula 22, 24 Yang et al.

      ## substract the measurement from Z value to get attack vector c -> c  = z - y
      attack.vector.c = Zk - Yk
      ## Multiply the attack vector c with the Kalman Gain, formula 37 Yang et al.: a = Kc 
      attack.vector.a = K%*%attack.vector.c
      ## prediction, formula 38 yang et al.
      Xkp.attack = Xkp - attack.vector.a
      ## estimate, formula formula 40 yang et al
      Yk.attack = Yk - attack.vector.c  
      Xk = Xkp + K%*%(Yk-H%*%Xkp) # the normal run, needed for Xkp
      Xkp.a = Xkp.attack
      Yk.a = Yk.attack
      Xk.a = Xkp.a + K%*%(Yk.a-H%*%Xkp.a)
      
      # create performance index Jk
      # Jk =  abs(estimated measurement vector - true vector of measurements) / 
      #       abs(noise (real) measurement vector - true vector of measurements)
      Jk =  abs(Xkp.a - Xk.a)  / abs(Yk.a - Xk.a)
      lonperformance = cbind(Jk[1], Xkp.a[1], Xk.a[1], Yk.a[1])
      latperformance = cbind(Jk[2], Xkp.a[2], Xk.a[2], Yk.a[2])
      altperformance = cbind(Jk[3], Xkp.a[3], Xk.a[3], Yk.a[3])
      graphdata[i,1]=i
      graphdata[i,2]=Xkp.a[1] # predicted long
      graphdata[i,3]=Xkp.a[2] # predicted lat
      graphdata[i,4]=Xkp.a[3] # predicted alt
      graphdata[i,5]=Xk.a[1] # Kalman filter state long
      graphdata[i,6]=Xk.a[2] # Kalman filter state lat
      graphdata[i,7]=Xk.a[3] # Kalman filter state alt
    }
    
  # 2 = wave-based (identical to attack model 1 but reversed direction innovation factor)
  } else if (attack.model == 2){
    if(innovationfactor[1]<0){
      ## anomaly detection threshold = lambda.max 
      Zk = abs(Xkp) + (lambda.max * matrix(c(sqrt(rho2[1]), sqrt(rho2[8]), sqrt(rho2[15]), 
                                        sqrt(rho2[22]), sqrt(rho2[29]), sqrt(rho2[36])),6,1))  # formula 22, 24 Yang et al.

      ## substract the measurement from Z value to get attack vector c -> c  = z - y
      attack.vector.c = Zk - Yk
      ## Multiply the attack vector c with the Kalman Gain, formula 37 Yang et al.: a = Kc 
      attack.vector.a = K%*%attack.vector.c
      ## prediction, formula 38 yang et al.
      Xkp.attack = Xkp - attack.vector.a
      ## estimate, formula formula 40 yang et al
      Yk.attack = Yk - attack.vector.c  
      Xk = Xkp + K%*%(Yk-H%*%Xkp) # the normal run, needed for Xkp
      Xkp.a = Xkp.attack
      Yk.a = Yk.attack
      Xk.a = Xkp.a + K%*%(Yk.a-H%*%Xkp.a)
      
      # create performance index Jk
      # Jk =  abs(estimated measurement vector - true vector of measurements) / 
      #       abs(noise (real) measurement vector - true vector of measurements)
      Jk =  abs(Xkp.a - Xk.a)  / abs(Yk.a - Xk.a)
      lonperformance = cbind(Jk[1], Xkp.a[1], Xk.a[1], Yk.a[1])
      latperformance = cbind(Jk[2], Xkp.a[2], Xk.a[2], Yk.a[2])
      altperformance = cbind(Jk[3], Xkp.a[3], Xk.a[3], Yk.a[3])
      graphdata[i,1]=i
      graphdata[i,2]=Xkp.a[1] # predicted long
      graphdata[i,3]=Xkp.a[2] # predicted lat
      graphdata[i,4]=Xkp.a[3] # predicted alt
      graphdata[i,5]=Xk.a[1] # Kalman filter state long
      graphdata[i,6]=Xk.a[2] # Kalman filter state lat
      graphdata[i,7]=Xk.a[3] # Kalman filter state alt
    } else if(innovationfactor[1]>=0){
      ## anomaly detection threshold = lambda.max 
      Zk = abs(Xkp) - (lambda.max * matrix(c(sqrt(rho2[1]), sqrt(rho2[8]), sqrt(rho2[15]), 
                                             sqrt(rho2[22]), sqrt(rho2[29]), sqrt(rho2[36])),6,1))  # formula 22, 24 Yang et al.
      ## substract the measurement from Z value to get attack vector c -> c  = z - y
      attack.vector.c = Zk - Yk
      ## Multiply the attack vector c with the Kalman Gain, formula 37 Yang et al.: a = Kc 
      attack.vector.a = K%*%attack.vector.c
      ## prediction, formula 38 yang et al.
      Xkp.attack = Xkp - attack.vector.a
      ## estimate, formula formula 40 yang et al
      Yk.attack = Yk - attack.vector.c  
      Xk = Xkp + K%*%(Yk-H%*%Xkp) # the normal run, needed for Xkp
      Xkp.a = Xkp.attack
      Yk.a = Yk.attack
      Xk.a = Xkp.a + K%*%(Yk.a-H%*%Xkp.a)
      
      # create performance index Jk
      # Jk =  abs(estimated measurement vector - true vector of measurements) / 
      #       abs(noise (real) measurement vector - true vector of measurements)
      Jk =  abs(Xkp.a - Xk.a)  / abs(Yk.a - Xk.a)
      lonperformance = cbind(Jk[1], Xkp.a[1], Xk.a[1], Yk.a[1])
      latperformance = cbind(Jk[2], Xkp.a[2], Xk.a[2], Yk.a[2])
      altperformance = cbind(Jk[3], Xkp.a[3], Xk.a[3], Yk.a[3])
      graphdata[i,1]=i
      graphdata[i,2]=Xkp.a[1] # predicted long
      graphdata[i,3]=Xkp.a[2] # predicted lat
      graphdata[i,4]=Xkp.a[3] # predicted alt
      graphdata[i,5]=Xk.a[1] # Kalman filter state long
      graphdata[i,6]=Xk.a[2] # Kalman filter state lat
      graphdata[i,7]=Xk.a[3] # Kalman filter state alt
    }
  
  # 3 = positive deviation
  } else if (attack.model == 3){
    ## anomaly detection threshold = lambda.max 
    Zk = Xkp + (lambda.max * matrix(c(sqrt(rho2[1]), sqrt(rho2[8]), sqrt(rho2[15]), 
                                           sqrt(rho2[22]), sqrt(rho2[29]), sqrt(rho2[36])),6,1))
    ## substract the measurement from Z value to get attack vector c -> c  = z - y
    attack.vector.c = Zk - Yk
    ## Multiply the attack vector c with the Kalman Gain, formula 37 Yang et al.: a = Kc 
    attack.vector.a = K%*%attack.vector.c
    ## prediction, formula 38 yang et al.
    Xkp.attack = Xkp - attack.vector.a
    ## estimate, formula formula 40 yang et al
    Yk.attack = Yk - attack.vector.c  
    Xk = Xkp + K%*%(Yk-H%*%Xkp) # the normal run, needed for Xkp
    Xkp.a = Xkp.attack
    Yk.a = Yk.attack
    Xk.a = Xkp.a + K%*%(Yk.a-H%*%Xkp.a)
    
    # create performance index Jk
    # Jk =  abs(estimated measurement vector - true vector of measurements) / 
    #       abs(noise (real) measurement vector - true vector of measurements)
    Jk =  abs(Xkp.a - Xk.a)  / abs(Yk.a - Xk.a)
    lonperformance = cbind(Jk[1], Xkp.a[1], Xk.a[1], Yk.a[1])
    latperformance = cbind(Jk[2], Xkp.a[2], Xk.a[2], Yk.a[2])
    altperformance = cbind(Jk[3], Xkp.a[3], Xk.a[3], Yk.a[3])
    graphdata[i,1]=i
    graphdata[i,2]=Xkp.a[1] # predicted long
    graphdata[i,3]=Xkp.a[2] # predicted lat
    graphdata[i,4]=Xkp.a[3] # predicted alt
    graphdata[i,5]=Xk.a[1] # Kalman filter state long
    graphdata[i,6]=Xk.a[2] # Kalman filter state lat
    graphdata[i,7]=Xk.a[3] # Kalman filter state alt
    
  # 4 = negative deviation
  } else if (attack.model == 4){
    ## anomaly detection threshold = lambda.max 
    Zk = Xkp - (lambda.max * matrix(c(sqrt(rho2[1]), sqrt(rho2[8]), sqrt(rho2[15]), 
                                      sqrt(rho2[22]), sqrt(rho2[29]), sqrt(rho2[36])),6,1))
    ## substract the measurement from Z value to get attack vector c -> c  = z - y
    attack.vector.c = Zk - Yk
    ## Multiply the attack vector c with the Kalman Gain, formula 37 Yang et al.: a = Kc 
    attack.vector.a = K%*%attack.vector.c
    ## prediction, formula 38 yang et al.
    Xkp.attack = Xkp - attack.vector.a
    ## estimate, formula formula 40 yang et al
    Yk.attack = Yk - attack.vector.c  
    Xk = Xkp + K%*%(Yk-H%*%Xkp) # the normal run, needed for Xkp
    Xkp.a = Xkp.attack
    Yk.a = Yk.attack
    Xk.a = Xkp.a + K%*%(Yk.a-H%*%Xkp.a)
    
    # create performance index Jk
    # Jk =  abs(estimated measurement vector - true vector of measurements) / 
    #       abs(noise (real) measurement vector - true vector of measurements)
    Jk =  abs(Xkp.a - Xk.a)  / abs(Yk.a - Xk.a)
    lonperformance = cbind(Jk[1], Xkp.a[1], Xk.a[1], Yk.a[1])
    latperformance = cbind(Jk[2], Xkp.a[2], Xk.a[2], Yk.a[2])
    altperformance = cbind(Jk[3], Xkp.a[3], Xk.a[3], Yk.a[3])
    graphdata[i,1]=i
    graphdata[i,2]=Xkp.a[1] # predicted long
    graphdata[i,3]=Xkp.a[2] # predicted lat
    graphdata[i,4]=Xkp.a[3] # predicted alt
    graphdata[i,5]=Xk.a[1] # Kalman filter state long
    graphdata[i,6]=Xk.a[2] # Kalman filter state lat
    graphdata[i,7]=Xk.a[3] # Kalman filter state alt
  } else {
    # do nothing
    print("no attack model specified")
    stop
  }

  ################ END INJECTION #################

  # 7. Pk = (I-KH)Pkp & Xk -> update process cov matrix
  Pk = (I - K*H) %*% Pkp

# write data for subsequent analyses
  write.table(lonperformance, "Jk_longitude.txt", append = T, col.names = F, row.names = F)
  write.table(latperformance, "Jk_latitude.txt", append = T, col.names = F, row.names = F)
  write.table(altperformance, "Jk_altitude.txt", append = T, col.names = F, row.names = F)
  write.table(K[1], "testK.txt", append = T, col.names = F, row.names = F)
  write.table(Pkp[1], "testPkp.txt", append = T, col.names = F, row.names = F)
  write.table(innovationfactor[1], "testresidualslongitude.txt", append = T, col.names = F, row.names = F)
  write.table(innovationfactor[2], "testresidualslatitude.txt", append = T, col.names = F, row.names = F)
  write.table(innovationfactor[3], "testresidualsaltitude.txt", append = T, col.names = F, row.names = F)
  write.table(normalized_innovationfactor[1], "normalized_if_lon.txt", append = T, col.names = F, row.names = F)
  write.table(normalized_innovationfactor[2], "normalized_if_lat.txt", append = T, col.names = F, row.names = F)
  write.table(normalized_innovationfactor[3], "normalized_if_alt.txt", append = T, col.names = F, row.names = F)
  ow <- options(warn = 2) 
} 

########### END LOOP ###########

## GRAPH ## 
# there will always be an extra prediction, which is unreliable because it 
# does not take the Nth measurement into account at Ykm, but takes the previous N-1 Ykm
data.select = cbind(seq(1,nrow(data),1),data,1)
  colnames(data.select)<-c("Time","Longitude","Lattitude","Altitude","Group")
predicted = cbind(seq(1,(nrow(graphdata)+2),1),
                  rbind(c(start.lon,start.lat,start.alt),graphdata[1:nrow(graphdata),c(2,3,4)],c(Xkp[1],Xkp[2],Xkp[3])),2) # predicted
  colnames(predicted)<-c("Time","Longitude","Lattitude","Altitude","Group")
kalman = cbind(seq(1,(nrow(graphdata)+1),1),
                  rbind(c(start.lon,start.lat,start.alt),graphdata[1:nrow(graphdata),c(5,6,7)]),3) # kalman
  colnames(kalman)<-c("Time","Longitude","Lattitude","Altitude","Group")

graphdata = rbind(data.select, predicted, kalman)
# zoom
"graphdata = rbind(data.select[c(1:100),], predicted[c(1:100),], kalman[c(1:100),])

# long lat first 100
plot3d(data.select[c(1:100),1], data.select[c(1:100),2], data.select[c(1:100),3])
plot3d(predicted[c(1:100),1], predicted[c(1:100),2], predicted[c(1:100),3])
plot3d(kalman[c(1:100),1], kalman[c(1:100),2], kalman[c(1:100),3])

plot3d(predicted[c(1:100),2], predicted[c(1:100),3], predicted[c(1:100),4])

plot3d(data.select[,2], data.select[,3], data.select[,4])
plot3d(predicted[,2], predicted[,3], predicted[,4])

"


plot3d(data.select[,2], data.select[,3], data.select[,4], col = "green", 
       xlab = "Longitude", ylab = "Latitude", zlab = "Altitude")

points3d(predicted[,2], predicted[,3], predicted[,4], col="red")
points3d(kalman[,2], kalman[,3], kalman[,4], col="blue")

#rgl.bbox(color=c("#333377","black"), emission="#3399FF",
#         specular="#3333FF", shininess=10, alpha=0.8, nticks = 3 ) 


# convert factor to numeric for convenience 
ngroups=nvar+1
# get the range for the x and y axis 
xrange <- c(0,max(graphdata[,1]))
yrange <- c(4,10)

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
legend(1, 85000, 1:ngroups, cex = 0.8, col=colors, lty = linetype, title = "Groups")
