rm(list=ls(all=TRUE))
##########################################
##                                      ##
##      KALMAN FILTER SIMULATION        ##
##      MP Roeling, miniproject 1       ##
##      Oxford University               ##   
##      May/June 2016                   ##
##                                      ##
##########################################

setwd("C:/Users/mproeling/Documents/Oxford/miniproject/Kalman/simulation/")

file.remove("Jk_position.txt")
file.remove("Jk_velocity.txt")
file.remove("testK.txt")
file.remove("testPkp.txt")
file.remove("testresidualsposition.txt")
file.remove("testresidualsvelocity.txt")
file.remove("testAx.txt")
file.remove("normalized_if_pos.txt")
file.remove("normalized_if_vel.txt")

attack.model = 1

##################################################################################
output = read.table("output.txt", sep="\t")
  colnames(output) = c("iteration", "timestep", "radar.distance")

# simulate noise
noise = rnorm(10000,0,250)
for(i in 1:nrow(output)){
  output[i,]$radar.distance = output[i,]$radar.distance + sample(noise,1)
}

# Rplot_linear_distance_time  
# plot(output$radar.distance, xlab = "Time (sec)", ylab = "Distance (meters)", main = "Raw data distance to radar over time")

# calculate velocity
# normally you would calculate the velocity on the fly, but here  we already have the data
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
  Zerror = 0 # observation errors in mechanism (eg. electronic delays)
  # simulate white noise
  white_noise_gaussian = rnorm(nrow(data), 0, sd(data[,2]))
  # attack model
  # 0 = no attack
  # 1 = maximum magnitude 
  # 2 = wave-based
  # 3 = positive deviation
  # 4 = negative deviation
  # anomaly detection threshold
  lambda.max = matrix(c(1.5,1.5),2,1) # anomaly detection threshold for distance and velocity

set.seed(100)
# for (i in 1:nrow(data)){
for (i in 1:100){
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

# 5. Yk = CYkm + Zm -> represents observed state
  # Zerror = matrix(c(sample(white_noise_gaussian, 1),
  #                   sample(white_noise_gaussian, 1)),nvar,1)
  Ykm = matrix(c(data[i+1,1],data[i+1,2]),nvar,1)
  Yk = C%*%Ykm+Zerror
  
  ############# DATA INJECTION ###############
  ## Yk = matrix with distance and velocity
  innovationfactor = Yk-H%*%Xkp   # formula 15 yang et al.
  rho2 = H%*%Pkp*t(H)+R           # formula 17 yang etal.
  normalized_innovationfactor = innovationfactor / matrix(c(sqrt(rho2)[1],sqrt(rho2)[4]),2,1) # formula 16 yang et al.
  
  ### BLOCK FOR ATTACK MODEL 3+4 ###
  ## anomaly detection threshold = lambda.max 
  Zk.lowerbound = Xkp + (lambda.max * matrix(c(sqrt(rho2)[1],sqrt(rho2)[4]),2,1)) # formula 21 Yang et al.
  Zk.upperbound = Xkp - (lambda.max * matrix(c(sqrt(rho2)[1],sqrt(rho2)[4]),2,1)) # formula 21 Yang et al.

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
    pperformance = cbind(Jk[1], Xkp[1], Xk[1], Yk[1])
    vperformance = cbind(Jk[2], Xkp[2], Xk[2], Yk[2])
    graphdata[i,1]=i
    graphdata[i,2]=Xkp[1] # predicted position
    graphdata[i,3]=Xkp[2] # predicted velocity
    graphdata[i,4]=Xk[1] # Kalman filter state position
    graphdata[i,5]=Xk[2] # Kalman filter state velocity
  # 1 = maximum magnitude
  } else if (attack.model == 1){   
    # determine whether innovation vector is positive or negative 
    # formula 23 Yang et al.
    if(innovationfactor[1]>=0){
      ## anomaly detection threshold = lambda.max 
      Zk = Xkp + (lambda.max * matrix(c(sqrt(rho2)[1],sqrt(rho2)[4]),2,1)) # formula 22, 24 Yang et al.
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
      graphdata[i,2]=Xkp.a[1] # predicted position
      graphdata[i,3]=Xkp.a[2] # predicted velocity
      graphdata[i,4]=Xk.a[1] # Kalman filter state position
      graphdata[i,5]=Xk.a[2] # Kalman filter state velocity
    } else if(innovationfactor[1]<0){
      ## anomaly detection threshold = lambda.max 
      Zk = Xkp - (lambda.max * matrix(c(sqrt(rho2)[1],sqrt(rho2)[4]),2,1)) # formula 22, 24 Yang et al.
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
      graphdata[i,2]=Xkp.a[1] # predicted position
      graphdata[i,3]=Xkp.a[2] # predicted velocity
      graphdata[i,4]=Xk.a[1] # Kalman filter state position
      graphdata[i,5]=Xk.a[2] # Kalman filter state velocity
    }
    
  # 2 = wave-based (identical to attack model 1 but reversed direction innovation factor)
  } else if (attack.model == 2){
    if(innovationfactor[1]<0){
      ## anomaly detection threshold = lambda.max 
      Zk = abs(Xkp) + (lambda.max * matrix(c(sqrt(rho2)[1],sqrt(rho2)[4]),2,1)) # formula 22, 24 Yang et al.
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
      graphdata[i,2]=Xkp.a[1] # predicted position
      graphdata[i,3]=Xkp.a[2] # predicted velocity
      graphdata[i,4]=Xk.a[1] # Kalman filter state position
      graphdata[i,5]=Xk.a[2] # Kalman filter state velocity
    } else if(innovationfactor[1]>=0){
      ## anomaly detection threshold = lambda.max 
      Zk = abs(Xkp) - (lambda.max * matrix(c(sqrt(rho2)[1],sqrt(rho2)[4]),2,1)) # formula 22, 24 Yang et al.
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
      graphdata[i,2]=Xkp.a[1] # predicted position
      graphdata[i,3]=Xkp.a[2] # predicted velocity
      graphdata[i,4]=Xk.a[1] # Kalman filter state position
      graphdata[i,5]=Xk.a[2] # Kalman filter state velocity
    }
    
  # 3 = positive deviation
  } else if (attack.model == 3){
    ## anomaly detection threshold = lambda.max 
    Zk = Xkp + (lambda.max * matrix(c(sqrt(rho2)[1],sqrt(rho2)[4]),2,1)) # formula 22, 24 Yang et al.
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
    graphdata[i,2]=Xkp.a[1] # predicted position
    graphdata[i,3]=Xkp.a[2] # predicted velocity
    graphdata[i,4]=Xk.a[1] # Kalman filter state position
    graphdata[i,5]=Xk.a[2] # Kalman filter state velocity
      
  # 4 = negative deviation
  } else if (attack.model == 4){
    ## anomaly detection threshold = lambda.max 
    Zk = Xkp - (lambda.max * matrix(c(sqrt(rho2)[1],sqrt(rho2)[4]),2,1)) # formula 22, 24 Yang et al.
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
    graphdata[i,2]=Xkp.a[1] # predicted position
    graphdata[i,3]=Xkp.a[2] # predicted velocity
    graphdata[i,4]=Xk.a[1] # Kalman filter state position
    graphdata[i,5]=Xk.a[2] # Kalman filter state velocity
    
  } else {
    # do nothing
    print("no attack model specified")
    stop
  }
################## END INJECTION #################

  # 7. Pk = (I-KH)Pkp & Xk -> update process cov matrix
  Pk = (I - K*H) %*% Pkp

# write data for subsequent analyses
  write.table(pperformance, "Jk_position.txt", append = T, col.names = F, row.names = F)
  write.table(vperformance, "Jk_velocity.txt", append = T, col.names = F, row.names = F)
  write.table(K[1], "testK.txt", append = T, col.names = F, row.names = F)
  write.table(Pkp[1], "testPkp.txt", append = T, col.names = F, row.names = F)
  write.table(innovationfactor[1], "testresidualsposition.txt", append = T, col.names = F, row.names = F)
  write.table(innovationfactor[2], "testresidualsvelocity.txt", append = T, col.names = F, row.names = F)
  write.table(Ax, "testAx.txt", append = T, col.names = F, row.names = F)
  write.table(normalized_innovationfactor[1], "normalized_if_pos.txt", append = T, col.names = F, row.names = F)
  write.table(normalized_innovationfactor[2], "normalized_if_vel.txt", append = T, col.names = F, row.names = F)
  ow <- options(warn = 2) 
} 

########### END LOOP ###########

## GRAPH ## 
# there will always be an extra prediction, which is unreliable because it 
# does not take the Nth measurement into account at Ykm, but takes the previous N-1 Ykm
data.select = cbind(seq(1,nrow(data),1),data,1)
  colnames(data.select)<-c("Time","Position","Velocity","Group")
predicted = cbind(seq(1,(nrow(graphdata)+2),1),rbind(c(X.0,Vx.0),graphdata[1:nrow(graphdata),c(2,3)],c(Xkp[1],Xkp[2])),2) # predicted
  colnames(predicted)<-c("Time","Position","Velocity","Group")
kalman = cbind(seq(1,(nrow(graphdata)+1),1),rbind(c(X.0,Vx.0),graphdata[1:nrow(graphdata),c(4,5)]),3) # kalman
  colnames(kalman)<-c("Time","Position","Velocity","Group")

# graphdata = rbind(data.select, predicted, kalman)
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
