rm(list=ls(all=TRUE))
##########################################
##                                      ##
##      KALMAN FILTER SIMULATION        ##
##      MP Roeling, miniproject 1       ##
##      Oxford University               ##   
##      April 2016                      ##
##                                      ##
##########################################

# covariance example [https://www.youtube.com/watch?v=9B5vEVjH2Pk]
nvar=3
ncases=5
A=matrix(c(c(90,90,60,30,30),              # math
           c(80,60,50,40,20),              # physics
           c(40,80,70,70,90)),ncases,nvar) # english
# unity matrix 
I=matrix(1,ncases,ncases)
# unity matrix * case scores = sum, *1/ncases = mean
A_distance=A-I%*%A*1/ncases
# covariance matrix ATA (AT = A transpose)
A_cov=t(A_distance)%*%A

# Kalman filter process
# H = transformation matrix, R = error in the observation, Y = observation matrix
# 1. Xkp = AXk-1 + BUk + Wk -> predicted state
# 2. Pkp = APk-1AT + Qk
# 3. process cov matrix
# 4. K = Pkp HT / KPkpHT + R
# 5. Yk = CYkm + Zm
# 6. Xk = Xkp + K[Y-HXkp]
# 7. Pk = (I-KH)Pkp & Xk

nvar=2
data=matrix(c(c(4000,4260,4550,4860,5110),c(280,282,285,286,290)),5,2)
graphdata=matrix(0,(nrow(data)-1),nvar*2+1)

# eigenlijk begin je pas data op te halen bij de 2e regel bij observed state schatting (stap 5) vs predicted
# die eerste regel is al gegeven.
for (i in 1:nrow(data)){
  # given
  Vx.0=280 ; X.0=4000
  Vy.0=120 ; Y.0=3000
  # initial conditions
  Ax=2
  Delta.T=1
  Vx=280
  Delta.X=25
  # process error in process covariance matrix
  Delta.Px=20
  Delta.Pvx=5
  # observation errors
  Delta.X=25
  Delta.Vx=6

# 1. Xkp = AXk-1 + BUk + Wk -> predicted state
  # A allows to update position and velocity
  A = diag(nvar) 
  A[row(A) == (col(A) - 1)]=Delta.T
  if(i==1){
    X = matrix(c(X.0,Vx.0),2,1) # first iteration take start values
  } else {
    X = Xk # second or more iterations take previous estimate
  }
  # B translates acceleration (Ax) into adjustment to position and velocity,
  B = matrix(c(1/2*Delta.T^2,Delta.T),nvar,1)
  # error (Wk) eq to 0
  Xkp = A%*%X + B%*%Ax + 0
  if(i==nrow(data)){break} # I let the loop continue to this point to get the extra prediction Xkp

# 2. Pkp.i = APk-1AT + Qk -> initial process cov matrix
  if(i==1){ # first iteration take start values
    Pkp.i = matrix(c(Delta.Px^2, Delta.Px*Delta.Pvx,
                   Delta.Px*Delta.Pvx, Delta.Pvx^2),2,2)
    # setting the cross terms to zero
    Pkp.i[row(Pkp.i) == (col(Pkp.i) - 1)]=0 # upper offdiagonal
    Pkp.i[row(Pkp.i) == (col(Pkp.i) + 1)]=0 # lower offdiagonal
  } else {
    Pkp.i = Pk # second or more iterations take previous estimate
  }
  
# 3. Pkp = APkp.iAT + QR >- Predicted / adjusted process cov matrix
  Qr = 0 # error in the process of calculating the process cov matrix
  Pkp = A%*%Pkp.i%*%t(A) + Qr
  # setting the cross terms to zero
  Pkp[row(Pkp) == (col(Pkp) - 1)]=0 # upper offdiagonal
  Pkp[row(Pkp) == (col(Pkp) + 1)]=0 # lower offdiagonal

# 4. K = Pkp HT / KPkpHT + R -> Calculate the Kalman Gain
  H = diag(nvar) # matrix to allows to change the format of Pkp
  R = matrix(c(Delta.X^2,0,0,Delta.Vx^2),2,2) # observation errors
  # Pkp%*%t(H) <- upper term
  # H%*%Pkp*t(H)+R <- lower term
  # dividing equals multiplication with inverse (solve function)
  K = Pkp%*%t(H)*solve(H%*%Pkp*t(H)+R)

# 5. Yk = CYkm + Zm -> represents observed state
  C = diag(nvar)
  Zk = 0 # observation errors in mechanism (eg. electronic delays)
  Ykm = matrix(c(data[i+1,1],data[i+1,2]),nvar,1)
  Yk = C%*%Ykm
  
# 6. Xk = Xkp + K[Y-HXkp] -> calculate current state
  Xk = Xkp + K%*%(Yk-H%*%Xkp)

# 7. Pk = (I-KH)Pkp & Xk -> update process cov matrix
  I = diag(nvar)
  Pk = (I - K*H) %*% Pkp

  graphdata[i,1]=i
  graphdata[i,2]=Xkp[1] # predicted position
  graphdata[i,3]=Xkp[2] # predicted velocity
  graphdata[i,4]=Xk[1] # Kalman filter state position
  graphdata[i,5]=Xk[2] # Kalman filter state velocity
}

## GRAPH ## 
# there will always be an extra prediction, which is unreliable because it 
# does not take the Nth measurement into account at Ykm, but takes the previous N-1 Ykm
graphdata=rbind(cbind(seq(1,nrow(data),1),data,1),
                cbind(seq(1,(nrow(graphdata)+2),1),rbind(c(X.0,Vx.0),graphdata[1:nrow(graphdata),c(2,3)],c(Xkp[1],Xkp[2])),2), # predicted
              cbind(seq(1,(nrow(graphdata)+1),1),rbind(c(X.0,Vx.0),graphdata[1:nrow(graphdata),c(4,5)]),3)) # kalman
colnames(graphdata)<-c("Time","Position","Velocity","Group")

# convert factor to numeric for convenience 
ngroups=nvar+1
# get the range for the x and y axis 
xrange <- c(0,max(graphdata[,1]))
yrange <- c(3800,max(graphdata[,2]))
# set up the plot
plot(xrange, yrange, type="n", xlab="Time (sec)",
     ylab="Position" ) 
colors <- rainbow(ngroups) 
linetype <- c(2,3,1) 
group1=graphdata[graphdata[,ncol(graphdata)]==1,]
group2=graphdata[graphdata[,ncol(graphdata)]==2,]
group3=graphdata[graphdata[,ncol(graphdata)]==3,]
lines(group1[,1],group1[,2], type="b", lwd=1.5, col=colors[1], lty=linetype[1])
lines(group2[,1],group2[,2], type="b", lwd=1.5, col=colors[2], lty=linetype[2])
lines(group3[,1],group3[,2], type="b", lwd=1.5, col=colors[3], lty=linetype[3])
# add a title and subtitle 
title("Observed (red), predicted (green), and Kalman (blue)")
# add a legend 
legend(xrange[1], yrange[2], 1:ngroups, cex=0.8, col=colors,
       lty=linetype, title="Groups")

############################################################






  
