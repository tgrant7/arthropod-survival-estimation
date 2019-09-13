#Simulate and analyze data for Scenario 4 and 5 in Grant et al.  

### Scenario 4 ################

#daily survival probabilities
S = c(0.53,0.68,0.90,0.94,0.97,0.96)

n=60
pupamax = 16

#number of eggs laid per day = B=L in this code

#normal distr with mean = 30 and SD = 15
#max 50 eggs laid per day
CF = 50/dnorm(30, mean = 30, sd = 15) #correction factor to make density at 30 = 50
L = round(CF*dnorm(1:60, mean = 30, sd = 15))
#L4.1=L
#L4.4=L
L4.7=L
L4.10=L

#max 25 eggs laid per day
CF = 25/dnorm(30, mean = 30, sd = 15) #correction factor to make density at 30 = 25
L = round(CF*dnorm(1:60, mean = 30, sd = 15))
#L4.2=L
#L4.5=L
L4.8=L
L4.11=L

#max 100 eggs laid per day
CF = 100/dnorm(30, mean = 30, sd = 15) #correction factor to make density at 30 = 100
L = round(CF*dnorm(1:60, mean = 30, sd = 15))
#L4.3=L
#L4.6=L
L4.9=L
L4.12=L

S = c(0.53,0.68,0.90,0.94,0.97,0.96)

MT = rep(21,n); length(MT) #for fig in ms

## ADD ####
ADD = matrix(nrow=n+pupamax, ncol=n); colnames(ADD) = 1:n
for (j in 1:n){
  for (i in (j+1):(n+1)) {
    ADD[i,j] = sum(MT[j:(i-1)])
  }
}
#set upper right off diagonal cells to 0 - set to 0 so following sum operations will work smoothly
for (j in 1:n) {
  for (i in 1:j){
    ADD[i,j] = 0    
  }
}
#set lower cells to 0 as well
for (j in 1:n) {
  for (i in (n+2):(n+pupamax)){
    ADD[i,j] = 0    
  }
}

#calculate stage durations for each cohort. 
#note that some numbers in the lower right of the matrix are not correct, but are not used further in the code - these cohorts didn't complete development during the time period of the study
#rows are stages, cols are day/cohort
#this matrix corresponds to S_jk and c_ik in Grant et al.
SDC = matrix(nrow = 7, ncol = n); colnames(SDC) = 1:n; rownames(SDC) = c("egg","first","second","third","fourth","fifth","pupa")
for(j in 1:n) {
  SDC[1,j] = 1
  SDC[2,j] = ifelse(ADD[j+2,j] > 45,3,
                    ifelse(ADD[j+3,j] > 45,4,
                           5))
  SDC[3,j] = ifelse(ADD[j+3,j] > 77.3,4,
                    ifelse(ADD[j+4,j] > 77.3,5,
                           ifelse(ADD[j+5,j] > 77.3,6,
                                  7)))
  SDC[4,j] = ifelse(ADD[j+4,j] > 105.1,5,
                    ifelse(ADD[j+5,j] > 105.1,6,
                           ifelse(ADD[j+6,j] > 105.1,7,
                                  ifelse(ADD[j+7,j] > 105.1,8,
                                         9))))
  SDC[5,j] = ifelse(ADD[j+5,j] > 129.6,6,
                    ifelse(ADD[j+6,j] > 129.6,7,
                           ifelse(ADD[j+7,j] > 129.6,8,
                                  ifelse(ADD[j+8,j] > 129.6,9,
                                         10))))
  SDC[6,j] = ifelse(ADD[j+6,j] > 165.3,7,
                    ifelse(ADD[j+7,j] > 165.3,8,
                           ifelse(ADD[j+8,j] > 165.3,9,
                                  ifelse(ADD[j+9,j] > 165.3,10,
                                         ifelse(ADD[j+10,j] > 165.3,11,
                                                ifelse(ADD[j+11,j] > 165.3,12,
                                                       13))))))
  SDC[7,j] = ifelse(ADD[j+9,j] > 231.9,10,
                    ifelse(ADD[j+10,j] > 231.9,11,
                           ifelse(ADD[j+11,j] > 231.9,12,
                                  ifelse(ADD[j+12,j] > 231.9,13,
                                         ifelse(ADD[j+13,j] > 231.9,14,
                                                ifelse(ADD[j+14,j] > 231.9,15,
                                                       ifelse(ADD[j+15,j] > 231.9,16,
                                                              17)))))))
}


#########  calculate mean stage durations  ###############
#vector of mean stage durations
MSD = vector(length = 6)
#cohorts that finished 5th instar
I5 = ADD[(n+1),] > 231.9
#stage durations of cohorts that reached pupation
SDC[7,I5]
#fifth instar stage durations for each cohort
SD5 = SDC[7,I5]-SDC[6,I5]
MSD[6] = mean(SD5)
#cohorts that finished 4th instar
I4 = ADD[(n+1),] > 165.3
MSD[5] = mean(SDC[6,I4]-SDC[5,I4])
#cohorts that finished 3rd instar
I3 = ADD[(n+1),] > 129.6
MSD[4] = mean(SDC[5,I3]-SDC[4,I3])
#cohorts that finished 2nd instar
I2 = ADD[(n+1),] > 105.1
MSD[3] = mean(SDC[4,I2]-SDC[3,I2])
#cohorts that finished 1st instar
I1 = ADD[(n+1),] > 77.3
MSD[2] = mean(SDC[3,I1]-SDC[2,I1])
#cohorts that hatched
I0 = ADD[(n+1),] > 45
MSD[1] = mean(SDC[2,I0]-SDC[1,I0])
MSD

#calculate stage survival probabilities and cumulative survival probability
S1 = S[1]^MSD[1] 
S2 = S[2]^MSD[2]
S3 = S[3]^MSD[3]
S4 = S[4]^MSD[4]
S5 = S[5]^MSD[5]
S6 = S[6]^MSD[6]
SC = S1*S2*S3*S4*S5*S6 #cumulative survival probability
SN = S1*S2*S3*S4*S5 #survival probability most similar to Nail et al. 2015 estimator
SC

#define final count matrix
CD = matrix(nrow = n, ncol = 6)

#matrix of eggs
M1 = matrix(nrow=n+pupamax, ncol=n); colnames(M1) = 1:n
#matrix of 1st instar,etc
M2 = matrix(nrow=n+pupamax, ncol=n); colnames(M2) = 1:n
M3 = matrix(nrow=n+pupamax, ncol=n); colnames(M3) = 1:n
M4 = matrix(nrow=n+pupamax, ncol=n); colnames(M4) = 1:n
M5 = matrix(nrow=n+pupamax, ncol=n); colnames(M5) = 1:n
M6 = matrix(nrow=n+pupamax, ncol=n); colnames(M6) = 1:n


#############################################
#this code is different than PEDD and NEDD, because why? I think I decided to try a different way
#or it helps when adding in individual variation?

#ALL COHORTS OF EGGS
#initialize each cohort
for (i in 1:n){
  M1[i,i] = L[i]
}

#vector of current pop
N = L

#j are cols/cohorts, i are rows
for (j in 1:n) {
  for (i in (j+1):n){
    if (ADD[i,j] < 45) {
      M1[i,j] = N[j]*S[1]
      N[j] = M1[i,j]
    }
  }
}
##############################################
#ALL COHORTS OF FIRST INSTARS
#k counter keeps first survival rate as S[1], becuase the assumption is that it was an egg up until it was found as 1st instar
k = 0
for (j in 1:n) {
  #reset index for each column
  k = 0
  for (i in (j+1):n){
    if (ADD[i,j] >= 45 && ADD[i,j] < 77.3) {
      k = k+1
      #first survival rate is egg survival rate, because that's the survival rate previous to first detection of 1st instar
      if (k == 1){
        M2[i,j] = N[j]*S[1]
        N[j] = M2[i,j]
      }
      #every subsequent survival rate is 1st instar rate
      if (k > 1){
        M2[i,j] = N[j]*S[2]
        N[j] = M2[i,j]
      }
    }
  }
}
##################################
#ALL COHORTS OF SECOND INSTARS
k = 0
for (j in 1:n) {
  #reset index for each column
  k = 0
  for (i in (j+1):n){
    if (ADD[i,j] >= 77.3 && ADD[i,j] < 105.1) {
      k = k+1
      #first survival rate is egg survival rate, because that's the survival rate previous to first detection of 1st instar
      if (k == 1){
        M3[i,j] = N[j]*S[2]
        N[j] = M3[i,j]
      }
      #every subsequent survival rate is 1st instar rate
      if (k > 1){
        M3[i,j] = N[j]*S[3]
        N[j] = M3[i,j]
      }
    }
  }
}
###############################
#ALL COHORTS OF THIRD INSTARS
k = 0
for (j in 1:n) {
  #reset index for each column
  k = 0
  for (i in (j+1):n){
    if (ADD[i,j] >= 105.1 && ADD[i,j] < 129.6) {
      k = k+1
      #first survival rate is egg survival rate, because that's the survival rate previous to first detection of 1st instar
      if (k == 1){
        M4[i,j] = N[j]*S[3]
        N[j] = M4[i,j]
      }
      #every subsequent survival rate is 1st instar rate
      if (k > 1){
        M4[i,j] = N[j]*S[4]
        N[j] = M4[i,j]
      }
    }
  }
}
###############################
#ALL COHORTS OF FOURTH INSTARS
k = 0
for (j in 1:n) {
  #reset index for each column
  k = 0
  for (i in (j+1):n){
    if (ADD[i,j] >= 129.6 && ADD[i,j] < 165.3) {
      k = k+1
      #first survival rate is egg survival rate, because that's the survival rate previous to first detection of 1st instar
      if (k == 1){
        M5[i,j] = N[j]*S[4]
        N[j] = M5[i,j]
      }
      #every subsequent survival rate is 1st instar rate
      if (k > 1){
        M5[i,j] = N[j]*S[5]
        N[j] = M5[i,j]
      }
    }
  }
}
###############################
#ALL COHORTS OF FIFTH INSTARS
k = 0
for (j in 1:n) {
  #reset index for each column
  k = 0
  for (i in (j+1):n){
    if (ADD[i,j] >= 165.3 && ADD[i,j] < 231.9) {
      k = k+1
      #first survival rate is egg survival rate, because that's the survival rate previous to first detection of 1st instar
      if (k == 1){
        M6[i,j] = N[j]*S[5]
        N[j] = M6[i,j]
      }
      #every subsequent survival rate is 1st instar rate
      if (k > 1){
        M6[i,j] = N[j]*S[6]
        N[j] = M6[i,j]
      }
    }
  }
}
###############################

#sum cohort contributions for final count matrix
#sum eggs
for (i in 1:n){
  CD[i,1] = sum(M1[i,], na.rm = TRUE)
}
#sum first instars
for (i in 1:n){
  CD[i,2] = sum(M2[i,], na.rm = TRUE)
}
#sum second instars
for (i in 1:n){
  CD[i,3] = sum(M3[i,], na.rm = TRUE)
}
#sum third instars
for (i in 1:n){
  CD[i,4] = sum(M4[i,], na.rm = TRUE)
}
#sum fourth instars
for (i in 1:n){
  CD[i,5] = sum(M5[i,], na.rm = TRUE)
}
#sum fifth instars
for (i in 1:n){
  CD[i,6] = sum(M6[i,], na.rm = TRUE)
}
CD

#simulate 100 scenarios with Poisson random variation

#add pupamax buffer to beginning of count matrix
Buffer = matrix(0, nrow = pupamax, ncol = 6)
NDDB = rbind(Buffer,CD)
NDDB

#scenario 1
#Simulate Y_ij as from a Poisson distribution with mean l_ij
#100 simulated count matrices
C4.3 = list()
for (i in 1:100){
  C4.3[[i]] = apply(NDDB, MARGIN = c(1,2), FUN = function(x) rpois(1,x))
}


#save as study design sim 1
SDSC1 = C4.1 #max 50 eggs per day
SDSC2 = C4.2 #max 25 eggs per day
SDSC3 = C4.3 #max 100 eggs per day

save(SDSC1,SDSC2,SDSC3, file = "SDSC1-3.RData") #subsampled in other R script for different sampling schemes
#renamed file in Windows Explorer to Scen4.RData








###############################  Simulate data for Scenario 5  ###################
library(msm) #for rtnorm function

#daily survival probabilities
S = c(0.53,0.68,0.90,0.94,0.97,0.96)

n = 60
pupamax = 16

#eggs laid per day - normal distr with mean = 30 and SD = 15
#max 50 eggs laid per day - SDSC.50
CF = 50/dnorm(30, mean = 30, sd = 15) #correction factor to make density at 30 = 50
L = round(CF*dnorm(1:60, mean = 30, sd = 15))
#L5.1 = L
#L5.4 = L
L5.7 = L

#max 25 eggs laid per day - SDSC.25
CF = 25/dnorm(30, mean = 30, sd = 15) #correction factor to make density at 30 = 25
L = round(CF*dnorm(1:60, mean = 30, sd = 15))
#L5.2 = L
#L5.5 = L
L5.8 = L

#max 100 eggs laid per day - SDSC.100
CF = 100/dnorm(30, mean = 30, sd = 15) #correction factor to make density at 30 = 100
L = round(CF*dnorm(1:60, mean = 30, sd = 15))
#L5.3 = L
#L5.6 = L
L5.9 = L

#lists for C matrices
Cs5.1 = list() #max 50
Cs5.2 = list() #max 25
Cs5.3 = list() #max 100

#large loop to create many simulated datasets - closing bracket at ~ L709
for (m in 1:100){

  #random survival probabilites based on Zalucki 1982 Table 2 field data. Rows are individuals, cols are stages
#these are day degrees that it takes each individual to transition to the next stage
SDI = matrix(nrow = sum(L), ncol = 6)
SDI[,1] = rtnorm(sum(L), mean = 44.5, sd = 4.32, lower = 44.5-2*4.32, upper = 44.5+2*4.32)
SDI[,2] = rtnorm(sum(L), mean = 30.9, sd = 3.97, lower = 30.9-2*3.97, upper = 30.9+2*3.97)
SDI[,3] = rtnorm(sum(L), mean = 27.2, sd = 3.71, lower = 27.2-2*3.71, upper = 27.2+2*3.71)
SDI[,4] = rtnorm(sum(L), mean = 25.6, sd = 2.88, lower = 25.6-2*2.88, upper = 25.6+2*2.88)
SDI[,5] = rtnorm(sum(L), mean = 33.8, sd = 4.09, lower = 33.8-2*4.09, upper = 33.8+2*4.09)
SDI[,6] = rtnorm(sum(L), mean = 58.5, sd = 3.89, lower = 58.5-2*3.89, upper = 58.5+2*3.89)

CH = sum(L != 0)

#define matrices to put results into
#matrix for egg stage - rows are individuals, columns are days
MI1 = matrix(nrow=sum(L), ncol=n+pupamax); colnames(MI1) = 1:(n+pupamax)
#matrix for 1st instar and so on
MI2 = matrix(nrow=sum(L), ncol=n+pupamax); colnames(MI2) = 1:(n+pupamax)
MI3 = matrix(nrow=sum(L), ncol=n+pupamax); colnames(MI3) = 1:(n+pupamax)
MI4 = matrix(nrow=sum(L), ncol=n+pupamax); colnames(MI4) = 1:(n+pupamax)
MI5 = matrix(nrow=sum(L), ncol=n+pupamax); colnames(MI5) = 1:(n+pupamax)
MI6 = matrix(nrow=sum(L), ncol=n+pupamax); colnames(MI6) = 1:(n+pupamax)

#matrix of individual ending I for each stage, to carry over to each stage from previous stage
IM = matrix(nrow=sum(L), ncol=6)
#matrix of individual ending i (day/col of MI) for each stage, to carry over to each stage from previous stage
iM = matrix(nrow=sum(L), ncol=6)

#final counts of each stage each day
C = matrix(nrow = n+pupamax, ncol = 6)


MT = rep(21,n+pupamax); length(MT) 
#standard ADD, SDC, and MSD code for temps:

#calculate accumulated day-degrees for each cohort - need to add pupamax to rows and cols - days and cohorts, respectively
ADD = matrix(nrow=n+pupamax, ncol=n+pupamax); colnames(ADD) = 1:(n+pupamax)
for (j in 1:n){
  for (i in (j+1):(n+1)) {
    ADD[i,j] = sum(MT[j:(i-1)])
  }
}
#set upper right off diagonal cells to 0 - these have to be set to 0 for JAGS (this code was originally in the model)
for (j in 1:n) {
  for (i in 1:j){
    ADD[i,j] = 0    
  }
}
#set lower cells to 1000 so doesn't keep while loop going
for (j in 1:n) {
  for (i in (n+2):(n+pupamax)){
    ADD[i,j] = 1000    
  }
}
#set right buffer cells to 1000 so doesn't keep while loop going
for (j in (n+1):(n+pupamax)) {
  for (i in 1:(n+pupamax)){
    ADD[i,j] = 1000    
  }
}

############ Eggs  ##################
Inds = 1
#k loops through cohorts
for (k in 1:CH){
  #set initial value to 1
  MI1[Inds:(Inds+L[k]-1),k] = 1
  #j loops through individuals of a cohort
  for (j in 1:L[k]){
    I=1
    i=1+k #i is the col/day index
    while (ADD[i,k] < SDI[Inds,1]) {
      R = runif(1)
      if (R > S[1]){
        I = 0
      }
      MI1[Inds,i] = I #kind of confusing because days are rows in MI1 but days are cols in ADD
      #save final I state
      IM[Inds,1] = I
      i=i+1
      #save final i index
      iM[Inds,1] = i
    }
    Inds = Inds + 1
  }
}


########  1st instar  ###############
Inds = 1
#k loops through cohorts
for (k in 1:CH){
  #j loops through individuals of a cohort
  for (j in 1:L[k]){
    #take ending I from previous stage
    I=IM[Inds,1]
    i=iM[Inds,1] #i is the col/day index
    #index for no of times it has been through while loop
    l = 1
    while (ADD[i,k] < (SDI[Inds,1]+SDI[Inds,2])) {
      R = runif(1)
      #survival before first detection as first instar
      if (l == 1) {
        if (R > S[1]){
          I = 0
        }
      }
      #survival for following days
      if (l > 1){
        if (R > S[2]){
          I = 0
        }
      }
      MI2[Inds,i] = I #kind of confusing because days are rows in MI1 but days are cols in ADD
      #save final I state
      IM[Inds,2] = I
      i=i+1
      iM[Inds,2] = i
      l=l+1
    }
    
    if (ADD[i,k] >= (SDI[Inds,1]+SDI[Inds,2])) { #this if loop takes care of IM and iM if 1st instar skipped because its warm
      IM[Inds,2] = I
      iM[Inds,2] = i
    }
    
    
    Inds = Inds + 1
  }
}

########  2nd instar  ###############
Inds = 1
#k loops through cohorts
for (k in 1:CH){
  #j loops through individuals of a cohort
  for (j in 1:L[k]){
    #take ending I from previous stage
    I=IM[Inds,2]
    i=iM[Inds,2] #i is the col/day index
    #index for no of times it has been through while loop
    l = 1
    while (ADD[i,k] < (SDI[Inds,1]+SDI[Inds,2]+SDI[Inds,3])) {
      R = runif(1)
      #survival before first detection as first instar
      if (l == 1) {
        if (R > S[2]){
          I = 0
        }
      }
      #survival for following days
      if (l > 1){
        if (R > S[3]){
          I = 0
        }
      }
      MI3[Inds,i] = I #kind of confusing because days are rows in MI1 but days are cols in ADD
      #save final I state
      IM[Inds,3] = I
      i=i+1
      iM[Inds,3] = i
      l=l+1
    }
    
    if (ADD[i,k] >= (SDI[Inds,1]+SDI[Inds,2]+SDI[Inds,3])) { #this if loop takes care of IM and iM if 1st instar skipped because its warm
      IM[Inds,3] = I
      iM[Inds,3] = i
    }
    
    
    Inds = Inds + 1
  }
}


########  3rd instar  ###############
Inds = 1
#k loops through cohorts
for (k in 1:CH){
  #j loops through individuals of a cohort
  for (j in 1:L[k]){
    #take ending I from previous stage
    I=IM[Inds,3]
    i=iM[Inds,3] #i is the col/day index
    #index for no of times it has been through while loop
    l = 1
    while (ADD[i,k] < (SDI[Inds,1]+SDI[Inds,2]+SDI[Inds,3]+SDI[Inds,4])) {
      R = runif(1)
      #survival before first detection as first instar
      if (l == 1) {
        if (R > S[3]){
          I = 0
        }
      }
      #survival for following days
      if (l > 1){
        if (R > S[4]){
          I = 0
        }
      }
      MI4[Inds,i] = I #kind of confusing because days are rows in MI1 but days are cols in ADD
      #save final I state
      IM[Inds,4] = I
      i=i+1
      iM[Inds,4] = i
      l=l+1
    }
    
    if (ADD[i,k] >= (SDI[Inds,1]+SDI[Inds,2]+SDI[Inds,3]+SDI[Inds,4])) { #this if loop takes care of IM and iM if 1st instar skipped because its warm
      IM[Inds,4] = I
      iM[Inds,4] = i
    }
    
    Inds = Inds + 1
  }
}

########  4th instar  ###############
Inds = 1
#k loops through cohorts
for (k in 1:CH){
  #j loops through individuals of a cohort
  for (j in 1:L[k]){
    #take ending I from previous stage
    I=IM[Inds,4]
    i=iM[Inds,4] #i is the col/day index
    #index for no of times it has been through while loop
    l = 1
    while (ADD[i,k] < (SDI[Inds,1]+SDI[Inds,2]+SDI[Inds,3]+SDI[Inds,4]+SDI[Inds,5])) {
      R = runif(1)
      #survival before first detection as first instar
      if (l == 1) {
        if (R > S[4]){
          I = 0
        }
      }
      #survival for following days
      if (l > 1){
        if (R > S[5]){
          I = 0
        }
      }
      MI5[Inds,i] = I #kind of confusing because days are rows in MI1 but days are cols in ADD
      #save final I state
      IM[Inds,5] = I
      i=i+1
      iM[Inds,5] = i
      l=l+1
    }
    
    if (ADD[i,k] >= (SDI[Inds,1]+SDI[Inds,2]+SDI[Inds,3]+SDI[Inds,4]+SDI[Inds,5])) { #this if loop takes care of IM and iM if 1st instar skipped because its warm
      IM[Inds,5] = I
      iM[Inds,5] = i
    }
    
    Inds = Inds + 1
  }
}

########  5th instar  ###############
Inds = 1
#k loops through cohorts
for (k in 1:CH){
  #j loops through individuals of a cohort
  for (j in 1:L[k]){
    #take ending I from previous stage
    I=IM[Inds,5]
    i=iM[Inds,5] #i is the col/day index
    #index for no of times it has been through while loop
    l = 1
    while (ADD[i,k] < (SDI[Inds,1]+SDI[Inds,2]+SDI[Inds,3]+SDI[Inds,4]+SDI[Inds,5]+SDI[Inds,6])) {
      R = runif(1)
      #survival before first detection as first instar
      if (l == 1) {
        if (R > S[5]){
          I = 0
        }
      }
      #survival for following days
      if (l > 1){
        if (R > S[6]){
          I = 0
        }
      }
      MI6[Inds,i] = I #kind of confusing because days are rows in MI1 but days are cols in ADD
      #save final I state
      IM[Inds,6] = I
      i=i+1
      iM[Inds,6] = i
      l=l+1
    }
    Inds = Inds + 1
  }
}


#sum columns to get eggs per day
#NS1 means summed N eggs
NS1 = vector(length = n+pupamax)
for (i in 1:n){
  NS1[i] = sum(MI1[,i], na.rm = TRUE)
}
NS1
NS2 = vector(length = n+pupamax)
for (i in 1:n){
  NS2[i] = sum(MI2[,i], na.rm = TRUE)
}
NS2
NS3 = vector(length = n+pupamax)
for (i in 1:n){
  NS3[i] = sum(MI3[,i], na.rm = TRUE)
}
NS3
NS4 = vector(length = n+pupamax)
for (i in 1:n){
  NS4[i] = sum(MI4[,i], na.rm = TRUE)
}
NS4
NS5 = vector(length = n+pupamax)
for (i in 1:n){
  NS5[i] = sum(MI5[,i], na.rm = TRUE)
}
NS5
NS6 = vector(length = n+pupamax)
for (i in 1:n){
  NS6[i] = sum(MI6[,i], na.rm = TRUE)
}
NS6

#put sums into count matrix in analysis format

C[,1] = NS1
C[,2] = NS2
C[,3] = NS3
C[,4] = NS4
C[,5] = NS5
C[,6] = NS6
C

C = C[1:n,]
C

Buffer = matrix(0, nrow = pupamax, ncol = 6)
C = rbind(Buffer,C)
C

#apply random Poisson variation and add to list of datasets
#Cs5.1[[m]] = apply(C, MARGIN = c(1,2), FUN = function(x) rpois(1,x))
#Cs5.2[[m]] = apply(C, MARGIN = c(1,2), FUN = function(x) rpois(1,x))
Cs5.3[[m]] = apply(C, MARGIN = c(1,2), FUN = function(x) rpois(1,x))

}

save(Cs5.1,Cs5.2,Cs5.3, file = "Scen5Counts.RData") #subsampled in other R script for different sampling schemes



######  Analyze Scenario 4 and Scenario 5  #######################
library(runjags)
runjags.options(force.summary=TRUE)

pupamax = 16
n = 60+pupamax #adding pupamax to n is not necessary, but doesn't hurt, for scenario 1. But pupamax is added to all sims for simplicity.  
#also note that in sim code I don't combine pupamax with n, though maybe I should have
MT = rep(21,n)

#calculate accumulated day-degrees for each cohort
ADD = matrix(nrow=n+pupamax, ncol=n); colnames(ADD) = 1:n
for (j in 1:n){
  for (i in (j+1):(n+1)) {
    ADD[i,j] = sum(MT[j:(i-1)])
  }
}
#set upper right off diagonal cells to 0 - these have to be set to 0 for JAGS (this code was originally in the model)
for (j in 1:n) {
  for (i in 1:j){
    ADD[i,j] = 0    
  }
}
#set lower cells to 0
for (j in 1:n) {
  for (i in (n+2):(n+pupamax)){
    ADD[i,j] = 0    
  }
}

#calculate stage durations for each cohort (this code was originally in the JAGS model, which is why it is a string of ifelse() statements instead of something more elegant)
SDC = matrix(nrow = 7, ncol = n); colnames(SDC) = 1:n; rownames(SDC) = c("egg","first","second","third","fourth","fifth","pupa")
for(j in 1:n) {
  SDC[1,j] = 1
  SDC[2,j] = ifelse(ADD[j+2,j] > 45,3,
                    ifelse(ADD[j+3,j] > 45,4,
                           5))
  SDC[3,j] = ifelse(ADD[j+3,j] > 77.3,4,
                    ifelse(ADD[j+4,j] > 77.3,5,
                           ifelse(ADD[j+5,j] > 77.3,6,
                                  7)))
  SDC[4,j] = ifelse(ADD[j+4,j] > 105.1,5,
                    ifelse(ADD[j+5,j] > 105.1,6,
                           ifelse(ADD[j+6,j] > 105.1,7,
                                  ifelse(ADD[j+7,j] > 105.1,8,
                                         9))))
  SDC[5,j] = ifelse(ADD[j+5,j] > 129.6,6,
                    ifelse(ADD[j+6,j] > 129.6,7,
                           ifelse(ADD[j+7,j] > 129.6,8,
                                  ifelse(ADD[j+8,j] > 129.6,9,
                                         10))))
  SDC[6,j] = ifelse(ADD[j+6,j] > 165.3,7,
                    ifelse(ADD[j+7,j] > 165.3,8,
                           ifelse(ADD[j+8,j] > 165.3,9,
                                  ifelse(ADD[j+9,j] > 165.3,10,
                                         ifelse(ADD[j+10,j] > 165.3,11,
                                                ifelse(ADD[j+11,j] > 165.3,12,
                                                       13))))))
  SDC[7,j] = ifelse(ADD[j+9,j] > 231.9,10,
                    ifelse(ADD[j+10,j] > 231.9,11,
                           ifelse(ADD[j+11,j] > 231.9,12,
                                  ifelse(ADD[j+12,j] > 231.9,13,
                                         ifelse(ADD[j+13,j] > 231.9,14,
                                                ifelse(ADD[j+14,j] > 231.9,15,
                                                       ifelse(ADD[j+15,j] > 231.9,16,
                                                              17)))))))
}



#########  calculate mean stage durations  ###############
#vector of mean stage durations
MSD = vector(length = 6)
#cohorts that finished 5th instar
I5 = ADD[(n+1),] > 231.9
#stage durations of cohorts that reached pupation
SDC[7,I5]
#fifth instar stage durations for each cohort
SD5 = SDC[7,I5]-SDC[6,I5]
MSD[6] = mean(SD5)
#cohorts that finished 4th instar
I4 = ADD[(n+1),] > 165.3
MSD[5] = mean(SDC[6,I4]-SDC[5,I4])
#cohorts that finished 3rd instar
I3 = ADD[(n+1),] > 129.6
MSD[4] = mean(SDC[5,I3]-SDC[4,I3])
#cohorts that finished 2nd instar
I2 = ADD[(n+1),] > 105.1
MSD[3] = mean(SDC[4,I2]-SDC[3,I2])
#cohorts that finished 1st instar
I1 = ADD[(n+1),] > 77.3
MSD[2] = mean(SDC[3,I1]-SDC[2,I1])
#cohorts that hatched
I0 = ADD[(n+1),] > 45
MSD[1] = mean(SDC[2,I0]-SDC[1,I0])
MSD

#calculate stage survival probabilities and cumulative survival probability
S1 = S[1]^MSD[1] 
S2 = S[2]^MSD[2]
S3 = S[3]^MSD[3]
S4 = S[4]^MSD[4]
S5 = S[5]^MSD[5]
S6 = S[6]^MSD[6]
SC = S1*S2*S3*S4*S5*S6 #cumulative survival probability
SN = S1*S2*S3*S4*S5 #survival probability most similar to Nail et al. 2015 estimator
SC


#specify initial values for three chains
inits = list(list(lambda2 = rep(50, n), S = rep(0.8,6)),
             list(lambda2 = rep(5, n), S = rep(0.5,6)),
             list(lambda2 = rep(100, n), S = rep(0.2,6)))

#list for results
R4.1 = list()
R4.2 = list()
R4.3 = list()
R5.1 = list()
R5.2 = list()
R5.3 = list()

#list for time
time = list()



#Scenario 4.1 - Daily
time[[1]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=C4.1[[i]])
    #call run.jags
    R4.1[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    i
  }
)


Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L4.1,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R4.1[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R4.1[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R4.1[[1]]$summary$quantile[,5]-R4.1[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R4.1[[i]]$summary$quantile[,5]-R4.1[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
write.csv(Res,"Scen4_1.csv")

save(R4.1,Res,file = "Scen4.1.Rdata")



#Scenario 4.2 - Daily
time[[2]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=C4.2[[i]])
    #call run.jags
    R4.2[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    i
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L4.2,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R4.2[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R4.2[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R4.2[[1]]$summary$quantile[,5]-R4.2[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R4.2[[i]]$summary$quantile[,5]-R4.2[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
write.csv(Res,"Scen4_2.csv")

save(R4.2,Res,file = "Scen4.2.Rdata")



#Scenario 4.3 - Daily
time[[3]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=C4.3[[i]])
    #call run.jags
    R4.3[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    i
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L4.3,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R4.3[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R4.3[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R4.3[[1]]$summary$quantile[,5]-R4.3[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R4.3[[i]]$summary$quantile[,5]-R4.3[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
write.csv(Res,"Scen4_3.csv")


save(R4.3,Res,file = "Scen4.3.Rdata")
load(file = "Scen4.1.Rdata")
save(time,file = "simtime.Rdata")
load(file = "simtime.scen4.Rdata")


#############  Scenario 5  #####################

#list for time in scen 5 sims
time2 = list()

#Scenario 5.1 - Daily
time2[[4]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=Cs5.1[[i]])
    #call run.jags
    R5.1[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(i)
  }
)


Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L5.1,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R5.1[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R5.1[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R5.1[[1]]$summary$quantile[,5]-R5.1[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R5.1[[i]]$summary$quantile[,5]-R5.1[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
write.csv(Res,"Scen5_1.csv")

save(R5.1,Res,file = "Scen5.1.Rdata")


#Scenario 5.2 - Daily
time2[[5]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=Cs5.2[[i]])
    #call run.jags
    R5.2[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(i)
  }
)


Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L5.2,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R5.2[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R5.2[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R5.2[[1]]$summary$quantile[,5]-R5.2[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R5.2[[i]]$summary$quantile[,5]-R5.2[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
write.csv(Res,"Scen5_2.csv")

save(R5.2,Res,file = "Scen5.2.Rdata")


#Scenario 5.3 - Daily
time2[[6]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=Cs5.3[[i]])
    #call run.jags
    R5.3[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(i)
  }
)


Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L5.3,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R5.3[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R5.3[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R5.3[[1]]$summary$quantile[,5]-R5.3[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R5.3[[i]]$summary$quantile[,5]-R5.3[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
write.csv(Res,"Scen5_3.csv")

save(R5.3,Res,file = "Scen5.3.Rdata")

############
#list of results - 1:6 are daily sampling, 1:3 are scen 4, 4:6 are scen 5, 

ResAll = list()
ResAll[[4]] = Res

save(ResAll,file = "ResAll.Rdata")


############## Other Sampling Schemes besides Daily ##########################

#load scenario 4 datasets
load(file = "SDSC1-3.RData")
SDSC1[[1]]

#subset the data

########################### twice weekly counts ##################################
#SDSC4 = 50 max eggs, twice weekly sampling, no individual variation
SDSC4 = SDSC1
SDSC4[[1]][1:16,]=NA
SDSC4[[1]][17:18,] = NA; SDSC4[[1]][20:21,] = NA; SDSC4[[1]][23:25,] = NA; SDSC4[[1]][27:28,] = NA; SDSC4[[1]][30:32,] = NA; SDSC4[[1]][34:35,] = NA; SDSC4[[1]][37:39,] = NA; 
SDSC4[[1]][41:42,] = NA; SDSC4[[1]][44:46,] = NA; SDSC4[[1]][48:49,] = NA; SDSC4[[1]][51:53,] = NA; SDSC4[[1]][55:56,] = NA; SDSC4[[1]][58:60,] = NA; SDSC4[[1]][62:63,] = NA; 
SDSC4[[1]][65:67,] = NA; SDSC4[[1]][69:70,] = NA; SDSC4[[1]][72:74,] = NA; SDSC4[[1]][76,] = NA
SDSC4[[1]]

for (i in 1:100){
  SDSC4[[i]][1:16,]=NA
  SDSC4[[i]][17:18,] = NA; SDSC4[[i]][20:21,] = NA; SDSC4[[i]][23:25,] = NA; SDSC4[[i]][27:28,] = NA; SDSC4[[i]][30:32,] = NA; SDSC4[[i]][34:35,] = NA; SDSC4[[i]][37:39,] = NA; 
  SDSC4[[i]][41:42,] = NA; SDSC4[[i]][44:46,] = NA; SDSC4[[i]][48:49,] = NA; SDSC4[[i]][51:53,] = NA; SDSC4[[i]][55:56,] = NA; SDSC4[[i]][58:60,] = NA; SDSC4[[i]][62:63,] = NA; 
  SDSC4[[i]][65:67,] = NA; SDSC4[[i]][69:70,] = NA; SDSC4[[i]][72:74,] = NA; SDSC4[[i]][76,] = NA
}
SDSC4[[4]]

#SDSC5 - 25 max eggs, twice weekly sampling, no individual variation
SDSC5 = SDSC2
for (i in 1:100){
  SDSC5[[i]][1:16,]=NA
  SDSC5[[i]][17:18,] = NA; SDSC5[[i]][20:21,] = NA; SDSC5[[i]][23:25,] = NA; SDSC5[[i]][27:28,] = NA; SDSC5[[i]][30:32,] = NA; SDSC5[[i]][34:35,] = NA; SDSC5[[i]][37:39,] = NA; 
  SDSC5[[i]][41:42,] = NA; SDSC5[[i]][44:46,] = NA; SDSC5[[i]][48:49,] = NA; SDSC5[[i]][51:53,] = NA; SDSC5[[i]][55:56,] = NA; SDSC5[[i]][58:60,] = NA; SDSC5[[i]][62:63,] = NA; 
  SDSC5[[i]][65:67,] = NA; SDSC5[[i]][69:70,] = NA; SDSC5[[i]][72:74,] = NA; SDSC5[[i]][76,] = NA
}
SDSC5[[4]]


#SCSC6 - 100 max eggs, twice weekly sampling, no individual variation
SDSC6 = SDSC3
for (i in 1:100){
  SDSC6[[i]][1:16,]=NA
  SDSC6[[i]][17:18,] = NA; SDSC6[[i]][20:21,] = NA; SDSC6[[i]][23:25,] = NA; SDSC6[[i]][27:28,] = NA; SDSC6[[i]][30:32,] = NA; SDSC6[[i]][34:35,] = NA; SDSC6[[i]][37:39,] = NA; 
  SDSC6[[i]][41:42,] = NA; SDSC6[[i]][44:46,] = NA; SDSC6[[i]][48:49,] = NA; SDSC6[[i]][51:53,] = NA; SDSC6[[i]][55:56,] = NA; SDSC6[[i]][58:60,] = NA; SDSC6[[i]][62:63,] = NA; 
  SDSC6[[i]][65:67,] = NA; SDSC6[[i]][69:70,] = NA; SDSC6[[i]][72:74,] = NA; SDSC6[[i]][76,] = NA
}
SDSC6[[4]]

#Cs5.4 - 50 max eggs, twice weekly sampling, with ind var
Cs5.4 = Cs5.1
for (i in 1:100){
  Cs5.4[[i]][1:16,]=NA
  Cs5.4[[i]][17:18,] = NA; Cs5.4[[i]][20:21,] = NA; Cs5.4[[i]][23:25,] = NA; Cs5.4[[i]][27:28,] = NA; Cs5.4[[i]][30:32,] = NA; Cs5.4[[i]][34:35,] = NA; Cs5.4[[i]][37:39,] = NA; 
  Cs5.4[[i]][41:42,] = NA; Cs5.4[[i]][44:46,] = NA; Cs5.4[[i]][48:49,] = NA; Cs5.4[[i]][51:53,] = NA; Cs5.4[[i]][55:56,] = NA; Cs5.4[[i]][58:60,] = NA; Cs5.4[[i]][62:63,] = NA; 
  Cs5.4[[i]][65:67,] = NA; Cs5.4[[i]][69:70,] = NA; Cs5.4[[i]][72:74,] = NA; Cs5.4[[i]][76,] = NA
}
Cs5.4[[4]]

#Cs5.5
Cs5.5 = Cs5.2
for (i in 1:100){
  Cs5.5[[i]][1:16,]=NA
  Cs5.5[[i]][17:18,] = NA; Cs5.5[[i]][20:21,] = NA; Cs5.5[[i]][23:25,] = NA; Cs5.5[[i]][27:28,] = NA; Cs5.5[[i]][30:32,] = NA; Cs5.5[[i]][34:35,] = NA; Cs5.5[[i]][37:39,] = NA; 
  Cs5.5[[i]][41:42,] = NA; Cs5.5[[i]][44:46,] = NA; Cs5.5[[i]][48:49,] = NA; Cs5.5[[i]][51:53,] = NA; Cs5.5[[i]][55:56,] = NA; Cs5.5[[i]][58:60,] = NA; Cs5.5[[i]][62:63,] = NA; 
  Cs5.5[[i]][65:67,] = NA; Cs5.5[[i]][69:70,] = NA; Cs5.5[[i]][72:74,] = NA; Cs5.5[[i]][76,] = NA
}
Cs5.5[[4]]

#Cs5.6
Cs5.6 = Cs5.3
for (i in 1:100){
  Cs5.6[[i]][1:16,]=NA
  Cs5.6[[i]][17:18,] = NA; Cs5.6[[i]][20:21,] = NA; Cs5.6[[i]][23:25,] = NA; Cs5.6[[i]][27:28,] = NA; Cs5.6[[i]][30:32,] = NA; Cs5.6[[i]][34:35,] = NA; Cs5.6[[i]][37:39,] = NA; 
  Cs5.6[[i]][41:42,] = NA; Cs5.6[[i]][44:46,] = NA; Cs5.6[[i]][48:49,] = NA; Cs5.6[[i]][51:53,] = NA; Cs5.6[[i]][55:56,] = NA; Cs5.6[[i]][58:60,] = NA; Cs5.6[[i]][62:63,] = NA; 
  Cs5.6[[i]][65:67,] = NA; Cs5.6[[i]][69:70,] = NA; Cs5.6[[i]][72:74,] = NA; Cs5.6[[i]][76,] = NA
}
Cs5.6[[4]]




##########  ANALYZE TWICE WEEKLY FOR SCEN 4-5 ################################################
#Rerun these to be safe
library(runjags)
runjags.options(force.summary=TRUE)

load(file = "ResAll.Rdata")

pupamax = 16
n = 60+pupamax #adding pupamax to n is not necessary, but doesn't hurt, for scenario 1. But pupamax is added to all sims for simplicity.  
#also note that in sim code I don't combine pupamax with n, though maybe I should have
MT = rep(21,n)

#calculate accumulated day-degrees for each cohort
ADD = matrix(nrow=n+pupamax, ncol=n); colnames(ADD) = 1:n
for (j in 1:n){
  for (i in (j+1):(n+1)) {
    ADD[i,j] = sum(MT[j:(i-1)])
  }
}
#set upper right off diagonal cells to 0 - these have to be set to 0 for JAGS (this code was originally in the model)
for (j in 1:n) {
  for (i in 1:j){
    ADD[i,j] = 0    
  }
}
#set lower cells to 0
for (j in 1:n) {
  for (i in (n+2):(n+pupamax)){
    ADD[i,j] = 0    
  }
}

#calculate stage durations for each cohort (this code was originally in the JAGS model, which is why it is a string of ifelse() statements instead of something more elegant)
SDC = matrix(nrow = 7, ncol = n); colnames(SDC) = 1:n; rownames(SDC) = c("egg","first","second","third","fourth","fifth","pupa")
for(j in 1:n) {
  SDC[1,j] = 1
  SDC[2,j] = ifelse(ADD[j+2,j] > 45,3,
                    ifelse(ADD[j+3,j] > 45,4,
                           5))
  SDC[3,j] = ifelse(ADD[j+3,j] > 77.3,4,
                    ifelse(ADD[j+4,j] > 77.3,5,
                           ifelse(ADD[j+5,j] > 77.3,6,
                                  7)))
  SDC[4,j] = ifelse(ADD[j+4,j] > 105.1,5,
                    ifelse(ADD[j+5,j] > 105.1,6,
                           ifelse(ADD[j+6,j] > 105.1,7,
                                  ifelse(ADD[j+7,j] > 105.1,8,
                                         9))))
  SDC[5,j] = ifelse(ADD[j+5,j] > 129.6,6,
                    ifelse(ADD[j+6,j] > 129.6,7,
                           ifelse(ADD[j+7,j] > 129.6,8,
                                  ifelse(ADD[j+8,j] > 129.6,9,
                                         10))))
  SDC[6,j] = ifelse(ADD[j+6,j] > 165.3,7,
                    ifelse(ADD[j+7,j] > 165.3,8,
                           ifelse(ADD[j+8,j] > 165.3,9,
                                  ifelse(ADD[j+9,j] > 165.3,10,
                                         ifelse(ADD[j+10,j] > 165.3,11,
                                                ifelse(ADD[j+11,j] > 165.3,12,
                                                       13))))))
  SDC[7,j] = ifelse(ADD[j+9,j] > 231.9,10,
                    ifelse(ADD[j+10,j] > 231.9,11,
                           ifelse(ADD[j+11,j] > 231.9,12,
                                  ifelse(ADD[j+12,j] > 231.9,13,
                                         ifelse(ADD[j+13,j] > 231.9,14,
                                                ifelse(ADD[j+14,j] > 231.9,15,
                                                       ifelse(ADD[j+15,j] > 231.9,16,
                                                              17)))))))
}



#########  calculate mean stage durations  ###############
#vector of mean stage durations
MSD = vector(length = 6)
#cohorts that finished 5th instar
I5 = ADD[(n+1),] > 231.9
#stage durations of cohorts that reached pupation
SDC[7,I5]
#fifth instar stage durations for each cohort
SD5 = SDC[7,I5]-SDC[6,I5]
MSD[6] = mean(SD5)
#cohorts that finished 4th instar
I4 = ADD[(n+1),] > 165.3
MSD[5] = mean(SDC[6,I4]-SDC[5,I4])
#cohorts that finished 3rd instar
I3 = ADD[(n+1),] > 129.6
MSD[4] = mean(SDC[5,I3]-SDC[4,I3])
#cohorts that finished 2nd instar
I2 = ADD[(n+1),] > 105.1
MSD[3] = mean(SDC[4,I2]-SDC[3,I2])
#cohorts that finished 1st instar
I1 = ADD[(n+1),] > 77.3
MSD[2] = mean(SDC[3,I1]-SDC[2,I1])
#cohorts that hatched
I0 = ADD[(n+1),] > 45
MSD[1] = mean(SDC[2,I0]-SDC[1,I0])
MSD

#calculate stage survival probabilities and cumulative survival probability
S1 = S[1]^MSD[1] 
S2 = S[2]^MSD[2]
S3 = S[3]^MSD[3]
S4 = S[4]^MSD[4]
S5 = S[5]^MSD[5]
S6 = S[6]^MSD[6]
SC = S1*S2*S3*S4*S5*S6 #cumulative survival probability
SN = S1*S2*S3*S4*S5 #survival probability most similar to Nail et al. 2015 estimator
SC


#specify initial values for three chains
inits = list(list(lambda2 = rep(50, n), S = rep(0.8,6)),
             list(lambda2 = rep(5, n), S = rep(0.5,6)),
             list(lambda2 = rep(100, n), S = rep(0.2,6)))

#list for results
R4.4 = list()
R4.5 = list()
R4.6 = list()
R5.4 = list()
R5.5 = list()
R5.6 = list()

#list for time
runtime = list()


#Scenario 4.4 - twice weekly
runtime[[7]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=SDSC4[[i]])
    #call run.jags
    R4.4[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(paste0("Simulation 4.4.",i," finished running."))
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L4.4,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R4.4[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R4.4[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R4.4[[1]]$summary$quantile[,5]-R4.4[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R4.4[[i]]$summary$quantile[,5]-R4.4[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
ResAll[[7]] = Res
write.csv(Res,"Scen4_4.csv")
save(R4.4,Res,file = "Scen4.4.Rdata")
rm(R4.4)

#Scen 4.5
runtime[[8]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=SDSC5[[i]])
    #call run.jags
    R4.5[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(paste0("Simulation ",i," finished running."))
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L4.5,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R4.5[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R4.5[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R4.5[[1]]$summary$quantile[,5]-R4.5[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R4.5[[i]]$summary$quantile[,5]-R4.5[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
ResAll[[8]] = Res
write.csv(Res,"Scen4.5.csv")
save(R4.5,Res,file = "Scen4.5.Rdata")
rm(R4.5)


#Scen 4.6
runtime[[9]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=SDSC6[[i]])
    #call run.jags
    R4.6[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(paste0("Simulation ",i," finished running."))
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L4.6,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R4.6[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R4.6[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R4.6[[1]]$summary$quantile[,5]-R4.6[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R4.6[[i]]$summary$quantile[,5]-R4.6[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
ResAll[[9]] = Res
write.csv(Res,"Scen4.6.csv")
save(R4.6,Res,file = "Scen4.6.Rdata")
rm(R4.6)


#Scen5.4
runtime[[10]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=Cs5.4[[i]])
    #call run.jags
    R5.4[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(paste0("Simulation ",i," finished running."))
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L5.4,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R5.4[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R5.4[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R5.4[[1]]$summary$quantile[,5]-R5.4[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R5.4[[i]]$summary$quantile[,5]-R5.4[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
ResAll[[10]] = Res
write.csv(Res,"Scen5.4.csv")
save(R5.4,Res,file = "Scen5.4.Rdata")
rm(R5.4)

#scen5.5
runtime[[11]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=Cs5.5[[i]])
    #call run.jags
    R5.5[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(paste0("Simulation ",i," finished running."))
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L5.5,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R5.5[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R5.5[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R5.5[[1]]$summary$quantile[,5]-R5.5[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R5.5[[i]]$summary$quantile[,5]-R5.5[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
ResAll[[11]] = Res
write.csv(Res,"Scen5.5.csv")
save(R5.5,Res,file = "Scen5.5.Rdata")
rm(R5.5)


#Scen5.6
runtime[[12]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=Cs5.6[[i]])
    #call run.jags
    R5.6[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(paste0("Simulation ",i," finished running."))
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L5.6,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R5.6[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R5.6[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R5.6[[1]]$summary$quantile[,5]-R5.6[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R5.6[[i]]$summary$quantile[,5]-R5.6[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
ResAll[[12]] = Res
write.csv(Res,"Scen5.6.csv")
save(R5.6,Res,file = "Scen5.6.Rdata")
rm(R5.6)

save(ResAll,file = "ResAll.Rdata")
save(runtime,file = "runtime.Rdata")



########################### once weekly counts ##################################

load("ResAll.Rdata")
load("runtime.Rdata")
runtime
load(file = "SDSC1-3.RData")
load("Scen5Counts.Rdata")
Cs5.1[[1]]

SDSC7 = SDSC1
for (i in 1:100){
  SDSC7[[i]][1:16,]=NA
  SDSC7[[i]][18:23,] = NA; SDSC7[[i]][25:30,] = NA; SDSC7[[i]][32:37,] = NA; SDSC7[[i]][39:44,] = NA; SDSC7[[i]][46:51,] = NA; 
  SDSC7[[i]][53:58,] = NA; SDSC7[[i]][60:65,] = NA; SDSC7[[i]][67:72,] = NA; SDSC7[[i]][74:76,] = NA; 
}
SDSC7[[4]]

SDSC8 = SDSC2
for (i in 1:100){
  SDSC8[[i]][1:16,]=NA
  SDSC8[[i]][18:23,] = NA; SDSC8[[i]][25:30,] = NA; SDSC8[[i]][32:37,] = NA; SDSC8[[i]][39:44,] = NA; SDSC8[[i]][46:51,] = NA; 
  SDSC8[[i]][53:58,] = NA; SDSC8[[i]][60:65,] = NA; SDSC8[[i]][67:72,] = NA; SDSC8[[i]][74:76,] = NA; 
}
SDSC8[[4]]

SDSC9 = SDSC3
for (i in 1:100){
  SDSC9[[i]][1:16,]=NA
  SDSC9[[i]][18:23,] = NA; SDSC9[[i]][25:30,] = NA; SDSC9[[i]][32:37,] = NA; SDSC9[[i]][39:44,] = NA; SDSC9[[i]][46:51,] = NA; 
  SDSC9[[i]][53:58,] = NA; SDSC9[[i]][60:65,] = NA; SDSC9[[i]][67:72,] = NA; SDSC9[[i]][74:76,] = NA; 
}
SDSC9[[4]]

Cs5.7 = Cs5.1
for (i in 1:100){
  Cs5.7[[i]][1:16,]=NA
  Cs5.7[[i]][18:23,] = NA; Cs5.7[[i]][25:30,] = NA; Cs5.7[[i]][32:37,] = NA; Cs5.7[[i]][39:44,] = NA; Cs5.7[[i]][46:51,] = NA; 
  Cs5.7[[i]][53:58,] = NA; Cs5.7[[i]][60:65,] = NA; Cs5.7[[i]][67:72,] = NA; Cs5.7[[i]][74:76,] = NA; 
}
Cs5.7[[4]]

Cs5.8 = Cs5.2
for (i in 1:100){
  Cs5.8[[i]][1:16,]=NA
  Cs5.8[[i]][18:23,] = NA; Cs5.8[[i]][25:30,] = NA; Cs5.8[[i]][32:37,] = NA; Cs5.8[[i]][39:44,] = NA; Cs5.8[[i]][46:51,] = NA; 
  Cs5.8[[i]][53:58,] = NA; Cs5.8[[i]][60:65,] = NA; Cs5.8[[i]][67:72,] = NA; Cs5.8[[i]][74:76,] = NA; 
}
Cs5.8[[4]]

Cs5.9 = Cs5.3
for (i in 1:100){
  Cs5.9[[i]][1:16,]=NA
  Cs5.9[[i]][18:23,] = NA; Cs5.9[[i]][25:30,] = NA; Cs5.9[[i]][32:37,] = NA; Cs5.9[[i]][39:44,] = NA; Cs5.9[[i]][46:51,] = NA; 
  Cs5.9[[i]][53:58,] = NA; Cs5.9[[i]][60:65,] = NA; Cs5.9[[i]][67:72,] = NA; Cs5.9[[i]][74:76,] = NA; 
}
Cs5.9[[4]]




##########  ANALYZE ONCE WEEKLY FOR SCEN 4-5 ################################################
#Rerun these to be safe
library(runjags)
runjags.options(force.summary=TRUE)

pupamax = 16
n = 60+pupamax #adding pupamax to n is not necessary, but doesn't hurt, for scenario 1. But pupamax is added to all sims for simplicity.  
#also note that in sim code I don't combine pupamax with n, though maybe I should have
MT = rep(21,n)

#calculate accumulated day-degrees for each cohort
ADD = matrix(nrow=n+pupamax, ncol=n); colnames(ADD) = 1:n
for (j in 1:n){
  for (i in (j+1):(n+1)) {
    ADD[i,j] = sum(MT[j:(i-1)])
  }
}
#set upper right off diagonal cells to 0 - these have to be set to 0 for JAGS (this code was originally in the model)
for (j in 1:n) {
  for (i in 1:j){
    ADD[i,j] = 0    
  }
}
#set lower cells to 0
for (j in 1:n) {
  for (i in (n+2):(n+pupamax)){
    ADD[i,j] = 0    
  }
}

#calculate stage durations for each cohort (this code was originally in the JAGS model, which is why it is a string of ifelse() statements instead of something more elegant)
SDC = matrix(nrow = 7, ncol = n); colnames(SDC) = 1:n; rownames(SDC) = c("egg","first","second","third","fourth","fifth","pupa")
for(j in 1:n) {
  SDC[1,j] = 1
  SDC[2,j] = ifelse(ADD[j+2,j] > 45,3,
                    ifelse(ADD[j+3,j] > 45,4,
                           5))
  SDC[3,j] = ifelse(ADD[j+3,j] > 77.3,4,
                    ifelse(ADD[j+4,j] > 77.3,5,
                           ifelse(ADD[j+5,j] > 77.3,6,
                                  7)))
  SDC[4,j] = ifelse(ADD[j+4,j] > 105.1,5,
                    ifelse(ADD[j+5,j] > 105.1,6,
                           ifelse(ADD[j+6,j] > 105.1,7,
                                  ifelse(ADD[j+7,j] > 105.1,8,
                                         9))))
  SDC[5,j] = ifelse(ADD[j+5,j] > 129.6,6,
                    ifelse(ADD[j+6,j] > 129.6,7,
                           ifelse(ADD[j+7,j] > 129.6,8,
                                  ifelse(ADD[j+8,j] > 129.6,9,
                                         10))))
  SDC[6,j] = ifelse(ADD[j+6,j] > 165.3,7,
                    ifelse(ADD[j+7,j] > 165.3,8,
                           ifelse(ADD[j+8,j] > 165.3,9,
                                  ifelse(ADD[j+9,j] > 165.3,10,
                                         ifelse(ADD[j+10,j] > 165.3,11,
                                                ifelse(ADD[j+11,j] > 165.3,12,
                                                       13))))))
  SDC[7,j] = ifelse(ADD[j+9,j] > 231.9,10,
                    ifelse(ADD[j+10,j] > 231.9,11,
                           ifelse(ADD[j+11,j] > 231.9,12,
                                  ifelse(ADD[j+12,j] > 231.9,13,
                                         ifelse(ADD[j+13,j] > 231.9,14,
                                                ifelse(ADD[j+14,j] > 231.9,15,
                                                       ifelse(ADD[j+15,j] > 231.9,16,
                                                              17)))))))
}



#########  calculate mean stage durations  ###############
#vector of mean stage durations
MSD = vector(length = 6)
#cohorts that finished 5th instar
I5 = ADD[(n+1),] > 231.9
#stage durations of cohorts that reached pupation
SDC[7,I5]
#fifth instar stage durations for each cohort
SD5 = SDC[7,I5]-SDC[6,I5]
MSD[6] = mean(SD5)
#cohorts that finished 4th instar
I4 = ADD[(n+1),] > 165.3
MSD[5] = mean(SDC[6,I4]-SDC[5,I4])
#cohorts that finished 3rd instar
I3 = ADD[(n+1),] > 129.6
MSD[4] = mean(SDC[5,I3]-SDC[4,I3])
#cohorts that finished 2nd instar
I2 = ADD[(n+1),] > 105.1
MSD[3] = mean(SDC[4,I2]-SDC[3,I2])
#cohorts that finished 1st instar
I1 = ADD[(n+1),] > 77.3
MSD[2] = mean(SDC[3,I1]-SDC[2,I1])
#cohorts that hatched
I0 = ADD[(n+1),] > 45
MSD[1] = mean(SDC[2,I0]-SDC[1,I0])
MSD

#calculate stage survival probabilities and cumulative survival probability
S1 = S[1]^MSD[1] 
S2 = S[2]^MSD[2]
S3 = S[3]^MSD[3]
S4 = S[4]^MSD[4]
S5 = S[5]^MSD[5]
S6 = S[6]^MSD[6]
SC = S1*S2*S3*S4*S5*S6 #cumulative survival probability
SN = S1*S2*S3*S4*S5 #survival probability most similar to Nail et al. 2015 estimator
SC


#specify initial values for three chains
inits = list(list(lambda2 = rep(50, n), S = rep(0.8,6)),
             list(lambda2 = rep(5, n), S = rep(0.5,6)),
             list(lambda2 = rep(100, n), S = rep(0.2,6)))

#list for results
R4.7 = list()
R4.8 = list()
R4.9 = list()
R5.7 = list()
R5.8 = list()
R5.9 = list()


#Scenario 4.7 - once weekly
runtime[[13]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=SDSC7[[i]])
    #call run.jags
    R4.7[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(paste0("Simulation 4.7.",i," finished running."))
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L4.7,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R4.7[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R4.7[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R4.7[[1]]$summary$quantile[,5]-R4.7[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R4.7[[i]]$summary$quantile[,5]-R4.7[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
ResAll[[13]] = Res
write.csv(Res,"Scen4.7.csv")
save(R4.7,Res,file = "Scen4.7.Rdata")
rm(R4.7)

#Scenario 4.8 - once weekly
runtime[[14]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=SDSC8[[i]])
    #call run.jags
    R4.8[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(paste0("Simulation 4.8.",i," finished running."))
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L4.8,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R4.8[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R4.8[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R4.8[[1]]$summary$quantile[,5]-R4.8[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R4.8[[i]]$summary$quantile[,5]-R4.8[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
ResAll[[14]] = Res
write.csv(Res,"Scen4.8.csv")
save(R4.8,Res,file = "Scen4.8.Rdata")
rm(R4.8)

#Scenario 4.9 - once weekly
runtime[[15]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=SDSC9[[i]])
    #call run.jags
    R4.9[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(paste0("Simulation 4.9.",i," finished running."))
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L4.9,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R4.9[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R4.9[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R4.9[[1]]$summary$quantile[,5]-R4.9[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R4.9[[i]]$summary$quantile[,5]-R4.9[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
ResAll[[15]] = Res
write.csv(Res,"Scen4.9.csv")
save(R4.9,Res,file = "Scen4.9.Rdata")
rm(R4.9)

#Scenario 5.7 - once weekly
runtime[[16]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=Cs5.7[[i]])
    #call run.jags
    R5.7[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(paste0("Simulation 5.7.",i," finished running."))
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L5.7,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R5.7[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R5.7[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R5.7[[1]]$summary$quantile[,5]-R5.7[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R5.7[[i]]$summary$quantile[,5]-R5.7[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
ResAll[[16]] = Res
write.csv(Res,"Scen5.7.csv")
save(R5.7,Res,file = "Scen5.7.Rdata")
rm(R5.7)

#Scenario 5.8 - once weekly
runtime[[17]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=Cs5.8[[i]])
    #call run.jags
    R5.8[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(paste0("Simulation 5.8.",i," finished running."))
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L5.8,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R5.8[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R5.8[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R5.8[[1]]$summary$quantile[,5]-R5.8[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R5.8[[i]]$summary$quantile[,5]-R5.8[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
ResAll[[17]] = Res
write.csv(Res,"Scen5.8.csv")
save(R5.8,Res,file = "Scen5.8.Rdata")
rm(R5.8)


#Scenario 5.9 - once weekly
runtime[[18]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=Cs5.9[[i]])
    #call run.jags
    R5.9[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(paste0("Simulation 5.9.",i," finished running."))
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L5.9,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R5.9[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R5.9[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R5.9[[1]]$summary$quantile[,5]-R5.9[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R5.9[[i]]$summary$quantile[,5]-R5.9[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
ResAll[[18]] = Res
write.csv(Res,"Scen5.9.csv")
save(R5.9,Res,file = "Scen5.9.Rdata")
rm(R5.9)

names(ResAll) = c("4.1","4.2","4.3","5.1","5.2","5.3","4.4","4.5","4.6","5.4","5.5","5.6","4.7","4.8","4.9","5.7","5.8","5.9")

save(ResAll,file = "ResAll2.Rdata")
save(runtime,file = "runtime2.Rdata")



######### Two Consecutive Days Every Week Sampling ###################################

#Scen 4.10
SDSC10 = SDSC1
for (i in 1:100){
  SDSC10[[i]][1:16,]=NA
  SDSC10[[i]][19:23,] = NA; SDSC10[[i]][26:30,] = NA; SDSC10[[i]][33:37,] = NA; SDSC10[[i]][40:44,] = NA; SDSC10[[i]][47:51,] = NA; 
  SDSC10[[i]][54:58,] = NA; SDSC10[[i]][61:65,] = NA; SDSC10[[i]][68:72,] = NA; SDSC10[[i]][75:76,] = NA; 
}
SDSC10[[5]]

#Scen 4.11
SDSC11 = SDSC2
for (i in 1:100){
  SDSC11[[i]][1:16,]=NA
  SDSC11[[i]][19:23,] = NA; SDSC11[[i]][26:30,] = NA; SDSC11[[i]][33:37,] = NA; SDSC11[[i]][40:44,] = NA; SDSC11[[i]][47:51,] = NA; 
  SDSC11[[i]][54:58,] = NA; SDSC11[[i]][61:65,] = NA; SDSC11[[i]][68:72,] = NA; SDSC11[[i]][75:76,] = NA; 
}
SDSC11[[4]]

#Scen 4.12
SDSC12 = SDSC3
for (i in 1:100){
  SDSC12[[i]][1:16,]=NA
  SDSC12[[i]][19:23,] = NA; SDSC12[[i]][26:30,] = NA; SDSC12[[i]][33:37,] = NA; SDSC12[[i]][40:44,] = NA; SDSC12[[i]][47:51,] = NA; 
  SDSC12[[i]][54:58,] = NA; SDSC12[[i]][61:65,] = NA; SDSC12[[i]][68:72,] = NA; SDSC12[[i]][75:76,] = NA; 
}
SDSC12[[4]]

#Scen 5.10
Cs5.10 = Cs5.1
for (i in 1:100){
  Cs5.10[[i]][1:16,]=NA
  Cs5.10[[i]][19:23,] = NA; Cs5.10[[i]][26:30,] = NA; Cs5.10[[i]][33:37,] = NA; Cs5.10[[i]][40:44,] = NA; Cs5.10[[i]][47:51,] = NA; 
  Cs5.10[[i]][54:58,] = NA; Cs5.10[[i]][61:65,] = NA; Cs5.10[[i]][68:72,] = NA; Cs5.10[[i]][75:76,] = NA; 
}
Cs5.10[[4]]

#Scen 5.11
Cs5.11 = Cs5.2
for (i in 1:100){
  Cs5.11[[i]][1:16,]=NA
  Cs5.11[[i]][19:23,] = NA; Cs5.11[[i]][26:30,] = NA; Cs5.11[[i]][33:37,] = NA; Cs5.11[[i]][40:44,] = NA; Cs5.11[[i]][47:51,] = NA; 
  Cs5.11[[i]][54:58,] = NA; Cs5.11[[i]][61:65,] = NA; Cs5.11[[i]][68:72,] = NA; Cs5.11[[i]][75:76,] = NA; 
}
Cs5.11[[4]]

#Scen 5.12
Cs5.12 = Cs5.3
for (i in 1:100){
  Cs5.12[[i]][1:16,]=NA
  Cs5.12[[i]][19:23,] = NA; Cs5.12[[i]][26:30,] = NA; Cs5.12[[i]][33:37,] = NA; Cs5.12[[i]][40:44,] = NA; Cs5.12[[i]][47:51,] = NA; 
  Cs5.12[[i]][54:58,] = NA; Cs5.12[[i]][61:65,] = NA; Cs5.12[[i]][68:72,] = NA; Cs5.12[[i]][75:76,] = NA; 
}
Cs5.12[[4]]


############## ANALYSIS OF TWO CONSECUTIVE DAYS ###########################

#list for results
R4.10 = list()
R4.11 = list()
R4.12 = list()
R5.10 = list()
R5.11 = list()
R5.12 = list()


#Scenario 4.10 - two consecutive days per week
runtime[[19]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=SDSC10[[i]])
    #call run.jags
    R4.10[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(paste0("Simulation 4.10.",i," finished running."))
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L4.10,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R4.10[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R4.10[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R4.10[[1]]$summary$quantile[,5]-R4.10[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R4.10[[i]]$summary$quantile[,5]-R4.10[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
ResAll[[19]] = Res
write.csv(Res,"Scen4.10.csv")
save(R4.10,Res,file = "Scen4.10.Rdata")
rm(R4.10)


#Scenario 4.11 - two consecutive days per week
runtime[[20]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=SDSC11[[i]])
    #call run.jags
    R4.11[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(paste0("Simulation 4.11.",i," finished running."))
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L4.11,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R4.11[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R4.11[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R4.11[[1]]$summary$quantile[,5]-R4.11[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R4.11[[i]]$summary$quantile[,5]-R4.11[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
ResAll[[20]] = Res
write.csv(Res,"Scen4.11.csv")
save(R4.11,Res,file = "Scen4.11.Rdata")
rm(R4.11)


#Scenario 4.12 - two consecutive days per week
runtime[[21]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=SDSC12[[i]])
    #call run.jags
    R4.12[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(paste0("Simulation 4.12.",i," finished running."))
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L4.12,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R4.12[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R4.12[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R4.12[[1]]$summary$quantile[,5]-R4.12[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R4.12[[i]]$summary$quantile[,5]-R4.12[[i]]$summary$quantile[,1])
}
CI[89,]#just SC CIs
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
ResAll[[21]] = Res
write.csv(Res,"Scen4.12.csv")
save(R4.12,Res,file = "Scen4.12.Rdata")
rm(R4.12)

#troubleshooting the large CIS for R4.12
CI[89,]#just SC CIs
#get psrf
psrf = R4.12[[1]]$psrf$psrf[89,1]
for (i in 2:100){
  psrf[i] = R4.12[[i]]$psrf$psrf[89,1]
}
psrf
gud = psrf < 1.10
sum(gud) #43 has psrf less than 1.10

R4.12S = R4.12[gud]

#Scenario 5.10 - two consecutive days per week
runtime[[22]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=Cs5.10[[i]])
    #call run.jags
    R5.10[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(paste0("Simulation 5.10.",i," finished running."))
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L5.10,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R5.10[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R5.10[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R5.10[[1]]$summary$quantile[,5]-R5.10[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R5.10[[i]]$summary$quantile[,5]-R5.10[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
ResAll[[22]] = Res
write.csv(Res,"Scen5.10.csv")
save(R5.10,Res,file = "Scen5.10.Rdata")
rm(R5.10)


#Scenario 5.11 - two consecutive days per week
runtime[[23]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=Cs5.11[[i]])
    #call run.jags
    R5.11[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(paste0("Simulation 5.11.",i," finished running."))
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L5.11,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R5.11[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R5.11[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R5.11[[1]]$summary$quantile[,5]-R5.11[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R5.11[[i]]$summary$quantile[,5]-R5.11[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
ResAll[[23]] = Res
write.csv(Res,"Scen5.11.csv")
save(R5.11,Res,file = "Scen5.11.Rdata")
rm(R5.11)



#Scenario 5.12 - two consecutive days per week
runtime[[24]]=system.time(
  for (i in 1:100){
    #specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=Cs5.12[[i]])
    #call run.jags
    R5.12[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(paste0("Simulation 5.12.",i," finished running."))
  }
)

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L5.12,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R5.12[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R5.12[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R5.12[[1]]$summary$quantile[,5]-R5.12[[1]]$summary$quantile[,1])
for (i in 2:100){
  CI = cbind(CI,R5.12[[i]]$summary$quantile[,5]-R5.12[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
ResAll[[24]] = Res
write.csv(Res,"Scen5.12.csv")
save(R5.12,Res,file = "Scen5.12.Rdata")
rm(R5.12)


names(ResAll) = c("4.1","4.2","4.3","5.1","5.2","5.3","4.4","4.5","4.6","5.4","5.5","5.6","4.7","4.8","4.9","5.7","5.8","5.9","4.10","4.11","4.12","5.10","5.11","5.12")

save(ResAll,file = "ResAll3.Rdata")
save(runtime,file = "runtime3.Rdata")



#fix ResAll - true values didn't get put in for a bunch of the sims.  checked estimates against .csv's and they are correct.  
#problem started at 4.7, 5.6 was fine.  I think I failed to calculate S, S1, ..., SC.  
#Also need to add 4.1-4.3

#5.6
ResAll[[12]]
#4.7
ResAll[[13]]
ResAll[[13]][,1] = c(rep(0,16),L4.7,S,S1,S2,S3,S4,S5,S6,SC)
#4.8
ResAll[[14]]
ResAll[[14]][,1] = c(rep(0,16),L4.8,S,S1,S2,S3,S4,S5,S6,SC)
#4.9
ResAll[[15]][,1] = c(rep(0,16),L4.9,S,S1,S2,S3,S4,S5,S6,SC)
ResAll[[15]]
#5.7
ResAll[[16]][,1] = c(rep(0,16),L5.7,S,S1,S2,S3,S4,S5,S6,SC)
ResAll[[16]]
#5.8
ResAll[[17]][,1] = c(rep(0,16),L5.8,S,S1,S2,S3,S4,S5,S6,SC)
ResAll[[17]]
#5.9
ResAll[[18]][,1] = c(rep(0,16),L5.9,S,S1,S2,S3,S4,S5,S6,SC)
ResAll[[18]]
#4.10
ResAll[[19]][,1] = c(rep(0,16),L4.7,S,S1,S2,S3,S4,S5,S6,SC)
ResAll[[19]]
#4.11
ResAll[[20]][,1] = c(rep(0,16),L4.8,S,S1,S2,S3,S4,S5,S6,SC)
ResAll[[20]]
#4.12
ResAll[[21]][,1] = c(rep(0,16),L4.9,S,S1,S2,S3,S4,S5,S6,SC)
ResAll[[21]]
#5.10
ResAll[[22]][,1] = c(rep(0,16),L5.7,S,S1,S2,S3,S4,S5,S6,SC)
ResAll[[22]]
#5.11
ResAll[[23]][,1] = c(rep(0,16),L5.8,S,S1,S2,S3,S4,S5,S6,SC)
ResAll[[23]]
#5.12
ResAll[[24]][,1] = c(rep(0,16),L5.9,S,S1,S2,S3,S4,S5,S6,SC)
ResAll[[24]]

#calc bias and percent bias
for (i in 13:24){
  ResAll[[i]][,5] = ResAll[[i]][,2]-ResAll[[i]][,1]
  ResAll[[i]][,6] = (ResAll[[i]][,5]/ResAll[[i]][,1])*100
}

save(ResAll,file = "ResAllFinal.Rdata")

#get 4.1-4.3
load("ResAllFinal.Rdata")
load("Scen4.1.Rdata")
ResAll[[1]] = Res
save(ResAll,file = "ResAllFinal.Rdata")
ResAll[[2]] = Res
save(ResAll,file = "ResAllFinal.Rdata")
ResAll[[3]] = Res
save(ResAll,file = "ResAllFinal.Rdata")

#I guess I need to put these into Excel to easily format the data for plotting :(
#write each list element to a separate .csv named after the list name
mapply(write.csv, ResAll, paste0(names(ResAll), ".csv"))


####################### FIGURES  ###########################################

#Scenario 4

#plot CI Width of SC
library(ggplot2)
plotSCCI = read.csv("SCCIWidth2.csv")
#levels for legend order
levels(plotSCCI$Sampling)
plotSCCI$Sampling = factor(plotSCCI$Sampling, levels = c("Daily","Two Days","2x Weekly","1x Weekly"))

p = ggplot(plotSCCI, aes(x=Lmax, y=CI.Width, linetype=Sampling)) + geom_line(size = 1)
p + ylab("Credible Interval Width") + xlab("Maximum Eggs Laid per Day") +
  theme(legend.position = c(0.85,0.8)) +
  #scale_linetype_manual(values = c("dotdash", "longdash", "solid", "dotted")) +
  #scale_linetype_manual(values = c("solid", "dotted", "longdash", "dotdash")) +
  scale_linetype_manual(values = c("solid", "dotdash","longdash","dotted")) +
  guides(linetype = guide_legend(override.aes = list(size = 0.5))) 


#plot bias in SC
plotSCBias = read.csv("SCBias2.csv")
plotSCBias$Sampling = factor(plotSCBias$Sampling, levels = c("Daily","Two Days","2x Weekly","1x Weekly"))

p = ggplot(plotSCBias, aes(x=Lmax, y=Bias, linetype=Sampling)) + geom_line(size = 1)
p + ylab("Bias") + xlab("Maximum Eggs Laid per Day") +
  theme(legend.position = c(0.85,0.7)) +
  #scale_linetype_manual(values = c("dotdash", "longdash", "solid", "dotted")) +
  #scale_linetype_manual(values = c("longdash", "dotdash", "solid", "dotted")) +
  scale_linetype_manual(values = c("solid", "dotdash","longdash","dotted")) +
  guides(linetype = guide_legend(override.aes = list(size = 0.5)))



#Scenario 5

#plot CI Width of SC
plotIndVarSD = read.csv("StudyDesignIndVarSimsPlot2.csv")
levels(plotIndVarSD$Sampling)
plotIndVarSD$Sampling = factor(plotIndVarSD$Sampling, levels = c("Daily","Two Days  ","2x Weekly","1x Weekly"))

p = ggplot(plotIndVarSD, aes(x=Lmax, y=MeanCIWidth, linetype=Sampling)) + geom_line(size = 1)
p + ylab("Credible Interval Width") + xlab("Maximum Eggs Laid per Day") +
  theme(legend.position = c(0.85,0.8)) +
  #scale_linetype_manual(values = c("dotdash", "longdash", "solid", "dotted")) +
  #scale_linetype_manual(values = c("longdash", "dotted", "solid", "dotdash")) +
  scale_linetype_manual(values = c("solid", "dotdash","longdash","dotted")) +
  guides(linetype = guide_legend(override.aes = list(size = 0.5)))

#plot bias in SC
p = ggplot(plotIndVarSD, aes(x=Lmax, y=MeanSCBias, linetype=Sampling)) + geom_line(size = 1)
p + ylab("Bias") + xlab("Maximum Eggs Laid per Day") +
  theme(legend.position = c(0.85,0.30)) +
  #scale_linetype_manual(values = c("dotdash", "longdash", "solid", "dotted")) +
  #scale_linetype_manual(values = c("dotted", "longdash", "dotdash", "solid")) +
  scale_linetype_manual(values = c("solid", "dotdash","longdash","dotted")) +
  guides(linetype = guide_legend(override.aes = list(size = 0.5)))

range(plotIndVarSD$MeanSCBias)


#code to subset 4.12 and 4.9 to get runs with psrf < 1.10
#code to troubleshoot the problem is under 4.12 above

#subset of 4.12 of only psrf < 1.10
R4.12S

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L4.12,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R4.12S[[1]]$summary$quantile[,3]
for (i in 2:43){
  P = cbind(P, R4.12S[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R4.12S[[1]]$summary$quantile[,5]-R4.12S[[1]]$summary$quantile[,1])
for (i in 2:43){
  CI = cbind(CI,R4.12S[[i]]$summary$quantile[,5]-R4.12S[[i]]$summary$quantile[,1])
}
CI[89,]#just SC CIs
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
write.csv(Res,"Scen4.12S.csv")
save(R4.12S,Res,file = "Scen4.12S.Rdata")


#subset 4.9
#get psrf
psrf = R4.9[[1]]$psrf$psrf[89,1]
for (i in 2:100){
  psrf[i] = R4.9[[i]]$psrf$psrf[89,1]
}
psrf
gud = psrf < 1.10
sum(gud) #28 has psrf less than 1.10
R4.9S = R4.9[gud]
R4.9[[38]]#bad
R4.9[[1]]#good

Res = data.frame(matrix(nrow = 16+60+13, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),L4.12,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R4.9S[[1]]$summary$quantile[,3]
for (i in 2:28){
  P = cbind(P, R4.9S[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R4.9S[[1]]$summary$quantile[,5]-R4.9S[[1]]$summary$quantile[,1])
for (i in 2:28){
  CI = cbind(CI,R4.9S[[i]]$summary$quantile[,5]-R4.9S[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
write.csv(Res,"Scen4.9S.csv")
save(R4.9S,Res,file = "Scen4.9S.Rdata")

