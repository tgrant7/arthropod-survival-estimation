##############  R Code for Grant et al.  #####################################################################
##############  Simulate and Analyze Monarch Butterfly Field Counts  ###############################

# Table of Contents
#1.  Calculate Expected Values of Counts - Line 14
#2.  Simulate counts with individual variation in stage duration - Line 275
#3.  Analysis of simulated counts - Line 655
#4.  Analysis of field data - Line 786

##################################################################################################



#################### 1. Simple Simulation of Counts for Given Input Variables #######################

##Input Variables##

#number of days in the simulation
n=50

#maximum possible time from oviposition to pupation in monarch butterflies
#used to determine the earliest day on which eggs could have been laid
#also used for a buffer at the end of ADD, which has no biological meaning, it's just for convenience
pupamax = 16

#B is the number of eggs laid per day.  Should be of length n
#example scenario from the paper
B = c(50,100,200,500,750,500,200,100,50,25);B[11:n] = 0  
#or simulate eggs laid as a normal distr with max eggs laid as n
#  SF = n/dnorm(n/2, mean = n/2, sd = n/4) #scale factor
#  B = round(SF*dnorm(1:n, mean = n/2, sd = n/4))
B

#daily survival probabilities
#Note that these are denoted as the Greek letter Pi in Grant et al.
S = c(0.53,0.68,0.90,0.94,0.97,0.96)

#mean daily temps - denoted as t in Grant et al.
MT = rep(21,n) #same temp every day for simulations


##Simulate Counts##

#ADD = Accumulated Day Degrees for each cohort over time.  D_jk in Grant et al.  Rows are days, cols are cohorts
#the buffer rows at the end of number pupamax is so that the indexing doesn't throw errors
#it is an inelegant solution but it works.  pupamax has no biological meaning as buffer at the end, it is just conveniently available
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


#PEDD = Proportion surviving each day for each cohort.  Rows are days, cols are cohorts.  This is p_ik*pi_j-k-bik from Grant et al.
PEDD = matrix(nrow=n+pupamax, ncol=n); colnames(PEDD) = 1:n
for(j in 1:n){
  PEDD[j,j] = 1
  for (i in 2:SDC[2,j]){
    PEDD[i+j-1,j] = (S[1])^(i-1)
  }
  for (i in (SDC[2,j]+1):SDC[3,j]){
    PEDD[i+j-1,j] = (S[1])^(SDC[2,j]-1)*(S[2])^(i-SDC[2,j])
  }
  for (i in (SDC[3,j]+1):SDC[4,j]){
    PEDD[i+j-1,j] = (S[1])^(SDC[2,j]-1)*(S[2])^(SDC[3,j]-SDC[2,j])*S[3]^(i-SDC[3,j])
  }
  for (i in (SDC[4,j]+1):SDC[5,j]){
    PEDD[i+j-1,j] = (S[1])^(SDC[2,j]-1)*(S[2])^(SDC[3,j]-SDC[2,j])*S[3]^(SDC[4,j]-SDC[3,j])*S[4]^(i-SDC[4,j])
  }
  for (i in (SDC[5,j]+1):SDC[6,j]){
    PEDD[i+j-1,j] = (S[1])^(SDC[2,j]-1)*(S[2])^(SDC[3,j]-SDC[2,j])*S[3]^(SDC[4,j]-SDC[3,j])*S[4]^(SDC[5,j]-SDC[4,j])*S[5]^(i-SDC[5,j])
  }
  for (i in (SDC[6,j]+1):SDC[7,j]){
    PEDD[i+j-1,j] = (S[1])^(SDC[2,j]-1)*(S[2])^(SDC[3,j]-SDC[2,j])*S[3]^(SDC[4,j]-SDC[3,j])*S[4]^(SDC[5,j]-SDC[4,j])*S[5]^(SDC[6,j]-SDC[5,j])*S[6]^(i-SDC[6,j])
  }
}
#set NA's to 0
for (j in 2:n) {
  for (i in 1:(j-1)) {
    PEDD[i,j] = 0
  }
}
for (j in 1:(n-1)) {
  for (i in (SDC[7,j]+j):(n+pupamax)) {
    PEDD[i,j] = 0
  }
}


#Number surviving each day of each cohort.  rows are days, cols are cohorts.  These are the cohort contributions that are summed to get l_ij of Grant et al.  
NDEDD = matrix(nrow = n+pupamax, ncol = n); colnames(NDEDD) = 1:n
for (j in 1:n) {
  for (i in j:(j+SDC[7,j]-1)){
    NDEDD[i,j] = B[j]*PEDD[i,j]
  }
}
for (j in 2:n) {
  for (i in 1:(j-1)) {
    NDEDD[i,j] = 0
  }
}
for (j in 1:(n-1)) {
  for (i in (j+SDC[7,j]):(n+pupamax)){
    NDEDD[i,j] = 0
  }
}

#ST - index indicating what stage each cohort is in each day.  rows are days, cols are cohorts.  99 means it has pupated
#this is essentially the indicator function I(s_kj = i) of Grant et al. 
ST = matrix(nrow=n, ncol=n); colnames(ST) = 1:n
for (j in 1:(n-1)) {
  for (i in (j+1):n) {
    ST[i,j] = ifelse(ADD[i,j] >= 1 && ADD[i,j] < 45, 1, 
                     ifelse(ADD[i,j] >= 45 && ADD[i,j] < 77.3, 2,
                            ifelse(ADD[i,j] >= 77.3 && ADD[i,j] < 105.1, 3,
                                   ifelse(ADD[i,j] >= 105.1 && ADD[i,j] < 129.6, 4,
                                          ifelse(ADD[i,j] >= 129.6 && ADD[i,j] < 165.3, 5,
                                                 ifelse(ADD[i,j] >= 165.3 && ADD[i,j] < 231.9, 6, 99))))))
  }
}
for (j in 1:n){
  for (i in j:j){
    ST[i,j] = 1
  }
}
for (j in 2:n){
  for (i in 1:(j-1)){
    ST[i,j] = 0
  }
}

#A matrix for each stage - copy numbers from NDEDD using ST to tell which numbers go into which matrix
#i.e., if a ST[4,5] = 1, then NDEDD[4,5] goes into matrix M1
#this makes it easier/possible to sum the population across cohorts in each stage on any particular day
M1 = matrix(nrow=n, ncol=n); colnames(M1) = 1:n
M2 = matrix(nrow=n, ncol=n); colnames(M2) = 1:n
M3 = matrix(nrow=n, ncol=n); colnames(M3) = 1:n
M4 = matrix(nrow=n, ncol=n); colnames(M4) = 1:n
M5 = matrix(nrow=n, ncol=n); colnames(M5) = 1:n
M6 = matrix(nrow=n, ncol=n); colnames(M6) = 1:n

for (j in 1:n){
  for (i in 1:n) {
    M1[i,j] = ifelse(ST[i,j]==1,NDEDD[i,j],0)
    M2[i,j] = ifelse(ST[i,j]==2,NDEDD[i,j],0)
    M3[i,j] = ifelse(ST[i,j]==3,NDEDD[i,j],0)
    M4[i,j] = ifelse(ST[i,j]==4,NDEDD[i,j],0)
    M5[i,j] = ifelse(ST[i,j]==5,NDEDD[i,j],0)
    M6[i,j] = ifelse(ST[i,j]==6,NDEDD[i,j],0)
  }
}

#NDD = sum populations from each cohort into a matrix of population of each stage each day.  rows are days, cols are stages.  
#This is l_ij of Grant et al. 
NDD = matrix(nrow=n,ncol=6); colnames(NDD) = c("E","I","II","III","IV","V")

for (i in 1:n) {
  NDD[i,1] =  sum(M1[i,])
}
for (i in 1:n) {
  NDD[i,2] =  sum(M2[i,])
}
for (i in 1:n) {
  NDD[i,3] =  sum(M3[i,])
}
for (i in 1:n) {
  NDD[i,4] =  sum(M4[i,])
}
for (i in 1:n) {
  NDD[i,5] =  sum(M5[i,])
}
for (i in 1:n) {
  NDD[i,6] =  sum(M6[i,])
}
NDD

#add pupamax buffer to beginning of count matrix
Buffer = matrix(0, nrow = pupamax, ncol = 6)
NDDB = rbind(Buffer,NDD)
NDDB

#scenario 1
#Simulate Y_ij as from a Poisson distribution with mean l_ij
#100 simulated count matrices
C = list()
for (i in 1:100){
  C[[i]] = apply(NDDB, MARGIN = c(1,2), FUN = function(x) rpois(1,x))
}



#detection probability scenarios - Scenario 3
  #relative det prob for each stage
  #1st scenario - 100% for eggs, but lower for 5th
  Det1 = c(1,0.95,0.9,0.85,0.8,0.75)
  #2nd scenario - 100% for 5th, but lower for eggs
  Det2 = c(0.75,0.8,0.85,0.9,0.95,1)

C = list()
CDet1 = list()
CDet2 = list()
for (i in 1:100){
  C[[i]] = apply(NDDB, MARGIN = c(1,2), FUN = function(x) rpois(1,x))

  #det prob calcs
  CDet1[[i]] = C[[i]]
  for(j in 1:6) {
    CDet1[[i]][,j] = C[[i]][,j]*Det1[j]
  }
  CDet1[[i]] = round(CDet1[[i]])

  CDet2[[i]] = C[[i]]
  for(j in 1:6) {
    CDet2[[i]][,j] = C[[i]][,j]*Det2[j]
  }
  CDet2[[i]] = round(CDet2[[i]])
}





######### 2. Simulated data with individual stochasticity in survival and individual variation in stage duration #######################
#this section uses much of the same code as the last section, with the addition of code for random stage durations and stochasticity in survival
library(msm) #for rtnorm function, a function for taking random numbers from a truncated normal distribution


############ Input Variables ##########################################

n=50
pupamax = 16
S = c(0.53,0.68,0.90,0.94,0.97,0.96)
B = c(50,100,200,500,750,500,200,100,50,25);B[11:n] = 0  #example scenario from Grant et al.

#number of cohorts, ie. number of days that eggs were laid
#as the code is written, simulated cohorts must be sequential, and days with no eggs must follow all the days with eggs
#i.e., no days with no eggs laid between days with eggs laid
#and the first day must have eggs laid
CH = sum(B != 0)

#MT needs to be of length n+pupamax instead of just n
MT = rep(21,n+pupamax)

#list of simulated datasets
Cs = list()
#large loop to create many simulated datasets - closing bracket at ~ L670
for (m in 1:100){


#random stage durations based on Zalucki 1982 Table 2 field data. Rows are individuals, cols are stages
#units are degree days that it takes each individual to transition to the next stage
SDI = matrix(nrow = sum(B), ncol = 6)
SDI[,1] = rtnorm(sum(B), mean = 44.5, sd = 4.32, lower = 44.5-2*4.32, upper = 44.5+2*4.32)
SDI[,2] = rtnorm(sum(B), mean = 30.9, sd = 3.97, lower = 30.9-2*3.97, upper = 30.9+2*3.97)
SDI[,3] = rtnorm(sum(B), mean = 27.2, sd = 3.71, lower = 27.2-2*3.71, upper = 27.2+2*3.71)
SDI[,4] = rtnorm(sum(B), mean = 25.6, sd = 2.88, lower = 25.6-2*2.88, upper = 25.6+2*2.88)
SDI[,5] = rtnorm(sum(B), mean = 33.8, sd = 4.09, lower = 33.8-2*4.09, upper = 33.8+2*4.09)
SDI[,6] = rtnorm(sum(B), mean = 58.5, sd = 3.89, lower = 58.5-2*3.89, upper = 58.5+2*3.89)

########################################################################

#define matrices to be filled from code below
#MI1-MI6 keep track of the survival history for each egg.  Each day a monarch is alive, it's coded as 1 in the matrix
#after it dies, it's coded as 0 in the matrix.  Each matrix is for a different stage - MI1 for eggs, MI2 for 1st instar, etc.
#state I (alive or dead) and day index i (the same day index is used for each matrix) is carried over between matrices with IM and iM
#i is last day an individual was in a stage; it's not the same for each individual because of varying stage duration

#rows are individual monarchs, columns are days
MI1 = matrix(nrow=sum(B), ncol=n+pupamax); colnames(MI1) = 1:(n+pupamax)
MI2 = matrix(nrow=sum(B), ncol=n+pupamax); colnames(MI2) = 1:(n+pupamax)
MI3 = matrix(nrow=sum(B), ncol=n+pupamax); colnames(MI3) = 1:(n+pupamax)
MI4 = matrix(nrow=sum(B), ncol=n+pupamax); colnames(MI4) = 1:(n+pupamax)
MI5 = matrix(nrow=sum(B), ncol=n+pupamax); colnames(MI5) = 1:(n+pupamax)
MI6 = matrix(nrow=sum(B), ncol=n+pupamax); colnames(MI6) = 1:(n+pupamax)

#matrix of I for each individual
IM = matrix(nrow=sum(B), ncol=6)
#matrix of i (day/col of MI) for each individual, to carry over to each stage from previous stage
iM = matrix(nrow=sum(B), ncol=6)

#final counts of each stage each day
C = matrix(nrow = n+pupamax, ncol = 6)


############### Calculations on Input Variables ########################

#ADD - added extra cols of number pupamax this sim and set value to 1000
#this forces while loops below stop when ADD > 231.9
ADD = matrix(nrow=n+pupamax, ncol=n+pupamax); colnames(ADD) = 1:(n+pupamax)
for (j in 1:n){
  for (i in (j+1):(n+1)) {
    ADD[i,j] = sum(MT[j:(i-1)])
  }
}
#set upper right off diagonal cells to 0
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

#Next simulate history of each individual monarch.  The code is set up to simulate by stage.
#each stage is a separate code block

############ Eggs  ##################
#initialize Inds index to 1 - this is the index for individual monarchs
Inds = 1
#k loops through cohorts
for (k in 1:CH){
  #set initial value to 1
  MI1[Inds:(Inds+B[k]-1),k] = 1
  #j loops through individuals of a cohort
  for (j in 1:B[k]){
    I=1
    i=1+k #i is the col/day index
    while (ADD[i,k] < SDI[Inds,1]) {
      R = runif(1) #choose a random number
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
  for (j in 1:B[k]){
    #take ending I from previous stage
    I=IM[Inds,1]
    i=iM[Inds,1] #i is the col/day index
    #index for no of times it has been through while loop
    l = 1
    while (ADD[i,k] < (SDI[Inds,1]+SDI[Inds,2])) {
      R = runif(1)
      #survival before first detection as first instar - in other words, the survival rate previous to the first time a 1st instar
      #larvae is found is assumed to be egg survival
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
  for (j in 1:B[k]){
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
    
    if (ADD[i,k] >= (SDI[Inds,1]+SDI[Inds,2]+SDI[Inds,3])) { #this if loop takes care of IM and iM if this instar is skipped because short stage dur
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
  for (j in 1:B[k]){
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
    
    if (ADD[i,k] >= (SDI[Inds,1]+SDI[Inds,2]+SDI[Inds,3]+SDI[Inds,4])) { #this if loop takes care of IM and iM if this instar skipped because its warm
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
  for (j in 1:B[k]){
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
    
    if (ADD[i,k] >= (SDI[Inds,1]+SDI[Inds,2]+SDI[Inds,3]+SDI[Inds,4]+SDI[Inds,5])) { #this if loop takes care of IM and iM if this instar skipped because its warm
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
  for (j in 1:B[k]){
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


#sum columns to get individual in each stage per day
#NS1 = summed N eggs = number of eggs in the population each day
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

#remove buffer rows from end
C = C[1:n,]
#To allow the model to estimate eggs laid before the first day of the study, 16 days are added before the first day of the study
#add buffer rows to beginning
Buffer = matrix(0, nrow = pupamax, ncol = 6)
C = rbind(Buffer,C)
C

#apply random Poisson variation and add to list of datasets
Cs[[m]] = apply(C, MARGIN = c(1,2), FUN = function(x) rpois(1,x))

}





#######################  3. Analysis of Simulated Counts  ############################
#JAGS
library(runjags)
#option to always calculate summary statistics
runjags.options(force.summary=TRUE)

#The input data is n, pupamax, MT, and C
#n and MT should be the same as for the simulated data 
#some of the below code may be repetition of above code, but is included here so this section can be run with only the below input variables

pupamax = 16
n = 50+pupamax #adding pupamax to n is not necessary, but doesn't hurt, for scenario 1. But pupamax is added to all sims for simplicity.  
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
R = list()
R2 = list()

#Scenario 2
system.time(
for (i in 1:100){
#specify data.  be sure to create ADD, SDC, and MSD.  C is the count matrix.  
data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=Cs[[i]])
#call run.jags
R[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
i
}
)
plot(R[[1]], vars = c("S","S1","S2","S3","S4","S5","S6","SC"))


#scenario 3
system.time(
  for (i in 1:100){
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=CDet1[[i]])
    R[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(i)
  }
)
system.time(
  for (i in 1:100){
    data = list(n=n, pupamax=pupamax, SDC = SDC, ADD = ADD, MSD = MSD, C=CDet2[[i]])
    R2[[i]] = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 3, data = data, inits = inits)
    print(i)
  }
)





#summary statistics
Res = data.frame(matrix(nrow = 79, ncol = 6)); colnames(Res) = c("True","MeanEst","SDofMeanEst","MeanCI","Bias","PercentBias")
Res[,1] = c(rep(0,16),B,S,S1,S2,S3,S4,S5,S6,SC)
#matrix of parameter ests (median of posterior)
P = R[[1]]$summary$quantile[,3]
for (i in 2:100){
  P = cbind(P, R[[i]]$summary$quantile[,3])
}
Res[,2] = apply(P, MARGIN = 1, mean)
Res[,3] = apply(P, 1, sd)
#ci width
CI = c(R[[1]]$summary$quantile[,5]-R[[1]]$summary$quantile[,1])
for (i in 2:100){
CI = cbind(CI,R[[i]]$summary$quantile[,5]-R[[i]]$summary$quantile[,1])
}
Res[,4] = apply(CI, 1, mean)
Res[,5] = Res[,2]-Res[,1]
Res[,6] = (Res[,5]/Res[,1])*100
write.csv(Res,"Scen3_2.csv")




#Simulation Scenario 1 results
R1 = R
Res1 = Res
save(R1,Res1,file = "Scen1.Rdata")

#Simulation Scenario 2 results
R2 = R
Res2 = Res
save(R2,Res2,file = "Scen2.Rdata")

#Scen 3, Case 1
R3.1 = R
Res3.1 = Res
save(R3.1,Res3.1,file = "Scen3.1.Rdata")

#Scen 3, Case 2
R3.2 = R2
Res3.2 = Res
save(R3.2,Res3.2,file = "Scen3.2.Rdata")

load(file = "C://Users//tgrant//Documents//Monarch Butterflies//Projects//Immature Survival Estimation//Sim Results Files 2//Scen3.1.Rdata")







####################### 4. Analysis of Field Data  #################################
#The major difference with real data is that it typically requires 16 leading days to be added to n
#this is because the model needs to be able to estimate eggs laid from any observations on the first day of field counts
#for example, if any 5th instars are observed on the first day, they were eggs up to 16 day previous to that, so the model must be defined for 16
#days prior to the first day.  The 16 leading days are not used to estimate MSD.  

#First dataset is Bitzer et al. 2016 data from 140th Street, collected in 2015
#data for the 3 other example datasets in Grant et al. follow at the end of this section

#maximum time to pupation,
pupamax = 16
#number of days in the study (48) + pupamax (16) = 64
n = 48 + pupamax
#mean temperature each day of the study period of length n days - from Ames WSW weather station
MT = c(21.9,21.9,20.8,17.5,18.1,19.7,22.2,24.7,17.2,16.4,18.9,20.3,23.3,26.9,28.3,26.4,
       24.4,23.65,26.7,26.7,24.15,23.6,20.8,20.55,21.95,23.6,26.1,24.4,24.75,22.5,21.95,21.95,22.75,23.6,25.25,23.6,20.85,21.7,21.4,23.85,24.45,24.15,
       22.5,21.7,22.2,21.7,24.4,23.6,23.6,24.15,19.75,15.8,17.2,20.25,20.85,17.5,15.25,15.85,16.65,18.05,18.9,18.35,18.85,22.8)
length(MT)

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


#########  calculate stage durations for each cohort and mean stage durations  ###############
#vector of mean stage durations
MSD = vector(length = 6)
#cohorts that finished 5th instar
I5 = ADD[(n+1),] > 231.9
#remove first 16 days before study started
I5[1:16] = FALSE
#stage durations of cohorts that reached pupation
SDC[7,I5]
#fifth instar stage durations for each cohort
SD5 = SDC[7,I5]-SDC[6,I5]
MSD[6] = mean(SD5); range(SD5)
#cohorts that finished 4th instar
I4 = ADD[(n+1),] > 165.3
I4[1:16] = FALSE
MSD[5] = mean(SDC[6,I4]-SDC[5,I4])
#cohorts that finished 3rd instar
I3 = ADD[(n+1),] > 129.6; I3[1:16] = FALSE
MSD[4] = mean(SDC[5,I3]-SDC[4,I3])
#cohorts that finished 2nd instar
I2 = ADD[(n+1),] > 105.1; I2[1:16] = FALSE
MSD[3] = mean(SDC[4,I2]-SDC[3,I2])
#cohorts that finished 1st instar
I1 = ADD[(n+1),] > 77.3; I1[1:16] = FALSE
MSD[2] = mean(SDC[3,I1]-SDC[2,I1])
#cohorts that hatched
I0 = ADD[(n+1),] > 45; I0[1:16] = FALSE
MSD[1] = mean(SDC[2,I0]-SDC[1,I0])
MSD

#Field Counts - counts were made approximately once per week for several weeks
#counts were taken on day 1 and the last count was taken on day 48
C = matrix(nrow = 48, ncol = 6)
C[1,] = c(3,0,0,1,0,0)
C[10,] = c(96,12,13,4,1,4)
C[20,] = c(138,8,8,3,0,2)
C[27,] = c(137,9,10,8,3,3)
C[36,] = c(65,3,5,1,2,2)
C[41,] = c(47,0,1,0,1,1)
C[48,] = c(27,1,0,0,1,0)

#To allow the model to estimate eggs laid before the first day of the study, 16 days are added before the first day of the study
Buffer = matrix(nrow = pupamax, ncol = 6)
CB = rbind(Buffer,C)

#run the model with runjags package
#data inputs for JAGS
data = list(n=n, pupamax=pupamax, MSD = MSD, SDC = SDC, ADD = ADD, C=CB)
#initial vlaues for lambda2 and S
inits = list(list(lambda2 = rep(200, n), S = rep(0.7,6)),
             list(lambda2 = rep(50, n), S = rep(0.8,6)),
             list(lambda2 = rep(10, n), S = rep(0.9,6)),
             list(lambda2 = rep(5, n), S = rep(0.2,6)),
             list(lambda2 = rep(75, n), S = rep(0.3,6)),
             list(lambda2 = rep(150, n), S = rep(0.4,6)),
             list(lambda2 = rep(100, n), S = rep(0.5,6)))
#runjags, monitoring daily survival, stage survival, cumulative survival (SC), and eggs laid B.  
out.T.140 = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
plot(out.T.140, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
#add some more samples
out.T.140.e = extend.jags(out.T.140, sample = 90000)
plot(out.T.140.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
write.csv(out.T.140.e$summaries, "out.T.csv")






#### 260th Stree dataset from Bitzer et al. 2016 data, collected in 2015 #############

pupamax = 16
n = 46 + pupamax
MT = c(21.9,21.9,20.8,17.5,18.1,19.7,22.2,24.7,17.2,16.4,18.9,20.3,23.3,26.9,28.3,26.4,
       24.4,23.65,26.7,26.7,24.15,23.6,20.8,20.55,21.95,23.6,26.1,24.4,24.75,22.5,21.95,21.95,22.75,23.6,25.25,23.6,20.85,21.7,21.4,23.85,24.45,24.15,
       22.5,21.7,22.2,21.7,24.4,23.6,23.6,24.15,19.75,15.8,17.2,20.25,20.85,17.5,15.25,15.85,16.65,18.05,18.9,18.35)
length(MT)

C = matrix(nrow = 46, ncol = 6)
C[1,] = c(37,2,4,0,0,0)
C[9,] = c(90,15,9,7,1,1)
C[17,] = c(135,19,5,5,3,7)
C[24,] = c(156,22,7,3,1,3)
C[31,] = c(121,5,2,3,1,3)
C[38,] = c(77,1,7,2,1,1)
C[46,] = c(11,1,0,2,3,2)


#### Trillium/Vanmeer dataset from Pitman et al. (2018)  ################################

n = 36+16
#maximum time to pupation, this shouldn't need changing, it makes indexing work
pupamax = 16
#mean temperature each day of the study period of length n days
MT = c(18.6,15,13,17.4,17.9,16.2,15.5,16.8,19.6,21.1,20.9,16.3,18.2,18.4,19.8,20.8,
       21.6,21.6,17.3,16.4,21,25.1,23.2,20.9,20.7,17.8,18.9,20.4,22.8,22.4,23.6,23.7,24.2,23.5,21.5,20,20.7,20.2,19.6,16.2,15.5,18.8,
       18.4,20,20,19.9,17,17.3,22.5,23.1,22.2,24.2)
length(MT)

C = matrix(nrow = 36, ncol = 6)
C[1,] = c(9,5,2,0,2,1)
C[8,] = c(59,8,5,1,0,2)
C[15,] = c(140,54,7,1,2,4)
C[22,] = c(133,60,27,13,6,3)
C[29,] = c(124,39,17,5,2,7)
C[36,] = c(51,22,16,2,5,0)


####  Mccracken/Quarter/Storm  dataset from Pitman et al. (2018) #######################

pupamax = 16
n = 36+pupamax
#mean temperature each day of the study period of length n days
MT = c(22.6,24.4,22,17,16.1,16.7,17.6,17.9,18.4,24.4,24.1,25.4,24.7,20.2,21.9,21.5,24.2,26.2,23.9,20.5,18.1,18,23.1,19.1,18.5,21.2,26.6,
       24.5,21.6,24.9,21.7,23,23.2,22.8,22.4,22,22,21.6,22.2,22.8,25.4,21.7,20.5,19.8,21.4,26,26.8,22,26.7,22.4,21.7,21.7)
length(MT)

C = matrix(nrow = 36, ncol = 6)
C[1,] = c(2,0,0,0,0,0)
C[8,] = c(1,0,0,1,0,0)
C[15,] = c(5,0,0,0,0,0)
C[22,] = c(29,16,3,0,0,1)
C[29,] = c(21,13,2,2,0,0)
C[36,] = c(25,2,2,3,5,3)







