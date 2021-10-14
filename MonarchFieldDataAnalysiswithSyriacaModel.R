#Analysis of Monarch Field Data for Grant et al.

#code below just copied directly from SimulateandAnalyzeMonarchCounts.R and then adapted
#Degree day numbers have been changed to new DD syriaca model. Should add them as input variable vector sometime. 

library(runjags)
#option to always calculate summary statistics
runjags.options(force.summary=TRUE)
runjags.options(jagspath = 'C:/Users/tgrant/AppData/Local/JAGS/JAGS-4.3.0/x64/bin')

#Example dataset is Bitzer et al. 2016 data from 140th Street, collected in 2015
#IA ROW1 in ms supp info

#maximum time or greater to pupation, this number serves as a buffer in some matrices/dataframes
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
  SDC[2,j] = ifelse(ADD[j+2,j] > 90.2,3,
                    ifelse(ADD[j+3,j] > 90.2,4,
                           5))
  SDC[3,j] = ifelse(ADD[j+3,j] > 136.8,4,
                    ifelse(ADD[j+4,j] > 136.8,5,
                           ifelse(ADD[j+5,j] > 136.8,6,
                                  7)))
  SDC[4,j] = ifelse(ADD[j+4,j] > 176.7,5,
                    ifelse(ADD[j+5,j] > 176.7,6,
                           ifelse(ADD[j+6,j] > 176.7,7,
                                  ifelse(ADD[j+7,j] > 176.7,8,
                                         9))))
  SDC[5,j] = ifelse(ADD[j+5,j] > 214.2,6,
                    ifelse(ADD[j+6,j] > 214.2,7,
                           ifelse(ADD[j+7,j] > 214.2,8,
                                  ifelse(ADD[j+8,j] > 214.2,9,
                                         10))))
  SDC[6,j] = ifelse(ADD[j+6,j] > 268.4,7,
                    ifelse(ADD[j+7,j] > 268.4,8,
                           ifelse(ADD[j+8,j] > 268.4,9,
                                  ifelse(ADD[j+9,j] > 268.4,10,
                                         ifelse(ADD[j+10,j] > 268.4,11,
                                                ifelse(ADD[j+11,j] > 268.4,12,
                                                       13))))))
  SDC[7,j] = ifelse(ADD[j+9,j] > 359.1,10,
                    ifelse(ADD[j+10,j] > 359.1,11,
                           ifelse(ADD[j+11,j] > 359.1,12,
                                  ifelse(ADD[j+12,j] > 359.1,13,
                                         ifelse(ADD[j+13,j] > 359.1,14,
                                                ifelse(ADD[j+14,j] > 359.1,15,
                                                       ifelse(ADD[j+15,j] > 359.1,16,
                                                              17)))))))
}


#########  calculate stage durations for each cohort and mean stage durations  ###############
#vector of mean stage durations
MSD = vector(length = 6)
#cohorts that finished 5th instar
I5 = ADD[(n+1),] > 359.1
#remove first 16 days before study started
I5[1:16] = FALSE
#stage durations of cohorts that reached pupation
SDC[7,I5]
#fifth instar stage durations for each cohort
SD5 = SDC[7,I5]-SDC[6,I5]
MSD[6] = mean(SD5); range(SD5)
#cohorts that finished 4th instar
I4 = ADD[(n+1),] > 268.4
I4[1:16] = FALSE
MSD[5] = mean(SDC[6,I4]-SDC[5,I4])
#cohorts that finished 3rd instar
I3 = ADD[(n+1),] > 214.2; I3[1:16] = FALSE
MSD[4] = mean(SDC[5,I3]-SDC[4,I3])
#cohorts that finished 2nd instar
I2 = ADD[(n+1),] > 176.7; I2[1:16] = FALSE
MSD[3] = mean(SDC[4,I2]-SDC[3,I2])
#cohorts that finished 1st instar
I1 = ADD[(n+1),] > 136.8; I1[1:16] = FALSE
MSD[2] = mean(SDC[3,I1]-SDC[2,I1])
#cohorts that hatched
I0 = ADD[(n+1),] > 90.2; I0[1:16] = FALSE
MSD[1] = mean(SDC[2,I0]-SDC[1,I0])
MSD #should this always be whole numbers?

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
#initial vlaues for lambda2 and S (S is denoted pi in the ms)
inits = list(list(lambda2 = 200*runif(n), S = runif(6)),
             list(lambda2 = 200*runif(n), S = runif(6)),
             list(lambda2 = 200*runif(n), S = runif(6)),
             list(lambda2 = 200*runif(n), S = runif(6)),
             list(lambda2 = 200*runif(n), S = runif(6)),
             list(lambda2 = 200*runif(n), S = runif(6)),
             list(lambda2 = 200*runif(n), S = runif(6)))

#runjags, monitoring daily survival, stage survival, and cumulative survival (SC).
out.T.140 = run.jags(model="MonarchModelSyriaca.bug", monitor = c("S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
plot(out.T.140, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
#add some more samples
out.T.140.e = extend.jags(out.T.140, sample = 190000)
plot(out.T.140.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
out.T.140.e2 = extend.jags(out.T.140.e, sample = 100000)
plot(out.T.140.e2, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
#save summaries
out.T.140.Sum = out.T.140.e2$summary
save(out.T.140.Sum, file = "FieldAnalysisSummaries.RData")




#### 260th Stree dataset from Bitzer et al. 2016 data, collected in 2015 #############
#IA ROW2 in ms supp info

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

#rerun ADD, SDC, MSD, then CB, data, and inits above for the data list

out.T.260 = run.jags(model="MonarchModelSyriaca.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
out.T.260.e = extend.jags(out.T.260, sample = 90000)
plot(out.T.260.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
out.T.260.e2 = extend.jags(out.T.260.e, sample = 100000)
plot(out.T.260.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))


