#Analysis of Monarch Field Data for Grant et al.

#code below just copied directly from SimulateandAnalyzeMonarchCounts.R and then adapted

library(runjags)
#option to always calculate summary statistics
runjags.options(force.summary=TRUE)

#First dataset is Bitzer et al. 2016 data from 140th Street, collected in 2015
#IA ROW1 in ms supp info

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
#initial vlaues for lambda2 and S (S is denoted pi in the ms)
inits = list(list(lambda2 = 200*runif(n), S = runif(6)),
             list(lambda2 = 200*runif(n), S = runif(6)),
             list(lambda2 = 200*runif(n), S = runif(6)),
             list(lambda2 = 200*runif(n), S = runif(6)),
             list(lambda2 = 200*runif(n), S = runif(6)),
             list(lambda2 = 200*runif(n), S = runif(6)),
             list(lambda2 = 200*runif(n), S = runif(6)))
#runjags, monitoring daily survival, stage survival, cumulative survival (SC), and eggs laid B.  
out.T.140 = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
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

#rerun ADD, CB, etc. above for the data list

out.T.260 = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
out.T.260.e = extend.jags(out.T.260, sample = 90000)
plot(out.T.260.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
out.T.260.e2 = extend.jags(out.T.260.e, sample = 100000)
plot(out.T.260.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))


####### Ontario sites ###############

########## AGR SITES #################

#Trillium - AGR1

n = 36+16
#maximum time to pupation, this shouldn't need changing, it makes indexing work
pupamax = 16
#mean temperature each day of the study period of length n days
MT = c(18.6,15,13,17.4,17.9,16.2,15.5,16.8,19.6,21.1,20.9,16.3,18.2,18.4,19.8,20.8,
       21.6,21.6,17.3,16.4,21,25.1,23.2,20.9,20.7,17.8,18.9,20.4,22.8,22.4,23.6,23.7,24.2,23.5,21.5,20,20.7,20.2,19.6,16.2,15.5,18.8,
       18.4,20,20,19.9,17,17.3,22.5,23.1,22.2,24.2)
length(MT)

C = matrix(nrow = 36, ncol = 6)
C[1,] = c(0,5,2,0,2,1)
C[8,] = c(47,5,1,1,0,2)
C[15,] = c(100,44,5,0,0,3)
C[22,] = c(52,31,19,8,5,2)
C[29,] = c(63,19,7,3,0,2)
C[36,] = c(26,13,5,1,4,0)

out.T.T = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
out.T.T.e = extend.jags(out.T.T, sample = 90000)
plot(out.T.T.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))

out.T.T.Sum = out.T.T.e$summaries
out.T.T.Sum2 = out.T.T.e$summary
save(out.T.T.Sum, out.T.T.Sum2, file = "SingleSiteSummaries1.RData")



#Vanmeer - AGR2
#same temps as for Trillium

C = matrix(nrow = 36, ncol = 6)
C[1,] = c(9,0,0,0,0,0)
C[8,] = c(12,3,4,0,0,0)
C[15,] = c(40,10,2,1,2,1)
C[22,] = c(81,29,8,5,1,1)
C[29,] = c(61,20,10,2,2,5)
C[36,] = c(25,9,11,1,1,0)

out.T.V = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
out.T.V.e = extend.jags(out.T.V, sample = 90000)
plot(out.T.V.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))

out.T.V.Sum = out.T.V.e$summaries
out.T.V.Sum2 = out.T.V.e$summary
save(out.T.T.Sum, out.T.T.Sum2, out.T.V.Sum, out.T.V.Sum2, file = "SingleSiteSummaries2.RData")



#Wooley - AGR3

n = 36+16
pupamax = 16
MT = c(17.9,16.2,15.5,16.8,19.6,21.1,20.9,16.3,18.2,18.4,19.8,20.8,21.6,21.6,17.3,16.4,21,25.1,23.2,20.9,20.7,17.8,18.9,20.4,22.8,22.4,
       23.6,23.7,24.2,23.5,21.5,20,20.7,20.2,19.6,16.2,15.5,18.8,18.4,20,20,19.9,17,17.3,22.5,23.1,22.2,24.2,23.1,23.7,19.8,17.3)
length(MT)

C = matrix(nrow = 36, ncol = 6)
C[1,] = c(13,0,0,0,0,1)
C[8,] = c(47,18,1,3,1,0)
C[15,] = c(67,5,12,7,10,0)
C[22,] = c(52,40,4,1,0,6)
C[29,] = c(19,14,18,11,1,5)
C[36,] = c(30,4,8,3,3,11)

out.T.W = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
plot(out.T.W, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
out.T.W.e = extend.jags(out.T.W, sample = 190000)
#out.T.W.e = extend.jags(out.T.W.e, sample = 100000)
plot(out.T.W.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))

out.T.W.Sum = out.T.W.e$summaries
out.T.W.Sum2 = out.T.W.e$summary
save(out.T.T.Sum, out.T.T.Sum2, out.T.V.Sum, out.T.V.Sum2, out.T.W.Sum, out.T.W.Sum2, file = "SingleSiteSummaries3.RData")



#Jan - AGR4

n = 36+16
pupamax = 16
MT = c(17,16.1,16.7,17.6,17.9,18.4,24.4,24.1,25.4,24.7,20.2,21.9,21.5,24.2,26.2,23.9,20.5,18.1,18,23.1,19.1,18.5,21.2,26.6,24.5,21.6,
       24.9,21.7,23,23.2,22.8,22.4,22,22,21.6,22.2,22.8,25.4,21.7,20.5,19.8,21.4,26,26.8,22,26.7,22.4,21.7,21.7,22,22.2,22.3)
length(MT)

C = matrix(nrow = 36, ncol = 6)
C[1,] = c(1,0,0,0,0,0)
C[8,] = c(7,0,0,0,0,0)
C[15,] = c(25,1,2,0,0,0)
C[22,] = c(13,2,3,0,0,0)
C[29,] = c(6,6,7,2,0,2)
C[36,] = c(3,2,1,0,1,3)
sum(C, na.rm = TRUE)

out.T.J = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
plot(out.T.J, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
out.T.J.e = extend.jags(out.T.J, sample = 90000)
plot(out.T.J.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))

out.T.J.Sum = out.T.J.e$summaries
out.T.J.Sum2 = out.T.J.e$summary
load("SingleSiteSummaries4.RData")
save(out.T.T.Sum, out.T.T.Sum2, out.T.V.Sum, out.T.V.Sum2, out.T.W.Sum, out.T.W.Sum2, out.T.J.Sum, out.T.J.Sum2,
     file = "SingleSiteSummaries4.RData")


#Chipps - AGR5 - not enough data - 15 data points total



#Middleton - AGR6

n = 36+16
pupamax = 16
MT = c(16.1,16.7,17.6,17.9,18.4,24.4,24.1,25.4,24.7,20.2,21.9,21.5,24.2,26.2,23.9,20.5,18.1,18,23.1,19.1,18.5,21.2,26.6,24.5,21.6,24.9,21.7,
       23,23.2,22.8,22.4,22,22,21.6,22.2,22.8,25.4,21.7,20.5,19.8,21.4,26,26.8,22,26.7,22.4,21.7,21.7,22,22.2,22.3,25)
length(MT)

C = matrix(nrow = 36, ncol = 6)
C[1,] = c(12,2,0,0,0,0)
C[8,] = c(8,0,4,1,0,0)
C[15,] = c(24,6,3,0,1,2)
C[22,] = c(69,17,12,2,1,1)
C[29,] = c(88,42,18,15,4,4)
C[36,] = c(83,17,18,10,9,9)
sum(C, na.rm = TRUE)

out.T.M = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
plot(out.T.M, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
out.T.M.e = extend.jags(out.T.M, sample = 90000)
plot(out.T.M.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))

out.T.M.Sum = out.T.M.e$summaries
out.T.M.Sum2 = out.T.M.e$summary
save(out.T.T.Sum, out.T.T.Sum2, out.T.V.Sum, out.T.V.Sum2, out.T.W.Sum, out.T.W.Sum2, out.T.J.Sum, out.T.J.Sum2, out.T.M.Sum, out.T.M.Sum2,
     file = "SingleSiteSummaries5.RData")



#Petker - AGR7

n = 36+16
pupamax = 16
MT = c(22,22.6,24.4,22,17,16.1,16.7,17.6,17.9,18.4,24.4,24.1,25.4,24.7,20.2,21.9,21.5,24.2,26.2,23.9,20.5,18.1,18,23.1,19.1,18.5,21.2,
       26.6,24.5,21.6,24.9,21.7,23,23.2,22.8,22.4,22,22,21.6,22.2,22.8,25.4,21.7,20.5,19.8,21.4,26,26.8,22,26.7,22.4,21.7)
length(MT)

C = matrix(nrow = 36, ncol = 6)
C[1,] = c(12,0,1,0,0,3)
C[8,] = c(5,0,4,0,0,2)
C[15,] = c(5,7,1,0,0,0)
C[22,] = c(99,0,1,1,2,0)
C[29,] = c(16,15,21,10,2,2)
C[36,] = c(16,4,3,3,10,10)
sum(C, na.rm = TRUE)

out.T.P = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
plot(out.T.P, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
out.T.P.e = extend.jags(out.T.P, sample = 90000)
plot(out.T.P.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))

#WOULDN'T CONVERGE - probably would with an informative prior, becuase of where S1 and S6 got stuck

out.T.P.Sum = out.T.P.e$summaries
out.T.P.Sum2 = out.T.P.e$summary
save(out.T.T.Sum, out.T.T.Sum2, out.T.V.Sum, out.T.V.Sum2, out.T.W.Sum, out.T.W.Sum2, out.T.J.Sum, out.T.J.Sum2, out.T.M.Sum, out.T.M.Sum2,
     out.T.P.Sum, out.T.P.Sum2, file = "SingleSiteSummaries6.RData")



##################### NAT SITES ###########################


#Anderson - NAT1

n = 37+16
pupamax = 16
MT = c(17.4,18.1,16.2,15.5,16.8,19.6,21.1,20.9,16.3,18.2,18.4,19.8,20.8,21.6,21.6,17.3,16.4,21,25.1,23.2,20.9,20.7,17.8,18.9,20.4,22.8,
       22.4,23.6,23.7,24.2,23.5,21.5,20,20.7,20.2,19.6,16.2,15.5,18.8,18.4,20,20,19.9,17,17.3,22.5,23.1,22.2,24.2,23.1,23.7,19.8,17.3)
length(MT)

C = matrix(nrow = 37, ncol = 6)
C[1,] = c(5,0,0,0,0,0)
C[8,] = c(8,1,0,0,0,1)
C[15,] = c(25,2,1,1,0,0)
C[22,] = c(38,11,2,0,0,0)
C[29,] = c(27,13,7,1,0,0)
C[37,] = c(31,10,6,5,2,1)

out.T.A = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
plot(out.T.A, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
out.T.A.e = extend.jags(out.T.A, sample = 90000)
plot(out.T.A.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))

out.T.A.Sum = out.T.A.e$summaries
out.T.A.Sum2 = out.T.A.e$summary
save(out.T.A.Sum, out.T.A.Sum2, file = "SingleSiteSummariesB1.RData")



#Helen - NAT2

n = 36+16
pupamax = 16
MT = c(16.2,15.5,16.8,19.6,21.1,20.9,16.3,18.2,18.4,19.8,20.8,21.6,21.6,17.3,16.4,21,25.1,23.2,20.9,20.7,17.8,18.9,20.4,22.8,22.4,23.6,
       23.7,24.2,23.5,21.5,20,20.7,20.2,19.6,16.2,15.5,18.8,18.4,20,20,19.9,17,17.3,22.5,23.1,22.2,24.2,23.1,23.7,19.8,17.3,17)
length(MT)

C = matrix(nrow = 36, ncol = 6)
C[1,] = c(20,0,2,0,0,0)
C[8,] = c(45,3,3,0,0,0)
C[15,] = c(22,8,3,1,3,2)
C[22,] = c(106,22,4,1,1,2)
C[29,] = c(39,20,19,5,1,1)
C[36,] = c(26,10,5,4,4,4)

out.T.H = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
plot(out.T.H, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
out.T.H.e = extend.jags(out.T.H, sample = 90000)
plot(out.T.H.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))

out.T.H.Sum = out.T.H.e$summaries
out.T.H.Sum2 = out.T.H.e$summary
save(out.T.A.Sum, out.T.A.Sum2, out.T.H.Sum, out.T.H.Sum2, file = "SingleSiteSummariesB1.RData") #saved as B1 instead of B2, oh well



#Sands - NAT3

n = 36+16
pupamax = 16
MT = c(17.4,18.1,16.2,15.5,16.8,19.6,21.1,20.9,16.3,18.2,18.4,19.8,20.8,21.6,21.6,17.3,16.4,21,25.1,23.2,20.9,20.7,17.8,18.9,20.4,22.8,22.4,
       23.6,23.7,24.2,23.5,21.5,20,20.7,20.2,19.6,16.2,15.5,18.8,18.4,20,20,19.9,17,17.3,22.5,23.1,22.2,24.2,23.1,23.7,19.8)
length(MT)

C = matrix(nrow = 36, ncol = 6)
C[1,] = c(3,0,0,0,0,0)
C[8,] = c(11,0,0,0,0,1)
C[15,] = c(16,0,1,0,0,0)
C[22,] = c(11,0,0,0,0,0)
C[29,] = c(30,3,3,1,0,0)
C[36,] = c(21,13,0,1,0,1)

out.T.S = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
plot(out.T.S, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
out.T.S.e = extend.jags(out.T.S, sample = 90000)
plot(out.T.S.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))

out.T.S.Sum = out.T.S.e$summaries
out.T.S.Sum2 = out.T.S.e$summary
save(out.T.A.Sum, out.T.A.Sum2, out.T.H.Sum, out.T.H.Sum2, out.T.S.Sum, out.T.S.Sum2, file = "SingleSiteSummariesB3.RData")



#Crane - NAT4

n = 36+16
pupamax = 16
MT = c(17,16.1,16.7,17.6,17.9,18.4,24.4,24.1,25.4,24.7,20.2,21.9,21.5,24.2,26.2,23.9,20.5,18.1,18,23.1,19.1,18.5,21.2,26.6,24.5,21.6,24.9,
       21.7,23,23.2,22.8,22.4,22,22,21.6,22.2,22.8,25.4,21.7,20.5,19.8,21.4,26,26.8,22,26.7,22.4,21.7,21.7,22,22.2,22.3)
length(MT)

C = matrix(nrow = 36, ncol = 6)
C[1,] = c(0,0,0,0,0,0)
C[8,] = c(1,0,0,0,0,0)
C[15,] = c(0,0,0,0,0,0)
C[22,] = c(40,5,0,0,0,0)
C[29,] = c(29,9,9,5,0,0)
C[36,] = c(24,6,8,1,1,4)

out.T.C = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
plot(out.T.C, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
out.T.C.e = extend.jags(out.T.C, sample = 90000)
plot(out.T.C.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))

out.T.C.Sum = out.T.C.e$summaries
out.T.C.Sum2 = out.T.C.e$summary
save(out.T.A.Sum, out.T.A.Sum2, out.T.H.Sum, out.T.H.Sum2, out.T.S.Sum, out.T.S.Sum2, out.T.C.Sum, out.T.C.Sum2, 
     file = "SingleSiteSummariesB4.RData")



#Falconer - NAT5

n = 36+16
pupamax = 16
MT = c(22,17,16.1,16.7,17.6,17.9,18.4,24.4,24.1,25.4,24.7,20.2,21.9,21.5,24.2,26.2,23.9,20.5,18.1,18,23.1,19.1,18.5,21.2,26.6,24.5,21.6,
       24.9,21.7,23,23.2,22.8,22.4,22,22,21.6,22.2,22.8,25.4,21.7,20.5,19.8,21.4,26,26.8,22,26.7,22.4,21.7,21.7,22,22.2)
length(MT)

C = matrix(nrow = 36, ncol = 6)
C[1,] = c(0,0,0,0,0,0)
C[8,] = c(4,0,0,0,0,0)
C[15,] = c(11,0,1,0,0,0)
C[22,] = c(6,6,0,0,0,0)
C[29,] = c(3,0,1,0,1,0)
C[36,] = c(15,0,0,2,0,1)

out.T.F = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
plot(out.T.F, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
out.T.F.e = extend.jags(out.T.F, sample = 90000)
plot(out.T.F.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))

out.T.F.Sum = out.T.F.e$summaries
out.T.F.Sum2 = out.T.F.e$summary
save(out.T.A.Sum, out.T.A.Sum2, out.T.H.Sum, out.T.H.Sum2, out.T.S.Sum, out.T.S.Sum2, out.T.C.Sum, out.T.C.Sum2, out.T.F.Sum, out.T.F.Sum2,
     file = "SingleSiteSummariesB5.RData")


#mackenzie - NAT6
#mccracken - NAT7
#pond - NAT8
#quarter - NAT9
#these 4 have no 5th instars observations



#storm - NAT10

n = 36+16
pupamax = 16
MT = c(22.6,24.4,22,17,16.1,16.7,17.6,17.9,18.4,24.4,24.1,25.4,24.7,20.2,21.9,21.5,24.2,26.2,23.9,20.5,18.1,18,23.1,19.1,18.5,21.2,26.6,24.5,
       21.6,24.9,21.7,23,23.2,22.8,22.4,22,22,21.6,22.2,22.8,25.4,21.7,20.5,19.8,21.4,26,26.8,22,26.7,22.4,21.7,21.7)
length(MT)

C = matrix(nrow = 36, ncol = 6)
C[1,] = c(2,0,0,0,0,0)
C[8,] = c(1,0,0,1,0,0)
C[15,] = c(5,0,0,0,0,0)
C[22,] = c(10,4,3,0,0,1)
C[29,] = c(0,8,2,1,0,0)
C[36,] = c(22,1,2,1,4,3)

out.T.ST = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
plot(out.T.ST, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
out.T.ST.e = extend.jags(out.T.ST, sample = 90000)
plot(out.T.ST.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))

out.T.ST.Sum = out.T.ST.e$summaries
out.T.ST.Sum2 = out.T.ST.e$summary
save(out.T.A.Sum, out.T.A.Sum2, out.T.H.Sum, out.T.H.Sum2, out.T.S.Sum, out.T.S.Sum2, out.T.C.Sum, out.T.C.Sum2, out.T.F.Sum, out.T.F.Sum2,
     out.T.ST.Sum, out.T.ST.Sum2, file = "SingleSiteSummariesB6.RData")



############# ROW SITES  #####################

#eight - ROW1

n = 36+16
pupamax = 16
MT = c(18.1,16.2,15.5,16.8,19.6,21.1,20.9,16.3,18.2,18.4,19.8,20.8,21.6,21.6,17.3,16.4,21,25.1,23.2,20.9,20.7,17.8,18.9,20.4,22.8,22.4,23.6,
       23.7,24.2,23.5,21.5,20,20.7,20.2,19.6,16.2,15.5,18.8,18.4,20,20,19.9,17,17.3,22.5,23.1,22.2,24.2,23.1,23.7,19.8,17.3)
length(MT)

C = matrix(nrow = 36, ncol = 6)
C[1,] = c(7,0,0,0,0,1)
C[8,] = c(4,2,1,0,0,0)
C[15,] = c(19,1,1,1,0,0)
C[22,] = c(17,5,2,0,0,1)
C[29,] = c(1,3,1,0,0,2)
C[36,] = c(5,1,0,0,0,0)

out.T.R = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
plot(out.T.R, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
out.T.R.e = extend.jags(out.T.R, sample = 90000)
plot(out.T.R.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))

out.T.R.Sum = out.T.R.e$summaries
out.T.R.Sum2 = out.T.R.e$summary
save(out.T.R.Sum, out.T.R.Sum2, file = "SingleSiteSummariesC1.RData")



#six - ROW2

n = 36+16
pupamax = 16
MT = c(13,17.4,18.1,16.2,15.5,16.8,19.6,21.1,20.9,16.3,18.2,18.4,19.8,20.8,21.6,21.6,17.3,16.4,21,25.1,23.2,20.9,20.7,17.8,18.9,20.4,22.8,
       22.4,23.6,23.7,24.2,23.5,21.5,20,20.7,20.2,19.6,16.2,15.5,18.8,18.4,20,20,19.9,17,17.3,22.5,23.1,22.2,24.2,23.1,23.7)
length(MT)

C = matrix(nrow = 36, ncol = 6)
C[1,] = c(0,0,0,0,0,0)
C[8,] = c(1,0,0,0,0,0)
C[15,] = c(10,0,0,0,0,0)
C[22,] = c(11,1,1,0,0,0)
C[29,] = c(5,2,0,1,0,0)
C[36,] = c(0,4,0,1,0,1)

out.T.SI = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
plot(out.T.SI, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
out.T.SI.e = extend.jags(out.T.SI, sample = 90000)
plot(out.T.SI.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))

out.T.SI.Sum = out.T.SI.e$summaries
out.T.SI.Sum2 = out.T.SI.e$summary
save(out.T.R.Sum, out.T.R.Sum2, out.T.SI.Sum, out.T.SI.Sum2, file = "SingleSiteSummariesC2.RData")



#three - ROW3 - no 5ths

#causeway - ROW4

n = 36+16
pupamax = 16
MT = c(24.4,22,17,16.1,16.7,17.6,17.9,18.4,24.4,24.1,25.4,24.7,20.2,21.9,21.5,24.2,26.2,23.9,20.5,18.1,18,23.1,19.1,18.5,21.2,26.6,24.5,21.6,
       24.9,21.7,23,23.2,22.8,22.4,22,22,21.6,22.2,22.8,25.4,21.7,20.5,19.8,21.4,26,26.8,22,26.7,22.4,21.7,21.7,22)
length(MT)

C = matrix(nrow = 36, ncol = 6)
C[1,] = c(11,1,1,0,0,0)
C[8,] = c(5,3,2,0,3,1)
C[15,] = c(3,2,0,0,0,3)
C[22,] = c(4,2,1,2,0,0)
C[29,] = c(2,10,3,0,0,1)
C[36,] = c(4,17,3,3,4,0)
sum(C, na.rm = TRUE)

out.T.CA = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
plot(out.T.CA, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
out.T.CA.e = extend.jags(out.T.CA, sample = 90000)
plot(out.T.CA.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))

out.T.CA.Sum = out.T.CA.e$summaries
out.T.CA.Sum2 = out.T.CA.e$summary
save(out.T.R.Sum, out.T.R.Sum2, out.T.SI.Sum, out.T.SI.Sum2, out.T.CA.Sum, out.T.CA.Sum2, file = "SingleSiteSummariesC3.RData")


#cornturn - ROW5
#crypt - ROW6
#neither have 5th instars


#peacock - ROW7
n = 36+16
pupamax = 16
MT = c(22.6,24.4,22,17,16.1,16.7,17.6,17.9,18.4,24.4,24.1,25.4,24.7,20.2,21.9,21.5,24.2,26.2,23.9,20.5,18.1,18,23.1,19.1,18.5,21.2,26.6,24.5,
       21.6,24.9,21.7,23,23.2,22.8,22.4,22,22,21.6,22.2,22.8,25.4,21.7,20.5,19.8,21.4,26,26.8,22,26.7,22.4,21.7,21.7)
length(MT)

C = matrix(nrow = 36, ncol = 6)
C[1,] = c(0,0,0,0,0,0)
C[8,] = c(0,0,0,0,0,0)
C[15,] = c(0,1,0,0,0,0)
C[22,] = c(15,3,0,0,0,0)
C[29,] = c(2,6,1,2,0,0)
C[36,] = c(0,0,2,2,1,1)
sum(C, na.rm = TRUE)

out.T.P = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
plot(out.T.P, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
out.T.P.e = extend.jags(out.T.P, sample = 90000)
plot(out.T.P.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))

out.T.P.Sum = out.T.P.e$summaries
out.T.P.Sum2 = out.T.P.e$summary
save(out.T.R.Sum, out.T.R.Sum2, out.T.SI.Sum, out.T.SI.Sum2, out.T.CA.Sum, out.T.CA.Sum2, out.T.P.Sum, out.T.P.Sum2, 
     file = "SingleSiteSummariesC4.RData")


#str2 - ROW8

n = 36+16
pupamax = 16
MT = c(22,17,16.1,16.7,17.6,17.9,18.4,24.4,24.1,25.4,24.7,20.2,21.9,21.5,24.2,26.2,23.9,20.5,18.1,18,23.1,19.1,18.5,21.2,26.6,24.5,21.6,
       24.9,21.7,23,23.2,22.8,22.4,22,22,21.6,22.2,22.8,25.4,21.7,20.5,19.8,21.4,26,26.8,22,26.7,22.4,21.7,21.7,22,22.2)
length(MT)

C = matrix(nrow = 36, ncol = 6)
C[1,] = c(3,0,0,0,0,0)
C[8,] = c(5,0,3,0,0,0)
C[15,] = c(10,3,0,0,0,1)
C[22,] = c(9,2,1,0,0,0)
C[29,] = c(15,5,3,0,0,1)
C[36,] = c(4,5,3,3,0,4)
sum(C, na.rm = TRUE)

out.T.STR = run.jags(model="MonarchModel.bug", monitor = c("B","S","S1","S2","S3","S4","S5","S6","SC"), n.chains = 7, data = data, inits = inits)
plot(out.T.STR, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))
out.T.STR.e = extend.jags(out.T.STR, sample = 90000)
plot(out.T.STR.e, vars = c("S","S1","S2","S3","S4","S5","S6","SC"))

out.T.STR.Sum = out.T.STR.e$summaries
out.T.STR.Sum2 = out.T.STR.e$summary
save(out.T.R.Sum, out.T.R.Sum2, out.T.SI.Sum, out.T.SI.Sum2, out.T.CA.Sum, out.T.CA.Sum2, out.T.P.Sum, out.T.P.Sum2, 
     out.T.STR.Sum, out.T.STR.Sum2, file = "SingleSiteSummariesC5.RData")



##################### Meta Analysis #########################

load("FieldAnalysisSummaries4.RData")
load("SingleSiteSummaries5.RData")
load("SingleSiteSummariesB6.RData")
load("SingleSiteSummariesC5.RData")

#dataframe with parm ests and 95% CI

#list of ests to loop through
#for some reason runjags has 2 summary types, both of which have different things
#forgot to get both types for out.T.140, so using the same type as out.T.140.sum
EstsSum2 = list(
                out.T.140.Sum, out.T.260.Sum2, #IA
                out.T.T.Sum2, out.T.V.Sum2, out.T.W.Sum2, out.T.J.Sum2, out.T.M.Sum2, #AGR
                out.T.A.Sum2, out.T.H.Sum2, out.T.S.Sum2, out.T.C.Sum2, out.T.F.Sum2, out.T.ST.Sum2, #NAT
                out.T.R.Sum2, out.T.SI.Sum2, out.T.CA.Sum2, out.T.P.Sum2, out.T.STR.Sum2 #ROW
                )
length(EstsSum2)

Reg = data.frame(matrix(NA, nrow = 18, ncol = 10))
colnames(Reg) = c("SiteNo","SiteName","Cat","S","LCL","UCL","d1","d2","d3")
Reg[,1] = 1:18
Reg[,2] = c("140th","260th","trillium","vanmeer","wooley","jan","middleton","anderson","helen","sands","crane","falconer","storm",
            "eight","six","causeway","peacock","str2")
Reg[,3] = c("IA","IA","AGR","AGR","AGR","AGR","AGR","NAT","NAT","NAT","NAT","NAT","NAT","ROW","ROW","ROW","ROW","ROW")
Reg[,7] = c(rep(1,13),rep(0,5))
Reg[,8] = c(rep(1,7),rep(0,11))
Reg[,9] = c(rep(1,2),rep(0,16))

for (i in 1:18){
  Reg[i,4] = EstsSum2[[i]]$quantiles[length(EstsSum2[[i]]$quantiles[,3]),3]
  Reg[i,5] = EstsSum2[[i]]$quantiles[length(EstsSum2[[i]]$quantiles[,3]),1]
  Reg[i,6] = EstsSum2[[i]]$quantiles[length(EstsSum2[[i]]$quantiles[,3]),5]
}

mean(Reg[,4]) #0.0143, same as est from model below, even though the weighting is slightly different
median(Reg[,4]) #0.009 - lower than mean because of skewed distr


#graph of all individual ests.  
library(ggplot2)

p = ggplot(Reg, aes(x = factor(Reg$SiteName, levels = Reg$SiteName), y = S, fill = Cat)) +
  geom_bar(stat = "identity") + xlab("Site") + ylab("Survival Probability") + labs(fill = "Type") +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width=0.2) +
  scale_x_discrete(labels = c("IAROW1","IAROW2","ONAGR1","ONAGR2","ONAGR3","ONAGR4","ONAGR6","ONAGR7","ONNAT1","ONNAT2","ONNAT3","ONNAT4","ONNAT5",
                              "ONNAT10","ONROW1","ONROW2","ONROW4","ONROW7","ONROW8")) +
  geom_text(aes(label = round(S, digits = 3)), vjust = -0.3) +
  theme(axis.title = element_text(size = 20))
p



###################  Meta Regression  ##################################

data = list(Pi = Reg$S, d1=Reg$d1, d2=Reg$d2, d3=Reg$d3)
#initial vlaues for lambda2 and S
inits = list(list(a = runif(1), b1 = runif(1), b2 = runif(1), b3 = runif(1)),
             list(a = runif(1), b1 = runif(1), b2 = runif(1), b3 = runif(1)),
             list(a = runif(1), b1 = runif(1), b2 = runif(1), b3 = runif(1)))
#runjags, monitoring daily survival, stage survival, cumulative survival (SC), and eggs laid B.  
out.reg = run.jags(model="SecondStage2.bug", monitor = c("a","b1","b2","b3","prec","predicted_mean_IA","predicted_mean_AGR",
                  "predicted_mean_NAT","predicted_mean_ROW","predicted_mean_ALL"), n.chains = 3, data = data, inits = inits)
out.reg
plot(out.reg)


#plot means
#into df
out.df = data.frame(out.reg$summaries)
means = out.df[6:10,]
means[,12] = c("IA ROW","ON AGR","ON NAT","ON ROW","ALL")
colnames(means)[12] = "Cat"
means[,13] = c("Right of Ways","Agricultural Borders","Natural Areas","Right of Ways","All")
colnames(means)[13] = "Cat2"
means[,14] = c("Iowa","Ontario","Ontario","Ontario","Iowa and Ontario")
colnames(means)[14] = "Cat3"
means$Cat4 = factor(means$Cat3, levels = c("Iowa","Ontario","Iowa and Ontario")) #sets order of x-axis in graph

#with 95% CI
p2 = ggplot(means, aes(x = factor(means$Cat, levels = means$Cat), y = means$Median)) +
  geom_bar(stat = "identity") + xlab("Category") + ylab("Mean Survival Probability") +
  geom_errorbar(aes(ymin = means$Lower95, ymax = means$Upper95), width = 0.2) +
  geom_text(aes(label = round(means$Median, digits = 4), vjust = -0.3))
p2

#with SD - used in ms
p2 = ggplot(means, aes(x = Cat2, y = means$Median)) +
  facet_wrap(~Cat4, strip.position = "bottom", scales = "free_x") +
  geom_bar(stat = "identity") + xlab("Landcover Type") + ylab("Mean Survival Probability") +
  geom_errorbar(aes(ymin = means$Median-means$SD, ymax = means$Median+means$SD), width = 0.2) +
  geom_text(aes(label = format(round(means$Median, digits = 3))), vjust = -0.5, size = 5, hjust = 1.5) +
  geom_hline(yintercept = 0) +
  theme(axis.title.x = element_text(size = 20, margin = margin(t=10)),
        axis.title.y = element_text(size = 20, margin = margin(r=10)),
        axis.text = element_text(size = 16),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        #panel.spacing = unit(0, "lines"),
        #strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size = 16))
p2

library(grid)
gtable::gtable_show_layout(gt)

gt = ggplot_gtable(ggplot_build(p2))
gt$widths[9] = 3*gt$widths[5]
grid.draw(gt)




########## Eggs Per Stem Analysis ###########################
#calculate monarchs per X stems, easier to understand

#vector of eggs per MW - IA, AGR_15, AGR_16, NAT_15, NAT_16, ROW_15, ROW_16 - from spreadsheet Grace Pitman sent, and Blader data
EP = c(0.433, 0.224, 0.051, 0.239, 0.055, 0.076, 0.015)
#vector of egg density - IA, AGR_15, AGR_16, NAT_15, NAT_16, ROW_15, ROW_16
meanallEP = c(EP[1],mean(EP[2:3]),mean(EP[4:5]),mean(EP[6:7]))
#number of stems - calculates adults produced for this number of stems
St = 1
St = 1/(EP[7]*out.df$Median[4+5]) #scaled so that smallest production is 1 adult per 5159.939 stems

PerStem = data.frame(matrix(nrow = 11, ncol = 7))
colnames(PerStem) = c("Adults","LCL","UCL","LowerSD","UpperSD","Cat")
PerStem$Cat = c("IA_15","AGR_15","AGR_16","NAT_15","NAT_16","ROW_15","ROW_16","AGR","NAT","ROW","ALL")

#Adults per St = 1000 stems
PerStem[1,1] = St*EP[1]*out.df$Median[1+5]
PerStem[2,1] = St*EP[2]*out.df$Median[2+5]
PerStem[3,1] = St*EP[3]*out.df$Median[2+5]
PerStem[4,1] = St*EP[4]*out.df$Median[3+5]
PerStem[5,1] = St*EP[5]*out.df$Median[3+5]
PerStem[6,1] = St*EP[6]*out.df$Median[4+5]
PerStem[7,1] = St*EP[7]*out.df$Median[4+5]
PerStem[8,1] = St*mean(EP[2:3])*out.df$Median[2+5]
PerStem[9,1] = St*mean(EP[4:5])*out.df$Median[3+5]
PerStem[10,1] = St*mean(EP[6:7])*out.df$Median[4+5]
PerStem[11,1] = St*mean(meanallEP)*out.df$Median[5+5]

###### LowerSD ########
PerStem[1,4] = St*EP[1]*(out.df$Median[1+5]-out.df$SD[1+5]) #IA_15
PerStem[2,4] = St*EP[2]*(out.df$Median[2+5]-out.df$SD[2+5]) #AGR_15
PerStem[3,4] = St*EP[3]*(out.df$Median[2+5]-out.df$SD[2+5]) #AGR_16
PerStem[4,4] = St*EP[4]*(out.df$Median[3+5]-out.df$SD[3+5]) #NAT_15
PerStem[5,4] = St*EP[5]*(out.df$Median[3+5]-out.df$SD[3+5]) #NAT_16
PerStem[6,4] = St*EP[6]*(out.df$Median[4+5]-out.df$SD[4+5]) #ROW_15
PerStem[7,4] = St*EP[7]*(out.df$Median[4+5]-out.df$SD[4+5]) #ROW_16
PerStem[8,4] = St*mean(EP[2:3])*(out.df$Median[2+5]-out.df$SD[2+5]) #AGR
PerStem[9,4] = St*mean(EP[4:5])*(out.df$Median[3+5]-out.df$SD[3+5]) #NAT
PerStem[10,4] = St*mean(EP[6:7])*(out.df$Median[4+5]-out.df$SD[4+5]) #ROW
PerStem[11,4] = St*mean(meanallEP)*(out.df$Median[5+5]-out.df$SD[5+5]) #All

###### UpperSD ########
PerStem[1,5] = St*EP[1]*(out.df$Median[1+5]+out.df$SD[1+5]) #IA_15
PerStem[2,5] = St*EP[2]*(out.df$Median[2+5]+out.df$SD[2+5]) #AGR_15
PerStem[3,5] = St*EP[3]*(out.df$Median[2+5]+out.df$SD[2+5]) #AGR_16
PerStem[4,5] = St*EP[4]*(out.df$Median[3+5]+out.df$SD[3+5]) #NAT_15
PerStem[5,5] = St*EP[5]*(out.df$Median[3+5]+out.df$SD[3+5]) #NAT_16
PerStem[6,5] = St*EP[6]*(out.df$Median[4+5]+out.df$SD[4+5]) #ROW_15
PerStem[7,5] = St*EP[7]*(out.df$Median[4+5]+out.df$SD[4+5]) #ROW_16
PerStem[8,5] = St*mean(EP[2:3])*(out.df$Median[2+5]+out.df$SD[2+5]) #AGR
PerStem[9,5] = St*mean(EP[4:5])*(out.df$Median[3+5]+out.df$SD[3+5]) #NAT
PerStem[10,5] = St*mean(EP[6:7])*(out.df$Median[4+5]+out.df$SD[4+5]) #ROW
PerStem[11,5] = St*mean(meanallEP)*(out.df$Median[5+5]+out.df$SD[5+5]) #All

#should make a column just for SD
St*EP[2]*(out.df$SD[2+5]) #SD AGR_15
St*EP[7]*(out.df$SD[4+5]) #SD ROW_16

#just graph yearxlandcover type - figure used in ms
PerStemB = PerStem[1:7,]
PerStemB$Year = c("2015","2015","2016","2015","2016","2015","2016")
PerStemB$Cat2 = factor(PerStemB$Cat, levels = c("IA_15","AGR_15","AGR_16","NAT_15","NAT_16","ROW_15","ROW_16"))
PerStemB = PerStemB[-c(2,3,7)]
PerStemB$Cat3 = c("Iowa","Ontario","Ontario","Ontario","Ontario","Ontario","Ontario")
PerStemB$Adults2 = as.character(round(PerStemB$Adults, digits = 1))
PerStemB$Adults2[7] = "1.0" #had to do this to get ending 0 to display, otherwise just displayed "1" and wouldn't get spacing right as data label
PerStemB$Adults2[1] = "24.0"
PerStemB$Cat4 = c("Right of Ways","Agricultural Borders","Agricultural Borders","Natural Areas","Natural Areas","Right of Ways","Right of Ways")


library(RColorBrewer)
display.brewer.all()

p3 = ggplot(PerStemB, aes(x = Cat2, y = Adults, fill = Year)) +
  facet_wrap(~Cat3, strip.position = "bottom", scales = "free_x")+
  geom_bar(stat = "identity") + xlab("Landcover Type") + ylab("Adults per 5,160 Milkweed Stems") +
  geom_errorbar(aes(ymin = PerStemB$LowerSD, ymax = PerStemB$UpperSD), width = 0.2) +
  geom_text(aes(label = Adults2, vjust = -0.5, size = 5, hjust = 2), show.legend = FALSE) +
  scale_x_discrete(labels = c("Right of Ways","Agricultural\nBorders","Agricultural\nBorders","Natural\nAreas","Natural\nAreas","Right of Ways","Right of Ways")) +
  geom_hline(yintercept = 0) +
  scale_fill_brewer(palette = "Set1") +
  #guides(fill = FALSE) +
  theme(axis.title.x = element_text(size = 20, margin = margin(t=10)),
        axis.title.y = element_text(size = 20, margin = margin(r=10)),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        legend.position = c(0.85,0.85),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        strip.placement = "outside",
        strip.text.x = element_text(size = 16)
  )

p3

gt = ggplot_gtable(ggplot_build(p3))
gt$widths[9] = 5*gt$widths[5]
grid.draw(gt)




#with SD - used in ms
p2 = ggplot(means, aes(x = Cat2, y = means$Median)) +
  facet_wrap(~Cat4, strip.position = "bottom", scales = "free_x") +
  geom_bar(stat = "identity") + xlab("Landcover Type") + ylab("Mean Survival Probability") +
  geom_errorbar(aes(ymin = means$Median-means$SD, ymax = means$Median+means$SD), width = 0.2) +
  geom_text(aes(label = format(round(means$Median, digits = 3))), vjust = -0.5, size = 5, hjust = 1.5) +
  geom_hline(yintercept = 0) +
  theme(axis.title.x = element_text(size = 20, margin = margin(t=10)),
        axis.title.y = element_text(size = 20, margin = margin(r=10)),
        axis.text = element_text(size = 16),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        #panel.spacing = unit(0, "lines"),
        #strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size = 16))
p2

library(grid)
gtable::gtable_show_layout(gt)

gt = ggplot_gtable(ggplot_build(p2))
gt$widths[9] = 3*gt$widths[5]
grid.draw(gt)







