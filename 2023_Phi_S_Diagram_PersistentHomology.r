################################################################################
# Method for the generation of a Phi_S diagram based persistent Homology
# Author: Juan G. Diaz Ochoa
# ID: 
# March - 2023
################################################################################

# Download libraries
library("ggplot2")
library("TDA")
library(readr)
local <- 0
if(local == 0){
  fn <- "http://archive.ics.uci.edu/ml/machine-learning-databases/00319/MHEALTHDATASET.zip"
  download.file(fn,destfile="zip")
  unzip("zip",list = TRUE)
  subj1 <- read.table(unzip("zip",files = "MHEALTHDATASET/mHealth_subject1.log"))
  subj2 <- read.table(unzip("zip",files = "MHEALTHDATASET/mHealth_subject2.log"))
  subj3 <- read.table(unzip("zip",files = "MHEALTHDATASET/mHealth_subject3.log"))
  subj4 <- read.table(unzip("zip",files = "MHEALTHDATASET/mHealth_subject4.log"))
  subj5 <- read.table(unzip("zip",files = "MHEALTHDATASET/mHealth_subject5.log"))
  subj6 <- read.table(unzip("zip",files = "MHEALTHDATASET/mHealth_subject6.log"))
  subj7 <- read.table(unzip("zip",files = "MHEALTHDATASET/mHealth_subject7.log"))
  subj8 <- read.table(unzip("zip",files = "MHEALTHDATASET/mHealth_subject8.log"))
}
# Synthetic random data - no coupled
if(local == 2){
  subj1 <- as.data.frame(cbind(runif(161280),runif(161280),runif(161280),runif(161280),runif(161280)))
  subj2 <- as.data.frame(cbind(runif(161280),runif(161280),runif(161280),runif(161280),runif(161280)))
  subj3 <- as.data.frame(cbind(runif(161280),runif(161280),runif(161280),runif(161280),runif(161280)))
  subj4 <- as.data.frame(cbind(runif(161280),runif(161280),runif(161280),runif(161280),runif(161280)))
  subj5 <- as.data.frame(cbind(runif(161280),runif(161280),runif(161280),runif(161280),runif(161280)))
  subj6 <- as.data.frame(cbind(runif(161280),runif(161280),runif(161280),runif(161280),runif(161280)))
  subj7 <- as.data.frame(cbind(runif(161280),runif(161280),runif(161280),runif(161280),runif(161280)))
  subj8 <- as.data.frame(cbind(runif(161280),runif(161280),runif(161280),runif(161280),runif(161280)))
}
# Synthetic data: different combination of semi-periodic time series - No real coupled model
if(local == 3){
  synthetic <- read_table2("RootDirectory/synthetic.txt")
  subj1 <- as.data.frame(cbind(synthetic$`0.0000000e+000`,synthetic$`0.0000000e+000_1`,synthetic$`0.0000000e+000_3`,synthetic$`0.0000000e+000_4`,synthetic$`0.0000000e+000_5`))
  subj2 <- as.data.frame(cbind(synthetic$`0.0000000e+000_1`,synthetic$`0.0000000e+000`,synthetic$`0.0000000e+000_2`,synthetic$`0.0000000e+000_3`,synthetic$`0.0000000e+000_4`))
  subj3 <- as.data.frame(cbind(synthetic$`0.0000000e+000_5`,synthetic$`0.0000000e+000_4`,synthetic$`0.0000000e+000_3`,synthetic$`0.0000000e+000_2`,synthetic$`0.0000000e+000_1`))
  subj4 <- as.data.frame(cbind(synthetic$`0.0000000e+000_2`,synthetic$`0.0000000e+000_3`,synthetic$`0.0000000e+000_4`,synthetic$`0.0000000e+000`,synthetic$`0.0000000e+000_5`))
  subj5 <- as.data.frame(cbind(synthetic$`0.0000000e+000_3`,synthetic$`0.0000000e+000_2`,synthetic$`0.0000000e+000_1`,synthetic$`0.0000000e+000`,synthetic$`0.0000000e+000_5`))
  subj6 <- as.data.frame(cbind(synthetic$`0.0000000e+000_4`,synthetic$`0.0000000e+000`,synthetic$`0.0000000e+000_1`,synthetic$`0.0000000e+000_2`,synthetic$`0.0000000e+000_3`))
  subj7 <- as.data.frame(cbind(synthetic$`0.0000000e+000_3`,synthetic$`0.0000000e+000_2`,synthetic$`0.0000000e+000_3`,synthetic$`0.0000000e+000`,synthetic$`0.0000000e+000_5`))
  subj8 <- as.data.frame(cbind(synthetic$`0.0000000e+000_5`,synthetic$`0.0000000e+000_4`,synthetic$`0.0000000e+000_3`,synthetic$`0.0000000e+000_2`,synthetic$`0.0000000e+000_1`))
}
# Define lists from the downloaded data
SP <- list()
SP[[1]] <- subj1
SP[[2]] <- subj2
SP[[3]] <- subj3
SP[[4]] <- subj4
SP[[5]] <- subj5
SP[[6]] <- subj6
SP[[7]] <- subj7
SP[[8]] <- subj8
NumberPatients <- 8
#Generate a grid for Analysis using homology persistence
Xlim <- c(-10, 10); Ylim <- c(-10, 10); by <- 0.03
Xseq <- seq(Xlim[1], Xlim[2], by = by)
Yseq <- seq(Ylim[1], Ylim[2], by = by)
h <- 0.02
Grid <- expand.grid(Xseq, Yseq)
# Select axis where the computation should be performed
axis <- 3
# Select the axis from the M_Health data where the evaluation should be performed
Column_ECG <- 4
#Define time periods for the evaluation of causality
TP <- c()
TP[1] <- 1
TP[2] <- 25000
#Second time period
TP[3] <- 20000
TP[4] <- 45000
# Third time period
TP[5] <- 40000#115000
TP[6] <- 60000#119000 

XP <- list()
XP_TP <- list()
XP_T <- list()
# introduce axis where this analysis should be performed
if(axis == 3){
  s_pa <- axis
  s_pb <- Column_ECG 
}
if(axis == 2){
  s_pa <- axis
  s_pb <- Column_ECG 
}
if(axis == 1){
  s_pa <- axis
  s_pb <- Column_ECG 
}
# Definition of the different time periods
# Iteration over all the individuals - axes are normalized
# XP is the normalized trajectory combining two different parameters

for(i in 1:NumberPatients){
  XP[[1]] <- cbind(SP[[i]][[s_pa]][TP[1]:TP[2]]/max(abs(SP[[i]][[s_pa]][TP[1]:TP[2]])),SP[[i]][[s_pb]][TP[1]:TP[2]]/max(abs(SP[[i]][[s_pb]][TP[1]:TP[2]]))) 
  XP[[2]] <- cbind(SP[[i]][[s_pa]][TP[3]:TP[4]]/max(abs(SP[[i]][[s_pa]][TP[3]:TP[4]])),SP[[i]][[s_pb]][TP[3]:TP[4]]/max(abs(SP[[i]][[s_pb]][TP[3]:TP[4]])))
  XP[[3]] <- cbind(SP[[i]][[s_pa]][TP[5]:TP[6]]/max(abs(SP[[i]][[s_pa]][TP[5]:TP[6]])),SP[[i]][[s_pb]][TP[5]:TP[6]]/max(abs(SP[[i]][[s_pb]][TP[5]:TP[6]])))
  
  XP_TP[[i]] <- list(XP[[1]],XP[[2]],XP[[3]])
}
# Define computation in parallel
library(doParallel)
library(foreach)
DiagPH <- list()
DiagPH_T <- list()
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)
# Compute Gamma^2
gamma <- 0
if(gamma ==  0){
  for(i in 1:NumberPatients){
    DiagPH[[1]] <- gridDiag(X = XP_TP[[i]][[1]], FUN = kde, h = h, lim = cbind(Xlim, Ylim),by = by, sublevel = FALSE, library = "Dionysus",printProgress = TRUE)
    DiagPH[[2]] <- gridDiag(X = XP_TP[[i]][[2]], FUN = kde, h = h, lim = cbind(Xlim, Ylim),by = by, sublevel = FALSE, library = "Dionysus",printProgress = TRUE)
    DiagPH[[3]] <- gridDiag(X = XP_TP[[i]][[3]], FUN = kde, h = h, lim = cbind(Xlim, Ylim),by = by, sublevel = FALSE, library = "Dionysus",printProgress = TRUE)
    DiagPH_T[[i]] <- list(DiagPH[[1]],DiagPH[[2]],DiagPH[[3]])
  }
}

if(gamma ==  1){
  for(i in 1:NumberPatients){
    DiagPH[[1]] <- gridDiag(X = XP_TP[[i]][[1]]*XP_TP[[i]][[1]], FUN = kde, h = h, lim = cbind(Xlim, Ylim),by = by, sublevel = FALSE, library = "Dionysus",printProgress = TRUE)
    DiagPH[[2]] <- gridDiag(X = XP_TP[[i]][[2]]*XP_TP[[i]][[2]], FUN = kde, h = h, lim = cbind(Xlim, Ylim),by = by, sublevel = FALSE, library = "Dionysus",printProgress = TRUE)
    DiagPH[[3]] <- gridDiag(X = XP_TP[[i]][[3]]*XP_TP[[i]][[3]], FUN = kde, h = h, lim = cbind(Xlim, Ylim),by = by, sublevel = FALSE, library = "Dionysus",printProgress = TRUE)
    DiagPH_T[[i]] <- list(DiagPH[[1]],DiagPH[[2]],DiagPH[[3]])
  }
}
# This threshold parameter filters all the bars that deploy noise
NL <- 0.05 #0.0005
Diff_HG <- list()
Diff_HGT <- list()
Diff_HGT2 <- list()
Perst_T1 <-list()
Perst_T2 <- list()
Perst_T3 <- list()
for(i in 1:NumberPatients){
  # The difference of the persistence bars will be computed for a single individual - analysis of time series for the single individual: for this reason i = j
  j <- i
  auxp_i <- 1
  auxp_j <- 2
  auxp_k <- 3
  
  # This one is the estimation of the persistence bars
  #Bar <- Bar/max(Bar)
  Bar <- DiagPH_T[[i]][[auxp_i]]$diagram[,2] - DiagPH_T[[i]][[auxp_i]]$diagram[,3]
  #Bar <- unique(ifelse(abs(Bar) < NL,Bar <- Bar,Bar<- NA))
  Bar1 <- DiagPH_T[[j]][[auxp_j]]$diagram[,2] - DiagPH_T[[j]][[auxp_j]]$diagram[,3]
  Bar1 <- Bar1/max(Bar1)
  #Bar1 <- unique(ifelse(abs(Bar) < NL,Bar <- Bar,Bar<- NA))
  #Bar1 <- unique(ifelse(abs(Bar1) < NL,Bar1 <- Bar1,Bar1<- NA))
  Bar2 <- DiagPH_T[[j]][[auxp_k]]$diagram[,2] - DiagPH_T[[j]][[auxp_k]]$diagram[,3]
  #Bar2 <- unique(ifelse(abs(Bar2) < NL,Bar2 <- Bar2,Bar1<- NA))  
  Bar2 <- Bar2/max(Bar2)
  #Different time periods same individual
  
  aux_1 <- (Bar[1:15] - Bar1[1:15]) 
  aux_2 <- (Bar1[1:15] - Bar2[1:15]) 
  aux_1 <- ifelse(aux_1 > 1, aux_1 <- NA, aux_1 <- aux_1)
  aux_2 <- ifelse(aux_2 > 1, aux_2 <- NA, aux_2 <- aux_2)  
  #Diff_HG[[j]] <- abs(aux_1[!is.na(aux_1)])
  #Diff_HG[[j]] <- (aux_1[!is.na(aux_1)])
  Diff_HGT[[i]] <- abs(aux_1[!is.na(aux_1)])
  Diff_HGT2[[i]] <- abs(aux_2[!is.na(aux_2)])
  Perst_T1[[i]] <- Bar[!is.na(Bar)]#[1:length()]
  Perst_T2[[i]] <- Bar1[!is.na(Bar1)]#[1:20]
  Perst_T3[[i]] <- Bar2[!is.na(Bar2)]#[1:20]
  # Estimation of the persistence bars
}
# Computation of the entropy
library("infotheo")
#library("DescTools")
#library("FNN")
aux <- 0
#library("philentropy")
library("entropy")
Geom_Phi <- c()
Mut_inf <- c()
Geom_Phi1 <- c()
Mut_inf1 <- c()
for(i in 1:8){
  #aux <- ifelse(Diff_HGT[[i]] > 1, aux <- NA, aux<- Diff_HGT[[i]])
  #aux2 <- ifelse(Diff_HGT2[[i]] >= 1, aux2 <- NA, aux2<- Diff_HGT2[[i]])
  P <-abs(Perst_T1[[i]])[1:20]   
  #P <- aux[!is.na(aux)][1:4]
  Q <- abs(Perst_T2[[i]])[1:20]
  Q1 <- abs(Perst_T3[[i]])[1:20]
  #Q <- aux2[!is.na(aux2)][1:4]
  #x <- rbind(P,Q)
  join_PQ <- cbind(P,Q)
  join_PQ1 <- cbind(Q,Q1)
  # Compute Kullback-Leiber Entropy
  Geom_Phi[i] <- KL.plugin(P, Q)
  Geom_Phi1[i] <- KL.plugin(Q, Q1)
  # Compute the mutual information - in order to get the limit of the entropies
  Mut_inf[i] <- entropy(P,Q)
  #Mut_inf[i] <- mi.plugin(join_PQ)
  Mut_inf1[i] <- entropy(Q,Q1)
  #Mut_inf1[i] <-  mi.plugin(join_PQ1)
}
NL <- 0.005
Diff_HG <- list()
Diff_HGT <- list()
for(i in 1:NumberPatients){
  for(j in 1:NumberPatients){
    if(i == j){
      auxp_i <- 1
      auxp_j <- 2
    } 
    else{
      auxp_i <- 2
      auxp_j <- 2
    }
    #Bar <- 0
    #Bar1 <- 0
    Bar <- DiagPH_T[[i]][[auxp_i]]$diagram[,2] - DiagPH_T[[i]][[auxp_i]]$diagram[,3]
    #Bar <- Bar/max(Bar)
    Bar <- unique(ifelse(abs(Bar) >= NL,Bar <- Bar,Bar<- NA))
    Bar1 <- DiagPH_T[[j]][[auxp_j]]$diagram[,2] - DiagPH_T[[j]][[auxp_j]]$diagram[,3]
    #Bar1 <- Bar1/max(Bar1)
    Bar1 <- unique(ifelse(abs(Bar1) >= NL,Bar1 <- Bar1,Bar1<- NA))
    #Different time periods same individual
    
    aux_1 <- (Bar[1:15] - Bar1[1:15]) 
    #Diff_HG[[j]] <- abs(aux_1[!is.na(aux_1)])
    Diff_HG[[j]] <- (aux_1[!is.na(aux_1)])
  }
  Diff_HGT[[i]] <- list(Diff_HG[[1]], Diff_HG[[2]], Diff_HG[[3]],Diff_HG[[4]],Diff_HG[[5]], Diff_HG[[6]], Diff_HG[[7]],Diff_HG[[8]])
}
library(entropy)
aux_x2 <- matrix()
for(k in 1: NumberPatients){
  
  p2 <- c(entropy(sqrt((Diff_HGT[[k]][[1]])**2)),entropy(sqrt((Diff_HGT[[k]][[2]])**2)),entropy(sqrt((Diff_HGT[[k]][[3]])**2)),entropy(sqrt((Diff_HGT[[k]][[4]])**2)),entropy(sqrt((Diff_HGT[[k]][[5]])**2)),entropy(sqrt((Diff_HGT[[k]][[6]])**2)),entropy(sqrt((Diff_HGT[[k]][[7]])**2)),entropy(sqrt((Diff_HGT[[k]][[8]])**2)))
  ifelse(k == 1, aux_x2 <- p2,aux_x2 <- cbind(aux_x2,p2))
}
aux_x3 <- matrix()
for(k in 1: NumberPatients){
  
  p3 <- c(entropy((Diff_HGT[[k]][[1]]**2)),entropy((Diff_HGT[[k]][[2]]**2)),entropy((Diff_HGT[[k]][[3]]**2)),entropy((Diff_HGT[[k]][[4]]**2)),entropy((Diff_HGT[[k]][[5]]**2)),entropy((Diff_HGT[[k]][[6]]**2)),entropy((Diff_HGT[[k]][[7]]**2)),entropy((Diff_HGT[[k]][[8]]**2)))
  ifelse(k == 1, aux_x3 <- p3,aux_x3 <- cbind(aux_x3,p3))
}
# Final visualization of the Phi_S diagram
library(ggplot2)
S_M <- c()
refP <- mean(Mut_inf)
for(i in 1:8){
  S_M[i] <- mean(aux_x3[,i][!is.na(aux_x3[,i])]) #mean(aux_x3[,i])
}
Pat_Nr <- c(1,2,3,4,5,6,7,8)
aux <- cbind(S_M,Geom_Phi,Pat_Nr)
aux1 <- cbind(S_M,Geom_Phi1,Pat_Nr)
Obs <- as.data.frame(rbind(aux,aux1))#as.data.frame(cbind(S_M,Geom_Phi,Pat_Nr))
p <- ggplot(Obs,aes(S_M,Geom_Phi))+  geom_point(aes(size = 1,col=Pat_Nr))
p + xlim(c(0,2.0))+ylim(c(0, 3.0)) + geom_abline(intercept = refP, slope = 0)+ labs(x=(expression(S(Gamma))), y=(expression(Phi[GP]))) #+geom_abline(intercept = 0, slope = 1)


