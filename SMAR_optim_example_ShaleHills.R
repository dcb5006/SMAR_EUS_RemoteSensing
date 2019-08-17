rm(list = ls())
library(grDevices); library(xts); library(plyr); library(hydroGOF); 
library(MASS); library(coda); library(DEoptim); library(sirad)
library(FME); library(fields); library(rgdal); library(imputeR)


SH_data <- read.csv('Moisture_data.csv')


Texture_data <- read.csv("SH_RealtimeSite_Data.csv", header = TRUE)



###############
# Model Data
###############


SM_dates <- SH_data[,1]
SM_dates <- ISOdate(substr(SM_dates,1,4),substr(SM_dates,6,7),substr(SM_dates,9,10))

# Get Time axes labels


Years <- unique(format(SM_dates, format = "%Y"))
Years <- ISOdate(Years,rep(1,length(Years)),rep(1,length(Years)))

Months <- unique(format(SM_dates, format = "%Y %m"))
Months <- ISOdate(substr(Months,1,4),substr(Months,6,7),rep(1,length(Months)))
Months_labels <- Months[seq(1,length(Months),by = 3)]



S_5cm <- as.numeric(SH_data[,16]) ## Optimize to satellite data
S_Root <- as.numeric(SH_data[,3])


## Deal with NA gaps at front/back of timeseries ##


starting_point <- min(which(!is.na(S_Root) & !is.na(S_5cm)))
ending_point <- max(which(!is.na(S_Root) & !is.na(S_5cm)))

S_5cm <- S_5cm[starting_point:ending_point]
S_Root <- S_Root[starting_point:ending_point]
SM_dates <- SM_dates[starting_point:ending_point]




# Just Focus on 2009 and 2010 for optimization

S_5cm <- S_5cm[which(format(SM_dates, "%Y") == "2009" | format(SM_dates, "%Y") == "2010")]
S_Root <- S_Root[which(format(SM_dates, "%Y") == "2009" | format(SM_dates, "%Y") == "2010")]
SM_dates <- SM_dates[which(format(SM_dates, "%Y") == "2009" | format(SM_dates, "%Y") == "2010")]


S_5cm <- S_5cm/max(S_5cm,na.rm = TRUE)


# Clean root zone moisture data

S_Root[which(S_Root < 0.2)] <- NA
S_Root[which(S_Root > 1)] <- 1


# Remove surface moisture NAs

S_Root <- S_Root[which(!is.na(S_5cm))]
SM_dates <- SM_dates[which(!is.na(S_5cm))]
S_5cm <- S_5cm[which(!is.na(S_5cm))]




########################
## SMAR Model (original)
########################

SMAR_predict_v1 <- function(p,x,RZSM1){ 

  Sfc_sl <- p[1]
  a <- p[2]
  b <- p[3]
  Sw_rz <- p[4]
  
  n <- length(x)
  
  #### METHODS ###### 
  
  
  RZSM <- c()
  Sy <- c()
  

  RZSM[1] <- RZSM1  
  
  
  for(i in 2:n){
    Sy[i] <- x[i] - Sfc_sl
    Sy[i] <- ifelse(Sy[i] < 0, 0, Sy[i])
    Sy[i] <- ifelse(Sy[i] > 1, 1, Sy[i])
      
    
    RZSM[i] <- Sw_rz + (RZSM[i - 1] - Sw_rz) * exp(-a) + ((1 - Sw_rz) * b * Sy[i])
    
      
  }
  
  for(p in 1:length(RZSM)){
    
    RZSM[p] <- ifelse(RZSM[p] > 1 & is.na(RZSM[p]) == FALSE, 1, RZSM[p])
  }
  
  return(RZSM)
}




#######################################
# Pre-Calibration (Levenberg-Marquardt)
#######################################



Obs <- data.frame(x = S_5cm, y = S_Root)


Residuals <- function(p,RZSM1){ 
  
  
  Sfc_sl <- p[1]
  a <- p[2]
  b <- p[3]
  Sw_rz <- p[4]
  
  n <- length(Obs$x)
  
  #### METHODS ###### 
  
  
  RZSM <- c()
  Sy <- c()
  
  
  ## CHANGE 4 HERE ##
  
  RZSM[1] <- RZSM1 
  
  
  for(i in 2:n){
    Sy[i] <- Obs$x[i] - Sfc_sl
    Sy[i] <- ifelse(Sy[i] < 0, 0, Sy[i])
    Sy[i] <- ifelse(Sy[i] > 1, 1, Sy[i])
    
    
    RZSM[i] <- Sw_rz + (RZSM[i - 1] - Sw_rz) * exp(-a) + ((1 - Sw_rz) * b * Sy[i])
    
    
  }
  
  
  for(P in 1:length(RZSM)){
    
    RZSM[P] <- ifelse(RZSM[P] > 1 & is.na(RZSM[P]) == FALSE, 1, RZSM[P])
  }
  
  ## CHANGE 5 HERE ##
  
  RES <- RZSM - Obs$y
  RES <- na.omit(RES)
  
  return(RES)
  
}



# Look at site texture for starting values

Texture_data[which(Texture_data[,1] == "Site15"),]


0.62   # FC start
((524.8/365/10)+2)/((1 - 0.34)*95*0.479) # a start
(0.501*5)/((1 - 0.34)*95*0.479) # b start
0.34 # WL start


P <- modFit(f = Residuals, p = c(0.62, 0.071, 0.083, 0.34),lower = c(0.0001,0.0001,0.0001,0.0001), upper = c(1,2,1,1),
            RZSM1 = 0.78)
sP <- summary(P); sP



Covar <- sP$cov.scaled 
s2prior <- sP$modVariance


set.seed(1234567)



MCMC_SMAR <- modMCMC(f = Residuals, p = P$par, niter = 18000, burninlength = 6000, ntrydr = 2, 
                     var0 = s2prior, wvar0 = 1, updatecov = 100, lower = c(0.0001,0.0001,0.0001,0.0001), 
                     upper = c(1,2,1,1), RZSM1 = 0.78)



MC_results <- as.mcmc(MCMC_SMAR$pars)




MCMC_SMAR2 <- modMCMC(f = Residuals, p = c(0.62, 0.071, 0.083, 0.34), niter = 32000, burninlength = 20000, ntrydr = 2, 
                     var0 = s2prior, wvar0 = 1, updatecov = 100, lower = c(0.0001,0.0001,0.0001,0.0001), 
                     upper = c(1,2,1,1), RZSM1 = 0.78)



MC_results2 <- as.mcmc(MCMC_SMAR2$pars)




MCMC_SMAR3 <- modMCMC(f = Residuals, p = 0.50*P$par, niter = 18000, burninlength = 6000, ntrydr = 2, 
                      var0 = s2prior, wvar0 = 1, updatecov = 100, lower = c(0.0000001,0.0001,0.0001,0.0001), 
                      upper = c(1,2,1,1), RZSM1 = 0.78)


MC_results3 <- as.mcmc(MCMC_SMAR3$pars)



# Gelman-Rubin test (convergence)

Sample1 <- as.data.frame(MCMC_SMAR$pars)
Sample2 <- as.data.frame(MCMC_SMAR2$pars)
Sample3 <- as.data.frame(MCMC_SMAR3$pars)


colnames(Sample1) <- c("FC","a","b","WL") 
colnames(Sample2) <- c("FC","a","b","WL")
colnames(Sample3) <- c("FC","a","b","WL")

mcmc_1 <- mcmc(Sample1)
mcmc_2 <- mcmc(Sample2)
mcmc_3 <- mcmc(Sample3)

mcmc_all <- mcmc.list(list(mcmc_1,mcmc_2,mcmc_3))

gelman.diag(mcmc_all)
Gelman_diag <- gelman.diag(mcmc_all)


# MCMC trace plots


png(paste("MCMC_Traceplots_","Site15",".png",sep = ""), width = 9, height = 6, 
    units = "in", pointsize = 12, res = 300)
par(mfrow = c(2, 2), mar = c(2,2,2,1) + 0.1, oma = c(4,4,0,0) + 0.1)

plot(Sample1[,1], ylim = c(0,1), type = "l", lwd = 1, col = "red", xlab = "Iteration", main = "Field Capacity", ylab = "")
lines(Sample2[,1], lwd = 1, col = "blue", xlab = "Iteration")
lines(Sample3[,1], lwd = 1, col = "green", xlab = "Iteration")
text(9000,0.76,paste("Gelman-Rubin = ",round(Gelman_diag$psrf[1,1],2)))

plot(Sample1[,2], ylim = c(0,2), type = "l", lwd = 1, col = "red", xlab = "Iteration", main = "Water Loss (a)", ylab = "")
lines(Sample2[,2], lwd = 1, col = "blue", xlab = "Iteration")
lines(Sample3[,2], lwd = 1, col = "green", xlab = "Iteration")
text(9000,0.67,paste("Gelman-Rubin = ",round(Gelman_diag$psrf[2,1],2)))

plot(Sample1[,3], ylim = c(0,1), type = "l", lwd = 1, col = "red", xlab = "Iteration", main = "Diffusivity (b)", ylab = "")
lines(Sample2[,3], lwd = 1, col = "blue", xlab = "Iteration")
lines(Sample3[,3], lwd = 1, col = "green", xlab = "Iteration")
text(9000,0.67,paste("Gelman-Rubin = ",round(Gelman_diag$psrf[3,1],2)))

plot(Sample1[,4], ylim = c(0,1), type = "l", lwd = 1, col = "red", xlab = "Iteration", main = "Wilting Level", ylab = "")
lines(Sample2[,4], lwd = 1, col = "blue", xlab = "Iteration")
lines(Sample3[,4], lwd = 1, col = "green", xlab = "Iteration")
text(9000,0.70,paste("Gelman-Rubin = ",round(Gelman_diag$psrf[4,1],2)))



title(ylab = "",
      xlab = "Iteration", outer = TRUE, line = 1.8, cex.lab = 2)

dev.off()




# Histograms of posteriors


png(paste("Posterior_Histograms_","Site15",".png",sep = ""), width = 6, height = 6, 
    units = "in", pointsize = 12, res = 300)
par(mfrow = c(2, 2), mar = c(2,2,2,1) + 0.1, oma = c(4,4,0,0) + 0.1)

FC_samples <- c(Sample1[,1],Sample2[,1],Sample3[,1])
hist(FC_samples,xlim = c(0,1), col = "white", main = "Field Capacity", prob = TRUE, border = "white")
polygon(c(0,0,1,1),c(0,3,3,0),border = NA, col = "light grey")
lines(density(FC_samples),col = "blue", lwd = 3)
lines(c(mean(FC_samples),mean(FC_samples)),c(0,1000),lty = 2, lwd = 1, col = "red")

a_samples <- c(Sample1[,2],Sample2[,2],Sample3[,2])
hist(a_samples,xlim = c(0,1), col = "white", border = "NA",prob = TRUE, main = "Water Loss (a)")
polygon(c(0,0,1,1),c(0,3,3,0),border = NA, col = "light grey")
lines(density(a_samples),col = "blue", lwd = 3)
lines(c(mean(a_samples),mean(a_samples)),c(0,1000),lty = 2, lwd = 1, col = "red")

b_samples <- c(Sample1[,3],Sample2[,3],Sample3[,3])
hist(b_samples,xlim = c(0,1), col = "white", border = NA, prob = TRUE, main = "Diffusion (b)")
polygon(c(0,0,1,1),c(0,3,3,0),border = NA, col = "light grey")
lines(density(b_samples),col = "blue", lwd = 3)
lines(c(mean(b_samples),mean(b_samples)),c(0,1000),lty = 2, lwd = 1, col = "red")

WL_samples <- c(Sample1[,4],Sample2[,4],Sample3[,4])
hist(WL_samples,xlim = c(0,1), col = "white", border = NA, prob = TRUE, main = "Wilting Level")
polygon(c(0,0,1,1),c(0,3,3,0),border = NA, col = "light grey")
lines(density(WL_samples),col = "blue", lwd = 3)
lines(c(mean(WL_samples),mean(WL_samples)),c(0,1000),lty = 2, lwd = 1, col = "red")


title(ylab = "Density",
      xlab = "Parameter", outer = TRUE, line = 1.8, cex.lab = 2)

dev.off()



SMAR_params_est <- c(mean(FC_samples), mean(a_samples),mean(b_samples),mean(WL_samples))




#####################
# Model Prediction
#####################



# Model fit diagnostics 

MCMC_fit <- SMAR_predict_v1(SMAR_params_est,x = S_5cm,RZSM1 = 0.78)

MCMC_R2 <- summary(lm(S_Root~MCMC_fit))$r.squared
MCMC_RMSE <- rmse(S_Root,MCMC_fit)
MCMC_RMSE_vol <- rmse(S_Root*0.479,MCMC_fit*0.479)  # Convert to volumetric with porosity at this site (0.479 cm3/cm3)
MCMC_ResSD <- sd(na.omit(MCMC_fit - S_Root))





