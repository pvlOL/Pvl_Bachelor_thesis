### Applying causal influence tools to lake community data
### By Paula van Laak
### Oldenburg, 19.08.2019

# All EDM analyses that are presented and discribed in the thesis can be found in this script. Preliminary analyses 
# are not part of this script. 
# The order of the script is meant to fit the order of the thesis and should therefore follow a distinct structure:

# set working directory
setwd()

# Load data
lakes <- read.csv("Iowa_lakes_data_corrected_2014.csv", header=T)
dist_geo <- read.csv("lake_dist.csv", header=T)
phyto <- read.csv("phyto_ID.csv")
phyto_fg <- read.csv("phyto_fg_ID.csv")
dist_env <- read.csv("lake_dist_env.csv", sep=",", dec=".")
zooplankton_data <- read.csv("Iowa zoop data with diversity FINAL.csv", sep=",", dec=".")

# load additional packages
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(rEDM)
library(stringr)
library(reshape2)
library(forcats)


#-----------------EDM: Biomasse zwischen den Seen pro funkt. Gruppe-----------------------------------
# EDM for prediction of biomass.

#-1.1-------------Data preparation--------------------------------------------------------------------

# sort phytoplankton data alphabetically
phyto_fg <- data.frame(phyto_fg[order(phyto_fg$lake_name, phyto_fg$sample_ID),])

# Get names of functional groups and add total biomass for PP and ZP
fg <- names(phyto_fg)[-c(1:2, 11:12)]
fg[c(8,9)] <- c("PPbiom", "ZPbiom")

# Extract environmental data and adjust names
env_dat <- lakes %>%
  rename(., "lake_name" = Lake_name, "sample_year" = Year, "sampling_round" = Round) 

# Build data frame with combines data
dat <- left_join(phyto_fg, env_dat) 
dat <- right_join(dat, rename(merge(unique(dat$lake_name), unique(dat$sample_ID)),"lake_name" = x, "sample_ID" = y))

# Adjust zooplankton data and compute total zooplankton biomass per sampling round
dat_zoo <- zooplankton_data %>%
  mutate("Copepoda_Nauplii.bio" = Copepoda_Nauplii.bio + Nauplii_Nauplii.bio) %>%
  select(c(3,6:7), ends_with(".bio", ignore.case = F), -Nauplii_Nauplii.bio) %>%
  gather(-c(1:3), key = "Species", value = "ZP_biom") %>%
  mutate(func_group = unlist(str_split(Species, "_",n = 2))[c(T,F)]) %>%
  # drei funktionelle Gruppen: Cladocera, Copepoda, Rotifera
  filter(func_group == c("Cladocera", "Copepoda", "Rotifera")) %>%
  group_by(lake_name, year_sampling, season_round) %>%
  summarise(ZPbiom = sum(ZP_biom)) %>%
  rename(., sample_year = year_sampling, sampling_round = season_round) %>%
  ungroup()

# Add ZP data to data frame and arrange alphabetically
dat <- left_join(dat, dat_zoo, by = c("lake_name", "sample_year", "sampling_round"))
dat <- arrange(dat, lake_name, sample_ID)

#----------------Functions----------------------------------------------------------------------------------
# Function for Normalization
normalize <- function(x, ...) {
  (x - mean(x, ..., na.rm = T))/sd(x, ..., na.rm = T)
}


### 1. Find dynamics with nonlinear behaviour----------------------------------------------------------------

# Choose relevant variables; no additional environmental parameters needed
vars <- fg

# Create list split by lakes and normalise data
composite_ts <- dat[, vars]
dat_by_lake <- split(composite_ts, dat$lake_name)
normalized_data <- lapply(dat_by_lake, function(df) sapply(df, normalize))
composite_ts <- cbind(sample_ID = dat$sample_ID, data.frame(do.call(rbind,normalized_data)))

# Get segments from data set for each lake
segments_end <- cumsum(sapply(dat_by_lake, NROW))
segments_begin <- c(1, segments_end[-length(segments_end)] + 1)
segments <- cbind(segments_begin, segments_end)
row.names(segments) <- names(segments_end)

# Create data frame for final results
dat.out <- data.frame(merge(unique(dat$lake_name), fg[1:9]), "rho" = NA)

# for-loop for application of S-mapping for each lake separately
for (j in 1:length(unique(dat$lake_name))) {  
  # Get indices from data set for the corresponding lake
  seg_A <- segments[unique(dat$lake_name)[j] == row.names(segments),]
  
  # Apply simplex projections to get best embedding dimension 
  simplex_out_AfromA <- lapply(names(composite_ts)[c(2:10)], function(var) {
    simplex(composite_ts[, c( var)], E = 1:10, lib = seg_A, pred = seg_A, silent = TRUE)
  })
  names(simplex_out_AfromA) <- names(composite_ts)[c(2:10)]
  
  #nur einmalig fuer Speichern des besten rho fuer dei Vorhersage innerhalb eines jeden Sees verwendet
  #dat.out[dat.out$x == unique(dat$lake_name)[j], 3] <- unlist(lapply(simplex_out_AfromA, function(var) {max(var$rho)}))
  
  # Save best E
  best_E_A <- lapply(simplex_out_AfromA, function(df) {if (is.finite(max(df$rho))) {df$E[which.max(df$rho)]}
    else c(0)})
  
  # Apply S-mapping to evaluate nonlinear dynamics
  smap_out_AfromA <- lapply(names(composite_ts)[c(2:10)], function(var) {
    s_map(composite_ts[, c(var)], E = unlist(best_E_A[var]), lib = seg_A,
          pred = seg_A, silent = TRUE)})
  names(smap_out_AfromA) <- names(simplex_out_AfromA)
  # Save highest prediction skill (rho)
  rho <- lapply(smap_out_AfromA, function(var) {var$rho[is.na(var$rho)] = 0; return(var$rho*(var$rho > 0))} )
  dat.out[dat.out$x == unique(dat$lake_name)[j], 3] <- unlist(lapply(rho, max))

  
  # Plot prediction skill rho as a function of theta for every functional group and save plot to later evaluate nonlinearity
  png(sprintf("Plots/NonLin_/%s.png", unique(dat$lake_name)[j]), width = 800, height = 700)
  par(mfrow= c(2,5))
  p <- lapply(names(smap_out_AfromA), function(var){
    plot(smap_out_AfromA[[var]]$theta, rho[[var]], main = unique(dat$lake_name)[j], ylab = var)
  })
  dev.off()
}

# Write .csv-file with results of prediction skill for each lake predicting itself
#write.csv(dat.out, "InLake_simplex.csv", row.names = F)

# Nonlinearity was analysed manually by looking at the plots. Increasing prediction skill with increasing theta 
# is an indicator of nonlinear dynamics beacause it shows that a nonlinear model fits better than a linear model. 
nonLin <- read.csv("NonLin_sampleID.csv", sep=",", dec=".")
# Add prediciton skill to data of nonlinear dynamics
InLake_simp <- read.csv("InLake_simplex_Rho_sampleID.csv", sep=",", dec=".")
InLake_pred <- full_join(rename(InLake_simp, Lake.name = x, func_group = y), 
                         gather(nonLin[,c(2:ncol(nonLin))], -Lake.name, key = "func_group", value = "nonLin"))

### 2-----Convergent-Cross-Mapping---------------------------------------------------------------------------

# 2.1-----CCM of environmental parameters and biomass of functional groups-----------------------------------

# Create list of all lakes
  # Choose which environmental parameters should be considered
vars <- c( fg, "Secchi", "TP","TN", "Temp", "Chl_a", "TSS", "VSS", "ISS")
composite_ts <- dat[, vars]
dat_by_lake <- split(composite_ts, dat$lake_name)
  # normalise data by lake and variable
normalized_data <- lapply(dat_by_lake, function(df) sapply(df, normalize))
  # Create data frame with normalised data
composite_ts <- cbind(sample_ID = dat$sample_ID, data.frame(do.call(rbind,normalized_data)))

# Get segments of lakes for attractor reconstruction; Indices of each lake are extracted
segments_end <- cumsum(sapply(dat_by_lake, NROW))
segments_begin <- c(1, segments_end[-length(segments_end)] + 1)
segments <- cbind(segments_begin, segments_end)
row.names(segments) <- names(segments_end)

# Create data frame for CCM results
env_vars <- names(composite_ts)[11:NCOL(composite_ts)]
dat.out.final <-left_join( 
  rename(merge(
    rename(merge(
      rename(merge(unique(dat$lake_name), fg), "lake_name" = x, "func_group" = y), env_vars),
      "env_var" = y),
    c("PPxmapENV", "ENVxmapPP")), "direction" = y), rename(InLake_pred[,-3], lake_name = Lake.name))

# Set library size and interval; Add a column for each library value
lib_sizes <- seq(10,40,by = 2)
dat.out.final[,as.character(lib_sizes)] <- NA

# for-loop for CCM of all lake-environmental parameter pair
for (i in 1:9) {  
  fg_z <- i+1           # Set index of functional group
  fg_n <- fg[i]         # Set name of functional group
  # Create a temporary data frame for results of this functional group; 
  # Only choose lakes with nonlinear dynamics (1-> strong inidcation of nonlinearity, 2-> weak indication of nonlinearity,
  #                                            0-> no nonlinearity)
  dat.out <- select(filter(dat.out.final, func_group == fg[i], (nonLin== 1 | nonLin == 2)), c(1,3:NCOL(dat.out.final)))

  for (j in 1:(dim(dat.out)[1]*0.5)) {  # Select each pair one by one (Only half of all pairs beacause each pair is tested
                                        # in both directions with each in a separate column)
    # select environmental parameter
    env <- as.character(dat.out$env_var[j])
    # Get index of corresponding lake within the dataset
    seg <- segments[dat.out$lake_name[j] == row.names(segments),]    
    
    # Apply simplex projection for best embedding dimension
    simplex_out <- lapply(c(fg_n, env), function(var) {
        simplex(composite_ts[,c("sample_ID", var)], E = 1:10, lib = seg, pred = seg, silent = TRUE)
    })
    names(simplex_out) <- c(fg_n, env)
    
    # Save best embedding dimension for both variables
    best_E_PP  <- simplex_out[[fg_n]][["E"]][which.max(simplex_out[[fg_n]][["rho"]])]
    best_E_env <- simplex_out[[env]][["E"]][which.max(simplex_out[[env]][["rho"]])]
    
    # Perform CCM in both directions with the corresponding best E
    ccm_PPfromEnv <- ccm(composite_ts[, c("sample_ID", fg_n, env)], E = best_E_env, lib = seg, pred = seg, 
                         lib_column = env, target_column = fg_n,  first_column_time = T, silent = T,lib_sizes = lib_sizes)
    ccm_EnvfromPP <- ccm(composite_ts[, c("sample_ID", fg_n, env)], E = best_E_PP, lib = seg, pred = seg,
                         lib_column = fg_n, target_column = env,  first_column_time = T, silent = T,lib_sizes = lib_sizes)
    # Get mean values of results
    ccm_PPfromEnv_mn <- ccm_means(ccm_PPfromEnv)
    ccm_EnvfromPP_mn <- ccm_means(ccm_EnvfromPP)    
    
    # Save results in temporary data frame
    dat.out[j,c(5:NCOL(dat.out))] <- ccm_EnvfromPP_mn$rho
    dat.out[j+(0.5*NROW(dat.out)), c(5:NCOL(dat.out))] <- ccm_PPfromEnv_mn$rho
  }
  # Transfer data from temporary data frame to final data frame
   dat.out.final[dat.out.final$func_group == fg_n &(dat.out.final$nonLin== 1 | dat.out.final$nonLin == 2), 
                 c(6:ncol(dat.out.final))] <- dat.out[,c(5:NCOL(dat.out))]
}

# Write .csv-file with final results
#write.csv(dat.out.final, "CCM_PP_ENV_sampleID.csv", row.names = F)

# It is essential that prediction skill increases with library size in order to identify causal relationships 
# among phytoplankton biomass and environmental parameters.
# It is expected that environmental parameters show a greater impact on pyhtoplankton biomass than vice versa. 
# Therefore, prediction skill should yield higher values for PPxmapENV.

# 2.2----------CCM of Phytoplankton and Zooplankton biomass wihtin lakes---------------------------------

# Get information on nonlinear dynamics
InLake_pred <- nonLin[,c("Lake.name", "PPbiom", "ZPbiom")]

# Create subset of lakes with nonlinear dynamics
samp_nonlin <- InLake_pred$Lake.name[InLake_pred$PPbiom %in% c(1,2) & InLake_pred$ZPbiom %in% c(1,2)]
Plank_dat <- dat[dat$lake_name %in% samp_nonlin, ]

# Select relevant variables
vars <- c("PPbiom", "ZPbiom")

# Create list split by lakes and normalise data
composite_ts <- Plank_dat[, vars]
dat_by_lake <- split(composite_ts, as.character(Plank_dat$lake_name))
normalized_data <- lapply(dat_by_lake, function(df) sapply(df, normalize))
composite_ts <- cbind(sample_ID = Plank_dat$sample_ID, data.frame(do.call(rbind,normalized_data)))

# Get indices of lakes from dataset
segments_end <- cumsum(sapply(dat_by_lake, NROW))
segments_begin <- c(1, segments_end[-length(segments_end)] + 1)
segments <- cbind(segments_begin, segments_end)
row.names(segments) <- names(segments_end)

# Create data frame for final results
dat.out <- rename(merge(unique(Plank_dat$lake_name), c("PPxmapZP", "ZPxmapPP")), "lake_name" = x, "direction" = y)

# Set library size and intervals; Add columns for each library value to data frame for final results
lib_sizes <- seq(10,40,by = 2)
dat.out[,as.character(lib_sizes)] <- NA

# for-loop for application of CCM 
for (j in 1:(dim(dat.out)[1]*0.5)) {  
  # Select corresponding data for this variable pair from dataset
  seg <- segments[dat.out$lake_name[j] == row.names(segments),]    
  
  # Apply simplex projection to get best embedding dimension
  simplex_out <- lapply(c("ZPbiom", "PPbiom"), function(var) {
    simplex(composite_ts[,c("sample_ID", var)], E = 1:10, lib = seg, pred = seg, silent = TRUE)
  })
  names(simplex_out) <- c("ZPbiom", "PPbiom")
  # Save best E for both variables
  best_E_PP <- simplex_out[["PPbiom"]][["E"]][which.max(simplex_out[["PPbiom"]][["rho"]])]
  best_E_ZP <- simplex_out[["ZPbiom"]][["E"]][which.max(simplex_out[["ZPbiom"]][["rho"]])]
  
  # Apply CCM with best E in both directions
  ccm_PPfromZP <- ccm(composite_ts[, c("sample_ID", "ZPbiom","PPbiom")], E = best_E_ZP, lib = seg, pred = seg,
                      lib_column = "ZPbiom", target_column = "PPbiom",  first_column_time = T, silent = T,lib_sizes = lib_sizes)
  ccm_ZPfromPP <- ccm(composite_ts[, c("sample_ID", "ZPbiom","PPbiom")], E = best_E_PP, lib = seg, pred = seg,
                      lib_column = "PPbiom", target_column = "ZPbiom",  first_column_time = T, silent = T,lib_sizes = lib_sizes)
  # Get means of all results
  ccm_PPfromZP_mn <- ccm_means(ccm_PPfromZP)
  ccm_ZPfromPP_mn <- ccm_means(ccm_ZPfromPP)    
  
  # Save results in final data frame
  dat.out[j,c(3:NCOL(dat.out))] <- ccm_ZPfromPP_mn$rho
  dat.out[j+(0.5*NROW(dat.out)), c(3:NCOL(dat.out))] <- ccm_PPfromZP_mn$rho
}

# Write .csv-file with all results 
#write.csv(dat.out, "CCM_PP_ZP_sampleID.csv", row.names = F)


#-3.1.-------------Univariate embedding---------------------------------------------------------------

# Choose relevant variables -> in the univariate model no environmental parameters are added
vars <- fg

# Create list of relevant data split by lakes and normalise data
composite_ts <- dat[, vars]
dat_by_lake <- split(composite_ts, dat$lake_name)
normalized_data <- lapply(dat_by_lake, function(df) sapply(df, normalize))
composite_ts <- cbind(sample_ID = dat$sample_ID, data.frame(do.call(rbind,normalized_data)))

# Get segments of each lake from data set
segments_end <- cumsum(sapply(dat_by_lake, NROW))
segments_begin <- c(1, segments_end[-length(segments_end)] + 1)
segments <- cbind(segments_begin, segments_end)
row.names(segments) <- names(segments_end)

# Create dat aframe for final results of simplex projection
dat.out <- full_join(dist_env, dist_geo); names(dat.out)[c(3,4)] <- c("env_dist", "geo_dist")
dat.out <- data.frame(dat.out, "best_E_B" = NA, "best_E_A" = NA, "rho_BfromA" = NA, "rho_AfromB" = NA, 
                      "rmse_BfromA" = NA, "rmse_AfromB" = NA)
dat.out.final <- rename(merge(fg, dat.out), "func_group" = x)

# Select subset of lakes with nonlinear dynamics
dat.out.final <-  
  select(filter(merge(merge(dat.out.final, rename(gather(nonLin[,-1], fg , key = "func_group", value = "nonLin2"), "lake2" = Lake.name)), 
                      rename(gather(nonLin[,-1], fg , key = "func_group", value = "nonLin1"), "lake1" = Lake.name)),
                nonLin1 != 0, nonLin2 != 0),
          -nonLin1, -nonLin2)

# for-loop for Application of simplex projection on all lake pairs
for (i in 1:length(fg)) { 
  fg_z <- i+1           # Select index of one functional group
  fg_n <- fg[i]         # Select name of functional group 
  env_add <- c()        # no additional environmental parameter
  
  # Create temporary data frame for results of just this functional group
  dat.out <- select(filter(dat.out.final, func_group == fg[i]), c(2:ncol(dat.out.final)))
  
  for (j in 1:dim(dat.out)[1]) {
    # Select Indices of the corresponding lakes from data set
    seg_A <- segments[dat.out$lake1[j] == row.names(segments),]    
    seg_B <- segments[dat.out$lake2[j] == row.names(segments),]    
    
    # Make predictions for each lake itself
    simplex_out_AfromA <- simplex(composite_ts[,c("sample_ID", fg_n)], E = 1:10, lib = seg_A, pred = seg_A, silent = TRUE)
    simplex_out_BfromB <- simplex(composite_ts[,c("sample_ID", fg_n)], E = 1:10, lib = seg_B, pred = seg_B, silent = TRUE)
    # get best embedding dimensions
    best_E_A <- simplex_out_AfromA$E[which.max(abs(simplex_out_AfromA$rho))]
    best_E_B <- simplex_out_BfromB$E[which.max(abs(simplex_out_BfromB$rho))]
    
    # Make predictions in both directions for both lakes respectively
    simplex_out_BfromA <- lapply(names(composite_ts)[c(fg_z,env_add)], function(var) {
      simplex(composite_ts[, c("sample_ID", var)], E = best_E_A, lib = seg_A, pred = seg_B)
    })
    simplex_out_AfromB <- lapply(names(composite_ts)[c(fg_z,env_add)], function(var) {
      simplex(composite_ts[, c("sample_ID", var)], E = best_E_B, lib = seg_B, pred = seg_A)
    })
    names(simplex_out_BfromA) <- names(composite_ts)[c(fg_z,env_add)]
    names(simplex_out_AfromB) <- names(composite_ts)[c(fg_z,env_add)]
    
    # Save the best rho      
    rho_BfromA <- sapply(simplex_out_BfromA, function(df) {max(df$rho)})
    rho_AfromB <- sapply(simplex_out_AfromB, function(df) {max(df$rho)})
    # Save corresponding error (root mean squared error)
    rmse_BfromA <- sapply(simplex_out_BfromA, function(df) {df$rmse[which.max(df$rho)]})
    rmse_AfromB <- sapply(simplex_out_AfromB, function(df) {df$rmse[which.max(df$rho)]})    
    
    # Set non-numeric values to NA
    if (is.numeric(best_E_A) == FALSE) {best_E_A <- NA}
    if (is.numeric(best_E_B) == FALSE)  {best_E_B <- NA}
    if (is.numeric(rho_AfromB) == FALSE) {rho_AfromB <- NA}
    if (is.numeric(rho_BfromA) == FALSE) {rho_BfromA <- NA}
    if (is.numeric(rmse_AfromB) == FALSE) {rmse_AfromB <- NA}
    if (is.numeric(rmse_BfromA) == FALSE) {rmse_BfromA <- NA}
    
    # Save results to final data frame
    dat.out[j,c(5:10)] <- c(best_E_B, best_E_A, rho_BfromA, rho_AfromB, rmse_BfromA, rmse_AfromB)
  }
  # alle Daten fuer die funktionelle Gruppe gesammelt in den finalen Datensatz uebertragen
  dat.out.final[dat.out.final$func_group == fg[i],c(2:ncol(dat.out.final))] <- dat.out
}

# Write .csv-file with results
#write.csv(dat.out.final, "simplex_allcomb_sampleID.csv", row.names = F)


#-3.2.-------------Multivariate embedding---------------------------------------------------------------

# Select relevant additional environmental parameters
vars <- c("sample_ID", fg, "Secchi", "TP","TN", "Temp", "Chl_a", "TSS")

# Create a list of all relevant data split by lakes and normalise data
composite_ts <- dat[, vars]
dat_by_lake <- split(composite_ts, dat$lake_name)
normalized_data <- lapply(dat_by_lake, function(df) sapply(df, normalize))
composite_ts <- cbind(sample_ID = dat$sample_ID, data.frame(do.call(rbind,normalized_data))[,-1])

# Get segments of each lake from dataset
segments_end <- cumsum(sapply(dat_by_lake, NROW))
segments_begin <- c(1, segments_end[-length(segments_end)] + 1)
segments <- cbind(segments_begin, segments_end)
row.names(segments) <- names(segments_end)

# Create data frame for final results
dat.out <- full_join(dist_env, dist_geo); names(dat.out)[c(3,4)] <- c("env_dist", "geo_dist")
dat.out <- data.frame(dat.out, "rho_BfromA" = NA, "rho_AfromB" = NA, "p_val_BfromA" = NA, "p_val_AfromB" = NA,
                      "rmse_BfromA" = NA, "rmse_AfromB" = NA, "rho_BfromA_multi" = NA, "rho_AfromB_multi" = NA, "p_val_BfromA_multi" = NA, "p_val_AfromB_multi" = NA,
                      "rmse_BfromA_multi" = NA, "rmse_AfromB_multi" = NA)
dat.out.final <- rename(merge(fg, dat.out), "func_group" = x)

# Select subset of lakes with nonlinear dynamics
dat.out.final <-  
  select(filter(merge(merge(dat.out.final, rename(gather(nonLin[,-1], fg , key = "func_group", value = "nonLin2"), "lake2" = Lake.name)), 
                      rename(gather(nonLin[,-1], fg , key = "func_group", value = "nonLin1"), "lake1" = Lake.name)),
                nonLin1 != 0, nonLin2 != 0),
         -nonLin1, -nonLin2)


# Select columns for prediction: It can be chosen between two or three time lags of phytoplankton biomass for attractor reconstruction
simple_columns <- c("fg_tm", "fg_tm1", "fg_tm2")
#simple_columns <- c("fg_tm", "fg_tm1", "fg_tm2", "fg_tm3")

# Option: "pred1" --> includes three time lage of one phytoplankton group, TP, TN, Temperature and Secchi depth
#pred_columns <- c("fg_tm", "fg_tm1", "fg_tm2", "fg_tm3", "TP_tm", "TP_tm1", "TP_tm2", "TN_tm3","TN_tm", "TN_tm1", "TN_tm2", "TN_tm3", "Temp_tm",
#                  "Temp_tm1", "Temp_tm2", "Temp_tm3", "St_tm", "St_tm1", "St_tm2", "St_tm3")

# Option: "pred2" --> includes three time lags of phytoplankton group and Secchi depth
#pred_columns <- c("fg_tm", "fg_tm1", "fg_tm2", "fg_tm3", "St_tm", "St_tm1", "St_tm2", "St_tm3")
# Option: "pred2_2" --> includes two time lags of phytoplankton group and Secchi depth
#pred_columns <- c("fg_tm", "fg_tm1", "fg_tm2", "St_tm", "St_tm1", "St_tm2")

# Option: "pred3" --> includes three time lags of phytoplankton group and TN
#pred_columns <- c("fg_tm", "fg_tm1", "fg_tm2", "fg_tm3", "TN_tm3","TN_tm", "TN_tm1", "TN_tm2", "TN_tm3")

# Variante: "pred 4" --> includes two time lags of phytoplankton group, TP and temperature
pred_columns <- c("fg_tm", "fg_tm1", "fg_tm2", "TP_tm", "TP_tm1", "TP_tm2", "Temp_tm", "Temp_tm1", "Temp_tm2")


for (i in 1:length(fg)) {
  fg_z <- i+1          
  fg_n <- fg[fg_z-1]
  
  # Manuall ycreate time lags of possible environmental parameters and the time series of the functional group
  block_data <- do.call(rbind, lapply(dat_by_lake, function(df) {
    n <- NROW(df)
    temp <- data.frame(Time = df$sample_ID)
    temp$fg_tm <- df[[fg_n]]
    temp$fg_tm1 <- c(NA, temp$fg_tm[-n])
    temp$fg_tm2 <- c(NA, temp$fg_tm1[-n])
    temp$fg_tm3 <- c(NA, temp$fg_tm2[-n])
    
    temp$TP_tm <- df$TP
    temp$TP_tm1 <- c(NA, temp$TP_tm[-n])
    temp$TP_tm2 <- c(NA, temp$TP_tm1[-n])
    temp$TP_tm3 <- c(NA, temp$TP_tm2[-n])
    
    temp$TN_tm <- df$TN
    temp$TN_tm1 <- c(NA, temp$TN_tm[-n])
    temp$TN_tm2 <- c(NA, temp$TN_tm1[-n])
    temp$TN_tm3 <- c(NA, temp$TN_tm2[-n])
    
    temp$Temp_tm <- df$Temp
    temp$Temp_tm1 <- c(NA, temp$Temp_tm[-n])
    temp$Temp_tm2 <- c(NA, temp$Temp_tm1[-n])
    temp$Temp_tm3 <- c(NA, temp$Temp_tm2[-n])
    
    temp$St_tm <- df$Secchi
    temp$St_tm1 <- c(NA, temp$St_tm[-n])
    temp$St_tm2 <- c(NA, temp$St_tm1[-n])
    temp$St_tm3 <- c(NA, temp$St_tm2[-n])
    
    return(temp)
  }))
  
  # Apply block_lnlp function
  
  # Create temporary data frame for results
  dat.out <- select(filter(dat.out.final, func_group == fg[i]), c(2:ncol(dat.out.final)))
  
  for (j in 1:nrow(dat.out)) {
      # Select corresponding segemnts from data set 
      segA <- segments[(dat.out[j,1] == row.names(segments)),]
      segB <- segments[(dat.out[j,2] == row.names(segments)),]
     
      # Make predictions from lake A on lake B applying the uni- and the multivariate approach
      simple_output <- block_lnlp(block_data, lib = segA, pred = segB,
                                  columns = simple_columns, target_column = 1, stats_only = FALSE, first_column_time = TRUE)
      
      multi_output <- block_lnlp(block_data, lib = segA, pred = segB, 
                                 columns = pred_columns, target_column = 1, stats_only = FALSE, first_column_time = TRUE)
      
      # Save results
      dat.out[j, c(5,7,9)] <- simple_output[,c("rho", "p_val", "rmse")]
      dat.out[j, c(11,13,15)] <- multi_output[,c("rho", "p_val", "rmse")]

      # Make predictions from lake B on lake A applying the uni- and the multivariate approach
      simple_output <- block_lnlp(block_data, lib = segB, pred = segA, 
                                  columns = simple_columns, target_column = 1, stats_only = FALSE, first_column_time = TRUE)
      
      multi_output <- block_lnlp(block_data, lib = segB, pred = segA, 
                                 columns = pred_columns, target_column = 1, stats_only = FALSE, first_column_time = TRUE)
      # Save results
      dat.out[j, c(6,8,10)] <- simple_output[,c("rho", "p_val", "rmse")]
      dat.out[j, c(12,14,16)] <- multi_output[,c("rho", "p_val", "rmse")]
  }
  # Save temporary data frame to final data frame for all results
  dat.out.final[dat.out.final$func_group==fg[i],c(2:NCOL(dat.out.final))] <- dat.out
}

# Write .csv-file with all results
#write.csv(dat.out.final, "EDM_multi_pred4_2lags.csv", row.names = F)

