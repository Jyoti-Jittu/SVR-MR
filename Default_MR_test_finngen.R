set.seed(100)

# load Libraries

library(TwoSampleMR)
library(dplyr)
library(purrr)
library(MRPRESSO)
library(ieugwasr)
library(MRInstruments) 
library(readr)
library(ggplot2)
library(MendelianRandomization)

# Load IVO function
source("Methods/t-value.R")

# Exposures =====
exposures_list_tValue <- list.files(path = "Methods/Exposures/")
exposures_tValue <- list()

for (i in seq_along(exposures_list_tValue)) {
  exposures_tValue[[i]] <- read_csv(paste0("Methods/Exposures/",
                                           file = exposures_list_tValue[i]))[-1]
  exposures_tValue[[i]] <- exposures_tValue[[i]][exposures_tValue[[i]]$pval.exposure<0.001,]
  exposures_tValue[[i]] <- t_value(exposures_tValue[[i]])
  gc()
}
names(exposures_tValue) <- exposures_list_tValue
gc()
i <- NULL

# Outcomes =====
Outcomes_list_tValue <- list.files(path = 
                                     "Methods/Outcomes/")
Outcomes_tValue <- list()

for (i in seq_along(Outcomes_list_tValue)) {
  Outcomes_tValue[[i]] <- read_csv(paste0("Methods/Outcomes/",
                                          file = Outcomes_list_tValue[i]))[-1]
  Outcomes_tValue[[i]] <- Outcomes_tValue[[i]][Outcomes_tValue[[i]]$pval.outcome<0.001,]
  Outcomes_tValue[[i]] <- t_value(Outcomes_tValue[[i]])
  gc()
}
names(Outcomes_tValue) <- Outcomes_list_tValue
gc()
i <- NULL

### Harmonisation ---> LD clumping ---> MR analysis ====
harmonized_set_tValue <- list()
name_harmo_tValue <- list()
clump_set_tValue <- list()
MR_set_tValue <- list()
m <- 0
i <- 0
j <- 0
gc()
for (j in seq_along(Outcomes_list_tValue)) {
  for (i in seq_along(exposures_list_tValue)) {
    if (unlist(map(strsplit(names(Outcomes_tValue[j]),"_"),1)) != unlist(map(strsplit(names(exposures_tValue[i]),"_"),1))){
      m <- m+1
      harmonized_set_tValue[[m]] <- harmonise_data(exposure_dat = exposures_tValue[[i]], outcome_dat = Outcomes_tValue[[j]])
      #name_harmo_tValue[[m]] <- paste0(exposures_list_tValue[i],"_", Outcomes_list_tValue[j])
      # LD clumping
      if (length(harmonized_set_tValue[[m]]$SNP)>1){
        clump_set_tValue[[m]] <- clump_data(harmonized_set_tValue[[m]])
        if (length(clump_set_tValue[[m]] > 1)){
          MR_set_tValue[[m]] <- mr(clump_set_tValue[[m]])
          name_harmo_tValue[[m]] <- unique(paste0(MR_set_tValue[[m]]$exposure,"-",MR_set_tValue[[m]]$outcome))
          print(name_harmo_tValue[[m]])
        }
      } else {
        clump_set_tValue[[m]] <- 0
        MR_set_tValue[[m]] <- 0
      }
      gc()
    }
    gc()
  }
}
names(harmonized_set_tValue) <- name_harmo_tValue
names(clump_set_tValue) <- name_harmo_tValue
names(MR_set_tValue) <- name_harmo_tValue
m <- 0
i <- 0
j <- 0
saveRDS(harmonized_set_tValue, "harmonized_set_t_value.RData")
saveRDS(clump_set_tValue, "clump_set_t_value.RData")
saveRDS(MR_set_tValue, "MR_set_t_value.RData")


gc()
