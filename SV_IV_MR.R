set.seed()
# Load required libraries
library(e1071)
library(MendelianRandomization)
library(caret)  # For data splitting
library(dplyr)  # For data manipulation
library(TwoSampleMR)
library(meta)
library(ggplot2)
source("Methods/ivw.R")

# Load the GWAS summary statistics data 
clump_set_t_value <- readRDS("clump_set_t_value.RData")
MR_set_Default <- readRDS("MR_set.RData") 


# Needed results
data <- NULL
mse <- list()
r2 <- list()
svr_model <- list()
name_list <- list()
clump_set_SV <- list()
clump_set_SV_MR <- list()
genes_SV <- list()
odds <- list()
odds_svm <- list()
odds_default <- list()
pleiotropy <- list()
mr_sim <- list()
columns= c("nsnp","b","se","pval") 
svm_result <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(svm_result) <- columns
columns_nsnp = c("Default_IV","Default_IV_MR","SV_IV","SV_IV_MR", "F_statistics_SV") 
nSNP_result <- data.frame(matrix(nrow = 0, ncol = length(columns_nsnp)))
colnames(nSNP_result) <- columns_nsnp

for (i in seq_along(clump_set_t_value)) {
  if (length(clump_set_t_value[[i]])>1){
    if(length(clump_set_t_value[[i]]$SNP)>3){
      data <- clump_set_t_value[[i]]
      data$pval.exposure <- -log10(data$pval.exposure)
      data$pval.outcome <- -log10(data$pval.outcome)
      # Use independent features as predictors (SE, P-value, and Allele Frequency)
      X <- data %>% select(beta.exposure, eaf.exposure, eaf.outcome, 
                           pval.outcome, pval.exposure,
                           se.outcome, se.exposure,
                           t_value.x, t_value.y)
      y <- data$beta.outcome 
      rownames(X) <- data$SNP
      print("++++++++++++++++++++++++++++++++++++++++++++++++")
      # SVR Model
      svr_model[[i]] <- svm(X, y, kernel = "linear",
                            ranges=list(epsilon=seq(0,1,0.1),gamma = seq(0.5,5,0.5),
                                        cost=seq(1,100,1)))
      # Predict on the test set
      y_pred <- predict(svr_model[[i]], X)
      # Calculate performance metrics (R^2 and MSE)
      r2[[i]] <- cor(y, y_pred)^2
      mse[[i]] <- mean((y - y_pred)^2)
      rownames(svr_model[[i]]$coefs) <-  rownames(svr_model[[i]]$SV)
      clump_set_SV[[i]] <- clump_set_t_value[[i]][clump_set_t_value[[i]]$SNP %in% 
                                                    rownames(svr_model[[i]]$SV),]
      genes_SV[[i]] <- clump_set_SV[[i]]$gene.exposure
      clump_set_SV_MR[[i]] <- mr(clump_set_SV[[i]])
      pleiotropy[[i]] <- mr_pleiotropy_test(clump_set_SV[[i]])
      mr_sim[[i]] <- mr_simex(clump_set_SV[[i]])
      odds[[i]] <- generate_odds_ratios(clump_set_SV_MR[[i]]
                                        [clump_set_SV_MR[[i]]$method=="Inverse variance weighted",])
      name_list[[i]] <- unique(paste0(clump_set_t_value[[i]]$exposure,"-",clump_set_t_value[[i]]$outcome))
      svm_ivw_MR <- ivw_svm(clump_set_SV[[i]]$beta.exposure, clump_set_SV[[i]]$beta.outcome,
                            clump_set_SV[[i]]$se.exposure, clump_set_SV[[i]]$se.outcome)
      svm_result[i,] <- cbind(svm_ivw_MR$nsnp, svm_ivw_MR$b, svm_ivw_MR$se, svm_ivw_MR$pval)
      row.names(svm_result)[i] <- name_list[[i]]
      odds_svm[[i]] <- generate_odds_ratios(svm_result[i,])
      print(name_list[[i]])
      # number of SNPs
      nSNP_result[i,] <- cbind(length(clump_set_t_value[[i]]$SNP), 
                               ifelse(length(MR_set_Default[[i]])>1, unique(MR_set_Default[[i]]$nsnp), 0), 
                               length(clump_set_SV[[i]]$SNP), 
                               odds_svm[[i]]$nsnp,
                               mean(clump_set_SV[[i]]$t_value.y**2))
      row.names(nSNP_result)[i] <- name_list[[i]]
      if(length(MR_set_Default[[i]])>1){
        odds_default[[i]] <- generate_odds_ratios(MR_set_Default[[i]]
                                                  [MR_set_Default[[i]]$method=="Inverse variance weighted",])
      }
    }
  } 
}
names(r2) <- name_list
names(mse) <- name_list
names(svr_model) <- name_list
names(clump_set_SV_MR) <- name_list
names(clump_set_SV) <- name_list
names(pleiotropy) <- name_list
names(mr_sim) <- name_list
names(genes_SV) <- name_list
names(odds) <- name_list
names(odds_default) <- names(MR_set_Default)
names(odds_svm) <- name_list

gc()

saveRDS(clump_set_SV_MR, "MR_set_SVR.RData")
saveRDS(clump_set_SV, "clump_set_SV.RData")
write.csv(svm_result, "FinnGen_SVM_result.csv", quote = F)

gc()