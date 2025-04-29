# SVR-MR (MR with Machine-learning)

## Required Packages
  1) e1071      # For SVM implementation
  2) caret      # For data splitting
  3) dplyr      # For data manipulation
  4) TwoSampleMR   # For MR analysis
  5) meta      # For meta-analysis and plotting
  6) readr     # For reading csv

## Default_MR_test_finngen.R
* To prepare data for SVM analysis and
* To perform default MR analysis
* Prepare Exposure and Outcome data based on trait of the study
    * Here FinnGen study datasets are used (Data is available on request to the FinnGen project investigators at https://www.finngen.fi/en/access_results)

## SV_IV_MR.R
* IVs are selected using SVM
* SVR-MR is implemented
