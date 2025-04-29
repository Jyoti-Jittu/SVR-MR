ivw_svm <- function (b_exp, b_out, se_exp, se_out, parameters = default_parameters()) 
{
  if (sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & 
          !is.na(se_out)) < 2) 
    return(list(b = NA, se = NA, pval = NA, nsnp = NA))
  # Incorporating Support Vectors
  ivw.res <- summary(stats::lm(b_out ~ -1 + b_exp, 
                               weights = (1/svr_model[[i]]$coefs^2) * (1/se_out^2)))
  b <- ivw.res$coef["b_exp", "Estimate"]
  se <- ivw.res$coef["b_exp", "Std. Error"]/min(1, ivw.res$sigma)
  pval <- 2 * pnorm(abs(b/se), lower.tail = FALSE)
  Q_df <- length(b_exp) - 1
  Q <- ivw.res$sigma^2 * Q_df
  Q_pval <- pchisq(Q, Q_df, lower.tail = FALSE)
  return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), 
              Q = Q, Q_df = Q_df, Q_pval = Q_pval))
}


