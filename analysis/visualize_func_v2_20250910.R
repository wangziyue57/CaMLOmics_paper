# generate figures for the checkpoints and overall RPS
# 1. time-dependent ROC and AUC
# 2. KM for combined score
# 3. KM for each probe and exp; 
# 4. violin/boxplot for probe and exp in 3 subgroups
# 5. forest plot for HR from cox model


# -------------------------------------------------------------
#  1. Function: Evaluate Cox Model + Time-Dependent AUC
#            with optional K-fold CV, Bootstrap, or .632 Bootstrap
# function: KM curve from the PRS
# -------------------------------------------------------------
eval_roc <- function(train_data,
                     test_data, 
                     time_col,
                     event_col,
                     predictors,
                     covariates,
                     times,
                     internal_validation = c("none", "cv", "boot", "boot632"),
                     cv_folds            = 5,
                     n_boot              = 50,
                     seed                = 123) {
  
  internal_validation <- match.arg(internal_validation)
  
  # Cox formula
  all_predictors <- c(predictors, covariates)
  cox_formula <- as.formula(
    paste0("Surv(", time_col, ", ", event_col, ") ~ ",
           paste(all_predictors, collapse = " + "))
  )
  
  # Fit Cox model
  cox_model <- coxph(cox_formula, data = train_data, x = TRUE, y = TRUE)
  
  # Add risk scores
  train_data$risk_score <- predict(cox_model, newdata = train_data, type = "lp")
  test_data$risk_score  <- predict(cox_model, newdata = test_data,  type = "lp")
  
  # Apparent ROC (train)
  roc_train_apparent <- timeROC(
    T      = train_data[[time_col]],
    delta  = train_data[[event_col]],
    marker = train_data$risk_score,
    cause  = 1,
    weighting = "marginal",
    times  = times,
    ROC    = TRUE,
    iid    = TRUE
  )
  CI_train_apparent <- confint(roc_train_apparent)
  
  # ROC (test)
  roc_test <- timeROC(
    T      = test_data[[time_col]],
    delta  = test_data[[event_col]],
    marker = test_data$risk_score,
    cause  = 1,
    weighting = "marginal",
    times  = times,
    ROC    = TRUE,
    iid    = TRUE
  )
  CI_test <- confint(roc_test)
  
  # -------------------------
  # Internal Validation
  # -------------------------
  cv_results <- boot_results <- boot632_results <- NULL
  roc_cv <- roc_oob <- roc_app_boot <- NULL
  set.seed(seed)
  
  ## --- K-fold CV
  if (internal_validation == "cv") {
    n <- nrow(train_data)
    folds_vec <- sample(rep(1:cv_folds, length.out = n))
    
    auc_matrix <- matrix(NA, nrow = cv_folds, ncol = length(times))
    roc_cv <- vector("list", cv_folds)
    
    for (f in 1:cv_folds) {
      test_idx <- which(folds_vec == f)
      cv_train <- train_data[-test_idx, ]
      cv_test  <- train_data[test_idx, ]
      
      cox_cv <- coxph(cox_formula, data = cv_train, x = TRUE, y = TRUE)
      cv_test$risk_score <- predict(cox_cv, newdata = cv_test, type = "lp")
      
      roc_cv[[f]] <- timeROC(
        T      = cv_test[[time_col]],
        delta  = cv_test[[event_col]],
        marker = cv_test$risk_score,
        cause  = 1,
        weighting = "cox",
        times  = times,
        ROC    = TRUE,
        iid    = FALSE 
      )
      auc_matrix[f, ] <- roc_cv[[f]]$AUC
    }
    
    auc_mean <- colMeans(auc_matrix, na.rm = TRUE)
    auc_sd   <- apply(auc_matrix, 2, sd, na.rm = TRUE)
    t_crit   <- qt(0.975, df = cv_folds - 1)
    
    ci_lower <- auc_mean - t_crit * auc_sd / sqrt(cv_folds)
    ci_upper <- auc_mean + t_crit * auc_sd / sqrt(cv_folds)
    
    cv_results <- data.frame(
      time_point  = times,
      AUC_mean    = auc_mean,
      AUC_sd      = auc_sd,
      CI_lower    = ci_lower,
      CI_upper    = ci_upper
    )
  }
  
  ## --- Bootstrap & .632 bootstrap
  if (internal_validation %in% c("boot", "boot632")) {
    n <- nrow(train_data)
    auc_oob_matrix <- matrix(NA, nrow = n_boot, ncol = length(times))
    auc_app_matrix <- matrix(NA, nrow = n_boot, ncol = length(times))
    
    roc_oob <- vector("list", n_boot)
    roc_app_boot <- vector("list", n_boot)
    
    for (b in seq_len(n_boot)) {
      boot_idx    <- sample(seq_len(n), size = n, replace = TRUE)
      boot_sample <- train_data[boot_idx, ]
      oob_mask    <- !(seq_len(n) %in% boot_idx)
      oob_data    <- train_data[oob_mask, ]
      
      cox_boot <- coxph(cox_formula, data = boot_sample, x = TRUE, y = TRUE)
      
      # OOB evaluation
      if (nrow(oob_data) > 0) {
        oob_data$risk_score <- predict(cox_boot, newdata = oob_data, type = "lp")
        roc_oob[[b]] <- timeROC(
          T      = oob_data[[time_col]],
          delta  = oob_data[[event_col]],
          marker = oob_data$risk_score,
          cause  = 1,
          weighting = "cox",
          times  = times,
          ROC    = TRUE,
          iid    = FALSE 
        )
        auc_oob_matrix[b, ] <- roc_oob[[b]]$AUC
      }
      
      # Apparent bootstrap (for .632)
      if (internal_validation == "boot632") {
        boot_sample$risk_score <- predict(cox_boot, newdata = boot_sample, type = "lp")
        roc_app_boot[[b]] <- timeROC(
          T      = boot_sample[[time_col]],
          delta  = boot_sample[[event_col]],
          marker = boot_sample$risk_score,
          cause  = 1,
          weighting = "cox",
          times  = times,
          ROC    = TRUE,
          iid    = FALSE 
        )
        auc_app_matrix[b, ] <- roc_app_boot[[b]]$AUC
      }
    }
    
    if (internal_validation == "boot") {
      auc_oob_mean <- colMeans(auc_oob_matrix, na.rm = TRUE)
      auc_ci_lower <- apply(auc_oob_matrix, 2, quantile, probs = 0.025, na.rm = TRUE)
      auc_ci_upper <- apply(auc_oob_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)
      
      boot_results <- data.frame(
        time_point   = times,
        AUC_oob_mean = auc_oob_mean,
        CI_lower     = auc_ci_lower,
        CI_upper     = auc_ci_upper
      )
    }
    
    if (internal_validation == "boot632") {
      auc_632_matrix <- (1 - 0.632)*auc_app_matrix + 0.632*auc_oob_matrix
      auc_632_mean <- colMeans(auc_632_matrix, na.rm = TRUE)
      auc_ci_lower <- apply(auc_632_matrix, 2, quantile, probs = 0.025, na.rm = TRUE)
      auc_ci_upper <- apply(auc_632_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)
      
      boot632_results <- data.frame(
        time_point   = times,
        AUC_632_mean = auc_632_mean,
        CI_lower     = auc_ci_lower,
        CI_upper     = auc_ci_upper
      )
    }
  }
  
  # -------------------------
  # Return
  # -------------------------
  return(list(
    cox_model          = cox_model,
    roc_train_apparent = roc_train_apparent,
    CI_train_apparent  = CI_train_apparent,
    roc_test           = roc_test,
    CI_test            = CI_test,
    roc_cv             = roc_cv,
    cv_results         = cv_results,
    roc_oob            = roc_oob,
    boot_results       = boot_results,
    roc_app_boot       = roc_app_boot,
    boot632_results    = boot632_results,
    train_data         = train_data,  # with risk_score
    test_data          = test_data    # with risk_score
  ))
}


plot_km <- function(dat, time_col, event_col,
                    cutoff = 0.5,
                    cutoff_type = c("median", "quantile"),
                    stratify_var = NULL,
                    palette = c("#E64B35FF", "#1F78A0")) {
  
  cutoff_type <- match.arg(cutoff_type)
  
  # 1. Determine cutoff value
  cutoff_value <- if (cutoff_type == "median") {
    median(dat$risk_score, na.rm = TRUE)
  } else {
    quantile(dat$risk_score, cutoff, na.rm = TRUE)
  }
  
  # 2. Assign risk groups
  dat$risk_group <- ifelse(dat$risk_score >= cutoff_value, "High risk", "Low risk")
  
  # 3. Build survival formula
  surv_formula <- as.formula(
    paste0("Surv(", time_col, ", ", event_col, ") ~ risk_group")
  )
  
  # 4. Stratified analysis
  if (!is.null(stratify_var)) {
    strat_levels <- unique(dat[[stratify_var]])
    plots <- list()
    
    for (lev in strat_levels) {
      dat_sub <- subset(dat, dat[[stratify_var]] == lev)
      fit <- surv_fit(surv_formula, data = dat_sub)
      
      plots[[as.character(lev)]] <- ggsurvplot(
        fit, data = dat_sub,
        conf.int = FALSE, risk.table = TRUE,
        xlab = "Time (days)", #xscale = 365.25,
        legend.title = "PRS",
        legend.labs = c("High", "Low"),
        surv.median.line = "hv",
        pval = TRUE,
        risk.table.y.text.col = TRUE, risk.table.y.text = FALSE,
        #title = paste(stratify_var, "=", lev),
        palette = palette
      )
    }
    return(plots)
    
  } else {
    fit <- surv_fit(surv_formula, data = dat)
    
    return(
      ggsurvplot(
        fit, data = dat,
        conf.int = FALSE, risk.table = TRUE,
        xlab = "Time (days)", #xscale = 365.25,
        legend.title = "PRS",
        legend.labs = c("High", "Low"),
        surv.median.line = "hv",
        pval = TRUE,
        risk.table.y.text.col = TRUE, risk.table.y.text = FALSE,
        #title = "Overall",
        palette = palette
      )
    )
  }
}




# -------------------------------------------------------------
#  2 Function: KM for single biomarker
#     chose different cutoff for each subgroup
# -------------------------------------------------------------
plot_km_single_subgroup <- function(dat, biomarker,
                                    cutoff, cutoff1, cutoff2, cutoff3,
                                    sample1, sample2, sample3){
  
  # Create groups based on the cutoff
  cutoff_value <- quantile(dat[[biomarker]], cutoff, na.rm = TRUE)
  dat <- dat %>%
    mutate(biomarker_group = ifelse(!!sym(biomarker) >= cutoff_value, "High", "Low"))
  
  # Subset data for each subgroup and compute risks 
  cutoff_value <- quantile(dat[sample1,biomarker], cutoff1, na.rm = TRUE)
  dat1 <- dat[sample1, ] %>%
    mutate(biomarker_group = ifelse(!!sym(biomarker) >= cutoff_value, "High", "Low"))
  
  cutoff_value <- quantile(dat[sample2,biomarker], cutoff2, na.rm = TRUE)
  dat2 <- dat[sample2, ] %>%
    mutate(biomarker_group = ifelse(!!sym(biomarker) >= cutoff_value, "High", "Low"))
  
  cutoff_value <- quantile(dat[sample3,biomarker], cutoff3, na.rm = TRUE)
  dat3 <- dat[sample3, ] %>%
    mutate(biomarker_group = ifelse(!!sym(biomarker) >= cutoff_value, "High", "Low"))
  
  
  # Fit the Kaplan-Meier survival model stratified by biomarker_group
  fit <- survfit(Surv(overall_survival, status) ~ biomarker_group, data = dat)
  fit1 <- survfit(Surv(overall_survival, status) ~ biomarker_group, data = dat1)
  fit2 <- survfit(Surv(overall_survival, status) ~ biomarker_group, data = dat2)
  fit3 <- survfit(Surv(overall_survival, status) ~ biomarker_group, data = dat3)
  
  # Plot the KM curve
  ggsurv <- list()
  ggsurv[[1]] <- ggsurvplot(fit, data = dat,
                            conf.int = FALSE, risk.table = TRUE,
                            xlab = "Time (days)", #xscale = 365.25,
                            #title = biomarker,
                            legend.title = biomarker,
                            legend.labs = c("High", "Low"),
                            surv.median.line = "hv",
                            pval = TRUE)
  
  ggsurv[[2]] <- ggsurvplot(fit1, data = dat1,
                            conf.int = FALSE, risk.table = TRUE,
                            xlab = "Time (days)", #xscale = 365.25,
                            #title = biomarker,
                            legend.title = biomarker,
                            legend.labs = c("High", "Low"),
                            surv.median.line = "hv",
                            pval = TRUE)
  
  ggsurv[[3]] <- ggsurvplot(fit2, data = dat2,
                            conf.int = FALSE, risk.table = TRUE,
                            xlab = "Time (days)", #xscale = 365.25,
                            #title = biomarker,
                            legend.title = biomarker,
                            legend.labs = c("High", "Low"),
                            surv.median.line = "hv",
                            pval = TRUE)
  
  ggsurv[[4]] <- ggsurvplot(fit3, data = dat3,
                            conf.int = FALSE, risk.table = TRUE,
                            xlab = "Time (days)", #xscale = 365.25,
                            #title = biomarker,
                            legend.title = biomarker,
                            legend.labs = c("High", "Low"),
                            surv.median.line = "hv",
                            pval = TRUE)
  
  return(ggsurv = ggsurv)
}


# -------------------------------------------------------------
#  3.1 Function: boxplot of cpg and mrna across 3 subgroups
# -------------------------------------------------------------
plot_boxplot_prob_mrna <- function(dat, probe, genes) {
  
  clean_theme <- theme_classic() +
    theme(
      plot.title   = element_text(size = 12),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 10)
    )
  
  # --- Boxplot for probe (always single) ---
  p_boxplot_probe <- ggplot(dat, aes(x = group, y = !!sym(probe), fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = probe, y = "Beta value", fill = "Group") +
    scale_fill_manual(labels = c("IDH-MT",
                                   "IDH-WT, mMGMT",
                                   "IDH-WT, unMGMT"),
                        values = c("MT" = "#1b9e77",
                                   "WT, Methy" = "#d95f02",
                                   "WT, Unmethy" = "#7570b3")) +
    clean_theme
  
  # --- Boxplots for each gene in the vector ---
  p_boxplot_genes <- lapply(genes, function(gene) {
    ggplot(dat, aes(x = group, y = !!sym(gene), fill = group)) +
      geom_boxplot(outlier.shape = NA) +
      labs(title = gene, y = "Expr value", fill = "Group") +
      scale_fill_manual(labels = c("IDH-MT",
                                   "IDH-WT, mMGMT",
                                   "IDH-WT, unMGMT"),
                        values = c("MT" = "#1b9e77",
                                   "WT, Methy" = "#d95f02",
                                   "WT, Unmethy" = "#7570b3")) +
      clean_theme
  })
  
  # --- Return as a list ---
  return(c(list(p_boxplot_probe), p_boxplot_genes))
}


# -------------------------------------------------------------
#  3.2 Function: correlation of cppg vs mrna
# -------------------------------------------------------------
plot_corr_prob_mrna <- function(dat, probe, anno_gene, med_gene){
  ## correlation
  probe_vals <- rlang::eval_tidy(sym(probe), dat)
  anno_gene_vals <- rlang::eval_tidy(sym(anno_gene), dat)
  med_gene_vals <- rlang::eval_tidy(sym(med_gene), dat)
  
  ct_anno <- cor.test(probe_vals, anno_gene_vals, method = "spearman")
  corr_value_anno <- round(ct_anno$estimate, 3)
  pvalue_anno <- round(ct_anno$p.value, 3)
  
  ct_med <- cor.test(probe_vals, med_gene_vals, method = "spearman")
  corr_value_med <- round(ct_med$estimate, 3)
  pvalue_med <- round(ct_med$p.value, 3)
  
  # Create the scatter plot with a linear regression line per group
  p_scatter_anno <- ggplot(dat, aes(x = !!sym(probe), y = !!sym(anno_gene), color = group)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(#title = paste0("r = ", corr_value_anno, ", p =", pvalue_anno),
      x = probe,
      y = anno_gene, color = "Group") +
    scale_color_manual(labels = c("IDH-MT",
                                 "IDH-WT, mMGMT",
                                 "IDH-WT, unMGMT"),
                      values = c("MT" = "#1b9e77",
                                 "WT, Methy" = "#d95f02",
                                 "WT, Unmethy" = "#7570b3")) +
    theme_classic()
  
  p_scatter_med <- ggplot(dat, aes(x = !!sym(probe), y = !!sym(med_gene), color = group)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(#title = paste0("r = ", corr_value_med, ", p =", pvalue_med),
      x = probe,
      y = med_gene, color = "Group") +
    scale_color_manual(labels = c("IDH-MT",
                                 "IDH-WT, mMGMT",
                                 "IDH-WT, unMGMT"),
                      values = c("MT" = "#1b9e77",
                                 "WT, Methy" = "#d95f02",
                                 "WT, Unmethy" = "#7570b3")) +
    theme_classic()
  
  return(list(p_scatter_anno, p_scatter_med))
  
}


# -------------------------------------------------------------
#  4. Function: forest plot of HR
# multivariate cox from whole sample and stratified sample
# -------------------------------------------------------------
plot_forest <- function(dat, 
                        time_col,
                        event_col,
                        probe, anno_gene, med_gene,
                        covariates,
                        sample){
  
  cox_formula <- as.formula(paste("Surv(", time_col, ", ", event_col, ") ~ ", 
                                  probe, "+", 
                                  paste(anno_gene, collapse = "+"), "+", 
                                  paste(med_gene, collapse = "+"), "+", 
                                  paste(covariates, collapse = " + ")))
  
  
  # only show one subgroup results
  # Fit the Cox models:
  cox_model <- coxph(cox_formula, data = dat[sample,])

  # Extract results
  cox_results <- summary(cox_model)
  
  # Create dataframe for plotting
  results <- data.frame(
    Variable = rownames(cox_results$coefficients),
    HR = exp(cox_results$coefficients[, "coef"]),
    lower = cox_results$conf.int[, "lower .95"],
    upper = cox_results$conf.int[, "upper .95"],
    p = cox_results$coefficients[, "Pr(>|z|)"]
  )

  # Format labels
  results <- results %>%
    mutate(
      HR_label = sprintf("%.2f (%.2f–%.2f)", HR, lower, upper),
      p_label = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
    )
  
  # Forest plot
  p_forest <- ggplot(results, aes(x = Variable, y = HR, ymin = lower, ymax = upper)) +
    geom_pointrange(shape=15,size = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_text(aes(label = sprintf("p = %.3f", p)), 
              hjust = -0.2, size = 3) +
    coord_flip() +
    theme_classic() +
    theme(
      axis.ticks.y= element_blank(),
      axis.line.y.left = element_blank(),
      axis.title.y= element_blank(),
      axis.text.y        = element_text(face = "bold",size = 12),
      axis.text.x= element_text(size = 11)) + 
    labs(
      x = "",
      y = "Hazard Ratio (95% CI)"
    ) 
  
  return(p_forest)
}

