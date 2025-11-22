
source(paste0(CODE,"visualize_func_v2_20250910.R"))


##################################
## load data
##################################
# omics data
dat.met.glioma.norm <- readRDS(file = paste0(DATA, "glioma_met_clean_norm_20250209.rds")) ###beta value
dat.exp.glioma.norm <- readRDS(file = paste0(DATA, "glioma_exp_clean_norm_20250209.rds")) ###normalized count

# clinical data
dat.clin.glioma.met <- readRDS(file = paste0(RESULT,"dat_met_clin_20250205.rds"))
dat.clin.glioma.exp <- readRDS(file = paste0(RESULT,"dat_exp_clin_20250205.rds"))

sample_id <- unique(c(intersect(dat.clin.glioma.met$sample, dat.clin.glioma.exp$sample)))
dat.clin.glioma.inte <- dat.clin.glioma.met %>%
  dplyr::select(sample) %>%
  full_join(dat.clin.glioma.exp, by="sample") %>%
  filter(sample %in% sample_id)

# combine age-scale and radiation
dat.clin.glioma.inte$treatments[[2]][2,]$treatment_or_therapy=="no"
dat.clin.glioma.inte <- dat.clin.glioma.inte %>%
  mutate(radiation = sapply(dat.clin.glioma.inte$treatments, function(x) x[2,]$treatment_or_therapy!="no")) %>%
  mutate(age_scale = as.vector(scale(age_at_diagnosis))) 
rownames(dat.clin.glioma.inte) <- dat.clin.glioma.inte$sample

# defina a new group variable by combination of IDH and MGMT
dat.clin.glioma.inte <- dat.clin.glioma.inte %>%
  mutate(group = ifelse(paper_IDH.status == "WT" & paper_MGMT.promoter.status == "Unmethylated", "WT, Unmethy",
                        ifelse(paper_IDH.status == "WT" & paper_MGMT.promoter.status == "Methylated", "WT, Methy",
                               "MT")))
sampleID.1 <- which(rownames(dat.clin.glioma.inte) %in% rownames(dat.clin.glioma.inte %>% filter(group == "WT, Unmethy" )))
sampleID.2 <- which(rownames(dat.clin.glioma.inte) %in% rownames(dat.clin.glioma.inte %>% filter(group == "WT, Methy")))
sampleID.3 <- which(rownames(dat.clin.glioma.inte) %in% rownames(dat.clin.glioma.inte %>% filter(group == "MT")))

# create the integrated dataset with methy, mrna and clinical factors
# normalize and standardize (?) methylation and gene expression
# converting DNA methylation probes beta values to m_values
dat.met.glioma.toUse <- dat.met.glioma.norm[,rownames(dat.clin.glioma.inte)]
dat.exp.glioma.toUse <- dat.exp.glioma.norm[,rownames(dat.clin.glioma.inte)]

## raw scale of omics data (beta value, and log-exp)
dat.integrate.raw <- data.frame(cancer = dat.clin.glioma.inte$cancer,
                                age_at_diagnosis = dat.clin.glioma.inte$age_at_diagnosis, 
                                age_scale = dat.clin.glioma.inte$age_scale,
                                gender = factor(dat.clin.glioma.inte$gender), 
                                grade = factor(dat.clin.glioma.inte$paper_Grade), 
                                overall_survival = dat.clin.glioma.inte$overall_survival,
                                status = dat.clin.glioma.inte$status,
                                IDH = dat.clin.glioma.inte$paper_IDH.status,
                                ch1p19q = dat.clin.glioma.inte$paper_X1p.19q.codeletion,
                                MGMT = dat.clin.glioma.inte$paper_MGMT.promoter.status,
                                group = dat.clin.glioma.inte$group,
                                radiation = dat.clin.glioma.inte$radiation,
                                t(rbind(dat.met.glioma.toUse, log1p(dat.exp.glioma.toUse))))

## scaled omics data
dat.met.glioma.toUse.mval.scale <- t(apply(dat.met.glioma.toUse, 1, function(x) {
  # Copy the row vector
  new_x <- log2(x/(1-x))
  # Scale within subgroup 1
  new_x[sampleID.1] <- as.numeric(scale(x[sampleID.1]))
  # Scale within subgroup 2
  new_x[sampleID.2] <- as.numeric(scale(x[sampleID.2]))
  # Scale within subgroup 3
  new_x[sampleID.3] <- as.numeric(scale(x[sampleID.3]))
  new_x
}))

dat.exp.glioma.toUse.norm.scale <- t(apply(dat.exp.glioma.toUse, 1, function(x) {
  # Copy the row vector
  new_x <- log1p(x)
  # Scale within subgroup 1
  new_x[sampleID.1] <- as.numeric(scale(x[sampleID.1]))
  # Scale within subgroup 2
  new_x[sampleID.2] <- as.numeric(scale(x[sampleID.2]))
  # Scale within subgroup 3
  new_x[sampleID.3] <- as.numeric(scale(x[sampleID.3]))
  new_x
}))

dat.integrate.scale <- data.frame(cancer = dat.clin.glioma.inte$cancer,
                                  age_at_diagnosis = dat.clin.glioma.inte$age_at_diagnosis, 
                                  age_scale = dat.clin.glioma.inte$age_scale,
                                  gender = factor(dat.clin.glioma.inte$gender), 
                                  grade = factor(dat.clin.glioma.inte$paper_Grade), 
                                  overall_survival = dat.clin.glioma.inte$overall_survival,
                                  status = dat.clin.glioma.inte$status,
                                  IDH = dat.clin.glioma.inte$paper_IDH.status,
                                  ch1p19q = dat.clin.glioma.inte$paper_X1p.19q.codeletion,
                                  MGMT = dat.clin.glioma.inte$paper_MGMT.promoter.status,
                                  group = dat.clin.glioma.inte$group,
                                  radiation = dat.clin.glioma.inte$radiation,
                                  t(rbind(dat.met.glioma.toUse.mval.scale, dat.exp.glioma.toUse.norm.scale)))



##################################
# helper function
##################################
fig_theme <- theme_classic() +
  theme(
    plot.title   = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    #axis.ticks = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10)
  )

plot_km_roc_stratify <- function(dat,
                                 time_col = "overall_survival",
                                 event_col = "status",
                                 predictors,
                                 covariates,
                                 treatment_col = "radiation",
                                 times,
                                 palette = c("#E64B35FF", "#1F78A0")) {
  
  # --- 1. Handle treatment variable automatically ---
  if (is.logical(dat[[treatment_col]])) {
    dat[[treatment_col]] <- factor(dat[[treatment_col]],
                                   levels = c(FALSE, TRUE),
                                   labels = c("Non-radiated", "Radiated"))
  } else if (is.numeric(dat[[treatment_col]])) {
    dat[[treatment_col]] <- factor(dat[[treatment_col]],
                                   levels = c(0, 1),
                                   labels = c("Non-radiated", "Radiated"))
  } else {
    dat[[treatment_col]] <- factor(dat[[treatment_col]])
  }
  
  # --- 2. Fit Cox model for PRS / risk score ---
  all_predictors <- c(predictors, covariates)
  cox_formula <- as.formula(
    paste0("Surv(", time_col, ", ", event_col, ") ~ ",
           paste(all_predictors, collapse = " + "))
  )
  cox_model <- coxph(cox_formula, data = dat, x = TRUE, y = TRUE)
  dat$risk_score <- predict(cox_model, newdata = dat, type = "lp")
  
  # --- 3. Define PRS groups by median ---
  dat$PRS_group <- ifelse(
    dat$risk_score >= median(dat$risk_score, na.rm = TRUE),
    "High", "Low"
  )
  
  # --- 4. Cox model for PRS × treatment interaction ---
  inter_formula <- as.formula(
    paste0("Surv(", time_col, ", ", event_col, ") ~ risk_score * ",
           treatment_col)
  )
  cox_inter <- coxph(inter_formula, data = dat)
  inter_p <- signif(coef(summary(cox_inter))[grep("risk_score:", rownames(coef(summary(cox_inter)))), "Pr(>|z|)"], 3)
  
  # --- 5. KM fits per treatment group ---
  fits <- dat %>%
    group_split(!!sym(treatment_col)) %>%
    setNames(levels(dat[[treatment_col]])) %>%
    lapply(function(df) {
      surv_fit(as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ PRS_group")),
               data = df)
    })
  
  # --- 6. KM plots for each treatment arm ---
  ggsurv <- mapply(function(sfit, arm) {
    ggsurvplot(
      sfit,
      data = dat[dat[[treatment_col]] == arm, ],
      conf.int = TRUE,
      risk.table = TRUE,
      pval = TRUE,
      risk.table.y.text.col = TRUE,
      risk.table.y.text = FALSE,
      xlab = "Time (days)",
      legend.labs = c("High", "Low"),
      legend.title = "PRS",
      surv.median.line = "hv",
      title = arm,
      palette = palette,
      ggtheme = fig_theme
    )
  }, fits, names(fits), SIMPLIFY = FALSE)
  
  # --- 7. Compute ROC per treatment group ---
  roc_train_apparent <- dat %>%
    group_split(!!sym(treatment_col)) %>%
    setNames(levels(dat[[treatment_col]])) %>%
    lapply(function(df) {
      timeROC(
        T = df[[time_col]],
        delta = df[[event_col]],
        marker = df$risk_score,
        cause = 1,
        weighting = "marginal",
        times = times,
        ROC = TRUE,
        iid = TRUE
      )
    })
  CI_train_apparent <- lapply(roc_train_apparent, function(x) confint(x))
  
  return(list(ggsurv=ggsurv,
              cox_model          = cox_model,
              cox_interaction    = cox_inter,
              roc_train_apparent = roc_train_apparent,
              CI_train_apparent  = CI_train_apparent))
}

roc_panel_aligned_stratify <- function(res,
                                       times,                     
                                       title,
                                       cols = c("#1B9E77", "#D95F02", "#7570B3"),
                                       mar  = c(5.0, 5.4, 2.2, 1.6),  # ⬅️ more top/bottom margin
                                       mgp  = c(2.2, 0.6, 0),
                                       tcl  = -0.3,
                                       cex_axis = 0.8,   # ≈ pt
                                       cex_lab  = 1,   # ≈ pt
                                       legend_xy = c(0.55, 0.14),      
                                       cex_leg = 1) {  # ≈ pt
  
  roc_list <- res$roc_train_apparent
  ci_list  <- res$CI_train_apparent
  if (!is.list(roc_list)) stop("res$roc_train_apparent must be a list.")
  if (!is.list(ci_list))  stop("res$CI_train_apparent must be a list.")
  
  roc_plots <- mapply(function(roc, ci, name) {
    
    ggplotify::as.ggplot(function() {
      op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
      par(mar = mar, mgp = mgp, tcl = tcl,
          cex.axis = cex_axis, cex.lab = cex_lab,
          font.lab = 1, font.axis = 1, family = "sans")  # plain text
      
      ## 1. Base ROC curve
      plot(roc, time = times[1], col = cols[1],
           title = FALSE, xlab = "1 – Specificity", ylab = "Sensitivity",
           lwd = 2.2)
      if (length(times) > 1) {
        for (i in 2:length(times)) {
          plot(roc, time = times[i], add = TRUE, col = cols[i],
               title = FALSE, lwd = 2.2)
        }
      }
      
      ## 2. Legend text
      idx <- match(times, roc$times)
      yrs <- round(times / 365)
      auc <- roc$AUC[idx]
      ci$CI_AUC[idx, 2] <- pmin(ci$CI_AUC[idx, 2], 100)
      lbl <- paste0(yrs, "yr, ",
                    sprintf("%.2f", auc), " (",
                    sprintf("%.2f", ci$CI_AUC[idx, 1] / 100), ", ",
                    sprintf("%.2f", ci$CI_AUC[idx, 2] / 100), ")")
      
      legend(legend_xy[1], legend_xy[2],
             legend = c("AUC (95% CI)", lbl),
             col    = c(NA, cols[seq_along(times)]),
             lty    = c(NA, rep(1, length(times))),
             lwd    = c(NA, rep(2.2, length(times))),
             xjust = 0, yjust = 0,
             bty = "n", cex = cex_leg, y.intersp = 1.05,
             x.intersp = 0.8, seg.len = 2.5, text.font = 1)
      
      ## 3. Panel title (plain)
      mtext(name, side = 3, adj = 0, line = 0.5, font = 1, cex = 1)
    }) +
      theme(plot.margin = margin(12, 5, 12, 5))  # outer padding
  },
  roc_list, ci_list, names(roc_list), SIMPLIFY = FALSE)
}


##################################
# cg00884093(CELSR1); all mediators
# CEACAM1,MSH2,MAP3K1,MAFB,LINC00671,AC131097.3
# train in MT
##################################
biomarkers <- c("cg00884093","CELSR1","CEACAM1","MSH2","MAP3K1","MAFB","LINC00671","AC131097.3")
km_mt <- plot_km_roc_stratify(
  dat = dat.integrate.scale %>%
    filter(group == "MT"),
  time_col = "overall_survival",
  event_col = "status",
  treatment_col = "radiation",     # logical TRUE/FALSE is fine now
  predictors = biomarkers,
  covariates = c("age_scale","gender","grade","MGMT"),
  times      = c(3*365, 5*365, 10*365))

##################################
# fig 5a. box plot
##################################
p_boxplot <- plot_boxplot_prob_mrna(dat.integrate.raw, 
                                    probe = "cg00884093",
                                    genes = c("CELSR1","CEACAM1","MSH2","MAP3K1","MAFB","LINC00671","AC131097.3"))
p_boxplot_all <- ggarrange(plotlist=p_boxplot, 
                           ncol = 4, nrow = 2, 
                           legend = "top", common.legend=TRUE)
p_boxplot_all <- ggarrange(plotlist=p_boxplot, 
                           ncol = 8, nrow = 1, 
                           legend = "top", common.legend=TRUE)

##################################
# fig 5b. forest plot
##################################
# Extract relevant statistics
model_summary <- summary(km_mt$cox_model)
results <- data.frame(
  term = rownames(model_summary$coefficients),  # Variable names
  estimate = exp(model_summary$coefficients[, "coef"]),  # Convert log HR to HR
  conf.low = model_summary$conf.int[, "lower .95"],  # Lower 95% CI
  conf.high = model_summary$conf.int[, "upper .95"],  # Upper 95% CI
  p.value = model_summary$coefficients[, "Pr(>|z|)"]  # p-values
)
results <- results %>%
  mutate(term = recode(term,
                       "age_scale" = "Age",
                       "gendermale" = "Male",
                       "gradeG3" = "Grade G3",
                       "gradeG4" = "Grade G4",
                       "MGMTUnmethylated" = "MGMT unmethylated",
                       "MAP3K1" = "MAP3K1",
                       "MAFB" = "MAFB",
                       "LINC00671" = "LINC00671",
                       "MSH2" = "MSH2",
                       "CEACAM1" = "CEACAM1",
                       "CELSR1" = "CELSR1",
                       "cg00884093" = "cg00884093",
                       "AC131097.3" = "AC131097.3"),
         term = factor(term,
                        levels = c(
                          "MGMT unmethylated",
                          "Grade G4",
                          "Grade G3",
                          "Male",
                          "Age",
                          "AC131097.3",
                          "LINC00671",
                          "MAFB",
                          "MAP3K1",
                          "MSH2",
                          "CELSR1",
                          "CEACAM1",
                          "cg00884093")),
         #term <- forcats::fct_rev(term),
         HR_CI = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high),
         p.value = sprintf("%.3f", p.value))


# Create the forest plot
# not show G4, since therea re only 3 samples and yield a large CI range
hr_breaks <- c(0, 0.5, 1, 2, 5) # choose ticks you like
# Compute rightward offset dynamically (to keep space but not stretch axis)
pval_offset <- max(results[-12,]$conf.high, na.rm = TRUE) * 1.25

p_forest <- ggplot(results[-12,], aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange(shape=15) +  # Points and error bars
  #geom_point(shape=15) + 
  #geom_errorbar() +  # Points and error bars
  geom_text(aes(label = HR_CI), vjust = -1, size = 3.5) +  # Add HR labels
  # # p-value column aligned to right side
  geom_text(
    aes(x = term, y = max(conf.high, na.rm = TRUE) * 1.2,  # position to the right
        label = ifelse(p.value < 0.001, "<0.001", p.value)),
    hjust = 0, size = 3.5
  ) +
  # Add column title "p-value" at top
  annotate("text", x = 12.5, y = pval_offset, 
           label = "p-value", hjust = 0.3, vjust = -1, fontface = "bold", size = 4) +
  coord_flip(clip = "off") +  # Flip coordinates for a vertical forest plot
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  # Reference line at HR = 1
  scale_y_log10(breaks = hr_breaks, labels = hr_breaks) +
  #scale_y_continuous(breaks = hr_breaks, labels = hr_breaks) +
  labs(
    title = "",
    x = "",
    y = "Hazard Ratio (HR), 95% CI"
  ) +
  theme_classic() +
  theme(
    axis.ticks.y= element_blank(),
    axis.line.y.left = element_blank(),
    axis.title.y= element_blank(),
    axis.title = element_text(size = 12),
    axis.text   = element_text(size = 10),
    plot.margin = margin(12, 50, 12, 5))  # ← increased right margin for p-values

  
##################################
# fig 5c. KM for RPS; stratify in radiation
##################################
# KM for PRS within radition vs. non radition
p_km_rad     <- km_mt$ggsurv$Radiated$plot + ggtitle("Radiated")
p_km_nonrad  <- km_mt$ggsurv$`Non-radiated`$plot + ggtitle("Non-radiated")

p_km <- ggarrange(p_km_rad, p_km_nonrad, 
                  ncol = 2, nrow = 1)


##################################
# fig 5e. ROC: stratify in radiation
##################################
p_roc_mt <- roc_panel_aligned_stratify(km_mt, times = c(3, 5, 10) * 365, 
                           legend_xy = c(0.2, 0))
p_timeROC <- ggarrange(p_roc_mt$Radiated, p_roc_mt$`Non-radiated`,
                       ncol = 2, nrow = 1)


##################################
# cg01878435(CD86); 
# H2BC11
# train in MT
##################################
biomarkers <- c("cg01878435","CD86","H2BC11")
km_mt <- plot_km_roc_stratify(
  dat = dat.integrate.scale %>%
    filter(group == "MT"),
  time_col = "overall_survival",
  event_col = "status",
  treatment_col = "radiation",     # logical TRUE/FALSE is fine now
  predictors = biomarkers,
  covariates = c("age_scale","gender","grade","MGMT"),
  times      = c(3*365, 5*365, 10*365))


##################################
# fig 5a. box plot
##################################
p_boxplot <- plot_boxplot_prob_mrna(dat.integrate.raw, 
                                    probe = "cg01878435",
                                    genes = c("CD86","H2BC11"))
p_boxplot_all <- ggarrange(plotlist=p_boxplot, 
                           ncol = 3, nrow = 1, 
                           legend = "top", common.legend=TRUE)

##################################
# fig 5b. forest plot
##################################
# Extract relevant statistics
model_summary <- summary(km_mt$cox_model)
results <- data.frame(
  term = rownames(model_summary$coefficients),  # Variable names
  estimate = exp(model_summary$coefficients[, "coef"]),  # Convert log HR to HR
  conf.low = model_summary$conf.int[, "lower .95"],  # Lower 95% CI
  conf.high = model_summary$conf.int[, "upper .95"],  # Upper 95% CI
  p.value = model_summary$coefficients[, "Pr(>|z|)"]  # p-values
)
results <- results %>%
  mutate(term = recode(term,
                       "age_scale" = "Age",
                       "gendermale" = "Male",
                       "gradeG3" = "Grade G3",
                       "gradeG4" = "Grade G4",
                       "MGMTUnmethylated" = "MGMT unmethylated"),
         term = factor(term,
                       levels = c(
                         "MGMT unmethylated",
                         "Grade G4",
                         "Grade G3",
                         "Male",
                         "Age",
                         "H2BC11",
                         "CD86",
                         "cg01878435")),
         #term <- forcats::fct_rev(term),
         HR_CI = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high),
         p.value = sprintf("%.3f", p.value))

# Create the forest plot
# not show G4, since therea re only 3 samples and yield a large CI range
hr_breaks <- c(0, 0.5, 1, 2, 5) # choose ticks you like
# Compute rightward offset dynamically (to keep space but not stretch axis)
pval_offset <- max(results[-7,]$conf.high, na.rm = TRUE) * 1.25

p_forest <- ggplot(results[-7,], aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange(shape=15) +  # Points and error bars
  #geom_point(shape=15) + 
  #geom_errorbar() +  # Points and error bars
  geom_text(aes(label = HR_CI), vjust = -1, size = 3.5) +  # Add HR labels
  # # p-value column aligned to right side
  geom_text(
    aes(x = term, y = max(conf.high, na.rm = TRUE) * 1.2,  # position to the right
        label = ifelse(p.value < 0.001, "<0.001", p.value)),
    hjust = 0, size = 3.5
  ) +
  # Add column title "p-value" at top
  annotate("text", x = 7.5, y = pval_offset, 
           label = "p-value", hjust = 0.3, vjust = -1, fontface = "bold", size = 4) +
  coord_flip(clip = "off") +  # Flip coordinates for a vertical forest plot
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  # Reference line at HR = 1
  scale_y_log10(breaks = hr_breaks, labels = hr_breaks) +
  #scale_y_continuous(breaks = hr_breaks, labels = hr_breaks) +
  labs(
    title = "",
    x = "",
    y = "Hazard Ratio (HR), 95% CI"
  ) +
  theme_classic() +
  theme(
    axis.ticks.y= element_blank(),
    axis.line.y.left = element_blank(),
    axis.title.y= element_blank(),
    axis.title = element_text(size = 12),
    axis.text   = element_text(size = 10),
    plot.margin = margin(12, 50, 12, 5))  # ← increased right margin for p-values


##################################
# fig 5c. KM for RPS; stratify
##################################
# KM for PRS within radition vs. non radition
p_km_rad     <- km_mt$ggsurv$Radiated$plot + ggtitle("Radiated")
p_km_nonrad  <- km_mt$ggsurv$`Non-radiated`$plot + ggtitle("Non-radiated")

p_km <- ggarrange(p_km_rad, p_km_nonrad, 
                  ncol = 2, nrow = 1)


##################################
# fig 5e. ROC: stratify in radiation
##################################
p_roc_mt <- roc_panel_aligned_stratify(km_mt, times = c(3, 5, 10) * 365, 
                                       legend_xy = c(0.2, 0))
p_timeROC <- ggarrange(p_roc_mt$Radiated, p_roc_mt$`Non-radiated`,
                       ncol = 2, nrow = 1)


