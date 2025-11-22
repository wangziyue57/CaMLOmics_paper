
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
# run cox on top hub genes (or all genes)
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

# helper function
plot_km_roc_v2 <- function(dat,
                        time_col = "overall_survival",
                        event_col = "status",
                        predictors,
                        covariates,
                        times,
                        palette = c("#E64B35FF", "#1F78A0")) {
  
  # 1. Define risk groups by median split
  # Fit Cox model
  all_predictors <- c(predictors, covariates)
  cox_formula <- as.formula(
    paste0("Surv(", time_col, ", ", event_col, ") ~ ",
           paste(all_predictors, collapse = " + "))
  )
  cox_model <- coxph(cox_formula, data = dat, x = TRUE, y = TRUE)
  dat$risk_score <- predict(cox_model, newdata = dat, type = "lp")
  
  group_sample <- ifelse(
    dat$risk_score >= quantile(dat$risk_score, 0.5),
    "High risk", "Low risk"
  )
  
  # 2. get HR for PRS
  cox_model_PRS <- coxph(
    as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ dat$risk_score")),
    data = dat,
    x = TRUE, y = TRUE
  )
  s <- summary(cox_model)
  HR    <- exp(s$coefficients[, "coef"])
  lower <- s$conf.int[, "lower .95"]
  upper <- s$conf.int[, "upper .95"]
  HR_CI <- sprintf("%.2f (%.2f, %.2f)", HR, lower, upper)
  
  # 3. Kaplan–Meier fit
  sfit <- surv_fit(
    as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ group_sample")),
    data = dat
  )
  
  # 4. Plot with ggsurvplot
  ggsurv <- ggsurvplot(
    sfit,
    data = dat,
    conf.int = TRUE,
    risk.table = TRUE,
    pval = TRUE,
    xlab = "Time (days)",
    legend.labs = c("High", "Low"),
    legend.title = "PRS",
    surv.median.line = "hv",
    risk.table.y.text.col = TRUE,
    risk.table.y.text = FALSE,
    palette = palette,
    ggtheme = fig_theme
  )
  
  # 5. Add HR + log-rank annotation
  # ggsurv$plot <- ggsurv$plot +
  #   annotate("text", x = 200, y = 0.1,
  #            label = paste0("HR (95% CI):\n ", HR_CI, "\n",
  #                           "Log-rank:\n ", surv_pvalue(sfit)$pval.txt),
  #            hjust = 0, size = 4)
  
  # ROC curve with appearant AUC on training whole data
  roc_train_apparent <- timeROC(
    T      = dat[[time_col]],
    delta  = dat[[event_col]],
    marker = dat$risk_score,
    cause  = 1,
    weighting = "marginal",
    times  = times,
    ROC    = TRUE,
    iid    = TRUE
  )
  CI_train_apparent <- confint(roc_train_apparent)
  
  return(list(ggsurv=ggsurv,
              risk_group = group_sample,
              cox_model          = cox_model,
              cox_model_PRS = cox_model_PRS,
              roc_train_apparent = roc_train_apparent,
              CI_train_apparent  = CI_train_apparent))
}


roc_panel_aligned <- function(res,
                              times,                     
                              title,
                              cols = c("#1B9E77", "#D95F02", "#7570B3"),
                              mar  = c(5.0, 5.4, 2.2, 1.6),  # ⬅️ more top/bottom margin
                              mgp  = c(2.2, 0.6, 0),
                              tcl  = -0.3,
                              cex_axis = 0.8,   # ≈ pt
                              cex_lab  = 1,   # ≈ pt
                              legend_xy = c(0.55, 0.14),      
                              cex_leg = 1) {
  
  ggplotify::as.ggplot(function() {
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mar = mar, mgp = mgp, tcl = tcl,
        cex.axis = cex_axis, cex.lab = cex_lab,
        font.lab = 1, font.axis = 1, family = "sans")  # plain text
    
    roc <- res$roc_train_apparent
    ci  <- res$CI_train_apparent
    
    # plot first curve
    plot(roc, time = times[1], col = cols[1],
         title = FALSE, xlab = "1-Specificity", ylab = "Sensitivity",
         lwd = 2.2)
    
    if (length(times) > 1) {
      for (i in 2:length(times)) {
        plot(roc, time = times[i], add = TRUE, col = cols[i], title = FALSE,
             lwd = 2.2)
      }
    }
    
    ## ---- legend text (CI divided by 100) ----
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
    mtext(title, side = 3, adj = 0, line = 0.5, font = 1, cex = 1)
  }) +
    theme(plot.margin = margin(12, 5, 12, 5)) 
}



##################################
# cg14217861(NOTCH1);
# CD274
##################################
biomarkers <- c("cg14217861", "NOTCH1", "CD274")
km_cd274 <- plot_km_roc_v2(dat.integrate.scale[sampleID.1,], time_col = "overall_survival", event_col = "status",
                             predictors = biomarkers,
                             covariates = c("age_scale","gender","grade"),
                             times      = c(1*365, 2*365, 3*365))


##################################
# fig 6a. correlation plot
##################################
p_corr <- plot_corr_prob_mrna(dat.integrate.raw, 
                              probe = "cg14217861",anno_gene = "NOTCH1",med_gene = "CD274")
p_corr_all_cd274 <- ggarrange(plotlist = p_corr, ncol = 2,  legend = "top", common.legend=TRUE)
p_corr_all_cd274 <- ggarrange(plotlist = p_corr, nrow = 2,  legend = "top", common.legend=TRUE)


##################################
# fig 6b. forest plot
##################################
# Extract relevant statistics
model_summary <- summary(km_cd274$cox_model)
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
                       "gradeG4" = "Grade G4"),
         term = factor(term,
                       levels = c(
                         "Grade G4",
                         "Grade G3",
                         "Male",
                         "Age",
                         "CD274",
                         "NOTCH1",
                         "cg14217861")),
         #term <- forcats::fct_rev(term),
         HR_CI = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high),
         p.value = sprintf("%.3f", p.value))

# Create the forest plot
# not show G4, since therea re only 3 samples and yield a large CI range
hr_breaks <- c(0, 0.5, 1, 2, 5, 10) # choose ticks you like
# Compute rightward offset dynamically (to keep space but not stretch axis)
pval_offset <- max(results$conf.high, na.rm = TRUE) * 1.25

p_forest_cd274 <- ggplot(results, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
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
# fig 6c. KM for RPS
##################################
p_km_cd274 <- km_cd274$ggsurv$plot + ggtitle("IDH-WT, unMGMT") + coord_cartesian(xlim = c(0, 2000)) + scale_x_continuous(breaks = seq(0, 2000, 400))

##################################
# fig 6d. ROC
##################################
p_timeROC_cd274 <- roc_panel_aligned(
  km_cd274, times = c(1, 2, 3) * 365, title = "IDH-WT, unMGMT",
  legend_xy = c(0.2, 0)
)


##################################
# cg10586317(CD276);	NEK6
# identified in screening by whole sample
# train in WT, methy
# only 2 people in G2, yield a large HR for grade
# combine G2 and G3
##################################
train_subgroup <- dat.integrate.scale %>%
  filter(group == "WT, Methy") %>%
  mutate(grade = as.character(grade)) %>%
  mutate(grade = ifelse(grade == "G2", "G3", grade)) %>%
  mutate(grade = factor(grade)) %>%
  mutate(grade = relevel(grade, ref = "G3"))

biomarkers <- c("cg10586317", "CD276", "NEK6")
km_cd276 <- plot_km_roc_v2(train_subgroup, time_col = "overall_survival", event_col = "status",
                        predictors = biomarkers,
                        covariates = c("age_scale","gender","grade"),
                        times      = c(1*365, 2*365, 3*365))


##################################
# fig 6a. correlation plot
##################################
p_corr <- plot_corr_prob_mrna(dat.integrate.raw, 
                              probe = "cg10586317",anno_gene = "CD276",med_gene = "NEK6")
p_corr_all_cd276 <- ggarrange(plotlist = p_corr, nrow = 2,  legend = "top", common.legend=TRUE)
p_corr_all_cd276 <- ggarrange(plotlist = p_corr, ncol = 2,  legend = "top", common.legend=TRUE)


##################################
# fig 6b. forest plot
##################################
# Extract relevant statistics
model_summary <- summary(km_cd276$cox_model)
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
                       "gradeG4" = "Grade G4"),
         term = factor(term,
                       levels = c(
                         "Grade G4",
                         "Grade G3",
                         "Male",
                         "Age",
                         "NEK6",
                         "CD276",
                         "cg10586317")),
         #term <- forcats::fct_rev(term),
         HR_CI = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high),
         p.value = sprintf("%.3f", p.value))

# Create the forest plot
# not show G4, since therea re only 3 samples and yield a large CI range
hr_breaks <- c(0, 0.5, 1, 2, 5) # choose ticks you like
# Compute rightward offset dynamically (to keep space but not stretch axis)
pval_offset <- max(results$conf.high, na.rm = TRUE) * 1.25

p_forest_cd276 <- ggplot(results, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
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
  annotate("text", x = 6.5, y = pval_offset, 
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
# fig 6c. KM for RPS
##################################
p_km_cd276 <- km_cd276$ggsurv$plot + ggtitle("IDH-WT, mMGMT")


##################################
# fig 6d. ROC
##################################
p_timeROC_cd276 <- roc_panel_aligned(
  km_cd276, times = c(1, 2, 3) * 365, title = "IDH-WT, mMGMT",
  legend_xy = c(0.2, 0)
)


##################################
# fig 6e. boxplot plot of CD274 and cpg in high vs. low
##################################
p_boxplot_genes <- function(gene) {
  ggplot(dat, aes(x = risk_group, y = !!sym(gene), fill = risk_group)) +
    geom_boxplot(outlier.shape = NA) +
    labs(title = gene, y = "Expr value", fill = "Group") +
    scale_fill_manual(values = c("High risk" = "#d95f02",
                                 "Low risk" = "#7570b3")) +
    theme_classic()+
    theme(
      plot.title   = element_text(size = 12),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 10)
    )
}

dat <- data.frame(CD274 = dat.integrate.raw[sampleID.1,]$CD274,
                  risk_group = km_cd274$risk_group)
p_boxplot_cd274 <- p_boxplot_genes("CD274")

dat <- data.frame(CD276 = dat.integrate.raw[sampleID.2,]$CD276,
                  risk_group = km_cd276$risk_group)
p_boxplot_cd276 <- p_boxplot_genes("CD276")

p_boxplot <- ggarrange(p_boxplot_cd274, p_boxplot_cd276, nrow = 2,  legend = "right", common.legend=TRUE)
p_boxplot <- ggarrange(p_boxplot_cd274, p_boxplot_cd276, ncol = 2,  legend = "right", common.legend=TRUE, labels = "K")




