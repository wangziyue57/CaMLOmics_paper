
source(paste0(CODE,"visualize_func_v2_20250910.R"))


##################################
# input data
##################################
# summary of mediation results for the identified prognostic biomarkers
env1 <- new.env()
load(file = paste0(paper, "res_prognosis_mediation_summary.RData"), envir = env1)
env2 <- new.env()
load(file = paste0(paper, "res_mediation_summary.RData"), envir = env2)

# pathway results from metascape
res_enrich_mt <- read_excel(paste0(RESULT, "metascape_result_mt_prog.xlsx"), sheet = "Enrichment")
res_enrich_wt_unmethy <- read_excel(paste0(RESULT, "metascape_result_wt_unmethy_prog.xlsx"), sheet = "Enrichment")
res_enrich_wt_methy <- read_excel(paste0(RESULT, "metascape_result_wt_methy_prog.xlsx"), sheet = "Enrichment")

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
# fig 4a. barplot of top enrichment results
##################################
plot_enrichment_bar <- function(df, top_n = 20,
                                y_col = "Description",
                                q_col = "LogP",
                                count_col = "InTerm_InList",
                                group_col = "GroupID",
                                truncate_width = 60,
                                glue_term_desc = TRUE) {
  
  df_plot <- df %>%
    # keep only “Summary” rows
    filter(grepl("Summary", .data[[group_col]])) %>%
    # build Count from “x/y”
    mutate(Count = ifelse(grepl("/", .data[[count_col]]),
                          as.numeric(sub("/.*", "", .data[[count_col]])),
                          NA_real_)) %>%
    # NEW: make a combined label “TERM: Description”
    mutate(
      .term_desc = if (glue_term_desc) paste0(.data[["Term"]], ": ", .data[["Description"]])
      else .data[[y_col]]
    ) %>%
    # rank by significance (more negative logP means stronger ⇒ sort ascending LogP)
    arrange(.data[[q_col]]) %>%
    slice_head(n = top_n) %>%
    # truncate long labels for cleaner panels
    mutate(.term_desc = ifelse(nchar(.term_desc) > truncate_width,
                               paste0(substr(.term_desc, 1, truncate_width), "..."),
                               .term_desc))
  
  ggplot(
    df_plot,
    aes(x = Count,
        y = reorder(.term_desc, .data[[q_col]], decreasing = TRUE),
        fill = .data[[q_col]])
  ) +
    geom_col(width = 0.7) +
    scale_fill_gradientn(
      colours = c("red", "purple", "blue"),
      name = expression(-log[10]~p)  # adjust to q if you're using Log(q-value)
    ) +
    theme_bw(base_size = 12) +
    labs(x = "Gene count", y = NULL) +
    theme(
      plot.title   = element_text(size = 12),
      axis.title = element_text(size = 12),
      axis.text.y  = element_text(size = 10, hjust = 1),
      axis.text.x  = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )
}


# Generate plots
p1 <- plot_enrichment_bar(res_enrich_mt, top_n = 15) + ggtitle("IDH-MT")
p2 <- plot_enrichment_bar(res_enrich_wt_unmethy, top_n = 15) + ggtitle("IDH-WT, unMGMT")
p3 <- plot_enrichment_bar(res_enrich_wt_methy, top_n = 15) + ggtitle("IDH-WT, mMGMT")

p_enrich <-  ggarrange(p1,p2,p3,
          ncol = 1, nrow = 3, align = "v",heights = c(2,1,1),
          legend = "right", common.legend = FALSE)


##################################
# fig 4b. 2D hub map: Scatterplot (x = connectivity, y = pathway frequency). 
# Label top genes; dashed lines mark medians to show “hub quadrant”.
##################################
# frequency of pathways
gene_freq_mt <- res_enrich_mt %>%
  # 1) split each “A/B/C…” geneID string into its own rows
  separate_rows(Symbols, sep = ",") %>%
  # 2) trim any stray whitespace
  mutate(Symbols = str_trim(Symbols)) %>%
  # 3) count how often each gene appears
  count(Symbols, name = "frequency") %>%
  # 4) sort descending
  arrange(desc(frequency))

gene_freq_wt_unmethy <- res_enrich_wt_unmethy %>%
  # 1) split each “A/B/C…” geneID string into its own rows
  separate_rows(Symbols, sep = ",") %>%
  # 2) trim any stray whitespace
  mutate(Symbols = str_trim(Symbols)) %>%
  # 3) count how often each gene appears
  count(Symbols, name = "frequency") %>%
  # 4) sort descending
  arrange(desc(frequency))

gene_freq_wt_methy <- res_enrich_wt_methy %>%
  # 1) split each “A/B/C…” geneID string into its own rows
  separate_rows(Symbols, sep = ",") %>%
  # 2) trim any stray whitespace
  mutate(Symbols = str_trim(Symbols)) %>%
  # 3) count how often each gene appears
  count(Symbols, name = "frequency") %>%
  # 4) sort descending
  arrange(desc(frequency))


# frequency of combination of meditors and pathways
gene_freq_all_mt <- gene_freq_mt %>%
  left_join(env2$plots_mt$mrna_data, by = join_by(Symbols == target_gene)) %>%
  mutate(n_cpg = ifelse(is.na(n_cpg), 1, n_cpg))

gene_freq_all_wt_unmethy <- gene_freq_wt_unmethy %>%
  left_join(env2$plots_wt_unmethy$mrna_data, by = join_by(Symbols == target_gene)) %>%
  mutate(n_cpg = ifelse(is.na(n_cpg), 1, n_cpg))

gene_freq_all_wt_methy <- gene_freq_wt_methy %>%
  left_join(env2$plots_wt_methy$mrna_data, by = join_by(Symbols == target_gene)) %>%
  mutate(n_cpg = ifelse(is.na(n_cpg), 1, n_cpg))


# ---- Function for 2D hub gene map ----
make_hubscore_plot <- function(df, subgroup_name,
                               top_n = 20,
                               plot_title = subgroup_name) {
  
  # calculate z-scores & hub score
  df_ranked <- df %>%
    mutate(
      z_freq = scale(frequency, center = TRUE, scale = TRUE)[,1],
      z_cpg  = scale(n_cpg, center = TRUE, scale = TRUE)[,1],
      hub_score = z_freq + z_cpg,
      subgroup = subgroup_name
    ) %>%
    arrange(desc(hub_score)) %>%
    mutate(rank = row_number())
  
  # scatter plot
  p <- ggplot(df_ranked, aes(x = n_cpg, y = frequency)) +
    geom_point(aes(size = hub_score, color = hub_score)) +
    ggrepel::geom_text_repel(aes(label = ifelse(rank <= top_n, Symbols, "")),
                             size = 3, max.overlaps = 15) +
    scale_color_viridis_c() +
    labs(x = "CpG–mRNA connectivity (#CpGs)",
         y = "Pathway frequency (#pathways)",
         color = "Hub score",
         size = "Hub score",
         title = plot_title) +
    theme_classic()
  
  return(list(ranked_table = df_ranked, plot = p))
}

# Run for three subgroups
res_mt <- make_hubscore_plot(gene_freq_all_mt,
                             subgroup_name = "IDH-MT")
res_wt_mMGMT <- make_hubscore_plot(gene_freq_all_wt_methy,
                                   subgroup_name = "IDH-WT, mMGMT")
res_wt_unMGMT <- make_hubscore_plot(gene_freq_all_wt_unmethy,
                                    subgroup_name = "IDH-WT, unMGMT")

# Combine 2D gene map side by side
p_hub_map <- ggarrange(res_mt$plot, res_wt_unMGMT$plot, res_wt_mMGMT$plot, 
                      ncol = 3, nrow = 1,
                      legend = "top", common.legend = TRUE)


##################################
# run cox on top hub genes (or all genes)
##################################
# helper function
fig_theme <- theme_classic() +
  theme(
    plot.title   = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    #axis.ticks = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10)
  )
plot_km_roc <- function(dat,
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
  s <- summary(cox_model_PRS)
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
    pval = FALSE,
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
  ggsurv$plot <- ggsurv$plot +
    annotate("text", x = 100, y = 0.1,
             label = paste0("HR (95% CI): ", HR_CI, "\n",
                            "Log-rank: ", surv_pvalue(sfit)$pval.txt),
             hjust = 0, size = 4)
  
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

# IDH-mt
biomarkers <- res_mt$ranked_table[c(1:20),]$Symbols
biomarkers <- sub("-", ".", biomarkers)
km_mt <- plot_km_roc(dat.integrate.scale[sampleID.3,], time_col = "overall_survival", event_col = "status",
                     predictors = biomarkers,
                     covariates = c("age_scale","gender","grade","MGMT"),
                     times      = c(3*365, 5*365, 10*365))         

# IDH-wt, mMGMT
biomarkers <- res_wt_mMGMT$ranked_table[c(1:20),]$Symbols
biomarkers <- sub("-", ".", biomarkers)[which(sub("-", ".", biomarkers) != "CHLSN")]
km_wt_methy <- plot_km_roc(dat.integrate.scale[sampleID.2,], time_col = "overall_survival", event_col = "status",
                           predictors = biomarkers,
                           covariates = c("age_scale","gender","grade"),
                           times      = c(1*365, 2*365, 3*365))

# IDH-wt, unMGMT
biomarkers <- res_wt_unMGMT$ranked_table[c(1:20),]$Symbols
biomarkers <- sub("-", ".", biomarkers)
km_wt_unmethy <- plot_km_roc(dat.integrate.scale[sampleID.1,], time_col = "overall_survival", event_col = "status",
                             predictors = biomarkers,
                             covariates = c("age_scale","gender","grade"),
                             times      = c(1*365, 2*365, 3*365))


##################################
# fig 4d. KM for refined PRS
##################################
# Extract the KM curves only (not risk tables)
p_km_mt        <- km_mt$ggsurv$plot + ggtitle("IDH-MT")
p_km_wt_methy  <- km_wt_methy$ggsurv$plot + ggtitle("IDH-WT, mMGMT")
p_km_wt_unmethy<- km_wt_unmethy$ggsurv$plot + ggtitle("IDH-WT, unMGMT") + coord_cartesian(xlim = c(0, 2000)) + scale_x_continuous(breaks = seq(0, 2000, 400))

# Combine side by side
final_km <- ggarrange(p_km_mt, p_km_wt_unmethy, p_km_wt_methy, 
                      ncol = 3, nrow = 1)


##################################
# fig 4e. ROC for PRS
##################################
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
                              cex_leg = 1) {  # ≈ pt
  
  ggplotify::as.ggplot(function() {
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    par(mar = mar, mgp = mgp, tcl = tcl,
        cex.axis = cex_axis, cex.lab = cex_lab,
        font.lab = 1, font.axis = 1, family = "sans")  # ⬅️ plain text, no bold
    
    roc <- res$roc_train_apparent
    ci  <- res$CI_train_apparent
    
    ## ---- first ROC curve ----
    plot(roc, time = times[1], col = cols[1],
         title = FALSE, xlab = "1 – Specificity", ylab = "Sensitivity",
         lwd = 2.2)
    
    ## ---- add more curves if present ----
    if (length(times) > 1) {
      for (i in 2:length(times)) {
        plot(roc, time = times[i], add = TRUE, col = cols[i],
             title = FALSE, lwd = 2.2)
      }
    }
    
    ## ---- legend text ----
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
    
    ## ---- title (plain) ----
    mtext(title, side = 3, adj = 0, line = 0.5, font = 1, cex = 1)
  }) +
    theme(plot.margin = margin(12, 5, 12, 5))
}

p_roc_mt <- roc_panel_aligned(
  km_mt, times = c(3, 5, 10) * 365, title = "IDH-MT",
  legend_xy = c(0.2, 0)  # nudge left/right as you like (0..1 scale)
)

p_roc_wt_methy <- roc_panel_aligned(
  km_wt_methy, times = c(1, 2, 3) * 365, title = "IDH-WT, mMGMT",
  legend_xy = c(0.2, 0)
)

p_roc_wt_unmethy <- roc_panel_aligned(
  km_wt_unmethy, times = c(1, 2, 3) * 365, title = "IDH-WT, unMGMT",
  legend_xy = c(0.2, 0)
)

# Combine side by side
final_roc <- ggarrange(p_roc_mt, p_roc_wt_unmethy, p_roc_wt_methy, 
                       ncol = 3, nrow = 1)

