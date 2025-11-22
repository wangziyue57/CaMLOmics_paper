
##################################
# input data
##################################
load(file = paste0(RESULT, "dat_integrate_clean_raw_20250210.RData"))

## cis mediation effect for all cpgs
res_cis_med_mt_mp <- read.csv(file = paste0(RESULT, "res_cis_mt_mp_20250221.csv"))
res_cis_med_wt_unmethy_mp <- read.csv(file = paste0(RESULT, "res_cis_wt_unmethy_mp_20250221.csv"))
res_cis_med_wt_methy_mp <- read.csv(file = paste0(RESULT, "res_cis_wt_methy_mp_filterAll_20250221.csv"))


## trans mediation effect for all cpgs
res_med_mt <- readRDS(file = paste0(RESULT, "res_mediation_glioma_IDHmt_20250205.rds"))
res_trans_med_mt_mp <- res_med_mt$res_med_sig %>%
  dplyr::select(prob_id, gene_id, alpha, beta, 'alpha_est*beta_est', DE, P_fdr_sobel) %>%
  rename(INE = 'alpha_est*beta_est', p_val = P_fdr_sobel) %>%
  group_by(prob_id) %>%
  mutate(total = sum(INE) + DE) %>%
  mutate(total_abs = sum(abs(INE)) + abs(DE)) %>%
  mutate(MP = INE/total) %>%
  mutate(MP_abs = abs(INE)/total_abs)

res_med_wt_unmethy <- readRDS(file = paste0(RESULT, "res_mediation_glioma_IDHwt_unmethy_20250205.rds"))
res_trans_med_wt_unmethy_mp <- res_med_wt_unmethy$res_med_sig %>%
  dplyr::select(prob_id, gene_id, alpha, beta, 'alpha_est*beta_est', DE, P_fdr_sobel) %>%
  rename(INE = 'alpha_est*beta_est', p_val = P_fdr_sobel) %>%
  group_by(prob_id) %>%
  mutate(total = sum(INE) + DE) %>%
  mutate(total_abs = sum(abs(INE)) + abs(DE)) %>%
  mutate(MP = INE/total) %>%
  mutate(MP_abs = abs(INE)/total_abs)

res_med_wt_methy <- readRDS(file = paste0(RESULT, "res_mediation_glioma_IDHwt_methy_20250125.rds")) ## contains CD276
res_trans_med_wt_methy_mp <- res_med_wt_methy$res_med_sig %>%
  dplyr::select(prob_id, gene_name, alpha, beta, 'alpha_est*beta_est', DE, P_fdr_sobel) %>%
  rename(INE = 'alpha_est*beta_est', 
         p_val = P_fdr_sobel,
         gene_id = gene_name) %>%
  group_by(prob_id) %>%
  mutate(total = sum(INE) + DE) %>%
  mutate(total_abs = sum(abs(INE)) + abs(DE)) %>%
  mutate(MP = INE/total) %>%
  mutate(MP_abs = abs(INE)/total_abs)


##################################
# helper function
##################################
# used for mediation results for all cpgs
generate_lollipop_plots <- function(data1, data2, top_n = 30) {
  # Step 1. Combine cis- and trans-
  scatter_cis <- data1 %>%
    mutate(source = "cis-") %>%
    rename(target_gene = gene_id)
  
  scatter_trans <- data2 %>%
    mutate(source = "trans-") %>%
    rename(target_gene = gene_id)
  
  scatter_df <- bind_rows(scatter_cis, scatter_trans) %>%
    filter(!is.na(target_gene) & target_gene != "") %>%   # removes NA genes
    distinct()
  
  # Step 2a. Counts (full, for summaries)
  count_cpg_full <- scatter_df %>%
    group_by(prob_id) %>%
    summarise(n_mRNA = n_distinct(target_gene), .groups = "drop")
  
  count_mrna_full <- scatter_df %>%
    group_by(target_gene) %>%
    summarise(n_cpg = n_distinct(prob_id), .groups = "drop")
  
  # Step 2b. Summaries (use all data)
  summary_cpg <- count_cpg_full %>%
    mutate(
      count_group = case_when(
        n_mRNA == 1 ~ "1 mRNA",
        n_mRNA == 2 ~ "2 mRNAs",
        n_mRNA > 2  ~ ">2 mRNAs"
      )
    ) %>%
    count(count_group, name = "n_prob_ids") %>%
    mutate(
      proportion = n_prob_ids / sum(n_prob_ids) * 100,
      count_group = factor(count_group, levels = c("1 mRNA", "2 mRNAs", ">2 mRNAs"))
    )
  
  summary_mrna <- count_mrna_full %>%
    mutate(
      count_group = case_when(
        n_cpg == 1 ~ "1 CpG",
        n_cpg == 2 ~ "2 CpGs",
        n_cpg > 2  ~ ">2 CpGs"
      )
    ) %>%
    count(count_group, name = "n_genes") %>%
    mutate(
      proportion = n_genes / sum(n_genes) * 100,
      count_group = factor(count_group, levels = c("1 CpG", "2 CpGs", ">2 CpGs"))
    )
  
  summary_type <- scatter_df %>%
    count(source, name = "n_pairs") %>%
    mutate(proportion = n_pairs / sum(n_pairs) * 100)
  
  # Step 3. Restrict to top N for plots
  count_cpg_top <- count_cpg_full %>%
    arrange(desc(n_mRNA)) %>%
    slice_head(n = top_n) %>%
    mutate(prob_id = factor(prob_id, levels = rev(prob_id)))
  
  count_mrna_top <- count_mrna_full %>%
    arrange(desc(n_cpg)) %>%
    slice_head(n = top_n) %>%
    mutate(target_gene = factor(target_gene, levels = rev(target_gene)))
  
  # --- Plots ---
  p_cpg <- ggplot(count_cpg_top, aes(x = n_mRNA, y = prob_id)) +
    geom_segment(aes(x = 0, xend = n_mRNA, y = prob_id, yend = prob_id),
                 color = "gray70") +
    geom_point(size = 3, color = "#1f77b4") +
    scale_x_continuous(expand = c(0,0)) +
    labs(x = "Number of mRNAs (mediators) per CpG", y = "CpG (top hits)") +
    theme_bw(base_size = 12) +
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
  
  p_mrna <- ggplot(count_mrna_top, aes(x = n_cpg, y = target_gene)) +
    geom_segment(aes(x = 0, xend = n_cpg, y = target_gene, yend = target_gene),
                 color = "gray70") +
    geom_point(size = 3, color = "#D95F02") +
    scale_x_continuous(expand = c(0,0)) +
    labs(x = "Number of CpGs mediated per gene", y = "mRNA (top hits)") +
    theme_bw(base_size = 12) +
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
  
  # --- Return ---
  list(
    scatter_df = scatter_df,
    cpg_lollipop  = p_cpg,
    mrna_lollipop = p_mrna,
    cpg_data      = count_cpg_full,
    mrna_data     = count_mrna_full,
    cpg_summary   = summary_cpg,
    mrna_summary  = summary_mrna,
    type_summary  = summary_type
  )
}



##################################
# save output
##################################
plots_mt <- generate_lollipop_plots(data1 = res_cis_med_mt_mp,
                                    data2 = res_trans_med_mt_mp,
                                    top_n = 30)

plots_wt_unmethy <- generate_lollipop_plots(data1 = res_cis_med_wt_unmethy_mp,
                                            data2 = res_trans_med_wt_unmethy_mp,
                                            top_n = 30)

plots_wt_methy <- generate_lollipop_plots(data1 = res_cis_med_wt_methy_mp,
                                          data2 = res_trans_med_wt_methy_mp,
                                          top_n = 30)

save(plots_mt, plots_wt_unmethy, plots_wt_methy,
     file = paste0(paper, "res_mediation_summary.RData"))


##################################
# fig 2a. lollipop plot of cpg vs. mrna
##################################
pA <-
  ggarrange((plots_mt$cpg_lollipop         + labs(subtitle = "IDH-MT")         +
               scale_x_continuous(limits = c(0, max(plots_mt$cpg_data$n_mRNA)+1), 
                                  #breaks = brk_cpg, 
                                  expand = c(0,0))),
            (plots_wt_unmethy$cpg_lollipop + labs(subtitle = "IDH-WT, unMGMT") +
               scale_x_continuous(limits = c(0, max(plots_wt_unmethy$cpg_data$n_mRNA)+1), 
                                  #breaks = brk_cpg, 
                                  expand = c(0,0))),
            (plots_wt_methy$cpg_lollipop   + labs(subtitle = "IDH-WT, mMGMT")   +
               scale_x_continuous(limits = c(0, max(plots_wt_methy$cpg_data$n_mRNA)+1), 
                                  #breaks = brk_cpg, 
                                  expand = c(0,0))),
            ncol = 1, nrow = 3,
            heights = c(1, 1, 1))


pB <-
  ggarrange((plots_mt$mrna_lollipop         + labs(subtitle = "IDH-MT")         +
               scale_x_continuous(limits = c(0, max(plots_mt$mrna_data$n_cpg)+1), 
                                  #breaks = brk_cpg, 
                                  expand = c(0,0))),
            (plots_wt_unmethy$mrna_lollipop + labs(subtitle = "IDH-WT, unMGMT") +
               scale_x_continuous(limits = c(0, max(plots_wt_unmethy$mrna_data$n_cpg)+1), 
                                  #breaks = brk_cpg, 
                                  expand = c(0,0))),
            (plots_wt_methy$mrna_lollipop   + labs(subtitle = "IDH-WT, mMGMT")   +
               scale_x_continuous(limits = c(0, max(plots_wt_methy$mrna_data$n_cpg)+1), 
                                  #breaks = brk_cpg, 
                                  expand = c(0,0))),
            ncol = 1, nrow = 3,
            heights = c(1, 1, 1))


##################################
# fig 2b. (Stacked) bar plot of number of mrnas (mediators) per cpg
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

# Combine CpG summaries
overall_cpg_summary <- bind_rows(
  plots_wt_methy$cpg_summary,
  plots_wt_unmethy$cpg_summary,
  plots_mt$cpg_summary
) %>%
  group_by(count_group) %>%
  summarise(n_prob_ids = sum(n_prob_ids), .groups = "drop") %>%
  mutate(proportion = n_prob_ids / sum(n_prob_ids) * 100)

p_cpg_summary <- ggplot(overall_cpg_summary,
                        aes(x = count_group, y = proportion, fill = count_group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(round(proportion, 1), "%")),
            vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("1 mRNA" = "#1b9e77",
                               "2 mRNAs" = "#d95f02",
                               ">2 mRNAs" = "#7570b3")) +
  labs(y = "Proportion of CpGs (%)", x = NULL,
       title = "") +
  ylim(0, 100) +
  fig_theme+
  theme(legend.position = "none")


##################################
# fig 2c. (Stacked) bar plot of number of cpgs per mrna
##################################
# Combine mRNA summaries
overall_mrna_summary <- bind_rows(
  plots_wt_methy$mrna_summary,
  plots_wt_unmethy$mrna_summary,
  plots_mt$mrna_summary
) %>%
  group_by(count_group) %>%
  summarise(n_genes = sum(n_genes), .groups = "drop") %>%
  mutate(proportion = n_genes / sum(n_genes) * 100)

p_mrna_summary <- ggplot(overall_mrna_summary, 
                         aes(x = count_group, y = proportion, fill = count_group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(round(proportion, 1), "%")),
            vjust = -0.5, size = 5, color = "black") +
  scale_fill_manual(values = c("1 CpG" = "#1b9e77",
                               "2 CpGs" = "#d95f02",
                               ">2 CpGs" = "#7570b3")) +
  labs(x = NULL, y = "Proportion of mRNAs (%)",
       title = "") +
  ylim(0, 100) +
  fig_theme +
  theme(legend.position = "none")


##################################
# fig 2d. (Stacked) bar plot of number of cis vs trans
##################################
# Combine cis/trans summaries
overall_type_summary <- bind_rows(
  plots_wt_methy$type_summary,
  plots_wt_unmethy$type_summary,
  plots_mt$type_summary
) %>%
  group_by(source) %>%
  summarise(n_pairs = sum(n_pairs), .groups = "drop") %>%
  mutate(proportion = n_pairs / sum(n_pairs) * 100)

p_type <- ggplot(overall_type_summary, aes(x = source, y = proportion, fill = source)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(round(proportion, 1), "%")),
            vjust = -0.5, size = 5, color = "black") +
  scale_fill_manual(values = c("cis-" = "#1b9e77", "trans-" = "#d95f02")) +
  labs(x = NULL, y = "Proportion of CpG–mRNA pairs (%)",
       title = "") +
  ylim(0, 100) +
  fig_theme +
  theme(legend.position = "none")


##################################
# fig 2e. (Stacked) bar plot of mediation effect
# maybe 4 levels (+/- for alpha and beta)
# or 2 levels (+/- for alpha*beta)
##################################
# Categorize based on alpha and beta signs
overall_mediation_effect <- bind_rows(plots_mt$scatter_df, 
                                      plots_wt_methy$scatter_df,
                                      plots_wt_unmethy$scatter_df) %>%
  distinct()

summary_effect <- overall_mediation_effect %>%
  mutate(effect_category = case_when(
    alpha > 0 & beta > 0 ~ "alpha>0, beta>0",
    alpha < 0 & beta < 0 ~ "alpha<0, beta<0",
    alpha > 0 & beta < 0 ~ "alpha>0, beta<0",
    alpha < 0 & beta > 0 ~ "alpha<0, beta>0",
    TRUE ~ "Other"
  )) %>%
  count(effect_category, name = "n_pairs") %>%
  mutate(proportion = n_pairs / sum(n_pairs) * 100,
         effect_category = factor(effect_category,
                                  levels = c("alpha>0, beta>0",
                                             "alpha<0, beta<0",
                                             "alpha>0, beta<0",
                                             "alpha<0, beta>0")))

# Bar plot of proportions
p_effect <- ggplot(summary_effect,
                   aes(x = effect_category, y = proportion, fill = effect_category)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(round(proportion, 1), "%")),
            vjust = -0.5, size = 5, color = "black") +
  scale_fill_manual(values = c("alpha>0, beta>0" = "#1b9e77",
                               "alpha<0, beta<0" = "#d95f02",
                               "alpha>0, beta<0" = "#7570b3",
                               "alpha<0, beta>0" = "#e7298a")) +
  labs(x = NULL, y = "Proportion of effect direction (%)",
       title = "") +
  ylim(0, 100) +
  fig_theme +
  theme(legend.position = "none")



##################################
# combine and output figure 2
##################################
# --- Panel C–F: overall summaries (already built above) ---
pC <- p_cpg_summary + ggtitle("Summary of mediators")
pD <- p_mrna_summary + ggtitle("Summary of CpGs")
pE <- p_type + ggtitle("Summary of cis-/trans- relation")
pF <- p_effect + ggtitle("Summary of mediation effect")

# --- Final layout ---
final_plot <-
  (pA | pB) /
  (pC | pD | pE | pF) +
  plot_layout(heights = c(4, 1)) +   # tweak if you want the first row shorter/ taller
  plot_annotation(tag_levels = "A") &   # <- tags A–F applied once, to the 6 top-level panels
  theme(plot.tag = element_text(size = 14, face = "bold"))
final_plot


