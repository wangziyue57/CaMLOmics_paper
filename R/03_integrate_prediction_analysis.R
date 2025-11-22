# run lten group lasso cox-ph model on stratified sample with MGMT methylation and IDH 
# other covariates: age, gender, grade


##################################
## load data
##################################
# omics data
load(file = paste0(RESULT, " filter_met_glioma_formediation_20250209.RData"))
load(file = paste0(RESULT, " filter_exp_glioma_formediation_20250209.RData"))
dat.exp.glioma.norm <- readRDS(file = paste0(DATA, "glioma_exp_clean_norm_20250209.rds"))
MET_gene_glioma <- read.csv(file = paste0(RESULT, "MET_gene_all_glioma_20250209.csv"))

# normalize and standardize (?) methylation and gene expression
# converting DNA methylation probes beta values to m_values
dat.met.glioma.mval.wt <- t(apply(dat.met.glioma.filter.wt, 1, function(x) log2(x/(1-x))))
dat.met.glioma.mval.scale.wt <- t(scale(t(dat.met.glioma.mval.wt)))

dat.met.glioma.mval.mt <- t(apply(dat.met.glioma.filter.mt, 1, function(x) log2(x/(1-x))))
dat.met.glioma.mval.scale.mt <- t(scale(t(dat.met.glioma.mval.mt)))

# normalized gene expression value
intersect(rownames(dat.exp.glioma.filter.wt),rownames(dat.exp.glioma.norm))
dat.exp.glioma.norm.wt <- dat.exp.glioma.norm[rownames(dat.exp.glioma.filter.wt), colnames(dat.exp.glioma.filter.wt)]
dat.exp.glioma.norm.scale.wt <- t(scale(t(log1p(dat.exp.glioma.norm.wt))))

dat.exp.glioma.norm.mt <- dat.exp.glioma.norm[rownames(dat.exp.glioma.filter.mt), colnames(dat.exp.glioma.filter.mt)]
dat.exp.glioma.norm.scale.mt <- t(scale(t(log1p(dat.exp.glioma.norm.mt))))


# clinical data
dat.clin.glioma.met <- readRDS(file = paste0(RESULT,"dat_met_clin_20250205.rds"))
dat.clin.glioma.exp <- readRDS(file = paste0(RESULT,"dat_exp_clin_20250205.rds"))

sample_id <- unique(c(intersect(dat.clin.glioma.met$sample, dat.clin.glioma.exp$sample)))
dat.clin.glioma.inte <- dat.clin.glioma.met %>%
  select(sample) %>%
  full_join(dat.clin.glioma.exp, by="sample") %>%
  filter(sample %in% sample_id)

# combine age-scale and radiation
dat.clin.glioma.inte$treatments[[2]][2,]$treatment_or_therapy=="no"
dat.clin.glioma.inte <- dat.clin.glioma.inte %>%
  mutate(radiation = sapply(dat.clin.glioma.inte$treatments, function(x) x[2,]$treatment_or_therapy=="no")) %>%
  mutate(age_scale = as.vector(scale(age_at_diagnosis))) 
rownames(dat.clin.glioma.inte) <- dat.clin.glioma.inte$sample

dat.clin.glioma.toUse <- dat.clin.glioma.inte %>%
  rename(IDH = paper_IDH.status,
         ch1p19q = paper_X1p.19q.codeletion,
         type = paper_IDH.codel.subtype,
         MGMT = paper_MGMT.promoter.status)


##################################
## train a latent group lasso Cox model
##################################
# latent group lasso survival function
source(paste0(CODE,"OverlapLasso.R"))

# stability selection for stable variables
run_lasso <- function(res_med, 
                      dat_met, # norm methylation dat
                      dat_exp, # norm expression dat
                      covariate, # vector with covariates name
                      dat_clin,
                      sampleID,
                      met_anno){
  
  # summary mediation res
  # normalize and standardize (?) methylation and gene expression
  # converting DNA methylation probes beta values to m_values
  met.med <- unique(res_med$res_med_sig$prob_id)
  dat_met <- dat_met[met.med,sampleID]
  
  # normalized gene expression value
  exp.med <- unique(res_med$res_med_sig$gene_id)
  dat_exp <- dat_exp[which(rownames(dat_exp) %in% exp.med),sampleID]
  
  #  run overlap group lasso
  dat.integrate <- rbind(dat_met, dat_exp)
  x <- cbind(age = dat_clin %>% filter(sample %in% sampleID) %>% pull(age_at_diagnosis), 
             gender = dat_clin %>% filter(sample %in% sampleID) %>% pull(gender) %>% as.factor(), 
             grade = dat_clin %>% filter(sample %in% sampleID) %>% pull(paper_Grade) %>% as.factor(),
             t(dat.integrate))
  y <- cbind(time = dat_clin %>% filter(sample %in% sampleID) %>% pull(overall_survival), 
             status = dat_clin %>% filter(sample %in% sampleID) %>% pull(status))
  
  # create index for each module
  group <- list()
  group[[1]] <- covariate #c("age","gender","grade")
  for (i in 2:(length(res_med$met.pathway)+1)) {
    group[[i]] <- res_med$met.pathway[[i-1]]
  }
  cvfit <- cv.grpsurvOverlap(X=x, y=y, group=group, nfolds=10, penalty = 'grLasso')
  pred_lp <- predict(cvfit, X = x, type = "link")
  
  # DNAm->mRNA prognostic pathways
  prog.pathway <- list()
  for(i in 1:length(names(coef(cvfit))[which(coef(cvfit) != 0)])){
    
    met.prog <- names(coef(cvfit))[which(coef(cvfit) != 0)][i]
    met.gene.prog <- met_anno %>%
      filter(probeID == met.prog) %>%
      pull(genesUniq)
    
    exp.prog <- res_med$res_med_sig %>%
      filter(prob_id == met.prog) %>%
      pull(gene_id)
    
    prog.pathway[[i]] <- c(met.prog, (paste(met.gene.prog, collapse = ", ")), (paste(exp.prog, collapse = ", ")))
  }
  prog.pathway <- t(do.call(cbind, prog.pathway)) %>%
    data.frame() %>%
    dplyr::rename(probe = X1, annotion.gene = X2, mediator.gene = X3)
  
  return(res_lasso = list(cvfit=cvfit, pred_lp=pred_lp, prog.pathway=prog.pathway, x=x, y=y))
}


##################################
## IDH WT + MGMT methy ##
## add age, gender, grade as covariate
##################################
res_med_wt_methy <- readRDS(file = paste0(RESULT, "res_mediation_glioma_IDHwt_methy_20250205.rds"))
sampleID <- dat.clin.glioma.toUse %>%
  filter(IDH == "WT" & MGMT == "Methylated") %>%
  pull(sample)

res_prog_wt_methy <- run_lasso(res_med = res_med_wt_methy,
                               dat_met=dat.met.glioma.mval.scale.wt, # norm methylation dat
                               dat_exp=dat.exp.glioma.norm.scale.wt, # norm expression dat
                               dat_clin=dat.clin.glioma.toUse,
                               covariate = c("paper_Grade", "age_scale", "gender"),
                               sampleID = sampleID,
                               met_anno=MET_gene_glioma)

# output mediation proportion for the prognostic pathways
res_prog_wt_methy_mp <- res_med_wt_methy$res_med_sig %>%
  filter(prob_id %in% res_prog_wt_methy$prog.pathway$probe) %>%
  left_join(res_prog_wt_methy$prog.pathway, join_by(prob_id == probe)) %>%
  dplyr::select(prob_id, annotion.gene, gene_id, alpha, beta, 'alpha_est*beta_est', DE) %>%
  rename(INE = 'alpha_est*beta_est') %>%group_by(prob_id) %>%
  mutate(total = sum(INE) + DE) %>%
  mutate(total_abs = sum(abs(INE)) + abs(DE)) %>%
  mutate(MP = INE/total) %>%
  mutate(MP_abs = abs(INE)/total_abs)


##################################
## IDH WT + MGMT unmethy ##
## add age, gender, grade as covariate
##################################
res_med_wt_unmethy <- readRDS(file = paste0(RESULT, "res_mediation_glioma_IDHwt_unmethy_20250205.rds"))
sampleID <- dat.clin.glioma.toUse %>%
  filter(IDH == "WT" & MGMT == "Unmethylated") %>%
  pull(sample)

res_prog_wt_unmethy <- run_lasso(res_med = res_med_wt_unmethy,
                                 dat_met=dat.met.glioma.mval.scale.wt, # norm methylation dat
                                 dat_exp=dat.exp.glioma.norm.scale.wt, # norm expression dat
                                 dat_clin=dat.clin.glioma.toUse,
                                 covariate = c("paper_Grade", "age_scale", "gender"),
                                 sampleID = sampleID,
                                 met_anno=MET_gene_glioma)

# output mediation proportion for the prognostic pathways
res_prog_wt_unmethy_mp <- res_med_wt_unmethy$res_med_sig %>%
  filter(prob_id %in% res_prog_wt_unmethy$prog.pathway$probe) %>%
  left_join(res_prog_wt_unmethy$prog.pathway, join_by(prob_id == probe)) %>%
  dplyr::select(prob_id, annotion.gene, gene_id, alpha, beta, 'alpha_est*beta_est', DE) %>%
  rename(INE = 'alpha_est*beta_est') %>%group_by(prob_id) %>%
  mutate(total = sum(INE) + DE) %>%
  mutate(total_abs = sum(abs(INE)) + abs(DE)) %>%
  mutate(MP = INE/total) %>%
  mutate(MP_abs = abs(INE)/total_abs)


##################################
## IDH MT ##
## add MGMT, ch1p19q, age, gender, grade as covariate
##################################
res_med_mt <- readRDS(file = paste0(RESULT, "res_mediation_glioma_IDHmt_20250205.rds"))
sampleID <- dat.clin.glioma.toUse %>%
  filter(IDH == "Mutant" & !is.na(MGMT) & !is.na(ch1p19q)) %>%
  pull(sample)

res_prog_mt <- run_lasso(res_med = res_med_mt,
                         dat_met=dat.met.glioma.mval.scale.mt, # norm methylation dat
                         dat_exp=dat.exp.glioma.norm.scale.mt, # norm expression dat
                         dat_clin=dat.clin.glioma.toUse,
                         covariate = c("paper_Grade", "age_scale", "gender","MGMT","ch1p19q"),
                         sampleID = sampleID,
                         met_anno=MET_gene_glioma)

# output mediation proportion for the prognostic pathways
res_prog_mt_mp <- res_med_mt$res_med_sig %>%
  filter(prob_id %in% res_prog_mt$prog.pathway$probe) %>%
  left_join(res_prog_mt$prog.pathway, join_by(prob_id == probe)) %>%
  dplyr::select(prob_id, annotion.gene, gene_id, alpha, beta, 'alpha_est*beta_est', DE) %>%
  rename(INE = 'alpha_est*beta_est') %>%group_by(prob_id) %>%
  mutate(total = sum(INE) + DE) %>%
  mutate(total_abs = sum(abs(INE)) + abs(DE)) %>%
  mutate(MP = INE/total) %>%
  mutate(MP_abs = abs(INE)/total_abs)

