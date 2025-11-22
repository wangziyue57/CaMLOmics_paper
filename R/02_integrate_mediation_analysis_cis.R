# run cis mediation analysis on stratified sample with MGMT methylation and IDH 
# other covariates: age, gender, grade

# IDH WT, MGMT methy
# IDH WT, MGMT unmethy
# IDH MT, MGMT as covariate

# cox ~ probe + cis-gene
# cis-gene ~ probe



##################################
## load data
##################################
# omics data
load(file = paste0(RESULT, " filter_met_glioma_formediation_20250209.RData"))
load(file = paste0(RESULT, " filter_exp_glioma_formediation_20250209.RData"))
dat.exp.glioma.norm <- readRDS(file = paste0(DATA, "glioma_exp_clean_norm_20250209.rds"))
MET_gene_glioma <- read.csv(file = paste0(RESULT, "MET_gene_all_glioma_20250209.csv")) ## top 10% variable cpgs

# normalize and standardize (?) methylation and gene expression
# converting DNA methylation probes beta values to m_values
dat.met.glioma.mval.wt <- t(apply(dat.met.glioma.filter.wt, 1, function(x) log2(x/(1-x))))
dat.met.glioma.mval.scale.wt <- t(scale(t(dat.met.glioma.mval.wt)))

dat.met.glioma.mval.mt <- t(apply(dat.met.glioma.filter.mt, 1, function(x) log2(x/(1-x))))
dat.met.glioma.mval.scale.mt <- t(scale(t(dat.met.glioma.mval.mt)))

# normalized gene expression value
intersect(rownames(dat.exp.glioma.filter.wt),rownames(dat.exp.glioma.norm))
dat.exp.glioma.norm.wt <- dat.exp.glioma.norm[, colnames(dat.exp.glioma.filter.wt)]
dat.exp.glioma.norm.scale.wt <- t(scale(t(log1p(dat.exp.glioma.norm.wt))))

dat.exp.glioma.norm.mt <- dat.exp.glioma.norm[, colnames(dat.exp.glioma.filter.mt)]
dat.exp.glioma.norm.scale.mt <- t(scale(t(log1p(dat.exp.glioma.norm.mt))))

# intersect(colnames(dat.met.glioma.filter.mt), colnames(dat.exp.glioma.filter.mt))
# intersect(colnames(dat.met.glioma.filter.wt), colnames(dat.exp.glioma.filter.wt))

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
## function to run cis mediation analysis and get cis DNAm-mRNA group
##################################
# estimation of alpha M~X
sis_alpha <- function(X, M, COV, p){
  s_alpha <- function(j){
    if (is.null(COV)) {
      MX <- data.frame(M = M[, j], X = X)
    } else {
      MX <- data.frame(M = M[, j], X = X, COV = COV)
    }
    fit <- glm(M ~., data = MX)
    s1 <- summary(fit)$cov.scaled[2,2]   #var for alpha
    s2 <- summary(fit)$coef[2]           #coefficients for alpha
    s3 <- summary(fit)$coef[2,4]         #p-value for alpha
    return(data.frame(s_var=s1, s_coef=s2, s_p=s3))
  }
  dat=data.frame(do.call(rbind, lapply(1:ncol(M), s_alpha)))
  alpha_sis <- t(dat)
  colnames(alpha_sis) = colnames(M)
  return(s_alpha=alpha_sis)
}
# cox mediation Y~X+M
hmcox <- function(X, Y, M, COV){
  n <- nrow(M)
  p <- ncol(M)
  
  # estimator of alpha
  alpha_s <- sis_alpha(X, M, COV, p)
  alpha_est <- alpha_s[2, ]   
  var_alpha <- alpha_s[1, ]
  P_alpha <- alpha_s[3, ]
  
  if (is.null(COV)) {
    YMX <- data.frame(Y = Y, M, X = X)
  } else {
    YMX <- data.frame(Y = Y, M, X = X, COV = COV)
  }
  cox_model <- survival::coxph(Y ~ ., data = YMX)
  
  # direct effect of X->Y given M and COV
  DE <- summary(cox_model)$coefficients[(ncol(M)+1), 1]
  
  # estimator of beta
  beta_est <- summary(cox_model)$coefficients[1: ncol(M)]     
  var_beta <- (summary(cox_model)$coefficients[1:ncol(M),3])^2
  P_beta <- summary(cox_model)$coefficients[1: ncol(M), 5]  
  
  # indirect effect = the estimator of alpha*beta
  ab_est <- alpha_est * beta_est   
  
  # var(alpha*beta)
  var_ab <- (alpha_est^2) * var_beta + var_alpha * (beta_est^2) + var_alpha * var_beta
  
  # confidence interval
  conf_low <- ab_est - 1.96 * sqrt(var_ab)
  conf_up <- ab_est + 1.96 * sqrt(var_ab)
  
  # sobel test for alpha*beta
  s.test <- abs(ab_est)/sqrt(var_ab)   #z-score for sobel test
  P_sobel <- 2 * (1-pnorm(s.test))     #p-value of sobel test
  
  # mediation proportion
  TE <- DE + ab_est
  MP = ab_est/TE
  
  results <- data.frame(alpha = alpha_est, beta = beta_est,
                        `alpha_est*beta_est` = ab_est, 
                        conf_low=conf_low, conf_up=conf_up,
                        p_val=P_sobel,
                        var_ab=var_ab, var_alpha=var_alpha, var_beta=var_beta,
                        DE, TE, MP, check.names = FALSE)
  
  return(results)
}
## mediation function
run_med_cis <- function(X = X, 
                    Y = Y, 
                    M = M, 
                    COV = Z,
                    met_anno = met_anno){
  
  fitMed <- list()
  for (i in 1:ncol(X)) {
    gene <- met_anno %>%
      filter(probeID == colnames(X)[i]) %>%
      pull(genesUniq)
      
    fitMed[[i]] <- hmcox(X = X[,i], 
                         Y = Y, 
                         M = M[,gene, drop=F], 
                         COV = Z)
    # combine gene name and probe id
    fitMed[[i]] <- fitMed[[i]] %>%
      mutate(prob_id = rep(colnames(X)[i], nrow(fitMed[[i]]))) %>%
      mutate(gene_id = gene)
  }
  res_med <- do.call(rbind, fitMed)
  
  # sig mediated results with p<0.05
  res_med_sig <- res_med %>%
    filter(p_val < 0.05)
  
  # build DNAm->mRNA pathways from mediation results
  met.pathway <- list()
  for(i in 1:length(unique(res_med_sig$prob_id))){
    met.med <- unique(res_med_sig$prob_id)[i]
    exp.med <- res_med_sig %>%
      filter(prob_id == met.med) %>%
      pull(gene_id)
    met.pathway[[i]] <- c(met.med, exp.med)
  }
  
  return(list(res_med=res_med, res_med_sig=res_med_sig, met.pathway=met.pathway))
}


##################################
## IDH WT##
## add age, gender, grade as covariate
##################################
sampleID <- dat.clin.glioma.toUse %>%
  filter(IDH == "WT") %>%
  pull(sample)

# prepare DNAm and mRNA
# keep probe with cis-mRNA, and the cis-mRNA has expression value
met_gene_wt <- MET_gene_glioma %>%
  filter(probeID %in% rownames(dat.met.glioma.mval.scale.wt)) %>%
  filter(!is.na(genesUniq)) %>%
  filter(genesUniq %in% rownames(dat.exp.glioma.norm.scale.wt))

X <- t(dat.met.glioma.mval.scale.wt[unique(met_gene_wt$probeID),sampleID]) # exposures
M <- t(dat.exp.glioma.norm.scale.wt[unique(met_gene_wt$genesUniq),sampleID]) # mediators
Z <- dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% 
  dplyr::select(c("paper_Grade", "age_scale", "gender")) # covariates
Y = Surv(dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% pull(overall_survival), 
         dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% pull(status))
res_cis_med_wt <- run_med_cis(X = X, Y = Y, M = M, COV = Z, met_anno = met_gene_wt)

res_cis_wt_mp <- res_cis_med_wt$res_med_sig %>%
  dplyr::select(prob_id, gene_id, alpha, beta, 'alpha_est*beta_est', DE, p_val) %>%
  rename(INE = 'alpha_est*beta_est') %>% 
  group_by(prob_id) %>%
  mutate(total = sum(INE) + DE) %>%
  mutate(total_abs = sum(abs(INE)) + abs(DE)) %>%
  mutate(MP = INE/total) %>%
  mutate(MP_abs = abs(INE)/total_abs)



##################################
## IDH MT ##
## add MGMT, age, gender, grade as covariate
##################################
sampleID <- dat.clin.glioma.toUse %>%
  filter(IDH == "Mutant" & !is.na(MGMT)) %>%
  pull(sample)

# prepare DNAm and mRNA
# keep probe with cis-mRNA, and the cis-mRNA has expression value
met_gene_mt <- MET_gene_glioma %>%
  filter(probeID %in% rownames(dat.met.glioma.mval.scale.mt)) %>%
  filter(!is.na(genesUniq)) %>%
  filter(genesUniq %in% rownames(dat.exp.glioma.norm.scale.mt))

X <- t(dat.met.glioma.mval.scale.mt[unique(met_gene_mt$probeID),sampleID]) # exposures
M <- t(dat.exp.glioma.norm.scale.mt[unique(met_gene_mt$genesUniq),sampleID]) # mediators
Z <- dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% 
  dplyr::select(c("paper_Grade", "age_scale", "gender","MGMT")) # covariates
Y = Surv(dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% pull(overall_survival), 
         dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% pull(status))
res_cis_med_mt <- run_med_cis(X = X, Y = Y, M = M, COV = Z, met_anno = met_gene_mt)

res_cis_mt_mp <- res_cis_med_mt$res_med_sig %>%
  dplyr::select(prob_id, gene_id, alpha, beta, 'alpha_est*beta_est', DE, p_val) %>%
  rename(INE = 'alpha_est*beta_est') %>% 
  group_by(prob_id) %>%
  mutate(total = sum(INE) + DE) %>%
  mutate(total_abs = sum(abs(INE)) + abs(DE)) %>%
  mutate(MP = INE/total) %>%
  mutate(MP_abs = abs(INE)/total_abs)


##################################
## IDH WT + MGMT unmethy ##
## add codeletion age, gender, grade as covariate
##################################
sampleID <- dat.clin.glioma.toUse %>%
  filter(IDH == "WT" & MGMT == "Unmethylated") %>%
  pull(sample)

# prepare DNAm and mRNA
# keep probe with cis-mRNA, and the cis-mRNA has expression value
met_gene_wt <- MET_gene_glioma %>%
  filter(probeID %in% rownames(dat.met.glioma.mval.scale.wt)) %>%
  filter(!is.na(genesUniq)) %>%
  filter(genesUniq %in% rownames(dat.exp.glioma.norm.scale.wt))

X <- t(dat.met.glioma.mval.scale.wt[unique(met_gene_wt$probeID),sampleID]) # exposures
M <- t(dat.exp.glioma.norm.scale.wt[unique(met_gene_wt$genesUniq),sampleID]) # mediators
Z <- dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% 
  dplyr::select(c("paper_Grade", "age_scale", "gender")) # covariates
Y = Surv(dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% pull(overall_survival), 
         dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% pull(status))
res_cis_med_wt_unmeth <- run_med_cis(X = X, Y = Y, M = M, COV = Z, met_anno = met_gene_wt)

res_cis_wt_unmethy_mp <- res_cis_med_wt_unmeth$res_med_sig %>%
  dplyr::select(prob_id, gene_id, alpha, beta, 'alpha_est*beta_est', DE, p_val) %>%
  rename(INE = 'alpha_est*beta_est') %>% 
  group_by(prob_id) %>%
  mutate(total = sum(INE) + DE) %>%
  mutate(total_abs = sum(abs(INE)) + abs(DE)) %>%
  mutate(MP = INE/total) %>%
  mutate(MP_abs = abs(INE)/total_abs)


##################################
## IDH WT + MGMT methy ##
## add age, gender, grade as covariate
## from whole sample filtering probe results
##################################
sampleID <- dat.clin.glioma.toUse %>%
  filter(IDH == "WT" & MGMT == "Methylated") %>%
  pull(sample)

# prepare DNAm and mRNA
# keep probe with cis-mRNA, and the cis-mRNA has expression value
met_gene_wt <- MET_gene_glioma %>%
  filter(probeID %in% rownames(dat.met.glioma.toUse.mval.scale)) %>%
  filter(!is.na(genesUniq)) %>%
  filter(genesUniq %in% rownames(dat.exp.glioma.norm.scale.wt))

X <- t(dat.met.glioma.toUse.mval.scale[unique(met_gene_wt$probeID),sampleID]) # exposures
M <- t(dat.exp.glioma.norm.scale.wt[unique(met_gene_wt$genesUniq),sampleID]) # mediators
Z <- dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% 
  select(c("paper_Grade", "age_scale", "gender")) # covariates
Y = Surv(dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% pull(overall_survival), 
         dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% pull(status))
res_cis_med_wt_meth <- run_med_cis(X = X, Y = Y, M = M, COV = Z, met_anno = met_gene_wt)

res_cis_wt_methy_mp <- res_cis_med_wt_meth$res_med_sig %>%
  dplyr::select(prob_id, gene_id, alpha, beta, 'alpha_est*beta_est', DE, p_val) %>%
  rename(INE = 'alpha_est*beta_est') %>% 
  group_by(prob_id) %>%
  mutate(total = sum(INE) + DE) %>%
  mutate(total_abs = sum(abs(INE)) + abs(DE)) %>%
  mutate(MP = INE/total) %>%
  mutate(MP_abs = abs(INE)/total_abs)
