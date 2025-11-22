# rerun mediation analysis on stratified sample with MGMT methylation and IDH
# other covariates: age, gender, grade

# IDH WT, MGMT methy
# IDH WT, MGMT unmethy
# IDH MT, MGMT as covariate


##################################
## load data
##################################
# omics data
load(file = paste0(RESULT, " filter_met_glioma_formediation_20250209.RData"))
load(file = paste0(RESULT, " filter_exp_glioma_formediation_20250209.RData"))
dat.exp.glioma.norm <- readRDS(file = paste0(DATA, "glioma_exp_clean_norm_20250209.rds"))
MET_gene_glioma <- read.csv(file = paste0(RESULT, "MET_gene_all_glioma_20250209.csv"))

# normalize and standardize methylation and gene expression
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
## function to run mediation analysis and get DNAm-mRNA group
##################################
# HIMA survival function
source(paste0(CODE, "HIMAsurvival.R"))
run_med <- function(X = X, 
                    Y = Y, 
                    M = M, 
                    COV = Z){
  
  fitMed <- list()
  for (i in 1:ncol(X)) {
    fitMed[[i]] <- hmas(X = X[,i], 
                        Y = Y, 
                        M = M, 
                        COV = Z, 
                        penalty = "MCP",
                        path = 'both',
                        topN = NULL,
                        # topN = ncol(M),
                        verbose = TRUE)
    
    # combine gene name and probe id
    fitMed[[i]] <- fitMed[[i]][[3]] %>%
      rownames_to_column("gene_id") %>%
      mutate(prob_id = rep(colnames(X)[i], nrow(fitMed[[i]][[3]]))) 
  }
  res_med <- do.call(rbind, fitMed)
  
  # sig mediated results with FDR<0.05
  res_med_sig <- res_med %>%
    filter(P_fdr_sobel < 0.05)
  
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
## IDH WT + MGMT methy ##
## add age, gender, grade as covariate
##################################
sampleID <- dat.clin.glioma.toUse %>%
  filter(IDH == "WT" & MGMT == "Methylated") %>%
  pull(sample)

X <- t(dat.met.glioma.mval.scale.wt[,sampleID]) # exposures
M <- t(dat.exp.glioma.norm.scale.wt[,sampleID]) # mediators
Z <- dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% 
  select(c("paper_Grade", "age_scale", "gender")) # covariates
Y = Surv(dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% pull(overall_survival), 
         dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% pull(status))
res_med_wt_methy <- run_med(X = X, Y = Y, M = M, COV = Z)


##################################
## IDH WT + MGMT unmethy ##
## add codeletion age, gender, grade as covariate
##################################
sampleID <- dat.clin.glioma.toUse %>%
  filter(IDH == "WT" & MGMT == "Unmethylated") %>%
  pull(sample)

X <- t(dat.met.glioma.mval.scale.wt[,sampleID]) # exposures
M <- t(dat.exp.glioma.norm.scale.wt[,sampleID]) # mediators
Z <- dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% 
  select(c("paper_Grade", "age_scale", "gender")) # covariates
Y = Surv(dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% pull(overall_survival), 
         dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% pull(status))
res_med_wt_unmethy <- run_med(X = X, Y = Y, M = M, COV = Z)


##################################
## IDH MT ##
## add MGMT, ch1p19q, age, gender, grade as covariate
##################################
sampleID <- dat.clin.glioma.toUse %>%
  filter(IDH == "Mutant" & !is.na(MGMT) & !is.na(ch1p19q)) %>%
  pull(sample)

X <- t(dat.met.glioma.mval.scale.mt[,sampleID]) # exposures
M <- t(dat.exp.glioma.norm.scale.mt[,sampleID]) # mediators
Z <- dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% 
  select(c("paper_Grade", "age_scale", "gender","MGMT","ch1p19q")) # covariates
Y = Surv(dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% pull(overall_survival), 
         dat.clin.glioma.toUse %>% filter(sample %in% sampleID) %>% pull(status))
res_med_mt <- run_med(X = X, Y = Y, M = M, COV = Z)
