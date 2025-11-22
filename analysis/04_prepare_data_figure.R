
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
