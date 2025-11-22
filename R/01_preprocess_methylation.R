# pre-process DNA methylation data
## top variance
## Cindex
## survival cutoff: 1yr, 2yr, 3yr and 5yr survival
## filtering for IDH WT vs MT separately


## input data
dat.clin.glioma <- readRDS(file = paste0(RESULT,"dat_met_clin_20250205.rds"))
dat.met.glioma.toUse <- readRDS(file = paste0(DATA, "glioma_met_clean.rds"))
dat.met.glioma <- as.data.frame(assay(dat.met.glioma.toUse))
dat.met.glioma <- dat.met.glioma[,rownames(dat.clin.glioma)]

## step 1: most variable CpGs (top 10%)
topVar <- featurefilter(dat.met.glioma, percentile = 10, method = "var", topN = 5)
met.filter1 <- rownames(topVar$filtered_data)
which(met.filter1=="cg10586317") ## CD276
dat.met.glioma <- dat.met.glioma[met.filter1,]
colnames(dat.met.glioma) <- dat.clin.glioma$sample

# step 2: differential methylation analysis
# fit the linear model (limma)
# converting beta values to m_values
mval <- t(apply(dat.met.glioma, 1, function(x) log2(x/(1-x))))
all(colnames(mval) %in% row.names(dat.clin.glioma))

# DMC analysis between survival group
DMC_analysis_survival <- function(data, mval){
  # IDH WT
  dat.clin.wt <- data %>% 
    filter(paper_IDH.status == "WT")
  
  # fit limma
  design <- model.matrix(~ survival + age_at_diagnosis + gender + paper_Grade + paper_MGMT.promoter.status, data = dat.clin.wt)
  fit.wt <- lmFit(mval[,rownames(dat.clin.wt)], design)
  fit.wt <- eBayes(fit.wt)
  
  # extracting significantly methylated probes with adj.p<0.1
  DMC.wt <- topTable(fit.wt, coef = 2, sort.by = "p",number = nrow(mval), adjust.method = "BH", p.value = 0.1)
  
  # IDH MT
  dat.clin.mt <- data %>% 
    filter(paper_IDH.status == "Mutant") 
  
  design <- model.matrix(~ survival + age_at_diagnosis + gender + paper_Grade + paper_MGMT.promoter.status, data = dat.clin.mt)
  fit.mt <- lmFit(mval[,rownames(dat.clin.mt)], design)
  fit.mt <- eBayes(fit.mt)
  
  DMC.mt <- topTable(fit.mt, coef = 2, sort.by = "p",number = nrow(mval), adjust.method = "BH", p.value = 0.1)
  
  rm(dat.clin.wt);rm(dat.clin.mt)
  
  return(list(DMC.wt = DMC.wt, DMC.mt = DMC.mt))
}
DMC_res_1yr <- DMC_analysis_survival(dat.clin.glioma.1yr, mval)
DMC_res_2yr <- DMC_analysis_survival(dat.clin.glioma.2yr, mval)
DMC_res_3yr <- DMC_analysis_survival(dat.clin.glioma.3yr, mval)
DMC_res_5yr <- DMC_analysis_survival(dat.clin.glioma.5yr, mval)

met.filter2.wt <- union(rownames(DMC_res_1yr$DMC.wt),
                        union(rownames(DMC_res_2yr$DMC.wt),
                              union(rownames(DMC_res_3yr$DMC.wt), rownames(DMC_res_5yr$DMC.wt))))
met.filter2.mt <- union(rownames(DMC_res_1yr$DMC.mt),
                        union(rownames(DMC_res_2yr$DMC.mt),
                              union(rownames(DMC_res_3yr$DMC.mt), rownames(DMC_res_5yr$DMC.mt))))

# step 3: C-index
# calculate for each feature i, using beta value
# 2 direction: hyper-methy (short survival time have higher beta values) and hypo-methy
# calculate cindex
dat.met.glioma.wt <- dat.met.glioma[met.filter1,rownames(dat.clin.glioma)[which(dat.clin.glioma$paper_IDH.status == "WT")]]
dat.clin.wt <- dat.clin.glioma %>% 
  filter(paper_IDH.status == "WT")
cindex_met_hyper.wt <- NULL
cindex_met_hypo.wt <- NULL
for (i in 1:nrow(dat.met.glioma.wt)) {
  time <- dat.clin.wt$overall_survival
  status <- dat.clin.wt$status
  S <- Surv(time, status)
  x <- t(dat.met.glioma.wt)[,i]
  
  cindex_met_hyper.wt[i] <- rcorr.cens(x, S,  outx=FALSE)[1]
  cindex_met_hypo.wt[i] <- 1 - cindex_met_hyper.wt[i]
}
cindex_met.wt <- pmax(cindex_met_hyper.wt, cindex_met_hypo.wt)
met.filter3.wt <- rownames(dat.met.glioma.wt)[which(cindex_met.wt >= quantile(cindex_met.wt, probs=0.75))]

dat.met.glioma.mt <- dat.met.glioma[met.filter1,rownames(dat.clin.glioma)[which(dat.clin.glioma$paper_IDH.status == "Mutant")]]
dat.clin.mt <- dat.clin.glioma %>% 
  filter(paper_IDH.status == "Mutant")
cindex_met_hyper.mt <- NULL
cindex_met_hypo.mt <- NULL
for (i in 1:nrow(dat.met.glioma.mt)) {
  time <- dat.clin.mt$overall_survival
  status <- dat.clin.mt$status
  S <- Surv(time, status)
  x <- t(dat.met.glioma.mt)[,i]
  
  cindex_met_hyper.mt[i] <- rcorr.cens(x, S,  outx=FALSE)[1]
  cindex_met_hypo.mt[i] <- 1 - cindex_met_hyper.mt[i]
}
cindex_met.mt <- pmax(cindex_met_hyper.mt, cindex_met_hypo.mt)
hist(cindex_met.mt)
quantile(cindex_met.mt)
quantile(cindex_met.mt, probs=0.8)
met.filter3.mt <- rownames(dat.met.glioma.mt)[which(cindex_met.mt >= quantile(cindex_met.mt, probs=0.8))]






