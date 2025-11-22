# pre-process gene expression data
## top expressions with less 0
## Cindex
## survival cutoff: 1yr, 2yr, 3yr and 5yr survival


## input data
dat.clin.glioma <- readRDS(file = paste0(RESULT,"dat_exp_clin_20250205.rds"))
dat.exp.glioma.toUse <- readRDS(file = paste0(DATA, "glioma_exp_clean.rds"))
dat.exp.glioma <- as.data.frame(assay(dat.exp.glioma.toUse))
dat.exp.glioma <- dat.exp.glioma[,rownames(dat.clin.glioma)]

## step 1: top expression values
keep <- rowSums(dat.exp.glioma >= 10) >= 2
dat.exp.glioma <- dat.exp.glioma[keep,]

# normalization
dds <- DESeqDataSetFromMatrix(countData = dat.exp.glioma, 
                              colData = dat.clin.glioma, 
                              design = ~ cancer)
dds <- estimateSizeFactors(dds)
dat.exp.glioma.norm <- counts(dds, normalized=TRUE)
colnames(dat.exp.glioma.norm) <- dat.clin.glioma$sample

# duplicate gene name for "ENSG00000158427.15" "ENSG00000269226.7" 
# search manually: ENSG00000158427.15->TMSB15B; ENSG00000269226.7->TMSB15C
which(rownames(dat.exp.glioma.norm) == "ENSG00000269226.7")
rownames(dat.exp.glioma.norm)[-29915] <- rowRanges(dat.exp.glioma.toUse)@elementMetadata@listData$gene_name[which(rowRanges(dat.exp.glioma.toUse)@elementMetadata@listData$gene_id %in% rownames(dat.exp.glioma)[-29915])]
rownames(dat.exp.glioma.norm)[29915] <- "TMSB15C"

## step 2: DE analysis
# DE analysis between survival group (same as DNA methylation data)
DEG_analysis_survival <- function(data, count){
  # IDH WT
  dat.clin.wt <- data %>% 
    filter(paper_IDH.status == "WT") %>%
    mutate(survival = as.factor(survival)) %>%
    mutate(gender = as.factor(gender))
  
  # fit DESeq2
  dds.wt <- DESeqDataSetFromMatrix(countData = count[,rownames(dat.clin.wt)],
                                   colData = dat.clin.wt,
                                   design= ~ survival + age_at_diagnosis + gender + paper_Grade + paper_MGMT.promoter.status)
  dds.wt <- DESeq(dds.wt)
  res.wt <- results(dds.wt, name="survival_TRUE_vs_FALSE", alpha = 0.1)
  
  # extracting DE genes with adj.p<0.1
  DEG.wt <- as.data.frame(subset(res.wt, padj < 0.1)) 
  
  # IDH MT
  dat.clin.mt <- data %>% 
    filter(paper_IDH.status == "Mutant")  %>%
    mutate(survival = as.factor(survival)) %>%
    mutate(gender = as.factor(gender))
  
  dds.mt <- DESeqDataSetFromMatrix(countData = count[,rownames(dat.clin.mt)],
                                   colData = dat.clin.mt,
                                   design= ~ survival + age_at_diagnosis + gender + paper_Grade + paper_MGMT.promoter.status)
  dds.mt <- DESeq(dds.mt)
  res.mt <- results(dds.mt, name="survival_TRUE_vs_FALSE", alpha = 0.1)
  DEG.mt <- as.data.frame(subset(res.mt, padj < 0.1)) 
  
  return(list(DEG.wt = DEG.wt, DEG.mt = DEG.mt))
}

DEG_res_1yr <- DEG_analysis_survival(dat.clin.glioma.1yr, dat.exp.glioma)
DEG_res_2yr <- DEG_analysis_survival(dat.clin.glioma.2yr, dat.exp.glioma)
DEG_res_3yr <- DEG_analysis_survival(dat.clin.glioma.3yr, dat.exp.glioma)
DEG_res_5yr <- DEG_analysis_survival(dat.clin.glioma.5yr, dat.exp.glioma)

exp.filter2.wt <- union(rownames(DEG_res_1yr$DEG.wt), 
                        union(rownames(DEG_res_2yr$DEG.wt), 
                              union(rownames(DEG_res_3yr$DEG.wt), rownames(DEG_res_5yr$DEG.wt))))
exp.filter2.mt <- union(rownames(DEG_res_1yr$DEG.mt), 
                        union(rownames(DEG_res_2yr$DEG.mt), 
                              union(rownames(DEG_res_3yr$DEG.mt), rownames(DEG_res_5yr$DEG.mt))))


## step 3: C-index
# calculate for each feature i, using beta value
# 2 direction: up-regulated (short survival time have higher beta values) and down-reg
dat.exp.glioma.wt <- dat.exp.glioma[,rownames(dat.clin.glioma)[which(dat.clin.glioma$paper_IDH.status == "WT")]]
dat.clin.wt <- dat.clin.glioma %>% 
  filter(paper_IDH.status == "WT")
cindex_exp_hyper.wt <- NULL
cindex_exp_hypo.wt <- NULL
for (i in 1:nrow(dat.exp.glioma.wt)) {
  time <- dat.clin.wt$overall_survival
  status <- dat.clin.wt$status
  S <- Surv(time, status)
  x <- t(dat.exp.glioma.wt)[,i]
  
  cindex_exp_hyper.wt[i] <- rcorr.cens(x, S,  outx=FALSE)[1]
  cindex_exp_hypo.wt[i] <- 1 - cindex_exp_hyper.wt[i]
}
cindex_exp.wt <- pmax(cindex_exp_hyper.wt, cindex_exp_hypo.wt)
hist(cindex_exp.wt)
quantile(cindex_exp.wt)
quantile(cindex_exp.wt, probs=0.9)
exp.filter3.wt <- rownames(dat.exp.glioma.wt)[which(cindex_exp.wt >= quantile(cindex_exp.wt, probs=0.9))]

dat.exp.glioma.mt <- dat.exp.glioma[,rownames(dat.clin.glioma)[which(dat.clin.glioma$paper_IDH.status == "Mutant")]]
dat.clin.mt <- dat.clin.glioma %>% 
  filter(paper_IDH.status == "Mutant")
cindex_exp_hyper.mt <- NULL
cindex_exp_hypo.mt <- NULL
for (i in 1:nrow(dat.exp.glioma.mt)) {
  time <- dat.clin.mt$overall_survival
  status <- dat.clin.mt$status
  S <- Surv(time, status)
  x <- t(dat.exp.glioma.mt)[,i]
  
  cindex_exp_hyper.mt[i] <- rcorr.cens(x, S,  outx=FALSE)[1]
  cindex_exp_hypo.mt[i] <- 1 - cindex_exp_hyper.mt[i]
}
cindex_exp.mt <- pmax(cindex_exp_hyper.mt, cindex_exp_hypo.mt)
hist(cindex_exp.mt)
quantile(cindex_exp.mt)
exp.filter3.mt <- rownames(dat.exp.glioma.mt)[which(cindex_exp.mt >= quantile(cindex_exp.mt, probs=0.9))]

