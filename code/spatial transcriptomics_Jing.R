# install.packages("devtools")
# devtools::install_github("analyxcompany/ForceAtlas2")
library(Seurat)
# library(ForceAtlas2)
# library(tidyverse)
# library(destiny)
# library(cccd)
library(stringr)
library(umap)
library(irlba) # for fast pca
library(vegan) # distance function
library(preprocessCore) # for quantile normalization
library(plyr)
memory.limit(size=54000)
output_dir <- "../data/output/Jing_code/"
if (!dir.exists(output_dir)){
  dir.create(output_dir)
}
# setwd("C:\\Users\\jtang\\Dropbox\\FIMM\\Olli\\Spatial transcriptomics")
rm(list=ls(all=TRUE))
# data <- readRDS("merged_samples.RDS") # 17677/17943 genes, 9129 cells
# data <- readRDS("merged_samples_old.RDS") # 17677/17943 genes, 9129 cells
data <- readRDS("../data/output/cellranger_agg/merged_samples.RDS") # 17677/17943 genes, 9129 cells. no big difference

# # Integration - does not work
# # https://satijalab.org/seurat/archive/v3.1/integration.html
# data.list = SplitObject(data, split.by = "Plate")
# for (i in names(data.list)) {
#   data.list[[i]] <- SCTransform(data.list[[i]], assay = "Spatial", verbose = FALSE)
# }
# 
# data.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
# data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = data.features)
# data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT",
#                                        anchor.features = data.features)
# data.integrated <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT")
# 
# data.integrated <- RunPCA(object = data.integrated, verbose = FALSE)
# data.integrated <- RunUMAP(object = data.integrated, dims = 1:30)
# 
# plots <- DimPlot(data.integrated, group.by = c("Sample"))
# plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 4, byrow = TRUE,
#                                                                      override.aes = list(size = 2.5)))

# # average count of a gene in a cell
# tmp = apply(data@assays$SCT@counts, 1, mean)
# summary(tmp)
# 
# # number of cells that are detecting 
# tmp1 = apply(data@assays$SCT@counts, 1, function(x) length(which(x > 1)))
# summary(tmp1)

# preprocessing
# names(data@meta.data)
# data = SCTransform(data, assay = "Spatial", verbose = F) # regress out patient or not? no need as the aim is to see the global structure

data2 = data@assays$SCT@scale.data # scale data or data or count?
# data2 = data@assays$SCT@data # scale data or data or count? data contains technical artifact, scale data may be more biologicall relevant

data2 = data.frame(data2)
dim(data2)

# quantile normalization
# data2 = normalize.quantiles(as.matrix(data2))

group1 = ifelse(str_detect(data$Sample, "marker"), "benign","cancer") # cancer or benign?
group2 = ifelse(str_detect(data$Sample, "BG"), "BG", "HC") # BG or HC?
group2[which(group1 == "benign")] = "benign"
group3 = ifelse(str_detect(data$Sample, "S_"), "S", "L") # S or L?
group3[which(group1 == "benign")] = "benign"
group4 = unlist(lapply(data$Sample, function(x) strsplit(x, "_")[[1]][2])) # patient id

# ------------------------------------------------------
# Q1 - whether cancer is different from benign samples?
# ------------------------------------------------------
P <- prcomp_irlba(data2, n = 5) # fast PCA, n PCs, bigger n more separate clusters, you can test with higher n
summary(P)
plot(P$rotation[,1], P$rotation[,2]) # PCA does not separate well

data2.umap = umap(P$rotation) # only work on the PCs
# head(data2.umap$layout, 3)
png(paste0(output_dir, "q1.png"))
plot(data2.umap$layout, col = factor(group1), pch = 16, xlim = c(-15,20), ylim = c(-15,15))
legend("topleft", legend = levels(factor(group1)), pch = 16, col = 1:length(factor(group1)))
dev.off()
# -----------------------------------------------
# Q2 - whether BG and HC are separated?
# -----------------------------------------------
png(paste0(output_dir, "q2.png"))
plot(data2.umap$layout, col = factor(group2), pch = 16, xlim = c(-15,20), ylim = c(-15,15))
legend("topleft", legend = levels(factor(group2)), pch = 16, col = 1:length(factor(group2)))
dev.off()
# ------------------------------------------------
# Q3 - whether S and L are separated?
# ------------------------------------------------
png(paste0(output_dir, "q3.png"))
plot(data2.umap$layout, col = factor(group3), pch = 16, xlim = c(-15,20), ylim = c(-15,15))
legend("topleft", legend = levels(factor(group3)), pch = 16, col = 1:length(factor(group3)))
dev.off()
# ------------------------------------------------
# Q4 - whether patients are separated?
# ------------------------------------------------
png(paste0(output_dir, "q4.png"))
plot(data2.umap$layout, col = factor(group4), pch = 16, xlim = c(-15,20), ylim = c(-15,15))
legend("topleft", legend = levels(factor(group4)), pch = 16, cex = 0.5, pt.cex = 1, col = 1:length(factor(group4)))
dev.off()

# ------------------------------------------------
# Q5 - whether HC predict survival better than BG?
# ------------------------------------------------
index1 = which(group2 == "HC")
index2 = which(group2 == "BG")

png(paste0(output_dir, "q5.png"), width = 10, height = 5, units = "in", res = 100)
par(mfrow = c(1, 2))
plot(data2.umap$layout[index1,], col = factor(group3[index1]), pch = 16, xlim = c(-15,20), ylim = c(-15,15), main = "HC cells")
legend("topleft", legend = levels(factor(group3[index1])), pch = 16, cex = 1, pt.cex = 1, col = 1:length(factor(group3[index1])))

plot(data2.umap$layout[index2,], col = factor(group3[index2]), pch = 16, xlim = c(-15,20), ylim = c(-15,15), main = "BG cells")
legend("topleft", legend = levels(factor(group3[index2])), pch = 16, cex = 1, pt.cex = 1, col = 1:length(factor(group3[index2])))
par(mfrow = c(1, 1))
dev.off()

# simple analysis - HC cells between S and L patients have higher distance than BG cells
data2.dist = vegdist(data2.umap$layout, method = "euclidean") # vegan library
data2.dist = as.matrix(data2.dist) # pairwise distance of all cells in the umap

# randomly take N cells in each patient group (L or S)
N = 1000
I = 1000
res = mat.or.vec(I, 3)
for(i in 1:I){
  cat("\r", i)
  index.L = sample(which(group3 == "L"), N)
  index.S = sample(which(group3 == "S"), N)
  # check how many are HC
  index.L.HC = index.L[which(group2[index.L] == "HC")] # HC cells for L 
  index.S.HC = index.S[which(group2[index.S] == "HC")] # HC cells for S
  
  index.L.BG = setdiff(index.L, index.L.HC) # BG cells for L
  index.S.BG = setdiff(index.S, index.S.HC) # BG cells for S
  
  # HC cells
  # within group distance
  dist1 = data2.dist[index.L.HC,index.L.HC]
  dist2 = data2.dist[index.S.HC,index.S.HC]
  
  # between group distance
  dist3 = data2.dist[index.L.HC,index.S.HC]
  dist4 = data2.dist[index.S.HC,index.L.HC]
  
  # t.test(c(dist1,dist2), c(dist3,dist4)) # should be significant
  
  # BG cells
  # within group distance
  dist5 = data2.dist[index.L.BG,index.L.BG]
  dist6 = data2.dist[index.S.BG,index.S.BG]
  
  # between group distance
  dist7 = data2.dist[index.L.BG,index.S.BG]
  dist8 = data2.dist[index.S.BG,index.L.BG]
  
  # t.test(c(dist5,dist6), c(dist7,dist8)) # most likely also significant
  
  # L/S cells within group distance
  dist9 = data2.dist[index.L,index.L]
  dist10 = data2.dist[index.S,index.S]
  
  # between group distance
  dist11 = data2.dist[index.L,index.S]
  dist12 = data2.dist[index.S,index.L]
  
  # t.test(c(dist9,dist10), c(dist11,dist12)) # most likely also significant
  
  
  # HC cells vs BG cells
  res[i, 1] = mean(c(dist3, dist4)) - mean(c(dist1, dist2))
  res[i, 2] = mean(c(dist7, dist8)) - mean(c(dist5, dist6))
  res[i, 3] = mean(c(dist11, dist12)) - mean(c(dist9, dist10))
}

t.test(res[,1], alternative = "greater") # between-group HC cells are more distanced than within-group HC cells
t.test(res[,2], alternative = "greater") # between-group BC cells are more distanced than within-group BC cells
t.test(res[,3], alternative = "greater") # between-group cells are more distanced than within-group cells
t.test(res[,1], res[,2], alternative = "greater") # HC cells separate L/S better than BG cells

delta = c(res[,1], res[,2])
group = c(rep("HC", 1000), rep("BG", 1000))
res.plot = data.frame(delta = delta, group = group)

png("q5_b.png", width = 5, height = 5, units = "in", res = 100)
ggplot(res.plot, aes(x = group, y = delta, color = group)) + geom_violin(trim = FALSE) +
  geom_jitter(position = position_jitter(0.2)) + labs(y = "L - d", x = "") + theme(legend.position = "none")
dev.off()

# -----------------------------------------------------
# Q6 averaging cells for each group and make a PCA plot
# -----------------------------------------------------
data3 = data.frame(t(data2))
data3$group1 = group1 # benign or cancer
data3$group2 = group2 # HC or BG
data3$group3 = group3 # S or L
data3$group4 = group4 # patient id
df2 <- ddply(data3, c("group4","group3","group2"), numcolwise(mean))

df2_HC = subset(df2, group2 == "HC")
df2_BG = subset(df2, group2 == "BG")
P1 <- prcomp_irlba(t(df2_HC[, -c(1:3)]), n = 3, scale. = T) 
summary(P1)
P2 <- prcomp_irlba(t(df2_BG[, -c(1:3)]), n = 3, scale. = T) 
summary(P2)

png("q5_c.png", width = 10, height = 5, units = "in", res = 100)
par(mfrow = c(1, 2))
plot(P1$rotation[,1], P1$rotation[,2], type = "n", main = "HC")
text(P1$rotation[,1], P1$rotation[,2], labels = df2_HC$group3, col = as.numeric(as.factor(df2_HC$group3)))

plot(P2$rotation[,1], P2$rotation[,2], type = "n", main = "BG")
text(P2$rotation[,1], P2$rotation[,2], labels = df2_BG$group3, col = as.numeric(as.factor(df2_BG$group3)))

par(mfrow = c(1, 1))
dev.off()


# ------------------------------------------------
# Q7 predicting the L/S state using HC/BG cells
# ------------------------------------------------
library(caret)
library(pROC)
library(h2o) # machine learning
library(Seurat)
library(stringr)
library(umap)
library(irlba) # for fast pca
library(vegan) # distance function
library(preprocessCore) # for quantile normalization
library(plyr)
memory.limit(size=54000)
setwd("C:\\Users\\jtang\\Dropbox\\FIMM\\Olli\\Spatial transcriptomics")
rm(list=ls(all=TRUE))
# data <- readRDS("merged_samples.RDS") # 17677/17943 genes, 9129 cells
# data <- readRDS("merged_samples_old.RDS") # 17677/17943 genes, 9129 cells
data <- readRDS("merged_samples_old_agg.RDS") # 17677/17943 genes, 9129 cells

# data = SCTransform(data, assay = "Spatial", vars.to.regress = "patient", verbose = F) # regress out patient
# data = SCTransform(data, assay = "Spatial", verbose = F)

data2 = data@assays$SCT@data
# data2 = data@assays$SCT@scale.data

data2 = data.frame(data2)
dim(data2)

# quantile normalization
# data2 = normalize.quantiles(as.matrix(data2))

group1 = ifelse(str_detect(data$Sample, "marker"), "benign","cancer") # cancer or benign?
group2 = ifelse(str_detect(data$Sample, "BG"), "BG", "HC") # BG or HC?
group2[which(group1 == "benign")] = "benign"
group3 = ifelse(str_detect(data$Sample, "S_"), "S", "L") # S or L?
group3[which(group1 == "benign")] = "benign"
group4 = unlist(lapply(data$Sample, function(x) strsplit(x, "_")[[1]][2])) # patient id

index1 = which(group2 == "HC")
index2 = which(group2 == "BG")

data_HC = data.frame(t(data2[,index1])) # 4166 cells, X features
data_BG = data.frame(t(data2[,index2])) # 4090 cells, X features

data_HC = cbind(data_HC, group = as.factor(group3[index1]), patient = as.factor(group4[index1]))
data_BG = cbind(data_BG, group = as.factor(group3[index2]), patient = as.factor(group4[index2]))

# ----------------------------------------------
# leave-one-patient-out
# ----------------------------------------------
patient = unique(group4)
lapply(patient, function(x) unique(group2[which(group4==x)])) # OVA38 does not have BG
patient = patient[patient != "OVA38"]
res = mat.or.vec(length(patient), 2)

h2o.init()
for(i in 1:length(patient)){
  # cat("\r", i)
  print(i)
  trainIndex = which(data_HC$patient != patient[i]) # the other patient as training data
  data_HC_train = data_HC[trainIndex, -dim(data_HC)[2]] # remove the patient id
  data_HC_test = data_HC[-trainIndex, -dim(data_HC)[2]]
  

  data_HC_train = as.h2o(data_HC_train)
  data_HC_test = as.h2o(data_HC_test)
  # fit_HC = h2o.glm(y = "group", training_frame = data_HC_train, family = "binomial", lambda_search = T) # alpha = 0.5 by default
  fit_HC = h2o.gbm(y = "group", training_frame = data_HC_train)
  
  pred_HC = h2o.predict(fit_HC, data_HC_test)
  
  pred = factor(as.vector(pred_HC[,"predict"]), levels = c("L","S"))
  true = factor(as.vector(data_HC_test$group), levels = c("L","S"))
  
  res[i, 1] = confusionMatrix(pred, true)$overall['Accuracy']
  h2o.removeAll()
  # h2o.shutdown(prompt = F)
  
  trainIndex = which(data_BG$patient != patient[i]) # the other patient as training data  
  data_BG_train = data_BG[trainIndex, -dim(data_BG)[2]]
  data_BG_test = data_BG[-trainIndex, -dim(data_BG)[2]]
  
  # h2o.init()
  data_BG_train = as.h2o(data_BG_train)
  data_BG_test = as.h2o(data_BG_test)
  # fit_BG = h2o.glm(y = "group", training_frame = data_BG_train, family = "binomial", lambda_search = T) # alpha = 0.5 by default
  fit_BG = h2o.gbm(y = "group", training_frame = data_BG_train)
  pred_BG = h2o.predict(fit_BG, data_BG_test)
  
  pred = factor(as.vector(pred_BG[,"predict"]), levels = c("L","S"))
  true = factor(as.vector(data_BG_test$group), levels = c("L","S"))
  
  res[i, 2] = confusionMatrix(pred, true)$overall['Accuracy']
  h2o.removeAll()
  print(res[i,])
  
}
h2o.shutdown(prompt = F)
colnames(res) = c("HC","BG")
rownames(res) = patient
res = data.frame(res)
res$PFI = unlist(lapply(patient, function(x) unique(group3[which(group4==x)])[1]))

library(openxlsx)
write.xlsx(res, "res_gbm.xlsx", rowNames = T)

# -------------------------------------------------------
# Q8: DEG
# -------------------------------------------------------
library(caret)
library(pROC)
library(h2o) # machine learning
library(Seurat)
library(stringr)
library(umap)
library(irlba) # for fast pca
library(vegan) # distance function
library(preprocessCore) # for quantile normalization
library(plyr)
library(sparseMatrixStats) # sparse matrix calculation
library(openxlsx)

library(glue)
library(ggsci)
library(survival)
library(ggfortify)
library(rlang)
library(survminer)
library(mice)
library(missForest)
library(rcompanion)
library(reshape2)
library(corrplot)
library(RColorBrewer)
library(rje)
library(Hmisc)
library(survminer)
library(readr)
library(tidyverse)
library(openxlsx)
library(survcomp)
library(reportROC)
library(gridExtra)
library(Hmisc)
library(vcdExtra)
require(ggplot2)
library(sparseMatrixStats) # sparse matrix calculation
memory.limit(size=54000)
setwd("C:\\Users\\jtang\\Dropbox\\FIMM\\Olli\\Spatial transcriptomics")
rm(list=ls(all=TRUE))
data <- readRDS("../data/output/cellranger_agg/merged_samples.RDS")
id_mappings = read.csv("../data/raw/gene_ensembl_symbol_table.csv")

data2 = data@assays$SCT@data
data2 = data.frame(data2)
dim(data2) # 17866 genes, 9129 cells

# group1 = ifelse(str_detect(data$Sample, "marker"), "benign","cancer") # cancer or benign?
# group2 = ifelse(str_detect(data$Sample, "BG"), "BG", "HC") # BG or HC?
# group2[which(group1 == "benign")] = "benign"
# group3 = ifelse(str_detect(data$Sample, "S_"), "S", "L") # S or L?
# group3[which(group1 == "benign")] = "benign"
# group4 = unlist(lapply(data$Sample, function(x) strsplit(x, "_")[[1]][2])) # patient id

# `%notin%` <- Negate(`%in%`)
# data3 = subset(data, subset = patient %notin% c("OVA38","OVA28","OVA10")) # OVA38 does not have BG data, OVA28 and OVA10 have poor prediction accuracy
# dim(data3)
data3 = data

# Load known BRCAness signatures
cov362 = read.csv("../../scRNA/lineage_tracing/2021.11.26_KURAMOCHI_AML_lineage_tracing/data/output/BRCAness/deg_COV362.csv") # JUN data
cov362[which(cov362$symbol == "BRCA1"),]

peng = read.xlsx("../data/raw/BRCAness_marker_Peng_PMID24553445_ensembl.xlsx") # Peng data
peng_supplementary = read.xlsx("../data/raw/Peng data.xlsx")
peng = merge(peng, peng_supplementary)
peng$ensembl = unlist(lapply(peng$ensembl, function(x) strsplit(x, "\\.")[[1]][1]))

# take all the samples
# res = FindMarkers(data3, ident.1 = colnames(data3)[which(data3$PFI == "Short")], 
#                   ident.2 = colnames(data3)[which(data3$PFI == "Long")], latent.vars = "patient", logfc.threshold = 0.25)
# 
# # sanity check
# res = data.frame(avg_log2FC = rowMeans2(data3@assays$SCT@data[, which(data3$PFI == "Short" & data3$Confidence == "High confidence")]) - 
#                      rowMeans2(data3@assays$SCT@data[, which(data3$PFI == "Long"& data3$Confidence == "High confidence")]))
# rownames(res) = rownames(data3)
# res["BRCA2",]
# 
# cor.test(res$avg_log2FC, cov362$avg_log2FC[match(rownames(res), cov362$symbol)])
# cor.test(res$avg_log2FC, peng$Gene.expression.of.BRCA1.shRNA.Control.shRNA.[match(rownames(res), peng$Gene.Symbol.)], method = "spearman")
# plot(2^res$avg_log2FC, peng$Gene.expression.of.BRCA1.shRNA.Control.shRNA.[match(rownames(res), peng$Gene.Symbol.)], xlim = c(0,2))
# 

# whether the HRD signature can predict patient prognosis?

# take only high confidence samples
res2 = FindMarkers(
  data3, ident.1 = colnames(data3)[which(data3$PFI == "Short" & data3$Confidence == "High confidence")], 
  ident.2 = colnames(data3)[which(data3$PFI == "Long" & data3$Confidence == "High confidence")], 
  latent.vars = "patient", logfc.threshold = 0.25, test.use = "negbinom")

res2$ensembl = id_mappings$ensembl[match(rownames(res2), id_mappings$symbol)]

# sanity check
res2_b = data.frame(
  avg_log2FC = rowMeans2(data3@assays$SCT@data[, which(data3$PFI == "Short" & data3$Confidence == "High confidence")]) - 
    rowMeans2(data3@assays$SCT@data[, which(data3$PFI == "Long"& data3$Confidence == "High confidence")])
  )
rownames(res2_b) = rownames(data3)
cor(res2$avg_log2FC, res2_b$avg_log2FC[match(rownames(res2), rownames(res2_b))]) # 0.95

cor.test(res2$avg_log2FC, cov362$avg_log2FC[match(rownames(res2), cov362$symbol)])
cor.test(res2$avg_log2FC, peng$Gene.expression.of.BRCA1.shRNA.Control.shRNA.[match(res2$ensembl, peng$ensembl)], method = "spearman") # 0.54

plot(res2$avg_log2FC, peng$Gene.expression.of.BRCA1.shRNA.Control.shRNA.[match(res2$ensembl, peng$ensembl)])

write.xlsx(res2, "../data/output/Jing_code/res2.xlsx", rowNames = T)

x = 2^res2$avg_log2FC
y =  peng$Gene.expression.of.BRCA1.shRNA.Control.shRNA.[match(res2$ensembl, peng$ensembl)]

png("../data/output/Jing_code/q8_a.png", width = 5, height = 5, units = "in", res = 100)
plot(x, y, 
     xlab = "Short/Long PFI", ylab = "HRD signature genes", xlim = c(0,3), type = "n", las = 1)
points(x[which(y<=1 & x <=1)], y[which(y <=1 & x <=1)], pch = 16)
text(x[which(y>0 & x>0)], y[which(y>0 & x >0)], labels = rownames(res2)[which(y>0 & x >0)])
abline(h = 1, lty = 2)
abline(v = 1, lty = 2)
dev.off()

# low confidence samples
res3 = FindMarkers(
  data3, 
  ident.1 = colnames(data3)[which(data3$PFI == "Short" & data3$Confidence == "Low confidence")],
  ident.2 = colnames(data3)[which(data3$PFI == "Long" & data3$Confidence == "Low confidence")],
  latent.vars = "patient",
  logfc.threshold = 0.25,
  test.use = "negbinom"
  )

res3$ensembl = id_mappings$ensembl[match(rownames(res3), id_mappings$symbol)]
cor.test(
  res3$avg_log2FC,
  cov362$avg_log2FC[match(rownames(res3), cov362$symbol)]
  )

cor.test(
  res3$avg_log2FC,
  peng$Gene.expression.of.BRCA1.shRNA.Control.shRNA.[match(res3$ensembl, peng$ensembl)],
  method = "spearman"
  ) # 0.46

res_all = merge(res2, res3, by = "row.names")

res_all$diff_log2FC <- res_all$avg_log2FC.x - res_all$avg_log2FC.y
write.xlsx(res_all, "../data/output/Jing_code/res_all.xlsx")

# load TCGA
setwd("C:\\Users\\jtang\\Dropbox\\FIMM\\Anna\\Data")
data_clinical_patient = read.delim("ov_tcga_pan_can_atlas_2018.tar//ov_tcga_pan_can_atlas_2018//data_clinical_patient.txt", skip = 4)
data_clinical_sample = read.delim("ov_tcga_pan_can_atlas_2018.tar//ov_tcga_pan_can_atlas_2018//data_clinical_sample.txt", skip = 4)
data_mrna = read.delim("ov_tcga_pan_can_atlas_2018.tar//ov_tcga_pan_can_atlas_2018//data_mrna_seq_v2_rsem.txt", check.names = F)
# data_mrna = read.delim("ov_tcga_pan_can_atlas_2018.tar//ov_tcga_pan_can_atlas_2018//data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt", check.names = F)
setwd("C:\\Users\\jtang\\Dropbox\\FIMM\\Olli\\Spatial transcriptomics")

identical(data_clinical_patient$PATIENT_ID, data_clinical_sample$PATIENT_ID) # data aligned already
table(data_clinical_sample$GRADE)
high_grade = grepl("G2|G3|G4", data_clinical_sample$GRADE)

data_clinical_sample = data_clinical_sample[high_grade, ]
data_clinical_patient = data_clinical_patient[high_grade, ]

patients = names(data_mrna)[-c(1:2)] # 300 patients
time = data_clinical_patient$PFS_MONTHS[match(patients, data_clinical_sample$SAMPLE_ID)]
status = data_clinical_patient$PFS_STATUS[match(patients, data_clinical_sample$SAMPLE_ID)]
table(status, exclude = NULL)
status[which(status == "")] = NA

length(which(is.na(time)==F)) # 271 has survival
length(which(is.na(status)==F)) # 272 has status

# gene = rownames(res2[order(res2$avg_log2FC, decreasing = T),])
# gene = rownames(res2_b)[res2_b$avg_log2FC > 0.2 | res2_b$avg_log2FC < -0.2]
gene2 = rownames(res2[order(res2$p_val_adj, decreasing = F),])[1:700] # from HC samples
gene3 = rownames(res3[order(res3$p_val_adj, decreasing = F),])[1:700] # from BG samples
library(ggvenn) # venn diagram
x = list(HC = gene2, LC = gene3)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 8, text_size = 6, show_percentage = F 
)


gene = gene2
res_gene = matrix(NA, nrow = length(gene), ncol = 2)
for(i in 1:length(gene)){
  cat('\r', i)
  index = which(data_mrna$Hugo_Symbol == gene[i])
  if(length(index) == 1){
    res1 = as.numeric(data_mrna[index, -c(1:2)])
    if(!all(is.na(res1))){
      group = ifelse(res1 > quantile(res1, 0.5), 1, 0)
      
      data = data.frame(time = time, status = as.numeric(as.factor(status)), group = group)
      fit <- survfit(Surv(time, status) ~ group, data = data)
      fit.coxph = coxph(Surv(time, status) ~ group, data = data)
      res_gene[i, 1] = summary(fit.coxph)$logtest[3]
    }
    # concordance(fit.coxph)$concordance
  }
  
  index_random = sample(1:dim(data_mrna)[1], 1)
  if(length(index_random) == 1){
    res1 = as.numeric(data_mrna[index_random, -c(1:2)])
    if(!all(is.na(res1))){
      group = ifelse(res1 > quantile(res1, 0.5), 1, 0)
      
      data = data.frame(time = time, status = as.numeric(as.factor(status)), group = group)
      fit <- survfit(Surv(time, status) ~ group, data = data)
      fit.coxph = coxph(Surv(time, status) ~ group, data = data)
      res_gene[i, 2] = summary(fit.coxph)$logtest[3]
    }
    # concordance(fit.coxph)$concordance
  }
}
rownames(res_gene) = gene
summary(res_gene)
ks.test(res_gene[,1], res_gene[,2])
t.test(res_gene[,1], res_gene[,2], paired = T)
which(res_gene[,1] < 0.003)

i = 25
res2[i,]
index = which(data_mrna$Hugo_Symbol == gene[i])
# index = which(data_mrna$Hugo_Symbol == "FAM43A")

res1 = as.numeric(data_mrna[index, -c(1:2)])
group = ifelse(res1 > quantile(res1, 0.5), 1, 0)
data = data.frame(time = time, status = as.numeric(as.factor(status)), group = group)
fit <- survfit(Surv(time, status) ~ group, data = data)
fit.coxph = coxph(Surv(time, status) ~ group, data = data)
summary(fit.coxph)

setwd("C:\\Users\\jtang\\Dropbox\\FIMM\\Olli\\Spatial transcriptomics")
png("q8_d_3.png", width = 7, height = 12, units = "in", res = 300, pointsize = 30)
ggsurvplot(fit, data = data, palette = c("green","red"), risk.table = T, legend.title = gene[i], legend.labs = c("low","high"), 
           xlab = "Time (months)", ylab = "PFS",
           ggtheme = theme_classic(base_size = 30), fontsize = 8) + guides(colour = guide_legend(nrow = 2))
dev.off()

