# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("ggfortify")
# BiocManager::install("limma")
 # BiocManager::install("edgeR")
# BiocManager::install(c("DEGreport", "scatterpie", "ggrepel", "pheatmap"))

library(edgeR)
library(DESeq2)
library(ggfortify)
library(limma)
library(DEGreport)
library(scatterpie)
library(ggrepel)
library(pheatmap)
library(nnet)
library(glmnet)
library(caret)
library(ggplot2)

######## 1 mount ########
counts = readRDS("data/KO_sample_counts.rds")
info = readRDS("data/KO_sample_info.rds")
dds = DESeqDataSetFromMatrix(countData = counts, 
                             colData = info, 
                             design = ~ Sex + Group)

autoplot(prcomp(t(normalized_counts)), data = info, 
         colour = 'Group', shape = 'Sex', 
         size = 3, label.repel = 1)

# plot PCA
normalized_counts = assay(vst(dds))
autoplot(prcomp(t(normalized_counts)), data = info, 
         colour = 'Group', shape = 'Sex', 
         size = 3, label.repel = 1)

adjusted_counts=removeBatchEffect(normalized_counts, dds$Sex)
autoplot(prcomp(t(adjusted_counts)), data = info, 
         colour = 'Group', shape = 'Sex', 
         size = 3, label.repel = 1)


dds_lrt = DESeq(dds, test="LRT", reduced = ~ Sex)

adjusted_counts_lrt=removeBatchEffect(assay(vst(dds_lrt)), dds$Sex)
res_LRT = results(dds_lrt)
sig_genes = rownames(subset(res_LRT, res_LRT$padj<0.05))
clusters <- degPatterns(adjusted_counts_lrt[sig_genes,], 
                        metadata = info,
                        time = "Group", col=NULL,
                        minc = )

dim(adjusted_counts_lrt[sig_genes,])


# ANOVA ####
######## 2 LC ########
counts = readRDS("data/LC_sample_counts.rds")
info = readRDS("data/LC_sample_info.rds")

cpm <- cpm(counts)

keep.exprs <- filterByExpr(counts, )
counts_keep.exprs <- counts[keep.exprs, ]
dim(counts)
dim(counts_keep.exprs)

dds = DESeqDataSetFromMatrix(countData = counts_keep.exprs, 
                             colData = info, 
                             design = ~ Group)

normalized_counts = assay(vst(dds))

pcaData = prcomp(t(normalized_counts))
vars <- (pcaData$sdev)^2  
props <- vars / sum(vars)    
props

percentVar <- round(props*100)
pcaData = as.data.frame(cbind(rownames(info), pcaData$x[,1], pcaData$x[,2], info))
colnames(pcaData)[1:3] = c("name", "PC1", "PC2")
ggplot(aes(x=PC1, y=PC2), data = pcaData) + xlim(-40, 40) +
  xlab(paste0("PC1: ", percentVar[1],"%")) +
  ylab(paste0("PC2: ", percentVar[2],"%")) + 
  scatterpie::geom_scatterpie(aes(x=PC1, y=PC2, r=3), 
                              data=pcaData, alpha=0.4, 
                              cols=c("Feature_A", "Feature_B", "Feature_C")) + 
  geom_text_repel(aes_string(label = "name", color = "Group"), size = 6,
                  data          = subset(pcaData, Group == "A"),
                  nudge_x       = -25 - subset(pcaData, Group == "A")$PC1,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 0) +
  geom_text_repel(aes_string(label = "name", color = "Group"), size = 6,
                  data          = subset(pcaData, Group == "B"),
                  nudge_x       = 30 - subset(pcaData, Group == "B")$PC1,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 0) +
  geom_text_repel(aes_string(label = "name", color = "Group"), size = 6,
                  data          = subset(pcaData, Group == "C"),
                  nudge_x       = -10 - subset(pcaData, Group == "C")$PC1,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 0) 

Lc <- cbind(info, t(normalized_counts))
ppt <- Lc[, 1:10]

######## 2.1 Only A vs BC ########
# 4770 9994 9998
# 389  957 1168 2129 3429 4247 4599 4647 5049 5161 5526 6297 7155 7335 7356 7772 9231
find_ANOVA <- normalized_counts[c(389, 957, 1168, 2129, 3429, 4247, 
                                  4599, 4647, 5049, 5161, 5526, 6297, 
                                  7155, 7335, 7356, 7772, 9231), ]
counts_cor <- cor(find_ANOVA)
pheatmap(counts_cor, annotation = info)



counts_cor <- cor(normalized_counts)
pheatmap(counts_cor, annotation = info)

######## 2.2 Only A vs B vs C ########
# BC
temp_group <- c(rep("B", 5), rep("C", 5))
group_a_bc <- rbind(temp_group, normalized_counts[, 6:15])
head(group_a_bc)

group_a_bc <- as.data.frame(t(group_a_bc))
group_a_bc[1:3, 1:3]
i <- 2

A_BC_p.value <- data.frame(Gene= colnames(group_a_bc)[-1])

i <- 2
for (i in 2:dim(group_a_bc)[2]){
  fit <- aov(as.numeric(as.character(group_a_bc[,i]))~ 
               group_a_bc$temp_group, data= group_a_bc)
  A_BC_p.value$p[i-1] <- summary(fit)[[1]][["Pr(>F)"]][[1]]
} 

# significant <- which(A_BC_p.value$p < (0.05/10000)) 

BC_extract <- A_BC_p.value[order(A_BC_p.value$p), ][1:20, ]
plot(BC_extract$p)

# 畫cor #
temp <- strsplit(as.character(BC_extract$Gene), "_")
BCsign <- as.numeric(lapply(temp, `[`, 2))
# 直接吃2.1的
A_BCsign <- c(389, 957, 1168, 2129, 3429, 4247, 
              4599, 4647, 5049, 5161, 5526, 6297, 
              7155, 7335, 7356, 7772, 9231)

find_ANOVA <- normalized_counts[c(A_BCsign, BCsign), ]
counts_cor <- cor(find_ANOVA)
pheatmap(counts_cor, annotation = info)

######## 2.3 PCA ########
A_BCsign <- c(389, 957, 1168, 2129, 3429, 4247, 
              4599, 4647, 5049, 5161, 5526, 6297, 
              7155, 7335, 7356, 7772, 9231)

pca <- prcomp(t(find_ANOVA))
vars <- (pca$sdev)^2  
props <- vars / sum(vars)    
props

find_ANOVA <- normalized_counts[c(A_BCsign, BCsign), ]
pcaData = prcomp(t(find_ANOVA))
percentVar <- round(props*100)
pcaData = as.data.frame(cbind(rownames(info), pcaData$x[,1], pcaData$x[,2], info))
colnames(pcaData)[1:3] = c("name", "PC1", "PC2")
ggplot(aes(x=PC1, y=PC2), data = pcaData) +
  xlab(paste0("PC1: ", percentVar[1],"%")) +
  ylab(paste0("PC2: ", percentVar[2],"%")) + 
  scatterpie::geom_scatterpie(aes(x=PC1, y=PC2, r=3), 
                              data=pcaData, alpha=0.4, 
                              cols=c("Feature_A", "Feature_B", "Feature_C")) + 
  geom_text_repel(aes_string(label = "name", color = "Group"), size = 6,
                  data          = subset(pcaData, Group == "A"),
                  nudge_x       = -25 - subset(pcaData, Group == "A")$PC1,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 0) +
  geom_text_repel(aes_string(label = "name", color = "Group"), size = 6,
                  data          = subset(pcaData, Group == "B"),
                  nudge_x       = 30 - subset(pcaData, Group == "B")$PC1,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 0) +
  geom_text_repel(aes_string(label = "name", color = "Group"), size = 6,
                  data          = subset(pcaData, Group == "C"),
                  nudge_x       = -10 - subset(pcaData, Group == "C")$PC1,
                  segment.size  = 0.2,
                  segment.color = "grey50",
                  direction     = "y",
                  hjust         = 0) 
#
pca <- prcomp(t(find_ANOVA)) # compute PCA on transposed KO_sample_counts
pca

# scree plot
plot(pca, type ="line", main="Scree Plot for KO_sample_counts") 

# 從pca中取出標準差(pca$sdev)後再平方，計算variance(特徵值)
vars <- (pca$sdev)^2  
vars

# 計算每個主成分的解釋比例 = 各個主成分的特徵值/總特徵值
props <- vars / sum(vars)    
props

# PCA
autoplot(prcomp(t(normalized_counts)), data = info, 
         colour = 'Group', 
         size = 3, label.repel = 1)

######## 2.4 Lasso ######## 
sele <- normalized_counts[c(A_BCsign, BCsign), ]
Lc <- cbind(info, t(sele))
Lc <- Lc[, -c(2,3,4)]

# Lc$Group <- ifelse(Lc$Group=="A", "A", "BC")
#### Matrix 
x <- model.matrix(Group~ .
                  , data = Lc)[, -1]

y <- as.numeric(Lc$Group)
# y <- ifelse(Lc$Group=="A", 1, 0)

ames_ridge <- glmnet(
  x = x,
  y = y,
  alpha = 1,
  standardize = FALSE,
  nlambda = 10000,
  family = "multinomial"
)

gene <- coef(ames_ridge, s=c(ames_ridge$lambda[1000]))
print(ames_ridge)[1000, 2]
gene[[1]]>0
which(gene[,2]>0)

plot(ames_ridge, xvar = "lambda")


######## 2.5 Multinomial ####
sele <- normalized_counts[c(A_BCsign, BCsign), ]
Lc <- cbind(info, t(sele))
Lc <- Lc[, -c(2,3,4)]

Lc$Group <- relevel(ml$prog, ref = "academic")
test <- multinom(Group ~ ., data = Lc)
summary(test)
######## 2.- MANOVA ########

library(dplyr)

TRY <- t(normalized_counts)
c_TRY <- cbind(TRY, info)

# 沒有inverse matrix, error
# lm_model <- lm(cbind(Feature_A, Feature_B, Feature_C) ~ . - Group,
#                data = c_TRY)
# fit_MANOVA <- Manova(lm_model)

model <- manova(cbind(Feature_A, Feature_B, Feature_C) ~ . - Group,
                data = c_TRY)
summary(model)
model 


# Multinomial ####
sele <- normalized_counts[c(A_BCsign, BCsign), ]
Lc <- cbind(info, t(sele))
Lc <- Lc[, -c(2,3,4)]

Lc$Group <- relevel(ml$prog, ref = "academic")
test <- multinom(Group ~ ., data = Lc)
summary(test)





######## 22 LC ########
counts = readRDS("~/Downloads/LC_sample_counts.rds")
info = readRDS("~/Downloads/LC_sample_info.rds")
dds = DESeqDataSetFromMatrix(countData = counts, 
                             colData = info, 
                             design = ~ Group)
normalized_counts = assay(vst(dds))
Lc <- cbind(info, t(normalized_counts))

trinary <- data.frame(group= c("A", "B", "C"))
i <- 6
for (i in 1:15) {
  start <- (i-1)*100 + 1
  end <- (i)*100
  trinary[start:end, 2:10001] <- Lc[i, 5:10004]
  
  endA <- start + Lc$Feature_A[i]-1
  endB <- endA + Lc$Feature_B[i]
  endC <- endB + Lc$Feature_C[i]
  
  if (Lc$Feature_A[i]!= 0) trinary[start:endA, 1] <- "A"
  if (Lc$Feature_B[i]!= 0) trinary[(endA+1):endB, 1] <- "B"
  if (Lc$Feature_C[i]!= 0) trinary[(endB+1):endC, 1] <- "C"
}



######## 22.1 Lasso ######## 
x <- model.matrix(group~ .
                  , data = trinary)[, -1]

y <- as.numeric(trinary$group)
y <- trinary$group
# y <- ifelse(Lc$Group=="A", 1, 0)

ames_ridge <- glmnet(
  x = x,
  y = y,
  alpha = 1,
  standardize = FALSE,
  nlambda = 1000,
  family = "multinomial",
  type.multinomial = "grouped"
)

TRY <- glmnet(
  x = x,
  y = y,
  alpha = 1,
  standardize = FALSE,
  family = "multinomial",
  type.multinomial = "grouped"
)

G1 <- as.data.frame(gene$`2`[1:10001])
print(ames_ridge)[1000, ]

plot(ames_ridge, xvar = "lambda")
plot(TRY, xvar = "lambda")

# CV
ames_ridge_cv <- cv.glmnet(
  x = x,
  y = y,
  alpha = 1,
  nlambda = 1000,
  family = "multinomial", 
  type.multinomial = "grouped",
  type.measure = "class"
)

plot(ames_ridge_cv)
best.lambda <- ames_ridge_cv$lambda.min
best.lambda <- ames_ridge_cv$lambda.1se

{
  par(mfcol = c(2, 2))
  plot(ames_ridge, xvar = "lambda") 
  abline(v = log(best.lambda), col = "red", lty = "dashed")
  par(mfcol = c(1, 1))
}

# selection
A <- as.numeric(unlist(coef(ames_ridge_cv, s = "lambda.min"))[[1]])[-1]
B <- as.numeric(unlist(coef(ames_ridge_cv, s = "lambda.min"))[[2]])[-1]
C <- as.numeric(unlist(coef(ames_ridge_cv, s = "lambda.min"))[[3]])[-1]
# x <- lapply(coef(ames_ridge_cv, s = "lambda.min")[[1]], as.numeric)

select.ind <- c(which(A != 0), which(B != 0), which(C != 0))
# select.ind = select.ind[-1]-1 # remove `Intercept` and 平移剩下的in
importantRNA <- unique(select.ind)
length(importantRNA)

######## 22.1 Lasso- Predict ####
# managenment data
sele <- normalized_counts[(importantRNA), ]
Lc <- cbind(info, t(normalized_counts))
Lc <- Lc[ ,c(1, importantRNA+4)]
Lc <- Lc[ ,-c(2,3,4)]
# new <- as.matrix(Lc[ ,(importantRNA+4)])
x <- model.matrix(Group~ .
                  , data = Lc)[, -1]

y <- Lc[, c(2,3,4)]

# Predict
TRY <- as.matrix(t(normalized_counts))
prd15 <- predict(ames_ridge_cv, s= "lambda.min", 
                 newx= x, type= "response")

options(digits = 3)
prd15_cbind <- as.data.frame(prd15[,,1]*100) %>% cbind(info)

######## 22.1 Lasso- type= deviance ######## 
x <- model.matrix(group~ .
                  , data = trinary)[, -1]

y <- as.numeric(trinary$group)
y <- trinary$group
# y <- ifelse(Lc$Group=="A", 1, 0)

ames_ridge <- glmnet(
  x = x,
  y = y,
  alpha = 1,
  standardize = FALSE,
  nlambda = 1000,
  family = "multinomial",
  type.multinomial = "grouped"
)

print(ames_ridge)[1000, ]

plot(ames_ridge, xvar = "lambda")

# CV
ames_auc_cv <- cv.glmnet(
  x = x,
  y = y,
  alpha = 1,
  nlambda = 1000,
  family = "multinomial", 
  type.multinomial = "grouped",
  type.measure = "auc"
)

plot(ames_auc_cv)
best.lambda_deviance <- ames_ridge_cv$lambda.min
best.lambda_deviance_1se <- ames_ridge_cv$lambda.1se

{
  par(mfcol = c(2, 2))
  plot(ames_ridge, xvar = "lambda") 
  abline(v = log(best.lambda_deviance), col = "red", lty = "dashed")
  par(mfcol = c(1, 1))
}

# selection
A <- as.numeric(unlist(coef(ames_ridge_cv, s = "lambda.min"))[[1]])[-1]
B <- as.numeric(unlist(coef(ames_ridge_cv, s = "lambda.min"))[[2]])[-1]
C <- as.numeric(unlist(coef(ames_ridge_cv, s = "lambda.min"))[[3]])[-1]
# x <- lapply(coef(ames_ridge_cv, s = "lambda.min")[[1]], as.numeric)

select.ind <- c(which(A != 0), which(B != 0), which(C != 0))
# select.ind = select.ind[-1]-1 # remove `Intercept` and 平移剩下的in
importantRNA_deviance <- unique(select.ind)
length(importantRNA_deviance)

######## 22.1 Lasso- type= deviance Predict ####
# managenment data
sele <- normalized_counts[(importantRNA_deviance), ]
Lc <- cbind(info, t(normalized_counts))
Lc <- Lc[ ,-c(2,3,4)]
# new <- as.matrix(Lc[ ,(importantRNA+4)])
x <- model.matrix(Group~ .
                  , data = Lc)[, -1]

y <- Lc[, c(2,3,4)]

# Predict
prd15_deviance <- predict(ames_auc_cv, s= "lambda.min", 
                 newx= x, type= "response")

options(digits = 3)
prd15_deviance_cbind <- as.data.frame(prd15_deviance[,,1]*100) %>% cbind(info)

######## 22.1 Lasso- CV ####
# fold
pre_CV <- info
pre_CV$A <- 0
pre_CV$B <- 0
pre_CV$C <- 0

RNA_CV <- data.frame(Lc= row.names(info))
pre_CV100 <- pre_CV
RNA_CV100 <- RNA_CV

Lc <- cbind(info, t(normalized_counts))
Lc <- Lc[ ,-c(2,3,4)]
Gx <- model.matrix(Group~ .
                  , data = Lc)[, -1]

i <-11
for (i in 5:15) {
  start <- (i-1)*100 + 1
  end <- (i)*100
  train <- trinary[-c(start:end), ]
  
  x <- model.matrix(group~ .
                    , data = train)[, -1]
  y <- train$group
  
  # Fit LASSO
  ames_ridge <- glmnet(
    x = x,
    y = y,
    alpha = 1,
    standardize = FALSE,
    # nlambda = 1000,
    family = "multinomial",
    type.multinomial = "grouped"
  )
  
  # cv.lasso to find the best lambda
  ames_ridge_cv <- cv.glmnet(
    x = x,
    y = y,
    alpha = 1,
    # nlambda = 1000,
    family = "multinomial", 
    type.multinomial = "grouped",
    type.measure = "deviance"
  )
  
  # selecte the RNA which is sign.
  select.ind <- c(which(A != 0), which(B != 0), which(C != 0))
  # select.ind = select.ind[-1]-1 # remove `Intercept` and 平移剩下的in
  important <- unique(select.ind)
  RNA_CV100[i, 1:length(important)+1] <- important
  
  # Predict
  prd <- predict(ames_ridge_cv, s= "lambda.min", 
                   newx= Gx[i:(i+1),], type= "response")
  
  options(digits = 3)
  pre_CV100$A[i] <- prd[1,1,1]*100
  pre_CV100$B[i] <- prd[1,2,1]*100
  pre_CV100$C[i] <- prd[1,3,1]*100
  print(i)

  }

######## 22.2 Multinomial CV####
muti_pre_CV <- info
muti_pre_CV$A <- 0
muti_pre_CV$B <- 0
muti_pre_CV$C <- 0

Lc <- cbind(info, t(normalized_counts))
Lc <- Lc[ ,-c(2,3,4)]
Gx <- model.matrix(Group~ .
                   , data = Lc)[, -1]

sele <- trinary[, c(1, c(A_BCsign, BCsign)+1)] #1500*15
pre <- as.data.frame(Gx[, c(A_BCsign, BCsign)])

i=1
for(i in 1:15){
  start <- (i-1)*100 + 1
  end <- (i)*100
  train_M <- sele[-c(start:end), ]
  test_M <- pre[i, ]
  model <- multinom(group ~ ., data = train_M)
  temp <- predict(model, newdata= test_M, type= "probs")
  
  muti_pre_CV$A[i] <- round(temp[1]*100)
  muti_pre_CV$B[i] <- round(temp[2]*100)
  muti_pre_CV$C[i] <- round(temp[3]*100)
}



######## 22.3 KM ####
find_ANOVA <- normalized_counts[importantRNA_deviance, ]
counts_cor <- cor(find_ANOVA)
pheatmap(counts_cor, annotation = info[, -1])







