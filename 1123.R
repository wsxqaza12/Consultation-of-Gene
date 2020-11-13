######## 1 Imput data ######## 
counts = readRDS("data/LC_sample_counts.rds")
info = readRDS("data/LC_sample_info.rds")
dds = DESeqDataSetFromMatrix(countData = counts, 
                             colData = info, 
                             design = ~ Group)
normalized_counts = assay(vst(dds))
Lc <- cbind(info, t(normalized_counts))

binary <- data.frame(group= c("A", "BC"))

for (i in 1:15) {
  start <- (i-1)*100 + 1
  end <- (i)*100
  binary[start:end, 2:10001] <- Lc[i, 5:10004]
  
  endA <- start + Lc$Feature_A[i]-1
  endB <- endA + Lc$Feature_B[i]
  endC <- endB + Lc$Feature_C[i]
  
  if (Lc$Feature_A[i]!= 0) binary[start:endA, 1] <- "A"
  if (Lc$Feature_B[i]!= 0) binary[(endA+1):endB, 1] <- "BC"
  if (Lc$Feature_C[i]!= 0) binary[(endB+1):endC, 1] <- "BC"
}


######## 2 Lasso- type= deviance ######## 
x <- model.matrix(group~ .
                  , data = binary)[, -1]

y <- as.numeric(binary$group)
y <- binary$group

ames_cv_bi <- cv.glmnet(
  x = x,
  y = y,
  alpha = 1,
  nlambda = 1000,
  family = "multinomial", 
  # type.multinomial = "grouped",
  type.measure = "deviance"
)

# plot(ames_cv_bi)
best.lambda_deviance <- ames_cv_bi$lambda.min
best.lambda_deviance_1se <- ames_cv_bi$lambda.1se

{
  par(mfcol = c(2, 2))
  plot(ames_ridge_bi, xvar = "lambda") 
  abline(v = log(best.lambda_deviance), col = "red", lty = "dashed")
  par(mfcol = c(1, 1))
}

# selection
A_bi <- as.numeric(unlist(coef(ames_cv_bi, s = "lambda.min"))[[1]])[-1]
B_bi <- as.numeric(unlist(coef(ames_cv_bi, s = "lambda.min"))[[2]])[-1]

# x <- lapply(coef(ames_cv_bi, s = "lambda.min")[[1]], as.numeric)

select.ind <- c(which(A_bi != 0), which(B_bi != 0))
# select.ind = select.ind[-1]-1 # remove `Intercept` and 平移剩下的in
importantRNA_deviance_bi <- unique(select.ind)
length(importantRNA_deviance_bi)

######## 3 Lasso- type= deviance Predict ####
# managenment data
sele <- normalized_counts[(importantRNA_deviance_bi), ]
Lc <- cbind(info, t(normalized_counts))
Lc <- Lc[ ,-c(2,3,4)]
# new <- as.matrix(Lc[ ,(importantRNA+4)])
x <- model.matrix(Group~ .
                  , data = Lc)[, -1]

y <- Lc[, c(2,3,4)]

# Predict
prd15_deviance_bi <- predict(ames_cv_bi, s= "lambda.min", 
                          newx= x, type= "response")

options(digits = 3)
prd15_deviance_cbind_bi <- as.data.frame(prd15_deviance_bi[,,1]*100) %>% cbind(info)

######## 4 Multinomial CV####
muti_pre_CV_bi <- info
muti_pre_CV_bi$A <- 0
muti_pre_CV_bi$BC <- 0


sele <- binary[, c(1, importantRNA_deviance_bi+1)] #1500*15
pre <- as.data.frame(Gx[, (importantRNA_deviance_bi)])
i=1

for(i in 1:15){
  start <- (i-1)*100 + 1
  end <- (i)*100
  train_M <- sele[-c(start:end), ]
  test_M <- pre[i, ]
  model <- multinom(group ~ ., data = train_M)
  temp <- predict(model, newdata= test_M, type= "probs")
  
  options(digits = 6)
  muti_pre_CV_bi$A[i] <- round(100- temp[1]*100)
  muti_pre_CV_bi$BC[i] <- round(temp[1]*100)
}

muti_pre_CV_bi %>% mutate(Feature_B = Feature_B+ Feature_C) %>%
  select(- Feature_C)
compair <- multinom(group ~ ., data = sele)
######## 5 KM ####
find_ANOVA <- normalized_counts[importantRNA_deviance, ]
counts_cor <- cor(find_ANOVA)
pheatmap(counts_cor, annotation = info[, -1])

######## 6 DirichletReg ####
library(DirichletReg)
Lc <- cbind(info, t(normalized_counts))
Di_data <- as.data.frame(Gx[, (importantRNA_deviance_bi)])

info = readRDS("data/LC_sample_info.rds")
precentLC <- DR_data(info[, 2:4]/100)
precentLC[,1] <- info[, 2]/100
precentLC[,2] <- info[, 3]/100
precentLC[,3] <- info[, 4]/100
precentLC

Di_data$prob <- DR_data(precentLC[, c('Feature_A', 'Feature_B', 'Feature_C')])

mod2 <- DirichReg(prob ~
                    gene_307 + gene_1485,
                  data = Di_data, model = "common")

predict(mod2, newdata = Di_data)




