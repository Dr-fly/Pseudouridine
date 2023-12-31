#####Drug Sensitivity####
#prediction
source("ffkkpharTCGA.R")
#x = expr, y = res, z = test
ffkkpharTCGA <- function(x, y){
  library (PharmacoGx)
  library (VennDiagram)
  library (RColorBrewer)
  library (R.utils)
  library (downloader)
  library (gdata)
  library (AnnotationDbi)
  library (hgu133a.db)
  library (hthgu133a.db)
  library (gridExtra)
  library (dplyr)
  
  library(oncoPredict)
  library(car)
  library(sva)
  library(ridge)
  calcPhenotype(trainingExprData = x,
                trainingPtype = y,
                testExprData = TCGAtrial,
                batchCorrect = 'eb',  #   "eb" for ComBat  
                powerTransformPhenotype = TRUE,
                removeLowVaryingGenes = 0.2,
                minNumSamples = 10, 
                printOutput = TRUE, 
                removeLowVaringGenesFrom = 'rawData' )
  
  a <- read.csv("./calcPhenotype_Output/DrugPredictions.csv", header = T, row.names = 1, sep = ",")
  loc <- paste("./calcPhenotype_Output/DrugPredictions.",deparse(substitute(x)),"TCGA.csv", sep = "")
  # loc <- paste("./calcPhenotype_Output/DrugPredictions.","GDSC2_Expr","TCGA.csv", sep = "")
  write.csv(a, file = loc)
  
  #screens drugs
  group <- read.csv("groupTCGA.csv", header = T,  sep = ",")
  sensit <- read.csv(loc,header = T,row.names = 1,sep = ",")
  # sensit <- read.csv("./calcPhenotype_Output/DrugPredictions.GDSC2_ExprTCGA.csv",header = T,row.names = 1,sep = ",")
  
  sensit$barcode <- rownames(sensit)
  mix <- merge(sensit,group,  by = "barcode", all.x = TRUE)
  rownames(mix) <- mix[,1]
  mix <- mix[,-1]
  
  # Create an empty table to store the results
  wilcox <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(wilcox) <- c("Drug", "W.statistic", "p.value", "direction")
  
  # Loop through each drug
  for (i in 1:(ncol(mix)-1)) {
    # Calculate the mean of the high and low groups
    high_mean <- mean(mix[mix$group == "high", i], na.rm = TRUE)
    low_mean <- mean(mix[mix$group == "low", i], na.rm = TRUE)
    
    # Determine the direction of high-low and assign a value
    if (high_mean > low_mean) {
      direction <- "resistant"
    } else {
      direction <- "sensitive"
    }
    
    result <- wilcox.test(mix[,i] ~ mix[,ncol(mix)], data = mix)
    
    # Store the results in a table
    wilcox[nrow(wilcox) + 1,] <- c(colnames(mix)[i], result$statistic, result$p.value, direction)
  }
  wilcox$p.value <- as.numeric(wilcox$p.value)
  wilcox_filter <- wilcox[wilcox$p.value < 0.05,]
}

ffkkpharTCGA(GDSC2_Expr, GDSC2_Res)
ffkkpharTCGA(GBM_Expr, GBM_Res)

#correlation
lnc_Expr <- read.csv("psilncTCGA.csv", header = T, row.names = 1)
rownames(lnc_Expr) <- lnc_Expr[,1]
lnc_Expr <- lnc_Expr[,13:ncol(lnc_Expr)]
lnc_Expr <- lnc_Expr[,-1]

TCGAlnc <- t(lnc_Expr)
TCGAlnc <- as.matrix(TCGAlnc)
TCGAsig <- TCGAlnc[c("DNAJC27.AS1","DNMBP.AS1","ZBTB20.AS4","GDNF.AS1"),]
TCGAsig <- t(TCGAsig)

#Extract predicted value
GDSC2_TCGAp <- read.csv("./calcPhenotype_Output/DrugPredictions.GDSC2_ExprTCGA.csv", header = T, row.names = 1)
GBM_TCGAp <- read.csv("./calcPhenotype_Output/DrugPredictions.GBM_ExprTCGA.csv", header = T, row.names = 1)
GDSC2_TCGAf <- GDSC2_TCGAp[rownames(GDSC2_TCGAp)%in%rownames(TCGAsig), colnames(GDSC2_TCGAp) %in% ab]
GBM_TCGAf <- GBM_TCGAp[rownames(GBM_TCGAp)%in%rownames(TCGAsig),colnames(GBM_TCGAp) %in% ab]

#Computational correlation
#GDSC TCGA
corGDSC2_TCGA <- data.frame(matrix(NA, nrow = 12, ncol = 4))
colnames(corGDSC2_TCGA) <- c("Drug", "lncRNA", "Spearman's rho", "p-value")
row0 <- 0
for (i in 1:4) {
  for (j in 1:3) {
    corGDSC2_TCGA[j+row0,1] <- colnames(TCGAsig)[i]
    corGDSC2_TCGA[j+row0,2] <- colnames(GDSC2_TCGAf)[j]
    corGDSC2_TCGA[j+row0,3] <- cor.test(TCGAsig[,i], GDSC2_TCGAf[,j], method = "spearman")$estimate
    corGDSC2_TCGA[j+row0,4] <- cor.test(TCGAsig[,i], GDSC2_TCGAf[,j], method = "spearman")$p.value
  }
  row0 <- i*3
}
corGDSC2_TCGAf <- corGDSC2_TCGA[corGDSC2_TCGA[,4]<0.05,]


# GBM TCGA
corGBM_TCGA <- data.frame(matrix(NA, nrow = 12, ncol = 4))
colnames(corGBM_TCGA) <- c("Drug", "lncRNA", "Spearman's rho", "p-value")
row0 <- 0
for (i in 1:4) {
  for (j in 1:3) {
    corGBM_TCGA[j+row0,1] <- colnames(TCGAsig)[i]
    corGBM_TCGA[j+row0,2] <- colnames(GBM_TCGAf)[j]
    corGBM_TCGA[j+row0,3] <- cor.test(TCGAsig[,i], GBM_TCGAf[,j], method = "spearman")$estimate
    corGBM_TCGA[j+row0,4] <- as.numeric(cor.test(TCGAsig[,i], GBM_TCGAf[,j], method = "spearman")$p.value)
  }
  row0 <- i*3
}
corGBM_TCGAf <- corGBM_TCGA[corGBM_TCGA[,4]<0.05,]



#There's only one drug left, Mitoxantrone.
ab <- "Mitoxantrone"

corGDSC2_TCGA <- data.frame(matrix(NA, nrow = 4, ncol = 3))
colnames(corGDSC2_TCGA) <- c("lncRNA", "Spearman's rho", "p-value")
row0 <- 0
for (i in 1:4) {
  corGDSC2_TCGA[i,1] <- colnames(TCGAsig)[i]
  # corGDSC2_TCGA[j+row0,2] <- colnames(GDSC2_TCGAf)[j]
  corGDSC2_TCGA[i,2] <- cor.test(TCGAsig[,i], GDSC2_TCGAf, method = "spearman")$estimate
  corGDSC2_TCGA[i,3] <- cor.test(TCGAsig[,i], GDSC2_TCGAf, method = "spearman")$p.value
}
corGDSC2_TCGAf <- corGDSC2_TCGA[corGDSC2_TCGA[,3]<0.05,]


corGBM_TCGA <- data.frame(matrix(NA, nrow = 4, ncol = 3))
colnames(corGBM_TCGA) <- c("lncRNA", "Spearman's rho", "p-value")
row0 <- 0
for (i in 1:4) {
  corGBM_TCGA[i,1] <- colnames(TCGAsig)[i]
  # corGDSC2_TCGA[j+row0,2] <- colnames(GDSC2_TCGAf)[j]
  corGBM_TCGA[i,2] <- cor.test(TCGAsig[,i], GBM_TCGAf, method = "spearman")$estimate
  corGBM_TCGA[i,3] <- cor.test(TCGAsig[,i], GBM_TCGAf, method = "spearman")$p.value
}
corGBM_TCGAf <- corGBM_TCGA[corGBM_TCGA[,3]<0.05,]



