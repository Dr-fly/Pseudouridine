#####Differential expression analysis with EdgeR#####
#Batch effect removal with SVA package
batchEffectRemoval <- function(dataset, dataRef) {
  require(sva)
  genesInCommon <-
    intersect(rownames(dataset), rownames(dataRef))
  dataset <- dataset[rownames(dataset)%in%genesInCommon, ]
  dataRef <- dataRef[rownames(dataRef)%in%genesInCommon, ]
  batch<-c(rep(1, length(dataset[1, ])), rep(2, length(dataRef[1, ])))
  combat_edata3 = ComBat(cbind(dataset, dataRef),
                         batch = batch,
                         ref.batch = 2)
  dataset <- combat_edata3[, 1:length(dataset[1, ])]
  unloadNamespace("sva")
  unloadNamespace("mgcv")#unload to avoid conflict with mlr3 packages
  dataset <- as.data.frame(dataset)
}
batchEffectRemoval_Seq <- function(dataset, dataRef) {
  require(sva)
  genesInCommon <-
    intersect(rownames(dataset), rownames(dataRef))
  dataset <- dataset[rownames(dataset)%in%genesInCommon, ]
  dataRef <- dataRef[rownames(dataRef)%in%genesInCommon, ]
  batch<-c(rep(1, length(dataset[1, ])), rep(2, length(dataRef[1, ])))
  dataset = ComBat_seq(cbind(dataset, dataRef),
                       batch = batch,full_mod=TRUE)
  unloadNamespace("sva")
  unloadNamespace("mgcv")
  dataset <- as.data.frame(dataset)
}
#Functions used to Convert rownames with ensebmle id to gene symbol
library(biomaRt)
ensembl <-
  useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
geneNM <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "transcript_biotype",
    "external_gene_name"
  ),
  filters    = "ensembl_gene_id",
  values     = rownames(dataset),
  mart       = ensembl
)
convertGeneName <- function(dataset) {
  temp <-
    merge(
      data.frame("ensembl_gene_id" = rownames(dataset)),
      data.frame(
        "ensembl_gene_id" = geneNM$ensembl_gene_id,
        geneNM$external_gene_name
      )
    )
  dataset <- dataset[order(rownames(dataset)), ]
  dataset <-
    dataset[rownames(dataset) %in% temp$ensembl_gene_id, ]
  temp <- temp[order(temp$ensembl_gene_id), ]
  keep <- !duplicated(temp$geneNM.external_gene_name)
  dataset <- dataset[keep, ]
  temp <- temp[keep, ]
  identical(rownames(dataset), temp$ensembl_gene_id)
  rownames(dataset) <- temp$geneNM.external_gene_name
  dataset <- dataset[!duplicated(rownames(dataset)), ]
}
#Differential expression analysis
library(edgeR)
y <- DGEList(counts = dt, group = group,remove.zeros = T)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
rm(dt)
gc()
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y,)
y <- estimateDisp(y, design, robust = T)
con <- makeContrasts(TCGA-GTEX, levels = design)
fit <- glmQLFit(y, design, robust = T)
qlf <- glmQLFTest(fit, contrast = con)
deres <- cbind(qlf$table, FDR = p.adjust(qlf$table$PValue,
                                         "fdr"))


#Volcano Plot
EnhancedVolcano::EnhancedVolcano(
  deres,
  lab = rownames(deres),
  x = 'logFC',
  y = 'FDR',
  subtitle = "Tumor Versus Normal"
)
#Heatmap
cal_z_score <- function(x) {
  (x - mean(x)) / sd(x)
}
data <-
  cbind(dataGTEXVST[rownames(dataGTEXVST) %in% ListOverlap,], dataTCGAVST[rownames(dataTCGAVST) %in% ListOverlap,])

mf1 <- data.frame(Type = ifelse(substr(colnames(data), 1, 4) == "TCGA", "Tumor", "Normal"))
rownames(mf1) <- colnames(data)
data <-
  t(apply(data, 1, cal_z_score))
mycolor <-
  colorRamp2(c(-2, 0, 2), c("#3B4992FF", "#EEEEEE", "#E64B35FF"))
pdf(
  file = "HMTCGAVersusGTEX.pdf",
  width = 9,
  height = 5.0,
  onefile = F
)
ht = ComplexHeatmap::pheatmap(
  data,
  annotation_col = mf1,
  show_rownames = T,
  show_colnames = F,
  annotation_names_row = F,
  cluster_rows = F,
  cluster_cols = F,
  fontsize_col = 5,
  color = mycolor,
  annotation_colors = ann_colors,
  name = "Z score"
)
print(ht)
dev.off()
#Lasso regression 
#Lasso-Cox Regression with 100 lambda
cv_fit_cox <-
  cv.glmnet(
    x = x.train,
    y = y.train,
    nlambda = 100,
    alpha = 1,
    nfolds = 50,
    family = "cox",
    type.measure = "C"
  )
#Plots
plot(cv_fit_cox)
plot(cv_fit_cox$glmnet.fit)
#StepWise Cox regression
StepCoxModel <-
  stepwiseCox(
    fml,
    data = lassodataTCGA,
    selection = "bidirection",
    select = "HQ",
    method = c("efron"),
    sle = 0.01,
    sls = 0.01
  )
#Foreast Plot
survminer::ggforest(
  FinalMultiCoxModel,
  data = lassodataTCGA,
  main = "Hazard ratio",
  fontsize = 1.0,
  refLabel = "1",
  noDigits = 2
)
#Corr Plot for risk score and risk score
plotRiskScore_Survival<-function(clinicalplotdata,name){
  pdf(
    file = paste0("CorrRisk-Surv-",name,".pdf"),
    width = 6,
    height = 5,
    onefile = F
  )
  plt<-ggplot(clinicalplotdata, aes(x = overall_survival, y = risk)) +
    geom_hline(yintercept = median(clinicalplotdata$risk), color="gray50",size=1,linetype = "dashed")+
    geom_point(
      aes(color = risktest),
      size = 1.4, 
      alpha = 0.8 # It's nice to add some transparency because there may be overlap.
    ) +
    # Use custom colors
    scale_color_manual(
      values = c(mypal[1],mypal[4]),
      name="Risk"
    ) +
    labs(
      x = "Survival Time (Days)",
      y = "Risk Score"
    )+
    scale_y_log10()+
    theme(
      legend.background = element_blank(),
      # This one removes the background behind each key in the legend
      legend.key = element_blank(),
      legend.title = element_text(face = "bold", size=14),
      legend.text = element_text(size = 12),
      # Customize title and subtitle font/size/color
      plot.title.position = "plot",
      # Adjust axis parameters such as size and color.
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 14),
      # Axis lines are now lighter than default
      axis.line = element_line(colour = "black"),
      # Only keep y-axis major grid lines, with a grey color and dashed type.
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      # Use a light color for the background of the plot and the panel.
      panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
      plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4")
    )
  print(plt)
  dev.off()
}
#Mutation WaterFall Plots
dtLow <-
  dt[substr(dt$Tumor_Sample_Barcode, 1, 16) %in% substr(survivalDataTCGA$barcode[survivalDataTCGA$risktest ==
                                                                                   "Low"], 1, 16), ]
dtHigh <-
  dt[substr(dt$Tumor_Sample_Barcode, 1, 16) %in% substr(survivalDataTCGA$barcode[survivalDataTCGA$risktest ==
                                                                                   "High"], 1, 16), ]
length(levels(as.factor(dtHigh$Tumor_Sample_Barcode)))
length(levels(as.factor(dtLow$Tumor_Sample_Barcode)))
duplicated(substr(dtHigh$Tumor_Sample_Barcode, 1, 16))
pdf(
  file = paste0("TCGAWaterfall-Low.pdf"),
  width = 13,
  height = 5
)
waterfall(
  dtLow,
  fileType = "MAF",
  mainRecurCutoff = 0.05,
  maxGenes = 15,
  main_geneLabSize = 14
)
dev.off()
pdf(
  file = paste0("TCGAWaterfall-High.pdf"),
  width = 13,
  height = 5
)
waterfall(
  dtHigh,
  fileType = "MAF",
  mainRecurCutoff = 0.05,
  maxGenes = 15,
  main_geneLabSize = 14
)
dev.off()