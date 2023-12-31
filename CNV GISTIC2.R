#####CNV GISTIC2#####
#step1 prepare data
#LGG
rm(list=ls())
LGGCNV <- read.table("TCGA-LGG.masked_cnv.tsv.gz", header = TRUE, sep = "\t")
LGG.sam <- substr(LGGCNV$sample,14,16) == "01A"
LGGCNV <- LGGCNV[LGG.sam,]
#GBM
GBMCNV <- read.table("TCGA-GBM.masked_cnv.tsv.gz", header = TRUE, sep = "\t")
GBM.sam <- substr(GBMCNV$sample,14,16) == "01A"
GBMCNV <- GBMCNV[GBM.sam,]

tumorCNV <- rbind(LGGCNV, GBMCNV)


riskScoreTCGA <- read.table("riskScoreTCGA.csv", header = TRUE, sep = ",")[,2:3]

highsample <- riskScoreTCGA[riskScoreTCGA$group == "high",1]
highCNV <- tumorCNV[tumorCNV$sample %in% highsample,-6]
highCNV$key <- result
highCNV <- highCNV[highCNV$key > 0 | is.na(highCNV$key),]



# step2 marker file data download and processing
library(data.table)
markersMatrix <- data.table::fread("snp6.na35.remap.hg38.subset.txt.gz",data.table = F)
str(markersMatrix)
fre <- markersMatrix$freqcnv == "FALSE"
table(fre)
markersMatrix <- markersMatrix[fre, 1:3]
colnames(markersMatrix) <- c("Marker_Name", "Chromosome", "Marker_Position")
#write.table(markersMatrix, file = "marker_file.txt", sep = "\t", row.names = F)


#step3 visualization
library(maftools)
highgistic <- readGistic(gisticAllLesionsFile = "./high1/all_lesions.conf_95.txt",
                         gisticAmpGenesFile = "./high1/amp_genes.conf_95.txt",
                         gisticDelGenesFile = "./high1/del_genes.conf_95.txt",
                         gisticScoresFile = "./high1/scores.gistic",
                         isTCGA = T)


getSampleSummary(highgistic)
genehigh <- getGeneSummary(highgistic)
getCytobandSummary(highgistic)
gisticChromPlot(gistic = highgistic, 
                y_lims = c(-0.15,0.9),
                fdrCutOff = 0.01, 
                markBands = "all",
                txtSize = 0.5)

legend("topright", legend = c("Del","Amp"), fill = c("blue", "red"), box.lwd = 0.05)

