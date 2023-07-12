#####Enrichment Functions#####
plotEnrichment<-function(deres,name){
  ##prepare genelist for enrichment analysis
  gene_list = deres$logFC
  names(gene_list) <- rownames(deres)
  gene_list <- sort(gene_list, decreasing = T)
  ##enrichment analysis
  library(clusterProfiler)
  goEnsembl <-
    read.gmt("E:/R/TCGA2.0/sourceData/c5.all.v7.5.1.symbols.gmt")
  gse <-
    GSEA(
      gene_list,
      TERM2GENE = goEnsembl,
      eps = 0,
      pAdjustMethod = "BH",
      nPermSimple = 10000
    )
  termNum <- 5
  kkUp = gse[gse$NES > 0,]
  showTerm = row.names(kkUp)[1:termNum]
  pdf(paste0("GseaUpEnrich-",name,".pdf"),
      width = 5,
      height = 4)
  p<-gseaplot2(
    gse,
    geneSetID = showTerm,
    base_size = 8,
    title = "Enriched in high risk group",
    color = mypal[1:5],
    ES_geom = "dot"
  )
  print(p)
  dev.off()
  rm(kkUp)
  kkDown = gse[gse$NES < 0,]
  showTerm = row.names(kkDown)[1:termNum]
  pdf(paste0("GseaDownEnrich-",name,".pdf"),
      width = 5,
      height = 4)
  p<-gseaplot2(
    gse,
    geneSetID = showTerm,
    base_size = 8,
    title = "Enriched in low risk group",
    color = mypal[1:5],
    ES_geom = "dot"
  )
  print(p)
  dev.off()
  rm(kkDown)
  ENTID <- geneNM %>%
    filter(geneNM$external_gene_name %in% rownames(deres)) %>%
    select(external_gene_name, entrezgene_id)
  ENTID <- ENTID[!duplicated(ENTID$external_gene_name), ]
  ENTID <-
    merge.data.frame(
      x = data.frame("external_gene_name" = rownames(deres)),
      y = ENTID,
      all.x = T
    )
  ENTID <- ENTID[order(ENTID$external_gene_name),]
  deres <- deres[order(rownames(deres)),]
  identical(ENTID$external_gene_name,rownames(deres))#
  deres$entrezgene_id<-ENTID$entrezgene_id#Add entrez geneID to deres dataframe
  #Select Top 2000 most differentially expressed genes
  genes <-
    deres$entrezgene_id[order(abs(deres$logFC), decreasing = T)][1:2000]
  genes <- genes[!duplicated(genes)]
  genes <- genes[!is.na(genes)]
  Universe <- deres$entrezgene_id[!duplicated(deres$entrezgene_id)]
  Universe <- Universe[!is.na(Universe)]
  p <- enrichGO(
    gene          = genes,
    universe      = as.character(Universe),
    keyType = "ENTREZID",
    OrgDb         = "org.Hs.eg.db",
    ont           = "ALL",
    pool = T,
    readable = T
  )
  pdf(paste0("Go-",name,".pdf"),
      width = 7,
      height = 6)
  p<-dotplot(p, split = "ONTOLOGY", showCategory = 5) + facet_grid(ONTOLOGY ~
                                                                     ., scale = "free")
  print(p)
  dev.off()
  p <-
    enrichKEGG(
      genes,
      organism = "hsa",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.1,
      use_internal_data = T
    )
  pdf(paste0("KEGG-",name,".pdf"),
      width = 6,
      height = 4)
  p<-dotplot(p)
  print(p)
  dev.off()
}