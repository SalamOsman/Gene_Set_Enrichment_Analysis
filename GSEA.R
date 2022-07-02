---
title: "GSEA"
author: "Abdus Salam Khan"
date: '2022-04-27'
output: html_document
---

# Loading the required libraries for GSE analysis.

library(tidyverse)
library(biomaRt)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cowplot)
library(ReactomePA)
library(ggnewscale)


# Importing the differential expressed genes list.

DGE_list <- read.csv("C:/Users/khan_/Desktop/GSEA/DEGs.csv", header = T, sep = ",")
head(DGE_list)


# Renaming column header in the data file.

colnames(DGE_list)[1] <- "Gene_symbol"
head(DGE_list)


Remember: You will be needing to filter the gene list based on adjusted p-value or padj < 0.05. In the current list, the genes have already been filtered by setting a cutoff. If you need to screen your genes list then, run the following step;
"DGE_list <- filter(DGE_list, DGE_list$padj < 0.05)"

# Annotating Gene symbols with Entrez gene ID in DGE_list file. For this step, we need to download Ensembl repository.

ensembl = useMart("ensembl")


# Look at the different organisms libraries has been stored in Ensembl repository. We'll be using "hsapiens_gene_ensembl" library because the gene symbols belong to human genome. It is Ok if some of the Gene_symbols misses Entrez IDs.

#listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)


# Retrieving the IDs column "entrez gene ids" from the Ensembl database by setiing an option to 'entrezgene_id' as a query while our gene list symbols were in 'external_gene_name' format. 
# Always look for the id format in the imported gene list you are currently exploring.   
and converting the to entrez gene ids.

IDs <- getBM(attributes=c('external_gene_name','entrezgene_id'), 
             filters = 'external_gene_name', values = DGE_list$Gene_symbol, mart = ensembl)


# Removing duplicate IDs from the files.

IDs <- IDs[!duplicated(IDs$external_gene_name), ]
DGE_list <- DGE_list[!duplicated(DGE_list$Gene_symbol), ]


# Renaming the column name for merging with the dataset.

names(IDs)[1] <- 'Gene_symbol'
DGE_list <- dplyr::inner_join(DGE_list, IDs, by="Gene_symbol")


# Removing rows with NA in IDs column.

DGE_list <- subset(DGE_list, !is.na(DGE_list$entrezgene_id))


# Preparing a gene lists file for enrichment analysis..

genelist <- cbind(DGE_list$entrezgene_id, DGE_list$log2FoldChange)
colnames(genelist) <- c("ID", "Log2FC")
genelist <- as.data.frame(genelist)
geneList <- as.numeric(genelist[,2])
names(geneList) <- as.character(genelist[,1])
geneList <- sort(geneList, decreasing = TRUE)
head(geneList)


This step is very crucial for setting a cutoff. I will choose value greater than absolute 1. You can play with this value to get an optimum results. If you need to know any further, then start learning "log2FC" of genes in transcriptomics profiling. 
# log2FC mean 2 fold increase in gene expression to regulate either in positive abd negative direction.  

gene <- names(geneList)[abs(geneList) > 1]


# GSEA using Mutational signature database (MsigDb) and getting H sapiens gene sets. The steps includes collecting from Msigdb for human, extraction of hallmarks info. and then retrieving gene sets for each hallmark.

msigdbr_species()
m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame
m_df_H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
head(m_df_H_t2g)
collections <- msigdbr_collections()
m_df_H <- m_df[m_df$gs_cat=="H",]
head(m_df_H)
msigdbr_t2g = m_df_H %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()


# Hallmarks enrichment analysis.

hallmark_enricher <- enricher(gene = gene, TERM2GENE = msigdbr_t2g)
res_dot <- dotplot(hallmark_enricher, showCategory=10, font.size = 8, title = "Dotplot of enriched hallmarks in my study")
res_bar <- barplot(hallmark_enricher, showCategory=10, font.size = 8, title = "Barplot of enriched hallmarks in my study")
res_bar
res_dot


# Plotting a enrichment network. 

cnetplot(hallmark_enricher, categorySize="pvalue", foldChange=geneList, cex_label_gene = 0.35, cex_label_category = 0.75, shadowtext = 'none', color_category = "gold")


# GO enrichment analysis. Look at the different option you have. You can gather information about this function by typing "help(gseGO)" in the command line. 

GO_gse <- gseGO(geneList=geneList, 
               ont ="ALL",#Specify if else; BP, MF, CC 
               keyType = "ENTREZID", 
               nPerm = 1000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Hs.eg.db, 
               pAdjustMethod = "none")

dotplot(GO_gse, showCategory=10, font.size = 8, split=".sign", title = "GO terms enriched in this study") + facet_grid(.~.sign)


# Reactome pathways enrichment analysis.

reactome_enriched <- enrichPathway(gene=gene, pvalueCutoff=0.05, readable=T)

bar_reactome <- barplot(reactome_enriched, showCategory=10, font.size = 10,  title = "Barplot of enriched Reactome terms in this study", col = 'pvalue')

dot_reactome <- dotplot(reactome_enriched, showCategory=10, font.size = 10, title = "Dotplot of enriched Reactome terms in this study")
bar_reactome
dot_reactome

