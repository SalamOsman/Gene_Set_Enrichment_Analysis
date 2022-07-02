# Gene Set Enrichment Analysis (GSEA)

GSEA is a short and more interesting method as compared to looking for the biological function of multiple genes. Your list of genes could be either derived from differential gene or miRNAs or DNA methylation expression entities which result in log2FC score. The value of log2FC determines the direction of expression in a sample study.

## Principles:
In GSEA, our list of multiple genes (look for the format in "DEGs.csv") were matched using knowledge-based method (using already trained/prepared datasets). In the current workflow, we used trained datasets from Molecular Signatures Database (MSigDB), Reactome Database, and Gene Ontology (GO) Database. The name of the main package in this workflow is "clusterProfiler". Alternative methods

## Steps:
(1) Make sure you have internet connections.
(2) Start loading the libraries and import your data that you prepared in the study.
(3) Also make sure that your IDs are in entrezgene id format, if not then make sure the name of your ids by typing "listAttributes(ensembl)"


## Intructions when using the scipt:
(1) Make sure you have internet connection.
(2) 


