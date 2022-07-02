# Gene Set Enrichment Analysis (GSEA)

GSEA is a short and more interesting method as compared to looking for the biological function of multiple genes. Your list of genes could be either derived from differential gene or miRNAs or DNA methylation expression entities which result in log2FC score. The value of log2FC determines the direction of expression in a sample study.

## Principles of operation:
In GSEA, our list of multiple genes (look for the format in "DEGs.csv") were matched using knowledge-based method (using already trained/prepared datasets). In the current workflow, we used trained datasets from Molecular Signatures Database (MSigDB), Reactome Database, and Gene Ontology (GO) Database. The name of the main package in this workflow is "clusterProfiler". Alternative methods are availbe in DAVID, g:Profiler servers etc. 

## Intructions when using the scipt:
  * Make sure you have internet connection.
  * You need to know the name of the id format of genes in your list. 
  * Apply the cutoff for log2FC and adjusted p value (adjP) accordingly.
  * Remember to use the proper host id when performing the enrichment analysis.

## Contact
abdussalam@precisionmedicine.pk




