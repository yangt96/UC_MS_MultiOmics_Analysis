Requirements | 运行环境

R version: 4.3.3

Key R packages:

limma, edgeR, DESeq2 (DEG)

WGCNA (network analysis)

clusterProfiler, enrichplot (GSEA/functional enrichment)

randomForest，XGB (machine learning)

Seurat, SingleR, Harmony (single-cell analysis)

Installation 

Clone the repository 

git clone https://github.com/yangt96/UC_MS_MultiOmics_Analysis.git
cd UC_MS_MultiOmics_Analysis


Install required R packages 

if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2","edgeR","clusterProfiler","enrichplot","WGCNA","SingleR","Seurat")
install.packages("randomForest"))


Usage 

Data preparation 

Download datasets from GEO (accession numbers listed in the manuscript)
Run DEG analysis 
Rscript scripts/DEG_analysis.R

Run WGCNA 

Rscript scripts/WGCNA_analysis.R

Run GSEA 

Rscript scripts/GSEA_analysis.R


Run Immune analysis 

Rscript scripts/Immune_infiltration.R


Run Machine Learning modeling 

Rscript scripts/ML_modeling.R


Run scRNA-seq analysis 单细胞分析

Rscript scripts/scRNAseq_analysis.R


Results (figures/tables) will be stored in the results/ folder.


Notes 

All datasets are publicly available from GEO/other databases.
This repository is made publicly available to ensure transparency and reproducibility of our analyses.
The repository will remain open.
If access is required beyond this date, please contact the corresponding author.
