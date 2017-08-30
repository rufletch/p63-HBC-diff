# p63-HBC-diff
Code and resources related to the olfactory epithelial HBC stem cell differentiation project (oeHBCdiff)

Below are the R scripts for analyzing the single-cell RNA-seq data from differentiating HBC stem cells of the olfactory epithelium, presented in the following manuscript:

Fletcher RB\*, Das D\*, Gadye L, Street KN, Baudhuin A, Wagner A, Cole MB, Flores Q, Choi YG, Yosef N, Purdom E, Dudoit S, Risso D, Ngai J. Deconstructing Olfactory Stem Cell Trajectories at Single Cell Resolution. Cell Stem Cell (2017). https://doi.org/10.1016/j.stem.2017.04.003 (\* co-first authors)

The data are available on GEO in [GSE95601](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95601).

The repository currently has scripts that take as input Expression Set data and perform a series of computations, interspersed with visualizations. First, the data are filtered for poor quality cells and less informative genes. The data are normalized, and biological contaminants and known doublets (based on co-expression of differentiated cell markers) are removed. Then, the data are re-filtered and re-normalized.

After filtering and normalization, we clustered the data using clusterExperiment, performed developmental ordering and inferred lineage trajectories and branching with slingshot. For each lineage, differentially expressed genes were identified and then clustered to reveal coordinated gene expression. We used Gene Set Enrichment Analysis to infer pathways regulating cell fates and transitions.

We created a number of visualizations based on clustering, experimental condition, and developmental order. We displayed coordinated and correlated differentially expressed genes including transcription factors, as well as a set of cell cycle genes and selected regulators of cell fate transitions along each lineage. The olfactory receptors (OR) and factors associated with OR regulation were plotted along the neuronal lineage. We also presented the top enriched gene sets for each cell cluster. 


### Create output directories and add to .gitignore
In project directory, run `mkdir -p output/{clust,data,gClust,romer,viz,DE,EDA}/oeHBCdiff`, and add new directories to `.gitignore`. Place the scripts in the 'scripts' directory and the initial eSet 'data' in the data directory.

### Filtering and Normalization
`oeHBCdiff_1_filt_norm.sh` performs the following analyses, by calling various R scripts (given in parentheses):

1. Filtering based on technical attributes (`oeHBCdiff_filtering.R`)
2. Normalization using SCONE (`oeHBCdiff_norm.R`)

`oeHBCdiff_2_renorm_clust_devO_DE.sh` performs the following analyses, by calling various R scripts (given in parentheses):

3. Create SummarizedExperiment object for desired normalization (`oeHBCdiff_makeSE.R`)
4. Identification of **biological** contaminants (`oeHBCdiff_exclude.R`)
4. Re-filtering after removal of contaminants (`oeHBCdiff_filtering.R`)
4. Re-normalization after removal of contaminants (`oeHBCdiff_norm.R`)
5. Create SummarizedExperiment object (`oeHBCdiff_makeSE.R`)

### Clustering, Developmental Ordering, Differential Expression
`oeHBCdiff_2_renorm_clust_devO_DE.sh` performs the following analyses, by calling various R scripts (given in parentheses):

6. Clustering using clusterExperiment (`oeHBCdiff_clust.R`)
8. Developmental ordering with slingshot (`oeHBCdiff_slingshot.Rmd`)
9. Differential gene expression using limma, along each lineage (`oeHBCdiff_de.Rmd`)
10. Clustering of differenitally expressed genes along each lineage (`oeHBCdiff_geneClustering.Rmd`)
11. Preparation of gene sets for Gene Set Enrichment Analysis (GSEA; `oeHBCdiff_GSEAprep.Rmd`)
11. GSEA based on cell clustering using limma romer (`oeHBCdiff_romerGSEA.R`)

### Visualizations
`oeHBCdiff_3_viz.sh` performs the following analyses, by calling various R scripts (given in parentheses):

12. Visualizations based on cell clustering (heatmap of marker genes, tSNE plots, PCA pairs plot, cluster & experimental condition bubble plots; `oeHBCdiff_clusterPlots.Rmd`)
12. Visualizations incorporating developmental ordering (3D-PCA plots, dot plots; `oeHBCdiff_devorderplots.Rmd`)
12. Heatmaps of cell cycle genes in the neuronal and sustentacular lineages (`oeHBCdiff_cellCycle.Rmd`)
12. Heatmaps of differentially expressed transcription factors by lineage (`oeHBCdiff_tf_hm.R`)
12. Transcription factor co-expression, network analysis, and visualizations (`oeHBCdiff_tf.Rmd`)
13. Heatmaps of gene clustering in developmental order (`oeHBCdiff_geneClustHeatmaps.Rmd`)
13. Plots of individual or pairs of genes in developmental order (`oeHBCdiff_genePlots.Rmd`)
14. Barplots of GSEA, showing top 100 enriched gene sets per cluster (`oeHBCdiff_GSEAplots.Rmd`)
15. Volcano plots of differentially expressed genes (`oeHBC_volcano.R`)
16. Olfactory Receptor (OR) gene and OR regulation associated gene expression plots (`oeHBCdiff_OR.R`)

### Wnt Pathway Visualization in Figure 6
PathVisio (Version 3.2.4, pathvisio.org) was used to display differential gene expression for Wnt pathway members expressed in the HBCs.  To reproduce this plot, download and install PathVisio and follow the instructions below. 

Download the Mm_Derby_Ensembl_85.bridge gene reference data via the downloads link at pathvisio.org. The Mm_Wnt_Signaling_Pathway_and_Pluripotency_WP723_89312.gpml pathway from the wikipathways_Mus_musculus_Curation-AnalysisCollection (available from the pathvisio.org download link) was used as a starting point. This pathway was modified to only include genes present in our data set after gene filtering. Then select genes were removed and added to focus the pathway on canonical Wnt signaling and to include relevant factors for our experiment. The genes in the pathway were colored by differential expression (log2FC) for the HBCs (cluster 1) relative to all other clusters. 

The input DE data (HBConeVsAllDE_WntSigPath.txt) and the modified pathway (WntPathway.gpml) are in the ref directory of this repository. To produce the diagram, load PathVisio; enter the File menu and open the WntPathway.gpml pathway file; then enter the Data menu and select "Select Gene Database" and load the Mm_Derby_Ensembl_85 data; enter the Data menu and select "Import expression data" and load HBConeVsAllDE_WntSigPath.txt. After the pathway is rendered, enter the Data menu, select "visualization options", select "Text Label", "Expression as color", and "logFC" and then modify the default color palette to match the final output (blue to red) presented in the paper.


### Dependencies/useful R packages:

- SCONE (normalization): https://github.com/YosefLab/scone
- clusterExperiment (clustering): https://github.com/epurdom/clusterExperiment
- slingshot (lineage trajectory algorithm): https://github.com/kstreet13/slingshot
