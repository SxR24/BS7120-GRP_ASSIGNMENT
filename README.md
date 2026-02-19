# BS7120-GRP_ASSIGNMENT
RNA-seq analyss group assignment 
# Figure 1B replication 
Figure 1B reconstructs lineage relationships in the human fetal neocortex (12–13 weeks postconception) by correlating single-cell transcriptomes with:

- Bulk RNA-seq data from cortical zones (VZ, iSVZ, oSVZ, CP)

- Bulk RNA-seq data from purified cell types (aRG, bRG, neurons)

The final output consists of two heatmap panels:

1. Normalized correlation to cortical zones

2. Normalized correlation to purified cell types

Both panels use identical cell ordering and Z-score normalization.
# Data Sources
The following GEO datasets are used:
1. Single-cell master matrix from GSE75140 
Contains both fetal and organoid cells , Only fetal 12–13 wpc cells are used for Figure 1B
2. Bulk cortical zones reference , GEO: GSE38805
File: GSE38805_human_FPKM.txt
Used columns: w13_VZ, w13_ISVZ, w13_OSVZ, w13_CP
3. Purified cell-type reference , GEO: GSE65000
File: GSE65000_hsa_fpkm_matrix.txt

# Reqiurements: 
install.packages(c("pheatmap", "gridExtra")) 
change the paths , they stored as variables 
