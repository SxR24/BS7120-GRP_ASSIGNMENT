library(shiny)

# Plot files

original_plots <- list(
  "Figure S1A" = c("Figure_S1A.rds"),
  "Figure 1B Heatmap A" = c("Figure_1B_heatmapa.rds"),
  "Figure 1B Heatmap B" = c("Figure_1B_heatmapb.rds"),
  "Figure 1C Heatmap" = c("Figure_1C_heatmap_plot.rds"),
  "Figure 1C Heatmap Part 2" = c("Figure_1C_part2_heatmap_plot.rds"),
  "Figure S1B Heatmap" = c("Figure_S1B_heatmap_plot.rds"),
  "Figure S1C Heatmap" = c("Figure_S1C_heatmap_plot.rds"),
  "Original Figure 2A" = c("OldpipelineFig2A.rds"),
  "Original Figure 2B.1" = c("OldpipelineFig2B1.rds"),
  "Original Figure 2B.2" = c("OldpipelineFig2B2.rds"),
  "Original Figure 2B.3" = c("OldpipelineFig2B3.rds"),
  "Original Figure 2B.4" = c("OldpipelineFig2B4.rds"),
  "Figure 3D" = c("Fig3D_plot.rds"),
  "Original Figure 3E" = c("OldpipelineFig3E.rds"),
  "Figure S3H" = c("FigS3H.rds"),
  "Organoid Cluster Comparison" = c("Organoidclustercomparison.rds"),
  "PCA Overlap Plot" = c("PCA_overlap_plot.rds"),
  "Figure 1C Comparison" = c("Fig1Ccomparison.rds"),
  "Figure S1A Comparison" = c("FigS1Acomparison.rds"),
  "Figure S1A Comparison 2" = c("FigS1Acomparison2.rds", "Fig1SAcomparison2.rds")
)

alternative_plots <- list(
  "Alternative Figure 3D" = c("Alt_Fig3D.rds", "Alt_Fig3D - Copy.rds"),
  "Alternative Figure S3H" = c("Alt_FigS3H.rds", "Alt_FigS3H - Copy.rds"),
  "Protein Prefilter" = c("Altpipelineprotienprefilter.rds", "Altpipelineproteinprefilter.rds"),
  "Protein Figure S1A" = c("AltFigure_S1A_protein.rds", "AltFigure_S1A_protein - Copy.rds"),
  "Protein Figure S1B Heatmap" = c(
    "AltFigure_S1B_heatmap_plot_protein.rds",
    "AltFigure_S1B_heatmap_plot_protein - Copy.rds",
    "AltFigure_S1B_heatmap_plot_protien.rds",
    "AltFigure_S1B_heatmap_plot_protien - Copy.rds"
  ),
  "Protein Figure S1C Heatmap" = c(
    "AltFigure_S1c_heatmap_plot_protein.rds",
    "AltFigure_S1c_heatmap_plot_protein - Copy.rds",
    "AltFigure_S1c_heatmap_plot_protien.rds",
    "AltFigure_S1c_heatmap_plot_protien - Copy.rds"
  ),
  "Protein Figure 3E" = c("AltprotienFig3E.rds"),
  "Protein Cell-Type Heatmap" = c(
    "Altpipeline_protien_fetalcelltypeheatmap.rds",
    "Altpipeline_protein_fetalcelltypeheatmap.rds"
  ),
  "Protein Cluster Expression" = c(
    "Altpipeline_protien_fetalclusterexpression.rds",
    "Altpipeline_protein_fetalclusterexpression.rds"
  ),
  "Protein Cluster-Zone Correlation" = c(
    "Altpipeline_protien_fetalclusterzonecorrelation.rds",
    "Altpipeline_protein_fetalclusterzonecorrelation.rds"
  ),
  "Protein Proto-Cluster" = c(
    "Altpipeline_protien_fetalprotocluster.rds",
    "Altpipeline_protein_fetalprotocluster.rds"
  ),
  "Protein Subcluster" = c(
    "Altpipeline_protien_fetalsubcluster.rds",
    "Altpipeline_protein_fetalsubcluster.rds"
  ),
  "Protein Subcluster Expression" = c(
    "Altpipeline_protien_fetalsubclusterexpression.rds",
    "Altpipeline_protein_fetalsubclusterexpression.rds"
  ),
  "Protein Subpopulation Cluster" = c(
    "Altpipeline_protien_fetalsubpoplationcluster.rds",
    "Altpipeline_protein_fetalsubpoplationcluster.rds"
  ),
  "Protein Zone Heatmap" = c(
    "Altpipeline_protien_fetalzoneheatmap.rds",
    "Altpipeline_protein_fetalzoneheatmap.rds"
  ),
  "Non-Coding Prefilter" = c("Altpipelinenoncodingprefilter.rds"),
  "Non-Coding Figure S1A" = c("AltFigure_S1A_noncoding.rds", "AltFigure_S1A_noncoding - Copy.rds"),
  "Non-Coding Figure S1B Heatmap" = c("AltFigure_S1B_heatmap_plot_noncoding.rds", "AltFigure_S1B_heatmap_plot_noncoding - Copy.rds"),
  "Non-Coding Figure S1C Heatmap" = c("AltFigure_S1c_heatmap_plot_noncoding.rds", "AltFigure_S1c_heatmap_plot_noncoding - Copy.rds"),
  "Non-Coding Figure 3D" = c("AltnoncodingFig3D.rds", "AltnoncodingFig3D - Copy.rds"),
  "Non-Coding Figure 3E" = c("AltnoncodingFig3E.rds", "AltnoncodingFig3E - Copy.rds"),
  "Non-Coding Figure S3H" = c("AltnoncodingFigS3H.rds", "AltnoncodingFigS3H - Copy.rds"),
  "Non-Coding Cell-Type Heatmap" = c("Altpipeline_noncoding_fetalcelltypeheatmap.rds"),
  "Non-Coding Cluster Expression" = c("Altpipeline_noncoding_fetalclusterexpression.rds"),
  "Non-Coding Cluster-Zone Correlation" = c("Altpipeline_noncoding_fetalclusterzonecorrelation.rds"),
  "Non-Coding Proto-Cluster" = c("Altpipeline_noncoding_fetalprotocluster.rds"),
  "Non-Coding Subcluster" = c("Altpipeline_noncoding_fetalsubcluster.rds"),
  "Non-Coding Subcluster Expression" = c("Altpipeline_noncoding_fetalsubclusterexpression.rds"),
  "Non-Coding Subpopulation Cluster" = c("Altpipeline_noncoding_fetalsubpoplationcluster.rds"),
  "Non-Coding Zone Heatmap" = c("Altpipeline_noncoding_fetalzoneheatmap.rds")
)

# Annotation

plot_notes <- list(
  "Figure S1A" = "PCA of 225 quality-filtered fetal single-cell transcriptomes (12 wpc chip 1, 12 wpc chip 2, and 13 wpc) reproduced the original NPC–neuron separation along PC1 (2.9% variance explained). Neurons formed a compact cluster at negative PC1 scores, while NPCs dispersed broadly at positive PC1 scores, reflecting the transcriptional heterogeneity inherent to cycling progenitor populations spanning multiple differentiation states. PC2 (1.6%) captured residual within-group variation without resolving additional subpopulations. Cells from all three samples were fully intermingled within both groups, confirming that PC1 separation is driven by cell identity rather than batch or gestational-age effects. ",
  "Figure 1B Heatmap A" = "Z-score normalised Spearman correlations between 220 fetal single-cell transcriptomes and bulk RNA-seq from four laser-microdissected cortical zones reveal two hierarchically separated groups: 123 cells with peak correlation to VZ, iSVZ, or oSVZ, representing progenitor and germinal-zone identities, and 97 cells with maximum correlation to the CP, reflecting post-mitotic neuronal identity. Cells from all three samples are distributed across both groups, confirming that zone assignment reflects cell identity rather than batch or gestational-age origin. ",
  "Figure 1B Heatmap B" = "Z-score normalised Spearman correlations against FACS-purified aRG, bRG, and neuron references cleanly partition the 220 cells into a progenitor group (48 aRG-like, 4 bRG-like) and a dominant neuronal group (168 cells). The strikingly low bRG representation is consistent with the developmental stage sampled , bRG expansion in the outer SVZ becomes prominent only after 13.5 weeks post-conception and anticipates the proportional underrepresentation of basal progenitors observed in subsequent clustering. ",
  "Figure 1C Heatmap" = "Hierarchical clustering of 220 fetal single-cell transcriptomes using 376 PCA-derived genes resolved seven transcriptionally distinct populations (AP1, AP2, BP1, BP2, N1, N2, N3). Apical progenitor clusters are separated by cell-cycle state, basal progenitor clusters by neurogenic commitment, and the three neuronal clusters by progressive maturation, trajectory corroborated by the cortical-zone sidebar, which places VZ-correlated cells at the progenitor end and CP-correlated cells at the neuronal end. Cluster non-randomness was confirmed by permutation testing against 1,000 random partitions (p < 0.001). ",
  "Figure 1C Heatmap Part 2" = " Selected marker genes confirm the identities of the seven fetal cell clusters defined in Figure 1C part 1. Progenitor-associated clusters show higher expression of radial glia and apical progenitor markers such as PAX6, SOX2, VIM, PROM1, and HES1, while intermediate basal progenitor states are supported by EOMES, NEUROD4, INSM1, ASPM, CCNB1, MKI67, and HES6. Neuronal clusters are distinguished by increased expression of differentiation and maturation markers including NEUROD6, BCL11B, MYT1L, TBR1, and MEF2C. This marker pattern supports a developmental ordering from apical progenitors through basal progenitor-like states to progressively maturing neurons. ",
  "Figure S1B Heatmap" = "Hierarchical clustering of fetal single-cell transcriptomes using PC1-associated gene sets resolves two dominant transcriptional programmes: genes enriched for cell cycle and forebrain development mark the progenitor cluster, while neuronal differentiation, cell adhesion, and vesicle-transport genes define the neuronal cluster.",
  "Figure S1C Heatmap" = "Hierarchical clustering using six interneuron marker genes (GAD1, DLX1, DLX2, DLX5, DLX6, ERBB4) identifies five cells with strong coordinated expression of the full marker panel, indicating a ventral telencephalic origin distinct from the dorsal neocortical lineage under study. The remaining cells show no appreciable interneuron marker expression, confirming the specificity of the separation. These five cells were excluded, retaining 220 dorsal neocortical cells for all subsequent clustering and lineage analyses. ",
  "Original Figure 2A" = "Pseudotemporal trajectory of 220 fetal neocortical cells reconstructed using Monocle 3 on ICA-derived coordinates, reproducing the AP→BP→Neuron lineage. The principal graph traces a continuous path from cycling apical progenitors (AP1, AP2) through transitional basal progenitors (BP1, BP2) to progressively maturing neurons (N1, N2, N3), with a visible side branch at the BP junction reflecting rare alternative neurogenic routes. Cell-type ordering along the trajectory is consistent with the hierarchical cluster identities established in Fig. 1C, confirming that the seven populations represent sequential rather than independent transcriptional states. ",
  "Original Figure 2B.1" = "The same pseudotemporal trajectory from figure 2A but coloured by bulk-zone assignment confirms spatial coherence of the lineage: VZ-correlated cells occupy the progenitor end, iSVZ and oSVZ cells populate the intermediate BP region, and CP-correlated cells concentrate at the neuronal end, spatially validating that pseudotime ordering recapitulates the known inside-out organisation of cortical neurogenesis",
  "Original Figure 2B.2" = "SOX2 expression is highest at the progenitor end of the trajectory and progressively extinguished along the path toward neurons, confirming its role as an apical radial glia maintenance factor that is downregulated upon neurogenic commitment. ",
  "Original Figure 2B.3" = "EOMES expression peaks at the AP-to-BP transition zone of the trajectory and is absent from both the earliest progenitors and mature neurons, reproducing its known function as a transient marker of basal progenitor specification and neurogenic commitment. ",
  "Original Figure 2B.4" = "MYT1L expression is absent in progenitors and rises sharply at the neuronal end of the trajectory, marking post-mitotic neuronal identity. Its restricted onset at the BP-to-neuron transition confirms that MYT1L induction coincides with terminal differentiation rather than intermediate progenitor states. ",
  "Figure 3D" = "t-SNE visualisation of 508 organoids partitioned by Seurat into 10 clusters (C1-C10), with cell shape encoding sample origin across nine time points and microdissected regions (33d-65d whole organoids and r1–r4 cortical regions). Clusters show partial temporal stratification earlier time points concentrate in distinct regions of the embedding while later and microdissected samples distribute across multiple clusters reflecting the progressive cellular diversification of the organoid over development. These 10 clusters form the basis for all subsequent cell-type classification and regional identity assignment. ",
  "Original Figure 3E" = "Eight marker genes overlaid onto the organoid t-SNE embedding assign regional and cell-type identities to the 10 clusters. FOXG1 and OTX2 show mutually exclusive distributions, distinguishing dorsal cortical from ventral forebrain identities within the same organoid. ASPM and LIN28A confirm NPC identity, MYT1L and NEUROD6 mark dorsal cortical neurons, RSPO2 localises the hem signalling centre, and DCN is exclusively restricted to the mesenchymal cluster. ",
  "Figure S3H" = "Heatmap of top differentially expressed genes across the 10 organoid clusters, grouped into six functional gene blocks: cycling cells, neuronal, dorsal forebrain, ventral forebrain, RSPO+, and mesenchyme. Each cluster displays a distinct expression signature confined largely to its assigned gene block, with cycling-cell markers (HMGB2, ASPM, TOP2A) enriched in progenitor clusters, neuronal markers (NEUROD6, MYT1L, BCL11A) restricted to neuronal clusters, and ECM genes (COL1A1, COL3A1, POSTN) exclusively marking the mesenchymal cluster. The clean block structure confirms that the 10 Seurat clusters capture biologically coherent cell populations spanning the full regional and lineage diversity of the cerebral organoid. ",
  "Organoid Cluster Comparison" = "Concordance of organoid cluster assignments between alternative pipelines and the original reproduction across six functional categories. Dorsal forebrain NPCs show the highest agreement (98.4-99.2%), followed by RSPO+ cells (95.5%) and mesenchymal cells (90.9-96.2%). Ventral forebrain neurons display the lowest concordance in both models (76.2-77.4%), consistent with their transcriptional similarity to other inhibitory neuronal states and relative rarity, while overall concordance exceeding 84% across all categories confirms that the organoid's regional and lineage architecture is reproducible across computational frameworks.",
  "PCA Overlap Plot" = "Percentage overlap between the PC-associated gene lists from each pipeline (original reproduction and alternative) and Camp et al.'s published Dataset S1 gene lists, computed across all nine PC loading directions. The original reproduction (green) achieves 73-96% overlap, confirming near-identical feature selection to the published study. Both alternative pipelines show substantially lower overlap (protein-coding 11–32%, non-coding 0-32%), reflecting their fundamentally different Seurat VST feature-selection strategy versus the original loading-threshold approach. ",
  "Figure 1C Comparison" = "Percentage overlap of PC gene lists with Camp et al.'s Dataset S1, evaluated specifically on the fetal cortical cell-classifying gene set. The original reproduction (green) maintains high concordance (83-96%) across all components. The protein-coding alternative (blue) achieves notably higher overlap than in the full dataset (73-98%), suggesting that restricting to the cortical subset recovers greater gene-level agreement with modern feature selection. The non-coding model (orange) shows moderate overlap only for PC2 components (42-50%) and zero overlap elsewhere, consistent with its expanded feature space diluting the cortical-specific signal. ",
  "Figure S1A Comparison" = "Per-category concordance of alternative pipeline cell assignments against the original reproduction across four classification levels. Broad categories are highly concordant NPCs (93.9-96.9%), neurons (98.8%), aRG (100%), and neurons by cell type (98.8%) confirming stable classification of abundant populations. Concordance decreases for rare transitional states, most notably BP2 (31.2-50.0%) and BP1 (40.0-62.5%), reflecting the inherent sensitivity of boundary-cell classification to feature selection differences rather than biological disagreement. ",
  "Figure S1A Comparison 2" = "Proportional cell-type distributions across all three pipelines show near-identical NPC/neuron ratios (~28-29% NPCs, ~71-72% neurons) and consistent zone and purified cell-type proportions. The non-coding model uniquely resolves a Cluster 4 (12.0% of subpopulations) capturing a transitional neuronal state distributed across N1 and N2 in the other pipelines, while all remaining subpopulation proportions remain broadly conserved, demonstrating that the overall fetal cortical cell-type architecture is robust across methodological frameworks. ",
  "Alternative Figure 3D" = "Reproduction of Figure 3D using our alternative pipeline. One fewer cluster is shown, though similar groupings can be seen. As with the original figure, shapes indicate the experiments, whilst colours indicate clusters.",
  "Alternative Figure S3H" = "Expression patterns of the identifying genes to distinguish cell type and regional identity. Top differentially expressed genes were identified via a Seurat ROC test, with average expression displayed. The top row highlights functional roles. ",
  "Protein Prefilter" = "Visualisation using the protein‑coding matrix, relating to the original protein-only gene subset, comprising of 20,121 genes. Distributions of library complexity (genes detected per cell), total reads per cell, and mitochondrial transcript fraction were examined. A lower bound of 2000 genes detected was used; cells showing fewer genes were removed. Cells with between 50,000 – 1,200,000 reads were kept. Cells with mitochondrial fraction >10% were removed due to being labelled as stressed or apoptotic.",
  "Protein Figure S1A" = "Visualisation using the protein‑coding matrix, relating to the original protein-only gene subset, comprising of 20,121 genes. Scripts adapted from the original pipeline were used; Principal component analysis was performed on the scaled data of the selected variable genes Cells were projected onto PC1 and PC2. Cells were classified with PC1 < –10 as neurons and those with PC1 > –10 as NPCs.",
  "Protein Figure S1B Heatmap" = "Visualisation using the protein‑coding matrix, relating to the original protein-only gene subset, comprising of 20,121 genes. As with the original Figure S1B, the hierarchical clustering heatmap shows per-cell correlation with PC1, with correlation shown in blue and anti-correlation shown in red. ",
  "Protein Figure S1C Heatmap" = "Visualisation using the protein‑coding matrix, relating to the original protein-only gene subset, comprising of 20,121 genes. Interneurons are identified using the same marker gene panel as in the original Fig. S1C (GAD1, DLX1, DLX2, DLX5, DLX6, ERBB4) ",
  "Protein Figure 3E" = "Visualisation using the protein‑coding matrix, relating to the original protein-only gene subset, comprising of 20,121 genes. Marker genes for each cluster, with expression shown via colour gradient. ",
  "Protein Cell-Type Heatmap" = "Visualisation using the protein‑coding matrix, relating to the original protein-only gene subset, comprising of 20,121 genes. Heatmap presents identification of broad cell types: Apical radial glia (aRG), basal radial glia (bRG) and neurons. The identification is based on correlation (via Spearman coefficient) after Z-score normalisation, with an external reference of FACS-purified cell types.",
  "Protein Cluster Expression" = "Visualisation using the protein‑coding matrix, relating to the original protein-only gene subset, comprising of 20,121 genes. UMAP of neural progenitor cells, with expression of cell‑type marker genes. These are the same list of genes used in the original analysis; colouring is based upon Z-score correlation with marker genes suggested in the original paper.",
  "Protein Cluster-Zone Correlation" = "Visualisation using the protein‑coding matrix, relating to the original protein-only gene subset, comprising of 20,121 genes. UMAP of neural progenitor cells, with expression of cortical zone; colouring is based upon Z-score correlation with bulk zone references: CP, iSVZ, oSVZ and VZ.",
  "Protein Proto-Cluster" = "Visualisation using the protein‑coding matrix, relating to the original protein-only gene subset, comprising of 20,121 genes. UMAP of neural progenitor cells, with 6 clusters. ",
  "Protein Subcluster" = "Visualisation using the protein‑coding matrix, relating to the original protein-only gene subset, comprising of 20,121 genes. Subset of such matrix is shown, showing 62 neural progenitor cells (NPCs); subset of those identified from Figure S1A PCA1 threshold. Nearest neighbour graph with k = 5, followed by clustering via the Leiden algorithm; 5 clusters are shown.",
  "Protein Subcluster Expression" = "Visualisation using the protein‑coding matrix, relating to the original protein-only gene subset, comprising of 20,121 genes. Subset of such matrix is shown, showing 62 neural progenitor cells (NPCs); subset of those identified from Figure S1A PCA1 threshold. UMAP of this subset is shown, with expression of cell‑type marker genes. These are the same list of genes used in the original analysis; colouring is based upon Z-score correlation with marker genes suggested in the original paper. ",
  "Protein Subpopulation Cluster" = "Visualisation using the protein‑coding matrix, relating to the original protein-only gene subset, comprising of 20,121 genes. Overlay of cell-types used in original organoid analysis (Apical progenitors, basal progenitors, neurons) for the identification of biological identities in the resulting clusters, from the UMAP of neural progenitor cells.",
  "Protein Zone Heatmap" = "Visualisation using the protein‑coding matrix, relating to the original protein-only gene subset, comprising of 20,121 genes. Heatmap of resulting clustering from UMAP of neural progenitor cells, where colouring is based upon per-cell Z-score correlation with bulk zone references: CP, iSVZ, oSVZ and VZ. This allows identification of further biological identities to each cluster. ",
  "Non-Coding Prefilter" = "Visualisation using the non-coding matrix, used to investigate whether original protein-coding-only subset remained optimal and only excludes non-coding genes with no annotation, comprising of 78,899 genes. Distributions of library complexity (genes detected per cell), total reads per cell, and mitochondrial transcript fraction were examined. A lower bound of 2000 genes detected was used; cells showing fewer genes were removed. Cells with between 50,000 – 1,200,000 reads were kept. Cells with mitochondrial fraction >10% were removed due to being labelled as stressed or apoptotic.",
  "Non-Coding Figure S1A" = "Visualisation using the non-coding matrix, used to investigate whether original protein-coding-only subset remained optimal and only excludes non-coding genes with no annotation, comprising of 78,899 genes. Scripts adapted from the original pipeline were used; Principal component analysis was performed on the scaled data of the selected variable genes Cells were projected onto PC1 and PC2. Cells were classified with PC1 < –10 as neurons and those with PC1 > –10 as NPCs. ",
  "Non-Coding Figure S1B Heatmap" = "Visualisation using the non-coding matrix, used to investigate whether original protein-coding-only subset remained optimal and only excludes non-coding genes with no annotation, comprising of 78,899 genes. As with the original Figure S1B, the hierarchical clustering heatmap shows per-cell correlation with PC1, with correlation shown in blue and anti-correlation shown in red.",
  "Non-Coding Figure S1C Heatmap" = "Visualisation using the non-coding matrix, used to investigate whether original protein-coding-only subset remained optimal and only excludes non-coding genes with no annotation, comprising of 78,899 genes. Interneurons are identified using the same marker gene panel as in the original Fig. S1C (GAD1, DLX1, DLX2, DLX5, DLX6, ERBB4).",
  "Non-Coding Figure 3D" = "Visualisation using the non-coding matrix, used to investigate whether original protein-coding-only subset remained optimal and only excludes non-coding genes with no annotation, comprising of 78,899 genes. Reproduction of Figure 3D using our alternative pipeline, via expanded gene matrix. Two fewer clusters are shown, though similar groupings can be seen. As with the original figure, shapes indicate the experiments, whilst colours indicate clusters.",
  "Non-Coding Figure 3E" = "Visualisation using the non-coding matrix, used to investigate whether original protein-coding-only subset remained optimal and only excludes non-coding genes with no annotation, comprising of 78,899 genes. Marker genes for each cluster, with expression shown via colour gradient. ",
  "Non-Coding Figure S3H" = "Visualisation using the non-coding matrix, used to investigate whether original protein-coding-only subset remained optimal and only excludes non-coding genes with no annotation, comprising of 78,899 genes. Expression patterns of the identifying genes to distinguish cell type and regional identity. Top differentially expressed genes were identified via a Seurat ROC test, with average expression displayed. The top row highlights functional roles.",
  "Non-Coding Cell-Type Heatmap" = "Visualisation using the non-coding matrix, used to investigate whether original protein-coding-only subset remained optimal and only excludes non-coding genes with no annotation, comprising of 78,899 genes. Heatmap presents identification of broad cell types: Apical radial glia (aRG), basal radial glia (bRG) and neurons. The identification is based on correlation (via Spearman coefficient) after Z-score normalisation, with an external reference of FACS-purified cell types.",
  "Non-Coding Cluster Expression" = "Visualisation using the non-coding matrix, used to investigate whether original protein-coding-only subset remained optimal and only excludes non-coding genes with no annotation, comprising of 78,899 genes. UMAP of neural progenitor cells, with expression of cell‑type marker genes. These are the same list of genes used in the original analysis; colouring is based upon Z-score correlation with marker genes suggested in the original paper.",
  "Non-Coding Cluster-Zone Correlation" = "Visualisation using the non-coding matrix, used to investigate whether original protein-coding-only subset remained optimal and only excludes non-coding genes with no annotation, comprising of 78,899 genes. UMAP of neural progenitor cells, with expression of cortical zone; colouring is based upon Z-score correlation with bulk zone references: CP, iSVZ, oSVZ and VZ.",
  "Non-Coding Proto-Cluster" = "Visualisation using the non-coding matrix, used to investigate whether original protein-coding-only subset remained optimal and only excludes non-coding genes with no annotation, comprising of 78,899 genes. UMAP of neural progenitor cells, with 6 clusters. Cluster 4 is shown to be isolated from the other clusters in the UMAP; it was therefore incorporated via manual comparison with top 30 differentially expressed genes. ",
  "Non-Coding Subcluster" = "Visualisation using the non-coding matrix, used to investigate whether original protein-coding-only subset remained optimal and only excludes non-coding genes with no annotation, comprising of 78,899 genes. Subset of such matrix is shown, showing 61 neural progenitor cells (NPCs); subset of those identified from Figure S1A PCA1 threshold. Nearest neighbour graph with k = 5, followed by clustering via the Leiden algorithm; 6 clusters are shown.",
  "Non-Coding Subcluster Expression" = "isualisation using the non-coding matrix, used to investigate whether original protein-coding-only subset remained optimal and only excludes non-coding genes with no annotation, comprising of 78,899 genes. Subset of such matrix is shown, showing 61 neural progenitor cells (NPCs); subset of those identified from Figure S1A PCA1 threshold. UMAP of this subset is shown, with expression of cell‑type marker genes. These are the same list of genes used in the original analysis; colouring is based upon Z-score correlation with marker genes suggested in the original paper. ",
  "Non-Coding Subpopulation Cluster" = "Visualisation using the non-coding matrix, used to investigate whether original protein-coding-only subset remained optimal and only excludes non-coding genes with no annotation, comprising of 78,899 genes. Overlay of cell-types used in original organoid analysis (Apical progenitors, basal progenitors, neurons) for the identification of biological identities in the resulting clusters, from the UMAP of neural progenitor cells. Note that group 4 appears to be an outlier and thus was integrated by manually comparing the top 30 differentially expressed genes with known neuronal maturation markers.",
  "Non-Coding Zone Heatmap" = "Visualisation using the non-coding matrix, used to investigate whether original protein-coding-only subset remained optimal and only excludes non-coding genes with no annotation, comprising of 78,899 genes. Heatmap of resulting clustering from UMAP of neural progenitor cells, where colouring is based upon per-cell Z-score correlation with bulk zone references: CP, iSVZ, oSVZ and VZ. This allows identification of further biological identities to each cluster. "
)

# Axis Inoformation



# Conclusion

conclusion_text <- "This website was built to bring the main parts of our project into one place. Instead of leaving the work spread across separate scripts, figures and intermediate files, it provides a clearer way of moving through and visualising the original pipeline and the alternative pipeline, as well as providing a level of interaction and the ability to make queries on specifc genes. This is important for a dataset like this one, where the analysis quickly becomes difficult to follow if everything is only presented as disconnected outputs. An important feature of the site is that it does not just display final figures; the original and alternative tabs let the user move through multiple plots in an organised way, while the interactive panel makes it possible to adjust the size of selected plots and inspect the cell assignment table alongside them. The gene query section provides another layer of analysis by allowing individual genes to be selected directly from the expression matrix, with a summary table, distribution plot and expression preview shown dynamically.Together, these parts make the website more than just a gallery of results. The website, therefore, acts as an accessible front-end for observing both the dataset and its analysis. At the same time, the website also reflects the limits of the work behind it. Not every part of the original study could be reproduced perfectly, and some sections are necessarily based on precomputed objects rather than being generated from scratch inside the app each time. A few outputs were also easier to anchor to the paper than others, especially where the published supplementary data were clearer. Considering these limits, the site still achieves its purpose: it makes the analyses easier to navigate, allows the user to compare the main strands of the project in one place, and gives a more accessible view of the Camp et al. (2015) dataset than a directory of separate plots alone."

# About
about_text <- "This website has been developed to facilitate the visual interpretation of the datasets originally produced and provided by Camp et al. (2015). This includes a reproduction of their analysis, including general scRNA-seq analysis of both cerebral organoid and foetal neocortex, as well as developmental single-cell trajectory analysis, using tools that the authors used. We also include an alternative pipeline, showing that the results of Camp et al. (2015) can be reproduced using varying methods. This website also contains tools for querying specific genes interactive plots. Finally, concluding remarks regarding our analysis can be found; each section is contained within a separate tab."
authors_text <- "Ali, Isa E\nKataria, Nitesh Ghansham\nShaban, Noor H.S.\nWroclawska, Weronika M.\nRamesh Kumar, Sohil Ananth\nHazzard, Jed "
reference_text <- "Camp, J.G. et al. (2015) 'Human cerebral organoids recapitulate gene expression programs of fetal neocortex development', Proceedings of the National Academy of Sciences of the United States of America, 112(51), pp. 15672–15677. Available at: https://doi.org/10.1073/pnas.1520760112 "

# Data files

cell_file <- if (file.exists("cell_assignments.txt")) "cell_assignments.txt" else "cell_assignments"
cell_data <- read.table(cell_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

gene_matrix_file <- if (file.exists("normalized_expression.rds")) {
  "normalized_expression.rds"
} else if (file.exists("/mnt/data/normalized_expression.rds")) {
  "/mnt/data/normalized_expression.rds"
} else {
  stop("Gene matrix file was not found.")
}

# User Interface UI

ui_base <- fluidPage(
  tags$head(
    tags$style(HTML("
      html, body { height: 100%; margin: 0; background: #f5f7fb; font-family: 'Segoe UI', Arial, sans-serif; color: #172033; }
      .container-fluid { padding-left: 0; padding-right: 0; }
      .topbar { min-height: 68px; background: linear-gradient(90deg, #02172f 0%, #041a3c 100%); display: flex; align-items: center; justify-content: space-between; padding: 0 22px; box-shadow: 0 2px 8px rgba(0,0,0,0.12); }
      .topbar-left { display: flex; align-items: center; gap: 20px; flex-wrap: wrap; }
      .brand { color: white; font-size: 18px; font-weight: 700; display: flex; align-items: center; gap: 12px; }
      .brand-icon { width: 26px; height: 26px; border-radius: 50%; border: 2px solid #1e90ff; display: inline-flex; align-items: center; justify-content: center; color: #1e90ff; font-weight: 700; }
      .top-nav { display: flex; align-items: center; gap: 10px; flex-wrap: wrap; }
      .nav-btn { background: transparent; color: #e7edf8; border: none; border-radius: 12px; padding: 12px 18px; font-size: 15px; font-weight: 600; cursor: pointer; }
      .nav-btn.active-tab { background: #153b79; color: white; }

      .app-shell { display: flex; height: calc(100vh - 68px); }
      .sidebar { width: 340px; min-width: 340px; background: #f8fafc; border-right: 1px solid #e2e8f0; padding: 18px 16px; box-sizing: border-box; }
      .sidebar-card { background: #ffffff; border-radius: 16px; box-shadow: 0 3px 14px rgba(15, 23, 42, 0.07); border: 1px solid #e5eaf2; padding: 16px; height: calc(100vh - 104px); display: flex; flex-direction: column; }
      .side-section-title { font-size: 12px; font-weight: 800; letter-spacing: 0.6px; color: #25406e; margin-bottom: 12px; text-transform: uppercase; }
      .plot-select-wrap { margin-top: 8px; }
      .plot-select-wrap .form-control { border-radius: 12px; border: 1px solid #d8e0ee; height: 46px; box-shadow: none; }
      .plot-list-wrap { margin-top: 14px; flex: 1 1 auto; min-height: 0; overflow-y: auto; border-top: 1px solid #eef2f7; padding-top: 10px; }
      .plot-list-button { width: 100%; text-align: left; border: none; background: transparent; border-radius: 14px; padding: 14px 14px; margin-bottom: 8px; color: #24344d; font-size: 15px; font-weight: 600; cursor: pointer; }
      .plot-list-button:hover { background: #eef4ff; }
      .plot-list-button.active-plot { background: #dbe6fb; color: #1c49b8; font-weight: 700; }

      .main-content { flex: 1; padding: 20px 22px; overflow: auto; }
      .content-header { display: flex; align-items: flex-start; justify-content: space-between; gap: 20px; margin-bottom: 16px; }
      .content-title { font-size: 28px; font-weight: 800; color: #172033; margin: 0 0 6px 0; }
      .breadcrumb-line { font-size: 15px; color: #5a6c8f; }
      .breadcrumb-line .crumb-accent { color: #2352d7; font-weight: 600; }

      .plot-panel { background: white; border-radius: 16px; border: 1px solid #e5eaf2; box-shadow: 0 3px 14px rgba(15, 23, 42, 0.07); padding: 8px; margin-bottom: 16px; }
      .annotation-panel { background: white; border-radius: 16px; border: 1px solid #e5eaf2; box-shadow: 0 3px 14px rgba(15, 23, 42, 0.07); padding: 22px; position: relative; margin-bottom: 16px; }
      .annotation-panel:before { content: ''; position: absolute; left: 0; top: 28px; bottom: 28px; width: 4px; background: #2563eb; border-radius: 3px; }
      .annotation-title { font-size: 16px; font-weight: 800; color: #172033; margin-bottom: 14px; padding-left: 18px; }
      .annotation-text { color: #4a5a75; font-size: 15px; line-height: 1.8; padding-left: 18px; white-space: pre-wrap; }

      .interactive-shell { padding: 20px 22px; }
      .interactive-card { background: white; border-radius: 16px; border: 1px solid #e5eaf2; box-shadow: 0 3px 14px rgba(15, 23, 42, 0.07); padding: 18px; margin-bottom: 16px; }
      .control-label { font-weight: 700; color: #334155; }
      .form-control { border-radius: 12px; border: 1px solid #d8e0ee; box-shadow: none; }
      .table { background: white; }

      .conclusion-box { background: white; border-radius: 16px; border: 1px solid #e5eaf2; box-shadow: 0 3px 14px rgba(15, 23, 42, 0.07); padding: 30px; min-height: 300px; }

      .axis-card, .gene-card { background: white; border-radius: 16px; border: 1px solid #e5eaf2; box-shadow: 0 3px 14px rgba(15, 23, 42, 0.07); padding: 18px; margin-top: 16px; }

      .axis-table-wrap { width: 100%; overflow-x: auto; }
      .axis-table { width: 100%; border-collapse: separate; border-spacing: 0; font-size: 15px; color: #24344d; background: #ffffff; border: 1px solid #e2e8f0; border-radius: 12px; overflow: hidden; }
      .axis-table thead th { background: #eff4fb; color: #10233f; font-weight: 700; text-align: left; padding: 14px 16px; border-bottom: 1px solid #dbe5f0; }
      .axis-table tbody td { padding: 14px 16px; border-bottom: 1px solid #edf2f7; vertical-align: top; }
      .axis-table tbody tr:last-child td { border-bottom: none; }
      .axis-table tbody tr:nth-child(even) { background: #fafcff; }
      .axis-col { width: 28%; font-weight: 700; color: #1c49b8; white-space: nowrap; }
      .desc-col { width: 72%; color: #334155; }

      .gene-layout { padding: 20px 22px; }
      .gene-page-grid { display: grid; grid-template-columns: 340px 1fr; gap: 18px; align-items: start; }
      .gene-control-card { background: white; border-radius: 16px; border: 1px solid #e5eaf2; box-shadow: 0 3px 14px rgba(15, 23, 42, 0.07); padding: 18px; }
      .gene-result-card { background: white; border-radius: 16px; border: 1px solid #e5eaf2; box-shadow: 0 3px 14px rgba(15, 23, 42, 0.07); padding: 18px; margin-bottom: 16px; }
      .small-muted { color: #64748b; font-size: 14px; line-height: 1.6; }
    "))
  ),
  
  div(
    class = "topbar",
    div(
      class = "topbar-left",
      div(class = "brand", "Group-A Dashboard: Analysis of Human Cerebral Organoids"),
      div(
        class = "top-nav",
        actionButton("nav_original", "Original Pipeline", class = "nav-btn"),
        actionButton("nav_alternative", "Alternative Pipeline", class = "nav-btn"),
        actionButton("nav_interactive", "Interactive Panel", class = "nav-btn"),
        actionButton("nav_gene", "Gene Query", class = "nav-btn"),
        actionButton("nav_conclusion", "Conclusion", class = "nav-btn"),
        actionButton("nav_about", "About", class = "nav-btn")
      )
    ),
    div()
  ),
  
  uiOutput("page_ui")
)

# Server 

server <- function(input, output, session) {
  
  plot_cache <- reactiveVal(new.env(parent = emptyenv()))
  gene_data_cache <- reactiveVal(NULL)
  current_page <- reactiveVal("original")
  
  observeEvent(input$nav_original, current_page("original"))
  observeEvent(input$nav_alternative, current_page("alternative"))
  observeEvent(input$nav_interactive, current_page("interactive"))
  observeEvent(input$nav_gene, current_page("gene"))
  observeEvent(input$nav_conclusion, current_page("conclusion"))
  observeEvent(input$nav_about, current_page("about"))
  
  output$orig_plot_list_ui <- renderUI({
    req(input$orig_plot_choice)
    tagList(lapply(names(original_plots), function(p) {
      btn_id <- paste0("orig_plot_btn_", gsub("[^A-Za-z0-9]", "_", p))
      cls <- if (identical(p, input$orig_plot_choice)) "plot-list-button active-plot" else "plot-list-button"
      actionButton(btn_id, label = p, class = cls)
    }))
  })
  
  output$alt_plot_list_ui <- renderUI({
    req(input$alt_plot_choice)
    tagList(lapply(names(alternative_plots), function(p) {
      btn_id <- paste0("alt_plot_btn_", gsub("[^A-Za-z0-9]", "_", p))
      cls <- if (identical(p, input$alt_plot_choice)) "plot-list-button active-plot" else "plot-list-button"
      actionButton(btn_id, label = p, class = cls)
    }))
  })
  
  observe({
    lapply(names(original_plots), function(p) {
      local({
        plot_name <- p
        btn_id <- paste0("orig_plot_btn_", gsub("[^A-Za-z0-9]", "_", plot_name))
        observeEvent(input[[btn_id]], {
          updateSelectInput(session, "orig_plot_choice", selected = plot_name)
        }, ignoreInit = TRUE)
      })
    })
  })
  
  observe({
    lapply(names(alternative_plots), function(p) {
      local({
        plot_name <- p
        btn_id <- paste0("alt_plot_btn_", gsub("[^A-Za-z0-9]", "_", plot_name))
        observeEvent(input[[btn_id]], {
          updateSelectInput(session, "alt_plot_choice", selected = plot_name)
        }, ignoreInit = TRUE)
      })
    })
  })
  
  output$page_ui <- renderUI({
    if (current_page() == "original") {
      div(
        class = "app-shell",
        div(
          class = "sidebar",
          div(
            class = "sidebar-card",
            div(class = "side-section-title", "Original Pipeline"),
            div(
              class = "plot-select-wrap",
              selectInput(
                "orig_plot_choice",
                NULL,
                choices = names(original_plots),
                selected = if (is.null(input$orig_plot_choice)) names(original_plots)[1] else input$orig_plot_choice,
                width = "100%"
              )
            ),
            div(class = "plot-list-wrap", uiOutput("orig_plot_list_ui"))
          )
        ),
        div(
          class = "main-content",
          div(
            class = "content-header",
            div(
              div(class = "content-title", textOutput("orig_selected_plot_title")),
              div(
                span(textOutput("orig_selected_plot_breadcrumb"))
              )
            )
          ),
          div(class = "plot-panel", plotOutput("orig_plot", height = "760px")),
          div(
            class = "annotation-panel",
            div(class = "annotation-title", "Annotation"),
            div(class = "annotation-text", textOutput("orig_annotation"))
          )
        )
      )
    } else if (current_page() == "alternative") {
      div(
        class = "app-shell",
        div(
          class = "sidebar",
          div(
            class = "sidebar-card",
            div(class = "side-section-title", "Alternative Pipeline"),
            div(
              class = "plot-select-wrap",
              selectInput(
                "alt_plot_choice",
                NULL,
                choices = names(alternative_plots),
                selected = if (is.null(input$alt_plot_choice)) names(alternative_plots)[1] else input$alt_plot_choice,
                width = "100%"
              )
            ),
            div(class = "plot-list-wrap", uiOutput("alt_plot_list_ui"))
          )
        ),
        div(
          class = "main-content",
          div(
            class = "content-header",
            div(
              div(class = "content-title", textOutput("alt_selected_plot_title")),
              div(
                
                span(textOutput("alt_selected_plot_breadcrumb"))
              )
            )
          ),
          div(class = "plot-panel", plotOutput("alt_plot", height = "760px")),
          div(
            class = "annotation-panel",
            div(class = "annotation-title", "Annotation or Caption"),
            div(class = "annotation-text", textOutput("alt_annotation"))
          )
        )
      )
    } else if (current_page() == "interactive") {
      div(
        class = "interactive-shell",
        div(
          class = "interactive-card",
          fluidRow(
            column(
              3,
              selectInput(
                "interactive_plot_choice",
                "Choose Plot",
                choices = c(
                  "PCA Overlap" = "pca",
                  "Organoid Cluster Comparison" = "org",
                  "Original Figure 3D" = "orig3d",
                  "Alternative Figure 3D" = "alt3d",
                  "Non Coding Figure 3D" = "nc3d"
                )
              ),
              sliderInput("interactive_plot_width", "Plot Width Percent", min = 50, max = 100, value = 100, step = 5),
              sliderInput("interactive_plot_height", "Plot Height Pixels", min = 500, max = 1200, value = 700, step = 50),
              numericInput("n_rows", "Rows of cell assignment table", value = 20, min = 5, max = nrow(cell_data), step = 5),
              checkboxInput("show_table", "Show cell assignment table", TRUE)
            ),
            column(
              9,
              uiOutput("interactive_plot_ui"),
              
            )
          )
        ),
        conditionalPanel(
          condition = "input.show_table == true",
          div(class = "interactive-card", h3("Cell Assignment Table"), tableOutput("cell_table"))
        )
      )
    } else if (current_page() == "gene") {
      div(
        class = "gene-layout",
        div(
          class = "content-header",
          div(
            div(class = "content-title", "Gene Query")
            
          )
        ),
        div(
          class = "gene-page-grid",
          div(
            class = "gene-control-card",
            h4("Select Gene"),
            p(class = "small-muted", "Choose a gene from the expression matrix dropdown."),
            uiOutput("gene_dropdown_ui"),
            br(),
          
          ),
          div(
            div(
              class = "gene-result-card",
              h3(textOutput("gene_result_title")),
              tableOutput("gene_summary_table")
            ),
            div(
              class = "gene-result-card",
              h3("Expression Distribution"),
              plotOutput("gene_expression_plot", height = "520px")
            ),
            div(
              class = "gene-result-card",
              h3("Preview of Expression Values"),
              tableOutput("gene_preview_table")
            )
          )
        )
      )
    } else if (current_page() == "about") {
      div(
        class = "interactive-shell",
        div(
          class = "conclusion-box",
          h2("About"),
          div(style = "white-space: pre-wrap; line-height: 1.8; color: #4a5a75;", about_text),
          h2("Authors"),
          div(style = "white-space: pre-wrap; line-height: 1.8; color: #4a5a75;", authors_text),
          h2("Reference"),
          div(style = "white-space: pre-wrap; line-height: 1.8; color: #4a5a75;", reference_text)
        )
      )
    } else {
      div(
        class = "interactive-shell",
        div(
          class = "conclusion-box",
          h2("Conclusion"),
          div(style = "white-space: pre-wrap; line-height: 1.8; color: #4a5a75;", conclusion_text)
        )
      )
    }
  })
  
  output$orig_selected_plot_title <- renderText({
    req(input$orig_plot_choice)
    input$orig_plot_choice
  })
  
  output$orig_selected_plot_breadcrumb <- renderText({
    req(input$orig_plot_choice)
    input$orig_plot_choice
  })
  
  output$alt_selected_plot_title <- renderText({
    req(input$alt_plot_choice)
    input$alt_plot_choice
  })
  
  output$alt_selected_plot_breadcrumb <- renderText({
    req(input$alt_plot_choice)
    input$alt_plot_choice
  })
  
  output$orig_plot <- renderPlot({
    req(input$orig_plot_choice)
    
    files <- original_plots[[input$orig_plot_choice]]
    env <- plot_cache()
    key <- paste0("original::", input$orig_plot_choice)
    
    if (!exists(key, envir = env, inherits = FALSE)) {
      found <- files[file.exists(files)]
      if (length(found) == 0) stop(paste("Missing plot file for", input$orig_plot_choice))
      assign(key, readRDS(found[1]), envir = env)
      plot_cache(env)
    }
    
    obj <- get(key, envir = env, inherits = FALSE)
    
    if ("recordedplot" %in% class(obj)) {
      replayPlot(obj)
    } else {
      print(obj)
    }
  }, res = 90, width = function() session$clientData$output_orig_plot_width, height = 760)
  
  output$alt_plot <- renderPlot({
    req(input$alt_plot_choice)
    
    files <- alternative_plots[[input$alt_plot_choice]]
    env <- plot_cache()
    key <- paste0("alternative::", input$alt_plot_choice)
    
    if (!exists(key, envir = env, inherits = FALSE)) {
      found <- files[file.exists(files)]
      if (length(found) == 0) stop(paste("Missing plot file for", input$alt_plot_choice))
      assign(key, readRDS(found[1]), envir = env)
      plot_cache(env)
    }
    
    obj <- get(key, envir = env, inherits = FALSE)
    
    if ("recordedplot" %in% class(obj)) {
      replayPlot(obj)
    } else {
      print(obj)
    }
  }, res = 90, width = function() session$clientData$output_alt_plot_width, height = 760)
  
  output$orig_annotation <- renderText({
    req(input$orig_plot_choice)
    if (is.null(plot_notes[[input$orig_plot_choice]])) "" else plot_notes[[input$orig_plot_choice]]
  })
  
  output$alt_annotation <- renderText({
    req(input$alt_plot_choice)
    if (is.null(plot_notes[[input$alt_plot_choice]])) "" else plot_notes[[input$alt_plot_choice]]
  })
  
  output$interactive_plot_ui <- renderUI({
    req(input$interactive_plot_width, input$interactive_plot_height)
    div(
      style = paste0("width:", input$interactive_plot_width, "%; margin: 0 auto;"),
      div(
        class = "plot-panel",
        plotOutput("interactive_plot", height = paste0(input$interactive_plot_height, "px"))
      )
    )
  })
  
  output$interactive_plot <- renderPlot({
    req(input$interactive_plot_choice)
    
    files <- switch(
      input$interactive_plot_choice,
      "pca" = c("PCA_overlap_plot.rds"),
      "org" = c("Organoidclustercomparison.rds"),
      "orig3d" = c("Fig3D_plot.rds"),
      "alt3d" = c("Alt_Fig3D.rds", "Alt_Fig3D - Copy.rds"),
      "nc3d" = c("AltnoncodingFig3D.rds", "AltnoncodingFig3D - Copy.rds")
    )
    
    key <- paste0("interactive::", input$interactive_plot_choice)
    env <- plot_cache()
    
    if (!exists(key, envir = env, inherits = FALSE)) {
      found <- files[file.exists(files)]
      if (length(found) == 0) stop("Missing interactive plot file")
      assign(key, readRDS(found[1]), envir = env)
      plot_cache(env)
    }
    
    obj <- get(key, envir = env, inherits = FALSE)
    
    if ("recordedplot" %in% class(obj)) {
      replayPlot(obj)
    } else {
      print(obj)
    }
  }, res = 110, width = function() session$clientData$output_interactive_plot_width, height = function() input$interactive_plot_height)

  
  output$cell_table <- renderTable({
    req(input$n_rows)
    head(cell_data, input$n_rows)
  }, striped = TRUE, bordered = TRUE, spacing = "s")
  
  output$gene_matrix_status <- renderText({
  })
  
  output$gene_dropdown_ui <- renderUI({
    obj <- gene_data_cache()
    
    if (is.null(obj)) {
      if (!file.exists(gene_matrix_file)) {
        return(paste("Gene matrix file not found:", gene_matrix_file))
      }
      
      obj <- readRDS(gene_matrix_file)
      gene_data_cache(obj)
    }
    
    if (is.matrix(obj) || inherits(obj, "Matrix")) {
      if (!is.null(rownames(obj))) {
        genes <- rownames(obj)
      } else {
        genes <- colnames(obj)
      }
    } else if (is.data.frame(obj)) {
      genes <- colnames(obj)
    } else {
      genes <- character(0)
    }
    
    genes <- sort(unique(genes))
    genes <- genes[!is.na(genes)]
    genes <- genes[genes != ""]
    
    selectInput(
      "selected_gene",
      "Select Gene",
      choices = genes,
      selected = if ("SOX2" %in% genes) "SOX2" else genes[1],
      width = "100%"
    )
  })
  
  selected_gene_result <- reactive({
    req(input$selected_gene)
    
    gene_name <- input$selected_gene
    obj <- gene_data_cache()
    
    if (is.null(obj)) {
      req(file.exists(gene_matrix_file))
      obj <- readRDS(gene_matrix_file)
      gene_data_cache(obj)
    }
    
    expression_values <- NULL
    
    if (is.matrix(obj) || inherits(obj, "Matrix")) {
      if (!is.null(rownames(obj)) && gene_name %in% rownames(obj)) {
        expression_values <- as.numeric(obj[gene_name, ])
      } else if (!is.null(colnames(obj)) && gene_name %in% colnames(obj)) {
        expression_values <- as.numeric(obj[, gene_name])
      }
    } else if (is.data.frame(obj)) {
      if (gene_name %in% colnames(obj)) {
        expression_values <- as.numeric(obj[[gene_name]])
      } else if (!is.null(rownames(obj)) && gene_name %in% rownames(obj)) {
        expression_values <- as.numeric(obj[gene_name, ])
      }
    }
    
    if (is.null(expression_values)) {
      list(found = FALSE, gene = gene_name)
    } else {
      expression_values <- expression_values[!is.na(expression_values)]
      list(found = TRUE, gene = gene_name, values = expression_values)
    }
  })
  
  output$gene_result_title <- renderText({
    res <- selected_gene_result()
    
    if (!isTRUE(res$found)) {
      paste("Gene not found:", res$gene)
    } else {
      paste("Gene:", res$gene)
    }
  })
  
  output$gene_summary_table <- renderTable({
    res <- selected_gene_result()
    
    if (!isTRUE(res$found)) {
      return(data.frame(Message = paste("The gene", res$gene, "was not found in the expression matrix.")))
    }
    
    values <- res$values
    
    data.frame(
      Metric = c("Gene", "Number of cells or samples", "Detected values above zero", "Mean expression", "Median expression", "Maximum expression"),
      Value = c(
        res$gene,
        length(values),
        sum(values > 0),
        round(mean(values), 4),
        round(median(values), 4),
        round(max(values), 4)
      ),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")
  
  output$gene_expression_plot <- renderPlot({
    res <- selected_gene_result()
    
    if (!isTRUE(res$found)) {
      plot.new()
      text(0.5, 0.5, "No valid gene selected.")
      return()
    }
    
    values <- res$values
    
    par(mar = c(5, 5, 4, 2))
    hist(
      values,
      breaks = 40,
      main = paste("Expression Distribution for", res$gene),
      xlab = "Expression value",
      ylab = "Number of cells or samples"
    )
  }, res = 110)
  
  output$gene_preview_table <- renderTable({
    res <- selected_gene_result()
    
    if (!isTRUE(res$found)) {
      return(data.frame(Message = "No expression values to show."))
    }
    
    values <- res$values
    n_preview <- min(20, length(values))
    
    data.frame(
      Index = seq_len(n_preview),
      Expression = round(values[seq_len(n_preview)], 4)
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")
  

}

ui <- tagList(ui_base)

shinyApp(ui, server)
