# Extended Data Figure 2: LPM Generation and characterization.
library(Seurat)
library(harmony)
library(Nebulosa)
library(patchwork)
library(fgsea)
library(SeuratDisk)
library(ggrepel)
library(GSVA)
library(ggpubr)
library(pbapply)
library(dplyr)


### load LPM
i <- 100
outFolder <- paste0("results/rds/P307earlyLPM_D", i)
setwd(outFolder)
load("P307earlyLPM.RData")

# Extended Figure 2e. UMAP visualization of single-cell RNA sequencing of LPM at day 6, showing clusters of possible lineages emerged from LPM.
library(RColorBrewer)
col <- c(brewer.pal(4, "Greens")[3:1], brewer.pal(9, "Blues")[5], brewer.pal(4, "Greens")[4], col2["Endothelia"])
names(col) <- 0:5

png(paste0("LPM_D100_noaxes_nolabel.png"), height = 1600, width = 1800, res = 300)
UMAPPlot(X, label = FALSE) +
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 0), axis.text.y = element_text(size = 0),
    axis.ticks.x = element_blank(), axis.ticks.y = element_blank()
  ) +
  xlab("") +
  ylab("") +
  scale_color_manual(values = col)
dev.off()
library(svglite)
svglite(paste0("LPM_D100_noaxes_nolabel.svg"), width = 6.75, height = 6)
UMAPPlot(X, label = FALSE) +
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 0), axis.text.y = element_text(size = 0),
    axis.ticks.x = element_blank(), axis.ticks.y = element_blank()
  ) +
  xlab("") +
  ylab("") +
  scale_color_manual(values = col)
dev.off()


# Extended Figure 2f. Heatmap of selected markers indicating different lineages. TBX1, a marker of cardiopharyngeal mesoderm, was not expressed (gray).
dge <- X

library(reshape2)
geneslist <- list(
  "LPM" = c("PRRX1", "TWIST1", "HAND1", "HAND2"),
  "Lineage" = c("HOXA10", "TBX3", "TBX4", "TBX5", "TBX6", "FOXF1", "TCF21", "TBX1"),
  "Hemangioblast" = c("TAL1", "ETV2", "LMO2", "KDR", "FLI1", "RUNX1", "SOX17")
)
genes <- unlist(geneslist)
print(genes[!(genes %in% rownames(dge))])
genes <- genes[which(genes %in% rownames(dge))]
genesgroup <- melt(geneslist)


# Heatmap by clusters
library(ComplexHeatmap)

### by cell types
centroid <- log(AverageExpression(X)$RNA + 1)
centroid.std <- (centroid - apply(centroid, 1, mean)) / apply(centroid, 1, sd)
centroid.std[is.nan(centroid.std)] <- NA

centroid["TBX1", ]
centroid.std["TBX1", ]


genes[!(genes %in% rownames(centroid))]
genes[!(genes %in% rownames(centroid.std))]
# Lineage  "TBX1"
length(which(genes %in% rownames(centroid.std))) # 17
genes2 <- genes[(genes %in% rownames(centroid.std))]

cc <- centroid[which(rownames(centroid) %in% genes), ]
cc.std <- centroid.std[which(rownames(centroid.std) %in% genes), ]
cc.std <- cc.std[genes2, ]
cc.std


png(paste0("LPM_Heatmap_celltype2_noclusterRow.png"), width = 1200, height = 1600, res = 300)
Heatmap(cc.std,
  name = "StdExp",
  row_split = genesgroup[, 2],
  cluster_rows = F,
  row_title_rot = 0,
  row_names_gp = gpar(fontsize = 8)
)
dev.off()
png(paste0("LPM_Heatmap_celltype2_nocluster.png"), width = 1200, height = 1400, res = 300)
Heatmap(cc.std,
  name = "StdExp",
  row_split = genesgroup[, 2],
  cluster_rows = F,
  cluster_columns = F,
  row_title_rot = 0,
  row_names_gp = gpar(fontsize = 8)
)
dev.off()

library(svglite)
svglite(paste0("LPM_Heatmap_celltype2_noclusterRow.svg"), width = 4.5, height = 6)
Heatmap(cc.std,
  name = "StdExp",
  row_split = genesgroup[, 2],
  cluster_rows = F,
  row_title_rot = 0,
  row_names_gp = gpar(fontsize = 8)
)
dev.off()
svglite(paste0("LPM_Heatmap_celltype2_nocluster.svg"), width = 4.5, height = 5.25)
Heatmap(cc.std,
  name = "StdExp",
  row_split = genesgroup[, 2],
  cluster_rows = F,
  cluster_columns = F,
  row_title_rot = 0,
  row_names_gp = gpar(fontsize = 8)
)
dev.off()



sessionInfo()
R version 4.2.3 (2023-03-15)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Rocky Linux 8.10 (Green Obsidian)

Matrix products: default

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] dplyr_1.1.1           pbapply_1.7-0         ggpubr_0.6.0         
 [4] GSVA_1.46.0           ggrepel_0.9.3         SeuratDisk_0.0.0.9020
 [7] fgsea_1.24.0          Nebulosa_1.8.0        patchwork_1.1.2      
[10] ggplot2_3.4.4         harmony_0.1.1         Rcpp_1.0.10          
[13] SeuratObject_4.1.3    Seurat_4.3.0         

loaded via a namespace (and not attached):
  [1] backports_1.4.1             fastmatch_1.1-3            
  [3] plyr_1.8.8                  igraph_1.4.1               
  [5] lazyeval_0.2.2              sp_1.6-0                   
  [7] GSEABase_1.60.0             splines_4.2.3              
  [9] BiocParallel_1.32.6         listenv_0.9.0              
 [11] scattermore_0.8             GenomeInfoDb_1.34.9        
 [13] digest_0.6.31               htmltools_0.5.5            
 [15] fansi_1.0.4                 memoise_2.0.1              
 [17] magrittr_2.0.3              ScaledMatrix_1.6.0         
 [19] tensor_1.5                  cluster_2.1.4              
 [21] ks_1.14.0                   ROCR_1.0-11                
 [23] Biostrings_2.66.0           globals_0.16.2             
 [25] annotate_1.76.0             matrixStats_0.63.0         
 [27] spatstat.sparse_3.0-1       colorspace_2.1-0           
 [29] blob_1.2.4                  crayon_1.5.2               
 [31] RCurl_1.98-1.12             jsonlite_1.8.4             
 [33] graph_1.76.0                progressr_0.13.0           
 [35] spatstat.data_3.0-1         survival_3.5-5             
 [37] zoo_1.8-11                  glue_1.6.2                 
 [39] polyclip_1.10-4             gtable_0.3.3               
 [41] zlibbioc_1.44.0             XVector_0.38.0             
 [43] leiden_0.4.3                DelayedArray_0.24.0        
 [45] car_3.1-2                   BiocSingular_1.14.0        
 [47] Rhdf5lib_1.20.0             future.apply_1.10.0        
 [49] SingleCellExperiment_1.20.1 HDF5Array_1.26.0           
 [51] BiocGenerics_0.44.0         abind_1.4-5                
 [53] scales_1.2.1                mvtnorm_1.1-3              
 [55] DBI_1.1.3                   rstatix_0.7.2              
 [57] spatstat.random_3.1-4       miniUI_0.1.1.1             
 [59] viridisLite_0.4.1           xtable_1.8-4               
 [61] reticulate_1.28             rsvd_1.0.5                 
 [63] bit_4.0.5                   mclust_6.0.0               
 [65] stats4_4.2.3                htmlwidgets_1.6.2          
 [67] httr_1.4.5                  RColorBrewer_1.1-3         
 [69] ellipsis_0.3.2              ica_1.0-3                  
 [71] pkgconfig_2.0.3             XML_3.99-0.14              
 [73] uwot_0.1.14                 deldir_1.0-6               
 [75] utf8_1.2.3                  tidyselect_1.2.0           
 [77] rlang_1.1.0                 reshape2_1.4.4             
 [55] DBI_1.1.3                   rstatix_0.7.2              
 [79] later_1.3.0                 AnnotationDbi_1.60.2       
 [81] cachem_1.0.7                munsell_0.5.0              
 [83] tools_4.2.3                 cli_3.6.1                  
 [85] RSQLite_2.3.0               generics_0.1.3             
 [87] broom_1.0.4                 ggridges_0.5.4             
 [89] stringr_1.5.0               fastmap_1.1.1              
 [91] goftest_1.2-3               bit64_4.0.5                
 [93] fitdistrplus_1.1-8          purrr_1.0.1                
 [95] RANN_2.6.1                  KEGGREST_1.38.0            
 [97] sparseMatrixStats_1.10.0    future_1.32.0              
 [99] nlme_3.1-162                mime_0.12                  
[101] pracma_2.4.2                hdf5r_1.3.8                
[103] compiler_4.2.3              plotly_4.10.1              
[105] png_0.1-8                   ggsignif_0.6.4             
[107] spatstat.utils_3.0-2        tibble_3.2.1               
[109] stringi_1.7.12              lattice_0.20-45            
[111] Matrix_1.5-3                vctrs_0.6.1                
[113] rhdf5filters_1.10.1         pillar_1.9.0               
[115] lifecycle_1.0.3             spatstat.geom_3.1-0        
[117] lmtest_0.9-40               RcppAnnoy_0.0.20           
[119] data.table_1.14.8           cowplot_1.1.1              
[121] bitops_1.0-7                irlba_2.3.5.1              
[123] httpuv_1.6.9                GenomicRanges_1.50.2       
[125] R6_2.5.1                    promises_1.2.0.1           
[127] KernSmooth_2.23-20          gridExtra_2.3              
[129] IRanges_2.32.0              parallelly_1.35.0          
[131] codetools_0.2-19            MASS_7.3-58.3              
[133] rhdf5_2.42.0                SummarizedExperiment_1.28.0
[135] withr_2.5.0                 sctransform_0.3.5          
[137] S4Vectors_0.36.2            GenomeInfoDbData_1.2.9     
[139] parallel_4.2.3              beachmat_2.14.0            
[141] grid_4.2.3                  tidyr_1.3.0                
[143] DelayedMatrixStats_1.20.0   carData_3.0-5              
[145] MatrixGenerics_1.10.0       Rtsne_0.16                 
[147] spatstat.explore_3.1-0      Biobase_2.58.0             
[149] shiny_1.7.4                