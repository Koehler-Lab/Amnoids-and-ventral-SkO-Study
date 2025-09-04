# Extended Data Fig. 9 | Molecular and mechanical features of epithelial–mesenchymal coupling in organoids.
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(pbapply)
library(Nebulosa)
library(ComplexHeatmap)
require(symphony)
source("rFunctions.R")

# load P239.3 with final cell type annotation and finalized sample labeling
X <- readRDS(file = "Annotation/P239.3_final.rds")
DefaultAssay(X) <- "RNA"
Idents(X) <- X$CT
P239.3 <- X


# Extended Figure 9a. Violin plots showing expression of YAP1 targets (ANXA1–3), the YAP1 negative regulator FRMD6, epithelial–mesenchymal transition markers (CCN1, CXCL14, CD24), and actomyosin 28 contractility gene KANK4 across multiple epithelial and mesenchymal populations.
# SEO (red) and cSkO (blue) groups are highlighted. (+ve): positive, (-ve): negative.
# violin plots comparing Aminion and Skin in Epidermis
epicelltypes <- c("AE", "BasalKC1", "cyAE", "EarlyKC", "BasalKC2", "AE/SE1", "SCP", "AE/SE2", "NP", "Cuticle-Cortex", "MEL", "HF-KC", "Periderm", "AE-Epiderm", "NCC", "Endothelial")
genes3 <- c("ANXA1", "ANXA2", "ANXA3", "FRMD6", "CCN1", "CXCL14", "CD24", "KANK4")

dge <- subset(dgeall, type %in% c("cSKO", "SE"))
dge$type <- factor(dge$type, levels = c("SE", "cSKO"))
dge1 <- subset(dge, CT %in% epicelltypes)
dge1$type <- factor(dge1$type, levels = c("SE", "cSKO"))

genes3[!(genes3 %in% rownames(dge1))]
png("P239.3TF-Vln_epi_genes3.png", width = 6667, height = 2000, res = 300)
VlnPlot(dge1, features = genes3, pt.size = -1, split.by = "type", ncol = 4, split = TRUE) & theme(legend.position = "right")
dev.off()
library(svglite)
svglite("P239.3TF-Vln_epi_genes3.svg", width = 25, height = 7.5)
VlnPlot(dge1, features = genes3, pt.size = -1, split.by = "type", ncol = 4, split = TRUE) & theme(legend.position = "right")
dev.off()
png("P239.3TF-Vln_epi_genes3_nolegend.png", width = 5600, height = 2000, res = 300)
VlnPlot(dge1, features = genes3, pt.size = -1, split.by = "type", ncol = 4, split = TRUE) & theme(legend.position = "none")
dev.off()
library(svglite)
svglite("P239.3TF-Vln_epi_genes3_nolegend.svg", width = 21, height = 7.5)
VlnPlot(dge1, features = genes3, pt.size = -1, split.by = "type", ncol = 4, split = TRUE) & theme(legend.position = "none")
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
[1] grid      stats     graphics  grDevices utils     datasets  methods  
[8] base     
other attached packages:
[1] ComplexHeatmap_2.14.0 Nebulosa_1.8.0        pbapply_1.7-0        
[4] patchwork_1.1.2       ggplot2_3.4.4         SeuratDisk_0.0.0.9020
[7] SeuratObject_4.1.3    Seurat_4.3.0         

loaded via a namespace (and not attached):
  [1] Rtsne_0.16                  colorspace_2.1-0           
  [3] rjson_0.2.21                deldir_1.0-6               
  [5] ellipsis_0.3.2              ggridges_0.5.4             
  [7] mclust_6.0.0                circlize_0.4.15            
  [9] XVector_0.38.0              GlobalOptions_0.1.2        
 [11] GenomicRanges_1.50.2        clue_0.3-64                
 [13] spatstat.data_3.0-1         leiden_0.4.3               
 [15] listenv_0.9.0               ggrepel_0.9.3              
 [17] bit64_4.0.5                 mvtnorm_1.1-3              
 [19] fansi_1.0.4                 codetools_0.2-19           
 [21] splines_4.2.3               doParallel_1.0.17          
 [23] polyclip_1.10-4             jsonlite_1.8.4             
 [25] ica_1.0-3                   cluster_2.1.4              
 [27] png_0.1-8                   uwot_0.1.14                
 [29] shiny_1.7.4                 sctransform_0.4.1          
 [31] spatstat.sparse_3.0-1       compiler_4.2.3             
 [33] httr_1.4.5                  Matrix_1.5-3               
 [35] fastmap_1.1.1               lazyeval_0.2.2             
 [37] cli_3.6.1                   later_1.3.0                
 [39] htmltools_0.5.5             tools_4.2.3                
 [41] igraph_2.0.3                GenomeInfoDbData_1.2.9     
 [43] gtable_0.3.3                glue_1.6.2                 
 [45] RANN_2.6.1                  reshape2_1.4.4             
 [47] dplyr_1.1.1                 Rcpp_1.0.10                
 [49] Biobase_2.58.0              scattermore_1.2            
 [51] vctrs_0.6.1                 spatstat.explore_3.1-0     
 [53] nlme_3.1-162                progressr_0.13.0           
 [55] iterators_1.0.14            lmtest_0.9-40              
 [57] spatstat.random_3.1-4       stringr_1.5.0              
 [59] globals_0.16.2              mime_0.12                  
 [61] miniUI_0.1.1.1              lifecycle_1.0.3            
 [63] irlba_2.3.5.1               goftest_1.2-3              
 [65] future_1.32.0               zlibbioc_1.44.0            
 [67] MASS_7.3-58.3               zoo_1.8-11                 
 [69] scales_1.2.1                promises_1.2.0.1           
 [71] MatrixGenerics_1.10.0       spatstat.utils_3.0-2       
 [73] parallel_4.2.3              SummarizedExperiment_1.28.0
 [75] RColorBrewer_1.1-3          SingleCellExperiment_1.20.1
 [77] reticulate_1.28             gridExtra_2.3              
 [79] stringi_1.7.12              S4Vectors_0.36.2           
 [81] foreach_1.5.2               BiocGenerics_0.44.0        
 [83] shape_1.4.6                 GenomeInfoDb_1.34.9        
 [85] bitops_1.0-7                rlang_1.1.4                
 [87] pkgconfig_2.0.3             matrixStats_0.63.0         
 [89] pracma_2.4.2                lattice_0.20-45            
 [91] ROCR_1.0-11                 purrr_1.0.1                
 [93] tensor_1.5                  ks_1.14.0                  
 [95] htmlwidgets_1.6.2           cowplot_1.1.1              
 [97] bit_4.0.5                   tidyselect_1.2.0           
 [99] parallelly_1.35.0           RcppAnnoy_0.0.20           
[101] plyr_1.8.8                  magrittr_2.0.3             
[103] R6_2.5.1                    IRanges_2.32.0             
[105] generics_0.1.3              DelayedArray_0.24.0        
[107] pillar_1.9.0                withr_2.5.0                
[109] fitdistrplus_1.1-8          RCurl_1.98-1.12            
[111] survival_3.5-5              abind_1.4-5                
[113] sp_1.6-0                    tibble_3.2.1               
[115] future.apply_1.10.0         crayon_1.5.2               
[117] hdf5r_1.3.8                 KernSmooth_2.23-20         
[119] utf8_1.2.3                  spatstat.geom_3.1-0        
[121] plotly_4.10.1               GetoptLong_1.0.5           
[123] data.table_1.14.8           digest_0.6.31              
[125] xtable_1.8-4                tidyr_1.3.0                
[127] httpuv_1.6.9                stats4_4.2.3               
[129] munsell_0.5.0               viridisLite_0.4.2     
