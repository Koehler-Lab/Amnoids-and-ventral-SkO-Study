# Figure 3: Ventralization of skin organoids by size control and assembly.
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(pbapply)
library(Nebulosa)
library(ComplexHeatmap)
source("rFunctions.R")


# load P307 annotated by P239.3 with cell types present in at least a certain percentage of the developmental stage
labels <- c("early", "mid", "late")
cutoverlap <- c(0.14, 0.1, 0.1) # 0.1 used for mid and late stages

i <- 3
label <- labels[i]
cut <- cutoverlap[i]
# use UMAP transferred from Project 239.3 - used this
P307 <- readRDS(file = paste0("P307", label, "_P239.3Annotated", cut, "pct.rds"))
Y <- P307
Ylate <- Y
Y2 <- Y

# decided to use col2 based on email from Phuong on 4/30/204
library(RColorBrewer)
col2 <- c(
  brewer.pal(9, "YlOrBr")[2:5],
  brewer.pal(9, "Oranges")[6],
  brewer.pal(7, "Reds")[-1],
  brewer.pal(9, "Blues"),
  brewer.pal(4, "GnBu")[3:4],
  brewer.pal(4, "Greens"),
  brewer.pal(5, "Purples"),
  brewer.pal(3, "PuRd")[3]
)
length(col2)
names(col2) <- levels(P239.3$CT)
matchColors2 <- col2

gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


# Figure 3g. UMAP of scRNA-seq data from day 128 vSkO reveals distinct dermal and epidermal clusters.
dge <- subset(Y2, sample2 %in% "d128 vSKO")
png(paste0("P307", label, "_P239.3Annotated", cut, "pct_d128vSKO_noaxes_biggerlabel.png"), height = 1600, width = 2000, res = 300)
UMAPPlot(dge, label = TRUE, label.size = 4.5, repel = TRUE) +
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 0), axis.text.y = element_text(size = 0),
    axis.ticks.x = element_blank(), axis.ticks.y = element_blank()
  ) +
  xlab("") +
  ylab("") +
  scale_color_manual(values = matchColors2[levels(Idents(Y2))])
dev.off()
svglite(paste0("P307", label, "_P239.3Annotated", cut, "pct_d128vSKO_noaxes_biggerlabel.svg"), width = 7.5, height = 6)
UMAPPlot(dge, label = TRUE, label.size = 4.5, repel = TRUE) +
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 0), axis.text.y = element_text(size = 0),
    axis.ticks.x = element_blank(), axis.ticks.y = element_blank()
  ) +
  xlab("") +
  ylab("") +
  scale_color_manual(values = matchColors2[levels(Idents(Y2))])
dev.off()


# Figure 3h. Heatmap of selected differentially expressed genes comparing vSkO and cSkO, highlighting ventral-specific transcriptional signatures in vSkOs.
# Heatmap of vSKO genes in all late stage samples by samples and by clusters
library(ComplexHeatmap)
genes <- gene <- c("ZSWIM6", "ANKUB1", "OOEP", "CFAP161", "LINC00486")
X <- Y2

### by samples
Idents(X) <- X$sample2
centroid <- log(AverageExpression(X)$RNA + 1)
centroid2.std <- (centroid - apply(centroid, 1, mean)) / apply(centroid, 1, sd)
centroid2.std <- na.omit(centroid2.std)
cc2.std <- centroid2.std[which(rownames(centroid2.std) %in% genes), ]
cc2.std <- cc2.std[genes2, ]

png(paste0("P307", label, "_P239.3Annotated", cut, "pct_vSKOgenes_Heatmap_sample1.png"), width = 800, height = 600, res = 300)
Heatmap(cc2.std,
  name = "StdExp",
  cluster_columns = F,
  row_title_rot = 0, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8)
)
dev.off()
png(paste0("P307", label, "_P239.3Annotated", cut, "pct_vSKOgenes_Heatmap_sample_cluster.png"), width = 800, height = 700, res = 300)
Heatmap(cc2.std,
  name = "StdExp",
  row_title_rot = 0, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8)
)
dev.off()
png(paste0("P307", label, "_P239.3Annotated", cut, "pct_vSKOgenes_Heatmap_sample_noCluster.png"), width = 700, height = 600, res = 300)
Heatmap(cc2.std,
  name = "StdExp",
  cluster_rows = F,
  cluster_columns = F,
  row_title_rot = 0, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8)
)
dev.off()

library(svglite)
svglite(paste0("P307", label, "_P239.3Annotated", cut, "pct_vSKOgenes_Heatmap_sample_noCluster.svg"), width = 3.0625, height = 2.625)
Heatmap(cc2.std,
  name = "StdExp",
  cluster_rows = F,
  cluster_columns = F,
  row_title_rot = 0, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8)
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
