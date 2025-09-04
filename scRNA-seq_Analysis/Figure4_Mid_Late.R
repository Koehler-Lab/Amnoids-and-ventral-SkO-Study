# Figure 4: Mesenchymal subtypes control amnion vs skin fates.
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


for (i in 2:3) {
  label <- labels[i]
  cut <- cutoverlap[i]
  P307 <- readRDS(file = paste0("P307", label, "_P239.3Annotated", cut, "pct.rds"))
  Y <- P307


  # genes of interest
  genes <- c("DLK1", "LEF1", "DKK1", "SERPINF1", "FLT1")

  # Figure 4b,d,h. Violin plots of DLK1, LEF1, DKK1, SERPINF1, FLT1 expression on day 28 and day100+ single-cell data from Amnioids versus cSkO fibroblasts.
  Amnioidfibroblasts <- c("earlyAmFibro", "AmFibro1", "AmFibro2")
  cSkOfibroblasts <- c("DCFibro", "FRZBFibro", "PRRXMesProgenitors", "DermFibro", "earlyFibro")
  celltypes <- c(Amnioidfibroblasts, cSkOfibroblasts)

  dge1 <- subset(Y, transferedCellType %in% celltypes)
  group <- as.character(dge1$transferedCellType)
  names(group) <- names(dge1$transferedCellType)

  group[which(group %in% Amnioidfibroblasts)] <- "Amnioid fibroblasts"
  group[which(group %in% cSkOfibroblasts)] <- "cSkO fibroblasts"
  table(group)
  dge1$group <- factor(group, levels = c("Amnioid fibroblasts", "cSkO fibroblasts"))
  Idents(dge1) <- dge1$group

  ncol <- 5
  p <- VlnPlot(dge1, features = genes, pt.size = -1, ncol = ncol) # ,split.by="group",split=TRUE
  plotlist <- patchwork::wrap_plots(p) + plot_layout(guides = "collect") &
    theme(legend.position = "none", plot.margin = margin(5, 5, 5, 10)) # , , down,left
  png(paste0("P307", label, "_P239.3Annotated", cut, "pct1_Vln_5genes_2groups.png"), width = 3000, height = 1200, res = 300)
  plotlist
  dev.off()
  library(svglite)
  svglite(paste0("P307", label, "_P239.3Annotated", cut, "pct1_Vln_5genes_2groups.svg"), height = 4.5, width = 11.25)
  plotlist
  dev.off()


  # Figure 4c,e,i. Feature plots of DLK1, LEF1, DKK1, SERPINF1, FLT1 expression on day 28 and day100+ single-cell data
  summary(FetchData(Y, vars = genes))
  maxexp <- max(FetchData(Y, vars = genes))
  maxexp # [1] 6.071071
  nrow <- 1
  ncol <- 5
  plots <- lapply(genes, function(gene) {
    FeaturePlot(Y, features = gene, pt.size = 0.5, cols = c("grey90", "#DB5545")) + NoLegend() + NoAxes() +
      ggtitle(bquote(italic(.(gene)))) +
      theme(plot.title = element_text(size = 28)) # , pt.size=0.5
  })
  plotlist <- cowplot::plot_grid(plotlist = plots, ncol = ncol)
  png(paste0("P307", label, "_P239.3Annotated", cut, "pct1_Feature_5genes.png"), res = 300, height = 850 * nrow, width = 800 * ncol)
  print(plotlist)
  dev.off()
  library(svglite)
  svglite(paste0("P307", label, "_P239.3Annotated", cut, "pct1_Feature_5genes.svg"), height = 3.1875 * nrow, width = 3 * ncol)
  print(plotlist)
  dev.off()


  plots <- lapply(genes, function(gene) {
    FeaturePlot(Y, features = gene, pt.size = 0.5, cols = c("grey90", "#DB5545")) + NoAxes() +
      ggtitle(bquote(italic(.(gene)))) +
      theme(plot.title = element_text(size = 28)) # , pt.size=0.5
  })
  plotlist <- cowplot::plot_grid(plotlist = plots, ncol = ncol)
  png(paste0("P307", label, "_P239.3Annotated", cut, "pct1_Feature_5genes_legend.png"), res = 300, height = 850 * nrow, width = 900 * ncol)
  print(plotlist)
  dev.off()



  plots <- lapply(genes, function(gene) {
    FeaturePlot(Y, features = gene, pt.size = 0.5) + NoLegend() + NoAxes() +
      scale_color_gradientn(colours = c("grey90", "#DB5545"), limits = c(0, maxexp)) +
      ggtitle(bquote(italic(.(gene)))) +
      theme(plot.title = element_text(size = 28)) #
  })
  plotlist <- cowplot::plot_grid(plotlist = plots, ncol = ncol)
  png(paste0("P307", label, "_P239.3Annotated", cut, "pct1_Feature_5genes_sameColorScale.png"), res = 300, height = 850 * nrow, width = 800 * ncol)
  print(plotlist)
  dev.off()

  plots <- lapply(genes, function(gene) {
    FeaturePlot(Y, features = gene, pt.size = 0.5) + NoAxes() +
      scale_color_gradientn(colours = c("grey90", "#DB5545"), limits = c(0, maxexp)) +
      ggtitle(bquote(italic(.(gene)))) +
      theme(plot.title = element_text(size = 28)) # , pt.size=0.5
  })
  plotlist <- patchwork::wrap_plots(plots, ncol = ncol) + plot_layout(guides = "collect") &
    theme(legend.position = "right")
  png(paste0("P307", label, "_P239.3Annotated", cut, "pct1_Feature_5genes_sameColorScale_legend.png"), res = 300, height = 850 * nrow, width = 850 * ncol)
  print(plotlist)
  dev.off()


  # For DKK1 and FLT1, make the cells expressing high levels of the marker more obvious
  # order cells expressing higher levels of the marker on top
  genes_to_saturate <- c("LEF1", "DKK1", "FLT1")
  plots <- lapply(genes, function(gene) {
    p <- FeaturePlot(Y,
      features = gene, pt.size = 0.5, cols = c("grey90", "#DB5545"),
      order = if (gene %in% genes_to_saturate) TRUE else FALSE # ensures high-expressing cells are plotted on top
    ) + NoAxes() +
      ggtitle(bquote(italic(.(gene)))) +
      theme(plot.title = element_text(size = 28), legend.position = "right")
    return(p)
  })
  #
  plotlist <- cowplot::plot_grid(plotlist = plots, ncol = ncol)
  png(paste0("P307", label, "_P239.3Annotated", cut, "pct1_Feature_5genes_order.png"), res = 300, height = 850 * nrow, width = 900 * ncol)
  print(plotlist)
  dev.off()


  # cap the color q80: all cells with expression in the top 20% (e.g., ≥ 80th percentile) are shown with the same maximal color, effectively capping the scale.
  genes_to_saturate <- c("LEF1", "DKK1", "FLT1")
  plots <- lapply(genes, function(gene) {
    p <- FeaturePlot(Y,
      features = gene, pt.size = 0.5, cols = c("grey90", "#DB5545"),
      max.cutoff = if (gene %in% genes_to_saturate) "q80" else NA # override default NA by q80 for certain markers
    ) + NoAxes() +
      ggtitle(bquote(italic(.(gene)))) +
      theme(plot.title = element_text(size = 28), legend.position = "right")
    return(p)
  })
  plotlist <- cowplot::plot_grid(plotlist = plots, ncol = ncol)
  png(paste0("P307", label, "_P239.3Annotated", cut, "pct1_Feature_5genes_top20pct.png"), res = 300, height = 850 * nrow, width = 900 * ncol)
  print(plotlist)
  dev.off()

  # cap the color q70: all cells with expression in the top 30% (e.g., ≥ 70th percentile) are shown with the same maximal color, effectively capping the scale.
  genes_to_saturate <- c("LEF1", "DKK1", "FLT1")
  plots <- lapply(genes, function(gene) {
    p <- FeaturePlot(Y,
      features = gene, pt.size = 0.5, cols = c("grey90", "#DB5545"),
      max.cutoff = if (gene %in% genes_to_saturate) "q70" else NA # override default NA by q80 for certain markers
    ) + NoAxes() +
      ggtitle(bquote(italic(.(gene)))) +
      theme(plot.title = element_text(size = 28), legend.position = "right")
    return(p)
  })
  plotlist <- cowplot::plot_grid(plotlist = plots, ncol = ncol)
  png(paste0("P307", label, "_P239.3Annotated", cut, "pct1_Feature_5genes_top30pct.png"), res = 300, height = 850 * nrow, width = 900 * ncol)
  print(plotlist)
  dev.off()
}


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