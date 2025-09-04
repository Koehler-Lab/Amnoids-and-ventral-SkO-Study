# Extended Data Figure 8: CellChat for cSKO vs SE + LPM and SE + LPM vs SE at late developmental stage
library(CellChat)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(pbapply)
library(Nebulosa)
library(ComplexHeatmap)
packageVersion("CellChat")
# [1] ‘2.1.2’

label <- "late"
output_dir <- "."
annotation_column <- "CTgroup2"
species <- "human"

sample_column <- NULL
annotation_selected <- NULL
group_column <- "group"
top_n <- 10


group_cmp <- list(
  c("d104 SE + LPM", "d104 cSKO + LPM"),
  c("d104 SE + LPM", "d128 vSKO"),
  c("d104 SE + LPM", "d133 cSKO"),
  c("d104 cSKO + LPM", "d128 vSKO"),
  c("d133 cSKO", "d104 cSKO + LPM"),
  c("d133 cSKO", "d128 vSKO")
)


source(paste("scDown/cellchat_utility_functions.R", sep = ""))
source(paste("scDown/run_cellchatV2.R", sep = ""))


for (cmp in group_cmp) {
  # Conditions for pair wise comparison
  cond_1 <- cmp[1]
  cond_2 <- cmp[2]
  if (cond_1 == cond_2) {
    warning("Skip the condition comparison since two conditions are identical: ", cond_1)
    next()
  }

  cond_1 <- cmp[1]
  cond_2 <- cmp[2]
  if (cond_1 == cond_2) {
    warning("Skip the condition comparison since two conditions are identical: ", cond_1)
    next()
  }

  seurat_obj_cond1 <- readRDS(paste0(output_dir, "/cellchat/rds/seurat_obj_", cond_1, ".rds"))
  cellchat_obj_cond1 <- readRDS(paste0(output_dir, "/cellchat/rds/cellchat_obj_", cond_1, ".rds"))
  # Condition 2
  seurat_obj_cond2 <- readRDS(paste0(output_dir, "/cellchat/rds/seurat_obj_", cond_2, ".rds"))
  cellchat_obj_cond2 <- readRDS(paste0(output_dir, "/cellchat/rds/cellchat_obj_", cond_2, ".rds"))


  dir_cellchat <- output_dir

  condition_col <- group_column
  condition_1 <- cond_1
  condition_2 <- cond_2

  if (!(identical(levels(cellchat_obj_cond1@idents), levels(cellchat_obj_cond2@idents)))) {
    message("Aligning cell types between objects.")
    aligned <- align_cell_labels(cellchat_obj_cond1, cellchat_obj_cond2)
    cellchat_obj_cond1 <- unlist(aligned)[[1]]
    cellchat_obj_cond2 <- unlist(aligned)[[2]]
  }

  # merge cellchat objects from 2 biological conditions, the objects being merged need to have the same cell type annotations
  object_list <- list()
  object_list[condition_1] <- cellchat_obj_cond1
  object_list[condition_2] <- cellchat_obj_cond2
  cellchat <- mergeCellChat(object_list, add.names = names(object_list))

  # record conditions and pathways in comparison
  cond_in_compare <- levels(cellchat@meta$datasets)
  pathways_to_compare <- top_pathways(X1 = cellchat_obj_cond1, X2 = cellchat_obj_cond2, top_n = top_n)
  message("Top pathways to compare ", condition_1, " and ", condition_2, " are calculated.")
  cat(pathways_to_compare, sep = ";\n")



  X <- cellchat


  # Extended Figure 8a. Network heatmap for NOTCH signaling showing overall signaling roles of cell clusters in day 100+ single-cell transcriptomic data.
  # graph specific pathways of interests side by side for visual comparison
  for (pathway in pathways_to_compare) {
    side_by_side_path_compr(dir_cellchat, X, object_list, cond_in_compare, pathway)
  }

  # Extended Figure 8b. Chord diagram comparing NOTCH ligands and receptors
  # between vSkO and cvSkO; DLK1 emerges as a key driver of signaling differences.
  # comparing cvSkO and vSkO NOTCH signaling against that of cSkO, reveal cell-type specific pathway contributions.
  # DEG by differential gene expression
  # define a positive dataset, i.e., the dataset with positive fold change against the other dataset
  pos.dataset <- cond_in_compare[2]
  features.name <- pos.dataset
  # Identify cells in clusters that are present in both groups
  group_clusters <- apply(table(X@idents[[3]], X@meta$datasets), 1, function(x) sum(x != 0))
  valid_clusters <- names(group_clusters)[which(group_clusters == 2)]
  # Subset your CellChat object to remove cells from the clusters not presenting in 2 groups
  X_filtered <- subsetCellChat(X, idents.use = valid_clusters)

  # perform differential expression analysis
  X_filtered <- identifyOverExpressedGenes(X_filtered, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 0.05)
  # map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
  net <- netMappingDEG(X_filtered, features.name = features.name)

  # visualize top up-regulated ligands and top down-regulated ligands with >1 logFC (2-fold change)
  for (logFC in c(0.1, 0.2)) { # used this for d104 cvSKO Vs. d133 cSKO # and used this for d128 vSKO vs. d104 cvSKO
    # extract the ligand-receptor pairs with upregulated ligands in LS
    net.up <- subsetCommunication(X_filtered, net = net, datasets = cond_in_compare[2], ligand.logFC = logFC, receptor.logFC = NULL)
    # extract the ligand-receptor pairs with upregulated ligands in NL, i.e.,downregulated in LS
    net.down <- subsetCommunication(X_filtered, net = net, datasets = cond_in_compare[1], ligand.logFC = -logFC, receptor.logFC = NULL)
    print(rbind(dim(net.up), dim(net.down)))
    net.up0 <- net.up
    net.down0 <- net.down
    up0.LRs.source.target <- c(
      paste0(net.up0$interaction_name, "_", net.up0$source, " -> ", net.up0$target, " (", names(X@net)[1], ")"),
      paste0(net.up0$interaction_name, "_", net.up0$source, " -> ", net.up0$target, " (", names(X@net)[2], ")")
    )
    down0.LRs.source.target <- c(
      paste0(net.down0$interaction_name, "_", net.down0$source, " -> ", net.down0$target, " (", names(X@net)[1], ")"),
      paste0(net.down0$interaction_name, "_", net.down0$source, " -> ", net.down0$target, " (", names(X@net)[2], ")")
    )
    write.csv(net.up, file = paste0(dir_cellchat, "/cellchat/csv/", cond_in_compare[2], "Vs", cond_in_compare[1], "_increased_signalingLR_diffExpession_", logFC, "logFC.csv", sep = ""))
    write.csv(net.down, file = paste0(dir_cellchat, "/cellchat/csv/", cond_in_compare[2], "Vs", cond_in_compare[1], "_decreased_signalingLR_diffExpession_", logFC, "logFC.csv", sep = ""))
    net.up0 <- read.csv(file = paste0(dir_cellchat, "/cellchat/csv/", cond_in_compare[2], "Vs", cond_in_compare[1], "_increased_signalingLR_diffExpession_", logFC, "logFC.csv", sep = ""))
    net.down0 <- read.csv(file = paste0(dir_cellchat, "/cellchat/csv/", cond_in_compare[2], "Vs", cond_in_compare[1], "_decreased_signalingLR_diffExpession_", logFC, "logFC.csv", sep = ""))

    # visualize top NOTCH signaling in bubble plots
    ### select NOTCH signaling ligands only
    signaling <- "NOTCH"
    grep("NOTCH", net.up0$receptor, value = T)
    grep("NOTCH", net.up0$ligand, value = T)
    ligand <- c("DLL1", "DLK1", "JAG1", "JAG2")
    receptor <- c("NOTCH1", "NOTCH2", "NOTCH3")
    LRs_uni <- c(ligand, receptor)
    net.up <- net.up0[which(net.up0$ligand %in% ligand), ]
    net.down <- net.down0[which(net.down0$ligand %in% ligand), ]
    print(rbind(dim(net.up), dim(net.down)))
    table(as.character(net.up$interaction_name))
    table(as.character(net.down$interaction_name))
    pairLR.use.up <- net.up[, "interaction_name", drop = F]
    pairLR.use.down <- net.down[, "interaction_name", drop = F]
    write.csv(net.up, file = paste0(dir_cellchat, "/cellchat/csv/", cond_in_compare[2], "Vs", cond_in_compare[1], "_", signaling, "increased_LR_diffExpession_", logFC, "logFCv2.csv", sep = ""))
    write.csv(net.down, file = paste0(dir_cellchat, "/cellchat/csv/", cond_in_compare[2], "Vs", cond_in_compare[1], "_", signaling, "decreased_LR_diffExpession_", logFC, "logFCv2.csv", sep = ""))

    # visualize top signaling in chord diagram
    png(paste(dir_cellchat, "/cellchat/csv/", cond_in_compare[2], "Vs", cond_in_compare[1], "_LR_diffExpession_", signaling, "_up", logFC, "logFC_chord1.png", sep = ""), height = 500 * 2.5, width = 700 * 2.8, res = 300, pointsize = 10)
    par(mar = c(1, 1, 1, 1), xpd = TRUE)
    netVisual_chord_gene(object_list[[2]], slot.name = "net", net = net.up, lab.cex = 0.8, legend.pos.x = 5, small.gap = 3.5, title.name = paste("Up-regulated", signaling, "signaling in", cond_in_compare[2], "for", cond_in_compare[2], "Vs.", cond_in_compare[1]))
    dev.off()
    png(paste(dir_cellchat, "/cellchat/csv/", cond_in_compare[2], "Vs", cond_in_compare[1], "_LR_diffExpession_", signaling, "_down", logFC, "logFC_chord1.png", sep = ""), height = 500 * 2.5, width = 700 * 2.8, res = 300, pointsize = 10)
    par(mar = c(1, 1, 1, 1), xpd = TRUE)
    netVisual_chord_gene(object_list[[1]], slot.name = "net", net = net.down, lab.cex = 0.8, legend.pos.x = 5, small.gap = 3.5, title.name = paste("Down-regulated", signaling, "signaling in", cond_in_compare[1], "for", cond_in_compare[2], "Vs.", cond_in_compare[1]))
    dev.off()

    library(svglite)
    svglite(paste(dir_cellchat, "/cellchat/csv/", cond_in_compare[2], "Vs", cond_in_compare[1], "_LR_diffExpession_", signaling, "_up", logFC, "logFC_chord1.svg", sep = ""), height = 2.1 * 2.5, width = 3 * 2.8)
    par(mar = c(1, 1, 1, 1), xpd = TRUE)
    netVisual_chord_gene(object_list[[2]], slot.name = "net", net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste("Up-regulated", signaling, "signaling in", cond_in_compare[2], "for", cond_in_compare[2], "Vs.", cond_in_compare[1]))
    dev.off()
    svglite(paste(dir_cellchat, "/cellchat/csv/", cond_in_compare[2], "Vs", cond_in_compare[1], "_LR_diffExpession_", signaling, "_down", logFC, "logFC_chord1.svg", sep = ""), height = 2.1 * 2.5, width = 3 * 2.8)
    par(mar = c(1, 1, 1, 1), xpd = TRUE)
    netVisual_chord_gene(object_list[[1]], slot.name = "net", net = net.down, lab.cex = 0.8, legend.pos.x = 5, small.gap = 3.5, title.name = paste("Down-regulated", signaling, "signaling in", cond_in_compare[1], "for", cond_in_compare[2], "Vs.", cond_in_compare[1]))
    dev.off()
  }
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
 [1] ggalluvial_0.12.5     NMF_0.26              cluster_2.1.4        
 [4] rngtools_1.5.2        registry_0.5-1        ComplexHeatmap_2.14.0
 [7] Nebulosa_1.8.0        pbapply_1.7-0         patchwork_1.1.2      
[10] SeuratDisk_0.0.0.9021 SeuratObject_4.1.3    Seurat_4.3.0         
[13] CellChat_2.1.2        Biobase_2.58.0        BiocGenerics_0.44.0  
[16] ggplot2_3.4.1         igraph_1.4.1          dplyr_1.1.4          
loaded via a namespace (and not attached):
  [1] utf8_1.2.3                  ks_1.14.0                  
  [3] spatstat.explore_3.1-0      reticulate_1.28            
  [5] tidyselect_1.2.0            htmlwidgets_1.6.2          
  [7] BiocParallel_1.32.6         Rtsne_0.16                 
  [9] munsell_0.5.0               codetools_0.2-19           
 [11] ica_1.0-3                   future_1.32.0              
 [13] miniUI_0.1.1.1              withr_2.5.0                
 [15] spatstat.random_3.1-4       colorspace_2.1-0           
 [17] progressr_0.13.0            stats4_4.2.3               
 [19] SingleCellExperiment_1.20.1 ROCR_1.0-11                
 [21] ggsignif_0.6.4              tensor_1.5                 
 [23] listenv_0.9.0               MatrixGenerics_1.10.0      
 [25] GenomeInfoDbData_1.2.9      polyclip_1.10-4            
 [27] bit64_4.0.5                 coda_0.19-4                
 [29] parallelly_1.35.0           vctrs_0.6.5                
 [31] generics_0.1.3              GenomeInfoDb_1.34.9        
 [33] R6_2.5.1                    doParallel_1.0.17          
 [35] clue_0.3-64                 hdf5r_1.3.8                
 [37] DelayedArray_0.24.0         bitops_1.0-7               
 [39] spatstat.utils_3.0-2        cachem_1.0.7               
 [41] promises_1.2.0.1            scales_1.2.1               
 [43] gtable_0.3.3                globals_0.16.2             
 [45] goftest_1.2-3               rlang_1.1.3                
 [47] systemfonts_1.0.4           GlobalOptions_0.1.2        
 [49] splines_4.2.3               rstatix_0.7.2              
 [51] lazyeval_0.2.2              spatstat.geom_3.1-0        
 [53] broom_1.0.4                 BiocManager_1.30.20        
 [55] reshape2_1.4.4              abind_1.4-5                
 [57] ggnetwork_0.5.12            backports_1.4.1            
 [59] httpuv_1.6.9                tools_4.2.3                
 [61] gridBase_0.4-7              statnet.common_4.8.0       
 [63] ellipsis_0.3.2              jquerylib_0.1.4            
 [65] RColorBrewer_1.1-3          ggridges_0.5.4             
 [67] Rcpp_1.0.10                 plyr_1.8.8                 
 [69] zlibbioc_1.44.0             RCurl_1.98-1.12            
 [71] purrr_1.0.1                 ggpubr_0.6.0               
 [73] deldir_1.0-6                GetoptLong_1.0.5           
 [75] cowplot_1.1.1               S4Vectors_0.36.2           
 [77] zoo_1.8-11                  SummarizedExperiment_1.28.0
 [79] ggrepel_0.9.3               magrittr_2.0.3             
 [81] data.table_1.14.8           RSpectra_0.16-1            
 [83] sna_2.7-1                   scattermore_0.8            
 [85] circlize_0.4.15             lmtest_0.9-40              
 [87] RANN_2.6.1                  mvtnorm_1.1-3              
 [89] fitdistrplus_1.1-8          matrixStats_0.63.0         
 [91] mime_0.12                   xtable_1.8-4               
 [93] mclust_6.0.0                IRanges_2.32.0             
 [95] gridExtra_2.3               shape_1.4.6                
 [97] compiler_4.2.3              tibble_3.2.1               
 [99] KernSmooth_2.23-20          crayon_1.5.2               
[101] htmltools_0.5.5             later_1.3.0                
[103] tidyr_1.3.0                 MASS_7.3-58.3              
[105] Matrix_1.5-3                car_3.1-2                  
[107] cli_3.6.1                   parallel_4.2.3             
[109] GenomicRanges_1.50.2        pkgconfig_2.0.3            
[111] sp_1.6-0                    plotly_4.10.1              
[113] spatstat.sparse_3.0-1       foreach_1.5.2              
[115] svglite_2.1.1               bslib_0.4.2                
[117] XVector_0.38.0              stringr_1.5.0              
[119] digest_0.6.31               pracma_2.4.2               
[121] sctransform_0.3.5           RcppAnnoy_0.0.20           
[123] spatstat.data_3.0-1         leiden_0.4.3               
[125] uwot_0.1.14                 shiny_1.7.4                
[127] rjson_0.2.21                lifecycle_1.0.3            
[129] nlme_3.1-162                jsonlite_1.8.4             
[131] carData_3.0-5               network_1.18.1             
[133] BiocNeighbors_1.16.0        viridisLite_0.4.1          
[135] fansi_1.0.4                 pillar_1.9.0               
[137] lattice_0.20-45             fastmap_1.1.1              
[139] httr_1.4.5                  survival_3.5-5             
[141] glue_1.6.2                  FNN_1.1.3.2                
[143] png_0.1-8                   iterators_1.0.14           
[145] bit_4.0.5                   presto_1.0.0               
[147] stringi_1.7.12              sass_0.4.5                 
[149] tidyverse_2.0.0             irlba_2.3.5.1              
[151] future.apply_1.10.0        