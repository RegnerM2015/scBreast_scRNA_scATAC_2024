################################################################################
# Matt Regner
# Franco Lab
# Ran in RStudio OnDemand environment, not Docker container
# Corrected dplyr namespace issue, confirm group by operations are working as intended
# Debug docker container with installation of homer and marge
################################################################################
# print(.libPaths())
# .libPaths(c("/usr/local/lib/R/site-library",
#             "/usr/local/lib/R/library",
#             "/home/regnerm/R/x86_64-pc-linux-gnu-library/4.1"))
# devtools::install_github('robertamezquita/marge', ref = 'master',force = TRUE)
# check_homer()

# old_path <- Sys.getenv("PATH")
# Sys.setenv(PATH = paste(old_path, "/home/regnerm/Software/homer/bin/homer", sep = ":"))
# options('homer_path' = '/home/regnerm/Software/homer/bin/homer')

library(ExperimentHub)
library(dplyr)
library(cowplot)
library(patchwork)
library(BSDA)
library(data.table)
library(ggplot2)
library(stringr)
library(marge)
check_homer()
library(purrr)
library(tidyr)
library(ArchR)
library(ComplexHeatmap)
library(circlize)

# Ran in RStudio OnDemand environment, not Docker container

homerHeatmaps <- function(peakPaths,
                          subtype,
                          outpath,
                          seed){
  
  regions <- vector(mode = "list", length = length(peakPaths))
  for(i in 1:length(peakPaths)){
    regions[[i]] <- fread(paste0(subtype,"_P2G_BED_files/",peakPaths[i]))
  }
  
  tbl <- tibble(id=gsub(pattern="_",replacement = "/",x=gsub("^.{0,16}", "", gsub('.{9}$', '', peakPaths))),
                regions=regions,
                path=paste0(outpath,"/homer_result_",gsub("^.{0,16}", "", gsub('.{9}$', '', peakPaths))))
  
  print("Starting Homer motif analysis for each peak set")
  
  pwalk(
    ## Varying parameters (regions, output directory)
    .l = list(x = tbl$regions, path = tbl$path),
    ## Function to find motifs across genome
    .f = find_motifs_genome,
    ## Constant parameters (motif search settings)
    genome = 'hg38',
    scan_size = 50,
    optimize_count = 2,
    only_known = TRUE,
    cores = 2, cache = 100,
    keep_minimal = TRUE, overwrite = TRUE
  )
  
  tbl_results <- tbl %>%
    mutate(known_results = map(path, read_known_results)) %>%
    select(-path, -regions) %>%
    unnest()
  
  tbl_results <- tbl_results[,colnames(tbl_results) %ni% c("motif_pwm")]
  
  fwrite(tbl_results,paste0(outpath,"/homer_motif_results_aggregated.csv"))
  
  # Heatmap of motif names
  tbl_motif_name <- tbl_results %>%
    dplyr::group_by(motif_name,id) %>%
    dplyr::summarise(top_log_p_value = max(log_p_value)) %>%
    ungroup() 
  
  res <- as.data.frame(pivot_wider(tbl_motif_name, names_from = "id", values_from = "top_log_p_value"))
  res <- res[!is.na(res$motif_name),]
  rownames(res) <- res[,1]
  res <- res[,-1]
  res <- res[rowSums(res) >0,]
  
  pdf(paste0(outpath,"/Heatmap_motif_names.pdf"),width=6,height = 8)
  ht <- Heatmap(data.matrix(res),
          col= colorRamp2(c(0, max(res)), c("#FFF5EB", "darkred")))
  draw(ht)
  dev.off()
  
  for( i in c(0.10,0.15,0.20)){
    pdf(paste0(outpath,"/Heatmap_motif_names-top_",i*100,"pct_most_variable.pdf"),width=6,height = 8)
    ht <- Heatmap(data.matrix(res[head(order(rowVars(data.matrix(res),useNames=TRUE), decreasing = TRUE), round(nrow(res)*i)),]),
                  col= colorRamp2(c(0, max(res)), c("#FFF5EB", "darkred")))
    draw(ht)
    dev.off()
  }

  # Heatmap of motif families 
  tbl_motif_family <- tbl_results %>%
    dplyr::group_by(motif_family,id) %>%
    dplyr::summarise(top_log_p_value = max(log_p_value)) %>%
    ungroup() 
  
  res <- as.data.frame(pivot_wider(tbl_motif_family, names_from = "id", values_from = "top_log_p_value"))
  res <- res[!is.na(res$motif_family),]
  rownames(res) <- res[,1]
  res <- res[,-1]
  res <- res[rowSums(res) >0,]
  
  pdf(paste0(outpath,"/Heatmap_motif_families.pdf"),width=6,height = 8)
  ht <- Heatmap(data.matrix(res),
          col= colorRamp2(c(0, max(res)), c("#FFF5EB", "darkred")))
  draw(ht)
  dev.off()
  
  for( i in c(0.10,0.15,0.20)){
    pdf(paste0(outpath,"/Heatmap_motif_families-top_",i*100,"pct_most_variable.pdf"),width=6,height = 8)
    ht <- Heatmap(data.matrix(res[head(order(rowVars(data.matrix(res),useNames=TRUE), decreasing = TRUE), round(nrow(res)*i)),]),
                  col= colorRamp2(c(0, max(res)), c("#FFF5EB", "darkred")))
    draw(ht)
    dev.off()
  }

  # Downsample peaks
  min <- min(unlist(lapply(regions, nrow)))
  for(i in 1:length(regions)){
    print(nrow(regions[[i]]))
    set.seed(seed)
    regions[[i]] <- regions[[i]][sample(1:nrow(regions[[i]]),min),]
  }
  
  tbl <- tibble(id=gsub(pattern="_",replacement = "/",x=gsub("^.{0,16}", "", gsub('.{9}$', '', peakPaths))),
                regions=regions,
                path=paste0(outpath,"/homer_result_downsample",gsub("^.{0,16}", "", gsub('.{9}$', '', peakPaths))))
  
  print("Starting Homer motif analysis for each peak set downsampled to lowest size")
  
  pwalk(
    ## Varying parameters (regions, output directory)
    .l = list(x = tbl$regions, path = tbl$path),
    ## Function to find motifs across genome
    .f = find_motifs_genome,
    ## Constant parameters (motif search settings)
    genome = 'hg38',
    scan_size = 50,
    optimize_count = 2,
    only_known = TRUE,
    cores = 2, cache = 100,
    keep_minimal = TRUE, overwrite = TRUE
  )
  
  tbl_results <- tbl %>%
    mutate(known_results = map(path, read_known_results)) %>%
    select(-path, -regions) %>%
    unnest()
  
  tbl_results <- tbl_results[,colnames(tbl_results) %ni% c("motif_pwm")]
  
  fwrite(tbl_results,paste0(outpath,"/homer_motif_results_aggregated-downsample.csv"))
  
  # Heatmap of motif names
  tbl_motif_name <- tbl_results %>%
    dplyr::group_by(motif_name, id) %>%
    dplyr::summarise(top_log_p_value = max(log_p_value)) %>%
    ungroup() 
  
  res <- as.data.frame(pivot_wider(tbl_motif_name, names_from = "id", values_from = "top_log_p_value"))
  res <- res[!is.na(res$motif_name),]
  rownames(res) <- res[,1]
  res <- res[,-1]
  res <- res[rowSums(res) >0,]
  
  pdf(paste0(outpath,"/Heatmap_motif_names-downsample.pdf"),width=6,height = 8)
  ht <- Heatmap(data.matrix(res),
          col= colorRamp2(c(0, max(res)), c("#FFF5EB", "darkred")))
  draw(ht)
  dev.off()
  
  for( i in c(0.10,0.15,0.20)){
    pdf(paste0(outpath,"/Heatmap_motif_names-downsample-top_",i*100,"pct_most_variable.pdf"),width=6,height = 8)
    ht <- Heatmap(data.matrix(res[head(order(rowVars(data.matrix(res),useNames=TRUE), decreasing = TRUE), round(nrow(res)*i)),]),
                  col= colorRamp2(c(0, max(res)), c("#FFF5EB", "darkred")))
    draw(ht)
    dev.off()
  }
  
  # Heatmap of motif families 
  tbl_motif_family <- tbl_results %>%
    dplyr::group_by(motif_family,id) %>%
    dplyr::summarise(top_log_p_value = max(log_p_value)) %>%
    ungroup() 
  
  res <- as.data.frame(pivot_wider(tbl_motif_family, names_from = "id", values_from = "top_log_p_value"))
  res <- res[!is.na(res$motif_family),]
  rownames(res) <- res[,1]
  res <- res[,-1]
  res <- res[rowSums(res) >0,]
  
  pdf(paste0(outpath,"/Heatmap_motif_families-downsample.pdf"),width=6,height = 8)
  ht <- Heatmap(data.matrix(res),
          col= colorRamp2(c(0, max(res)), c("#FFF5EB", "darkred")))
  draw(ht)
  dev.off()
  
  for( i in c(0.10,0.15,0.20)){
    pdf(paste0(outpath,"/Heatmap_motif_families-downsample-top_",i*100,"pct_most_variable.pdf"),width=6,height = 8)
    ht <- Heatmap(data.matrix(res[head(order(rowVars(data.matrix(res),useNames=TRUE), decreasing = TRUE), round(nrow(res)*i)),]),
                  col= colorRamp2(c(0, max(res)), c("#FFF5EB", "darkred")))
    draw(ht)
    dev.off()
  }
  
}

# Basal
homerHeatmaps(peakPaths=list.files(path="Basal_P2G_BED_files")[-c(6,9)],
              subtype="Basal",
              outpath="Basal_P2G_homer_motif_outputs",
              seed=123)

homerHeatmaps(peakPaths=list.files(path="Basal_P2G_BED_files")[c(6,9)],
              subtype="Basal",
              outpath="Basal_P2G_homer_motif_outputs-upgenes",
              seed=123)

# Luminal
homerHeatmaps(peakPaths=list.files(path="Luminal_P2G_BED_files")[-c(6,9)],
              subtype="Luminal",
              outpath="Luminal_P2G_homer_motif_outputs",
              seed=123)

homerHeatmaps(peakPaths=list.files(path="Luminal_P2G_BED_files")[c(6,9)],
              subtype="Luminal",
              outpath="Luminal_P2G_homer_motif_outputs-upgenes",
              seed=123)
