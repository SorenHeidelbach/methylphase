#!/usr/bin/env Rscript
pacman::p_load(
  "data.table",
  "dbscan",
  "patchwork"
)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggside)
  library(glue)
  library(readr)
  library(tidyr)
  library(purrr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (interactive() & length(args) == 0) {
  args <- c(
    "../development/methylation_phasing/ZymoFecal_bin2.474_contig_1330/per_red.tsv",
    "../development/methylation_phasing/ZymoFecal_bin2.474_contig_1330/aggregate.tsv",
    "output_ZymoFecal_bin2.474_contig_1330"
  )
}
if (interactive() & length(args) == 0) {
  args <- c(
    "../development/methylation_phasing/ZymoFecal_bin2.474_contig_1330/per_red.tsv",
    "../development/methylation_phasing/ZymoFecal_bin2.474_contig_1330/aggregate.tsv",
    "output_ZymoFecal_bin2.474_contig_1330"
  )
}

if (length(args) < 3) {
  stop("usage: analyze_methylation.R <per_read.tsv> <aggregate.tsv> <output_dir>")
}

per_read_path <- args[[1]]
aggregate_path <- args[[2]]
out_dir <- args[[3]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


motif_order <- function(df) {
  df %>%
    count(motif, sort = TRUE) %>%
    pull(motif)
}

save_plot <- function(plot, filename, width = 10, height = 6) {
  ggsave(
    filename = file.path(out_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

per_reads <- fread(per_read_path)
aggregate_reads <- fread(aggregate_path)

motif_levels <- motif_order(per_reads)
per_reads <- per_reads %>% mutate(motif = factor(motif, levels = motif_levels))
aggregate_reads <- aggregate_reads %>% mutate(motif = factor(motif, levels = motif_levels))


aggregate_reads_wide <- aggregate_reads %>% select(read_id, motif, mean_probability) %>% 
  pivot_wider(
  , names_from = motif, values_from = mean_probability
) %>%  na.omit()


clust_res <- hdbscan(aggregate_reads_wide %>% select(-read_id), minPts = 5)
aggregate_reads_wide <- aggregate_reads_wide %>%  mutate(
  cluster = clust_res$cluster,
  cluster_prob = clust_res$membership_prob,
)

ggplot(aggregate_reads_wide) +
  aes(
    x = GATC_6mA_2,
    y = GCWGC_5mC_2,
    fill = as.factor(cluster),
    alpha = cluster_prob
  ) +
  geom_point(shape = 21, size = 3) +
  theme_bw() +
  labs(
    x = "G6mATC, mean of probabilities in read",
    y = "G5mCWGC, mean of probabilities in read"
  )








temp <- fread("//wsl$/nix/home/shei/methylation_phasing/development/methylation_phasing/results/Escherichia-coli_DSMZ18039_split/read_clustering_raw.tsv")
ggplot(temp) +
  aes(
    x = TGATC_6mA_3,
    y = GCWGC_5mC_2,
    color = as.factor(cluster_id)
  ) +
  geom_point()



plot_raw <- temp  %>% pivot_longer(cols = -c("read_id", "cluster_id"))  %>%
  mutate(
    read_id = factor(read_id)
  ) %>% 
  ggplot() +
  aes(
    x = name,
    y = read_id,
    fill = value
  ) +
  geom_tile() +
  geom_point(
    inherit.aes = FALSE,
    x = 0.45,
    shape = 18,
    aes(col = as.factor(cluster_id), y = read_id),
    size = 2
  ) +
  scale_fill_gradientn(colors = rev(c("#FFC20A" , "#0C7BDC"))) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = "Motifs identified in bin",
    y = "Read id",
    col = "Cluster",
    fill = "Mean methylation \nprobabilities",
    title = "Raw"
  )

temp <- fread("//wsl$/nix//home/shei/methylation_phasing/development/methylation_phasing/results/Escherichia-coli_DSMZ18039_split/read_clustering.tsv")
ggplot(temp) +
  aes(
    x = GATC_6mA_1,
    y = CCWGG_5mC_1,
    col = as.factor(cluster_id)
  ) +
  geom_point(
    alpha = 0.05
  )



plot_impute <- temp  %>% pivot_longer(cols = -c("read_id", "cluster_id"))  %>%
  mutate(
    read_id = factor(read_id, levels = temp$read_id[y_order])
  ) %>% 
  mutate(
    cluster_x = "cluster"
  ) %>% 
  ggplot() +
  aes(
    x = name,
    y = read_id,
    fill = value
  ) +
  geom_tile() +
  geom_point(
    inherit.aes = FALSE,
    x = 0.45,
    shape = 18,
    aes(col = as.factor(cluster_id), y = read_id),
    size = 2
  ) +
  scale_fill_gradientn(colors = rev(c("#FFC20A" , "#0C7BDC"))) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = "Motifs identified in bin",
    y = "Read id",
    col = "Cluster",
    fill = "Mean methylation \nprobabilities",
    title = "Imputed"
  )


 plot_impute + plot_raw +
   plot_layout(guides = 'collect', axes = "collect")

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 longitudinal <- fread("../development/methylation_phasing/Escherichia-coli_DSMZ18039_longitudinal/longitudinal.tsv")
 
 
 
( ggplot(longitudinal) +
   aes(
     x = block_start,
     y = polarization
   ) +
   geom_point()) /
 (ggplot(longitudinal) +
   aes(
     x = block_start,
     y = methylation_variance
   ) +
   geom_point())
 
 
 ggplot(longitudinal) +
   aes(
     x = block_start
   ) +
   geom_line(
     aes( y = fully_methylated_reads + mixed_reads),
     color = "blue3"
   ) +
   geom_line(
     aes(y = fully_unmethylated_reads),
     color = "red3"
   )   
 
 
 ggplot(longitudinal) +
   aes(
     x = block_start
   ) +
   geom_line(
     aes( y = mean_methylation),
     color = "blue3"
   ) 
 
 
 
 