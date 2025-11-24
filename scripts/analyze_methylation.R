#!/usr/bin/env Rscript

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
    "../development/methylation_phasing/ZymoFecal_bin2.474_contig_1330_split/per_red.tsv",
    "../development/methylation_phasing/ZymoFecal_bin2.474_contig_1330_split/aggregate.tsv",
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

load_per_reads <- function(path) {
  read_tsv(
    path,
    col_types = cols(
      read_id = col_character(),
      contig = col_character(),
      motif = col_character(),
      motif_start = col_double(),
      motif_position = col_double(),
      read_position = col_double(),
      ref_position = col_character(),
      probability = col_character(),
      methylated = col_character(),
      mod_label = col_character()
    )
  ) %>%
    mutate(
      probability = as.numeric(probability),
      ref_position = suppressWarnings(as.numeric(ref_position))
    )
}

load_aggregate_reads <- function(path) {
  read_tsv(
    path,
    col_types = cols(
      read_id = col_character(),
      motif = col_character(),
      call_count = col_double(),
      mean_probability = col_double()
    )
  )
}

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

per_reads <- load_per_reads(per_read_path)
aggregate_reads <- load_aggregate_reads(aggregate_path)

motif_levels <- motif_order(per_reads)
per_reads <- per_reads %>% mutate(motif = factor(motif, levels = motif_levels))
aggregate_reads <- aggregate_reads %>% mutate(motif = factor(motif, levels = motif_levels))

# Histogram of per-site probabilities -----------------------------------------------------------
per_site_hist <- per_reads %>%
  filter(!is.na(probability)) %>%
  mutate(
    , mapped = !is.na(ref_position)
  ) %>% 
  ggplot(aes(probability, fill = mapped)) +
  geom_histogram(binwidth = 0.05, alpha = 0.8) +
  facet_wrap(~motif, scales = "free_y") +
  labs(
    title = "Per-site methylation probability",
    x = "Probability",
    y = "Count"
  ) +
  theme_minimal()

save_plot(per_site_hist, "per_site_probability_distribution.png")

# Histogram of per-read means -------------------------------------------------------------------
per_read_hist <- aggregate_reads %>%
  mutate(
    sites = case_when(
      call_count < 5 ~ "<5",
      call_count <= 10 ~ "5-10",
      call_count <= 20 ~ "11-20",
      call_count > 20 ~ ">20"
    )
  ) %>% 
  ggplot(aes(mean_probability, fill = sites)) +
  geom_histogram(binwidth = 0.05, alpha = 0.8) +
  facet_wrap(~motif, scales = "free_y") +
  labs(
    title = "Per-read mean methylation probability",
    x = "Mean probability",
    y = "Count"
  ) +
  theme_minimal()

save_plot(per_read_hist, "per_read_mean_probability.png")

# Contingency tables ----------------------------------------------------------------------------
threshold <- 0.5

contingency <- per_reads %>%
  mutate(
    methylated = case_when(
      is.na(probability) ~ "missing",
      probability >= threshold ~ "high",
      TRUE ~ "low"
    )
  ) %>%
  count(motif, methylated, name = "site_count") %>%
  group_by(motif) %>%
  mutate(fraction = site_count / sum(site_count)) %>%
  ungroup()

write_tsv(contingency, file.path(out_dir, "motif_methylation_contingency.tsv"))

contingency_plot <- contingency %>%
  mutate(methylated = factor(methylated, levels = c("high", "low", "missing"))) %>%
  ggplot(aes(x = motif, y = fraction, fill = methylated)) +
  geom_col(position = "stack") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = glue("Per-site contingency (threshold = {threshold})"),
    x = "Motif",
    y = "Fraction",
    fill = "Category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(contingency_plot, "motif_methylation_contingency.png", width = 9, height = 5)

read_summary <- aggregate_reads %>%
  mutate(category = ifelse(mean_probability >= threshold, "mean >= threshold", "mean < threshold")) %>%
  count(motif, category, name = "read_count") %>%
  group_by(motif) %>%
  mutate(fraction = read_count / sum(read_count)) %>%
  ungroup()

write_tsv(read_summary, file.path(out_dir, "motif_read_contingency.tsv"))

read_contingency_plot <- read_summary %>%
  ggplot(aes(x = motif, y = fraction, fill = category)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = glue("Per-read contingency (threshold = {threshold})"),
    x = "Motif",
    y = "Fraction",
    fill = "Category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(read_contingency_plot, "per_read_contingency.png", width = 9, height = 5)

# Reference chunk distribution ------------------------------------------------------------------
max_ref <- suppressWarnings(max(per_reads$ref_position, na.rm = TRUE))
chunk_size <- if (is.finite(max_ref) && max_ref > 0) {
  max(3000, round(max_ref / 25, -3))
} else {
  10000
}

chunk_plot <- per_reads %>%
  filter(!is.na(probability), !is.na(ref_position)) %>%
  mutate(
    chunk = floor(ref_position / chunk_size) * chunk_size,
    chunk_label = factor(
      glue("{chunk / 1000}-{(chunk + chunk_size) / 1000} kb"),
      levels = glue("{sort(unique(chunk)) / 1000}-{(sort(unique(chunk)) + chunk_size) / 1000} kb")
    )
  ) %>%
  ggplot(aes(x = chunk_label, y = probability, fill = motif, gorup = motif)) +
  #geom_boxplot(outlier.alpha = 0.2, outlier.size = 0.5) +
  geom_violin(
    scale = "width", 
    draw_quantiles = c(0.25, 0.5, 0.75), 
    quantile.linewidth = c(1, 2, 1)/2,
    quantile.colour = c("gray30", "gray0", "gray30")) +
  facet_grid(strand ~contig, scales = "free_x") +
  labs(
    title = glue("Per-site methylation across {chunk_size/1000} kb reference chunks"),
    x = "Reference chunk",
    y = "Probability",
    fill = "Motif"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
chunk_plot
save_plot(chunk_plot, "chunk_methylation_boxplot.png", width = 12, height = 6)



chunk_aggregate_plot <- per_reads %>%
  filter(!is.na(probability), !is.na(ref_position)) %>%
  mutate(
    chunk = floor(ref_position / chunk_size) * chunk_size,
    chunk_label = factor(
      glue("{chunk / 1000}-{(chunk + chunk_size) / 1000} kb"),
      levels = glue("{sort(unique(chunk)) / 1000}-{(sort(unique(chunk)) + chunk_size) / 1000} kb")
    )
  ) %>%
  group_by(
    chunk, contig, chunk_label, motif
  ) %>% 
  summarise(
    quantiles = quantile(probability, c(c(1, 2.5, 5, 7.5, 9)/10)),
    quantile_lab = as.factor(c(1, 2.5, 5, 7.5, 9)/10)
  ) %>% 
  ggplot(aes(x = chunk_label, y = quantiles, fill = quantile_lab, group = quantile_lab, color = quantile_lab )) +
  geom_line(aes(col = quantile_lab)) +
  geom_point(size = 5) +
  facet_wrap(motif~contig, scales = "free_x") +
  labs(
    title = glue("Per-site methylation across {chunk_size/1000} kb reference chunks"),
    x = "Reference chunk",
    y = "Probability",
    fill = "Motif"
  ) +
  scale_color_viridis_d(direction = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

chunk_aggregate_plot
save_plot(chunk_aggregate_plot, "chunk_methylation_quantile_plot.png", width = 12, height = 6)
# Pairwise motif relationships ------------------------------------------------------------------
generate_pair_plots <- function(agg_df, limit = 3) {
  motifs <- agg_df %>%
    count(motif, sort = TRUE) %>%
    head(limit + 1) %>%
    pull(motif)

  if (length(motifs) < 2) {
    return(invisible(NULL))
  }

  wide <- agg_df %>%
    select(read_id, motif, mean_probability) %>%
    distinct() %>%
    pivot_wider(names_from = motif, values_from = mean_probability)

  combos <- combn(motifs, 2, simplify = FALSE)
  
  plots <- walk(combos, function(pair) {
    x_col <- pair[[1]] %>%  as.character()
    print(x_col)
    y_col <- pair[[2]] %>%  as.character()
    plot_df <- wide %>%
      select(all_of(c(x_col, y_col))) %>%
      drop_na() 

    if (nrow(plot_df) == 0) {
      return()
    }

    p <- ggplot(plot_df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
      geom_point(alpha = 0.2, size = 0.8, color = "#984ea3") +
      geom_density_2d_filled(alpha = 0.3) +
      geom_xsidedensity(fill = "#984ea3", alpha = 0.3) +
      geom_ysidedensity(fill = "#984ea3", alpha = 0.3) +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1)) +
      theme_minimal() +
      labs(
        title = glue("Mean methylation relationship: {x_col} vs {y_col}"),
        x = glue("Mean probability ({x_col})"),
        y = glue("Mean probability ({y_col})")
      )

    filename <- glue("pair_{sanitize(pair[[1]])}_vs_{sanitize(pair[[2]])}.png")
    save_plot(p, filename, width = 7, height = 6)
    return(p)
  })
  return(plots)
}

sanitize <- function(label) {
  str_replace_all(label, "[^A-Za-z0-9]+", "_")
}

plots <- generate_pair_plots(aggregate_reads, limit = 4)

cat(glue("Analysis complete.\nOutputs written to {out_dir}\n"))





