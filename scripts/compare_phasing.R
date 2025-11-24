#!/usr/bin/env Rscript

suppressWarnings({
  args <- commandArgs(trailingOnly = TRUE)
})

print_usage <- function() {
  cat(
    "Usage: Rscript compare_phasing.R --output-prefix <prefix> --input label=path ...\n",
    "  Supports Floria haplosets files and split-read clustering TSVs.\n",
    "Example:\n",
    "  Rscript compare_phasing.R --output-prefix results/compare \\\n",
    "    --input floria_imputed=path/to/contig.haplosets \\\n",
    "    --input floria_raw=path/to/contig.haplosets \\\n",
    "    --input split_reads=path/to/read_clustering.tsv\n",
    sep = ""
  )
}

if (length(args) == 0) {
  print_usage()
  quit(status = 1)
}

input_files <- list()
output_prefix <- NULL
i <- 1
while (i <= length(args)) {
  token <- args[[i]]
  if (token == "--output-prefix") {
    if (i == length(args)) {
      stop("missing value after --output-prefix")
    }
    output_prefix <- args[[i + 1]]
    i <- i + 2
  } else if (token == "--input") {
    if (i == length(args)) {
      stop("missing value after --input")
    }
    spec <- args[[i + 1]]
    parts <- strsplit(spec, "=", fixed = TRUE)[[1]]
    if (length(parts) != 2) {
      stop("expected label=path for --input, got: ", spec)
    }
    label <- parts[[1]]
    path <- parts[[2]]
    input_files[[label]] <- path
    i <- i + 2
  } else if (token %in% c("-h", "--help")) {
    print_usage()
    quit(status = 0)
  } else {
    stop("unknown argument: ", token)
  }
}

maybe_prompt <- function(title, message, default = "") {
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
    value <- rstudioapi::showPrompt(title = title, message = message, default = default)
    if (is.null(value)) "" else value
  } else {
    cat(message, if (nzchar(default)) paste0(" [", default, "]") else "", ": ", sep = "")
    value <- readline()
    if (value == "") default else value
  }
}

maybe_select_file <- function(caption) {
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
    path <- rstudioapi::selectFile(caption = caption)
    if (is.null(path)) "" else path
  } else {
    cat(caption, ": ")
    path <- readline()
    path
  }
}

if ((is.null(output_prefix) || length(input_files) < 2) && interactive()) {
  message("Entering interactive input mode (provide at least two datasets).")
  if (is.null(output_prefix) || output_prefix == "") {
    output_prefix <- maybe_prompt("Output prefix", "Enter output prefix", "results/phasing_compare")
  }
  while (length(input_files) < 2) {
    label <- maybe_prompt("Dataset label", sprintf("Enter label for dataset #%d (blank to finish)", length(input_files) + 1), "")
    if (label == "") break
    path <- maybe_select_file(sprintf("Select file for dataset '%s'", label))
    if (path != "") {
      input_files[[label]] <- path
      message("Added dataset '", label, "' -> ", path)
    } else {
      message("No file selected for dataset '", label, "'.")
    }
  }
}

if (is.null(output_prefix)) {
  stop("missing --output-prefix")
}
if (length(input_files) < 2) {
  stop("provide at least two --input label=path arguments")
}

read_haploset <- function(path, label) {
  if (!file.exists(path)) {
    stop("file does not exist: ", path)
  }
  lines <- readLines(path)
  if (length(lines) == 0) {
    return(data.frame())
  }
  current_hap <- NA_character_
  rows <- list()
  for (line in lines) {
    if (line == "" || startsWith(line, "#")) {
      next
    }
    if (startsWith(line, ">")) {
      token <- substring(line, 2)
      hap_id <- strsplit(token, "\\s+")[[1]][1]
      hap_id <- sub("\\..*$", "", hap_id)
      current_hap <- hap_id
      next
    }
    parts <- strsplit(trimws(line), "\\s+")[[1]]
    if (length(parts) < 1 || is.na(current_hap)) {
      next
    }
    read_id <- parts[[1]]
    start <- if (length(parts) >= 2) as.numeric(parts[[2]]) else NA_real_
    end <- if (length(parts) >= 3) as.numeric(parts[[3]]) else NA_real_
    rows[[length(rows) + 1]] <- data.frame(
      read_id = read_id,
      hap = current_hap,
      start = start,
      end = end,
      dataset = label,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

read_split_reads <- function(path, label) {
  if (!file.exists(path)) {
    stop("file does not exist: ", path)
  }
  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  if (!("read_id" %in% colnames(df)) || !("cluster_id" %in% colnames(df))) {
    stop("split-read TSV missing required columns read_id and cluster_id: ", path)
  }
  df <- df[!is.na(df$cluster_id) & df$cluster_id != -1, , drop = FALSE]
  data.frame(
    read_id = df$read_id,
    hap = paste0("cluster_", df$cluster_id),
    start = NA_real_,
    end = NA_real_,
    dataset = label,
    stringsAsFactors = FALSE
  )
}

datasets <- lapply(names(input_files), function(label) {
  path <- input_files[[label]]
  if (grepl("\\.haplosets$", path, ignore.case = TRUE)) {
    read_haploset(path, label)
  } else if (grepl("read_clustering\\.tsv$", path, ignore.case = TRUE)) {
    read_split_reads(path, label)
  } else {
    warning("Unrecognized file type for ", path, "; attempting haploset parser.")
    read_haploset(path, label)
  }
})
names(datasets) <- names(input_files)

assignments <- lapply(names(datasets), function(label) {
  df <- datasets[[label]]
  if (nrow(df) == 0) {
    return(data.frame(read_id = character(), stringsAsFactors = FALSE))
  }
  unique_df <- df[!duplicated(df$read_id), c("read_id", "hap")]
  colnames(unique_df) <- c("read_id", label)
  unique_df
})

combined <- Reduce(function(x, y) merge(x, y, by = "read_id", all = TRUE), assignments)

combined_path <- paste0(output_prefix, "_read_assignments.tsv")
write.table(
  combined,
  file = combined_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pairwise_summary <- data.frame(
  dataset_a = character(),
  dataset_b = character(),
  matching = integer(),
  mismatching = integer(),
  only_a = integer(),
  only_b = integer(),
  shared = integer(),
  optimal_matching = double(),
  optimal_fraction = double(),
  mismatch_annotated = integer(),
  stringsAsFactors = FALSE
)

mapping_rows <- data.frame(
  dataset_a = character(),
  hap_a = character(),
  dataset_b = character(),
  hap_b = character(),
  shared_reads = integer(),
  stringsAsFactors = FALSE
)

compute_optimal_mapping <- function(contingency) {
  n_rows <- nrow(contingency)
  n_cols <- ncol(contingency)
  shared <- sum(contingency)
  if (shared == 0) {
    return(list(total = 0, mapping = data.frame()))
  }
  if (!requireNamespace("clue", quietly = TRUE)) {
    warning("Package 'clue' not installed; optimal mappings computed via identity only.")
    total <- sum(diag(contingency[, intersect(colnames(contingency), rownames(contingency)), drop = FALSE]))
    mapping <- data.frame()
    return(list(total = total, mapping = mapping))
  }
  size <- max(n_rows, n_cols)
  padded <- matrix(0, nrow = size, ncol = size)
  padded[seq_len(n_rows), seq_len(n_cols)] <- contingency
  assignment <- clue::solve_LSAP(padded, maximum = TRUE)
  total <- sum(padded[cbind(seq_len(size), assignment)])
  map_rows <- list()
  row_names <- rownames(contingency)
  col_names <- colnames(contingency)
  for (i in seq_len(n_rows)) {
    j <- assignment[i]
    if (j <= n_cols) {
      count <- contingency[i, j]
      if (count > 0) {
        map_rows[[length(map_rows) + 1]] <- data.frame(
          hap_a = row_names[i],
          hap_b = col_names[j],
          shared_reads = count,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  mapping <- if (length(map_rows) == 0) data.frame() else do.call(rbind, map_rows)
  list(total = total, mapping = mapping)
}

labels <- names(datasets)
if (length(labels) >= 2) {
  combos <- combn(labels, 2, simplify = FALSE)
  for (pair in combos) {
    a <- pair[[1]]
    b <- pair[[2]]
    df <- merge(
      datasets[[a]][, c("read_id", "hap")],
      datasets[[b]][, c("read_id", "hap")],
      by = "read_id",
      all = TRUE,
      suffixes = c("_a", "_b")
    )
    shared_mask <- !is.na(df$hap_a) & !is.na(df$hap_b)
    match_count <- sum(shared_mask & df$hap_a == df$hap_b)
  mismatch_count_raw <- sum(shared_mask & df$hap_a != df$hap_b)
  only_a <- sum(!is.na(df$hap_a) & is.na(df$hap_b))
  only_b <- sum(is.na(df$hap_a) & !is.na(df$hap_b))
    shared_reads <- match_count + mismatch_count_raw
    shared_df <- df[!is.na(df$hap_a) & !is.na(df$hap_b), , drop = FALSE]
    optimal_matching <- NA_real_
    optimal_fraction <- NA_real_
    matched_reads <- match_count
    mismatch_annotated <- mismatch_count_raw
    if (shared_reads > 0) {
      contingency <- table(shared_df$hap_a, shared_df$hap_b)
      opt <- compute_optimal_mapping(contingency)
      optimal_matching <- opt$total
      optimal_fraction <- if (shared_reads > 0) opt$total / shared_reads else NA_real_
      if (nrow(opt$mapping) > 0) {
        opt_map <- cbind(
          dataset_a = a,
          opt$mapping,
          dataset_b = b
        )
        opt_map <- opt_map[, c("dataset_a", "hap_a", "dataset_b", "hap_b", "shared_reads")]
        mapping_rows <- rbind(mapping_rows, opt_map)
      }
    }
    final_match_count <- matched_reads + ifelse(is.na(optimal_matching), 0, optimal_matching)
    mismatch_count <- max(shared_reads - final_match_count, 0)
    similarity <- if (shared_reads > 0) final_match_count / shared_reads else NA_real_
    pairwise_summary <- rbind(
      pairwise_summary,
      data.frame(
        dataset_a = a,
        dataset_b = b,
        matching = final_match_count,
        mismatching = mismatch_count,
        only_a = only_a,
        only_b = only_b,
        shared = shared_reads,
        optimal_matching = optimal_matching,
        optimal_fraction = optimal_fraction,
        mismatch_annotated = mismatch_annotated,
        similarity = similarity,
        stringsAsFactors = FALSE
      )
    )
  }
}

summary_path <- paste0(output_prefix, "_pairwise_summary.tsv")
write.table(
  pairwise_summary,
  file = summary_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

mapping_path <- NULL
if (nrow(mapping_rows) > 0) {
  mapping_path <- paste0(output_prefix, "_pairwise_optimal_mapping.tsv")
  write.table(
    mapping_rows,
    file = mapping_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

plot_dir <- dirname(output_prefix)
if (!dir.exists(plot_dir) && nzchar(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
}

plot_base <- paste0(output_prefix, "_plots")

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  warning("Package ggplot2 not installed; skipping pairwise similarity heatmap.")
} else if (nrow(pairwise_summary) > 0) {
  suppressPackageStartupMessages(library(ggplot2))
  heatmap_plot <- ggplot(pairwise_summary, aes(x = dataset_a, y = dataset_b, fill = similarity)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", similarity)), color = "black", size = 3) +
    scale_fill_gradient(low = "#f7fbff", high = "#08306b", na.value = "grey90", limits = c(0, 1)) +
    theme_minimal() +
    labs(title = "Pairwise Similarity (matching / total shared reads)",
         x = "Dataset A",
         y = "Dataset B",
         fill = "Similarity")

  heatmap_path <- paste0(plot_base, "_similarity_heatmap.png")
  ggsave(heatmap_path, heatmap_plot, width = 6, height = 5, dpi = 300)
}

if (!requireNamespace("UpSetR", quietly = TRUE)) {
  warning("Package UpSetR not installed; skipping UpSet plot.")
} else if (nrow(combined) > 0) {
  suppressPackageStartupMessages(library(UpSetR))
  long_assignments <- data.frame(
    read_id = character(),
    dataset = character(),
    hap = character(),
    stringsAsFactors = FALSE
  )
  for (label in names(input_files)) {
    df <- combined[, c("read_id", label), drop = FALSE]
    colnames(df) <- c("read_id", "hap")
    df <- df[!is.na(df$hap), , drop = FALSE]
    if (nrow(df) > 0) {
      df$dataset <- label
      long_assignments <- rbind(long_assignments, df[, c("read_id", "dataset", "hap")])
    }
  }
  if (nrow(long_assignments) > 0) {
    long_assignments$cluster_id <- paste(long_assignments$dataset, long_assignments$hap, sep = ":")
    cluster_ids <- long_assignments$cluster_id[order(-table(long_assignments$cluster_id)[long_assignments$cluster_id])]
    cluster_ids <- unique(cluster_ids)
    dataset_cols <- setNames(rainbow(length(input_files)), names(input_files))
    set_colors <- sapply(cluster_ids, function(cid) {
      dataset <- strsplit(cid, ":", fixed = TRUE)[[1]][1]
      dataset_cols[[dataset]]
    })
    names(set_colors) <- cluster_ids

    indicator <- as.data.frame(matrix(0, nrow = nrow(combined), ncol = length(cluster_ids)))
    colnames(indicator) <- cluster_ids
    indicator$read_id <- combined$read_id
    for (cid in cluster_ids) {
      reads <- long_assignments$read_id[long_assignments$cluster_id == cid]
      indicator[[cid]] <- as.integer(indicator$read_id %in% reads)
    }
    upset_data <- indicator
    upset_plot <- upset(
      upset_data,
      sets = cluster_ids,
      nsets = length(cluster_ids),
      nintersects = 20,
      keep.order = TRUE,
      order.by = "freq",
      decreasing = c(TRUE, TRUE),
      sets.bar.color = set_colors
    )
    upset_path <- paste0(plot_base, "_upset.png")
    png(upset_path, width = 1600, height = 1000, res = 150)
    print(upset_plot)
    dev.off()
  }
}

cat(
  "Comparison written to:\n  ", combined_path, "\n  ", summary_path, "\n",
  if (!is.null(mapping_path)) paste0("  ", mapping_path, "\n") else "",
  if (exists("heatmap_path")) paste0("  ", heatmap_path, "\n") else "",
  if (exists("upset_path")) paste0("  ", upset_path, "\n") else "",
  sep = ""
)
