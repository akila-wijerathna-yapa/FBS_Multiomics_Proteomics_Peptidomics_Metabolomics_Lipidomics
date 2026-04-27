#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
  library(mixOmics)
  library(MOFA2)
  library(glue)
  library(limma)
  library(ggrepel)
  library(circlize)
  library(WGCNA)
  library(igraph)
  library(gprofiler2)
  library(dynamicTreeCut)
})

# ---- User-configurable parameters -----------------------------------------
n_pca_components <- 5
n_diablo_components <- 2
n_mofa_factors <- 5
max_keepX <- 50
default_plot_height_increment <- 0
heatmap_plot_height_increment <- 2
pathway_plot_height_increment <- 2

# Colour palette requested by the user
grade_colors <- c(
  "High-grade" = "#EF3F37",
  "Medium-grade" = "#FBAF41",
  "Low-grade" = "#262161"
)

dataset_display_names <- c(
  lipidomics = "Lipidomics",
  metabolomics = "Metabolomics",
  peptidomics = "Peptidomics",
  proteomics = "Proteomics"
)

dataset_palette <- c(
  lipidomics = "#4C72B0",
  metabolomics = "#55A868",
  peptidomics = "#8172B2",
  proteomics = "#C44E52"
)

# ---- Paths ----------------------------------------------------------------
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", cmd_args[grep(file_arg, cmd_args)])
  if (length(script_path) > 0) {
    return(dirname(normalizePath(script_path)))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()
workspace_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
dataset_dir <- workspace_dir
results_root <- file.path(dirname(workspace_dir), "multiomics_results")

dir_create(results_root, recurse = TRUE)

analysis_dirs <- c(
  "PCA",
  "DIABLO",
  "MOFA2",
  "Differential",
  "Correlation",
  "Networks",
  "WGCNA",
  "Circos",
  "Pathways",
  file.path("Heatmaps", "Expression"),
  file.path("Heatmaps", "Differential"),
  file.path("Heatmaps", "Coexpression"),
  "MetaboAnalyst"
)
purrr::walk(analysis_dirs, ~dir_create(path(results_root, .x), recurse = TRUE))

# ---- Helper functions ------------------------------------------------------
standardize_sample_id <- function(x) {
  clean <- tolower(x)
  clean <- gsub("[^a-z0-9]", "_", clean)

  grade <- case_when(
    str_detect(clean, "high") ~ "High-grade",
    str_detect(clean, "medium") ~ "Medium-grade",
    str_detect(clean, "low") ~ "Low-grade",
    TRUE ~ NA_character_
  )

  replicate <- readr::parse_number(clean)
  replicate <- if_else(is.na(replicate), NA_real_, replicate)

  if (any(is.na(grade))) {
    stop("Unable to match one or more sample names to a grade: ",
         paste(x[is.na(grade)], collapse = ", "))
  }

  if (any(is.na(replicate))) {
  replicate <- stats::ave(replicate, grade, FUN = function(idx) {
      idx[is.na(idx)] <- seq_len(sum(is.na(idx)))
      idx
    })
  }

  standardized <- glue("{grade}_{as.integer(replicate)}")

  if (any(duplicated(standardized))) {
    dupes <- standardized[duplicated(standardized)]
    stop("Duplicate sample identifiers generated: ",
         paste(unique(dupes), collapse = ", "))
  }

  standardized
}

rename_sample_columns <- function(df, id_cols) {
  sample_cols <- setdiff(names(df), id_cols)
  standardized <- standardize_sample_id(sample_cols)
  colnames(df)[match(sample_cols, names(df))] <- standardized
  df
}

make_unique_ids <- function(ids) {
  if (anyDuplicated(ids)) {
    make.unique(ids, sep = "__dup")
  } else {
    ids
  }
}

shorten_labels <- function(labels, max_length = 30) {
  if (is.null(labels)) {
    return(labels)
  }

  truncated <- labels
  too_long <- !is.na(truncated) & stringr::str_length(truncated) > max_length
  truncated[too_long] <- paste0(stringr::str_sub(truncated[too_long], 1, max_length - 3), "...")

  placeholder <- "__NA_PLACEHOLDER__"
  safe_vals <- ifelse(is.na(truncated) | truncated == "", placeholder, truncated)
  safe_vals <- make.unique(safe_vals, sep = " ")
  safe_vals[safe_vals == placeholder] <- NA_character_

  safe_vals
}

shorten_metabolite_labels <- function(labels, max_words = 3, max_length = 28) {
  if (is.null(labels)) {
    return(labels)
  }

  trimmed <- stringr::str_trim(labels)
  leading_words <- stringr::word(trimmed, 1, max_words)
  leading_words <- ifelse(is.na(leading_words) | leading_words == "", trimmed, leading_words)
  shorten_labels(leading_words, max_length = max_length)
}

matrix_from_df <- function(df, feature_col) {
  df[[feature_col]] <- make_unique_ids(df[[feature_col]])
  mat <- df |> column_to_rownames(feature_col) |> as.matrix()
  mode(mat) <- "numeric"
  mat
}

drop_zero_variance_features <- function(mat) {
  if (nrow(mat) == 0) {
    return(mat)
  }
  sds <- apply(mat, 1, function(row_vals) stats::sd(row_vals, na.rm = TRUE))
  keep <- sds > 0
  mat[keep, , drop = FALSE]
}

save_plot <- function(plot_obj, file_stub, output_dir, width = 7, height = 5, height_increment = default_plot_height_increment) {
  dir_create(output_dir, recurse = TRUE)
  png_path <- file.path(output_dir, paste0(file_stub, ".png"))
  svg_path <- file.path(output_dir, paste0(file_stub, ".svg"))
  effective_height <- height + height_increment
  ggsave(filename = png_path, plot = plot_obj, width = width, height = effective_height, dpi = 300)
  ggsave(filename = svg_path, plot = plot_obj, width = width, height = effective_height, dpi = 300)
  invisible(list(png = png_path, svg = svg_path))
  svg_path <- file.path(output_dir, paste0(file_stub, ".svg"))
  effective_height <- height + height_increment
  ggsave(filename = png_path, plot = plot_obj, width = width, height = effective_height, dpi = 300)
  ggsave(filename = svg_path, plot = plot_obj, width = width, height = effective_height, dpi = 300)
  invisible(list(png = png_path, svg = svg_path))
}

save_base_plot <- function(draw_fn, file_stub, output_dir, width = 7, height = 7, height_increment = default_plot_height_increment) {
  dir_create(output_dir, recurse = TRUE)
  png_path <- file.path(output_dir, paste0(file_stub, ".png"))
  svg_path <- file.path(output_dir, paste0(file_stub, ".svg"))
  effective_height <- height + height_increment
  grDevices::png(filename = png_path, width = width, height = effective_height, units = "in", res = 300)
  draw_fn()
  grDevices::dev.off()
  grDevices::svg(filename = svg_path, width = width, height = effective_height)
  draw_fn()
  grDevices::dev.off()
  invisible(list(png = png_path, svg = svg_path))
}

geom_label_repel <- function(..., box.padding = 0.5, point.padding = 0.3, max.overlaps = Inf) {
  ggrepel::geom_text_repel(
    ...,
    box.padding = box.padding,
    point.padding = point.padding,
    max.overlaps = max.overlaps,
    min.segment.length = 0
  )
}

extract_primary_id <- function(ids) {
  ids <- ifelse(is.na(ids), "", ids)
  primary <- str_extract(ids, "^[^;,\\s]+")
  str_trim(primary)
}

annotate_protein_groups <- function(protein_groups) {
  if (length(protein_groups) == 0) {
    return(tibble())
  }

  primary_ids <- extract_primary_id(protein_groups)
  lookup_tbl <- tibble(
    protein_group = protein_groups,
    primary_id = primary_ids
  ) |> distinct(protein_group, .keep_all = TRUE)

  unique_ids <- unique(primary_ids)
  unique_ids <- unique_ids[!is.na(unique_ids) & unique_ids != ""]
  mapping <- tibble()

  if (length(unique_ids) > 0) {
    mapping <- tryCatch(
      gprofiler2::gconvert(
        query = unique_ids,
        organism = "hsapiens",
        filter_na = TRUE
      ),
      error = function(e) {
        warning(glue("Protein annotation lookup failed: {e$message}"))
        tibble()
      }
    )
  }

  if (nrow(mapping) > 0) {
    mapping <- mapping |> distinct(input, .keep_all = TRUE)
  }

  annotated <- lookup_tbl |> 
    left_join(mapping |> dplyr::select(input, name, description), by = c("primary_id" = "input")) |> 
    mutate(
      gene_symbol = name,
      protein_description = description,
      protein_label = coalesce(gene_symbol, primary_id)
    ) |> 
    dplyr::select(protein_group, primary_id, gene_symbol, protein_label, protein_description)

  annotated
}

ensure_feature_label_column <- function(feature_tbl) {
  if (is.null(feature_tbl)) {
    return(feature_tbl)
  }
  if (!"feature_label" %in% colnames(feature_tbl)) {
    feature_tbl <- feature_tbl |> mutate(feature_label = feature_id)
  } else {
    feature_tbl <- feature_tbl |> mutate(
      feature_label = if_else(is.na(feature_label) | feature_label == "", feature_id, feature_label)
    )
  }
  feature_tbl
}

augment_feature_display_column <- function(feature_tbl) {
  if (is.null(feature_tbl) || nrow(feature_tbl) == 0) {
    return(feature_tbl)
  }

  feature_tbl <- ensure_feature_label_column(feature_tbl)

  feature_tbl <- feature_tbl |> 
    mutate(
      feature_label = coalesce(feature_label, feature_id),
      feature_label = if_else(is.na(feature_label) | feature_label == "", feature_id, feature_label)
    )

  if (!"feature_display" %in% colnames(feature_tbl)) {
    feature_tbl <- feature_tbl |> mutate(feature_display = feature_label)
  }

  feature_tbl <- feature_tbl |> 
    mutate(feature_display = coalesce(feature_display, feature_label)) |> 
    mutate(feature_display = if_else(feature_display == "", feature_label, feature_display)) |> 
    mutate(base_display = feature_display) |> 
    group_by(base_display) |> 
    mutate(
      feature_display = dplyr::case_when(
        base_display == "" ~ feature_label,
        n() == 1 ~ base_display,
        TRUE ~ glue("{base_display} #{row_number()}")
      )
    ) |> 
    ungroup() |> 
    dplyr::select(-base_display)

  feature_tbl
}

strip_bracket_suffix <- function(x) {
  if (is.null(x)) {
    return(x)
  }
  stringr::str_trim(stringr::str_remove(x, "\\s*\\[[^\\]]*\\]$"))
}

remove_no_ms2_prefix <- function(x) {
  if (is.null(x)) {
    return(x)
  }
  x <- stringr::str_replace(x, "^\\s*no\\s+MS2:\\s*", "")
  x <- stringr::str_replace(x, "^\\s*w/o\\s+MS2:?\\s*", "")
  x
}

dataset_label <- function(dataset_name) {
  dataset_vec <- as.character(dataset_name)
  purrr::map_chr(dataset_vec, function(name) {
    label <- dataset_display_names[[name]]
    if (is.null(label)) {
      str_to_title(name)
    } else {
      label
    }
  })
}

safe_mixomics_plot <- function(file_stub, output_dir, plot_fn, width = 7, height = 7, height_increment = default_plot_height_increment) {
  tryCatch({
    save_base_plot(
      function() {
        suppressWarnings(plot_fn())
      },
      file_stub,
      output_dir,
      width = width,
      height = height,
      height_increment = height_increment
    )
  }, error = function(e) {
    warning(glue("Failed to create {file_stub}: {e$message}"))
    invisible(NULL)
  })
}

create_volcano_plot <- function(contrast_df, dataset_name, contrast_name, output_dir,
                                logfc_threshold = 1, adjp_threshold = 0.05,
                                up_label = NULL, down_label = NULL) {
  if (nrow(contrast_df) == 0) return(invisible(NULL))
  
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(glue)
  library(ggtext)
  library(scales)
  
  safe_contrast <- str_replace_all(tolower(contrast_name), "[^a-z0-9]+", "_")
  pretty_contrast <- str_replace_all(contrast_name, "_vs_", " vs ") |> str_replace_all("_", " ")
  
  if (!"feature_display" %in% colnames(contrast_df)) {
    contrast_df <- contrast_df |> mutate(feature_display = NA_character_)
  }
  
  contrast_df <- contrast_df |>
    mutate(
      feature_label   = coalesce(feature_label, feature_id),
      feature_display = if_else(is.na(feature_display) | feature_display == "", feature_label, feature_display)
    )
  
  up_caption   <- if (!is.null(up_label)   && !is.na(up_label))   glue("+{up_label}")   else "Up"
  down_caption <- if (!is.null(down_label) && !is.na(down_label)) glue("+{down_label}") else "Down"
  
  colour_values <- c("#EF3F37", "#262161", "#b3b3b3")
  names(colour_values) <- c(up_caption, down_caption, "Not significant")
  
  plot_df <- contrast_df |>
    mutate(
      neg_log_p = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
      regulation = dplyr::case_when(
        is.na(adj.P.Val) ~ "Not significant",
        adj.P.Val < adjp_threshold & logFC >= logfc_threshold  ~ up_caption,
        adj.P.Val < adjp_threshold & logFC <= -logfc_threshold ~ down_caption,
        TRUE ~ "Not significant"
      ),
      regulation = factor(regulation, levels = c(up_caption, down_caption, "Not significant"), ordered = TRUE)
    )
  
  top_labels <- plot_df |> arrange(adj.P.Val, dplyr::desc(abs(logFC))) |> slice_head(n = 10)
  
  volcano <- ggplot(plot_df, aes(x = logFC, y = neg_log_p, colour = regulation)) +
    geom_point(alpha = 0.7, size = 1.8) +
    scale_colour_manual(values = colour_values, drop = FALSE,
                        breaks = c(up_caption, down_caption)) +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold),
               linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(adjp_threshold),
               linetype = "dashed", colour = "grey40") +
    geom_text_repel(
      data = top_labels,
      aes(label = feature_display),
      size = 3,
      max.overlaps = 20,
      colour = "black",
      box.padding = 0.3,
      point.padding = 0.2,
      min.segment.length = 0
    ) +
    annotate("text", x = min(plot_df$logFC, na.rm = TRUE), y = -log10(adjp_threshold),
             label = "P Adj. = 0.05", hjust = 0, vjust = -0.6, size = 3) +
    labs(
      title = glue("{dataset_label(dataset_name)} volcano"),
      subtitle = pretty_contrast,
      x = "<span style='background-color:#FFEB3B;padding:2px 6px;'>Fold change log2</span>",
      y = "P value (-log10)",
      colour = NULL
    ) +
    theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
      legend.position = "right",
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = alpha("white", 0), colour = NA),
      legend.box.background = element_blank(),
      legend.key = element_rect(fill = alpha("white", 0), colour = NA),
      legend.text = element_text(size = 9),
      plot.title = element_text(face = "bold"),
      axis.title.x = element_markdown(),
      plot.margin = margin(5, 40, 5, 5)
    ) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  save_plot(volcano, glue("{dataset_name}_volcano_{safe_contrast}"), output_dir, width = 8, height = 5)
}




run_limma_differential <- function(mat, dataset_name, sample_metadata, feature_metadata, output_dir) {
  if (is.null(mat) || nrow(mat) < 3 || ncol(mat) < 3) {
    warning(glue("Skipping differential analysis for {dataset_name}: insufficient data."))
    return(NULL)
  }

  design <- model.matrix(~ 0 + group, data = sample_metadata)
  design_raw_names <- gsub("^group", "", colnames(design))
  design_sanitized <- make.names(design_raw_names)
  colnames(design) <- design_sanitized

  name_lookup <- setNames(design_sanitized, design_raw_names)
  required_groups <- c("High-grade", "Medium-grade", "Low-grade")
  missing_groups <- setdiff(required_groups, names(name_lookup))
  if (length(missing_groups) > 0) {
    warning(glue("Skipping differential analysis for {dataset_name}: missing groups {paste(missing_groups, collapse = ', ')}."))
    return(NULL)
  }

  if (qr(design)$rank < ncol(design)) {
    warning(glue("Design matrix for {dataset_name} is not full rank; skipping differential analysis."))
    return(NULL)
  }

  fit <- limma::lmFit(mat, design)

  contrast_defs <- c(
    High_vs_Low = glue("{name_lookup['High-grade']} - {name_lookup['Low-grade']}") ,
    High_vs_Medium = glue("{name_lookup['High-grade']} - {name_lookup['Medium-grade']}") ,
    Medium_vs_Low = glue("{name_lookup['Medium-grade']} - {name_lookup['Low-grade']}")
  )

  contrast_meta <- tibble(
    contrast = names(contrast_defs),
    numerator = c("High-grade", "High-grade", "Medium-grade"),
    denominator = c("Low-grade", "Medium-grade", "Low-grade")
  )

  contrast_matrix <- limma::makeContrasts(contrasts = contrast_defs, levels = design)

  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)

  results_tbl <- purrr::map_dfr(colnames(contrast_matrix), function(contrast_name) {
    limma::topTable(fit2, coef = contrast_name, number = Inf, sort.by = "P") |> 
      tibble::rownames_to_column("feature_id") |> 
      as_tibble() |> 
      mutate(contrast = contrast_name)
  })

  if (!is.null(feature_metadata)) {
    join_cols <- intersect(colnames(feature_metadata), c(
      "feature_id", "feature_label", "feature_display", "name", "name_raw", "protein_group",
      "peptide_sequence", "primary_id", "gene_symbol", "protein_label",
      "protein_description", "mz", "rt", "adduct", "formula", "mode"
    ))
    if (length(join_cols) > 1) {
  results_tbl <- results_tbl |> left_join(feature_metadata |> dplyr::select(all_of(join_cols)), by = "feature_id")
    } else {
      results_tbl <- results_tbl |> left_join(feature_metadata, by = "feature_id")
    }
  }

  if (!"feature_display" %in% colnames(results_tbl)) {
    results_tbl <- results_tbl |> mutate(feature_display = NA_character_)
  }

  results_tbl <- results_tbl |> 
    mutate(
      dataset = dataset_name,
      feature_label = coalesce(feature_label, feature_id),
      feature_display = if_else(is.na(feature_display) | feature_display == "", feature_label, feature_display)
    )

  write_tibble(results_tbl, glue("{dataset_name}_differential_results"), output_dir)

  purrr::walk(unique(results_tbl$contrast), function(contrast_name) {
    contrast_df <- results_tbl |> filter(.data$contrast == contrast_name)
    contrast_info <- contrast_meta |> filter(.data$contrast == contrast_name)
    up_group <- if (nrow(contrast_info) > 0) contrast_info$numerator else NA_character_
    down_group <- if (nrow(contrast_info) > 0) contrast_info$denominator else NA_character_
    create_volcano_plot(
      contrast_df,
      dataset_name,
      contrast_name,
      output_dir,
      up_label = up_group,
      down_label = down_group
    )
  })

  results_tbl
}

plot_sample_correlations <- function(mat, dataset_name, sample_metadata, output_dir) {
  if (is.null(mat) || ncol(mat) < 2) {
    return(invisible(NULL))
  }

  sample_labels <- sample_metadata |> 
    mutate(sample_label = if_else(is.na(sample_label) | sample_label == "", sample_id, sample_label)) |> 
    pull(sample_label)
  names(sample_labels) <- sample_metadata$sample_id
  sample_order <- sample_metadata$sample_id

  cor_mat <- stats::cor(mat, use = "pairwise.complete.obs")
  cor_df <- tibble::as_tibble(cor_mat, rownames = "sample_id_1") |> 
    pivot_longer(-sample_id_1, names_to = "sample_id_2", values_to = "correlation") |> 
    mutate(
      sample_id_1 = factor(sample_id_1, levels = sample_order),
      sample_id_2 = factor(sample_id_2, levels = sample_order),
      label_1 = sample_labels[as.character(sample_id_1)],
      label_2 = sample_labels[as.character(sample_id_2)]
    )

  correlation_plot <- ggplot(cor_df, aes(x = sample_id_1, y = sample_id_2, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-1, 1), low = "#262161", mid = "#ffffff", high = "#EF3F37", midpoint = 0) +
    labs(
      title = glue("{dataset_label(dataset_name)} sample correlations"),
      x = NULL,
      y = NULL,
      fill = "Pearson r"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = sample_labels[sample_order]) +
    scale_y_discrete(labels = sample_labels[sample_order])

  save_plot(
    correlation_plot,
    glue("{dataset_name}_sample_correlation"),
    output_dir,
    width = 7,
    height = 6,
    height_increment = heatmap_plot_height_increment
  )

  write_tibble(cor_df, glue("{dataset_name}_sample_correlations"), output_dir)
}

plot_cross_dataset_correlations <- function(assays, sample_metadata, output_dir) {
  dataset_summaries <- purrr::imap_dfr(assays, function(mat, dataset_name) {
    if (is.null(mat) || nrow(mat) == 0) {
      return(tibble())
    }

    summary_vals <- colMeans(mat, na.rm = TRUE)
    tibble(
      dataset = dataset_name,
      sample_id = names(summary_vals),
      mean_abundance = as.numeric(summary_vals)
    )
  })

  dataset_summaries <- dataset_summaries |> filter(!is.na(sample_id))
  if (nrow(dataset_summaries) == 0) {
    return(invisible(NULL))
  }

  dataset_wide <- dataset_summaries |> 
    pivot_wider(names_from = sample_id, values_from = mean_abundance)

  if (ncol(dataset_wide) <= 2) {
    return(invisible(NULL))
  }

  dataset_order <- dataset_wide$dataset
  dataset_labels <- setNames(purrr::map_chr(dataset_order, dataset_label), dataset_order)

  summary_mat <- dataset_wide |> 
    tibble::column_to_rownames("dataset") |> 
    as.matrix()

  cor_mat <- stats::cor(t(summary_mat), use = "pairwise.complete.obs")

  cor_df <- cor_mat |> 
    as_tibble(rownames = "dataset_1") |> 
    pivot_longer(-dataset_1, names_to = "dataset_2", values_to = "correlation") |> 
    mutate(
      dataset_1 = factor(dataset_1, levels = dataset_order),
      dataset_2 = factor(dataset_2, levels = dataset_order),
      label_1 = dataset_labels[as.character(dataset_1)],
      label_2 = dataset_labels[as.character(dataset_2)]
    )

  correlation_plot <- ggplot(cor_df, aes(x = dataset_1, y = dataset_2, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-1, 1), low = "#262161", mid = "#ffffff", high = "#EF3F37", midpoint = 0) +
    labs(
      title = "Cross-omics sample correlations",
      x = NULL,
      y = NULL,
      fill = "Pearson r"
    ) +
    theme_minimal() +
    scale_x_discrete(labels = dataset_labels[dataset_order]) +
    scale_y_discrete(labels = dataset_labels[dataset_order])

  save_plot(
    correlation_plot,
    "cross_omics_sample_correlation",
    output_dir,
    width = 6,
    height = 5,
    height_increment = heatmap_plot_height_increment
  )

  write_tibble(cor_df, "cross_omics_sample_correlations", output_dir)
}

plot_expression_heatmap <- function(mat, dataset_name, sample_metadata, feature_metadata, output_dir, top_n = 50) {
  if (is.null(mat) || nrow(mat) == 0) {
    return(invisible(NULL))
  }

  variances <- apply(mat, 1, stats::var, na.rm = TRUE)
  variances <- variances[is.finite(variances)]
  if (length(variances) == 0) {
    return(invisible(NULL))
  }

  selected <- names(sort(variances, decreasing = TRUE))[seq_len(min(top_n, length(variances)))]
  sample_order <- sample_metadata$sample_id
  selected_mat <- mat[selected, sample_order, drop = FALSE]
  scaled_mat <- t(scale(t(selected_mat)))
  scaled_mat[is.na(scaled_mat)] <- 0

  label_lookup <- NULL
  if (!is.null(feature_metadata) && "feature_display" %in% colnames(feature_metadata)) {
    label_lookup <- setNames(feature_metadata$feature_display, feature_metadata$feature_id)
  }

  feature_order <- selected
  feature_display_order <- purrr::map_chr(feature_order, function(fid) {
    if (is.null(label_lookup)) {
      return(fid)
    }
    val <- label_lookup[[fid]]
    if (is.null(val) || is.na(val) || val == "") fid else val
  })

  sample_labels <- setNames(sample_metadata$sample_label, sample_metadata$sample_id)
  missing_samples <- is.na(sample_labels) | sample_labels == ""
  sample_labels[missing_samples] <- names(sample_labels)[missing_samples]

  heatmap_df <- as_tibble(scaled_mat, rownames = "feature_id") |> 
    pivot_longer(-feature_id, names_to = "sample_id", values_to = "zscore") |> 
    mutate(
      feature_display = purrr::map_chr(feature_id, function(fid) {
        if (is.null(label_lookup)) {
          return(fid)
        }
        val <- label_lookup[[fid]]
        if (is.null(val) || is.na(val) || val == "") fid else val
      }),
      feature_display = factor(feature_display, levels = rev(unique(feature_display_order))),
      sample_id = factor(sample_id, levels = sample_order),
      sample_label = sample_labels[as.character(sample_id)]
    )

  plot_height <- max(6, 4 + 0.15 * length(feature_order))

  expression_plot <- ggplot(heatmap_df, aes(x = sample_id, y = feature_display, fill = zscore)) +
    geom_tile() +
    scale_fill_gradient2(low = "#262161", mid = "#ffffff", high = "#EF3F37", midpoint = 0) +
    labs(
      title = glue("{dataset_label(dataset_name)} expression heatmap"),
      x = "Sample",
      y = "Feature",
      fill = "Z-score"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = sample_labels[sample_order])

  save_plot(
    expression_plot,
    glue("{dataset_name}_expression_heatmap"),
    output_dir,
    width = 8,
    height = plot_height,
    height_increment = heatmap_plot_height_increment
  )

  write_tibble(heatmap_df |> mutate(feature_display = as.character(feature_display), sample_id = as.character(sample_id)), glue("{dataset_name}_expression_heatmap"), output_dir)
}

plot_differential_heatmap <- function(diff_tbl, dataset_name, output_dir, top_n = 50) {
  if (is.null(diff_tbl) || nrow(diff_tbl) == 0) {
    return(invisible(NULL))
  }

  ranked_features <- diff_tbl |> 
    filter(!is.na(adj.P.Val)) |> 
    group_by(feature_id) |> 
    summarise(best_adj_p = min(adj.P.Val, na.rm = TRUE), .groups = "drop") |> 
    arrange(best_adj_p) |> 
    pull(feature_id) |> 
    unique() |> 
    head(top_n)

  if (length(ranked_features) == 0) {
    return(invisible(NULL))
  }

  label_lookup <- diff_tbl |> 
    mutate(feature_display = coalesce(feature_display, feature_label, feature_id)) |> 
    distinct(feature_id, feature_display) |> 
    tibble::deframe()

  feature_display_levels <- purrr::map_chr(ranked_features, function(fid) {
    val <- label_lookup[[fid]]
    if (is.null(val) || is.na(val) || val == "") fid else val
  })

  contrast_levels <- unique(diff_tbl$contrast)

  plot_df <- diff_tbl |> 
    filter(feature_id %in% ranked_features) |> 
    mutate(
      feature_display = coalesce(feature_display, feature_label, feature_id),
      feature_display = if_else(is.na(feature_display) | feature_display == "", feature_id, feature_display),
      feature_display = factor(feature_display, levels = rev(unique(feature_display_levels))),
      contrast = factor(contrast, levels = contrast_levels)
    )

  plot_height <- max(6, 4 + 0.2 * length(ranked_features))

  diff_plot <- ggplot(plot_df, aes(x = contrast, y = feature_display, fill = logFC)) +
    geom_tile() +
    scale_fill_gradient2(low = "#262161", mid = "#ffffff", high = "#EF3F37", midpoint = 0) +
    labs(
      title = glue("{dataset_label(dataset_name)} differential heatmap"),
      x = "Contrast",
      y = "Feature",
      fill = "log2 FC"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_plot(
    diff_plot,
    glue("{dataset_name}_differential_heatmap"),
    output_dir,
    width = 8,
    height = plot_height,
    height_increment = heatmap_plot_height_increment
  )

  write_tibble(
    plot_df |> mutate(
      feature_display = as.character(feature_display),
      contrast = as.character(contrast)
    ),
    glue("{dataset_name}_differential_heatmap"),
    output_dir
  )
}

plot_coexpression_heatmap <- function(mat, dataset_name, output_dir, feature_metadata = NULL, top_n = 60) {
  if (is.null(mat) || nrow(mat) < 3) {
    return(invisible(NULL))
  }

  variances <- apply(mat, 1, stats::var, na.rm = TRUE)
  variances <- variances[is.finite(variances)]
  if (length(variances) == 0) {
    return(invisible(NULL))
  }

  selected <- names(sort(variances, decreasing = TRUE))[seq_len(min(top_n, length(variances)))]
  subset_mat <- mat[selected, , drop = FALSE]

  cor_mat <- stats::cor(t(subset_mat), use = "pairwise.complete.obs")
  if (!is.matrix(cor_mat) || nrow(cor_mat) < 2) {
    return(invisible(NULL))
  }

  label_lookup <- NULL
  if (!is.null(feature_metadata) && "feature_display" %in% colnames(feature_metadata)) {
    label_lookup <- setNames(feature_metadata$feature_display, feature_metadata$feature_id)
  }

  feature_levels <- rownames(cor_mat)
  feature_display_levels <- purrr::map_chr(feature_levels, function(fid) {
    if (is.null(label_lookup)) {
      return(fid)
    }
    val <- label_lookup[[fid]]
    if (is.null(val) || is.na(val) || val == "") fid else val
  })
  feature_display_levels <- unique(feature_display_levels)

  cor_df <- cor_mat |> 
    as_tibble(rownames = "feature_id_1") |> 
    pivot_longer(-feature_id_1, names_to = "feature_id_2", values_to = "correlation") |> 
    mutate(
      feature_display_1 = purrr::map_chr(feature_id_1, function(fid) {
        val <- if (is.null(label_lookup)) NULL else label_lookup[[fid]]
        if (is.null(val) || is.na(val) || val == "") fid else val
      }),
      feature_display_2 = purrr::map_chr(feature_id_2, function(fid) {
        val <- if (is.null(label_lookup)) NULL else label_lookup[[fid]]
        if (is.null(val) || is.na(val) || val == "") fid else val
      }),
      feature_display_1 = factor(feature_display_1, levels = feature_display_levels),
      feature_display_2 = factor(feature_display_2, levels = feature_display_levels)
    )

  plot_height <- max(6, 4 + 0.15 * length(feature_levels))

  coexpr_plot <- ggplot(cor_df, aes(x = feature_display_1, y = feature_display_2, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-1, 1), low = "#262161", mid = "#ffffff", high = "#EF3F37", midpoint = 0) +
    labs(
      title = glue("{dataset_label(dataset_name)} coexpression heatmap"),
      x = "Feature",
      y = "Feature",
      fill = "Pearson r"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_plot(
    coexpr_plot,
    glue("{dataset_name}_coexpression_heatmap"),
    output_dir,
    width = 7,
    height = plot_height,
    height_increment = heatmap_plot_height_increment
  )

  write_tibble(cor_df |> mutate(
    feature_display_1 = as.character(feature_display_1),
    feature_display_2 = as.character(feature_display_2)
  ), glue("{dataset_name}_coexpression_heatmap"), output_dir)
}

analyze_multiomics_coexpression <- function(assays, feature_metadata_list, sample_metadata,
                                            heatmap_output_dir, network_output_dir,
                                            top_n_per_dataset = 40, cor_threshold = 0.7,
                                            max_edges = 200, min_edges = 30, fallback_quantile = 0.85) {
  if (length(assays) < 2) {
    return(invisible(NULL))
  }

  sample_order <- sample_metadata$sample_id

  dataset_selections <- purrr::imap(assays, function(mat, dataset_name) {
    if (is.null(mat) || nrow(mat) < 3) {
      return(NULL)
    }

    variances <- apply(mat, 1, stats::var, na.rm = TRUE)
    variances <- variances[is.finite(variances)]
    if (length(variances) == 0) {
      return(NULL)
    }

    selected_ids <- names(sort(variances, decreasing = TRUE))[seq_len(min(top_n_per_dataset, length(variances)))]
    if (length(selected_ids) < 2) {
      return(NULL)
    }

    subset_mat <- mat[selected_ids, sample_order, drop = FALSE]

    feature_meta <- feature_metadata_list[[dataset_name]]
    display_lookup <- NULL
    if (!is.null(feature_meta) && "feature_display" %in% colnames(feature_meta)) {
      display_lookup <- setNames(feature_meta$feature_display, feature_meta$feature_id)
    }

    dataset_label_value <- dataset_label(dataset_name)

    info_tbl <- tibble(
      dataset = dataset_name,
      dataset_label = dataset_label_value,
      feature_id = selected_ids,
      feature_row = paste(dataset_name, selected_ids, sep = "::"),
      feature_display = purrr::map_chr(selected_ids, function(fid) {
        if (is.null(display_lookup)) {
          return(fid)
        }
        val <- display_lookup[[fid]]
        if (is.null(val) || is.na(val) || val == "") fid else val
      }),
      feature_axis_label = glue("{dataset_label_value} | {feature_display}")
    )

    rownames(subset_mat) <- info_tbl$feature_row
    list(mat = subset_mat, info = info_tbl)
  }) |> purrr::compact()

  if (length(dataset_selections) < 2) {
    warning("Multi-omics coexpression skipped: insufficient data across datasets.")
    return(invisible(NULL))
  }

  combined_mat <- purrr::map(dataset_selections, "mat") |> do.call(what = rbind)
  if (nrow(combined_mat) < 4) {
    warning("Multi-omics coexpression skipped: not enough combined features.")
    return(invisible(NULL))
  }

  combined_info <- purrr::map_dfr(dataset_selections, "info")

  cor_mat <- stats::cor(t(combined_mat), use = "pairwise.complete.obs")
  if (!is.matrix(cor_mat) || nrow(cor_mat) < 2) {
    warning("Multi-omics coexpression skipped: correlation matrix too small.")
    return(invisible(NULL))
  }

  feature_levels <- combined_info$feature_row
  feature_lookup <- combined_info |> dplyr::select(feature_row, dataset, dataset_label, feature_display, feature_axis_label)

  cor_df <- cor_mat |> 
    as_tibble(rownames = "feature_row_1") |> 
    pivot_longer(-feature_row_1, names_to = "feature_row_2", values_to = "correlation") |> 
    left_join(feature_lookup, by = c("feature_row_1" = "feature_row")) |> 
    rename(
      dataset_1 = dataset,
      dataset_label_1 = dataset_label,
      feature_display_1 = feature_display,
      feature_axis_label_1 = feature_axis_label
    ) |> 
    left_join(feature_lookup, by = c("feature_row_2" = "feature_row")) |> 
    rename(
      dataset_2 = dataset,
      dataset_label_2 = dataset_label,
      feature_display_2 = feature_display,
      feature_axis_label_2 = feature_axis_label
    )

  cor_df <- cor_df |> 
    mutate(
      feature_axis_label_1 = coalesce(feature_axis_label_1, feature_row_1),
      feature_axis_label_2 = coalesce(feature_axis_label_2, feature_row_2),
      feature_axis_label_1 = factor(feature_axis_label_1, levels = feature_lookup$feature_axis_label[match(feature_levels, feature_lookup$feature_row)]),
      feature_axis_label_2 = factor(feature_axis_label_2, levels = feature_lookup$feature_axis_label[match(feature_levels, feature_lookup$feature_row)])
    )

  heatmap_plot <- ggplot(cor_df, aes(x = feature_axis_label_1, y = feature_axis_label_2, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-1, 1), low = "#262161", mid = "#ffffff", high = "#EF3F37", midpoint = 0) +
    labs(
      title = "Multi-omics coexpression heatmap",
      x = "Feature",
      y = "Feature",
      fill = "Pearson r"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_plot(
    heatmap_plot,
    "multiomics_coexpression_heatmap",
    heatmap_output_dir,
    width = 20,
    height = max(8, 4 + 0.1 * nrow(combined_mat)),
    height_increment = heatmap_plot_height_increment
  )

  write_tibble(
    cor_df |> mutate(
      feature_display_1 = as.character(feature_display_1),
      feature_display_2 = as.character(feature_display_2),
      feature_axis_label_1 = as.character(feature_axis_label_1),
      feature_axis_label_2 = as.character(feature_axis_label_2)
    ),
    "multiomics_coexpression_correlations",
    heatmap_output_dir
  )

  edges_ranked <- cor_df |> 
    filter(feature_row_1 < feature_row_2) |> 
    filter(dataset_1 != dataset_2) |> 
    filter(!is.na(correlation)) |> 
    mutate(abs_corr = abs(correlation)) |> 
    arrange(desc(abs_corr))

  if (nrow(edges_ranked) == 0) {
    warning("No multi-omics coexpression pairs available after filtering.")
    return(invisible(NULL))
  }

  dynamic_threshold <- cor_threshold
  edges_df <- edges_ranked |> 
    filter(abs_corr >= dynamic_threshold)

  if (nrow(edges_df) < min_edges) {
    alt_threshold <- stats::quantile(edges_ranked$abs_corr, probs = fallback_quantile, na.rm = TRUE)
    if (is.finite(alt_threshold) && alt_threshold < dynamic_threshold) {
      dynamic_threshold <- max(alt_threshold, 0.4)
      edges_df <- edges_ranked |> filter(abs_corr >= dynamic_threshold)
    }
  }

  if (nrow(edges_df) < min_edges) {
    edges_df <- edges_ranked |> slice_head(n = min(max_edges, max(min_edges, nrow(edges_ranked))))
  } else {
    edges_df <- edges_df |> slice_head(n = max_edges)
  }

  if (nrow(edges_df) == 0) {
    warning("No multi-omics coexpression edges retained for network plotting.")
    return(invisible(NULL))
  }

  message(glue::glue(
    "Multi-omics coexpression network: using |r| >= {format(round(dynamic_threshold, 3), nsmall = 3)} to plot {nrow(edges_df)} edges (of {nrow(edges_ranked)} candidates)."
  ))

  edges_export <- edges_df |> 
    transmute(
      dataset_1,
      feature_id_1 = str_remove(feature_row_1, "^[^:]+::"),
      feature_display_1 = as.character(feature_display_1),
      feature_axis_label_1 = feature_axis_label_1,
      dataset_2,
      feature_id_2 = str_remove(feature_row_2, "^[^:]+::"),
      feature_display_2 = as.character(feature_display_2),
      feature_axis_label_2 = feature_axis_label_2,
      correlation
    )

  write_tibble(edges_export, "multiomics_coexpression_edges", network_output_dir)
  write_tibble(
    tibble(
      applied_threshold = dynamic_threshold,
      total_candidate_edges = nrow(edges_ranked),
      plotted_edges = nrow(edges_df)
    ),
    "multiomics_coexpression_edge_summary",
    network_output_dir
  )

  g <- igraph::graph_from_data_frame(edges_export, directed = FALSE)

  if (igraph::gorder(g) < 2) {
    warning("Multi-omics coexpression network skipped: insufficient nodes after filtering.")
    return(invisible(NULL))
  }

  vertex_df <- combined_info |> 
    mutate(
      vertex = feature_row,
      feature_display = coalesce(feature_display, feature_id),
      feature_axis_label = feature_axis_label,
      color = dataset_palette[dataset]
    )

  igraph::V(g)$color <- vertex_df$color[match(igraph::V(g)$name, vertex_df$vertex)]
  igraph::V(g)$dataset <- vertex_df$dataset[match(igraph::V(g)$name, vertex_df$vertex)]
  igraph::V(g)$label <- vertex_df$feature_display[match(igraph::V(g)$name, vertex_df$vertex)]

  vertex_degree <- igraph::degree(g)
  label_threshold <- sort(unique(vertex_degree), decreasing = TRUE)
  if (length(label_threshold) >= 15) {
    min_degree_for_label <- label_threshold[pmin(15, length(label_threshold))]
    show_label <- vertex_degree >= min_degree_for_label
  } else {
    show_label <- vertex_degree > 0
  }

  igraph::V(g)$label <- ifelse(show_label, igraph::V(g)$label, "")
  igraph::V(g)$size <- scales::rescale(vertex_degree, to = c(6, 18))
  igraph::V(g)$size[is.na(igraph::V(g)$size)] <- 6
  igraph::E(g)$weight <- edges_export$correlation[match(seq_len(igraph::ecount(g)), seq_len(igraph::ecount(g)))]
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, grade_colors[["High-grade"]], grade_colors[["Low-grade"]])
  edge_weight_abs <- abs(igraph::E(g)$weight)
  edge_weight_abs[!is.finite(edge_weight_abs)] <- 0

  save_base_plot(
    function() {
      layout <- igraph::layout_with_fr(g, weights = edge_weight_abs)
      igraph::plot.igraph(
        g,
        layout = layout,
        vertex.color = igraph::V(g)$color,
        vertex.label = igraph::V(g)$label,
        vertex.label.cex = 0.7,
        vertex.label.color = "black",
        vertex.size = igraph::V(g)$size,
        edge.width = scales::rescale(edge_weight_abs, to = c(0.5, 3)),
        edge.color = igraph::E(g)$color
      )
      legend(
        "topright",
        legend = names(dataset_palette),
        col = dataset_palette,
        pch = 16,
        bty = "n",
        cex = 0.8,
        title = "Dataset"
      )
      legend(
        "bottomleft",
        legend = c("Positive correlation", "Negative correlation"),
        col = c(grade_colors[["High-grade"]], grade_colors[["Low-grade"]]),
        lwd = 2,
        bty = "n",
        cex = 0.8,
        title = "Edge sign"
      )
    },
    "multiomics_coexpression_network",
    network_output_dir,
    width = 9,
    height = 7
  )

  invisible(list(correlations = cor_df, edges = edges_export))
}

prepare_metaboanalyst_exports <- function(assays, feature_metadata_list, sample_metadata, output_dir) {
  if (length(assays) == 0) {
    return(invisible(NULL))
  }

  dir_create(output_dir, recurse = TRUE)

  sample_export <- sample_metadata |> 
    mutate(
      sample_label = if_else(is.na(sample_label) | sample_label == "", sample_id, sample_label)
    )
  write_tibble(sample_export, "metaboanalyst_sample_metadata", output_dir)

  expression_long <- purrr::imap_dfr(assays, function(mat, dataset_name) {
    if (is.null(mat) || nrow(mat) == 0) {
      return(NULL)
    }

    feature_meta <- feature_metadata_list[[dataset_name]]
    display_lookup <- NULL
    if (!is.null(feature_meta) && "feature_display" %in% colnames(feature_meta)) {
      display_lookup <- setNames(feature_meta$feature_display, feature_meta$feature_id)
    }

    tibble(
      dataset = dataset_name,
      dataset_label = dataset_label(dataset_name),
      feature_id = rownames(mat),
      feature_display = purrr::map_chr(rownames(mat), function(fid) {
        if (is.null(display_lookup)) {
          return(fid)
        }
        val <- display_lookup[[fid]]
        if (is.null(val) || is.na(val) || val == "") fid else val
      })
    ) |> 
      bind_cols(as_tibble(mat)) |> 
      pivot_longer(cols = -c(dataset, dataset_label, feature_id, feature_display), names_to = "sample_id", values_to = "intensity")
  }) |> 
    mutate(sample_id = factor(sample_id, levels = sample_metadata$sample_id)) |> 
    arrange(dataset, feature_id, sample_id)

  if (nrow(expression_long) == 0) {
    return(invisible(NULL))
  }

  write_tibble(expression_long, "metaboanalyst_expression_long", output_dir)

  purrr::iwalk(assays, function(mat, dataset_name) {
    if (is.null(mat) || nrow(mat) == 0) {
      return(NULL)
    }

    feature_meta <- feature_metadata_list[[dataset_name]]
    display_lookup <- NULL
    if (!is.null(feature_meta) && "feature_display" %in% colnames(feature_meta)) {
      display_lookup <- setNames(feature_meta$feature_display, feature_meta$feature_id)
    }

    export_df <- as_tibble(mat, rownames = "feature_id") |> 
      mutate(
        feature_display = purrr::map_chr(feature_id, function(fid) {
          if (is.null(display_lookup)) {
            return(fid)
          }
          val <- display_lookup[[fid]]
          if (is.null(val) || is.na(val) || val == "") fid else val
        })
      ) |> 
      relocate(feature_display, .after = feature_id)

    write_tibble(export_df, glue("{dataset_name}_metaboanalyst_matrix"), output_dir)
  })

  combined_matrix <- expression_long |> 
    mutate(feature_key = paste(dataset, feature_id, sep = "::")) |> 
    dplyr::select(feature_key, sample_id, intensity) |> 
    pivot_wider(names_from = sample_id, values_from = intensity) |> 
    arrange(feature_key)

  write_tibble(combined_matrix, "multiomics_metaboanalyst_matrix", output_dir)

  feature_annotations <- expression_long |> 
    distinct(dataset, dataset_label, feature_id, feature_display) |> 
    arrange(dataset, feature_id)

  write_tibble(feature_annotations, "multiomics_metaboanalyst_feature_annotations", output_dir)

  invisible(list(
    expression_long = expression_long,
    sample_metadata = sample_export,
    combined_matrix = combined_matrix,
    feature_annotations = feature_annotations
  ))
}

impute_matrix_column_means <- function(mat) {
  if (!is.matrix(mat) || ncol(mat) == 0) {
    return(mat)
  }

  for (j in seq_len(ncol(mat))) {
    col_vals <- mat[, j]
    if (all(is.na(col_vals))) {
      mat[, j] <- 0
    } else {
      mean_val <- mean(col_vals, na.rm = TRUE)
      mat[is.na(col_vals), j] <- mean_val
    }
  }

  mat
}

plot_metaboanalyst_overview <- function(expression_long, sample_metadata, output_dir, top_n_heatmap = 40) {
  if (is.null(expression_long) || nrow(expression_long) == 0) {
    return(invisible(NULL))
  }

  dir_create(output_dir, recurse = TRUE)

  sample_info <- sample_metadata |> 
    mutate(sample_label = if_else(is.na(sample_label) | sample_label == "", sample_id, sample_label))

  sample_levels <- sample_info$sample_id
  sample_label_map <- setNames(sample_info$sample_label, sample_info$sample_id)
  sample_label_levels <- unique(unname(sample_label_map[sample_levels]))

  expression_augmented <- expression_long |> 
    mutate(
      feature_axis_label = glue("{dataset_label} | {feature_display}"),
      sample_id = factor(sample_id, levels = sample_levels)
    )

  expression_wide <- expression_augmented |> 
    dplyr::select(sample_id, feature_axis_label, intensity) |> 
    pivot_wider(names_from = feature_axis_label, values_from = intensity) |> 
    arrange(sample_id)

  if (ncol(expression_wide) <= 1) {
    return(invisible(NULL))
  }

  expression_matrix <- expression_wide |> 
    dplyr::select(-sample_id) |> 
    as.matrix()

  rownames(expression_matrix) <- as.character(expression_wide$sample_id)

  if (ncol(expression_matrix) == 0 || nrow(expression_matrix) < 3) {
    return(invisible(NULL))
  }

  col_vars <- apply(expression_matrix, 2, stats::var, na.rm = TRUE)
  keep_cols <- names(col_vars[is.finite(col_vars) & col_vars > 0])

  if (length(keep_cols) < 2) {
    return(invisible(NULL))
  }

  expression_matrix <- expression_matrix[, keep_cols, drop = FALSE]
  col_vars <- col_vars[keep_cols]

  expression_matrix <- impute_matrix_column_means(expression_matrix)
  scaled_matrix <- scale(expression_matrix, center = TRUE, scale = TRUE)
  scaled_matrix[!is.finite(scaled_matrix)] <- 0

  pca_res <- tryCatch({
    stats::prcomp(scaled_matrix, center = FALSE, scale. = FALSE)
  }, error = function(e) {
    NULL
  })

  if (!is.null(pca_res) && ncol(pca_res$x) >= 2) {
    pca_scores <- pca_res$x[, 1:2, drop = FALSE]
    colnames(pca_scores) <- c("PC1", "PC2")
    pca_df <- as_tibble(pca_scores, rownames = "sample_id")
  } else if (!is.null(pca_res) && ncol(pca_res$x) == 1) {
    scores <- pca_res$x
    pca_df <- tibble(
      sample_id = rownames(scores),
      PC1 = scores[, 1],
      PC2 = 0
    )
  } else {
    pca_df <- tibble()
  }

  if (nrow(pca_df) > 0) {
    pca_df <- pca_df |> 
      left_join(sample_info |> dplyr::select(sample_id, group, sample_label), by = "sample_id")

    pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = group)) +
      geom_point(size = 3) +
      geom_label_repel(aes(label = sample_label), size = 3, max.overlaps = 20, show.legend = FALSE) +
      labs(
        title = "MetaboAnalyst-style PCA",
        x = "PC1",
        y = "PC2",
        colour = "Group"
      ) +
      theme_minimal()

    group_levels <- unique(stats::na.omit(pca_df$group))
    if (length(group_levels) > 0 && all(group_levels %in% names(grade_colors))) {
      pca_plot <- pca_plot + scale_colour_manual(values = grade_colors, drop = FALSE)
    }

    save_plot(pca_plot, "metaboanalyst_pca", output_dir, width = 7, height = 5)

    write_tibble(pca_df, "metaboanalyst_pca_scores", output_dir)
  }

  top_n <- min(top_n_heatmap, length(col_vars))
  top_features <- names(sort(col_vars, decreasing = TRUE))[seq_len(top_n)]

  heatmap_matrix <- scaled_matrix[, top_features, drop = FALSE]

  if (ncol(heatmap_matrix) == 0) {
    return(invisible(NULL))
  }

  feature_levels <- rev(top_features)

  heatmap_df <- as_tibble(heatmap_matrix, rownames = "sample_id") |> 
    pivot_longer(-sample_id, names_to = "feature_axis_label", values_to = "zscore") |> 
    mutate(
  sample_label = sample_label_map[as.character(sample_id)],
      sample_label = factor(sample_label, levels = sample_label_levels),
      feature_axis_label = as.character(feature_axis_label)
    )

  feature_lookup <- expression_augmented |> 
    distinct(feature_axis_label, dataset, dataset_label, feature_display) |> 
    mutate(feature_axis_label = as.character(feature_axis_label))

  heatmap_df <- heatmap_df |> 
    left_join(feature_lookup, by = "feature_axis_label") |> 
    mutate(feature_axis_label = factor(feature_axis_label, levels = feature_levels))

  heatmap_plot <- ggplot(heatmap_df, aes(x = sample_label, y = feature_axis_label, fill = zscore)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-3, 3), low = "#262161", mid = "#ffffff", high = "#EF3F37", midpoint = 0, oob = scales::squish) +
    labs(
      title = "MetaboAnalyst-style heatmap",
      x = "Sample",
      y = "Feature",
      fill = "Z-score"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks = element_blank()
    )

  save_plot(
    heatmap_plot,
    "metaboanalyst_heatmap",
    output_dir,
    width = 9,
    height = max(6, 4 + 0.18 * length(feature_levels)),
    height_increment = heatmap_plot_height_increment
  )

  write_tibble(
    heatmap_df |> mutate(
      sample_label = as.character(sample_label),
      feature_axis_label = as.character(feature_axis_label)
    ),
    "metaboanalyst_heatmap_zscores",
    output_dir
  )

  invisible(list(pca = pca_df, heatmap = heatmap_df))
}

build_coexpression_network <- function(mat, dataset_name, output_dir, feature_metadata = NULL, top_n = 100, cor_threshold = 0.8) {
  if (is.null(mat) || nrow(mat) < 3 || ncol(mat) < 3) {
    warning(glue("Skipping coexpression network for {dataset_name}: insufficient data."))
    return(invisible(NULL))
  }

  variances <- apply(mat, 1, stats::var, na.rm = TRUE)
  variances <- variances[is.finite(variances)]
  if (length(variances) == 0) {
    warning(glue("Unable to compute variances for {dataset_name}; skipping coexpression network."))
    return(invisible(NULL))
  }

  top_features <- names(sort(variances, decreasing = TRUE))[seq_len(min(top_n, length(variances)))]
  mat_subset <- mat[top_features, , drop = FALSE]

  cor_mat <- stats::cor(t(mat_subset), use = "pairwise.complete.obs")
  cor_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  colnames(cor_df) <- c("feature_id_1", "feature_id_2", "correlation")

  edges_df <- cor_df |> 
    filter(feature_id_1 < feature_id_2) |> 
    filter(!is.na(correlation)) |> 
    filter(abs(correlation) >= cor_threshold)

  label_lookup <- NULL
  if (!is.null(feature_metadata)) {
    label_lookup <- feature_metadata |> 
      dplyr::select(feature_id, feature_label, feature_display) |> 
      distinct()
  }

  if (!is.null(label_lookup)) {
    edges_df <- edges_df |> 
      left_join(label_lookup, by = c("feature_id_1" = "feature_id")) |> 
      rename(feature_label_1 = feature_label, feature_display_1 = feature_display) |> 
      left_join(label_lookup, by = c("feature_id_2" = "feature_id")) |> 
      rename(feature_label_2 = feature_label, feature_display_2 = feature_display) |> 
      mutate(
        feature_display_1 = coalesce(feature_display_1, feature_label_1, feature_id_1),
        feature_display_2 = coalesce(feature_display_2, feature_label_2, feature_id_2)
      )
  }

  write_tibble(edges_df, glue("{dataset_name}_coexpression_edges"), output_dir)

  if (nrow(edges_df) == 0) {
    warning(glue("No edges passed the correlation threshold for {dataset_name}."))
    return(invisible(NULL))
  }

  g <- igraph::graph_from_data_frame(edges_df, directed = FALSE)

  node_color <- dataset_palette[[dataset_name]]
  if (is.null(node_color)) {
    node_color <- "#4C72B0"
  }

  vertex_labels <- igraph::V(g)$name
  if (!is.null(label_lookup)) {
    label_map <- setNames(coalesce(label_lookup$feature_display, label_lookup$feature_label, label_lookup$feature_id), label_lookup$feature_id)
    vertex_labels <- label_map[igraph::V(g)$name]
    vertex_labels[is.na(vertex_labels)] <- igraph::V(g)$name[is.na(vertex_labels)]
  }

  vertex_degree <- igraph::degree(g)
  label_threshold <- sort(unique(vertex_degree), decreasing = TRUE)
  if (length(label_threshold) >= 10) {
    min_degree_for_label <- label_threshold[pmin(10, length(label_threshold))]
    show_label <- vertex_degree >= min_degree_for_label
  } else {
    show_label <- vertex_degree > 0
  }

  igraph::V(g)$label <- ifelse(show_label, vertex_labels, "")
  igraph::V(g)$label.cex <- 0.7
  igraph::V(g)$label.color <- "black"
  igraph::V(g)$size <- scales::rescale(vertex_degree, to = c(6, 18))
  igraph::V(g)$size[is.na(igraph::V(g)$size)] <- 6
  igraph::V(g)$color <- node_color

  igraph::E(g)$color <- ifelse(igraph::E(g)$correlation > 0, grade_colors[["High-grade"]], grade_colors[["Low-grade"]])

  save_base_plot(
    function() {
      edge_weights <- scales::rescale(abs(igraph::E(g)$correlation), to = c(0.5, 3))
      igraph::plot.igraph(
        g,
        vertex.label = igraph::V(g)$label,
        vertex.label.cex = igraph::V(g)$label.cex,
        vertex.label.color = igraph::V(g)$label.color,
        vertex.size = igraph::V(g)$size,
        vertex.color = igraph::V(g)$color,
        edge.width = edge_weights,
        layout = igraph::layout_with_fr(g),
        edge.color = igraph::E(g)$color
      )
      legend(
        "topright",
        legend = c("Positive", "Negative"),
        col = c(grade_colors[["High-grade"]], grade_colors[["Low-grade"]]),
        lwd = 2,
        bty = "n",
        cex = 0.8
      )
    },
    glue("{dataset_name}_coexpression_network"),
    output_dir,
    width = 7,
    height = 7
  )
}

plot_loading_heatmap <- function(loadings_df, dataset_name, output_dir, top_features_per_component = 15) {
  if (nrow(loadings_df) == 0) {
    return(invisible(NULL))
  }

  if (!"feature_label" %in% colnames(loadings_df)) {
    loadings_df <- loadings_df |> mutate(feature_label = feature_id)
  }

  if (!"feature_display" %in% colnames(loadings_df)) {
    loadings_df <- loadings_df |> mutate(feature_display = feature_label)
  }

  component_cols <- grep("^(PC|Comp)", colnames(loadings_df), value = TRUE)
  if (length(component_cols) == 0) {
    component_cols <- setdiff(colnames(loadings_df), c("feature_id", "feature_label", "feature_display"))
  }

  loadings_long <- loadings_df |> 
    dplyr::select(feature_id, feature_label, feature_display, all_of(component_cols)) |> 
    pivot_longer(cols = all_of(component_cols), names_to = "component", values_to = "loading")

  top_loadings <- loadings_long |> 
    group_by(component) |> 
    arrange(dplyr::desc(abs(loading))) |> 
    mutate(rank = dplyr::row_number()) |> 
    filter(rank <= top_features_per_component) |> 
    mutate(feature_plot = factor(feature_display, levels = rev(unique(feature_display)))) |> 
    ungroup() |> 
    dplyr::select(-rank)

  loading_plot <- ggplot(top_loadings, aes(x = component, y = feature_plot, fill = loading)) +
    geom_tile() +
    scale_fill_gradient2(low = "#262161", mid = "#ffffff", high = "#EF3F37", midpoint = 0) +
    labs(
      title = glue("{dataset_label(dataset_name)} PCA loadings"),
      x = "Principal component",
      y = "Feature",
      fill = "Loading"
    ) +
  theme_minimal()

  save_plot(
    loading_plot,
    glue("{dataset_name}_pca_loading_heatmap"),
    output_dir,
    width = 7,
    height = 6,
    height_increment = heatmap_plot_height_increment
  )
}

run_wgcna_analysis <- function(mat, dataset_name, output_dir, sample_metadata, feature_metadata = NULL, top_features = 1000) {
  if (is.null(mat) || nrow(mat) < 50 || ncol(mat) < 4) {
    warning(glue("Skipping WGCNA for {dataset_name}: requires at least 50 features and 4 samples."))
    return(invisible(NULL))
  }

  variances <- apply(mat, 1, stats::var, na.rm = TRUE)
  variances <- variances[order(variances, decreasing = TRUE)]
  selected_features <- names(variances)[seq_len(min(top_features, length(variances)))]
  mat_subset <- mat[selected_features, , drop = FALSE]
  dat_expr <- t(mat_subset)

  WGCNA::allowWGCNAThreads()
  gsg <- WGCNA::goodSamplesGenes(dat_expr, verbose = 3)
  if (!gsg$allOK) {
    dat_expr <- dat_expr[gsg$goodSamples, gsg$goodGenes]
  }

  powers <- c(1:10, seq(12, 20, 2))
  sft <- WGCNA::pickSoftThreshold(dat_expr, powerVector = powers, verbose = 5)
  soft_power <- sft$powerEstimate
  if (is.na(soft_power)) {
    soft_power <- 6
  }

  adjacency <- WGCNA::adjacency(dat_expr, power = soft_power)
  tom <- WGCNA::TOMsimilarity(adjacency)
  diss_tom <- 1 - tom
  gene_tree <- hclust(as.dist(diss_tom), method = "average")

  module_labels <- dynamicTreeCut::cutreeDynamic(
    dendro = gene_tree,
    distM = diss_tom,
    deepSplit = 2,
    pamRespectsDendro = FALSE,
    minClusterSize = 20
  )
  module_colors <- WGCNA::labels2colors(module_labels)

  save_base_plot(
    function() {
      WGCNA::plotDendroAndColors(
        gene_tree,
        module_colors,
        "Module colors",
        dendroLabels = FALSE,
        hang = 0.03,
        main = glue("{dataset_label(dataset_name)} WGCNA")
      )
    },
    glue("{dataset_name}_wgcna_dendrogram"),
    output_dir,
    width = 7,
    height = 5
  )

  module_assignments <- tibble(
    feature_id = colnames(dat_expr),
    module = module_colors
  )

  if (!is.null(feature_metadata)) {
    cols_to_keep <- setdiff(colnames(feature_metadata), "dataset")
    module_assignments <- module_assignments |> 
      left_join(feature_metadata |> dplyr::select(any_of(cols_to_keep)), by = "feature_id")
  }

  if ("feature_label" %in% colnames(module_assignments)) {
    module_assignments <- module_assignments |> mutate(feature_label = coalesce(feature_label, feature_id))
  } else {
    module_assignments <- module_assignments |> mutate(feature_label = feature_id)
  }

  if (!"feature_display" %in% colnames(module_assignments)) {
    module_assignments <- module_assignments |> mutate(feature_display = feature_label)
  }

  module_assignments <- module_assignments |> 
    mutate(feature_display = if_else(is.na(feature_display) | feature_display == "", feature_label, feature_display))

  module_assignments <- module_assignments |> mutate(dataset = dataset_name)

  module_eigengenes <- WGCNA::moduleEigengenes(dat_expr, colors = module_colors)$eigengenes |> 
    tibble::rownames_to_column("sample_id") |> 
    as_tibble()

  module_eigengenes <- module_eigengenes |> 
  left_join(sample_metadata |> dplyr::select(sample_id, any_of(c("group", "color"))), by = "sample_id") |> 
    mutate(dataset = dataset_name)

  write_tibble(module_assignments, glue("{dataset_name}_wgcna_modules"), output_dir)
  write_tibble(module_eigengenes, glue("{dataset_name}_wgcna_eigengenes"), output_dir)

  invisible(list(modules = module_assignments, eigengenes = module_eigengenes))
}

draw_diablo_circos <- function(diablo_model, output_dir, cutoff = 0.6) {
  n_components <- ncol(diablo_model$variates[[1]])
  if (n_components == 0) {
    warning("No DIABLO components available for circos plotting.")
    return(invisible(NULL))
  }

  block_order <- names(diablo_model$X)
  block_colors <- dataset_palette[block_order]
  if (length(block_colors) != length(block_order)) {
    block_colors <- rep("#4C72B0", length(block_order))
  }
  missing_blocks <- is.na(block_colors)
  if (any(missing_blocks)) {
    block_colors[missing_blocks] <- "#4C72B0"
  }

  y_levels <- levels(diablo_model$Y)
  y_colors <- grade_colors[y_levels]
  if (length(y_colors) != length(y_levels)) {
    y_colors <- rep("#555555", length(y_levels))
  }
  missing_y <- is.na(y_colors)
  if (any(missing_y)) {
    y_colors[missing_y] <- "#555555"
  }

  label_lookup <- purrr::imap(diablo_model$loadings, function(loading_mat, block_name) {
    original_labels <- rownames(loading_mat)
    shorten_labels(original_labels, max_length = 32)
  })

  comp_seq <- seq_len(n_components)
  purrr::walk(comp_seq, function(comp_idx) {
    safe_mixomics_plot(
      file_stub = glue("diablo_circos_comp{comp_idx}"),
      output_dir = output_dir,
      width = 15,
      height = 15,
      plot_fn = function() {
        circlize::circos.clear() # Ensure enlarged tracks/labels start from a clean slate
        circlize::circos.par(track.margin = c(0.01, 0.18), cell.padding = c(0, 0, 0, 0))
        mixOmics::circosPlot(
          diablo_model,
          comp = comp_idx,
          cutoff = cutoff,
          line = TRUE,
          size.variables = 1.85,
          size.labels = 3,
          size.legend = 1.65,
          linkWidth = c(1.0, 6.5),
          color.blocks = block_colors,
          color.Y = y_colors,
          legend = TRUE,
          var.names = label_lookup,
          trackHeight = 0.28
        )
      }
    )
  })
}

plot_diablo_component_correlations <- function(diablo_scores_wide, output_dir) {
  if (nrow(diablo_scores_wide) == 0) {
    return(invisible(NULL))
  }

  comp_cols <- grep("^comp", colnames(diablo_scores_wide), value = TRUE)
  if (length(comp_cols) == 0) {
    return(invisible(NULL))
  }

  dataset_levels <- unique(diablo_scores_wide$dataset)
  dataset_labels <- setNames(purrr::map_chr(dataset_levels, dataset_label), dataset_levels)

  purrr::walk(comp_cols, function(comp_col) {
    comp_df <- diablo_scores_wide |> 
      dplyr::select(sample_id, dataset, all_of(comp_col)) |> 
      pivot_wider(names_from = dataset, values_from = all_of(comp_col))

    numeric_cols <- setdiff(colnames(comp_df), "sample_id")
    numeric_cols <- intersect(numeric_cols, dataset_levels)
    if (length(numeric_cols) < 2) {
      return()
    }

    comp_mat <- comp_df |> dplyr::select(all_of(numeric_cols)) |> as.matrix()
    colnames(comp_mat) <- numeric_cols
    rownames(comp_mat) <- comp_df$sample_id
    cor_mat <- stats::cor(comp_mat, use = "pairwise.complete.obs")
    cor_tbl <- tidyr::expand_grid(
      dataset_1 = numeric_cols,
      dataset_2 = numeric_cols
    ) |> 
      mutate(correlation = as.vector(cor_mat)) |> 
      mutate(
        dataset_1 = factor(dataset_1, levels = dataset_levels),
        dataset_2 = factor(dataset_2, levels = dataset_levels),
        label_1 = dataset_labels[as.character(dataset_1)],
        label_2 = dataset_labels[as.character(dataset_2)]
      )

    heatmap_plot <- ggplot(cor_tbl, aes(x = dataset_1, y = dataset_2, fill = correlation)) +
      geom_tile() +
      scale_fill_gradient2(limits = c(-1, 1), low = "#262161", mid = "#ffffff", high = "#EF3F37", midpoint = 0) +
      labs(
        title = glue("DIABLO cross-omics correlations ({comp_col})"),
        x = "Dataset",
        y = "Dataset",
        fill = "Pearson r"
      ) +
      theme_minimal()

    heatmap_plot <- heatmap_plot +
      scale_x_discrete(labels = dataset_labels[dataset_levels]) +
      scale_y_discrete(labels = dataset_labels[dataset_levels])

    safe_stub <- str_replace_all(comp_col, "[^a-z0-9]+", "_")
    save_plot(
      heatmap_plot,
      glue("diablo_cross_omics_{safe_stub}"),
      output_dir,
      width = 6,
      height = 5,
      height_increment = heatmap_plot_height_increment
    )
  })

  all_scores <- diablo_scores_wide |> 
    pivot_longer(cols = starts_with("comp"), names_to = "component", values_to = "score") |> 
    mutate(
      component_index = readr::parse_number(component),
      dataset_component = glue("{dataset_label(dataset)} • Comp {component_index}"),
      dataset_component_id = glue("{dataset}_{component}")
    ) |> 
    dplyr::select(sample_id, dataset_component_id, dataset_component, score)

  dataset_component_lookup <- all_scores |> 
    dplyr::select(dataset_component_id, dataset_component) |> 
    distinct()

  all_scores_wide <- all_scores |> 
    dplyr::select(sample_id, dataset_component_id, score) |> 
    pivot_wider(names_from = dataset_component_id, values_from = score)

  numeric_cols <- setdiff(colnames(all_scores_wide), "sample_id")
  if (length(numeric_cols) < 2) {
    return(invisible(NULL))
  }

  all_mat <- all_scores_wide |> dplyr::select(all_of(numeric_cols)) |> as.matrix()
  colnames(all_mat) <- numeric_cols
  rownames(all_mat) <- all_scores_wide$sample_id
  cor_mat_all <- stats::cor(all_mat, use = "pairwise.complete.obs")
  cor_tbl_all <- tidyr::expand_grid(
    feature_1 = numeric_cols,
    feature_2 = numeric_cols
  ) |> 
    mutate(correlation = as.vector(cor_mat_all))

  name_lookup <- setNames(dataset_component_lookup$dataset_component, dataset_component_lookup$dataset_component_id)

  feature_levels <- sort(unique(cor_tbl_all$feature_1))
  cor_tbl_all <- cor_tbl_all |> 
    mutate(
      feature_1 = factor(feature_1, levels = feature_levels),
      feature_2 = factor(feature_2, levels = feature_levels),
      label_1 = name_lookup[as.character(feature_1)],
      label_2 = name_lookup[as.character(feature_2)]
    )

  heatmap_all <- ggplot(cor_tbl_all, aes(x = feature_1, y = feature_2, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-1, 1), low = "#262161", mid = "#ffffff", high = "#EF3F37", midpoint = 0) +
    labs(
      title = "DIABLO cross-omics correlations (all components)",
      x = "Dataset & component",
      y = "Dataset & component",
      fill = "Pearson r"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(labels = name_lookup[feature_levels]) +
    scale_y_discrete(labels = name_lookup[feature_levels])

  save_plot(
    heatmap_all,
    "diablo_cross_omics_all_components",
    output_dir,
    width = 8,
    height = 6,
    height_increment = heatmap_plot_height_increment
  )
}

plot_diablo_loadings <- function(diablo_model, output_dir, top_n = 20, feature_lookup = NULL) {
  loadings_list <- diablo_model$loadings
  if (length(loadings_list) == 0) {
    return(invisible(NULL))
  }

  purrr::iwalk(loadings_list, function(loading_mat, dataset_name) {
    if (!is.null(feature_lookup) && !(dataset_name %in% names(feature_lookup))) {
      return()
    }

    if (is.null(loading_mat) || nrow(loading_mat) == 0) {
      return()
    }

    mapping_df <- NULL
    if (!is.null(feature_lookup) && !is.null(feature_lookup[[dataset_name]])) {
      mapping_df <- feature_lookup[[dataset_name]]
    }

    loading_tbl <- as_tibble(loading_mat, rownames = "feature_display") |> 
      pivot_longer(-feature_display, names_to = "component", values_to = "loading") |> 
      mutate(component = readr::parse_number(component))

    if (!is.null(mapping_df) && nrow(mapping_df) > 0) {
      loading_tbl <- loading_tbl |> left_join(mapping_df, by = "feature_display")
    }

    if (!"feature_id" %in% colnames(loading_tbl)) {
      loading_tbl <- loading_tbl |> mutate(feature_id = feature_display)
    }

    if (!"feature_label" %in% colnames(loading_tbl)) {
      loading_tbl <- loading_tbl |> mutate(feature_label = feature_display)
    }

    loading_tbl <- loading_tbl |> 
      mutate(feature_display = if_else(is.na(feature_display) | feature_display == "", feature_label, feature_display))

    purrr::walk(unique(loading_tbl$component), function(comp_idx) {
      plot_df <- loading_tbl |> 
        filter(component == comp_idx) |> 
        arrange(dplyr::desc(abs(loading))) |> 
        head(top_n) |> 
        mutate(feature_display = factor(feature_display, levels = rev(unique(feature_display))))

      if (nrow(plot_df) == 0) {
        return()
      }

      loading_plot <- ggplot(plot_df, aes(x = feature_display, y = loading, fill = loading > 0)) +
        geom_col(show.legend = FALSE) +
        coord_flip() +
        scale_fill_manual(values = c("TRUE" = "#EF3F37", "FALSE" = "#262161")) +
        labs(
          title = glue("{dataset_label(dataset_name)} loadings (Component {comp_idx})"),
          x = "Feature",
          y = "Loading"
        ) +
        theme_minimal()

      save_plot(
        loading_plot,
        glue("diablo_loadings_{dataset_name}_comp{comp_idx}"),
        output_dir,
        width = 7,
        height = 5
      )
    })
  })
}

run_diablo_performance <- function(diablo_model, output_dir, nrepeat = 25) {
  if (is.null(diablo_model)) {
    return(invisible(NULL))
  }

  class_counts <- table(diablo_model$Y)
  if (length(class_counts) == 0) {
    return(invisible(NULL))
  }

  max_folds <- suppressWarnings(min(5, max(2, min(class_counts))))
  if (!is.finite(max_folds) || max_folds < 2) {
    warning("Unable to run DIABLO performance evaluation: insufficient samples per class.")
    return(invisible(NULL))
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- NULL
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
  }
  set.seed(123)

  perf_res <- tryCatch(
    mixOmics::perf(
      diablo_model,
      validation = "Mfold",
      folds = max_folds,
      nrepeat = nrepeat,
      dist = "centroids.dist",
      progressBar = FALSE
    ),
    error = function(e) {
      warning(glue("DIABLO performance evaluation failed: {e$message}"))
      NULL
    }
  )

  if (had_seed && !is.null(old_seed)) {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  } else if (!had_seed && exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  if (is.null(perf_res)) {
    return(invisible(NULL))
  }

  saveRDS(perf_res, file.path(output_dir, "diablo_performance.rds"))

  extract_error_metrics <- function(err_obj, metric_name) {
    if (is.null(err_obj)) {
      return(tibble())
    }

    err_mat <- tryCatch(as.matrix(err_obj), error = function(...) NULL)
    if (is.null(err_mat)) {
      return(tibble())
    }

    if (nrow(err_mat) == 0) {
      return(tibble())
    }

    tibble(
      component = seq_len(nrow(err_mat)),
      mean_error = rowMeans(err_mat, na.rm = TRUE),
      sd_error = apply(err_mat, 1, stats::sd, na.rm = TRUE),
      metric = metric_name
    )
  }

  error_tables <- purrr::imap_dfr(perf_res$error.rate, extract_error_metrics)
  if (nrow(error_tables) > 0) {
    write_tibble(error_tables, "diablo_performance_error_rates", output_dir)
  }

  if (!is.null(perf_res$AUC)) {
    auc_tbl <- perf_res$AUC
    auc_df <- tryCatch(
      {
        auc_mat <- as.matrix(auc_tbl)
        tibble(
          component = seq_len(nrow(auc_mat)),
          mean_auc = rowMeans(auc_mat, na.rm = TRUE),
          sd_auc = apply(auc_mat, 1, stats::sd, na.rm = TRUE)
        )
      },
      error = function(...) tibble()
    )
    if (nrow(auc_df) > 0) {
      write_tibble(auc_df, "diablo_performance_auc", output_dir)
    }
  }

  safe_mixomics_plot(
    file_stub = "diablo_performance_plot",
    output_dir = output_dir,
    width = 7,
    height = 5,
    plot_fn = function() {
      mixOmics::plot(perf_res, overlay = TRUE, sd = TRUE)
    }
  )

  invisible(perf_res)
}

plot_diablo_cim <- function(diablo_model, output_dir) {
  if (is.null(diablo_model)) {
    return(invisible(NULL))
  }

  n_components <- ncol(diablo_model$variates[[1]])
  if (is.null(n_components) || n_components == 0) {
    return(invisible(NULL))
  }

  purrr::walk(seq_len(n_components), function(comp_idx) {
    safe_mixomics_plot(
      file_stub = glue("diablo_cim_comp{comp_idx}"),
      output_dir = output_dir,
      width = 8,
      height = 8,
      plot_fn = function() {
        mixOmics::cim(diablo_model, comp = comp_idx)
      }
    )
  })
}

plot_diablo_networks <- function(diablo_model, output_dir, threshold = 0.6) {
  if (is.null(diablo_model)) {
    return(invisible(NULL))
  }

  n_components <- ncol(diablo_model$variates[[1]])
  if (is.null(n_components) || n_components == 0) {
    return(invisible(NULL))
  }

  block_order <- names(diablo_model$X)
  node_colors <- dataset_palette[block_order]
  if (all(is.na(node_colors))) {
    node_colors <- rep("#4C72B0", length(block_order))
  } else {
    node_colors[is.na(node_colors)] <- "#4C72B0"
  }
  node_colors <- unname(node_colors)

  purrr::walk(seq_len(n_components), function(comp_idx) {
    safe_mixomics_plot(
      file_stub = glue("diablo_network_comp{comp_idx}"),
      output_dir = output_dir,
      width = 7,
      height = 7,
      plot_fn = function() {
        mixOmics::network(
          diablo_model,
          comp = comp_idx,
          threshold = threshold,
          color.node = node_colors,
          vertex.label.cex = 0.7
        )
      }
    )
  })
}

plot_mofa_variance_overview <- function(variance_total_df, output_dir) {
  if (is.null(variance_total_df) || nrow(variance_total_df) == 0) {
    return(invisible(NULL))
  }

  plot_df <- variance_total_df |> 
    mutate(
      dataset_label = dataset_label(dataset),
      dataset_factor = factor(dataset, levels = names(dataset_palette)),
      group = if_else(is.na(group) | group == "", "All samples", group)
    )

  variance_plot <- ggplot(plot_df, aes(x = dataset_label, y = variance_explained, fill = dataset_factor)) +
    geom_col(show.legend = FALSE) +
    facet_wrap(~ group, scales = "free_y") +
    scale_fill_manual(values = dataset_palette, guide = "none") +
    labs(
      title = "MOFA2 variance explained by view",
      x = "Dataset",
      y = "Variance explained"
    ) +
    theme_minimal()

  save_plot(
    variance_plot,
    "mofa_variance_overview",
    output_dir,
    width = 7,
    height = 5
  )
}

plot_mofa_variance_heatmap <- function(variance_per_factor_df, output_dir) {
  if (is.null(variance_per_factor_df) || nrow(variance_per_factor_df) == 0) {
    return(invisible(NULL))
  }

  plot_df <- variance_per_factor_df |> 
    mutate(
      dataset_label = dataset_label(dataset),
      group = if_else(is.na(group) | group == "", "All samples", group),
      factor = factor(factor, levels = sort(unique(factor)))
    )

  heatmap_plot <- ggplot(plot_df, aes(x = dataset_label, y = factor, fill = variance_explained)) +
    geom_tile() +
    scale_fill_gradient(low = "#ffffff", high = "#EF3F37") +
    facet_wrap(~ group) +
    labs(
      title = "MOFA2 variance explained per factor",
      x = "Dataset",
      y = "Factor",
      fill = "Variance"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_plot(
    heatmap_plot,
    "mofa_variance_per_factor",
    output_dir,
    width = 8,
    height = 6,
    height_increment = heatmap_plot_height_increment
  )
}

plot_mofa_factor_boxplots <- function(factors_df, sample_metadata, output_dir) {
  if (is.null(factors_df) || nrow(factors_df) == 0) {
    return(invisible(NULL))
  }

  factor_cols <- grep("^Factor", colnames(factors_df), value = TRUE)
  if (length(factor_cols) == 0) {
    return(invisible(NULL))
  }

  plot_df <- factors_df |> 
    dplyr::select(sample_id, all_of(factor_cols)) |> 
    pivot_longer(-sample_id, names_to = "factor", values_to = "score") |> 
    mutate(factor = factor(factor, levels = sort(unique(factor)))) |> 
    left_join(sample_metadata, by = "sample_id")

  if (!"group" %in% colnames(plot_df)) {
    return(invisible(NULL))
  }

  plot_df <- plot_df |> mutate(group = factor(group, levels = names(grade_colors)))

  box_plot <- ggplot(plot_df, aes(x = group, y = score, colour = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
    facet_wrap(~ factor, scales = "free_y") +
    scale_colour_manual(values = grade_colors, drop = FALSE) +
    labs(
      title = "MOFA2 factor scores by group",
      x = "Group",
      y = "Factor value",
      colour = "Group"
    ) +
    theme_minimal()

  save_plot(
    box_plot,
    "mofa_factor_boxplots",
    output_dir,
    width = 8,
    height = 6
  )
}

plot_mofa_factor_heatmap <- function(factors_df, sample_metadata, output_dir) {
  if (is.null(factors_df) || nrow(factors_df) == 0) {
    return(invisible(NULL))
  }

  factor_cols <- grep("^Factor", colnames(factors_df), value = TRUE)
  if (length(factor_cols) == 0) {
    return(invisible(NULL))
  }

  sample_order <- sample_metadata |> 
    arrange(factor(group, levels = names(grade_colors)), replicate) |> 
    pull(sample_id)

  sample_labels <- setNames(sample_metadata$sample_label, sample_metadata$sample_id)
  missing_labels <- is.na(sample_labels) | sample_labels == ""
  sample_labels[missing_labels] <- names(sample_labels)[missing_labels]

  heatmap_df <- factors_df |> 
    dplyr::select(sample_id, all_of(factor_cols)) |> 
    pivot_longer(-sample_id, names_to = "factor", values_to = "score") |> 
    mutate(
      factor = factor(factor, levels = sort(unique(factor))),
      sample_id = factor(sample_id, levels = sample_order)
    ) |> 
    left_join(sample_metadata |> dplyr::select(sample_id, group), by = "sample_id")

  if (nrow(heatmap_df) == 0) {
    return(invisible(NULL))
  }

  heatmap_plot <- ggplot(heatmap_df, aes(x = sample_id, y = factor, fill = score)) +
    geom_tile() +
    scale_fill_gradient2(low = "#262161", mid = "#ffffff", high = "#EF3F37", midpoint = 0) +
    scale_x_discrete(labels = sample_labels[sample_order]) +
    labs(
      title = "MOFA2 factor heatmap",
      x = "Sample",
      y = "Factor",
      fill = "Score"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  plot_height <- max(4, 3 + 0.35 * length(factor_cols))

  save_plot(
    heatmap_plot,
    "mofa_factor_heatmap",
    output_dir,
    width = 8,
    height = plot_height,
    height_increment = heatmap_plot_height_increment
  )
}

plot_mofa_top_weights <- function(weights_df, feature_metadata_list, output_dir, top_n = 20) {
  if (is.null(weights_df) || nrow(weights_df) == 0) {
    return(invisible(NULL))
  }

  feature_lookup <- purrr::imap_dfr(feature_metadata_list, function(meta_df, dataset_name) {
    if (is.null(meta_df) || nrow(meta_df) == 0) {
      return(tibble())
    }

    meta_df |> 
      dplyr::select(feature_id, any_of(c("feature_label", "feature_display"))) |> 
      mutate(dataset = dataset_name)
  })

  weights_prepared <- weights_df |> 
    mutate(
      factor_index = readr::parse_number(factor),
      factor_label = glue("Factor {factor_index}"),
      dataset_label = dataset_label(dataset)
    ) |> 
    left_join(feature_lookup, by = c("feature_id", "dataset")) |> 
    mutate(
      feature_label = coalesce(feature_label, feature_id),
      feature_display = coalesce(feature_display, feature_label, feature_id)
    )

  top_weights <- weights_prepared |> 
    mutate(abs_weight = abs(weight)) |> 
    group_by(dataset, factor_index) |> 
    slice_max(order_by = abs_weight, n = top_n, with_ties = FALSE) |> 
    ungroup()

  if (nrow(top_weights) == 0) {
    return(invisible(NULL))
  }

  write_tibble(top_weights, "mofa_top_weights", output_dir)

  top_weights |> 
    group_by(dataset, factor_index) |> 
    group_walk(function(df, key) {
      if (nrow(df) == 0) {
        return()
      }

      dataset_name <- key$dataset[[1]]
      factor_idx <- key$factor_index[[1]]

      plot_df <- df |> 
        arrange(abs_weight) |> 
        mutate(feature_display = forcats::fct_rev(forcats::fct_inorder(feature_display)))

      weight_plot <- ggplot(plot_df, aes(x = feature_display, y = weight, fill = weight > 0)) +
        geom_col(show.legend = FALSE) +
        coord_flip() +
        scale_fill_manual(values = c("TRUE" = "#EF3F37", "FALSE" = "#262161")) +
        labs(
          title = glue("{dataset_label(dataset_name)} top loadings (Factor {factor_idx})"),
          x = "Feature",
          y = "Weight"
        ) +
        theme_minimal()

      file_stub <- glue("mofa_top_weights_{dataset_name}_factor{factor_idx}")
      save_plot(weight_plot, file_stub, output_dir, width = 7, height = 5)
    })
}

summarise_mofa_factor_covariates <- function(factors_df, sample_metadata, output_dir) {
  if (is.null(factors_df) || nrow(factors_df) == 0) {
    return(invisible(NULL))
  }

  factor_cols <- grep("^Factor", colnames(factors_df), value = TRUE)
  if (length(factor_cols) == 0) {
    return(invisible(NULL))
  }

  long_df <- factors_df |> 
    dplyr::select(sample_id, all_of(factor_cols)) |> 
    pivot_longer(-sample_id, names_to = "factor", values_to = "score") |> 
    left_join(sample_metadata, by = "sample_id")

  if (!"group" %in% colnames(long_df)) {
    return(invisible(NULL))
  }

  group_stats <- long_df |> 
    group_by(factor, group) |> 
    summarise(
      mean_score = mean(score, na.rm = TRUE),
      sd_score = stats::sd(score, na.rm = TRUE),
      .groups = "drop"
    )

  write_tibble(group_stats, "mofa_factor_group_stats", output_dir)

  anova_stats <- long_df |> 
    group_by(factor) |> 
    group_map(function(df, key) {
      if (length(unique(df$group)) < 2) {
        return(tibble(factor = key$factor, p_value = NA_real_, eta_squared = NA_real_))
      }

      model <- tryCatch(stats::aov(score ~ group, data = df), error = function(...) NULL)
      if (is.null(model)) {
        return(tibble(factor = key$factor, p_value = NA_real_, eta_squared = NA_real_))
      }

      smry <- summary(model)[[1]]
      p_value <- suppressWarnings(smry[["Pr(>F)"]][1])
      ss_between <- smry[["Sum Sq"]][1]
      ss_total <- sum(smry[["Sum Sq"]], na.rm = TRUE)
      eta_sq <- ifelse(is.finite(ss_between) && is.finite(ss_total) && ss_total > 0, ss_between / ss_total, NA_real_)

      tibble(factor = key$factor, p_value = p_value, eta_squared = eta_sq)
    }) |> 
    dplyr::bind_rows()

  replicate_stats <- long_df |> 
    group_by(factor) |> 
    summarise(
      replicate_cor = tryCatch(stats::cor(score, replicate, use = "pairwise.complete.obs"), error = function(...) NA_real_),
      .groups = "drop"
    )

  covariate_stats <- anova_stats |> left_join(replicate_stats, by = "factor")
  write_tibble(covariate_stats, "mofa_factor_covariate_stats", output_dir)
}

run_pathway_enrichment <- function(diff_tbl, feature_metadata, dataset_name, output_dir, top_n = 150,
                                   id_candidates = c(
                                     "protein_group", "gene_symbol", "primary_id", "protein_label",
                                     "feature_label", "feature_display", "name", "name_raw", "feature_id"
                                   )) {
  if (is.null(diff_tbl) || nrow(diff_tbl) == 0) {
    warning(glue("Skipping pathway analysis for {dataset_name}: no differential results."))
    return(invisible(NULL))
  }

  if (is.null(feature_metadata) || !"feature_id" %in% colnames(feature_metadata)) {
    warning(glue("Skipping pathway analysis for {dataset_name}: missing feature metadata."))
    return(invisible(NULL))
  }

  available_ids <- id_candidates[id_candidates %in% colnames(feature_metadata)]
  if (length(available_ids) == 0) {
    warning(glue("Skipping pathway analysis for {dataset_name}: no identifier columns found."))
    return(invisible(NULL))
  }

  identifier_tbl <- feature_metadata |> 
    dplyr::select(feature_id, any_of(available_ids)) |> 
    pivot_longer(-feature_id, names_to = "id_type", values_to = "pathway_id") |> 
  mutate(pathway_id = stringr::str_replace_all(pathway_id, "\\|", ";")) |> 
    tidyr::separate_rows(pathway_id, sep = "[;,]\\s*") |> 
    mutate(
      pathway_id = extract_primary_id(pathway_id),
      pathway_id = strip_bracket_suffix(pathway_id),
      pathway_id = remove_no_ms2_prefix(pathway_id),
      pathway_id = str_trim(pathway_id)
    ) |> 
    filter(!is.na(pathway_id), pathway_id != "") |> 
    distinct()

  if (nrow(identifier_tbl) == 0) {
    warning(glue("Skipping pathway analysis for {dataset_name}: identifier columns empty after cleaning."))
    return(invisible(NULL))
  }

  purrr::walk(unique(diff_tbl$contrast), function(contrast_name) {
    contrast_stub <- str_replace_all(tolower(contrast_name), "[^a-z0-9]+", "_")

    top_features <- diff_tbl |> 
      filter(.data$contrast == contrast_name) |> 
      arrange(adj.P.Val) |> 
      pull(feature_id) |> 
      unique() |> 
      head(top_n)

    if (length(top_features) == 0) {
      warning(glue("No features available for pathway analysis in {dataset_name} ({contrast_name})."))
      return()
    }

    gene_list <- identifier_tbl |> 
      filter(feature_id %in% top_features) |> 
      pull(pathway_id) |> 
      unique()

    if (length(gene_list) < 3) {
      warning(glue("Not enough identifiers for pathway analysis in {dataset_name} ({contrast_name})."))
      return()
    }

    enrich_res <- tryCatch(
      gprofiler2::gost(
        query = gene_list,
        organism = "hsapiens",
        significant = TRUE,
        correction_method = "fdr"
      ),
      error = function(e) {
        warning(glue("Pathway analysis failed for {dataset_name} ({contrast_name}): {e$message}"))
        NULL
      }
    )

    if (is.null(enrich_res) || is.null(enrich_res$result) || nrow(enrich_res$result) == 0) {
      warning(glue("No enriched pathways detected for {dataset_name} ({contrast_name})."))
      return()
    }

    enrich_tbl <- as_tibble(enrich_res$result) |> 
      mutate(
        dataset = dataset_name,
        contrast = contrast_name,
        adj_p = pmax(p_value, .Machine$double.xmin),
        neg_log10_p = -log10(adj_p),
        term_name_short = str_trunc(term_name, 70)
      )

    write_tibble(enrich_tbl, glue("{dataset_name}_pathway_{contrast_stub}"), output_dir)

    ontology_sets <- list(
      GO = c("GO:BP", "GO:CC", "GO:MF"),
      KEGG = "KEGG",
      Reactome = "REAC",
      Hallmark = "HALLMARK"
    )

    contrast_title <- str_replace_all(contrast_name, "_", " ")

    purrr::iwalk(ontology_sets, function(sources, ontology_label) {
      subset_tbl <- enrich_tbl |> filter(source %in% sources)
      if (nrow(subset_tbl) == 0) {
        return()
      }

      subset_tbl <- subset_tbl |> arrange(p_value)
      write_tibble(subset_tbl, glue("{dataset_name}_{contrast_stub}_{str_to_lower(ontology_label)}_pathways"), output_dir)

      plot_df <- subset_tbl |> 
        slice_head(n = 20) |> 
        mutate(
          term_label = str_trunc(term_name, 60),
          term_label = factor(term_label, levels = rev(unique(term_label)))
        )

      plot_height <- max(4, 4 + 0.25 * nrow(plot_df))

      lollipop <- ggplot(plot_df, aes(x = neg_log10_p, y = term_label)) +
        geom_segment(aes(x = 0, xend = neg_log10_p, y = term_label, yend = term_label), colour = "#b3b3b3") +
        geom_point(aes(size = intersection_size, colour = neg_log10_p)) +
        scale_colour_gradient(low = "#262161", high = "#EF3F37") +
        scale_size_continuous(range = c(2, 6)) +
        labs(
          title = glue("{dataset_label(dataset_name)} {contrast_title} {ontology_label} pathways"),
          x = "-log10 adjusted p-value",
          y = NULL,
          size = "Overlap",
          colour = "-log10 adj. p"
        ) +
        theme_minimal()

      save_plot(
        lollipop,
        glue("{dataset_name}_{contrast_stub}_{str_to_lower(ontology_label)}_pathways"),
        output_dir,
        width = 8,
        height = plot_height,
        height_increment = pathway_plot_height_increment
      )
    })
  })
}
write_tibble <- function(tb, file_stub, output_dir) {
  dir_create(output_dir, recurse = TRUE)
  out_path <- file.path(output_dir, paste0(file_stub, ".csv"))
  readr::write_csv(tb, out_path)
  out_path
}

# ---- Data ingestion --------------------------------------------------------
sample_name_mappings <- list()

lipidomics_raw <- read_csv(
  file.path(dataset_dir, "lipiomics", "Lipidom_Xf_log2CSS.filtered.cleaned.csv"),
  show_col_types = FALSE
)

lipidomics_sample_cols_raw <- setdiff(colnames(lipidomics_raw), "feature")

lipidomics_data <- lipidomics_raw |> 
  rename(feature_id = feature) |> 
  rename_sample_columns(id_cols = "feature_id")

lipidomics_sample_cols_std <- setdiff(colnames(lipidomics_data), "feature_id")

if (length(lipidomics_sample_cols_raw) != length(lipidomics_sample_cols_std)) {
  stop(glue("Lipidomics sample column mismatch: raw {length(lipidomics_sample_cols_raw)} vs standardized {length(lipidomics_sample_cols_std)}"))
}

lipidomics_features <- read_csv(
  file.path(dataset_dir, "lipiomics", "Lipidom_feature_metadata.cleaned.noUNK.csv"),
  show_col_types = FALSE
) |> 
  rename(feature_id = feature) |> 
  mutate(
    feature_label = coalesce(name, name_raw, feature_id),
    feature_label = strip_bracket_suffix(feature_label),
    feature_label = remove_no_ms2_prefix(feature_label),
    feature_label = str_trim(feature_label),
    feature_display = feature_label
  )

sample_name_mappings[["lipidomics"]] <- tibble(
  dataset = "lipidomics",
  original = lipidomics_sample_cols_raw,
  standardized = lipidomics_sample_cols_std
)

metabolomics_raw <- read_csv(
  file.path(dataset_dir, "metabolomics", "Metabo_X_log2CSS.cleaned.csv"),
  show_col_types = FALSE
)

metabolomics_sample_cols_raw <- setdiff(colnames(metabolomics_raw), "feature")
metabolomics_data <- metabolomics_raw |> 
  rename(feature_id = feature) |> 
  rename_sample_columns(id_cols = "feature_id")

metabolomics_sample_cols_std <- setdiff(colnames(metabolomics_data), "feature_id")

if (length(metabolomics_sample_cols_raw) != length(metabolomics_sample_cols_std)) {
  stop(glue("Metabolomics sample column mismatch: raw {length(metabolomics_sample_cols_raw)} vs standardized {length(metabolomics_sample_cols_std)}"))
}

metabolomics_features <- read_csv(
  file.path(dataset_dir, "metabolomics", "Metabo_feature_metadata.csv"),
  show_col_types = FALSE
) |> 
  rename(feature_id = feature) |> 
  mutate(
    feature_label = coalesce(name, name_raw, feature_id),
    feature_label = strip_bracket_suffix(feature_label),
    feature_label = remove_no_ms2_prefix(feature_label),
    feature_label = str_trim(feature_label),
    feature_display = shorten_metabolite_labels(feature_label)
  )

sample_name_mappings[["metabolomics"]] <- tibble(
  dataset = "metabolomics",
  original = metabolomics_sample_cols_raw,
  standardized = metabolomics_sample_cols_std
)

peptidomics_raw <- read_csv(
  file.path(dataset_dir, "Peptidomics", "Peptidomics_Imputed_WithProteinIDs.csv"),
  show_col_types = FALSE
)

names(peptidomics_raw)[1:2] <- c("peptide_sequence", "protein_group")

peptidomics_prepared <- peptidomics_raw |> 
  mutate(feature_id = paste(protein_group, peptide_sequence, sep = "|")) |> 
  dplyr::select(feature_id, protein_group, peptide_sequence, dplyr::everything())

peptidomics_sample_cols_raw <- setdiff(colnames(peptidomics_prepared), c("feature_id", "protein_group", "peptide_sequence"))

peptidomics_data <- peptidomics_prepared |> 
  rename_sample_columns(id_cols = c("feature_id", "protein_group", "peptide_sequence"))

peptidomics_annotations <- annotate_protein_groups(peptidomics_data$protein_group)

if (nrow(peptidomics_annotations) > 0) {
  peptidomics_data <- peptidomics_data |> left_join(peptidomics_annotations, by = "protein_group")
}

peptidomics_data <- peptidomics_data |> 
  mutate(
    feature_label = dplyr::case_when(
      !is.na(protein_label) & protein_label != "" ~ protein_label,
      !is.na(gene_symbol) & gene_symbol != "" ~ gene_symbol,
      !is.na(primary_id) & primary_id != "" ~ primary_id,
      TRUE ~ peptide_sequence
    ),
    feature_display = feature_label
  )

peptidomics_metadata_cols <- c(
  "feature_id", "protein_group", "peptide_sequence", "primary_id",
  "gene_symbol", "protein_label", "protein_description", "feature_label",
  "feature_display"
)

peptidomics_sample_cols_std <- setdiff(colnames(peptidomics_data), peptidomics_metadata_cols)

if (length(peptidomics_sample_cols_raw) != length(peptidomics_sample_cols_std)) {
  stop(glue("Peptidomics sample column mismatch: raw {length(peptidomics_sample_cols_raw)} vs standardized {length(peptidomics_sample_cols_std)}"))
}

peptidomics_assay_df <- peptidomics_data |> 
  dplyr::select(feature_id, dplyr::all_of(peptidomics_sample_cols_std))

peptidomics_features <- peptidomics_data |> 
  dplyr::select(dplyr::any_of(peptidomics_metadata_cols))

sample_name_mappings[["peptidomics"]] <- tibble(
  dataset = "peptidomics",
  original = peptidomics_sample_cols_raw,
  standardized = peptidomics_sample_cols_std
)

proteomics_raw <- read_csv(
  file.path(dataset_dir, "Proteomics", "Proteomics_imputed_WithProteinIDs.csv"),
  show_col_types = FALSE
)

names(proteomics_raw)[1] <- "protein_group"

proteomics_prepared <- proteomics_raw |> 
  rename(peptide_sequence = Stripped.Sequence) |> 
  mutate(feature_id = paste(protein_group, peptide_sequence, sep = "|")) |> 
  dplyr::select(feature_id, protein_group, peptide_sequence, dplyr::everything())

proteomics_sample_cols_raw <- setdiff(colnames(proteomics_prepared), c("feature_id", "protein_group", "peptide_sequence"))

proteomics_data <- proteomics_prepared |> 
  rename_sample_columns(id_cols = c("feature_id", "protein_group", "peptide_sequence"))

proteomics_annotations <- annotate_protein_groups(proteomics_data$protein_group)

if (nrow(proteomics_annotations) > 0) {
  proteomics_data <- proteomics_data |> left_join(proteomics_annotations, by = "protein_group")
}

proteomics_data <- proteomics_data |> 
  mutate(
    feature_label = dplyr::case_when(
      !is.na(protein_label) & protein_label != "" ~ protein_label,
      TRUE ~ protein_group
    ),
    feature_display = feature_label
  )

proteomics_metadata_cols <- c(
  "feature_id", "protein_group", "peptide_sequence", "primary_id",
  "gene_symbol", "protein_label", "protein_description", "feature_label",
  "feature_display"
)

proteomics_sample_cols_std <- setdiff(colnames(proteomics_data), proteomics_metadata_cols)

if (length(proteomics_sample_cols_raw) != length(proteomics_sample_cols_std)) {
  stop(glue("Proteomics sample column mismatch: raw {length(proteomics_sample_cols_raw)} vs standardized {length(proteomics_sample_cols_std)}"))
}

proteomics_assay_df <- proteomics_data |> 
  dplyr::select(feature_id, dplyr::all_of(proteomics_sample_cols_std))

proteomics_features <- proteomics_data |> 
  dplyr::select(dplyr::any_of(proteomics_metadata_cols))

sample_name_mappings[["proteomics"]] <- tibble(
  dataset = "proteomics",
  original = proteomics_sample_cols_raw,
  standardized = proteomics_sample_cols_std
)

sample_name_mapping_tbl <- dplyr::bind_rows(sample_name_mappings)
write_tibble(sample_name_mapping_tbl, "sample_name_mapping", results_root)

# ---- Sample metadata -------------------------------------------------------
common_samples <- Reduce(intersect, list(
  lipidomics_sample_cols_std,
  metabolomics_sample_cols_std,
  peptidomics_sample_cols_std,
  proteomics_sample_cols_std
))

expected_order <- c(
  sprintf("High-grade_%d", 1:3),
  sprintf("Medium-grade_%d", 1:3),
  sprintf("Low-grade_%d", 1:3)
)

ordered_samples <- expected_order[expected_order %in% common_samples]
if (length(ordered_samples) != length(common_samples)) {
  extra_samples <- setdiff(common_samples, ordered_samples)
  if (length(extra_samples) > 0) {
    ordered_samples <- c(ordered_samples, sort(extra_samples))
  }
}

sample_metadata <- tibble(
  sample_id = ordered_samples,
  group = str_remove(sample_id, "_[0-9]+$"),
  replicate = as.integer(str_extract(sample_id, "[0-9]+$")),
  sample_label = if_else(is.na(replicate), group, glue("{group} #{replicate}")),
  color = grade_colors[group]
)

if (any(is.na(sample_metadata$color))) {
  missing_groups <- unique(sample_metadata$group[is.na(sample_metadata$color)])
  stop("Missing colour definitions for groups: ", paste(missing_groups, collapse = ", "))
}

write_tibble(sample_metadata, "sample_metadata", results_root)

# ---- Assay matrices --------------------------------------------------------
lipidomics_mat <- matrix_from_df(lipidomics_data, "feature_id")[, ordered_samples, drop = FALSE]
metabolomics_mat <- matrix_from_df(metabolomics_data, "feature_id")[, ordered_samples, drop = FALSE]
peptidomics_mat <- matrix_from_df(peptidomics_assay_df, "feature_id")[, ordered_samples, drop = FALSE]
proteomics_mat <- matrix_from_df(proteomics_assay_df, "feature_id")[, ordered_samples, drop = FALSE]

assays <- list(
  lipidomics = lipidomics_mat,
  metabolomics = metabolomics_mat,
  peptidomics = peptidomics_mat,
  proteomics = proteomics_mat
)

assays_features <- list(
  lipidomics = lipidomics_features,
  metabolomics = metabolomics_features,
  peptidomics = peptidomics_features,
  proteomics = proteomics_features
)

for (nm in names(assays)) {
  filtered_mat <- drop_zero_variance_features(assays[[nm]])
  removed <- nrow(assays[[nm]]) - nrow(filtered_mat)
  assays[[nm]] <- filtered_mat
  if (!is.null(assays_features[[nm]])) {
    assays_features[[nm]] <- assays_features[[nm]] |> 
      dplyr::filter(feature_id %in% rownames(filtered_mat))
  }
  if (removed > 0) {
    message(glue("Removed {removed} zero-variance feature(s) from {nm}."))
  }
}

assays_features <- purrr::map(assays_features, augment_feature_display_column)

feature_label_lookup <- purrr::imap_dfr(assays_features, function(df, dataset_name) {
  if (is.null(df) || nrow(df) == 0) {
    return(tibble())
  }
  df |> 
    mutate(dataset = dataset_name) |> 
    dplyr::select(dataset, feature_id, feature_label, feature_display)
})

stopifnot(all(purrr::map_int(assays, ~ncol(.x)) == nrow(sample_metadata)))

feature_summary <- tibble(
  dataset = names(assays),
  n_features = purrr::map_int(assays, nrow),
  n_samples = purrr::map_int(assays, ncol)
)
write_tibble(feature_summary, "feature_summary", results_root)

expression_heatmap_dir <- file.path(results_root, "Heatmaps", "Expression")
purrr::walk(names(assays), function(dataset_name) {
  plot_expression_heatmap(
    mat = assays[[dataset_name]],
    dataset_name = dataset_name,
    sample_metadata = sample_metadata,
    feature_metadata = assays_features[[dataset_name]],
    output_dir = expression_heatmap_dir
  )
})

diablo_feature_lookup <- purrr::imap(assays, function(mat, dataset_name) {
  feature_ids <- rownames(mat)
  if (length(feature_ids) == 0) {
    return(tibble(feature_id = character(), feature_label = character(), feature_display = character()))
  }

  meta_df <- assays_features[[dataset_name]]
  if (is.null(meta_df) || nrow(meta_df) == 0) {
    return(tibble(
      feature_id = feature_ids,
      feature_label = feature_ids,
      feature_display = feature_ids
    ))
  }

  lookup_df <- meta_df |> 
    dplyr::select(feature_id, feature_label, feature_display) |> 
    distinct()

  missing_ids <- setdiff(feature_ids, lookup_df$feature_id)
  if (length(missing_ids) > 0) {
    lookup_df <- dplyr::bind_rows(
      lookup_df,
      tibble(
        feature_id = missing_ids,
        feature_label = missing_ids,
        feature_display = missing_ids
      )
    )
  }

  lookup_df |> 
    mutate(
      feature_label = coalesce(feature_label, feature_id),
      feature_display = if_else(is.na(feature_display) | feature_display == "", feature_label, feature_display)
    )
})

# ---- Principal Component Analysis -----------------------------------------
run_pca <- function(mat, assay_name) {
  sample_matrix <- t(mat)
  zero_var <- which(apply(sample_matrix, 2, function(col) stats::sd(col, na.rm = TRUE)) == 0)
  if (length(zero_var) > 0) {
    sample_matrix <- sample_matrix[, -zero_var, drop = FALSE]
  }
  if (ncol(sample_matrix) < 2 || nrow(sample_matrix) < 2) {
    warning(glue("Skipping PCA for {assay_name}: insufficient variability after filtering."))
    return(invisible(NULL))
  }
  max_rank <- min(n_pca_components, ncol(sample_matrix), nrow(sample_matrix))
  pca_fit <- prcomp(sample_matrix, center = TRUE, scale. = TRUE, rank. = max_rank)

  scores_df <- as_tibble(pca_fit$x, rownames = "sample_id") |> 
    mutate(sample_id = factor(sample_id, levels = sample_metadata$sample_id)) |> 
    left_join(sample_metadata, by = "sample_id")

  loadings_df <- as_tibble(pca_fit$rotation, rownames = "feature_id")

  feature_meta <- assays_features[[assay_name]]
  if (!is.null(feature_meta)) {
    loadings_df <- loadings_df |> 
  left_join(feature_meta |> dplyr::select(feature_id, feature_label), by = "feature_id")
  } else {
    loadings_df <- loadings_df |> mutate(feature_label = feature_id)
  }

  variance_df <- tibble(
    component = paste0("PC", seq_along(pca_fit$sdev)),
    variance_explained = (pca_fit$sdev ^ 2) / sum(pca_fit$sdev ^ 2)
  )

  write_tibble(scores_df, glue("{assay_name}_pca_scores"), file.path(results_root, "PCA"))
  write_tibble(loadings_df, glue("{assay_name}_pca_loadings"), file.path(results_root, "PCA"))
  write_tibble(variance_df, glue("{assay_name}_pca_variance"), file.path(results_root, "PCA"))

  if (nrow(variance_df) >= 2) {
    subtitle_text <- glue(
      "PC1 ({scales::percent(variance_df$variance_explained[1], accuracy = 0.1)}) vs PC2 ({scales::percent(variance_df$variance_explained[2], accuracy = 0.1)})"
    )
  } else {
    subtitle_text <- glue("PC1 ({scales::percent(variance_df$variance_explained[1], accuracy = 0.1)})")
  }

  p <- ggplot(scores_df, aes(x = PC1, y = PC2, colour = group)) +
    geom_point(size = 3) +
    scale_colour_manual(values = grade_colors, drop = FALSE) +
    coord_equal() +
    labs(
      title = glue("{dataset_label(assay_name)} PCA"),
      subtitle = subtitle_text,
      x = "PC1",
      y = "PC2",
      colour = "Group"
    ) +
    theme_minimal()

  save_plot(p, glue("{assay_name}_pca_PC1_PC2"), file.path(results_root, "PCA"))

  if (all(c("PC1", "PC2") %in% colnames(loadings_df))) {
    top_loadings <- loadings_df |> 
      dplyr::select(feature_id, feature_label, PC1, PC2) |> 
      mutate(loading_norm = sqrt(PC1^2 + PC2^2)) |> 
      arrange(dplyr::desc(loading_norm)) |> 
      slice_head(n = 15)

    biplot_loadings <- top_loadings |> 
      mutate(
        PC1_scaled = scales::rescale(PC1, to = range(scores_df$PC1, na.rm = TRUE)),
        PC2_scaled = scales::rescale(PC2, to = range(scores_df$PC2, na.rm = TRUE))
      )

    biplot <- ggplot(scores_df, aes(x = PC1, y = PC2, colour = group)) +
      geom_point(size = 3) +
      geom_segment(
        data = biplot_loadings,
        aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled),
        inherit.aes = FALSE,
        arrow = grid::arrow(length = grid::unit(0.1, "inches")),
        colour = "grey40"
      ) +
      geom_label_repel(
        data = biplot_loadings,
        aes(x = PC1_scaled, y = PC2_scaled, label = feature_label),
        inherit.aes = FALSE,
        size = 3,
        colour = "grey20"
      ) +
      scale_colour_manual(values = grade_colors, drop = FALSE) +
      coord_equal() +
      labs(
        title = glue("{dataset_label(assay_name)} PCA biplot"),
        x = "PC1",
        y = "PC2",
        colour = "Group"
      ) +
      theme_minimal()

    save_plot(biplot, glue("{assay_name}_pca_biplot"), file.path(results_root, "PCA"), width = 7, height = 5)
  }

  plot_loading_heatmap(loadings_df, assay_name, file.path(results_root, "PCA"))
}

purrr::walk2(assays, names(assays), run_pca)

# ---- Differential analysis -------------------------------------------------
differential_dir <- file.path(results_root, "Differential")
differential_results <- purrr::imap(assays, function(mat, dataset_name) {
  run_limma_differential(mat, dataset_name, sample_metadata, assays_features[[dataset_name]], differential_dir)
})

differential_heatmap_dir <- file.path(results_root, "Heatmaps", "Differential")
purrr::iwalk(differential_results, function(diff_tbl, dataset_name) {
  plot_differential_heatmap(diff_tbl, dataset_name, differential_heatmap_dir)
})

# ---- Correlations and coexpression ----------------------------------------
correlation_dir <- file.path(results_root, "Correlation")
purrr::walk(names(assays), function(dataset_name) {
  plot_sample_correlations(assays[[dataset_name]], dataset_name, sample_metadata, correlation_dir)
})

plot_cross_dataset_correlations(assays, sample_metadata, correlation_dir)

network_dir <- file.path(results_root, "Networks")
purrr::walk(names(assays), function(dataset_name) {
  build_coexpression_network(assays[[dataset_name]], dataset_name, network_dir, assays_features[[dataset_name]])
  plot_coexpression_heatmap(
    mat = assays[[dataset_name]],
    dataset_name = dataset_name,
    output_dir = file.path(results_root, "Heatmaps", "Coexpression"),
    feature_metadata = assays_features[[dataset_name]]
  )
})

analyze_multiomics_coexpression(
  assays = assays,
  feature_metadata_list = assays_features,
  sample_metadata = sample_metadata,
  heatmap_output_dir = file.path(results_root, "Heatmaps", "Coexpression"),
  network_output_dir = network_dir,
  top_n_per_dataset = 60,
  cor_threshold = 0.6,
  max_edges = 300,
  min_edges = 45,
  fallback_quantile = 0.8
)

# ---- WGCNA -----------------------------------------------------------------
wgcna_dir <- file.path(results_root, "WGCNA")
purrr::walk(names(assays), function(dataset_name) {
  run_wgcna_analysis(
    mat = assays[[dataset_name]],
    dataset_name = dataset_name,
    output_dir = wgcna_dir,
    sample_metadata = sample_metadata,
    feature_metadata = assays_features[[dataset_name]]
  )
})

# ---- Pathway enrichment ----------------------------------------------------
pathway_dir <- file.path(results_root, "Pathways")
purrr::iwalk(differential_results, function(diff_tbl, dataset_name) {
  run_pathway_enrichment(diff_tbl, assays_features[[dataset_name]], dataset_name, pathway_dir)
})

# ---- DIABLO / mixOmics -----------------------------------------------------
diablo_input <- purrr::map(assays, ~ t(.x) |> as.data.frame())
names(diablo_input) <- names(assays)
for (nm in names(diablo_input)) {
  rownames(diablo_input[[nm]]) <- sample_metadata$sample_id

  mapping_df <- diablo_feature_lookup[[nm]]
  if (!is.null(mapping_df) && nrow(mapping_df) > 0) {
    display_vec <- setNames(mapping_df$feature_display, mapping_df$feature_id)
    current_cols <- colnames(diablo_input[[nm]])
    new_cols <- display_vec[current_cols]
    replace_idx <- is.na(new_cols) | new_cols == ""
    if (any(replace_idx)) {
      new_cols[replace_idx] <- current_cols[replace_idx]
    }
    if (anyDuplicated(new_cols)) {
      new_cols <- make.unique(new_cols, sep = " ")
    }
    colnames(diablo_input[[nm]]) <- new_cols

    id_to_display <- setNames(new_cols, current_cols)
    diablo_feature_lookup[[nm]] <- diablo_feature_lookup[[nm]] |> 
      mutate(feature_display = coalesce(id_to_display[feature_id], feature_display))
  }
}

keepX_list <- purrr::map(diablo_input, ~ rep(min(ncol(.x), max_keepX), n_diablo_components))
design_matrix <- matrix(1, nrow = length(diablo_input), ncol = length(diablo_input))
diag(design_matrix) <- 0

diablo_model <- block.splsda(
  X = diablo_input,
  Y = factor(sample_metadata$group, levels = names(grade_colors)),
  ncomp = n_diablo_components,
  keepX = keepX_list,
  design = design_matrix
)

diablo_plot_obj <- diablo_model
diablo_plot_obj$Y <- factor(sample_metadata$group, levels = names(grade_colors))

dir_create(file.path(results_root, "DIABLO"), recurse = TRUE)
saveRDS(diablo_model, file.path(results_root, "DIABLO", "diablo_model.rds"))

diablo_scores_wide <- purrr::imap_dfr(diablo_model$variates, function(mat, view) {
  if (!(view %in% names(diablo_input))) {
    return(tibble())
  }

  as_tibble(mat, rownames = "sample_id") |> 
    mutate(dataset = view, dataset_label = dataset_label(view)) |> 
    left_join(sample_metadata, by = "sample_id")
})

diablo_scores_long <- diablo_scores_wide |> 
  pivot_longer(cols = starts_with("comp"), names_to = "component", values_to = "score") |> 
  mutate(component_index = readr::parse_number(component))

write_tibble(diablo_scores_wide, "diablo_scores_wide", file.path(results_root, "DIABLO"))
write_tibble(diablo_scores_long, "diablo_scores_long", file.path(results_root, "DIABLO"))

selected_features <- purrr::imap(diablo_input, function(.x, view) {
  mapping_df <- diablo_feature_lookup[[view]]
  display_to_id <- NULL
  display_to_label <- NULL
  if (!is.null(mapping_df) && nrow(mapping_df) > 0) {
    display_to_id <- setNames(mapping_df$feature_id, mapping_df$feature_display)
    display_to_label <- setNames(mapping_df$feature_label, mapping_df$feature_display)
  }

  purrr::map_dfr(seq_len(n_diablo_components), function(comp) {
    vars <- selectVar(diablo_model, comp = comp, block = view)
    feature_display <- vars$name
    feature_id <- if (!is.null(display_to_id)) unname(display_to_id[feature_display]) else NA_character_
    feature_label <- if (!is.null(display_to_label)) unname(display_to_label[feature_display]) else NA_character_

    tibble(
      feature_display = feature_display,
      feature_id = coalesce(feature_id, feature_display),
      feature_label = coalesce(feature_label, feature_display),
      loading = vars$value,
      component = comp,
      dataset = view
    )
  })
}) |> dplyr::bind_rows()

selected_features <- selected_features |> mutate(dataset_label = dataset_label(dataset))

write_tibble(selected_features, "diablo_selected_features", file.path(results_root, "DIABLO"))

if (n_diablo_components >= 2 && all(c("comp1", "comp2") %in% colnames(diablo_scores_wide))) {
  diablo_plot <- ggplot(diablo_scores_wide, aes(x = comp1, y = comp2, colour = group)) +
    geom_point(size = 3) +
    facet_wrap(~dataset) +
    scale_colour_manual(values = grade_colors, drop = FALSE) +
    labs(
      title = "DIABLO component space",
      subtitle = "Component 1 vs Component 2 across omics blocks",
      x = "Component 1",
      y = "Component 2",
      colour = "Group"
    ) +
    theme_minimal()

  save_plot(diablo_plot, "diablo_component_scatter", file.path(results_root, "DIABLO"))
}

group_factor <- factor(sample_metadata$group, levels = names(grade_colors))
names(group_factor) <- sample_metadata$sample_id

if (n_diablo_components >= 2) {
  safe_mixomics_plot(
    file_stub = "diablo_individuals_comp1_vs_comp2",
    output_dir = file.path(results_root, "DIABLO"),
    width = 7,
    height = 7,
    plot_fn = function() {
      mixOmics::plotIndiv(
        diablo_plot_obj,
        comp = c(1, 2),
        group = group_factor,
        legend = TRUE,
        ind.names = FALSE,
        title = "DIABLO individuals (components 1 & 2)",
        col = grade_colors[levels(group_factor)]
      )
    }
  )

  safe_mixomics_plot(
    file_stub = "diablo_arrow_comp1_vs_comp2",
    output_dir = file.path(results_root, "DIABLO"),
    width = 7,
    height = 7,
    plot_fn = function() {
      mixOmics::plotArrow(
        diablo_plot_obj,
        comp = c(1, 2),
        group = group_factor,
        ind.names = FALSE,
        legend = TRUE,
        title = "DIABLO arrow plot"
      )
    }
  )

  safe_mixomics_plot(
    file_stub = "diablo_variable_correlation_comp1_comp2",
    output_dir = file.path(results_root, "DIABLO"),
    width = 7,
    height = 7,
    plot_fn = function() {
      mixOmics::plotVar(
        diablo_model,
        comp = c(1, 2),
        cutoff = 0.6,
        var.names = FALSE,
        legend = TRUE,
        title = "DIABLO variable correlations"
      )
    }
  )
}

plot_diablo_component_correlations(diablo_scores_wide, file.path(results_root, "Correlation"))
plot_diablo_loadings(diablo_model, file.path(results_root, "DIABLO"), feature_lookup = diablo_feature_lookup)
draw_diablo_circos(diablo_model, file.path(results_root, "Circos"))
run_diablo_performance(diablo_model, file.path(results_root, "DIABLO"))
plot_diablo_cim(diablo_model, file.path(results_root, "DIABLO"))
plot_diablo_networks(diablo_model, file.path(results_root, "DIABLO"))

# ---- MOFA2 -----------------------------------------------------------------
mofa_views <- purrr::map(assays, ~ .x)
names(mofa_views) <- names(assays)

mofa_obj <- create_mofa(mofa_views)
mofa_sample_metadata <- sample_metadata |> 
  transmute(
    sample = sample_id,
    group = "group1",
    grade = .data$group,
    replicate = replicate,
    color = color
  )

samples_metadata(mofa_obj) <- as.data.frame(mofa_sample_metadata)

model_opts <- get_default_model_options(mofa_obj)
model_opts$num_factors <- n_mofa_factors

data_opts <- get_default_data_options(mofa_obj)
training_opts <- get_default_training_options(mofa_obj)
training_opts$verbose <- TRUE

mofa_prepared <- prepare_mofa(
  mofa_obj,
  model_options = model_opts,
  training_options = training_opts,
  data_options = data_opts
)

mofa_model <- run_mofa(mofa_prepared, use_basilisk = TRUE)

dir_create(file.path(results_root, "MOFA2"), recurse = TRUE)
saveRDS(mofa_model, file.path(results_root, "MOFA2", "mofa_model.rds"))

factors_df <- get_factors(mofa_model, factors = "all", as.data.frame = TRUE) |> 
  as_tibble() |> 
  rename(sample_id = sample)

weights_list <- get_weights(mofa_model, factors = "all", as.data.frame = FALSE)
weights_df <- purrr::imap_dfr(weights_list, function(view_mat, view_name) {
  as_tibble(view_mat, rownames = "feature_id") |> 
    pivot_longer(cols = -feature_id, names_to = "factor", values_to = "weight") |> 
    mutate(dataset = view_name)
})

variance_explained <- get_variance_explained(mofa_model)

variance_total_df <- purrr::imap_dfr(
  variance_explained$r2_total,
  function(vec, grp) {
    tibble(
      group = grp,
      dataset = names(vec),
      variance_explained = as.numeric(vec)
    )
  }
)

variance_per_factor_df <- purrr::imap_dfr(
  variance_explained$r2_per_factor,
  function(mat, grp) {
    as_tibble(mat, rownames = "factor") |> 
      pivot_longer(
        cols = -factor,
        names_to = "dataset",
        values_to = "variance_explained"
      ) |> 
      mutate(group = grp)
  }
)

write_tibble(factors_df, "mofa_factors", file.path(results_root, "MOFA2"))
write_tibble(weights_df, "mofa_feature_weights", file.path(results_root, "MOFA2"))
write_tibble(variance_total_df, "mofa_variance_explained_total", file.path(results_root, "MOFA2"))
write_tibble(variance_per_factor_df, "mofa_variance_explained_per_factor", file.path(results_root, "MOFA2"))

if (n_mofa_factors >= 2) {
  factor_plot <- plot_factors(
    mofa_model,
    factors = 1:2,
    color_by = "grade",
    dot_size = 3
  ) + scale_colour_manual(values = grade_colors)

  save_plot(factor_plot, "mofa_factors_1_2", file.path(results_root, "MOFA2"))
}

plot_mofa_variance_overview(variance_total_df, file.path(results_root, "MOFA2"))
plot_mofa_variance_heatmap(variance_per_factor_df, file.path(results_root, "MOFA2"))
plot_mofa_factor_boxplots(factors_df, sample_metadata, file.path(results_root, "MOFA2"))
plot_mofa_factor_heatmap(factors_df, sample_metadata, file.path(results_root, "MOFA2"))
plot_mofa_top_weights(weights_df, assays_features, file.path(results_root, "MOFA2"))
summarise_mofa_factor_covariates(factors_df, sample_metadata, file.path(results_root, "MOFA2"))

# ---- MetaboAnalyst exports ------------------------------------------------
metaboanalyst_exports <- prepare_metaboanalyst_exports(
  assays = assays,
  feature_metadata_list = assays_features,
  sample_metadata = sample_metadata,
  output_dir = file.path(results_root, "MetaboAnalyst")
)

if (!is.null(metaboanalyst_exports) && !is.null(metaboanalyst_exports$expression_long)) {
  plot_metaboanalyst_overview(
    expression_long = metaboanalyst_exports$expression_long,
    sample_metadata = sample_metadata,
    output_dir = file.path(results_root, "MetaboAnalyst")
  )
}

message("Multi-omics analysis pipeline completed. Results saved to: ", results_root)
