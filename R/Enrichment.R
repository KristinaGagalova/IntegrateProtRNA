## Enrichment.R -- shared fgsea helpers, ported from norinXcadenza-shared
## (analysis/rnaseq/wheat/Enrichment.R). Kept behaviourally identical to the
## upstream version so gene-set membership rules (specific/generic/mixed)
## match the existing enrichment write-up exactly; only sourced here, not
## modified.

check_id_overlap_annotations <- function(

    ranked_df,
    annotation_df,
    rank_gene_col = "Majority.protein.IDs",
    assign_gene_col = "IDENTIFIER",
    term_col = "NAME",
    not_assigned_terms = c(
      "not assigned.not annotated",
      "not assigned.annotated"
    ),
    generic_patterns = c(
      "^uncharacterised",
      "^uncharacterized"
    ),
    id_sep = ";",
    remove_cds_suffix = TRUE
) {
  # dependencies
  require(dplyr)
  require(tidyr)

  # helper: generic term check
  is_generic_term <- function(x, patterns = generic_patterns) {
    x <- tolower(trimws(x))
    ifelse(
      is.na(x) | x == "",
      FALSE,
      Reduce(`|`, lapply(patterns, grepl, x = x, ignore.case = TRUE))
    )
  }

  # helper: normalize annotation term
  normalize_term <- function(x) {
    tolower(trimws(x))
  }

  # ----------------------------
  # 1. Prepare ranked table
  # ----------------------------
  rank_tmp <- ranked_df %>%
    mutate(
      row_id = row_number(),
      original_ids = .data[[rank_gene_col]]
    ) %>%
    select(row_id, original_ids)

  rank_long <- rank_tmp %>%
    separate_rows(original_ids, sep = id_sep) %>%
    mutate(
      original_ids = tolower(trimws(original_ids))
    )

  if (remove_cds_suffix) {
    rank_long <- rank_long %>%
      mutate(
        original_ids = gsub("\\.cds[0-9]+$", "", original_ids)
      )
  }

  # ----------------------------
  # 2. Prepare annotation table
  # ----------------------------
  annot_tmp <- annotation_df %>%
    transmute(
      gene_id = tolower(trimws(.data[[assign_gene_col]])),
      term = trimws(.data[[term_col]])
    ) %>%
    distinct() %>%
    mutate(
      term_norm = normalize_term(term),
      is_non_annotated_term = is.na(term_norm) |
        term_norm == "" |
        term_norm %in% tolower(not_assigned_terms),
      is_generic = ifelse(is_non_annotated_term, FALSE, is_generic_term(term)),
      is_specific = !is_non_annotated_term & !is_generic
    )

  # ----------------------------
  # 3. Join ranked IDs to annotations
  # ----------------------------
  joined <- rank_long %>%
    left_join(annot_tmp, by = c("original_ids" = "gene_id")) %>%
    mutate(
      is_missing_join = is.na(term),
      is_non_annotated = ifelse(is_missing_join, TRUE, is_non_annotated_term),
      is_generic = ifelse(is.na(is_generic), FALSE, is_generic),
      is_specific = ifelse(is.na(is_specific), FALSE, is_specific)
    )

  # ----------------------------
  # 4. Summarise per original row
  # ----------------------------
  summary_df <- joined %>%
    group_by(row_id) %>%
    summarise(
      n_ids = n(),

      n_non_annotated_raw = sum(is_non_annotated, na.rm = TRUE),
      n_generic = sum(is_generic, na.rm = TRUE),
      n_specific = sum(is_specific, na.rm = TRUE),

      n_unique_generic_terms = n_distinct(term[is_generic], na.rm = TRUE),
      n_unique_specific_terms = n_distinct(term[is_specific], na.rm = TRUE),

      generic_terms = paste(sort(unique(na.omit(term[is_generic]))), collapse = " | "),
      specific_terms = paste(sort(unique(na.omit(term[is_specific]))), collapse = " | "),
      non_annotated_terms = paste(
        sort(unique(na.omit(term[is_non_annotated & !is_missing_join]))),
        collapse = " | "
      ),

      has_annotated = (n_unique_generic_terms + n_unique_specific_terms) > 0,

      n_non_annotated = if_else(has_annotated, 0L, n_non_annotated_raw),

      status = case_when(
        !has_annotated ~ "non.annotated",

        # IMPORTANT: prioritize multiple specific terms
        n_unique_specific_terms > 1 ~ "mixed_terms",

        # then check specific + generic
        n_unique_specific_terms > 0 & n_unique_generic_terms > 0 ~ "mixed_generic",

        n_unique_specific_terms == 1 ~ "specific",

        n_unique_specific_terms == 0 & n_unique_generic_terms > 0 ~ "generic",

        TRUE ~ "other"
      ),

      keep_for_fgsea = case_when(
        status %in% c("specific", "generic", "mixed_generic") ~ TRUE,
        TRUE ~ FALSE
      ),

      .groups = "drop"
    )

  # ----------------------------
  # 5. Attach back to original table
  # ----------------------------
  checked_df <- ranked_df %>%
    mutate(row_id = row_number()) %>%
    left_join(summary_df, by = "row_id")

  # return everything useful
  list(
    checked_df = checked_df,
    summary_df = summary_df,
    joined_df = joined,
    annotation_df_processed = annot_tmp
  )
}


make_fgsea_annotation_table <- function(
    checked_df,
    id_col = "Majority.protein.IDs"
) {
  checked_df %>%
    mutate(
      term_for_fgsea = case_when(
        status == "generic" ~ generic_terms,
        status == "mixed_generic" ~ specific_terms,
        status == "specific" ~ specific_terms,
        TRUE ~ NA_character_
      )
    ) %>%
    select(
      row_id,
      gene_ids = all_of(id_col),
      status,
      keep_for_fgsea,
      term_for_fgsea,
      generic_terms,
      specific_terms,
      non_annotated_terms
    )
}


plot_fgsea_bubble <- function(fgsea_res,
                              top_n = 20,
                              plot_title = NULL,
                              padj_cutoff = 0.05,
                              wrap_width = 55,
                              debug = TRUE) {

  parse_pathway_label <- function(x,
                                  sep = "\\.",
                                  connector = " - ",
                                  keep_first = TRUE,
                                  keep_last_n = 2,
                                  max_length = 60,
                                  truncate_token = " ... ") {

    truncate_middle <- function(text,
                                max_length = 40,
                                token = " ... ") {

      if (nchar(text) <= max_length) {
        return(text)
      }

      keep_chars <- floor((max_length - nchar(token)) / 2)

      start_part <- substr(text, 1, keep_chars)

      end_part <- substr(
        text,
        nchar(text) - keep_chars + 1,
        nchar(text)
      )

      paste0(start_part, token, end_part)
    }

    sapply(x, function(term) {

      parts <- unlist(strsplit(term, sep))
      parts <- trimws(parts)
      parts <- parts[parts != ""]

      if (length(parts) == 0) {
        return(NA_character_)
      }

      if (length(parts) <= keep_last_n + 1) {

        label <- paste(parts, collapse = connector)

      } else {

        first_part <- parts[1]

        last_parts <- tail(parts, keep_last_n)

        label <- paste(
          c(first_part, last_parts),
          collapse = connector
        )
      }

      truncate_middle(
        label,
        max_length = max_length,
        token = truncate_token
      )
    })
  }

  df <- as.data.frame(fgsea_res)

  if (debug) {
    message("\n==============================")
    message("Plot: ", plot_title)
    message("Initial fgsea rows: ", nrow(df))
    message("Columns: ", paste(colnames(df), collapse = ", "))
  }

  if (nrow(df) == 0) {
    message("Stopping: fgsea result has 0 rows")
    return(NULL)
  }

  required_cols <- c("pathway", "padj", "NES", "size")

  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  df <- df %>%
    mutate(
      pathway_original = as.character(pathway),
      pathway = parse_pathway_label(pathway_original),
      #pathway = stringr::str_wrap(pathway, width = wrap_width),
      neglog10_padj = -log10(padj),
      direction = ifelse(NES >= 0, "Up", "Down")
    )

  if (debug) {
    message("Rows before padj filter: ", nrow(df))

    if (!is.null(padj_cutoff)) {
      message("padj cutoff: ", padj_cutoff)
      message("Rows passing padj cutoff: ", sum(df$padj <= padj_cutoff, na.rm = TRUE))
      message("padj summary:")
      print(summary(df$padj))
    } else {
      message("padj cutoff: NULL, no padj filtering")
    }
  }

  if (!is.null(padj_cutoff)) {
    df <- df %>%
      filter(!is.na(padj), padj <= padj_cutoff)
  }

  if (nrow(df) == 0) {
    message("Stopping: no rows left after padj filter for plot: ", plot_title)
    return(NULL)
  }

  df <- df %>%
    arrange(padj) %>%
    slice_head(n = top_n)

  if (debug) {
    message("Rows after top_n filter: ", nrow(df))
  }

  dup_pathways <- df$pathway[duplicated(df$pathway)]

  if (debug && length(dup_pathways) > 0) {
    message("\nDuplicated pathway labels detected in plot: ", plot_title)
    print(unique(dup_pathways))

    message("\nOriginal pathways causing duplicated labels:")
    print(
      df %>%
        filter(pathway %in% unique(dup_pathways)) %>%
        select(pathway_original, pathway, NES, padj, size)
    )
  }

  df <- df %>%
    mutate(
      pathway_label = pathway,
      pathway_plot = make.unique(pathway),
      pathway_plot = factor(pathway_plot, levels = pathway_plot[order(NES)])
    )

  if (debug) {
    message("Building plot now.")
  }

  p1 <- ggplot(df, aes(x = NES, y = pathway_plot)) +

    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.45,
      color = "grey40"
    ) +

    geom_point(
      aes(size = size, fill = neglog10_padj),
      shape = 21,
      color = "black",
      stroke = 0.35,
      alpha = 0.9
    ) +

    scale_y_discrete(
      labels = setNames(df$pathway_label, df$pathway_plot)
    ) +

    scale_size_continuous(
      range = c(3, 10),
      name = "Pathway size"
    ) +

    scale_fill_viridis_c(
      option = "magma",
      name = expression(-log[10](padj))
    ) +

    labs(
      title = plot_title,
      x = "Normalized enrichment score (NES)",
      y = NULL
    ) +

    theme_bw(base_size = 14) +

    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(size = 13, color = "black"),
      axis.text.y = element_text(size = 15, color = "black", lineheight = 1.05),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12, color = "black"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.6)
    )

  return(p1)
}
