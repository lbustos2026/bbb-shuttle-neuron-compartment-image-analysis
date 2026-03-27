# ================================================================================== #
# Figure 8 - NeuN+ Neuron mean Cy5-signal in Subcellular Compartment Analysis 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# This code is usable for segmented NeuN+ Neuron Somas using the IN Carta image analysis
# software. This code is using subcellular object data from segmented neurons and  
# the mean Cy5 signal in those objects (subcellular compartments).
# *** Key Details: 
#     - This script normalizes each compartment to Group A (naive) mean
#     - Detects small high-intensity right-tail populations by group (makes QC plots)
#     - Writes per-well & group-level summary tables
#     - Creates combined compartment overlay plots as seen in Figure 8D
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# LBustos - Updated 2026 :)
# ================================================================================= #

# ============================================================================ #
# --------------------------------- Setup ------------------------------------ #
# ============================================================================ #

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(rlang)
  library(ggplot2)
  library(scales)
  library(grid)
  library(tibble)
})


project_dir <- "." # <-----------enter your file path with the Data folder here
raw_root <- file.path(project_dir, "Data")
output_root <- file.path(project_dir, "Outputs")
combined_output_dir <- file.path(output_root, "Compartment_Analysis")
plots_dir <- file.path(combined_output_dir, "Plots")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(combined_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

well_label_filter <- NULL  # example: "A - 1"

meta_cols <- c("Plate ID", "WELL LABEL", "Row", "Column", "FOV", "Z", "T", "Object ID")

# ==================================================================================== #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~ Labels, Palettes & Data Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
# ==================================================================================== #

group_levels <- c(
  "Naive",
  "1hr Post Dose 1",
  "3hrs Post Dose 1",
  "24hrs Post Dose 1",
  "1wk Post Dose 1",
  "1wk Post Dose 2"
)

# Group info - these are the group names for the imaging experiments
group_map <- c(
  A = "Naive",
  G = "1hr Post Dose 1",
  H = "3hrs Post Dose 1",
  I = "24hrs Post Dose 1",
  J = "1wk Post Dose 1",
  K = "1wk Post Dose 2"
)

compartment_specs <- list(
  EE = list(
    short_name   = "EE",
    display_name = "EEA1+ Compartment",
    input_dir    = file.path(raw_root, "NeuNEEA1Cy5-HPRT"),
    output_dir   = file.path(output_root, "EE"),
    object_label = "Early Endosome",
    legend_name  = "Early Endosome",
    file_pattern = "Early Endosome_ObjectData.*\\.csv$",
    glob_stub    = "*Early Endosome_ObjectData*.csv",
    target_col   = "Early Endosome Mean Intensity - Local Background 3 : None (2) - Cy5",
    fuzzy_regex  = "Early\\s*Endosome\\s*Mean\\s*Intensity.*Local\\s*Background.*3.*Cy5",
    value_col    = "EE_mean_intensity3_bg_sub",
    norm_col     = "EE_mean_intensity3_bg_sub_normA",
    unit_name    = "obj"
  ),
  RE = list(
    short_name   = "RE",
    display_name = "RAB11+ Compartment",
    input_dir    = file.path(raw_root, "NeuNRAB11Cy5-HPRT"),
    output_dir   = file.path(output_root, "RE"),
    object_label = "RecycEndo",
    legend_name  = "RecycEndo",
    file_pattern = "RecycEndo_ObjectData.*\\.csv$",
    glob_stub    = "*RecycEndo_ObjectData*.csv",
    target_col   = "RecycEndo Mean Intensity - Local Background 3 : None (2) - Cy5",
    fuzzy_regex  = "RecycEndo\\s*Mean\\s*Intensity.*Local\\s*Background.*3.*Cy5",
    value_col    = "RE_mean_intensity3_bg_sub",
    norm_col     = "RE_mean_intensity3_bg_sub_normA",
    unit_name    = "obj"
  ),
  Lyso = list(
    short_name   = "Lyso",
    display_name = "LAMP1+ Compartment",
    input_dir    = file.path(raw_root, "NeuNLAMP1Cy5-HPRT"),
    output_dir   = file.path(output_root, "Lyso"),
    object_label = "Lysosome",
    legend_name  = "Lysosome",
    file_pattern = "Lyso_ObjectData.*\\.csv$",
    glob_stub    = "*Lyso_ObjectData*.csv",
    target_col   = "Lyso Mean Intensity - Local Background 3 : None (2) - Cy5",
    fuzzy_regex  = "Lyso\\s*Mean\\s*Intensity.*Local\\s*Background.*3.*Cy5",
    value_col    = "Lyso_mean_intensity3_bg_sub",
    norm_col     = "Lyso_mean_intensity3_bg_sub_normA",
    unit_name    = "lyso"
  )
)

comp_pal <- c(
  "EEA1+ Compartment"  = "#94A323",
  "RAB11+ Compartment" = "#F3A546",
  "LAMP1+ Compartment" = "#CF3E53"
)

# =================================================================================== #
# ----------------------- (Lots of) Helper Functions ------------------------------- #
# =================================================================================== #

pretty_dose_labels <- function(i) {
  group_levels[i] |>
    stringr::str_replace("1hr Post Dose 1",  "1hr Post\nDose 1") |>
    stringr::str_replace("3hrs Post Dose 1", "3hrs Post\nDose 1") |>
    stringr::str_replace("24hrs Post Dose 1","24hrs Post\nDose 1") |>
    stringr::str_replace("1wk Post Dose 1",  "1wk Post\nDose 1") |>
    stringr::str_replace("1wk Post Dose 2",  "1wk Post\nDose 2")
}

find_bimodal_tail_cutoff <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 200) return(NA_real_)
  
  lx <- log10(x)
  d  <- density(lx, n = 512)
  ys <- d$y
  xs <- d$x
  
  peaks <- which(diff(sign(diff(ys))) == -2)
  if (length(peaks) < 2) return(NA_real_)
  
  main_peak_idx <- peaks[which.max(ys[peaks])]
  tail_candidates <- peaks[peaks > main_peak_idx & ys[peaks] < ys[main_peak_idx]]
  if (length(tail_candidates) == 0) return(NA_real_)
  
  tail_peak_idx <- tail_candidates[1]
  valley_range <- which(xs > xs[main_peak_idx] & xs < xs[tail_peak_idx])
  if (length(valley_range) == 0) return(NA_real_)
  
  valley_idx <- valley_range[which.min(ys[valley_range])]
  10 ^ xs[valley_idx]
}

get_input_files <- function(spec) {
  stopifnot(dir.exists(spec$input_dir))
  
  files <- list.files(
    path = spec$input_dir,
    pattern = spec$file_pattern,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  files <- unique(c(
    files,
    Sys.glob(file.path(spec$input_dir, "*", spec$glob_stub)),
    Sys.glob(file.path(spec$input_dir, "*", toupper(spec$glob_stub)))
  ))
  
  files
}

read_one_object_file <- function(fp, spec) {
  grp <- basename(dirname(fp))
  
  df <- tryCatch(
    readr::read_csv(fp, show_col_types = FALSE),
    error = function(e) {
      warning(sprintf("Failed to read %s: %s", fp, e$message))
      return(NULL)
    }
  )
  if (is.null(df)) return(NULL)
  
  if (!("WELL LABEL" %in% names(df))) {
    well_hit <- names(df)[stringr::str_detect(
      names(df), stringr::regex("^WELL\\s*LAB", ignore_case = TRUE)
    )]
    if (length(well_hit)) {
      df <- dplyr::rename(df, `WELL LABEL` = !!rlang::sym(well_hit[1]))
    }
  }
  
  value_col_found <- if (spec$target_col %in% names(df)) {
    spec$target_col
  } else {
    cand <- names(df)[stringr::str_detect(
      names(df),
      stringr::regex(spec$fuzzy_regex, ignore_case = TRUE)
    )]
    if (length(cand) == 0) {
      warning(sprintf("No matching Cy5 BG-subtracted mean intensity column found in %s", basename(fp)))
      return(NULL)
    }
    cand[1]
  }
  
  out <- df %>%
    dplyr::select(dplyr::any_of(c(meta_cols, value_col_found))) %>%
    dplyr::rename(!!spec$value_col := !!rlang::sym(value_col_found)) %>%
    dplyr::mutate(
      !!spec$value_col := suppressWarnings(as.numeric(.data[[spec$value_col]])),
      .group       = grp,
      .source_file = basename(fp),
      .source_path = fp,
      Compartment  = spec$display_name
    )
  
  out
}

write_plot_set <- function(raw_pos_norm, raw_clean_norm, spec) {
  dir.create(spec$output_dir, recursive = TRUE, showWarnings = FALSE)
  
  norm_col <- spec$norm_col
  
  p_hist_overlay <- ggplot(raw_pos_norm, aes(x = .data[[norm_col]], fill = tail_flag)) +
    geom_histogram(bins = 150, position = "identity", alpha = 0.6) +
    scale_x_log10() +
    facet_wrap(~ .group, scales = "free_y") +
    scale_fill_manual(
      values = c("gray70", "red"),
      labels = c("Main distribution", "Small high-intensity tail"),
      name = spec$legend_name
    ) +
    theme_bw() +
    labs(
      title = paste0("Per-object ", spec$object_label, " Cy5 (normalized to A) - bimodal tail labeling"),
      x = paste0(spec$object_label, " Cy5 mean intensity (bg-sub, fold of A, log10)"),
      y = "Number of objects"
    )
  
  p_violin_overlay <- ggplot(raw_pos_norm, aes(x = .group, y = .data[[norm_col]], fill = tail_flag)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    scale_y_log10() +
    scale_fill_manual(
      values = c("gray70", "red"),
      labels = c("Main distribution", "Small high-intensity tail"),
      name = spec$legend_name
    ) +
    theme_bw() +
    labs(
      title = paste0(spec$object_label, ": Cy5 (normalized to A) - bimodal tail labeling"),
      x = "Group",
      y = paste0(spec$object_label, " Cy5 mean intensity (bg-sub, fold of A, log10)")
    )
  
  p_hist_after <- ggplot(raw_clean_norm, aes(x = .data[[norm_col]])) +
    geom_histogram(bins = 100) +
    scale_x_log10() +
    facet_wrap(~ .group, scales = "free_y") +
    theme_bw() +
    labs(
      title = paste0("Per-object ", spec$object_label, " Cy5 (normalized to A, log10) - small tail removed"),
      x = paste0(spec$object_label, " Cy5 mean intensity (bg-sub, fold of A, log10)"),
      y = "Number of objects"
    )
  
  p_violin_after <- ggplot(raw_clean_norm, aes(x = .group, y = .data[[norm_col]])) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
    scale_y_log10() +
    theme_bw() +
    labs(
      title = paste0(spec$object_label, ": Cy5 (normalized to A, log10) - small tail removed"),
      x = "Group",
      y = paste0(spec$object_label, " Cy5 mean intensity (bg-sub, fold of A, log10)")
    )
  
  print(p_hist_overlay)
  print(p_violin_overlay)
  print(p_hist_after)
  print(p_violin_after)
  
  ggsave(file.path(spec$output_dir, paste0(spec$short_name, "_Cy5_hist_overlay.png")),
         p_hist_overlay, width = 8, height = 4, dpi = 300)
  ggsave(file.path(spec$output_dir, paste0(spec$short_name, "_Cy5_violin_overlay.png")),
         p_violin_overlay, width = 8, height = 4, dpi = 300)
  ggsave(file.path(spec$output_dir, paste0(spec$short_name, "_Cy5_hist_tail_removed.png")),
         p_hist_after, width = 8, height = 4, dpi = 300)
  ggsave(file.path(spec$output_dir, paste0(spec$short_name, "_Cy5_violin_tail_removed.png")),
         p_violin_after, width = 5, height = 4, dpi = 300)
}

summarize_compartment <- function(raw_bimodal, spec) {
  norm_col <- spec$norm_col
  
  total_name <- paste0("total_", spec$unit_name)
  tail_name  <- paste0("tail_",  spec$unit_name)
  main_name  <- paste0("main_",  spec$unit_name)
  
  by_well <- raw_bimodal %>%
    group_by(.group, `WELL LABEL`, Compartment) %>%
    summarise(
      !!total_name := n(),
      !!tail_name  := sum(tail_flag, na.rm = TRUE),
      !!main_name  := sum(!tail_flag, na.rm = TRUE),
      
      mean_normA_main   = mean(.data[[norm_col]][!tail_flag], na.rm = TRUE),
      median_normA_main = median(.data[[norm_col]][!tail_flag], na.rm = TRUE),
      sd_normA_main     = sd(.data[[norm_col]][!tail_flag], na.rm = TRUE),
      
      mean_normA_tail   = mean(.data[[norm_col]][tail_flag], na.rm = TRUE),
      median_normA_tail = median(.data[[norm_col]][tail_flag], na.rm = TRUE),
      sd_normA_tail     = sd(.data[[norm_col]][tail_flag], na.rm = TRUE),
      
      mean_normA_all = (
        .data[[main_name]] * mean_normA_main +
          .data[[tail_name]] * mean_normA_tail
      ) / (.data[[main_name]] + .data[[tail_name]]),
      
      .groups = "drop"
    ) %>%
    mutate(
      across(starts_with("mean_"),   ~ replace_na(.x, 0)),
      across(starts_with("median_"), ~ replace_na(.x, 0)),
      across(starts_with("sd_"),     ~ replace_na(.x, 0)),
      mean_normA_all = replace_na(mean_normA_all, 0)
    )
  
  qc_summary <- by_well %>%
    group_by(.group, Compartment) %>%
    summarise(
      n_wells = n(),
      !!total_name := sum(.data[[total_name]], na.rm = TRUE),
      !!tail_name  := sum(.data[[tail_name]],  na.rm = TRUE),
      !!main_name  := sum(.data[[main_name]],  na.rm = TRUE),
      pct_tail = 100 * .data[[tail_name]] / .data[[total_name]],
      
      mean_normA_main   = mean(mean_normA_main, na.rm = TRUE),
      sd_normA_main     = if_else(sum(!is.na(mean_normA_main)) >= 2,
                                  sd(mean_normA_main, na.rm = TRUE), 0),
      median_normA_main = median(mean_normA_main, na.rm = TRUE),
      
      mean_normA_tail   = mean(mean_normA_tail, na.rm = TRUE),
      sd_normA_tail     = if_else(sum(!is.na(mean_normA_tail)) >= 2,
                                  sd(mean_normA_tail, na.rm = TRUE), 0),
      median_normA_tail = median(mean_normA_tail, na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    arrange(.group)
  
  list(by_well = by_well, qc_summary = qc_summary)
}

run_compartment <- function(spec) {
  message("\n=========(beep)====(beep)======(boop)===========")
  message("Running: ", spec$short_name)
  message("=========(beep)====(beep)======(boop)===========")
  
  files <- get_input_files(spec)
  if (!length(files)) {
    stop("No matching files found for ", spec$short_name, " in ", spec$input_dir)
  }
  
  message("Found ", length(files), " file(s).")
  
  raw <- purrr::map_dfr(files, read_one_object_file, spec = spec)
  if (!nrow(raw)) stop("No rows loaded for ", spec$short_name)
  
  if (!is.null(well_label_filter)) {
    raw <- raw %>% filter(`WELL LABEL` == well_label_filter)
  }
  
  raw <- raw %>% drop_na(.data[[spec$value_col]])
  
  A_mean <- raw %>%
    filter(.group == "A") %>%
    summarise(A_mean = mean(.data[[spec$value_col]], na.rm = TRUE)) %>%
    pull(A_mean)
  
  if (length(A_mean) == 0 || is.na(A_mean) || A_mean == 0) {
    stop("Could not compute valid Group A mean for ", spec$short_name)
  }
  
  raw <- raw %>%
    mutate(!!spec$norm_col := .data[[spec$value_col]] / A_mean)
  
  cutoffs <- raw %>%
    group_by(.group) %>%
    summarise(
      tail_cutoff_normA = find_bimodal_tail_cutoff(.data[[spec$norm_col]]),
      .groups = "drop"
    )
  
  print(cutoffs)
  
  raw_bimodal <- raw %>%
    left_join(cutoffs, by = ".group") %>%
    mutate(
      tail_flag = !is.na(tail_cutoff_normA) &
        .data[[spec$norm_col]] > tail_cutoff_normA
    )
  
  print(table(raw_bimodal$.group, raw_bimodal$tail_flag))
  
  raw_pos_norm <- raw_bimodal %>%
    filter(.data[[spec$norm_col]] > 0)
  
  raw_clean_norm <- raw_pos_norm %>%
    filter(!tail_flag)
  
  write_plot_set(raw_pos_norm, raw_clean_norm, spec)
  
  summaries <- summarize_compartment(raw_bimodal, spec)
  
  readr::write_csv(
    summaries$by_well,
    file.path(spec$output_dir, paste0(spec$short_name, "_Cy5_byWell_main_vs_tail.csv"))
  )
  
  readr::write_csv(
    summaries$qc_summary,
    file.path(spec$output_dir, paste0(spec$short_name, "_Cy5_groupSummary_main_vs_tail.csv"))
  )
  
  list(
    raw = raw,
    raw_bimodal = raw_bimodal,
    by_well = summaries$by_well,
    qc_summary = summaries$qc_summary
  )
}

make_combined_compartment_plots <- function(bywell_all) {
  long_dat <- bywell_all %>%
    mutate(
      GroupLabel = recode(as.character(.group), !!!group_map, .default = as.character(.group)),
      GroupLabel = factor(GroupLabel, levels = group_levels),
      Compartment = factor(
        Compartment,
        levels = c("EEA1+ Compartment", "RAB11+ Compartment", "LAMP1+ Compartment")
      ),
      Value = mean_normA_main
    ) %>%
    filter(!is.na(Value))
  
  summary_df <- long_dat %>%
    group_by(GroupLabel, Compartment) %>%
    summarise(
      mean = mean(Value, na.rm = TRUE),
      sd   = sd(Value, na.rm = TRUE),
      n    = n(),
      .groups = "drop"
    )
  
  global_y_min <- min(summary_df$mean - summary_df$sd, na.rm = TRUE)
  global_y_max <- max(summary_df$mean + summary_df$sd, na.rm = TRUE)
  y_margin <- 0.05 * (global_y_max - global_y_min)
  global_y_limits <- c(global_y_min - y_margin, global_y_max + y_margin)
  
  pos <- position_jitter(width = 0.05, height = 0, seed = 42)
  
  compartments_to_plot <- levels(summary_df$Compartment)
  compartments_to_plot <- compartments_to_plot[compartments_to_plot %in% summary_df$Compartment]
  
  for (compartment_name in compartments_to_plot) {
    sum_df_comp <- summary_df %>%
      filter(Compartment == compartment_name) %>%
      mutate(x_num = as.numeric(GroupLabel))
    
    long_dat_comp <- long_dat %>%
      filter(Compartment == compartment_name) %>%
      mutate(x_num = as.numeric(GroupLabel))
    
    this_col <- comp_pal[as.character(compartment_name)]
    
    x  <- sum_df_comp$x_num
    xs <- seq(min(x), max(x), length.out = 200)
    
    lo_mean <- loess(mean ~ x, data = sum_df_comp, span = 0.8)
    lo_lo   <- loess(I(mean - sd) ~ x, data = sum_df_comp, span = 0.8)
    lo_hi   <- loess(I(mean + sd) ~ x, data = sum_df_comp, span = 0.8)
    
    smooth_df <- data.frame(
      x_num  = xs,
      mean_s = predict(lo_mean, xs),
      lo_s   = predict(lo_lo, xs),
      hi_s   = predict(lo_hi, xs)
    )
    
    p_comp <- ggplot() +
      geom_ribbon(
        data = smooth_df,
        aes(x = x_num, ymin = lo_s, ymax = hi_s),
        fill = scales::alpha(this_col, 0.20),
        color = NA
      ) +
      geom_line(
        data = smooth_df,
        aes(x = x_num, y = mean_s),
        color = this_col,
        linewidth = 0.7
      ) +
      geom_point(
        data = sum_df_comp,
        aes(x = x_num, y = mean),
        color = this_col,
        size = 1.5,
        shape = 15
      ) +
      geom_point(
        data = long_dat_comp,
        aes(x = x_num, y = Value),
        position = pos,
        size = 1.3,
        alpha = 0.55,
        color = this_col
      ) +
      geom_point(
        data = long_dat_comp,
        aes(x = x_num, y = Value),
        position = pos,
        shape = 21,
        fill = NA,
        colour = "black",
        stroke = 0.3,
        size = 1.75
      ) +
      scale_x_continuous(
        breaks = seq_along(group_levels),
        labels = pretty_dose_labels,
        expand = expansion(mult = c(0.03, 0.08))
      ) +
      scale_y_continuous(
        limits = global_y_limits,
        labels = scales::label_number(big.mark = ","),
        breaks = scales::breaks_pretty(9)
      ) +
      labs(x = NULL, y = NULL, title = NULL) +
      theme_classic(base_size = 12) +
      theme(
        panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.position    = "none",
        axis.text.x        = element_text(size = 9, face = "bold", color = "black"),
        axis.text.y        = element_text(size = 8, face = "bold", color = "black"),
        plot.margin        = margin(5, 5, 5, 5),
        axis.line.x.bottom = element_line(linewidth = 0.25),
        axis.line.y.left   = element_line(linewidth = 0.25),
        axis.ticks         = element_line(linewidth = 0.5)
      )
    
    print(p_comp)
    
    short_name <- gsub("\\+?\\s*Compartment", "", compartment_name)
    short_name <- gsub("\\s+", "", short_name)
    
    ggsave(
      filename = file.path(plots_dir, paste0("siRNAComp_meanBGsub_norm_", short_name, ".tif")),
      plot = p_comp, width = 7, height = 2.5, dpi = 1200
    )
    
    ggsave(
      filename = file.path(plots_dir, paste0("siRNAComp_meanBGsub_norm_", short_name, "_big.tif")),
      plot = p_comp, width = 7, height = 4, dpi = 1200
    )
  }
  
  make_smoothed_overlay <- function(comp_subset, file_stub) {
    sum_sub <- summary_df %>%
      filter(Compartment %in% comp_subset) %>%
      mutate(x_num = as.numeric(GroupLabel))
    
    long_sub <- long_dat %>%
      filter(Compartment %in% comp_subset) %>%
      mutate(x_num = as.numeric(GroupLabel))
    
    pal_sub <- comp_pal[comp_subset]
    
    smooth_df_all <- sum_sub %>%
      group_by(Compartment) %>%
      group_split() %>%
      map_dfr(function(df) {
        this_comp <- df$Compartment[1]
        x  <- df$x_num
        xs <- seq(min(x), max(x), length.out = 200)
        
        lo_mean <- loess(mean ~ x, data = df, span = 0.8)
        lo_lo   <- loess(I(mean - sd) ~ x, data = df, span = 0.8)
        lo_hi   <- loess(I(mean + sd) ~ x, data = df, span = 0.8)
        
        tibble(
          Compartment = this_comp,
          x_num  = xs,
          mean_s = predict(lo_mean, xs),
          lo_s   = predict(lo_lo, xs),
          hi_s   = predict(lo_hi, xs)
        )
      })
    
    p <- ggplot() +
      geom_ribbon(
        data = smooth_df_all,
        aes(x = x_num, ymin = lo_s, ymax = hi_s, fill = Compartment),
        alpha = 0.20,
        color = NA
      ) +
      geom_line(
        data = smooth_df_all,
        aes(x = x_num, y = mean_s, color = Compartment),
        linewidth = 0.7
      ) +
      geom_point(
        data = sum_sub,
        aes(x = x_num, y = mean, color = Compartment),
        size = 1.5,
        shape = 15
      ) +
      geom_point(
        data = long_sub,
        aes(x = x_num, y = Value, color = Compartment),
        position = pos,
        size = 1.3,
        alpha = 0.55
      ) +
      geom_point(
        data = long_sub,
        aes(x = x_num, y = Value),
        position = pos,
        shape = 21,
        fill = NA,
        colour = "black",
        stroke = 0.3,
        size = 1.75
      ) +
      scale_x_continuous(
        breaks = seq_along(group_levels),
        labels = pretty_dose_labels,
        expand = expansion(mult = c(0.03, 0.05))
      ) +
      scale_y_continuous(
        limits = global_y_limits,
        labels = scales::label_number(big.mark = ","),
        breaks = scales::breaks_pretty(9)
      ) +
      scale_colour_manual(values = pal_sub, breaks = names(pal_sub), guide = guide_legend(title = NULL)) +
      scale_fill_manual(values = scales::alpha(pal_sub, 0.20), breaks = names(pal_sub), guide = guide_legend(title = NULL)) +
      labs(x = NULL, y = NULL, title = NULL) +
      theme_classic(base_size = 12) +
      theme(
        panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.position    = "top",
        legend.title       = element_blank(),
        legend.text        = element_text(size = 8, face = "bold"),
        legend.key.width   = unit(10, "pt"),
        legend.key.height  = unit(10, "pt"),
        legend.margin      = margin(t = 0, r = 0, b = -5, l = 0),
        axis.text.x        = element_text(size = 9, face = "bold", color = "black"),
        axis.text.y        = element_text(size = 8, face = "bold", color = "black"),
        plot.margin        = margin(5, 5, 5, 5),
        axis.line.x.bottom = element_line(linewidth = 0.25),
        axis.line.y.left   = element_line(linewidth = 0.25),
        axis.ticks         = element_line(linewidth = 0.5)
      )
    
    print(p)
    
    ggsave(
      filename = file.path(plots_dir, paste0(file_stub, ".tif")),
      plot = p, width = 7, height = 2.5, dpi = 1200
    )
    
    ggsave(
      filename = file.path(plots_dir, paste0(file_stub, "_big.tif")),
      plot = p, width = 7, height = 4, dpi = 1200
    )
  }
  
  make_smoothed_overlay(
    comp_subset = c("EEA1+ Compartment"),
    file_stub   = "Fig8_meanBGsub_norm_smoothedOverlay_EEA1"
  )
  
  make_smoothed_overlay(
    comp_subset = c("EEA1+ Compartment", "RAB11+ Compartment"),
    file_stub   = "Fig8_meanBGsub_norm_smoothedOverlay_EEA1_RAB11"
  )
  
  make_smoothed_overlay(
    comp_subset = c("EEA1+ Compartment", "RAB11+ Compartment", "LAMP1+ Compartment"),
    file_stub   = "Fig8_meanBGsub_norm_smoothedOverlay_all3"
  )
  
  invisible(list(long_dat = long_dat, summary_df = summary_df))
}

# ============================================================================= #
# ~~~~~~~~~~~~~~~~~~ Run EVERYTHING for ALL Compartments ~~~~~~~~~~~~~~~~~~~~~~ #
# ============================================================================= #

results <- purrr::map(compartment_specs, run_compartment)

bywell_all <- purrr::map_dfr(results, "by_well")
groupsummary_all <- purrr::map_dfr(results, "qc_summary")

readr::write_csv(
  bywell_all,
  file.path(combined_output_dir, "all_compartments_byWell_main_vs_tail.csv")
)

readr::write_csv(
  groupsummary_all,
  file.path(combined_output_dir, "all_compartments_groupSummary_main_vs_tail.csv")
)

make_combined_compartment_plots(bywell_all)

message("\nYay! Analysis and Plotting complete.")
