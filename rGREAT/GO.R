# -------------------------------
# 0. Install/load required packages
# -------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("rGREAT", quietly = TRUE))
  BiocManager::install("rGREAT")
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!requireNamespace("stringr", quietly = TRUE))
  install.packages("stringr")
if (!requireNamespace("GenomicRanges", quietly = TRUE))
  BiocManager::install("GenomicRanges")
if (!requireNamespace("IRanges", quietly = TRUE))
  BiocManager::install("IRanges")

library(rGREAT)
library(dplyr)
library(ggplot2)
library(stringr)
library(GenomicRanges)
library(IRanges)

# -------------------------------
# 1. Configuration
# -------------------------------
setwd("/Users/oukanyou/Desktop/mapping_results")

peak_files <- list(
  human_specific = "human_specific_peaks_conservative.narrowPeak",
  mouse_specific = "mouse_specific_peaks_conservative.narrowPeak",
  shared         = "shared_peaks_conservative.narrowPeak",
  mouse_to_human = "mouse_to_human_conservative.narrowPeak"
)

species_map <- list(
  human_specific = "hg38",
  mouse_specific = "mm10",
  shared         = "hg38",
  mouse_to_human = "hg38"
)

# Chromosome size limits — used to filter out-of-bounds peaks
# Prevents GREAT error when peaks exceed actual chrom lengths
CHROM_SIZES <- list(
  hg38 = c(
    chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
    chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
    chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
    chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
    chr17=83257441,  chr18=80373285,  chr19=58617616,  chr20=64444167,
    chr21=46709983,  chr22=50818468,  chrX=156040895,  chrY=57227415,
    chrM=16569
  ),
  mm10 = c(
    chr1=195471971,  chr2=182113224,  chr3=160039680,  chr4=156508116,
    chr5=151834684,  chr6=149736546,  chr7=145441459,  chr8=129401213,
    chr9=124595110,  chr10=130694993, chr11=122082543, chr12=120129022,
    chr13=120421639, chr14=124902244, chr15=104043685, chr16=98207768,
    chr17=94987271,  chr18=90702639,  chr19=61431566,  chrX=171031299,
    chrY=91744698,   chrM=16299
  )
)

# -------------------------------
# 2. Helper: read narrowPeak → GRanges with bounds filtering
# -------------------------------
read_narrowpeak <- function(filepath, genome) {
  df <- read.table(filepath, header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "")
  colnames(df)[1:3] <- c("chr", "start", "end")
  
  # 1. Keep standard chromosomes only
  df <- df[grepl("^chr([0-9]+|X|Y|M)$", df$chr), ]
  
  # 2. Filter peaks that exceed known chromosome boundaries
  chrom_sizes <- CHROM_SIZES[[genome]]
  if (!is.null(chrom_sizes)) {
    before <- nrow(df)
    df <- df[
      !is.na(chrom_sizes[df$chr]) &          # chr exists in size table
        df$end <= chrom_sizes[df$chr],          # end within chromosome
    ]
    removed <- before - nrow(df)
    if (removed > 0)
      cat(sprintf("  [Filter] Removed %d out-of-bounds peaks for %s\n",
                  removed, genome))
  }
  
  if (nrow(df) == 0) return(NULL)
  
  GRanges(
    seqnames = df$chr,
    ranges   = IRanges(start = df$start + 1,   # BED 0-based → 1-based
                       end   = df$end)
  )
}

# -------------------------------
# 3. Helper: extract GO sub-tables from GREAT v4 "GO" category
#    GREAT v4 returns one "GO" category containing BP, MF, CC together.
#    getEnrichmentTables() returns a list; inspect its names to split them.
# -------------------------------
extract_go_subtables <- function(job, cats) {
  go_bp <- NULL; go_mf <- NULL; go_cc <- NULL
  
  # ---- GREAT v4: single "GO" category ----
  if ("GO" %in% cats) {
    tbl_list <- tryCatch(
      getEnrichmentTables(job, category = "GO"),
      error = function(e) { warning("GO fetch error: ", e$message); NULL }
    )
    if (!is.null(tbl_list)) {
      cat("  GO sub-tables available:", paste(names(tbl_list), collapse=", "), "\n")
      # Names are typically "GO Biological Process", "GO Molecular Function",
      # "GO Cellular Component" in GREAT v4
      for (nm in names(tbl_list)) {
        if (grepl("Biological|BP", nm, ignore.case = TRUE))  go_bp <- tbl_list[[nm]]
        if (grepl("Molecular|MF",  nm, ignore.case = TRUE))  go_mf <- tbl_list[[nm]]
        if (grepl("Cellular|CC",   nm, ignore.case = TRUE))  go_cc <- tbl_list[[nm]]
      }
      # Fallback: if only one table returned, assign to BP
      if (is.null(go_bp) && is.null(go_mf) && is.null(go_cc) &&
          length(tbl_list) >= 1) {
        go_bp <- tbl_list[[1]]
        cat("  Warning: Could not identify GO sub-type; assigned first table to BP\n")
      }
    }
  }
  
  # ---- GREAT v3 fallback: separate "GO:BP", "GO:MF", "GO:CC" categories ----
  safe_get <- function(category) {
    if (!(category %in% cats)) return(NULL)
    tryCatch({
      tbl <- getEnrichmentTables(job, category = category)
      tbl[[1]]
    }, error = function(e) { warning("Fetch error (", category, "): ", e$message); NULL })
  }
  if (is.null(go_bp)) go_bp <- safe_get("GO:BP")
  if (is.null(go_mf)) go_mf <- safe_get("GO:MF")
  if (is.null(go_cc)) go_cc <- safe_get("GO:CC")
  if (is.null(go_bp)) go_bp <- safe_get("GO Biological Process")
  if (is.null(go_mf)) go_mf <- safe_get("GO Molecular Function")
  if (is.null(go_cc)) go_cc <- safe_get("GO Cellular Component")
  
  list(BP = go_bp, MF = go_mf, CC = go_cc)
}

# -------------------------------
# 4. Plotting helper
# -------------------------------
plot_top_terms <- function(tbl, ontology, prefix) {
  if (is.null(tbl) || nrow(tbl) == 0) return(invisible(NULL))
  
  # Detect p-value column (varies by GREAT version)
  pval_col <- intersect(
    colnames(tbl),
    c("Hyper_Raw_PValue", "p_value", "pvalue",
      "Binom_Raw_PValue", "BinomP", "HyperP")
  )
  if (length(pval_col) == 0) {
    # Last resort: pick any column whose name contains "pval" or "pvalue"
    pval_col <- grep("p.?val", colnames(tbl), ignore.case = TRUE, value = TRUE)
  }
  if (length(pval_col) == 0) {
    cat("  No p-value column found for", ontology, "— skipping plot\n")
    cat("  Available columns:", paste(colnames(tbl), collapse=", "), "\n")
    return(invisible(NULL))
  }
  pval_col <- pval_col[1]
  
  name_col <- intersect(colnames(tbl), c("name", "Term", "description", "ID"))
  if (length(name_col) == 0) name_col <- colnames(tbl)[1]
  name_col <- name_col[1]
  
  top <- tbl %>%
    filter(!is.na(.data[[pval_col]]), .data[[pval_col]] < 0.05) %>%
    arrange(.data[[pval_col]]) %>%
    head(20) %>%
    mutate(
      log10p   = -log10(.data[[pval_col]]),
      term_lab = str_wrap(.data[[name_col]], width = 45)
    )
  
  if (nrow(top) == 0) {
    cat("  No significant terms (p<0.05) for", ontology, "\n")
    return(invisible(NULL))
  }
  
  p <- ggplot(top, aes(x = log10p, y = reorder(term_lab, log10p))) +
    geom_col(fill = "#2C7BB6", alpha = 0.85) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               colour = "red", linewidth = 0.5) +
    labs(
      title    = paste0(prefix, "  |  ", ontology),
      subtitle = paste0("Top ", nrow(top), " significant terms (p < 0.05)"),
      x        = expression(-log[10](p-value)),
      y        = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title         = element_text(face = "bold"),
      axis.text.y        = element_text(size = 8),
      panel.grid.major.y = element_blank()
    )
  
  out_png <- paste0(prefix, "_", gsub("[ /]", "_", ontology), "_top20.png")
  ggsave(out_png, plot = p, width = 10, height = 7, dpi = 150)
  cat("  Saved plot:", out_png, "\n")
}

# -------------------------------
# 5. Core GREAT function
# -------------------------------
run_great <- function(filepath, species, prefix) {
  cat("\n========================================\n")
  cat("Running GREAT for:", filepath, "| Genome:", species, "\n")
  cat("========================================\n")
  
  if (!file.exists(filepath)) {
    warning("File not found, skipping: ", filepath)
    return(NULL)
  }
  
  # Read & filter peaks
  gr <- tryCatch(
    read_narrowpeak(filepath, species),
    error = function(e) {
      warning("Failed to read: ", filepath, "\n  ", e$message)
      NULL
    }
  )
  if (is.null(gr) || length(gr) == 0) {
    warning("No valid peaks after filtering: ", filepath)
    return(NULL)
  }
  cat("  Peaks after filtering:", length(gr), "\n")
  
  # Submit to GREAT
  job <- tryCatch(
    submitGreatJob(
      gr,
      species          = species,
      rule             = "basalPlusExt",
      adv_upstream     = 5.0,
      adv_downstream   = 1.0,
      adv_span         = 1000.0,
      request_interval = 30,
      help             = FALSE        # suppress version message
    ),
    error = function(e) {
      warning("submitGreatJob failed: ", filepath, "\n  ", e$message)
      NULL
    }
  )
  if (is.null(job)) return(NULL)
  
  # Available categories
  cats <- availableCategories(job)
  cat("  Available categories:", paste(cats, collapse = ", "), "\n")
  
  # Extract GO tables
  go_tables <- extract_go_subtables(job, cats)
  go_bp <- go_tables$BP
  go_mf <- go_tables$MF
  go_cc <- go_tables$CC
  
  # Report what we got
  cat(sprintf("  GO BP: %s terms | GO MF: %s terms | GO CC: %s terms\n",
              if (!is.null(go_bp)) nrow(go_bp) else "NULL",
              if (!is.null(go_mf)) nrow(go_mf) else "NULL",
              if (!is.null(go_cc)) nrow(go_cc) else "NULL"))
  
  # Save CSVs
  save_csv <- function(tbl, suffix) {
    if (!is.null(tbl) && nrow(tbl) > 0) {
      out <- paste0(prefix, "_", suffix, ".csv")
      write.csv(tbl, out, row.names = FALSE)
      cat("  Saved:", out, "\n")
    }
  }
  save_csv(go_bp, "GO_BP")
  save_csv(go_mf, "GO_MF")
  save_csv(go_cc, "GO_CC")
  
  # Plots
  plot_top_terms(go_bp, "GO Biological Process", prefix)
  plot_top_terms(go_mf, "GO Molecular Function", prefix)
  plot_top_terms(go_cc, "GO Cellular Component",  prefix)
  
  list(BP = go_bp, MF = go_mf, CC = go_cc, job = job)
}

# -------------------------------
# 6. Run for all peak sets
# -------------------------------
results <- list()
for (name in names(peak_files)) {
  results[[name]] <- run_great(
    filepath = peak_files[[name]],
    species  = species_map[[name]],
    prefix   = paste0("GREAT_", name)
  )
}

# -------------------------------
# 7. Summary
# -------------------------------
cat("\n========== SUMMARY ==========\n")
for (name in names(results)) {
  res <- results[[name]]
  if (is.null(res)) {
    cat(sprintf("  %-20s : FAILED\n", name))
  } else {
    bp_n <- if (!is.null(res$BP)) nrow(res$BP) else 0
    mf_n <- if (!is.null(res$MF)) nrow(res$MF) else 0
    cc_n <- if (!is.null(res$CC)) nrow(res$CC) else 0
    cat(sprintf("  %-20s | BP: %4d | MF: %4d | CC: %4d terms\n",
                name, bp_n, mf_n, cc_n))
  }
}
cat("==============================\n")