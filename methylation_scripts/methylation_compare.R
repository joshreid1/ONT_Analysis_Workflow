#!/usr/bin/env Rscript

# -----------------------------
# 1. Load libraries + CLI args
# -----------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(biomaRt)
  library(Rsamtools)
  library(GenomicRanges)
  library(IRanges)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(parallel)
  library(pheatmap)
  library(gprofiler2)
})

# Define CLI options
option_list <- list(
  make_option(c("--sample_a"), type = "character", help = "Sample A .bed.gz file (modkit output)"),
  make_option(c("--sample_b"), type = "character", help = "Sample B .bed.gz file (modkit output)"),
  make_option(c("--name_a"), type = "character", default = "SampleA", help = "Sample A name"),
  make_option(c("--name_b"), type = "character", default = "SampleB", help = "Sample B name"),
  make_option(c("--coverage_threshold"), type = "integer", default = 0, help = "Minimum coverage filter"),
  make_option(c("--region_bed"), type = "character", default = NULL, help = "Optional BED file of regions (e.g. CpG islands)"),
  make_option(c("--gene_list"), type = "character", help = "Path to TSV gene list or 'all_protein_coding_genes'"),
  make_option(c("--add_promotor"), action = "store_true", default = TRUE, help = "Include 2000bp upstream of gene TSS (strand-aware). [default: %default]"),
  make_option(c("--outdir"), type = "character", default = ".", help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Assign args to variables
sample_a_file <- opt$sample_a
sample_b_file <- opt$sample_b
sample_a_name <- opt$name_a
sample_b_name <- opt$name_b
coverage_threshold <- opt$coverage_threshold
region_bed_file <- opt$region_bed
gene_list_path <- opt$gene_list
add_promotor <- opt$add_promotor
outdir <- normalizePath(opt$outdir, mustWork = FALSE)

# The following example overriding values is disabled by default. Uncomment for local debugging only.
#TEST = TRUE
if (TEST) {
  sample_a_file <- "/vast/projects/reidj-project/nanopore/methylation/38703_genome_wide_modkit.bed.gz"
  sample_b_file <- "/vast/projects/reidj-project/nanopore/methylation/38704_genome_wide_modkit.bed.gz"
  sample_a_name <- '38703'
  sample_b_name <- '38704'
  
  sample_a_file <- "/vast/projects/reidj-project/nanopore/methylation/T1648-13_genome_wide_modkit.bed.gz"
  sample_b_file <- "/vast/projects/reidj-project/nanopore/methylation/8454-9_genome_wide_modkit.bed.gz"
  sample_a_name <- 'T1648'
  sample_b_name <- '8454'
  
  coverage_threshold <- 10
  region_bed_file <- "/vast/projects/reidj-project/Genes4Epilepsy/v2025-03/EpilepsyGenes_v2025-03.bed"
  gene_list_path <- "/vast/projects/reidj-project/Genes4Epilepsy/v2025-03/EpilepsyGenes_NamesOnly_v2025-03.tsv"
  gene_list_path <- 'all_protein_coding_genes'
  add_promotor <- TRUE
  outdir <- "/vast/scratch/users/reid.j/methylation"
}

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# parameter validation
if (is.null(sample_a_file) || sample_a_file == "" || is.null(sample_b_file) || sample_b_file == "") {
  stop("--sample_a and --sample_b are required")
}
if (!file.exists(sample_a_file)) {
  stop("Sample A path not found: ", sample_a_file)
}
if (!file.exists(sample_b_file)) {
  stop("Sample B path not found: ", sample_b_file)
}
if (!is.null(region_bed_file) && region_bed_file != "" && !file.exists(region_bed_file)) {
  stop("region_bed file path not found: ", region_bed_file)
}

# bedMethyl columns:
# 1. "chrom" (chromosome name)
# 2. "start" (start position of the methylation site)
# 3. "end" (end position of the methylation site)
# 4. "mod_base_motif" (modification type e.g. m = methylcytosine, h = hydroxymethylcytosine)   
# 5. "score" (equal to Nvalid_cov, but included for compatibility with standard BED format)
# 6. "strand" ("+" or "-", or "." when combined from both strands)
# 7. "start_compat" (start position for compatibility)
# 8. "end_compat" (end position for compatibility)
# 9. "color" (always "255,0,0" for compatibility with BED format)
# 10. "Nvalid_cov" (valid coverage, Nvalid_cov = Nmod + Nother_mod + Ncanonical)
# 11. "percent_modified" (Nmod / Nvalid_cov * 100, percentage of reads showing the modification)
# 12. "Nmod" (number of reads showing the modification, e.g. methylation)
# 13. "Ncanonical" (number of reads showing the canonical base, e.g. unmethylated)
# 14. "Nother_mod" (number of reads showing other modifications, e.g. hydroxymethylation)
# 15. "Ndelete"  (number of reads with a deletion at the site)
# 16. "Nfail" (number of reads where the probability of the call was below the threshold for a confident call, i.e. no-call)
# 17. "Ndiff" (number of reads with a base other than the reference or the called modification, e.g. a SNP or indel at the site) 
# 18. "Nnocall" (Number of reads with correct canonical base but without a modification call, i.e. Ncanonical reads that were not confidently called as modified or unmethylated)

# -----------------------------
# 2. Function Defintions
# -----------------------------

# Utility to log messages with timestamp
log_msg <- function(...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, paste0(...)))
}

extract_methyl_region <- function(bed_path, chr, start, end, coverage_threshold = 0) {
  region <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
  tabix_file <- TabixFile(bed_path)
  lines <- scanTabix(tabix_file, param = region)[[1]]
  
  if (length(lines) == 0) {
    warning("No methylation data found in specified region.")
    return(data.frame())
  }
  
  df <- read.table(text = paste(lines, collapse = "\n"), header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c(
    "chrom", "start", "end", "mod_base_motif", "score", "strand",
    "start_compat", "end_compat", "color", "Nvalid_cov",
    "percent_modified", "Nmod", "Ncanonical", "Nother_mod",
    "Ndelete", "Nfail", "Ndiff", "Nnocall"
  )
  df <- df[df$Nvalid_cov >= coverage_threshold, ]
  return(df)
}

read_bed_regions <- function(bed_path) {
  df <- read.table(bed_path, header = FALSE, stringsAsFactors = FALSE)
  colnames(df)[1:3] <- c("chrom", "start", "end")
  df
}

mean_methylation_per_region <- function(bed_df, sample_bed_file, coverage_threshold = 0) {
  total_regions <- nrow(bed_df)
  log_msg("Starting mean methylation calculation for ", total_regions, " regions (coverage >= ", coverage_threshold, ").")
  means <- numeric(total_regions)
  
  # Progress checkpoints at every 10%
  pct_steps <- seq(0.1, 1, by = 0.1)
  next_pct_idx <- 1
  
  for (i in seq_len(total_regions)) {
    region <- bed_df[i, ]
    df_region <- extract_methyl_region(
      sample_bed_file,
      region$chrom,
      region$start,
      region$end,
      coverage_threshold
    )
    means[i] <- if (nrow(df_region) > 0) mean(df_region$percent_modified, na.rm = TRUE) else NA
    
    progress <- i / total_regions
    if (next_pct_idx <= length(pct_steps) && progress >= pct_steps[next_pct_idx]) {
      log_msg(sprintf("Processed %d/%d (%.0f%%) regions...", i, total_regions, pct_steps[next_pct_idx] * 100))
      next_pct_idx <- next_pct_idx + 1
    }
  }
  log_msg("Completed mean methylation for all regions.")
  data.frame(bed_df, mean_methylation = means)
}

methylation_linear_correlation <- function(
    sample_a_file, sample_b_file,
    coverage_threshold = 0,
    sample_a_name = sample_a_name,
    sample_b_name = sample_b_name,
    plot_main = TRUE,
    chr = NULL,
    start = NULL,
    end = NULL,
    region_bed_file = NULL,
    n_cores = 1
) {
  normalize_chr <- function(x) {
    x <- as.character(x)
    if (length(x) == 0 || is.na(x[1])) return(character(0))
    if (!startsWith(x[1], "chr")) paste0("chr", x) else x
  }

  if (!is.null(region_bed_file)) {
    log_msg("Reading regions from ", region_bed_file)
    bed_regions <- read_bed_regions(region_bed_file)
    log_msg("Loaded ", nrow(bed_regions), " regions from BED file.")

    # Parallelized mean methylation per region
    mean_meth_parallel <- function(bed_df, sample_file) {
      mclapply(seq_len(nrow(bed_df)), function(i) {
        region <- bed_df[i, ]
        df <- extract_methyl_region(sample_file, region$chrom, region$start, region$end, coverage_threshold)
        mean_methyl <- if (nrow(df) > 0) mean(df$percent_modified, na.rm = TRUE) else NA
        return(mean_methyl)
      }, mc.cores = n_cores) %>% unlist()
    }

    bed_regions$mean_methylation_a <- mean_meth_parallel(bed_regions, sample_a_file)
    bed_regions$mean_methylation_b <- mean_meth_parallel(bed_regions, sample_b_file)

    merged_means <- bed_regions
    merged_means <- merged_means[!is.na(merged_means$mean_methylation_a) & !is.na(merged_means$mean_methylation_b), ]
    log_msg("Merged region means: ", nrow(merged_means), " overlapping regions with data.")

    if (nrow(merged_means) == 0) stop("No overlapping regions with methylation data.")

    # Pearson correlation and R²
    cor_val <- cor(merged_means$mean_methylation_a, merged_means$mean_methylation_b, method = "pearson")
    r2 <- cor_val^2
    log_msg(sprintf("Pearson R² for %s vs %s: %.4f", sample_a_name, sample_b_name, r2))

    # Write to TSV
    region_id <- gsub("\\.bed.*$", "", basename(region_bed_file))
    region_id <- gsub("[^A-Za-z0-9_\\-]", "_", region_id)
    outfile <- sprintf("mean_methylation_linear_%s_vs_%s_%s.tsv", sample_a_name, sample_b_name, region_id)

    write.table(
      merged_means,
      file = outfile,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    log_msg("Saved linear region comparison to ", outfile)

    # Plot
    if (plot_main) {
      plotfile <- file.path(outdir, sprintf("linear_comparison_%s_vs_%s_%s.png", sample_a_name, sample_b_name, region_id))
      png(plotfile, width = 1800, height = 1600, res = 300)
      plot_data <- merged_means %>%
      filter(
        !is.na(mean_methylation_a), !is.na(mean_methylation_b),
        mean_methylation_a >= 0, mean_methylation_a <= 100,
        mean_methylation_b >= 0, mean_methylation_b <= 100
      )

       plot(
        plot_data$mean_methylation_a,
        plot_data$mean_methylation_b,
        pch = 20,
        main = sprintf("Mean Methylation per Region\n(R² = %.3f)", r2),
        xlab = sample_a_name,
        ylab = sample_b_name,
        xlim = c(0, 100),
        ylim = c(0, 100),
        asp = 1
      )

      abline(a = 0, b = 1, col = "red", lty = 2)
      dev.off()
      log_msg("Saved plot to ", plotfile)
    }

    return(list(
      r2 = r2,
      merged_means = merged_means,
      tsv_path = outfile
    ))
  }
}

compare_methylation_clusters_for_genes <- function(
  gene_names, sample_a_file, sample_b_file, sample_a_name, sample_b_name,
  coverage_threshold = 5, delta_threshold = 25, min_cluster_size = 5, no_plot = FALSE,
  fdr_threshold = 0.05, max_gap = 100,
  add_promotor = TRUE,
  batch_num = NULL,
  n_batches = NULL
) {
  if (!requireNamespace("grDevices", quietly = TRUE)) {
    stop("Package 'grDevices' required for adjustcolor")
  }

  all_clusters <- list()
  total_genes <- length(gene_names)

  for (i in seq_along(gene_names)) {
    gene_name <- gene_names[i]
    if (i %% ceiling(total_genes / 5) == 0 || i == total_genes) {
      pct <- round(i / total_genes * 100)
      bnum <- ifelse(is.null(batch_num), NA, batch_num)
      nbatches <- ifelse(is.null(n_batches), NA, n_batches)
      message(sprintf(
        "[%s] Progress: %d%% (%d of %d) — Gene: %s — Batch %d/%d",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        pct, i, total_genes, gene_name, bnum, nbatches
      ))
    }

    tx <- gene_metadata %>%
      filter(hgnc_symbol == gene_name) %>%
      filter(transcript_is_canonical == 1, chromosome_name %in% valid_chromosomes) %>%
      slice(1)

    if (nrow(tx)==0) next
    chr <- paste0("chr", tx$chromosome_name); start <- tx$start_position; end <- tx$end_position
    if (isTRUE(add_promotor)) {
      if (tx$strand == 1) {
        start <- max(1, start - 2000)
      } else {
        end <- end + 2000
      }
    }
    df_a <- extract_methyl_region(sample_a_file,chr,start,end,coverage_threshold)
    df_b <- extract_methyl_region(sample_b_file,chr,start,end,coverage_threshold)
    if (nrow(df_a)==0 || nrow(df_b)==0) next
    merged <- merge(df_a, df_b, by=c("chrom","start","end","strand"), suffixes=c("_a","_b")) %>%
      arrange(start) %>%
      mutate(
        delta_percent_modified = percent_modified_a - percent_modified_b,
        abs_delta = abs(delta_percent_modified)
      )

    if (nrow(merged) < 2) next  # ⛑️ Prevent crash on empty data

    # Cluster by distance: group CpGs that are within max_gap bp of each other
    merged$cluster_id <- NA_integer_
    cluster_id <- 1
    count <- 1
    merged$cluster_id[1] <- cluster_id
    for (j in 2:nrow(merged)) {
      if ((merged$start[j] - merged$start[j-1]) <= max_gap) {
        merged$cluster_id[j] <- cluster_id
        count <- count + 1
      } else {
        if (count < min_cluster_size) {
          merged$cluster_id[merged$cluster_id == cluster_id] <- NA_integer_
        }
        cluster_id <- cluster_id + 1
        merged$cluster_id[j] <- cluster_id
        count <- 1
      }
    }

    # Final cluster size check
    if (count < min_cluster_size) {
      merged$cluster_id[merged$cluster_id == cluster_id] <- NA_integer_
    }

    cluster_info <- merged %>% filter(!is.na(cluster_id))
    if (nrow(cluster_info)==0) next

    # Summarise clusters and apply Fisher's exact test
    cluster_regions <- cluster_info %>% group_by(cluster_id) %>% summarise(
      gene = gene_name,
      chrom = first(chrom),
      start = min(start), end = max(end), n_sites = n(),
      mean_percent_modified_a = mean(percent_modified_a),
      mean_percent_modified_b = mean(percent_modified_b),
      mean_delta = mean(delta_percent_modified),
      meth_a = sum(Nmod_a), unmeth_a = sum(Ncanonical_a),
      meth_b = sum(Nmod_b), unmeth_b = sum(Ncanonical_b), .groups="drop"
    ) %>% rowwise() %>% mutate(
      fisher_p = fisher.test(matrix(c(meth_a,unmeth_a,meth_b,unmeth_b), nrow=2))$p.value
    ) %>% ungroup()

    all_clusters[[gene_name]] <- cluster_regions
  }

  combined <- bind_rows(all_clusters)
  combined <- combined %>% select(-cluster_id)

  # Rename mean columns
  names(combined)[names(combined)=="mean_percent_modified_a"] <- paste0("mean_percent_modified_", sample_a_name)
  names(combined)[names(combined)=="mean_percent_modified_b"] <- paste0("mean_percent_modified_", sample_b_name)
  names(combined)[names(combined) == "meth_a"] <- paste0("meth_", sample_a_name)
  names(combined)[names(combined) == "unmeth_a"] <- paste0("unmeth_", sample_a_name)
  names(combined)[names(combined) == "meth_b"] <- paste0("meth_", sample_b_name)
  names(combined)[names(combined) == "unmeth_b"] <- paste0("unmeth_", sample_b_name)

return(combined)
}

# Initialize Ensembl BioMart (run this once)
mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# -----------------------------

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
message(sprintf("[%s] Running on %d cores", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), n_cores))

# -----------------------------
# 3. Prepare gene list
# -----------------------------
log_msg("Loading gene list...")

#gene_list <- NULL

if (exists("gene_list_path") && !is.null(gene_list_path) && gene_list_path != "") {
  if (gene_list_path == "all_protein_coding_genes") {
    message("Retrieving all protein-coding genes from Ensembl...")
    mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    gene_df <- biomaRt::getBM(
      attributes = c("hgnc_symbol", "chromosome_name"),
      filters = "biotype",
      values = "protein_coding",
      mart = mart
    )
    gene_list <- unique(gene_df$hgnc_symbol[
      gene_df$hgnc_symbol != "" &
        !(gene_df$chromosome_name %in% c("MT", "M", "Y"))
    ])
  } else {
    message("Reading gene list from file: ", gene_list_path)
    gene_list <- read.table(gene_list_path, header = FALSE, stringsAsFactors = FALSE)[[1]]
  }
} else {
  message("No gene list provided — skipping gene filtering.")
}

log_msg("Gene list loaded: ", length(gene_list), " genes")

valid_chromosomes <- c(as.character(1:22), "X", "Y")
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_metadata <- getBM(
   attributes = c("hgnc_symbol","chromosome_name","start_position","end_position","strand","transcript_is_canonical"),
   filters = "hgnc_symbol",
   values = gene_list,
   mart = mart
 ) %>% filter(transcript_is_canonical==1, chromosome_name %in% valid_chromosomes)

log_msg("Gene metadata loaded")

parallel_compare_methylation <- function(
    gene_names, sample_a_file, sample_b_file, sample_a_name, sample_b_name,
    coverage_threshold = 5, delta_threshold = 25, min_cluster_size = 10,
    fdr_threshold = 0.05, max_gap = 100, no_plot = TRUE,
    add_promotor = TRUE,
    n_cores = 20,
    log_file = file.path(outdir, "parallel_compare_methylation.log")
) {
  
  log_msg <- function(...) {
    msg <- paste0(Sys.time(), " - ", paste0(..., collapse = ""), "\n")
    cat(msg)
    cat(msg, file = log_file, append = TRUE)
  }
  
  log_msg("Starting parallel_compare_methylation with ", length(gene_names), " genes and ", n_cores, " cores.")
  
  gene_batches <- split(gene_names, cut(seq_along(gene_names), n_cores, labels = FALSE))
  
  results_list <- parallel::mclapply(seq_along(gene_batches), function(batch_num) {
    batch <- gene_batches[[batch_num]]
    batch_size <- length(batch)
    chunk_size <- ceiling(batch_size * 0.05)  # 5% chunk size
    chunks <- split(batch, ceiling(seq_along(batch) / chunk_size))
    
    batch_results <- list()
    genes_processed <- 0
    
    for (i in seq_along(chunks)) {
      chunk_genes <- chunks[[i]]
      
      # Call your function on the chunk subset of genes
      chunk_result <- compare_methylation_clusters_for_genes(
        gene_names = chunk_genes,
        sample_a_file = sample_a_file,
        sample_b_file = sample_b_file,
        sample_a_name = sample_a_name,
        sample_b_name = sample_b_name,
        coverage_threshold = coverage_threshold,
        delta_threshold = delta_threshold,
        min_cluster_size = min_cluster_size,
        fdr_threshold = fdr_threshold,
        max_gap = max_gap,
        add_promotor = add_promotor,
        no_plot = no_plot,
        batch_num = batch_num,
        n_batches = n_cores
      )
      
      batch_results[[i]] <- chunk_result
      genes_processed <- genes_processed + length(chunk_genes)
      
      progress_pct <- round(genes_processed / batch_size * 100)
      log_msg("[Batch ", batch_num, "] Processed ", genes_processed, "/", batch_size, " genes (", progress_pct, "%)")
    }
    
    log_msg("[Batch ", batch_num, "] Done processing all genes.")
    
    # Combine all chunk results from this batch (assuming data frames)
    dplyr::bind_rows(batch_results)
    
  }, mc.cores = n_cores)
  
  # Filter out NULL results (if any chunk or batch returned NULL)
  results_list <- Filter(Negate(is.null), results_list)
  
  all_combined <- dplyr::bind_rows(results_list)
  
  # Remove duplicate clusters (same chrom, start, end) - keep first gene encountered
  all_combined <- all_combined %>%
    dplyr::group_by(chrom, start, end) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  log_msg("Adjusting p-values for ", nrow(all_combined), " total clusters")
  
  all_combined$p_adj <- p.adjust(all_combined$fisher_p, method = "BH")
  
  all_significant <- all_combined %>%
    dplyr::filter(p_adj <= fdr_threshold & abs(mean_delta) >= delta_threshold)
  
  log_msg("Finished parallel_compare_methylation: ", nrow(all_significant), " significant clusters found.")
  
  return(list(
    all_clusters         = all_combined,
    significant_clusters = all_significant
  ))
}

#Main functions running..
linear <- methylation_linear_correlation(
   sample_a_file = sample_a_file,
   sample_b_file = sample_b_file,
   sample_a_name = sample_a_name,
   sample_b_name = sample_b_name,
   coverage_threshold = coverage_threshold,
   region_bed_file = region_bed_file,
   plot_main = TRUE,
   n_cores = n_cores
 )

res <- parallel_compare_methylation(
  gene_names = gene_list,
  sample_a_file = sample_a_file,
  sample_b_file = sample_b_file,
  sample_a_name = sample_a_name,
  sample_b_name = sample_b_name,
  coverage_threshold = coverage_threshold,
  delta_threshold = 25,
  min_cluster_size = 10,
  add_promotor = add_promotor,
  n_cores = n_cores
)

# Extract tables
result_all <- res$all_clusters
result_sig <- res$significant_clusters

outfile <- file.path(outdir, paste0("methylation_clusters_", sample_a_name, "_vs_", sample_b_name, ".tsv"))
sig_outfile <- file.path(outdir, paste0("sig_methylation_clusters_", sample_a_name, "_vs_", sample_b_name, ".tsv"))

write.table(
  result_all,
  file = outfile,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Only write significant results if any are present
if (nrow(result_sig) > 0) {
  write.table(
    result_sig,
    file = sig_outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
} else {
  message("No significant clusters found — skipping write of ", sig_outfile)
}

# ---- Volcano Plot ----

# Define plot paths
plotfile_annotated <- file.path(outdir, paste0("volcano_plot_", sample_a_name, "_vs_", sample_b_name, ".png"))
plotfile_unannotated <- file.path(outdir, paste0("volcano_plot_", sample_a_name, "_vs_", sample_b_name, "_unannotated.png"))

# Label top hits from all tested clusters
top_hits <- result_all %>%
  filter(p_adj < 1e-5 & abs(mean_delta) >= 25)

# Determine dynamic x-axis limits
delta_max <- if (nrow(result_all) > 0 && any(!is.na(result_all$mean_delta))) {
  max(abs(result_all$mean_delta), na.rm = TRUE)
} else {
  50
}
x_limit <- if (delta_max <= 50) c(-50, 50) else c(-1.05, 1.05) * delta_max

# --- Annotated Volcano Plot ---
png(plotfile_annotated, width = 2400, height = 2000, res = 300)

ggplot(result_all, aes(x = mean_delta, y = -log10(p_adj))) +
  geom_point(aes(color = p_adj < 0.05 & abs(mean_delta) >= 25), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) +
  geom_text_repel(
    data = top_hits,
    aes(label = gene),
    size = 2.5,
    max.overlaps = 100,
    box.padding = 0.2,
    point.padding = 0.2,
    segment.color = 'grey50'
  ) +
  geom_vline(xintercept = c(-25, 25), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +
  scale_x_continuous(limits = x_limit) +
  labs(
    title = paste("Volcano Plot:", sample_a_name, "vs", sample_b_name),
    x = "Mean Δ Methylation (%)",
    y = expression(-log[10](FDR)),
    color = "Significant"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )

dev.off()

# --- Unannotated Volcano Plot (No Gene Labels) ---
png(plotfile_unannotated, width = 2400, height = 2000, res = 300)

ggplot(result_all, aes(x = mean_delta, y = -log10(p_adj))) +
  geom_point(aes(color = p_adj < 0.05 & abs(mean_delta) >= 25), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick")) +
  geom_vline(xintercept = c(-25, 25), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +
  scale_x_continuous(limits = x_limit) +
  labs(
    title = paste("Volcano Plot:", sample_a_name, "vs", sample_b_name),
    x = "Mean Δ Methylation (%)",
    y = expression(-log[10](FDR)),
    color = "Significant"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 12),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9)
  )

dev.off()

# # Save summary result (number of clusters per gene)
# output_txt <- sprintf("cluster_summary_%s_vs_%s.txt", sample_a_name, sample_b_name)
# save_lines <- sapply(names(results_list), function(g) {
#   n <- tryCatch(nrow(results_list[[g]]), error = function(e) 0)
#   sprintf("%s\t%d", g, n)
# })
# writeLines(save_lines, con = output_txt)
# log_msg("Saved cluster result summary to ", output_txt)

# # -----------------------------
# # 5. Linear methylation comparison
# # -----------------------------
# plot_file <- sprintf("linear_comparison_%s_vs_%s.png", sample_a_name, sample_b_name)
# png(plot_file, width = 1000, height = 800)
# result_linear <- methylation_linear_correlation(
#   sample_a_file = sample_a_file,
#   sample_b_file = sample_b_file,
#   coverage_threshold = coverage_threshold,
#   sample_a_name = sample_a_name,
#   sample_b_name = sample_b_name,
#   region_bed_file = region_bed_file,
#   plot_main = TRUE
# )
# dev.off()
# log_msg("Saved linear methylation comparison plot to ", plot_file)

infile <- file.path(outdir, basename(sig_outfile))
heatmap_file <- file.path(outdir, "test_heatmap.png")

if (!file.exists(infile)) {
  warning("Significant output file not found, skipping heatmap and enrichment steps: ", infile)
} else {
  # Read file
  df <- read.delim(infile, header = TRUE, stringsAsFactors = FALSE)

# Filter for strong differential methylation (|delta| ≥ 25)
df_filtered <- df %>% filter(abs(mean_delta) >= 30)

# Create matrix: rows = genes, columns = samples (use dynamically generated column names)
mean_col_a <- paste0("mean_percent_modified_", sample_a_name)
mean_col_b <- paste0("mean_percent_modified_", sample_b_name)

mat <- df_filtered %>%
  select(all_of(c(mean_col_a, mean_col_b))) %>%
  as.matrix()

# Add gene names as row labels
rownames(mat) <- df_filtered$gene
colnames(mat) <- c(paste0("Sample_", sample_a_name), paste0("Sample_", sample_b_name))

# Transpose matrix: samples on x-axis
mat_t <- t(mat)

# Colour scale: 0 (hypomethylated) → 100 (hypermethylated)
heat_cols <- colorRampPalette(c("navy", "white", "firebrick3"))(101)

library(grid)
library(gtable)

# Set column names
colnames(mat) <- c("Proband", "Other") 

# Open PNG device with high resolution
png(heatmap_file, width = 5000, height = 6000, res = 300, type = "cairo")

# Create heatmap and capture the gtable object
p <- pheatmap(
  mat,
  color         = heat_cols,
  breaks        = seq(0, 100, length.out = length(heat_cols) + 1),
  cluster_rows  = TRUE,       # samples — cluster rows
  cluster_cols  = FALSE,      # genes — no clustering columns
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize      = 20,
  fontsize_col  = 30,         # gene labels font size
  fontsize_row  = 16,         # sample labels font size
  cellwidth     = 200,        # width of sample columns
  cellheight    = 20,         # height of gene rows
  border_color  = NA,
  main          = "Differential Methylation Heatmap (|Δ%| ≥ 30)",
  legend        = TRUE,
  legend.width  = unit(2, "cm"),   # wider color key bar
  legend.height = unit(18.75, "cm"), # ~1/4 of 7500 px height at 300 dpi (7500 px / 300 dpi = 25 inch; 1/4 * 25 in = 6.25 in = 18.75 cm)
  angle_col     = 45
)

# Increase right margin space for legend separation
p$gtable$widths[[ncol(p$gtable)]] <- unit(3, "cm")

# Identify legend grob’s row in gtable layout
legend_pos <- which(grepl("legend", p$gtable$layout$name))

# Get the current row of the legend grob (usually all legend parts share the same row)
legend_row <- unique(p$gtable$layout[legend_pos, "t"])

# Insert an empty row above the legend to push it down
p$gtable <- gtable_add_rows(p$gtable, unit(2, "cm"), pos = legend_row - 1)

# Draw modified heatmap with updated spacing
grid.newpage()
grid.draw(p$gtable)

# Close device
dev.off()


gene_list <- unique(df_filtered$gene)

# Step 1: Split significant genes by methylation direction
sig_up_genes <- df_filtered %>%
  filter(mean_delta > 0, p_adj < 0.05) %>%
  pull(gene) %>%
  unique()

sig_down_genes <- df_filtered %>%
  filter(mean_delta < 0, p_adj < 0.05) %>%
  pull(gene) %>%
  unique()

# Step 2: Run g:Profiler enrichment for each group
enrich_up <- gost(
  query = sig_up_genes,
  organism = "hsapiens",
  sources = c("GO:BP", "KEGG", "REAC"),
  correction_method = "fdr",
  significant = FALSE
)$result

enrich_down <- gost(
  query = sig_down_genes,
  organism = "hsapiens",
  sources = c("GO:BP", "KEGG", "REAC"),
  correction_method = "fdr",
  significant = FALSE
)$result

# Step 3: Plot top pathways for each group

## Upregulated methylation (Sample 1 > Sample 2)
top_up <- enrich_up %>%
  arrange(p_value) %>%
  slice_head(n = 50) %>%
  mutate(term_name = paste0(term_name, " (", source, ")"),
         term_name = factor(term_name, levels = rev(term_name)))

p_up <- ggplot(top_up, aes(x = -log10(p_value), y = term_name)) +
  geom_point(aes(size = intersection_size), color = "#D95F02") +
  labs(
    title = "Pathways: Hypermethylated in Sample 1",
    x = "-log10(FDR-adjusted p-value)", y = NULL,
    size = "Gene Count"
  ) +
  theme_minimal(base_size = 12)

## Downregulated methylation (Sample 1 < Sample 2)
top_down <- enrich_down %>%
  arrange(p_value) %>%
  slice_head(n = 50) %>%
  mutate(term_name = paste0(term_name, " (", source, ")"),
         term_name = factor(term_name, levels = rev(term_name)))

p_down <- ggplot(top_down, aes(x = -log10(p_value), y = term_name)) +
  geom_point(aes(size = intersection_size), color = "#1B9E77") +
  labs(
    title = "Pathways: Hypomethylated in Sample 1",
    x = "-log10(FDR-adjusted p-value)", y = NULL,
    size = "Gene Count"
  ) +
  theme_minimal(base_size = 12)

# Display the plots (or use ggsave to export them)
print(p_up)
print(p_down)
}
dev.off()

# Input file (use dynamic sample names - this is a fallback block for manual testing)
infile <- file.path(outdir, paste0("sig_methylation_clusters_", sample_a_name, "_vs_", sample_b_name, ".tsv"))

# Read data
df <- read.delim(infile, header = TRUE, stringsAsFactors = FALSE)

# Filter for strong differential methylation (|mean_delta| >= 30)
df_filtered <- df %>% filter(abs(mean_delta) >= 30)

# Genes with higher methylation in Sample1 (positive delta)
genes_delta_pos <- unique(df_filtered$gene[df_filtered$mean_delta > 0])

# Genes with higher methylation in Sample2 (negative delta)
genes_delta_neg <- unique(df_filtered$gene[df_filtered$mean_delta < 0])

cat("Genes with higher methylation in Sample1:", length(genes_delta_pos), "\n")
cat("Genes with higher methylation in Sample2:", length(genes_delta_neg), "\n")

# Run enrichment analyses, include all results
res_pos <- gost(query = genes_delta_pos,
                organism = "hsapiens",
                sources = c("GO:BP", "KEGG", "REAC"),
                correction_method = "fdr",
                significant = FALSE)

res_neg <- gost(query = genes_delta_neg,
                organism = "hsapiens",
                sources = c("GO:BP", "KEGG", "REAC"),
                correction_method = "fdr",
                significant = FALSE)


# Extract enrichment results without 'intersection'
if (!is.null(res_pos$result) && nrow(res_pos$result) > 0) {
  df_pos <- dplyr::select(res_pos$result, term_id, term_name, source, p_value) 
  df_pos$sample <- "Higher methylation in Sample1"
} else {
  df_pos <- data.frame()
}

if (!is.null(res_neg$result) && nrow(res_neg$result) > 0) {
  df_neg <- dplyr::select(res_neg$result, term_id, term_name, source, p_value)
  df_neg$sample <- "Higher methylation in Sample2"
} else {
  df_neg <- data.frame()
}

# Rename p_value columns before join
colnames(df_pos)[colnames(df_pos) == "p_value"] <- "p_value_pos"
colnames(df_neg)[colnames(df_neg) == "p_value"] <- "p_value_neg"

# Join enrichment results
combined <- full_join(df_pos, df_neg, by = c("term_id", "term_name", "source"))

# Prepare data for plotting: pivot longer and filter out NA p_values
plot_df <- combined %>%
  pivot_longer(cols = c(p_value_pos, p_value_neg), names_to = "sample_pval", values_to = "p_value") %>%
  mutate(sample = ifelse(sample_pval == "p_value_pos", "Higher methylation in Sample1", "Higher methylation in Sample2")) %>%
  filter(!is.na(p_value))

# Calculate -log10 p-value for X-axis and circle size
plot_df$logp <- -log10(plot_df$p_value)

# Cap logp values for visualization
plot_df$logp <- ifelse(plot_df$logp > 20, 20, plot_df$logp)

# Select top 20 pathways by minimum p-value across samples
top_terms <- combined %>%
  mutate(min_p = pmin(p_value_pos, p_value_neg, na.rm = TRUE)) %>%
  arrange(min_p) %>%
  slice_head(n = 20) %>%
  pull(term_id)

plot_df <- plot_df %>% filter(term_id %in% top_terms)

# Plot
ggplot(plot_df, aes(x = logp, y = reorder(term_name, logp),
                    color = sample, size = logp)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Higher methylation in Sample1" = "blue", "Higher methylation in Sample2" = "red")) +
  labs(x = "-log10(adjusted p-value)",
       y = "Pathway",
       color = "Methylation direction",
       size = "Enrichment significance") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9))


names(res_pos$result)
names(res_neg$result)


####

# File path (use dynamic sample names)
infile <- file.path(outdir, paste0("sig_methylation_clusters_", sample_a_name, "_vs_", sample_b_name, ".tsv"))

# Read data
df <- read.delim(infile, header = TRUE, stringsAsFactors = FALSE)

# Filter for strong differential methylation (|mean_delta| >= 30)
#df_filtered <- df %>% filter(abs(mean_delta) >= 30)
df_filtered <- df

# Genes with higher methylation in Sample1 (positive delta)
genes_delta_pos <- unique(df_filtered$gene[df_filtered$mean_delta > 0])

# Genes with higher methylation in Sample2 (negative delta)
genes_delta_neg <- unique(df_filtered$gene[df_filtered$mean_delta < 0])

cat("Genes with higher methylation in Sample1:", length(genes_delta_pos), "\n")
cat("Genes with higher methylation in Sample2:", length(genes_delta_neg), "\n")

# Run enrichment analyses, include all results (significant = FALSE)
res_pos <- gost(query = genes_delta_pos,
                organism = "hsapiens",
                sources = c("GO:BP", "KEGG", "REAC"),
                correction_method = "fdr",
                significant = FALSE)

res_neg <- gost(query = genes_delta_neg,
                organism = "hsapiens",
                sources = c("GO:BP", "KEGG", "REAC"),
                correction_method = "fdr",
                significant = FALSE)

# Extract enrichment results explicitly (pipe-free)
if (!is.null(res_pos$result) && nrow(res_pos$result) > 0) {
    df_pos <- res_pos$result %>%
      dplyr::select(term_id, term_name, source, p_value, dplyr::any_of("intersection"))
    if (!"intersection" %in% colnames(df_pos)) df_pos$intersection <- NA_character_
    df_pos$sample <- "Higher methylation in Sample1"
  } else {
    df_pos <- data.frame(term_id=character(), term_name=character(), source=character(), p_value=numeric(), intersection=character(), sample=character(), stringsAsFactors = FALSE)
  }

if (!is.null(res_neg$result) && nrow(res_neg$result) > 0) {
    df_neg <- res_neg$result %>%
      dplyr::select(term_id, term_name, source, p_value, dplyr::any_of("intersection"))
    if (!"intersection" %in% colnames(df_neg)) df_neg$intersection <- NA_character_
    df_neg$sample <- "Higher methylation in Sample2"
  } else {
    df_neg <- data.frame(term_id=character(), term_name=character(), source=character(), p_value=numeric(), intersection=character(), sample=character(), stringsAsFactors = FALSE)
}
    
# Function to compute mean_delta per pathway from intersected genes
compute_mean_delta <- function(df_enrich, df_input) {
  if (nrow(df_enrich) == 0) return(df_enrich)

  df_enrich$mean_delta_pathway <- NA_real_

  for (i in seq_len(nrow(df_enrich))) {
    intersection_value <- df_enrich$intersection[i]
    if (is.na(intersection_value) || intersection_value == "") next

    genes_in_pathway <- unlist(strsplit(intersection_value, ","))
    genes_in_pathway <- trimws(genes_in_pathway)
    if (length(genes_in_pathway) == 0) next

    deltas <- df_input$mean_delta[df_input$gene %in% genes_in_pathway]
    if (length(deltas) > 0) {
      df_enrich$mean_delta_pathway[i] <- mean(deltas, na.rm=TRUE)
    }
  }
  df_enrich
}

# Compute mean_delta_pathway for pos and neg
df_pos <- compute_mean_delta(df_pos, df_filtered)
df_neg <- compute_mean_delta(df_neg, df_filtered)

# Rename mean_delta_pathway columns before joining
colnames(df_pos)[colnames(df_pos) == "mean_delta_pathway"] <- "mean_delta_pathway_pos"
colnames(df_neg)[colnames(df_neg) == "mean_delta_pathway"] <- "mean_delta_pathway_neg"

# Join enrichment results
combined <- full_join(df_pos, df_neg, by = c("term_id", "term_name", "source"), suffix = c("_pos", "_neg"))

# Prepare data for plotting: pivot longer and filter out NA p_values
plot_df <- combined %>%
  pivot_longer(cols = c(p_value_pos, p_value_neg), names_to = "sample_pval", values_to = "p_value") %>%
  mutate(sample = ifelse(sample_pval == "p_value_pos", "Higher methylation in Sample1", "Higher methylation in Sample2")) %>%
  filter(!is.na(p_value))

# Initialize mean_delta_pathway column
plot_df$mean_delta_pathway <- NA_real_

# Assign mean_delta_pathway from combined results
plot_df$mean_delta_pathway[plot_df$sample == "Higher methylation in Sample1"] <-
  combined$mean_delta_pathway_pos[match(plot_df$term_id[plot_df$sample == "Higher methylation in Sample1"], combined$term_id)]

plot_df$mean_delta_pathway[plot_df$sample == "Higher methylation in Sample2"] <-
  combined$mean_delta_pathway_neg[match(plot_df$term_id[plot_df$sample == "Higher methylation in Sample2"], combined$term_id)]

# Cap mean_delta_pathway values for visualization
plot_df$mean_delta_pathway <- ifelse(plot_df$mean_delta_pathway > 50, 50,
                                     ifelse(plot_df$mean_delta_pathway < -50, -50, plot_df$mean_delta_pathway))

# Calculate -log10 p-value and cap max
plot_df$logp <- -log10(plot_df$p_value)
plot_df$logp <- ifelse(plot_df$logp > 20, 20, plot_df$logp)

# Select top 20 pathways by minimum p-value across samples
top_terms <- combined %>%
  mutate(min_p = pmin(p_value_pos, p_value_neg, na.rm = TRUE)) %>%
  arrange(min_p) %>%
  slice_head(n = 20) %>%
  pull(term_id)

plot_df <- plot_df %>% filter(term_id %in% top_terms)

# Plot
dev.off()

ggplot(plot_df, aes(x = mean_delta_pathway, y = reorder(term_name, abs(mean_delta_pathway)),
                    color = sample, size = logp)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Higher methylation in Sample1" = "blue", "Higher methylation in Sample2" = "red")) +
  labs(x = "Average methylation difference (%) in pathway genes",
       y = "Pathway",
       color = "Methylation direction",
       size = "-log10(adj p-value)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 9))
}