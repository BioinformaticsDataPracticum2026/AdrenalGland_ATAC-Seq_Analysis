if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("rGREAT", quietly = TRUE))
  BiocManager::install("rGREAT")
if (!requireNamespace("ChIPpeakAnno", quietly = TRUE))
  BiocManager::install("ChIPpeakAnno")

library(rGREAT)
library(ChIPpeakAnno)
library(tidyverse)

go_analysis <- function(peaks_file, species, conservation) {
  species <- match.arg(species, choices = c("mouse", "human"))
  conservation <- match.arg(conservation, choices = c("conserved", "specific"))

  gr <- toGRanges(peaks_file, format = "BED")

  if (species == "human") {
    ref <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
  } else if (species == "mouse") {
    ref <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
  }

  res <- great(gr, "GO:BP", ref)
  tb <- getEnrichmentTable(res)

  filename <- paste(species, conservation, "great_results.tsv", sep = "_")
  write.table(tb, file = filename, sep = "\t")

  return(filename)
}

plot_go <- function(file) {
  res <- read.table(file, header = TRUE, sep = "\t")
  
  top_terms <- res %>% arrange(p_adjust) %>% slice_head(n = 20)

  filename <- paste(tools::file_path_sans_ext(file))

  ggplot(top_terms, aes(x = -log10(p_adjust), y = reorder(description, -log10(p_adjust)))) + geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "-log10(adjusted p-value)", y = "GO Biological Process", title = filename) + theme_bw()

  ggsave(file.path("rGREAT/plots", paste(filename, ".png")))
}


liftover_files <- list.files(
  path = "liftOver_pipeline/peak_sets",
  pattern = "\\.bed$",
  recursive = TRUE,
  full.names = TRUE
)
halper_files <- list.files(
  path = "HALPER_pipeline/mapping/results",
  pattern = "\\.narrowPeak$",
  recursive = TRUE,
  full.names = TRUE
)

print(length(liftover_files))
print(length(halper_files))
all_files <- c(liftover_files, halper_files)

# loop through all_files, extract species/conservation status for inputs, run go_analysis and plot_go
