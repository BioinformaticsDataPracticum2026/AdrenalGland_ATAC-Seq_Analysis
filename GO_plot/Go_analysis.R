# load library
library(rGREAT)
library(GenomicRanges)

# convert narrowPeak into GRanges
shared_region <- read.table("../shared_peaks_conservative.narrowPeak", header = FALSE)
human_specific <- read.table("../human_specific_peaks_conservative.narrowPeak", header = FALSE)
mouse_specific <- read.table("../mouse_specific_peaks_conservative.narrowPeak", header = FALSE)

# just need the first three cols in the narrowPeak
shared <- GRanges(
  seqnames = shared_region$V1,
  ranges   = IRanges(start = shared_region$V2 + 1, end = shared_region$V3)
)

human <- GRanges(
  seqnames = human_specific$V1,
  ranges   = IRanges(start = human_specific$V2 + 1, end = human_specific$V3)
)

mouse <- GRanges(
  seqnames = mouse_specific$V1,
  ranges   = IRanges(start = mouse_specific$V2 + 1, end = mouse_specific$V3)
)

# Use human gene annotation database
job1 <- great(shared, "BP", "TxDb.Hsapiens.UCSC.hg19.knownGene")
job2 <- great(human, "BP", "TxDb.Hsapiens.UCSC.hg19.knownGene")
job3 <- great(mouse, "BP", "TxDb.Hsapiens.UCSC.hg19.knownGene")

# get the enrichment table
go_bp1 <- getEnrichmentTables(job1)
head(go_bp1)

go_bp2 <- getEnrichmentTables(job2)
head(go_bp2)

go_bp3 <- getEnrichmentTables(job3)
head(go_bp3)

# load library for ploting
library(ggplot2)

# for shared region
top1 <- head(go_bp1, 10)  # top 10 terms
# give biological process name
top1$Term <- substr(top1$description, 1, 50)

# plot
ggplot(top1, aes(x = reorder(Term, fold_enrichment), y = fold_enrichment)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  xlab("") +
  ylab("Fold Enrichment") +
  ggtitle("GO BP Enrichment - Shared_peaks") +
  theme_bw()

ggsave("GO_Job1.png", width = 8, height = 6)

# for human specific
top2 <- head(go_bp2, 10)
top2$Term <- substr(top2$description, 1, 50)


ggplot(top2, aes(x = reorder(Term, fold_enrichment), y = fold_enrichment)) +
  geom_bar(stat = "identity", fill = "forestgreen") +
  coord_flip() +
  xlab("") +
  ylab("Fold Enrichment") +
  ggtitle("GO BP Enrichment - Human-specific") +
  theme_bw()
ggsave("GO_Job2.png", width = 8, height = 6)

# for mouse specific
top3 <- head(go_bp3, 10)
top3$Term <- substr(top3$description, 1, 50)

ggplot(top3, aes(x = reorder(Term, fold_enrichment), y = fold_enrichment)) +
  geom_bar(stat = "identity", fill = "tomato") +
  coord_flip() +
  xlab("") +
  ylab("Fold Enrichment") +
  ggtitle("GO BP Enrichment - Mouse-specific") +
  theme_bw()

ggsave("GO_Job3.png", width = 8, height = 6)

# bubble plot for shared region
ggplot(top1, aes(x = fold_enrichment, 
                 y = reorder(Term, fold_enrichment), 
                 size = observed_region_hits,        # size represents number of regions
                 color = -log10(p_adjust))) +       # color represents significance
  geom_point() +
  scale_color_gradient(low = "lightblue", high = "red") +
  xlab("Fold Enrichment") +
  ylab("") +
  ggtitle("Top 10 GO_BP terms - Shared peaks") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))
ggsave("bubble_Job1.png", width = 8, height = 6)

# bubble plot for human specific
ggplot(top2, aes(x = fold_enrichment, 
                 y = reorder(Term, fold_enrichment), 
                 size = observed_region_hits,        # size represents number of regions
                 color = -log10(p_adjust))) +       # color represents significance
  geom_point() +
  scale_color_gradient(low = "lightblue", high = "red") +
  xlab("Fold Enrichment") +
  ylab("") +
  ggtitle("Top 10 GO_BP terms - Human-sepcific") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))
ggsave("bubble_Job2.png", width = 8, height = 6)


# bubble plot for mouse specific
ggplot(top3, aes(x = fold_enrichment, 
                 y = reorder(Term, fold_enrichment), 
                 size = observed_region_hits,        # size represents number of regions
                 color = -log10(p_adjust))) +       # color represents significance
  geom_point() +
  scale_color_gradient(low = "lightblue", high = "red") +
  xlab("Fold Enrichment") +
  ylab("") +
  ggtitle("Top 10 GO_BP terms - Mouse-specific") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))
ggsave("bubble_Job3.png", width = 8, height = 6)
