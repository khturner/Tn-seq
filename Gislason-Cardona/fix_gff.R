library(tidyverse)
library(seqinr)
library(fuzzyjoin)

reference_fasta <- "k56.fasta"
reference_gff <- "k56.gff"

reference_genome <- read.fasta(reference_fasta)
template_sizes <- data.frame(template = names(reference_genome),
                             template_length = unlist(lapply(reference_genome, length)),
                             stringsAsFactors = F) %>% tbl_df

genome_features <- read_tsv(reference_gff, comment = "#",
                            col_names = c("seqname", "source", "feature", "start", "end", "score", "strand",
                                          "frame", "attribute"))

renamed_genomed_features <- template_sizes %>%
  mutate(end = cumsum(template_length), start = 1 + end - template_length) %>%
  select(template, start, end) %>%
  # interval_inner_join(sample_n(genome_features, 100)) %>%
  interval_inner_join(genome_features) %>%
  mutate(adjstart = start.y - start.x + 1, adjend = end.y - start.x + 1) %>%
  select(template, source, feature, start = adjstart, end = adjend, score, strand, frame, attribute)

system("grep '^#' k56.gff > k56.bycontig.gff", T)
write_tsv(renamed_genomed_features, "k56.bycontig.gff", append = T)
