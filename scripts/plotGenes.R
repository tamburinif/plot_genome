# input files
# gff_file <- snakemake@input[[1]]
# names_file <- snakemake@input[[2]]

library(gggenes)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(dplyr)

# gggenes documentation
# https://github.com/wilkox/gggenes

# input gff file
gff_file <- snakemake@input[[1]]

# read input gff
gff <- read.gff(gff_file)
colnames(gff)[9] <- "gene"

# read names file
# names <- read.table(names_file, sep = '\t', header = TRUE)

# rename samples
# gff$seqid <- names[match(gff$seqid, names$Sample), "Label"]

# add length column
gff$length <- (gff$end - gff$start)

## Gene reordering breaks if a cluster isn't present in every genome and if genes / samples are not sorted

# sort gff file
gff <- gff[order(gff$seqid, gff$start),]

## anchor with respect to longest gene - find the start of the longest gene which will be the new zero

# find longest gene
longest <- gff[rev(order(gff$length)),][1,"gene"]

anchors <- filter(gff, gene == longest)
anchors <- aggregate(start ~ seqid, data = anchors, min)
ends <- aggregate(end ~ seqid, data = gff, max)

# fix strand
gff$strand <- ifelse(gff$strand == '+', 1, 0)

## find core genome -- genes present in > p fraction of genomes

# minimum fraction of genomes
p <- 0.5

# number of samples
n_samp <- length(unique(gff$seqid))

# count unique genes
unique_genes <- unique(gff[, c("seqid", "gene")])
counts <- aggregate(seqid ~ gene, data = unique_genes, length)

# find core and non-core genes
non_core <- unique(filter(counts, seqid <= p*n_samp)$gene)
core <- unique(filter(counts, seqid > p*n_samp)$gene)

# color core genome
myCols <- colorRampPalette(brewer.pal(11, "Spectral"))
myPal <- myCols(length(core))
myPal <- sample(myPal)
core_df <- data.frame("gene" = core, "color" = myPal)
non_core_df <- data.frame("gene" = non_core, "color" = "gray90")
colors <- rbind(core_df, non_core_df)

# add color column
gff$color <- colors[match(gff$gene, colors$gene), "color"]

## sort!!
# gff$seqid <- factor(gff$seqid, levels = rev(sort(gff$seqid)))
# gff <- gff[order(gff$gene),]

df <- gff[, c(9, 11)]
df$gene_num <- gsub("Cluster ", '', df$gene)
colors <- unique(df[order(df$gene_num), ])$color

# set sample ordering
# gff$seqid <- factor(gff$seqid, levels = rev(names$Label))

gff$gene_num <- gsub("Cluster ", '', gff$gene)

# plot
ggplot2::ggplot(gff, ggplot2::aes(xmin = start, xmax = end, y = seqid, fill = gene, forward = strand, label = gene_num)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  geom_gene_label(align = "left") +
  # ggplot2::facet_wrap(~ seqid, scales = "free", ncol = 1) +
  # ggplot2::scale_fill_manual(values = myPal, guide = guide_legend(ncol = 4)) +
  # scale_fill_manual(values = as.character(gff$color)) +
  ggplot2::scale_fill_manual(values = as.character(colors)) +
  theme_bw() +
  theme(
    # title = element_text(size = 14),
    strip.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    legend.position = "none"
  ) +
  ylab("Genome") +
  xlab("Kilobases")

# save plot
ggsave(snakemake@input[[1]], device = "pdf", width = 16, height = 12, units = "in", useDingbats=FALSE)

# average number of genes
# mean(count(gff, "seqid")$freq)


