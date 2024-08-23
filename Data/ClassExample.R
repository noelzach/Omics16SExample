library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(Biostrings)
library(ggpubr)
library(decontam)
library(metagenomeSeq)
library(indicspecies)

packageVersion("phyloseq")
packageVersion("decontam")
packageVersion("metagenomeSeq")
packageVersion("decontam")
packageVersion("tidyverse")


install.packages("ggplot2")
install_version("ggplot2", version = "0.9.1", repos = "http://cran.us.r-project.org")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
BiocManager::install("metagenomeSeq")
BiocManager::install("decontam")
BiocManager::install("Biostrings")

# color blind pallet
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# load in our files from the HPC for input into phyloseq. 

# Taxonomy 
tax.path <- file.choose() # "/Users/zacharynoel/Library/CloudStorage/Box-Box/Auburn/Teaching/Omics in Agriculture/Lab 10 - Microbiome 2/16s_otus_sintax_tax_rdp.csv"
tax <- read.delim(tax.path, header = T, row.names = 1, sep = ",")
head(tax)
tax$OTU <- rownames(tax)
TAX.bacteria <- phyloseq::tax_table(as.matrix(tax))

#OTU table 
table.path <- file.choose() # "/Users/zacharynoel/Library/CloudStorage/Box-Box/Auburn/Teaching/Omics in Agriculture/Lab 10 - Microbiome 2/otu_table_16s_UPARSE.csv"
table <- read.csv(table.path)
rownames(table) <- table$OTU
table <- table[,-1]
OTU.bacteria <- phyloseq::otu_table(table, taxa_are_rows = TRUE)

#Meta data 
samples.path <- file.choose() # "/Users/zacharynoel/Library/CloudStorage/Box-Box/Auburn/Teaching/Omics in Agriculture/Lab 10 - Microbiome 2/soil_prok_map.csv"
samples <- read.csv(samples.path)
rownames(samples) <- samples$SampleID #row names must match OTU table headers
SAMP.bacteria <- phyloseq::sample_data(samples)

#fasta file
fasta.path <- file.choose() # "/Users/zacharynoel/Library/CloudStorage/Box-Box/Auburn/Teaching/Omics in Agriculture/Lab 10 - Microbiome 2/otus_16s.fasta"
FASTA.bacteria <- Biostrings::readDNAStringSet(fasta.path, format="fasta", seek.first.rec=TRUE, use.names=TRUE)

#Phylogeny
tree.path <- file.choose() # "/Users/zacharynoel/Library/CloudStorage/Box-Box/Auburn/Teaching/Omics in Agriculture/Lab 10 - Microbiome 2/otus_phylogeny.nwk"
tree <- phyloseq::read_tree(tree.path)

phyloseq.start <- phyloseq(SAMP.bacteria, TAX.bacteria, OTU.bacteria, FASTA.bacteria, tree)

phyloseq.start@otu_table # the OTU table
phyloseq.start@tax_table # the taxonomy table
phyloseq.start@sam_data # the metadata
phyloseq.start@refseq # the sequences 
phyloseq.start@phy_tree # the tree

# removing chloroplast or taxa not assigned at the domain level
physeq.no.chloro <- phyloseq.start %>% subset_taxa(Class!= "Chloroplast" & Domain!= "unidentified")

## DECONTAMINATE
#Use the full dataset to call contaminants, then remove them, if they exist in the non plant OTU dataset
sample_data(physeq.no.chloro)$is.neg <- sample_data(physeq.no.chloro)$SampleorControl == "Control Sample"
contamdf.prev <- isContaminant(physeq.no.chloro, method="prevalence", neg="is.neg", threshold = 0.1, normalize = TRUE)
badTaxa <- rownames(contamdf.prev[contamdf.prev$contaminant == TRUE,])

print(badTaxa)

# transform data to presence absence
ps.pa <- transform_sample_counts(physeq.no.chloro, function(abund) 1*(abund>0))

# making a dataframe for both negative and positive smaples.
ps.pa.neg <- prune_samples(sample_data(ps.pa)$SampleorControl == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$SampleorControl == "True Sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
decontaminated <- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + 
  geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") + 
  theme_classic() + 
  scale_color_manual(values = c(cbbPalette[[1]], cbbPalette[[2]]))

#Take out the contaminants 
goodTaxa <- setdiff(taxa_names(physeq.no.chloro), badTaxa)
physeq.clean <- prune_taxa(goodTaxa, physeq.no.chloro)

# Now that we have removed the contaminants using the negative controls lets get rid of those samples
physeq.clean.samples <- physeq.clean %>% 
  subset_samples(SampleorControl == "True Sample") %>%
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

# Now lets look at the read distribution per sample and decide if we need to get rid of some samples because of low sequence depth
# a good general rule of thumb is samples below 1000 reads could be eliminated, although this isn't a hard rule, you can remove at 10,000 or more if you want
sort(sample_sums(physeq.clean.samples), decreasing = T) # read distribution

# how many total reads are we now working with?
sum(sample_sums(physeq.clean.samples))

# What is our mean and median read depth per sample? 
mean(sample_sums(physeq.clean.samples))

median(sample_sums(physeq.clean.samples))

#Lets make a histogram of the read distribution and put the median read depth
read.depths <- data.frame(sample_sums(physeq.clean.samples))
colnames(read.depths) <- "read.depth"
read.depth.plot <- ggplot(read.depths, aes(read.depth)) +
  geom_histogram(fill = cbbPalette[[3]], color = "black") + 
  geom_vline(xintercept = median(sample_sums(physeq.clean.samples)), linetype = "dashed") + 
  theme_classic() + 
  xlab("Read Depth")


# Rarefaction analysis 
sam.data <- data.frame(physeq.clean.samples@sam_data)
bOTU.table <- otu_table(physeq.clean.samples) %>%
  as.data.frame() %>%
  as.matrix()

raremax <- min(rowSums(t(bOTU.table)))
rare.fun <- rarecurve(t(bOTU.table), step = 1000, sample = raremax, tidy = T)

bac.rare.curve.extract2 <- left_join(sam.data, rare.fun, by = c("SampleID" = "Site"))

bac.rare <- ggplot(bac.rare.curve.extract2, aes(x = Sample, y = Species, group = SampleID, color = Management)) + 
  #geom_point() +
  geom_line() + 
  xlab("Reads") + 
  ylab("Number of OTUs") +
  ggtitle("Prokaryote") +
  theme_classic() + 
  geom_vline(xintercept = median(sample_sums(physeq.clean.samples)), linetype = "dashed") +
  scale_color_manual(values = cbbPalette)


# Stick many figures together for a publication ready figure
SuplementalFig1 <- ggarrange(bac.rare, read.depth.plot, decontaminated, nrow = 1, labels = c("a", "b", "c"))
ggsave("SupplementalFig1.pdf", dpi=300, width = 40, height = 10, units = "cm")

# Normalize Sampling reads based on proportions
physeq.prop <- transform_sample_counts(physeq.clean.samples, function(x) x/sum(x)) # normalize using proportions

# Normalize Sampling reads based on cumulative sum scaling (CSS normalization)
MGS <- phyloseq_to_metagenomeSeq(physeq.clean.samples)
p <- metagenomeSeq::cumNormStatFast(MGS)
MGS <- metagenomeSeq::cumNorm(MGS, p =p)
metagenomeSeq::normFactors(MGS) # exports the normalized factors for each sample
norm.bacteria <- metagenomeSeq::MRcounts(MGS, norm = T)
norm.bacteria.OTU <- phyloseq::otu_table(norm.bacteria, taxa_are_rows = TRUE)

physeq.css <- phyloseq::phyloseq(norm.bacteria.OTU, SAMP.bacteria, TAX.bacteria, FASTA.bacteria, tree)

# Save RDS files of each type of normalized reads for easy loading in the future and reproducibility 
# Save an object to a file
saveRDS(physeq.css, file = "Bacteria_Soil_CSSnorm_102021.rds")
saveRDS(physeq.clean.samples, file = "Bacteria_Soil_nonnorm_102021.rds")

# This is how to read in an RDS file
# just uncomment the line below to read it
physeq.css <- readRDS(file = "Bacteria_Soil_CSSnorm_102021.rds")
physeq.clean.samples <- readRDS(file = "Bacteria_Soil_nonnorm_102021.rds")


# Getting comfortable 
# most abundant OTUs
most.abundant <- data.frame(sort(rowSums(physeq.clean.samples@otu_table), decreasing = TRUE))
head(most.abundant)
physeq.clean.samples@refseq$BOTU_2
tax[tax$OTU == "BOTU_2",]

# Occupancy (i.e., number of samples the bacteria were observed)
presence.absence <- transform_sample_counts(physeq.clean.samples, function(abund) 1*(abund>0))
presence.absence@otu_table
OTU.occupancy <- data.frame(sort(rowSums(presence.absence@otu_table)))
colnames(OTU.occupancy) <- "Occupancy"
OTU.occupancy$OTU <- rownames(OTU.occupancy)
OTU.occupancy$percent.occupancy <- (100*OTU.occupancy$Occupancy/max(OTU.occupancy$Occupancy))

physeq.clean.samples@otu_table

# ALPHA DIVERSITY
physeq.clean.samples@sam_data$shannon <- estimate_richness(physeq.clean.samples, measures=c("Shannon"))$Shannon
physeq.clean.samples@sam_data$invsimpson <- estimate_richness(physeq.clean.samples, measures=c("InvSimpson"))$InvSimpson
physeq.clean.samples@sam_data$richness <- estimate_richness(physeq.clean.samples, measures=c("Observed"))$Observed
physeq.clean.samples@sam_data$even <- physeq.clean.samples@sam_data$shannon/log(physeq.clean.samples@sam_data$richness)

sample.data.fungi <- data.frame(physeq.clean.samples@sam_data)

# Richness over time
richness.time <- ggplot(sample.data.fungi, aes(x = Collection, y = richness)) + 
  geom_boxplot() +
  geom_jitter() + 
  ylab("Richness") +
  stat_compare_means(method = "anova") + 
  xlab("")+
  theme_classic() 

# Richness by treatment
richness.management <- ggplot(sample.data.fungi, aes(x = Management, y = richness)) + 
  geom_boxplot() +
  geom_jitter() + 
  ylab("Richness") + 
  stat_compare_means(method = "anova") + 
  xlab("")+
  theme_classic() 

# Richness by treatment*collection interaction
richness.management.time <- ggplot(sample.data.fungi, aes(x = Collection, y = richness)) + 
  geom_boxplot() +
  geom_jitter() + 
  ylab("Richness") + 
  stat_compare_means(method = "anova") + 
  xlab("")+
  theme_classic() +
  facet_wrap(~Management)

# Shannon diversity by management
shannon.management <- ggplot(sample.data.fungi, aes(x = Management, y = shannon)) + 
  geom_boxplot() +
  geom_jitter() + 
  ylab("Shannon") + 
  stat_compare_means(method = "anova") + 
  xlab("")+
  theme_classic() 

# BETA DIVERSITY
# Principle coordinates analysis with Bray-Curtis distances
ordination.pcoa <- ordinate(physeq.css, "PCoA", "bray") # calculate the resemblance and ordinate using PCoA
ordination.pcoa$vectors # positions of your points on the pcoa graph
ordination.pcoa$values #values to calculate the variance explained on each axis (dimension)

pcoa <- plot_ordination(physeq.css, ordination = ordination.pcoa, type = "samples", color = "Management", shape = "Collection") +
  theme_classic() + 
  scale_color_manual(values = cbbPalette)
pcoa

pcoa.data <- pcoa$data # taking the data to make a fancy plot

ggplot() + 
  geom_point(data = pcoa.data, aes(x = Axis.1, y = Axis.2, shape = Collection, fill = Management), alpha = 0.8, size = 2) +
  theme_bw() +
  ylab("PCoA2 (13.7%)") + 
  xlab("PCoA1 (18.0%)") +
  scale_fill_manual(values=cbbPalette) +
  #stat_ellipse(data = global.pcoa.obj1.data, aes(x = Axis.1, y = Axis.2, group = Compartment), type = "norm", linetype = 2) +
  scale_shape_manual(values=c(21, 24, 23)) +
  guides(fill=guide_legend(override.aes=list(shape=21))) 

# Principle coordinates analysis with Weighted-unifrac distances
physeq.css@phy_tree
unifrac <- UniFrac(physeq.css, weighted = FALSE)
ordination.unifrac.pcoa <- ordinate(physeq.css, "PCoA", distance = unifrac) # calculate the resemblance and ordinate using PCoA
ordination.unifrac.pcoa$vectors # positions of your points on the pcoa graph
ordination.unifrac.pcoa$values #values to calculate the variance explained on each axis (dimension)

pcoa.unifrac <- plot_ordination(physeq.css, ordination = ordination.unifrac.pcoa, type = "samples", color = "Management") +
  theme_classic() + 
  scale_color_manual(values = cbbPalette) + 
  stat_ellipse()
pcoa.unifrac

# PERMANOVA - testing for differences in centroids
prok.dist.bray = phyloseq::distance(physeq.css, "bray") # create bray-curtis distance matrix
adonis2(unifrac~Collection*Management, as(sample_data(physeq.css), "data.frame")) #Are there significant changes?
# The R2 is the variation due to different factors 
# pvalue is the pvalue for the factor

anova(betadisper(unifrac, sample.data.fungi$Management, type = "centroid"))

### Indicator species analysis ### 
#### taxonomy table ####
tax.bac <- physeq.clean.samples %>%
  tax_table() %>%
  as.data.frame()

indicator.management <- indicspecies::multipatt(as.data.frame(t(physeq.css@otu_table)), cluster = physeq.css@sam_data$Management)
summary(indicator.management)

indicators <- indicator.management$sign
indicators$OTU <- rownames(indicators)

indicators2 <- left_join(indicators, tax.bac, by = "OTU")
indicators2$sig <- ifelse(indicators2$p.value <= 0.01 & indicators2$stat >= 0.5, TRUE, FALSE)

indicators2 %>%
  subset(index == 2 & sig == TRUE) %>%
  ggplot(aes(x = 1, fill = order)) +
  geom_col(position = "stack")

ggplot(indicators2, aes(x = stat, y = -log10(p.value), shape = sig)) + 
  geom_point()

### Differential Abundance analysis ###### 

diffabund_conv_notill <- physeq.clean.samples %>%
  phyloseq::subset_samples(Management %in% c("Conventional", "No-Till")) %>%
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

#BiocManager::install("DESeq2")
library(DESeq2)
library(ggrepel)
ibm.cbb <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")
tol.cbb <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")


diagdds = phyloseq_to_deseq2(diffabund_conv_notill, ~Management)

diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 1
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(diffabund_conv_notill)[rownames(sigtab), ], "matrix"))


volcano <- ggplot(sigtab, aes(x = log2FoldChange, y = -log10(padj), color = Phylum)) +
  geom_point() +
  geom_text_repel(data = sigtab[sigtab$padj <= 0.001 & abs(sigtab$log2FoldChange) >= 3,],
                  aes(label = Label), size = 1.5) +
  theme_classic() +
  scale_shape_manual(values = c(20,24), name = "p <= 0.05") + 
  geom_vline(xintercept = 0, linetype = "dashed")

volcano


## ANCOM-BC ###
library(ANCOMBC)



# takes a while to finish
out = ancombc2(data = diffabund_conv_notill, 
               assay_name = NULL,
               tax_level = "Class",
               p_adj_method = "holm", 
               prv_cut = 0.50, 
               fix_formula = "Management",
               group = "Management", 
               struc_zero = TRUE, 
               neg_lb = TRUE, 
               alpha = 0.05, 
               global = TRUE, 
               n_cl = 1, verbose = TRUE)

####### diff abundance test #####
diff.abund <- out$res

diff.abund.plot <- ggplot(diff.abund, aes(x = `lfc_ManagementNo-Till`, y = -log10(`q_ManagementNo-Till`), shape = `diff_ManagementNo-Till`)) + 
  geom_point() +
  theme_classic() + 
  geom_vline(xintercept = 0, lty = "dotted") + 
  scale_fill_manual(values = c(cbbPalette, ibm.cbb, tol.cbb)) + 
  xlab("log fold change") + 
  ylab("-log10(pvalue)") +
  geom_text_repel(data = diff.abund[diff.abund$`diff_ManagementNo-Till` == TRUE,],
                  aes(label = taxon), size = 4)
diff.abund.plot


















