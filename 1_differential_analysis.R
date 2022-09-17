# Load libraries
library(phyloseq)
library(tidyverse)

# Load metadata and move samples names to row names
filepath = "/home/leandro/Downloads/SraRunTable (7).txt"
metadata <- read.table(filepath, header = TRUE)
metadata <- metadata %>% remove_rownames() %>% column_to_rownames(var = "Run")
write.table(metadata,
            file = "/home/leandro/Downloads/SraRunTable(7).tsv",
            sep = '\t',
            row.names = TRUE)

# Subset metadata according to sex
meta_men <- metadata[metadata$host_sex == "male",]
meta_women <- metadata[metadata$host_sex == "female",]

# Calculate standard deviation before considering as covariates
meta_men$age <- 
  (meta_men$Host_Age - mean(meta_men$Host_Age) ) / sd(meta_men$Host_Age)
meta_men$bmi <- 
  (meta_men$host_body_mass_index - mean(meta_men$host_body_mass_index) ) / sd(meta_men$host_body_mass_index)
meta_men$disease <- factor(meta_men$Host_disease)

meta_women$age <- 
  (meta_women$Host_Age - mean(meta_women$Host_Age) ) / sd(meta_women$Host_Age)
meta_women$bmi <- 
  (meta_women$host_body_mass_index - mean(meta_women$host_body_mass_index) ) / sd(meta_women$host_body_mass_index)
meta_women$disease <- factor(meta_women$Host_disease)

# Load feature table to R
OTU_table = read.csv("/home/leandro/Documents/Microbiome_data_analysis/feature-table.txt",
                     header = TRUE, 
                     sep = '\t')

# Subset feature table according to sex
OTU_men <- OTU_table[colnames(OTU_table)[-c(1)] %in% meta_men$Run]
OTU_men <- OTU_men %>% remove_rownames() %>% column_to_rownames(var = "OTUID")
meta_men <- meta_men %>% remove_rownames() %>% column_to_rownames(var = "Run")

OTU_women <- OTU_table[colnames(OTU_table)[-c(1)] %in% meta_women$Run]
OTU_women <- cbind(OTU_table$OTUID, OTU_women)
OTU_women <- OTU_women %>% remove_rownames() %>% column_to_rownames(var = "OTU_table$OTUID")
OTU_women <- OTU_women + 1
meta_women <- meta_women %>% remove_rownames() %>% column_to_rownames(var = "Run")
write.table(meta_women,
            sep = '\t',
            file = "meta_women.tsv",
            row.names = TRUE)


# Create phyloseq object with OTU table and metadata
OTU_men <- otu_table(OTU_men,
                       taxa_are_rows = TRUE)
OTU_men <- phyloseq(OTU_men,
                      sample_data(meta_men))

OTU_women <- otu_table(OTU_women,
                     taxa_are_rows = TRUE)
OTU_women <- phyloseq(OTU_women,
                    sample_data(meta_men))


# Add taxa taxonomy information to phyloseq object
library(qiime2R)

taxa <- qza_to_phyloseq(taxonomy = "/home/leandro/Documents/Microbiome_data_analysis/taxonomy.qza")
tax_table(OTU_men) <- taxa
tax_table(OTU_women) <- taxa

# Load DESeq2 library
library(DESeq2)

# Convert phyloseq objects to DESeq2 objects
dds_men <- phyloseq_to_deseq2(OTU_men,
                              ~ age + bmi + disease)
dds_women <- phyloseq_to_deseq2(OTU_women,
                                ~ age + bmi + disease)

# Calculate dds
dds_men = DESeq(dds_men, test="Wald", fitType="parametric")
dds_women = DESeq(dds_women, test="Wald", fitType="parametric")

# Filter results to padj < 0.01
res_men = results(dds_men, cooksCutoff = FALSE)
alpha = 0.01
sigtab_men = res_men[which(res_men$padj < alpha), ]
sigtab_men = cbind(as(sigtab_men, "data.frame"), as(OTU_men@tax_table[rownames(sigtab_men), ], "matrix"))
sigtab_men <- as.data.frame(sigtab_men)
head(sigtab_men)


write.table(sigtab_men,
            file = "sigtab_men.csv",
            sep = ",",
            row.names = TRUE)


# Plot results for men
library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab_men$log2FoldChange, sigtab_men$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_men$Phylum = factor(as.character(sigtab_men$Phylum), levels=names(x))

# Genus order
x = tapply(sigtab_men$log2FoldChange, sigtab_men$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab_men$Genus = factor(as.character(sigtab_men$Genus), levels=names(x))
ggplot(sigtab_men, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Filter results to padj < 0.01
res_women = results(dds_women, cooksCutoff = FALSE)
sigtab_women = res_women[which(res_women$padj < alpha), ]
sigtab_women = cbind(as(sigtab_women, "data.frame"), as(OTU_women@tax_table[rownames(sigtab_women), ], "matrix"))
sigtab_women <- as.data.frame(sigtab_women)
head(sigtab_women)

write.table(sigtab_women,
            file = "sigtab_women.csv",
            sep = ",",
            row.names = TRUE)


# Plot results for women
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab_women$log2FoldChange, sigtab_women$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_women$Phylum = factor(as.character(sigtab_women$Phylum), levels=names(x))

# Genus order
x = tapply(sigtab_women$log2FoldChange, sigtab_women$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab_women$Genus = factor(as.character(sigtab_women$Genus), levels=names(x))
ggplot(sigtab_women, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
