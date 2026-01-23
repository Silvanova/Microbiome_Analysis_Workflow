#==============================================================================#
#                                                                              #
#  Prepare OTU tables, manage auxiliary data, and convert raw files to .RData  #
#                                                                              #
#  Written by Karelle Rheault: karh@ign.ku.dk                                  #
#                                                                              #
#==============================================================================#

# Set working directory (modify the path below to match your own)
setwd("some_folder_path/my_working_directory")

# Load packages
source("scripts/001_required_packages.R")

# Load auxiliary data ----------------------------------------------------------

# Auxiliary data may contain information about the sites where the samples were
# taken (metadata), chemical analyses, and other types of environmental data.

# Load an Excel sheet (.xlsx)
aux <- read.xlsx("raw data/aux_data.xlsx")

# This auxiliary data file contains site metadata and soil physicochemical analyses.

# 16S --------------------------------------------------------------------------

## Taxonomy table ####

# Load BLAST results
BLAST_raw <- read.delim("raw data/16S_BLAST.txt")

# Verify that the data loaded properly: e.g., check that numeric value should have
# class numeric or integer
summary(BLAST_raw)

# Extract taxonomy from BLAST results
taxonomy_raw <- BLAST_raw %>%
  # Split columns based on delimiter (;)
  separate_wider_delim(cols = X1st_hit, delim = ";",
                       names = c("kingdom", "phylum", "class", "order", "family",
                                 "genus", "species", "x1", "x2", "x3", "x4", "x5",
                                 "x6", "x7", "x8", "x9", "x10", "x11", "x12"),
                       too_few = "align_start") %>%
  # Split kingdom column containing a BLAST ID, based on delimiter (space)
  separate_wider_delim(cols = kingdom, delim = " ", names = c("BLAST_ID", "kingdom"),
                       too_few = "align_end") %>%
  # Filter out OTUs without significant similarity found
  filter(kingdom != "No_significant_similarity_found" & !is.na(kingdom)) %>%
  # Rename qseqid column to OTU
  as.data.frame() %>% dplyr::rename(OTU = qseqid) %>%
  # OTU names must be moved to row names
  `row.names<-`(.$OTU)

# Extract sequences from taxonomy raw
refseq_raw <- taxonomy_raw$query_seq %>% `names<-`(taxonomy_raw$OTU)


## OTU table ####

# Were your samples part of a larger sequencing run containing other samples?
# If your samples were mixed with multiple other samples in a larger library,
# you will need to filter your OTU table to keep only your own samples. You may
# then find OTUs present in none of your samples (they were present in the other
# samples you removed). Follow the procedure below to make sure to remove
# superfluous data before you proceed to the next scripts.

# Load OTU table
otu_table_raw <- read.delim("raw data/16S_table.txt")

# The first column contains OTU names, the rest are the count data per sample.
# Before proceeding to the next steps, make sure that the sample identifiers
# in your auxiliary file match those in the OTU table (column names).

# Add sample identifier (Soil_ID) as row names to aux data
row.names(aux) <- aux$Soil_ID

# Clean the otu_table_raw to keep only relevant samples. Make sure the selected
# columns match your aux data. Potentially missing samples may have not amplified
# properly and were removed during bioinformatic processing. Ask your lab tech.
otu_table_raw <- otu_table_raw %>%
  # Keep only samples from aux data
  select(OTU, contains(rownames(aux))) %>%
  # OTU names must be moved to row names
  as.data.frame() %>% `row.names<-`(.$OTU) %>% select(!OTU)


## Format to phyloseq ####

# Combine data into a phyloseq-class object
otu_table_raw <- phyloseq(otu_table(otu_table_raw, taxa_are_rows = T),
                          sample_data(aux),
                          tax_table(as.matrix(taxonomy_raw)),
                          DNAStringSet(refseq_raw))

# Remove all OTUs with sum == 0
ntaxa(otu_table_raw) # 38500 before filtering
otu_table_raw <- prune_taxa(taxa_sums(otu_table_raw) > 0, otu_table_raw)
ntaxa(otu_table_raw) # 25412 (34% of the OTUs are absent from this subset of samples)


## Save as RData ####

save(BLAST_raw, file = "RData/16S/BLAST_raw.RData")
save(taxonomy_raw, file = "RData/16S/taxonomy_raw.RData")
save(otu_table_raw, file = "RData/16S/otu_table_raw.RData")
save(refseq_raw, file = "RData/16S/refseq_raw.RData")
save(aux, file = "RData/auxfile.RData")



# ITS --------------------------------------------------------------------------

## Taxonomy table ####

# Load BLAST results
BLAST_raw <- read.delim("raw data/ITS_BLAST.txt")

# Verify that the data loaded properly: e.g., check that numeric value should have
# class numeric or integer
summary(BLAST_raw)

# Extract taxonomy from BLAST results
taxonomy_raw <- BLAST_raw %>%
  # Split columns based on delimiter (;)
  separate_wider_delim(cols = X1st_hit, delim = ";",
                       names = c("BLAST_ID","kingdom", "phylum", "class",
                                 "order", "family", "genus", "species", "x1",
                                 "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9",
                                 "x10", "x11", "x12"),
                       too_few = "align_start") %>%
  # Filter out OTUs without significant similarity found
  filter(kingdom != "No_significant_similarity_found" & !is.na(kingdom)) %>%
  # Rename qseqid column to OTU
  as.data.frame() %>% dplyr::rename(OTU = qseqid) %>%
  # OTU names must be moved to row names
  `row.names<-`(.$OTU)

# Extract sequences from taxonomy raw
refseq_raw <- taxonomy_raw$query_seq %>% `names<-`(taxonomy_raw$OTU)


## OTU table ####

# Were your samples part of a larger sequencing run containing other samples?
# If your samples were mixed with multiple other samples in a larger library,
# you will need to filter your OTU table to keep only your own samples. You may
# then find OTUs present in none of your samples (they were present in the other
# samples you removed). Follow the procedure below to make sure to remove
# superfluous data before you proceed to the next scripts.

# Load OTU table
otu_table_raw <- read.delim("raw data/ITS_table.txt")

# The first column contains OTU names, the rest are the count data per sample.
# Before proceeding to the next steps, make sure that the sample identifiers
# in your auxiliary file match those in the OTU table (column names).

# Add sample identifier (Library_ID) as row names to aux data
row.names(aux) <- aux$Library_ID

# Clean the otu_table_raw to keep only relevant samples. Make sure the selected
# columns match your aux data. Potentially missing samples may have not amplified
# properly and were removed during bioinformatic processing. Ask your lab tech.
otu_table_raw <- otu_table_raw %>%
  # Keep only samples from aux data
  select(OTU, contains(rownames(aux))) %>%
  # OTU names must be moved to row names
  as.data.frame() %>% `row.names<-`(.$OTU) %>% select(!OTU)


## Format to phyloseq ####

# Combine data into a phyloseq-class object
otu_table_raw <- phyloseq(otu_table(otu_table_raw, taxa_are_rows = T),
                          sample_data(aux),
                          tax_table(as.matrix(taxonomy_raw)),
                          DNAStringSet(refseq_raw))

# Remove all OTUs with sum == 0
ntaxa(otu_table_raw) # 27699 before filtering
otu_table_raw <- prune_taxa(taxa_sums(otu_table_raw) > 0, otu_table_raw)
ntaxa(otu_table_raw) # 16467 (41% of the OTUs are absent from this subset of samples)


## Save as RData ####

save(BLAST_raw, file = "RData/ITS/BLAST_raw.RData")
save(taxonomy_raw, file = "RData/ITS/taxonomy_raw.RData")
save(otu_table_raw, file = "RData/ITS/otu_table_raw.RData")
save(refseq_raw, file = "RData/ITS/refseq_raw.RData")
save(aux, file = "RData/auxfile.RData")



# CO1 --------------------------------------------------------------------------

## Taxonomy table ####

# Load BLAST results
BLAST_raw <- read.delim("raw data/CO1_BLAST.txt")

# Verify that the data loaded properly: e.g., check that numeric value should have
# class numeric or integer
summary(BLAST_raw)

# Extract taxonomy from BLAST results
taxonomy_raw <- BLAST_raw %>%
  # Split columns based on delimiter (;)
  separate_wider_delim(cols = X1st_hit, delim = ";",
                       names = c("BLAST_ID","kingdom", "phylum", "class",
                                 "order", "family", "genus", "species", "x1",
                                 "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9",
                                 "x10", "x11", "x12"),
                       too_few = "align_start") %>%
  # Filter out OTUs without significant similarity found
  filter(kingdom != "No_significant_similarity_found" & !is.na(kingdom)) %>%
  # Rename qseqid column to OTU
  as.data.frame() %>% dplyr::rename(OTU = qseqid) %>%
  # OTU names must be moved to row names
  `row.names<-`(.$OTU)

# Extract sequences from taxonomy raw
refseq_raw <- taxonomy_raw$query_seq %>% `names<-`(taxonomy_raw$OTU)


## OTU table ####

# Were your samples part of a larger sequencing run containing other samples?
# If your samples were mixed with multiple other samples in a larger library,
# you will need to filter your OTU table to keep only your own samples. You may
# then find OTUs present in none of your samples (they were present in the other
# samples you removed). Follow the procedure below to make sure to remove
# superfluous data before you proceed to the next scripts.

# Load OTU table
otu_table_raw <- read.delim("raw data/CO1_table.txt")

# The first column contains OTU names, the rest are the count data per sample.
# Before proceeding to the next steps, make sure that the sample identifiers
# in your auxiliary file match those in the OTU table (column names).

# Add sample identifier (Soil_ID) as row names to aux data
row.names(aux) <- aux$Soil_ID

# Clean the otu_table_raw to keep only relevant samples. Make sure the selected
# columns match your aux data. Potentially missing samples may have not amplified
# properly and were removed during bioinformatic processing. Ask your lab tech.
otu_table_raw <- otu_table_raw %>%
  # Keep only samples from aux data
  select(OTU, contains(rownames(aux))) %>%
  # OTU names must be moved to row names
  as.data.frame() %>% `row.names<-`(.$OTU) %>% select(!OTU)


## Format to phyloseq ####

# Combine data into a phyloseq-class object
otu_table_raw <- phyloseq(otu_table(otu_table_raw, taxa_are_rows = T),
                          sample_data(aux),
                          tax_table(as.matrix(taxonomy_raw)),
                          DNAStringSet(refseq_raw))

# Remove all OTUs with sum == 0
ntaxa(otu_table_raw) # 33027 before filtering
otu_table_raw <- prune_taxa(taxa_sums(otu_table_raw) > 0, otu_table_raw)
ntaxa(otu_table_raw) # 22511 (32% of the OTUs are absent from this subset of samples)


## Save as RData ####

save(BLAST_raw, file = "RData/CO1/BLAST_raw.RData")
save(taxonomy_raw, file = "RData/CO1/taxonomy_raw.RData")
save(otu_table_raw, file = "RData/CO1/otu_table_raw.RData")
save(refseq_raw, file = "RData/CO1/refseq_raw.RData")
save(aux, file = "RData/auxfile.RData")
