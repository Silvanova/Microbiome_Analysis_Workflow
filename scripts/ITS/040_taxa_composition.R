#==============================================================================#
#                                                                              #
#  Taxonomic composition: dominant OTUs, taxa plots, and fungal guilds and     #
#  lifestyles with FungalTraits                                                #
#                                                                              #
#  Modified from Sofia Gomes: s.i.gomes@biology.leidenuniv.nl                  #
#  Modified by Davide Barsotti: db@ign.ku.dk                                   #
#                                                                              #
#==============================================================================#

# Set working directory (modify the path below to match your own)
setwd("some_folder_path/my_working_directory")

# Load packages and homemade functions
source("scripts/001_required_packages.R")

# The figures of the ITS scripts are saved in a separate folder
dir.create("figures/ITS", showWarnings = FALSE)

# Load the cleaned OTU table from previous scripts
load("RData/ITS/otu_BLAST_clean.RData")

# Set seeds: for reproducibility, set seeds (random number generation) at the
# start of this script to always obtain the same results with rarefaction
set.seed(100)


# Rarefaction and relative abundances ------------------------------------------

# Taxa plots (stacked bar plots) show which taxa dominate the communities and
# how the composition changes between sites. They are descriptive: they do not
# test whether the differences are significant (see the beta diversity and
# indicator species scripts for statistical tests).

# As in 010_alpha_diversity.R, we start from the data rarefied to 1011 reads
# (same seed, so the same rarefied table). All samples then have the same
# number of reads.
otu_table_rarefy <- rarefy_even_depth(otu_BLAST_clean, sample.size = 1011,
                                      replace = FALSE, rngseed = 100)
nsamples(otu_table_rarefy) # 87 samples
ntaxa(otu_table_rarefy)    # 5153 OTUs

# Taxa plots show relative abundances: the proportion of the reads of a sample
# that belong to an OTU or a taxon.
otu_table_rel <- transform_sample_counts(otu_table_rarefy, function(x) x / sum(x))

# Check: the relative abundances of each sample should sum to 1
head(sample_sums(otu_table_rel))

## Long format table ####

# To summarise and plot the relative abundances with dplyr and ggplot2, we
# convert the phyloseq object to a long format table (one row per OTU and per
# sample) with phyloseq::psmelt(). The tax_table contains 37 columns (taxonomy
# and BLAST results). To keep the long table small, we only keep the taxonomic
# ranks.
tax_table(otu_table_rel) <- tax_table(otu_table_rel)[, c("phylum", "class",
                                                         "order", "family",
                                                         "genus", "species")]

rel_long <- psmelt(otu_table_rel) %>%
  select(OTU, Sample, SiteID, Land_use, Abundance,
         phylum, class, order, family, genus, species) %>%
  # OTUs without a confident assignment at a given rank were changed to NA in
  # 005_BLAST_filter.R. We rename them "Unassigned" to show them in the plots.
  mutate(across(phylum:species, ~ if_else(is.na(.x) | .x == "NA",
                                          "Unassigned", .x)))

# Order the sites by land use for the plots
rel_long$SiteID <- factor(rel_long$SiteID,
                          levels = unique(rel_long$SiteID[order(rel_long$Land_use,
                                                                rel_long$SiteID)]))

# NOTE: psmelt() keeps the OTUs with an abundance of 0. This is important when
# calculating means: an OTU absent from a sample must count as 0 in the mean,
# not be ignored.


# Dominant OTUs ----------------------------------------------------------------

# Which OTUs dominate the fungal communities of each site? We calculate the
# mean relative abundance of each OTU across the samples of a site.
otu_mean <- rel_long %>%
  group_by(SiteID, Land_use, OTU, phylum, class, order, family, genus, species) %>%
  summarise(mean_relab = mean(Abundance), .groups = "drop")

# The 10 most abundant OTUs in the organic agriculture site...
otu_mean %>%
  filter(SiteID == "Org_agri") %>%
  slice_max(mean_relab, n = 10) %>%
  select(mean_relab, phylum, order, family, genus, species)

# ... and in the old oak forest (reference site)
otu_mean %>%
  filter(SiteID == "Oak_ref") %>%
  slice_max(mean_relab, n = 10) %>%
  select(mean_relab, phylum, order, family, genus, species)

# The dominant OTUs differ between the two sites. In the agricultural site, the
# most abundant OTUs belong to typical soil saprotrophs and yeasts (e.g.,
# Plectosphaerellaceae, Mortierella, Solicoccozyma). In the old oak forest,
# ectomycorrhizal genera such as Russula and Inocybe, which live in symbiosis
# with the roots of oaks, are among the dominant OTUs. Even the most abundant
# OTUs represent only a few percent of the reads, and many of them cannot be
# assigned to a genus.

# Try other sites by changing the SiteID in the code above.


# Taxa plots -------------------------------------------------------------------

# There are too many taxa to give each of them a colour, especially at low
# taxonomic ranks (e.g., genus). We therefore show only the most abundant
# taxa, and group all other taxa into "Other". OTUs without a confident
# assignment at the chosen rank are shown as "Unassigned".

# How many taxa are there at each rank?
rel_long %>% summarise(across(phylum:genus, n_distinct))

# The steps are the same for every taxonomic rank, so we write them in a
# function. taxa_site() summarises the relative abundances at a chosen rank
# and returns the mean relative abundance of each taxon per site:
#   - data: the long format table
#   - rank: the taxonomic rank to summarise (e.g., "phylum")
#   - n_taxa: the number of most abundant taxa to show
taxa_site <- function(data, rank, n_taxa) {
  # 1. Sum the relative abundances of all OTUs of a taxon, in each sample
  taxa_sample <- data %>%
    group_by(Sample, SiteID, Land_use, taxon = .data[[rank]]) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")

  # 2. Find the n most abundant taxa (mean across all samples)
  top_taxa <- taxa_sample %>%
    filter(taxon != "Unassigned") %>%
    group_by(taxon) %>%
    summarise(mean_relab = mean(Abundance)) %>%
    slice_max(mean_relab, n = n_taxa) %>%
    pull(taxon)

  # 3. Group all other taxa into "Other", then calculate the mean per site
  taxa_sample %>%
    mutate(taxon = if_else(taxon %in% c(top_taxa, "Unassigned"), taxon, "Other")) %>%
    group_by(Sample, SiteID, Land_use, taxon) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    group_by(SiteID, Land_use, taxon) %>%
    summarise(Abundance = mean(Abundance), .groups = "drop") %>%
    # Order the taxa from most to least abundant, followed by Other and Unassigned
    mutate(taxon = factor(taxon, levels = c(top_taxa, "Other", "Unassigned")))
}

# plot_taxa() draws the stacked bar plot. Colours are interpolated from the
# RColorBrewer palette "Set3" (12 colours), so that any number of taxa can be
# shown. "Other" and "Unassigned" are shown in grey.
#   - taxa_df: the table returned by taxa_site()
#   - legend_title: the title of the legend (e.g., "Phylum")
plot_taxa <- function(taxa_df, legend_title) {
  n_top <- nlevels(taxa_df$taxon) - 2
  colours <- c(colorRampPalette(brewer.pal(12, "Set3"))(n_top),
               "grey60", "grey85")
  names(colours) <- levels(taxa_df$taxon)

  ggplot(taxa_df, aes(x = SiteID, y = Abundance, fill = taxon)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = colours, name = legend_title) +
    # One panel per land use, with the panel width proportional to the
    # number of sites
    facet_grid(~ Land_use, scales = "free_x", space = "free_x",
               labeller = label_wrap_gen(width = 12)) +
    ylab("Mean relative abundance") +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

## Phylum ####

phylum_site <- taxa_site(rel_long, rank = "phylum", n_taxa = 10)

# Check: the relative abundances of each site should sum to 1
phylum_site %>% group_by(SiteID) %>% summarise(total = sum(Abundance))

p_phylum <- plot_taxa(phylum_site, "Phylum")
p_phylum
ggsave("figures/ITS/040_taxa_phylum.png", p_phylum, width = 9, height = 5)

# Ascomycota and Basidiomycota dominate in all sites, as is common in soils.
# Ascomycota dominate in the agricultural sites and the grassland (46-62% of
# the reads), Basidiomycota in the oak forests (~60%), where many of them are
# ectomycorrhizal (see the guilds below). Mortierellomycota are more abundant
# in the agricultural sites (14-21%) than in the forests (6-8%).

## Genus ####

genus_site <- taxa_site(rel_long, rank = "genus", n_taxa = 20)

p_genus <- plot_taxa(genus_site, "Genus")
p_genus
ggsave("figures/ITS/040_taxa_genus.png", p_genus, width = 9, height = 6)

# At the genus level, 17% (Oak_20) to 50% (Conv_agri) of the reads are
# unassigned: many soil fungi
# have not been described, or their ITS sequences are too different from the
# reference sequences to be assigned to a genus with confidence (see the
# rank-level thresholds in 005_BLAST_filter.R). The 20 most abundant genera
# represent only a small part of the community.

# NOTE: It is common to remove the unassigned reads and to rescale the top
# taxa to 100%. This makes the plot easier to read, but it hides how little of
# the community is actually shown. If you do so, always state it in the legend.

## Variation between samples ####

# The site means hide the variation between the samples of the same site.
# microViz::comp_barplot() draws one bar per sample. It calculates the relative
# abundances itself, so we give it the rarefied counts.

# Keep only the taxonomic ranks in the tax_table (comp_barplot() would treat
# the BLAST columns as taxonomic ranks)
otu_table_ranks <- otu_table_rarefy
tax_table(otu_table_ranks) <- tax_table(otu_table_ranks)[, c("kingdom", "phylum",
                                                             "class", "order",
                                                             "family", "genus")]

# tax_fix() replaces the unassigned ranks (NA) with the name of the closest
# assigned rank above it (e.g., "Fungi Kingdom"), which comp_barplot()
# requires. Each panel is a site, each bar is a sample.
p_samples <- otu_table_ranks %>%
  tax_fix() %>%
  comp_barplot(tax_level = "phylum", n_taxa = 10, facet_by = "SiteID",
               bar_width = 0.9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p_samples
ggsave("figures/ITS/040_taxa_phylum_samples.png", p_samples, width = 10, height = 6)

# The samples of the same site have a similar composition, but the variation
# between samples can be large: in Oak_28, for example, Basidiomycota make up
# 35% to 90% of the reads of a sample.


# Fungal guilds and lifestyles -------------------------------------------------

# Fungi have very different ecological roles in soil: decomposers
# (saprotrophs), mycorrhizal symbionts (e.g., ectomycorrhizal fungi associated
# with oaks, arbuscular mycorrhizal fungi associated with grasses and most
# crops), and pathogens of plants or animals. Grouping fungi by their
# lifestyle or trophic mode (how they obtain nutrients) helps to interpret the
# changes in the communities.

# Trophic modes:
#   - Saprotroph: feeds on dead organic matter
#   - Symbiotroph: obtains nutrients by exchanging resources with a host
#     (e.g., mycorrhizal fungi, lichens, endophytes)
#   - Pathotroph: obtains nutrients by harming a host (pathogens, parasites)

# The most up-to-date database of fungal lifestyles is FungalTraits (Polme et
# al. 2020, Fungal Diversity,
# https://link.springer.com/article/10.1007/s13225-020-00466-2). It gives the
# primary and secondary lifestyle of more than 10,000 fungal genera.
# Keep in mind:
#   1. The ecology of most fungi is only known from observational studies.
#   2. Many fungi switch between strategies, e.g., an endophyte can become a
#      pathogen or a saprotroph.
#   3. The assignments are POTENTIAL functions: the actual function can vary
#      with the environment, the geography or the host.
#   4. The assignments are made at the genus level: most functional traits of
#      fungi are conserved within genera. Only OTUs assigned to a genus can be
#      assigned to a lifestyle, and in soil a large part of the fungal OTUs
#      cannot be assigned to a genus.

## Load the FungalTraits database ####

# FungalTraits is distributed as a table (supplementary material of the paper).
# We use version 1.2 (16 Dec 2020), as distributed with the R package microeco,
# saved as a .csv file in the raw data folder. One row per genus.
fungal_traits <- read.csv("raw data/FungalTraits_1.2.csv")
dim(fungal_traits) # 10771 genera, 24 columns

# The main columns are the genus (GENUS), its primary and secondary lifestyles,
# and more specific traits (e.g., decay type, ectomycorrhizal exploration
# type). Which primary lifestyles are there?
sort(table(fungal_traits$primary_lifestyle), decreasing = TRUE)

# A few genus names appear twice in the table, for homonymous genera in
# different groups of fungi. None of them is in our data, but a duplicated
# genus would duplicate OTUs in the join below. We keep the first row of each
# genus.
fungal_traits %>% dplyr::count(GENUS) %>% filter(n > 1)
fungal_traits <- fungal_traits %>% distinct(GENUS, .keep_all = TRUE)

# The lifestyle of many genera is not known: their primary_lifestyle is empty
# (""). We replace the empty values by NA, so that these genera are treated
# like the genera that are not in the database.
sum(fungal_traits$primary_lifestyle == "")
fungal_traits <- fungal_traits %>%
  mutate(primary_lifestyle = na_if(primary_lifestyle, ""))

## Trophic modes ####

# FungalTraits gives detailed lifestyles, which we group into the three
# trophic modes. Lifestyles that do not fit any of the three (e.g.,
# "unspecified", "bacterivorous") are grouped into "Other". Genera without a
# lifestyle (NA) stay NA.
lifestyle_to_trophic <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    grepl("saprotroph|sooty_mold", x) ~ "Saprotroph",
    grepl("mycorrhizal|endophyte|lichenized|epiphyte|symbiont|symbiotroph", x) ~ "Symbiotroph",
    grepl("pathogen|parasite|pathotroph", x) ~ "Pathotroph",
    TRUE ~ "Other")
}

# Check the grouping: which trophic mode was given to each lifestyle?
fungal_traits %>%
  filter(!is.na(primary_lifestyle)) %>%
  distinct(primary_lifestyle) %>%
  mutate(trophic_mode = lifestyle_to_trophic(primary_lifestyle)) %>%
  arrange(trophic_mode)

## Add the lifestyles to the taxonomy table ####

# Extract the taxonomy table, and join the FungalTraits columns that we need
# by genus. Unassigned genera (NA) do not match any genus of the database.
# NOTE: The genus names must be written the same way in both tables (here,
# e.g., "Russula"). Some taxonomy tables add prefixes such as "g__Russula":
# remove them before the join.
tax_traits <- as.data.frame(tax_table(otu_table_rarefy)) %>%
  select(taxa_OTU, phylum, class, order, family, genus, species) %>%
  left_join(fungal_traits %>%
              select(GENUS, primary_lifestyle, Secondary_lifestyle),
            by = c("genus" = "GENUS"))

# Check: the number of rows has not changed (one row per OTU)
nrow(tax_traits) == ntaxa(otu_table_rarefy)

tax_traits <- tax_traits %>%
  mutate(
    # Genera without a secondary lifestyle ("") become "none"
    Secondary_lifestyle = case_when(is.na(primary_lifestyle) ~ NA_character_,
                                    Secondary_lifestyle == "" ~ "none",
                                    TRUE ~ Secondary_lifestyle),
    # Primary trophic mode, from the primary lifestyle
    primary_trophic = lifestyle_to_trophic(primary_lifestyle),
    # A genus can have a second lifestyle with another trophic mode. For each
    # trophic mode, we note whether the primary OR the secondary lifestyle
    # belongs to it.
    Symbiotroph = if_else(lifestyle_to_trophic(primary_lifestyle) == "Symbiotroph" |
                            lifestyle_to_trophic(Secondary_lifestyle) == "Symbiotroph",
                          "yes", "no"),
    Pathotroph = if_else(lifestyle_to_trophic(primary_lifestyle) == "Pathotroph" |
                           lifestyle_to_trophic(Secondary_lifestyle) == "Pathotroph",
                         "yes", "no"),
    Saprotroph = if_else(lifestyle_to_trophic(primary_lifestyle) == "Saprotroph" |
                           lifestyle_to_trophic(Secondary_lifestyle) == "Saprotroph",
                         "yes", "no"),
    # OTUs with traits (genus found in FungalTraits)
    traits = if_else(is.na(primary_lifestyle), "no", "yes"))

# How many OTUs have traits?
table(tax_traits$traits)
# Percentage of OTUs without traits
mean(tax_traits$traits == "no") * 100

# Only 2044 of the 5153 OTUs (40%) could be assigned to a lifestyle: the others
# are not assigned to a genus, or the lifestyle of their genus is unknown.

# How many of the OTUs with traits can act as each trophic mode? An OTU can be
# counted in two trophic modes (primary and secondary lifestyles).
tax_traits %>%
  filter(traits == "yes") %>%
  summarise(Symb = sum(Symbiotroph == "yes"),
            Pat = sum(Pathotroph == "yes"),
            Sap = sum(Saprotroph == "yes"))

## Relative abundance of the trophic modes and lifestyles ####

# Join the traits to the long format table. OTUs without traits are labelled
# "No traits", so that the bars show which part of the community could be
# assigned.
traits_long <- rel_long %>%
  select(OTU, Sample, SiteID, Land_use, Abundance) %>%
  left_join(tax_traits %>%
              select(taxa_OTU, primary_trophic, primary_lifestyle,
                     Secondary_lifestyle),
            by = c("OTU" = "taxa_OTU")) %>%
  mutate(across(primary_trophic:Secondary_lifestyle,
                ~ if_else(is.na(.x), "No traits", .x)))

# As in the tutorial, we average the samples per land use. The Oak land use
# pools the samples of the six oak forests of different ages (see the site
# plots above for the individual sites).
# trait_landuse() sums the relative abundances per trait in each sample, then
# calculates the mean per land use (in %):
#   - data: the long format table with traits
#   - trait: the column to summarise (e.g., "primary_trophic")
trait_landuse <- function(data, trait) {
  trait_df <- data %>%
    group_by(Sample, Land_use, trait = .data[[trait]]) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    group_by(Land_use, trait) %>%
    summarise(Mean_Rel_Abundance = mean(Abundance) * 100, .groups = "drop")

  # Order the traits from most to least abundant, with "No traits" last
  trait_order <- trait_df %>%
    filter(trait != "No traits") %>%
    group_by(trait) %>%
    summarise(total = sum(Mean_Rel_Abundance)) %>%
    arrange(desc(total)) %>%
    pull(trait)
  trait_df %>%
    mutate(trait = factor(trait, levels = c(trait_order, "No traits")))
}

# plot_traits() draws the stacked bar plot, with "No traits" in grey
#   - trait_df: the table returned by trait_landuse()
#   - legend_title: the title of the legend
plot_traits <- function(trait_df, legend_title) {
  n_traits <- nlevels(trait_df$trait) - 1
  colours <- c(colorRampPalette(brewer.pal(12, "Set3"))(n_traits), "grey85")
  names(colours) <- levels(trait_df$trait)

  ggplot(trait_df, aes(x = Land_use, y = Mean_Rel_Abundance, fill = trait)) +
    geom_col(position = "stack", width = 0.7) +
    scale_fill_manual(values = colours, name = legend_title) +
    labs(x = "Land use", y = "Mean relative abundance (%)") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

### Primary trophic mode ####

trophic_landuse <- trait_landuse(traits_long, "primary_trophic")
p_trophic <- plot_traits(trophic_landuse, "Primary trophic mode")
p_trophic
ggsave("figures/ITS/040_guilds_trophic_mode.png", p_trophic, width = 7, height = 5)

# 21% (Oak) to 50% (conventional agriculture) of the reads belong to OTUs
# without traits. Among the others, saprotrophs dominate in the agricultural
# sites and the grassland (41-50% of the reads), where symbiotrophs are almost
# absent (~1%). In the oak forests, symbiotrophs are the largest group (47% in
# the oak plantations, 17% in the old reference forest). Pathotrophs represent
# 5-11% of the reads in all land uses.

### Primary lifestyle ####

primary_landuse <- trait_landuse(traits_long, "primary_lifestyle")
p_primary <- plot_traits(primary_landuse, "Primary lifestyle")
p_primary
ggsave("figures/ITS/040_guilds_primary_lifestyle.png", p_primary, width = 8, height = 5)

# The primary lifestyles show which symbionts these are: almost all are
# ectomycorrhizal fungi, which live with the roots of the oaks, and are absent
# from the open land uses (< 0.1%). Soil saprotrophs dominate in the
# agricultural sites and the grassland. Arbuscular mycorrhizal fungi
# (Glomeromycota) are rare everywhere (< 1%), although grasses and most crops
# are their hosts: they are poorly amplified by the common ITS primers, so
# their abundance is underestimated. Use specific primers (e.g., for the 18S
# rRNA gene) to study them.

### Secondary lifestyle ####

# Fungi can have more than one lifestyle. The database also describes a
# secondary lifestyle ("none" if the genus has only one lifestyle). In the
# agricultural sites, many saprotrophs are also root-associated or foliar
# endophytes, i.e. they can live inside living plants.
secondary_landuse <- trait_landuse(traits_long, "Secondary_lifestyle")
p_secondary <- plot_traits(secondary_landuse, "Secondary lifestyle")
p_secondary
ggsave("figures/ITS/040_guilds_secondary_lifestyle.png", p_secondary, width = 8, height = 5)

# NOTE: As for alpha diversity, the land uses are not replicated (one site per
# land use, except for the oak forests). Use these plots to describe patterns,
# not to test differences between land uses with the samples as replicates.
