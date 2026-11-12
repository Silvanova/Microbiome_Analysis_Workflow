#==============================================================================#
#                                                                              #
#  Taxonomy and functional traits: taxa plots of bacterial communities and     #
#  guild assignment of fungi                                                   #
#                                                                              #
#  Written by Davide Barsotti: db@ign.ku.dk                                    #
#  Modified from Sofia Gomes: s.i.gomes@biology.leidenuniv.nl                  #
#                                                                              #
#==============================================================================#

# Set working directory (modify the path below to match your own)
setwd("some_folder_path/my_working_directory")

# Load packages and homemade functions
source("scripts/001_required_packages.R")

# Load the cleaned bacterial OTU table from previous scripts
load("RData/16S/otu_BLAST_clean.RData")


# Taxonomic composition of bacteria --------------------------------------------

# Taxa plots (stacked bar plots) show which taxa dominate the communities and
# how the composition changes between samples or groups of samples. They are
# descriptive: they do not test whether the differences are significant (see
# the beta diversity and differential abundance scripts for statistical tests).

## Relative abundance ####

# Taxa plots show the relative abundance of each taxon, i.e., the proportion of
# the reads of a sample that belong to this taxon. Relative abundances are
# calculated per sample, so they are not affected by differences in library
# size between samples. For this reason, the data do not need to be rarefied
# first (rarefaction would give almost the same proportions, but discard reads).

otu_table_rel <- transform_sample_counts(otu_BLAST_clean, function(x) x / sum(x))

# Check: the relative abundances of each sample should sum to 1
head(sample_sums(otu_table_rel))

## Long format table ####

# To summarise and plot the relative abundances with dplyr and ggplot2, we
# convert the phyloseq object to a long format table (one row per OTU and per
# sample) with phyloseq::psmelt(). The tax_table contains 37 columns (taxonomy
# and BLAST results). To keep the long table small, we only keep the taxonomic
# ranks that we need.
tax_table(otu_table_rel) <- tax_table(otu_table_rel)[, c("phylum", "class",
                                                         "family", "genus")]

rel_long <- psmelt(otu_table_rel) %>%
  select(OTU, Sample, SiteID, Land_use, Abundance, phylum, class, family, genus) %>%
  # OTUs without a confident assignment at a given rank were changed to "NA"
  # in 005_BLAST_filter.R. We rename them "Unassigned" to show them in the plots.
  mutate(across(phylum:genus, ~ if_else(is.na(.x) | .x == "NA", "Unassigned", .x)))

# Order the sites by land use for the plots
rel_long$SiteID <- factor(rel_long$SiteID,
                          levels = unique(rel_long$SiteID[order(rel_long$Land_use,
                                                                rel_long$SiteID)]))

# NOTE: psmelt() keeps the OTUs with an abundance of 0. This is important when
# calculating means: a taxon absent from a sample must count as 0 in the mean,
# not be ignored.

## Dominant OTUs ####

# Which OTUs dominate the bacterial communities of each land use? We calculate
# the mean relative abundance of each OTU across the samples of a land use.
otu_mean <- rel_long %>%
  group_by(Land_use, OTU, phylum, class, family, genus) %>%
  summarise(mean_relab = mean(Abundance), .groups = "drop")

# Show the 5 most abundant OTUs of each land use
otu_mean %>%
  group_by(Land_use) %>%
  slice_max(mean_relab, n = 5) %>%
  select(Land_use, mean_relab, phylum, family, genus) %>%
  print(n = 25)

# Do the same OTUs dominate in all land uses? In this example, one OTU of
# Candidatus Udaeobacter (Verrucomicrobiota) is the most abundant in all land
# uses, with 2.5% (Old Ref) to 7% (Conventional agriculture) of the reads.
# Bacillus (Firmicutes) is among the dominant OTUs in the agricultural sites
# and the grassland, but not in the forests. Even the most abundant OTUs
# represent only a few percent of the reads: soil bacterial communities are
# very diverse. Note that many of the dominant OTUs cannot be assigned to a
# genus.

## Taxa plots ####

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
  colours <- c(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_top),
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

### Phylum ####

phylum_site <- taxa_site(rel_long, rank = "phylum", n_taxa = 10)

# Check: the relative abundances of each site should sum to 1
phylum_site %>% group_by(SiteID) %>% summarise(total = sum(Abundance))

plot_taxa(phylum_site, "Phylum")

# Proteobacteria, Actinobacteriota and Acidobacteriota are the most abundant
# phyla in all sites, as is common in soils. Compare the agricultural sites
# and the grassland with the oak forests: which phyla differ the most?

### Class ####

class_site <- taxa_site(rel_long, rank = "class", n_taxa = 15)
plot_taxa(class_site, "Class")

### Genus ####

genus_site <- taxa_site(rel_long, rank = "genus", n_taxa = 20)
plot_taxa(genus_site, "Genus")

# At the genus level, most reads are unassigned: many soil bacteria have not
# been cultured and described, and their 16S sequences are too different from
# the reference sequences to be assigned to a genus with confidence (see the
# rank-level thresholds in 005_BLAST_filter.R). The 20 most abundant genera
# represent only a small part of the community.

# NOTE: It is common to remove the unassigned reads and to rescale the top
# taxa to 100%. This makes the plot easier to read, but it hides how little of
# the community is actually shown. If you do so, always state it in the legend.

## Variation between samples ####

# The site means hide the variation between the samples of the same site.
# microViz::comp_barplot() draws one bar per sample and calculates the relative
# abundances itself, so it uses the count data.

# Keep only the taxonomic ranks in the tax_table (comp_barplot() would treat
# the BLAST columns as taxonomic ranks)
otu_table_ranks <- otu_BLAST_clean
tax_table(otu_table_ranks) <- tax_table(otu_table_ranks)[, c("kingdom", "phylum",
                                                             "class", "order",
                                                             "family", "genus")]

# tax_fix() replaces the unassigned ranks ("NA") with the name of the closest
# assigned rank above it (e.g., "Bacteria Kingdom"), which comp_barplot()
# requires. Each panel is a site, each bar is a sample.
otu_table_ranks %>%
  tax_fix(unknowns = "NA") %>%
  comp_barplot(tax_level = "phylum", n_taxa = 10, facet_by = "SiteID",
               bar_width = 0.9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# The samples of the same site have a similar composition, but the variation
# between samples can be large for some phyla. Notice, for example, the higher
# abundance of Firmicutes in the agricultural sites.


# Functional groups in bacteria ------------------------------------------------

# Bacterial functional groups are covered in a separate practical.


# Guilds in fungi --------------------------------------------------------------

# Fungi have very different ecological roles in soil: decomposers
# (saprotrophs), mycorrhizal symbionts (e.g., ectomycorrhizal fungi associated
# with oaks, arbuscular mycorrhizal fungi associated with grasses and most
# crops), and pathogens of plants or animals. Grouping fungi by guild (their
# ecological role) or trophic mode (how they obtain nutrients) helps to
# interpret changes in the communities.

# Trophic modes:
#   - Saprotroph: feeds on dead organic matter
#   - Symbiotroph: obtains nutrients by exchanging resources with a host
#     (e.g., mycorrhizal fungi, lichens, endophytes)
#   - Pathotroph: obtains nutrients by harming a host (pathogens, parasites)

# Guild and trophic mode are assigned from the taxonomy, usually at the genus
# level, using a database of the ecology of described fungi. Keep in mind:
#   1. The ecology of most fungi is only known from observational studies.
#   2. Many fungi switch between strategies, e.g., an endophyte can become a
#      pathogen or a saprotroph. A genus can therefore have several guilds.
#   3. The assignments are POTENTIAL functions: the actual function can vary
#      with the environment, the geography or the host.
#   4. Only OTUs assigned to a genus can be assigned to a guild. In soil, a
#      large part of the fungal OTUs cannot be assigned to a genus.

# Databases for fungal guilds:
#   - FUNGuild (Nguyen et al. 2016, Fungal Ecology): guilds and trophic modes
#     of genera and species, with a confidence ranking.
#   - FungalTraits (Polme et al. 2020, Fungal Diversity,
#     https://link.springer.com/article/10.1007/s13225-020-00466-2): primary
#     and secondary lifestyles of 10,210 genera and 92,623 species hypotheses
#     (UNITE). It is distributed as a table (supplementary material of the
#     paper).
#   - fungaltraits R package (Zanne et al. 2020, Biological Reviews): compiles
#     trait data from many sources, including the guilds and trophic modes of
#     FUNGuild (columns ending in "_fg"). We use this package here, because it
#     can be loaded directly in R.

## Prepare the fungal OTU table ####

# Load the raw ITS OTU table (from 002_data_preparation.R)
load("RData/ITS/otu_table_raw.RData")

# NOTE: In this workflow, the cleaning steps (scripts 003 to 006) were only
# applied to the 16S data. For this practical, we apply a few quick filters to
# the ITS data. For your own data, follow the complete workflow first.

# Keep the samples only (remove negative and positive controls)
otu_its <- subset_samples(otu_table_raw, sample_type == "sample")

# Keep Fungi only (the ITS primers also amplify other eukaryotes)
otu_its <- subset_taxa(otu_its, kingdom == "Fungi")

# Remove singletons (and the OTUs that were only present in the controls)
otu_its <- prune_taxa(taxa_sums(otu_its) > 1, otu_its)

# Remove samples with a low number of fungal reads
sort(sample_sums(otu_its))[1:5]
otu_its <- prune_samples(sample_sums(otu_its) >= 1000, otu_its)

# Count OTUs and samples before and after filtering
ntaxa(otu_table_raw) # 16691
ntaxa(otu_its) # 5827 (35% of the OTUs are Fungi, non-singletons)
nsamples(otu_its) # 87 (2 samples with < 1000 fungal reads removed)

## Genus-level taxonomy ####

# Extract the taxonomy table. In the ITS taxonomy, unassigned ranks are "."
its_tax <- as.data.frame(tax_table(otu_its)) %>%
  mutate(pident = as.numeric(pident))

# As for the 16S data (005_BLAST_filter.R), we only keep the genus assignment
# when the BLAST identity is high enough. For the fungal ITS region, a
# threshold of 94.3% is commonly used for the genus level (Vu et al. 2019,
# Studies in Mycology).
its_tax <- its_tax %>%
  mutate(genus = if_else(genus == "." | pident < 94.3, NA_character_, genus))

# How many OTUs are assigned to a genus?
sum(!is.na(its_tax$genus)) # 2079 out of 5827 OTUs (36%)

## Load the fungaltraits database ####

# The fungaltraits package is installed from GitHub (only needed once)
if (!require("fungaltraits", quietly = TRUE)) {
  if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
  remotes::install_github("ropenscilabs/datastorr")
  remotes::install_github("traitecoevo/fungaltraits")
}
library(fungaltraits)

# Retrieve the database. It is downloaded the first time, which may take a
# while.
dbFT <- fungal_traits()

# The database contains one row per record (a species in a study), and more
# than 100 columns of traits. Most traits are only known for a few species.
dim(dbFT)
head(dbFT)

# We only need the genus, guild and trophic mode. Select the columns by name
# rather than by number: the column order can change between versions of the
# database.
ft_guild <- dbFT %>%
  select(Genus, guild_fg, trophic_mode_fg) %>%
  # Keep only the genera present in our data
  filter(Genus %in% its_tax$genus) %>%
  # Remove the records without a guild assignment
  na.omit() %>%
  filter(!guild_fg %in% c("NULL", "NA", "")) %>%
  # The same genus appears in many records: remove the duplicates
  distinct()

head(ft_guild)

# How many of our genera were found in the database?
n_distinct(its_tax$genus, na.rm = TRUE)
n_distinct(ft_guild$Genus)

## Genera with several guilds ####

# Some genera have several records with different guilds or trophic modes.
# Let's see which ones:
ft_guild %>% dplyr::count(Genus) %>% filter(n > 1)

# For example, Fusarium:
ft_guild %>% filter(Genus == "Fusarium")

# If we join this table to our taxonomy table, every OTU of these genera is
# duplicated, once per record. The duplicated OTUs would be counted several
# times in all subsequent analyses.
nrow(its_tax)
nrow(left_join(its_tax, ft_guild, by = c("genus" = "Genus"),
               relationship = "many-to-many"))

# NOTE: Without relationship = "many-to-many", left_join() gives a warning
# when rows are duplicated. Do not ignore this warning in your own joins.

# To keep one row per genus, we combine all the guilds and trophic modes of a
# genus into one record. E.g., for Fusarium: "Plant Pathogen" +
# "Endophyte-Plant Pathogen-Wood Saprotroph" becomes
# "Endophyte-Plant Pathogen-Wood Saprotroph". This keeps all the potential
# functions of a genus.

# combine_guilds() splits the guilds of all records (separated by "-"), removes
# the duplicates and pastes them back together in alphabetical order. FUNGuild
# sometimes surrounds a guild with "|", which we remove first.
combine_guilds <- function(x) {
  x %>%
    gsub("\\|", "", .) %>%
    strsplit("-") %>%
    unlist() %>%
    trimws() %>%
    unique() %>%
    sort() %>%
    paste(collapse = "-")
}

ft_genus <- ft_guild %>%
  group_by(Genus) %>%
  summarise(guild = combine_guilds(guild_fg),
            trophic_mode = combine_guilds(trophic_mode_fg))

# Check: one row per genus
ft_genus %>% filter(Genus == "Fusarium")
nrow(ft_genus) == n_distinct(ft_genus$Genus)

# NOTE: This is one option among others. You could also keep only the most
# frequent or the most confident assignment of each genus (FUNGuild gives a
# confidence ranking in the column confidence_fg). Choose the option that fits
# your research question, and report it.

## Add the guilds to the taxonomy table ####

its_guild <- its_tax %>%
  select(taxa_OTU, phylum, class, order, family, genus, x1) %>%
  left_join(ft_genus, by = c("genus" = "Genus"))

# Check: the number of rows has not changed
nrow(its_tax) == nrow(its_guild)

# For the plots, we simplify the trophic modes: OTUs with several trophic
# modes are grouped into "Multiple", and we distinguish two types of
# unassigned OTUs.
its_guild <- its_guild %>%
  mutate(trophic_group = case_when(
    is.na(genus) ~ "Unknown genus",
    is.na(trophic_mode) ~ "Genus not in database",
    grepl("-", trophic_mode) ~ "Multiple",
    TRUE ~ trophic_mode))

# How many OTUs are in each group?
table(its_guild$trophic_group)

# Put the new taxonomy table back into the phyloseq object
row.names(its_guild) <- its_guild$taxa_OTU
otu_its_guild <- otu_its
tax_table(otu_its_guild) <- tax_table(as.matrix(its_guild))

## Cross-check with the UNITE EcM lineages ####

# The ITS taxonomy table contains a column x1, with the ectomycorrhizal (EcM)
# lineage of each OTU from UNITE ("non-EcM" for fungi that are not EcM). This
# is an independent assignment that we can compare with the database: OTUs
# with a guild containing "Ectomycorrhizal" should belong to an EcM lineage.
table(guild_EcM = grepl("Ectomycorrhizal", its_guild$guild),
      UNITE_lineage = case_when(its_guild$x1 == "non-EcM" ~ "non-EcM",
                                its_guild$x1 %in% c(".", "unverified") ~ "unknown",
                                TRUE ~ "EcM lineage"))

# The two sources do not always agree: UNITE assigns EcM lineages to OTUs that
# are not assigned to a genus, and some genera contain both EcM and non-EcM
# species. Disagreements show the limits of genus-level guild assignment.

## Trophic modes per site ####

# Relative abundances, as for the bacteria
otu_its_rel <- transform_sample_counts(otu_its_guild, function(x) x / sum(x))
tax_table(otu_its_rel) <- tax_table(otu_its_rel)[, c("genus", "guild",
                                                     "trophic_group")]

its_long <- psmelt(otu_its_rel) %>%
  select(OTU, Sample, SiteID, Land_use, Abundance, genus, guild, trophic_group)

# Order the sites by land use for the plots
its_long$SiteID <- factor(its_long$SiteID,
                          levels = unique(its_long$SiteID[order(its_long$Land_use,
                                                                its_long$SiteID)]))

# What proportion of the reads could be assigned to a trophic mode?
its_long %>%
  group_by(Sample, assigned = !trophic_group %in% c("Unknown genus",
                                                    "Genus not in database")) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  filter(assigned) %>%
  summarise(mean = mean(Abundance), min = min(Abundance), max = max(Abundance))

# We can reuse the functions written for the bacteria. taxa_site() treats the
# group "Unassigned" separately, so we combine the two unassigned groups.
trophic_site <- its_long %>%
  mutate(trophic_group = if_else(trophic_group %in% c("Unknown genus",
                                                      "Genus not in database"),
                                 "Unassigned", trophic_group)) %>%
  taxa_site(rank = "trophic_group", n_taxa = 10)

plot_taxa(trophic_site, "Trophic mode")

## Key guilds per site ####

# Some guilds are of particular interest for the comparison of land uses.
# For example, we expect more ectomycorrhizal fungi in the oak forests (oaks
# are EcM hosts), and more arbuscular mycorrhizal fungi in the grassland and
# agricultural sites (grasses and most crops are AM hosts).

# For each sample, we calculate the relative abundance of the OTUs that have
# each guild among their potential guilds. An OTU with several guilds is
# counted in each of its guilds, so the guilds do not sum to 100%.
guilds_key <- c("Ectomycorrhizal", "Arbuscular Mycorrhizal",
                "Plant Pathogen", "Endophyte")

guild_sample <- lapply(guilds_key, function(g) {
  its_long %>%
    group_by(Sample, SiteID, Land_use) %>%
    summarise(Abundance = sum(Abundance[grepl(g, guild)]), .groups = "drop") %>%
    mutate(guild = g)
}) %>%
  bind_rows()

# Plot the relative abundance of each guild per site
ggplot(guild_sample, aes(x = SiteID, y = Abundance * 100, fill = Land_use)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 1, alpha = 0.5, width = 0.1) +
  facet_wrap(~ guild, scales = "free_y") +
  ylab("Relative abundance (% of fungal reads)") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# NOTE: Arbuscular mycorrhizal fungi (Glomeromycota) are poorly amplified by
# the common ITS primers, so their abundance is usually underestimated. Use
# specific primers (e.g., for the 18S or 28S rRNA genes) to study them.

# As for alpha diversity, the land uses are not replicated (one site per land
# use, except for the oak forests). Compare the sites and use the plots to
# describe patterns, but do not test differences between land uses with the
# samples as replicates.


# Save fungal guild tables -----------------------------------------------------

# Taxonomy table with guilds and trophic modes
save(its_guild, file = "RData/ITS/its_guild.RData")

# Fungal OTU table (phyloseq) with guilds and trophic modes in the tax_table
save(otu_its_guild, file = "RData/ITS/otu_its_guild.RData")

# The taxonomy table can also be saved as a .csv file, e.g., to share it
# write.csv(its_guild, "its_guild.csv", row.names = FALSE)
