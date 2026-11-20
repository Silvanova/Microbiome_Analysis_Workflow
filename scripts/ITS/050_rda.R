#==============================================================================#
#                                                                              #
#  Redundancy analysis (RDA) and variation partitioning                        #
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
# start of this script to always obtain the same results with rarefaction and
# permutation tests
set.seed(100)

# Order the land uses from agriculture to the old forest (for the plots)
sample_data(otu_BLAST_clean)$Land_use <- factor(
  sample_data(otu_BLAST_clean)$Land_use,
  levels = c("Conventional agriculture", "Organic agriculture", "Grassland",
             "Oak", "Old Ref"))

# Colours for each site and each land use (as in 020_beta_diversity.R)
site_col <- c("Conv_agri" = "bisque", "Org_agri" = "bisque3",
              "Grass" = "orange",
              "Oak_20" = "darkolivegreen1", "Oak_28" = "darkolivegreen2",
              "Oak_33" = "darkolivegreen3", "Oak_44" = "palegreen3",
              "Oak_51" = "darkolivegreen4", "Oak_54" = "darkolivegreen",
              "Oak_ref" = "turquoise4")
land_use_col <- c("Conventional agriculture" = "bisque",
                  "Organic agriculture" = "bisque3",
                  "Grassland" = "orange",
                  "Oak" = "darkolivegreen2",
                  "Old Ref" = "turquoise4")


# Data preparation -------------------------------------------------------------

# Beta diversity (020_beta_diversity.R) showed that the fungal communities
# differ between land uses and sites. Which soil properties are related to
# these differences? A redundancy analysis (RDA) relates the community
# composition (many response variables: the OTUs) to a set of explanatory
# variables (here, the soil properties measured in each sample).

## Rarefied data and Hellinger transformation ####

# As in 010_alpha_diversity.R and 020_beta_diversity.R, we start from the
# data rarefied to 1011 reads (same seed, so the same rarefied table).
otu_table_rarefy <- rarefy_even_depth(otu_BLAST_clean, sample.size = 1011,
                                      replace = FALSE, rngseed = 100)
nsamples(otu_table_rarefy) # 87 samples
ntaxa(otu_table_rarefy)    # 5153 OTUs

# The RDA is a linear method: it fits a linear regression of each OTU on the
# explanatory variables. It is not suited to raw counts with many zeros, so we
# use the Hellinger transformation (square root of the relative abundances).
# An RDA on Hellinger-transformed data preserves the Hellinger distance
# between samples, and is a good choice for community data (Legendre &
# Gallagher, 2001).
otu_table_hel <- transform_sample_counts(otu_table_rarefy,
                                         function(x) sqrt(x / sum(x)))

## Extract the OTU table and the metadata ####

# vegan does not work with phyloseq objects: it needs the OTU table as a
# matrix with the samples as rows and the OTUs as columns. In our phyloseq
# object the OTUs are the rows (taxa_are_rows = TRUE), so the table is
# transposed.
taxa_are_rows(otu_table_hel)
otu_hel <- t(as(otu_table(otu_table_hel), "matrix"))
dim(otu_hel) # 87 samples x 5153 OTUs

# NOTE: phyloseq::ordinate(otu_table_hel, method = "RDA", formula = ~ pH + ...)
# also runs vegan::rda() on a phyloseq object, but the model selection and
# some tests below need to refit the model, which only works reliably when
# vegan gets the matrix and the data frame directly.

# Metadata of the samples
meta_rare <- data.frame(sample_data(otu_table_hel))

# Check that the samples are in the same order in the matrix and the metadata
all(rownames(otu_hel) == rownames(meta_rare))

# The 18 soil properties, from bulk density (BD) to phosphorus (P):
#   BD: bulk density; SWCg: gravimetric soil water content; pH; PO4:
#   phosphate; B, Na, Mg, Al, K, Ca, Mn, Fe: exchangeable elements; TN, TC:
#   total nitrogen and carbon; CN: C/N ratio; NH4, NO3: ammonium, nitrate;
#   P: phosphorus.
soil_vars <- c("BD", "SWCg", "pH", "PO4", "B", "Na", "Mg", "Al", "K", "Ca",
               "Mn", "Fe", "TN", "TC", "CN", "NH4", "NO3", "P")
meta_soil <- meta_rare[, soil_vars]
summary(meta_soil)

# No missing values (rda() would remove the samples with missing values)
sum(is.na(meta_soil))


# Redundancy analysis (RDA) ----------------------------------------------------

## Full model ####

# The formula otu_hel ~ . uses all the columns of meta_soil as explanatory
# variables.
mod_rda <- rda(otu_hel ~ ., data = meta_soil)
mod_rda

# NOTE: the tutorial adds distance = "bray" to rda(), but rda() has no
# distance argument and ignores it: an RDA always uses Euclidean distances,
# here between Hellinger-transformed samples (the Hellinger distance). To
# constrain an ordination of another distance (e.g. Bray-Curtis), use a
# distance-based RDA: dbrda(otu_hel ~ ., data = meta_soil, distance = "bray").

# The "Inertia" is the total variance of the OTU table. The constrained part
# (0.354 of 0.757, 47%) is the variance explained by the soil properties, the
# unconstrained part (53%) is the residual variance. The model has 18
# constrained axes (RDA1 to RDA18, one per explanatory variable) and 68
# unconstrained axes (PC1 to PC68). The first axis explains much more than the
# following ones.
summary(mod_rda)$cont$importance[, 1:5]

## Adjusted R2 ####

# The R2 always increases when variables are added, even random ones. The
# adjusted R2 corrects for the number of explanatory variables (Peres-Neto et
# al., 2006): it is the value to report and compare between models.
RsquareAdj(mod_rda)

# With 18 variables, R2 = 0.47 but the adjusted R2 is only 0.33.

## Permutation tests ####

# Is the model significant? anova.cca() permutes the samples to test whether
# the soil properties explain more variance than expected by chance.
set.seed(100)
anova.cca(mod_rda, permutations = 999)

# The model is significant (F = 3.32, p = 0.001).

# Which variables are significant? With by = "term", the variables are tested
# in sequence: each variable only gets the variance not explained by the
# variables before it. This takes about 30 seconds with 99 permutations (as
# in the tutorial; the smallest possible p-value is then 0.01).
set.seed(100)
anova.cca(mod_rda, by = "term", permutations = 99)

# Almost all the variables are significant. The results depend on the order
# of the variables in the formula: pH comes third and explains a large part
# of the variance (F = 13.7). Variables placed after correlated variables get
# little variance. The tutorial removes boron (B) and manganese (Mn), the
# variables that were not significant in its test.

## Collinearity ####

# Many soil properties are correlated with each other. For example, total
# carbon and total nitrogen:
cor(meta_soil$TC, meta_soil$TN) # 0.95

# The variance inflation factor (VIF) measures how much each variable is
# explained by the others. Values above 10 indicate strong collinearity: the
# effects of these variables cannot be separated, and their arrows in the RDA
# plot are unreliable.
sort(vif.cca(mod_rda))

# TC (206) and TN (231) are almost redundant, and Ca (52), pH (32) and SWCg
# (16) also have high VIFs. Removing B and Mn does not solve this.

## Model selection ####

# Instead of removing variables by hand, we select a smaller set of variables
# with a forward selection (Blanchet et al., 2008). ordiR2step() starts from
# an empty model (otu_hel ~ 1) and adds, at each step, the variable that
# increases the adjusted R2 the most, if it is significant. It stops when no
# variable is significant, or when the adjusted R2 of the full model (0.33) is
# reached, to avoid selecting too many variables.
#
# This takes about 10 minutes with 999 permutations at each step. Do not
# reduce the number of permutations: several variables have p-values close to
# 0.05, so with fewer permutations the selection becomes unstable and keeps
# adding collinear variables (with 199 permutations it even stops with an
# error from the singular value decomposition, 'dgesdd').
mod_rda0 <- rda(otu_hel ~ 1, data = meta_soil)

set.seed(100)
mod_sel <- ordiR2step(mod_rda0, scope = formula(mod_rda),
                      permutations = 999, trace = FALSE)

# Variables selected, in order, with the adjusted R2 after each step
mod_sel$anova

# Save the selection table, so that you do not have to run the selection
# again (it is also used by the website)
rda_selection <- mod_sel$anova
save(rda_selection, file = "RData/ITS/rda_selection.RData")

# pH is selected first: alone, it explains 11.5% of the variance (adjusted
# R2), as expected from the PERMANOVA in 020_beta_diversity.R. Then come the
# C/N ratio, nitrate, ammonium and phosphorus, and a few exchangeable
# elements (see the table above). TC, TN and SWCg are not selected.

# The selected model
mod_sel
RsquareAdj(mod_sel)

# The selected variables explain less variance than the full model, but with
# fewer variables (adjusted R2 0.29, instead of 0.33 with 18 variables).

# Collinearity in the selected model
sort(vif.cca(mod_sel))

# All the VIFs are now around 10 or below.

# Test the selected model
set.seed(100)
anova.cca(mod_sel, permutations = 999)


# Plot the RDA -----------------------------------------------------------------

# vegan can plot the RDA directly (plot(mod_sel)), but ggplot2 gives more
# control. We extract the scores of the samples (sites in vegan's
# vocabulary) and of the explanatory variables (biplot arrows) on the first
# two axes. With scaling = 2, the angles between the arrows reflect the
# correlations between the explanatory variables.
site_scores <- as.data.frame(scores(mod_sel, display = "sites",
                                    choices = c(1, 2), scaling = 2))
# Add the metadata (same order of the samples)
site_scores <- cbind(site_scores, meta_rare[, c("SiteID", "Land_use")])

env_scores <- as.data.frame(scores(mod_sel, display = "bp",
                                   choices = c(1, 2), scaling = 2))
env_scores$Variable <- rownames(env_scores)

# Variance explained by RDA1 and RDA2 (of the total variance), for the axis
# labels
rda_eig <- summary(mod_sel)$cont$importance
rda1_pct <- round(rda_eig[2, 1] * 100, 1)
rda2_pct <- round(rda_eig[2, 2] * 100, 1)

# The two plots use the same code, so we write it in a function.
# plot_rda() draws the samples on the first two RDA axes and the arrows of the
# explanatory variables:
#   - colour_by: the metadata column used for the colours of the samples
#   - colours: the colours for each value of that column
#   - legend: the title of the legend
plot_rda <- function(colour_by, colours, legend) {
  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
    geom_point(data = site_scores,
               aes(x = RDA1, y = RDA2, colour = .data[[colour_by]]),
               size = 3, alpha = 0.85) +
    geom_segment(data = env_scores,
                 aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                 arrow = arrow(length = unit(0.25, "cm")),
                 colour = "black", linewidth = 0.3) +
    geom_text(data = env_scores,
              aes(x = RDA1 * 1.1, y = RDA2 * 1.1, label = Variable),
              colour = "black", size = 3.5) +
    scale_colour_manual(values = colours, breaks = names(colours)) +
    labs(x = paste0("RDA1 (", rda1_pct, "%)"),
         y = paste0("RDA2 (", rda2_pct, "%)"),
         colour = legend) +
    theme_bw() +
    theme(panel.grid = element_blank())
}

## Coloured by land use ####

p_rda_land_use <- plot_rda("Land_use", land_use_col, "Land use")
p_rda_land_use
ggsave("figures/ITS/050_rda_land_use.png", p_rda_land_use,
       width = 8, height = 6)

## Coloured by site ####

p_rda_site <- plot_rda("SiteID", site_col, "Site")
p_rda_site
ggsave("figures/ITS/050_rda_site.png", p_rda_site, width = 8, height = 6)

# The first axis (RDA1, 14.6% of the total variance) separates the
# agricultural sites (left) from the forests (right), as the first axis of the
# PCoA in 020_beta_diversity.R. The agricultural soils have a higher pH and
# bulk density (BD), the forest soils more ammonium (NH4) and a higher C/N
# ratio. The second axis (RDA2, 5.6%) separates the grassland, the reference
# forest (Oak_ref) and some plantations (top, with more magnesium, nitrate and
# sodium) from the other plantations, mainly Oak_20, Oak_28 and Oak_44
# (bottom). The grassland is also associated with calcium. The plantations do
# not line up by age along the axes.
#
# NOTE: an arrow points in the direction in which the variable increases, and
# its length shows how strongly it is correlated with the axes. Samples far
# along an arrow have high values of that variable. Correlation is not
# causation: land use changes many soil properties at the same time, so the
# RDA cannot tell which of them actually affect the fungi.


# Variation partitioning -------------------------------------------------------

# Variation partitioning divides the variation explained by an RDA into the
# parts explained by groups of explanatory variables: the part explained
# uniquely by each group, and the parts shared between groups (Borcard et al.,
# 1992). As in the tutorial, we divide the variables into:
#   - X1, organic matter properties: water content, total C and N, C/N ratio,
#     bulk density;
#   - X2, soil chemistry: pH, phosphate and the exchangeable elements;
#   - X3, the age of the forest (0 for the agricultural sites).
mod_varpart <- varpart(otu_hel,
                       ~ SWCg + TC + TN + CN + BD,
                       ~ pH + PO4 + Na + Mg + Al + K + Ca + Fe + NH4 + NO3 + P,
                       ~ Age,
                       data = meta_rare)
mod_varpart

# The fractions are adjusted R2:
#   - [a], [b], [c]: variation explained uniquely by X1, X2, X3;
#   - [d], [e], [f], [g]: variation shared between groups (it cannot be
#     attributed to one of them because they are correlated);
#   - [h]: residuals, the variation not explained.
#
# Together, the three groups explain 33% of the variation (adjusted R2 of
# "All"), and 67% remains unexplained. Soil chemistry explains the largest
# unique part (15.5%), organic matter 4.8%, and age only 1.4%. 8.3% is shared
# between organic matter and chemistry.

# Shared fractions allocated equally to the groups
summary(mod_varpart)

## Plot ####

png("figures/ITS/050_varpart.png", width = 6, height = 6, units = "in",
    res = 300)
plot(mod_varpart, Xnames = c("Organic matter", "Chemistry", "Age"),
     bg = c("tan", "skyblue", "darkolivegreen3"), alpha = 100, digits = 2,
     id.size = 0.8)
dev.off()

plot(mod_varpart, Xnames = c("Organic matter", "Chemistry", "Age"),
     bg = c("tan", "skyblue", "darkolivegreen3"), alpha = 100, digits = 2,
     id.size = 0.8)

# Negative fractions (if any) are not shown in the plot: they can occur with
# adjusted R2 and are interpreted as zero.

## Test the unique fractions ####

# Only the "testable" fractions can be tested, with an RDA in which the other
# groups are partialled out with Condition(). Is the unique effect of age [c]
# significant?
set.seed(100)
anova.cca(rda(otu_hel ~ Age +
                Condition(SWCg + TC + TN + CN + BD +
                            pH + PO4 + Na + Mg + Al + K + Ca + Fe + NH4 +
                            NO3 + P),
              data = meta_rare),
          permutations = 999)

# Age explains a small but significant part of the variation once the soil
# properties are accounted for (F = 2.46, p = 0.001).


# Exercise ---------------------------------------------------------------------

# Run the PCoA (020_beta_diversity.R) and the RDA only for the oak plantations
# and the old reference forest. Do the same soil properties explain the
# composition of the forest communities?
#
# Hints:
#   - Keep only the samples of the land uses "Oak" and "Old Ref", e.g. with
#     subset_samples() on otu_table_hel, then extract the OTU matrix and the
#     metadata again.
#   - Remove the OTUs with no reads in the remaining samples
#     (prune_taxa(taxa_sums(...) > 0, ...)).
#   - Are the soil properties still collinear? Is Age selected by
#     ordiR2step() if you add it to the scope?
