# KARMA Clean Code Script
# Jordan Sims
# July 2026

# This script contains the code needed to re-produce all the analyses from the KARMA paper.

# INPUTS

# KARMA_V1-filt-tree.rds: a phyloseq object containing the biofilm bacteria and bacterioplankton 16S dataset that has been filtered to remove unassigned taxa, mitochondrial and chloroplast DNA, and taxa found in the negative controls. Negative controls and samples with too few reads have already been removed from the dataset.

# KARMA_environmental_data.csv: CSV file containing all environmental data collected at each site, including metadata, temperature metrics, pH, salinity, nutrient concentrations, benthic community data, coral community data, and macroalgae community data.

# KARMA.xxx.csv: CSV file containing an output produced by the iCAMP function. These outputs can be produced by following the analysis found in this script. Because the analysis takes significant time and computational resources, the outputs have also been provided. See the iCAMP documentation for details on each file.

#---------Read in required packages and data---------

# Packages
library('phyloseq'); packageVersion('phyloseq') # 1.52.0
library('tidyverse'); packageVersion('tidyverse') # 2.0.0
library('vegan'); packageVersion('vegan') # 2.7.2
library('geosphere'); packageVersion('geosphere') # 1.6.8
library('EcolUtils'); packageVersion('EcolUtils') # 0.1
library('iCAMP'); packageVersion('iCAMP') # 1.5.2
library('Maaslin2'); packageVersion('Maaslin2') # 1.22.0

# Custom functions

# Runs assumption tests for continuous variables with two groups to determine which statistical test to use
test_2groups <- function(df, index, group) {
  ifelse(
    shapiro.test(df[[index]])$p.value <= 0.05,
    "Data are not normally distributed. You should use a Mann-Whitney U test, wilcox.test(y ~ x, data)",
    ifelse(
      rstatix::levene_test(df, df[[index]] ~ as.factor(df[[group]])) <= 0.05,
      "Variances are different between groups. You should use a Welch's t-test, t.test(y ~ x, data)",
      "Variances are similar between groups. You should use a Student's t-test, t.test(y ~ x, data, var.equal = T)"
    )
  )
}

# Builds a data frame of Shannon and Simpson diversity values and sample metadata from a phyloseq object
adiv_df <- function(ps) {
  estimate_richness(ps, measures = c("Shannon", "Simpson")) %>% 
    rownames_to_column(var = "Sample") %>% 
    merge(rownames_to_column(data.frame(sample_data(ps)), var = "Sample"), by = "Sample")
}

# Runs assumption tests for continuous variables with three or more groups to determine which statistical test to use
test_3groups <- function(df, index, group) {
  
  model <- lm(df[[index]] ~ df[[group]], data = df)
  
  ifelse(
    shapiro.test(residuals(model))$p.value <= 0.05,
    "The residuals are not normally distributed. You should use a Kruskal-Wallis test, kruskal.test(y ~ x, data)",
    ifelse(
      sum(check_normal_group(df, index, group)$p.value <= 0.05) >= 1,
      "The samples are not normally distributed. You should use a Kruskal-Wallis test, kruskal.test(y ~ x, data)",
      ifelse(
        rstatix::levene_test(df, df[[index]] ~ as.factor(df[[group]]))$p <= 0.05,
        "Residuals and samples are normally distributed, but variances are different between groups. You should use Welch's ANOVA, welch_anova_test(data, y ~ x)",
        "Residuals and samples are normally distributed, and variances are similar between groups. You should use ANOVA, aov(y ~ x, data)"
      )
    )
  )
}

# Helper function for test_3groups
check_normal_group <- function(df, index, group) {
  
  group_list <- vector(mode = "list", length = length(unique(df[[group]])))
  
  names(group_list) <- unique(df[[group]])
  
  for (i in unique(df[[group]])) {
    group_list[[i]] <- 
      rstatix::shapiro_test(
        df[df[[group]] == i,][,index])
  }
  
  output <-
    do.call(rbind, group_list) %>% 
    select(-variable) %>% 
    mutate(Group = unique(df[[group]]))
}

# Data

ps <- readRDS("data/KARMA_V1-filt-tree.rds")

# Split ps object into water, all biofilms, CCA only, and carbonate only

water <- subset_samples(ps, SampleType == "Water")
water <- prune_taxa(taxa_sums(water) > 0, water)

bio <- subset_samples(ps, SampleType == "Biofilm")
bio <- prune_taxa(taxa_sums(bio) > 0, bio)

CCA <- subset_samples(ps, Substrate == "CCA")
CCA <- prune_taxa(taxa_sums(CCA) > 0, CCA)

carb <- subset_samples(ps, Substrate == "Carbonate")
carb <- prune_taxa(taxa_sums(carb) > 0, carb)

#---------Environmental conditions among sites---------

# Get site-level environmental data

env <- 
  read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, AvgSiteTemp:Silicate) %>%
  filter(!is.na(Site)) %>% 
  distinct() %>% 
  column_to_rownames(var = "Site")

# Standardize data by centering and scaling
env_stand <- decostand(env, method = "standardize")

# Calculate euclidean distance to get environmental distance matrix
env_dist <- vegdist(env_stand, method = "euclidian")

# Read in lat-long data
geo <- 
  read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, Longitude, Latitude) %>% 
  column_to_rownames(var = "Site")

# Calculate geographic distance using Haversine distance
geo_dist <- 
  distm(geo, fun = distHaversine) %>% 
  as.dist()

# Run Mantel test

# Do environmental conditions co-vary with geographic distance?

mantel(env_dist, geo_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.0298, p = 0.3896

# Coral community vs geography

# Get coral community data

coral_cov <-
  read_csv("data/KARMA_environmental_data.csv") %>%
  select(Site = Name, ACER:SSID) %>% 
  column_to_rownames(var = "Site") %>% 
  data.frame()

# Hellinger transform

coral.hel <- 
  decostand(coral_cov, method = "hellinger") %>% 
  data.frame()

# Calculate bray curtis distances

coral_dist <- vegdist(coral.hel, method = "bray")

# Run Mantel test

# Does coral community composition co-vary with geographic distance?

mantel(coral_dist, geo_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.2628, p = 0.0214

# Algae community vs geography

# Get algae community data

alg_cov <-
  read_csv("data/KARMA_environmental_data.csv") %>%
  select(Site = Name, Amphiroa:TurfAlgae) %>% 
  column_to_rownames(var = "Site") %>% 
  data.frame()

# Hellinger transform

alg.hel <- 
  decostand(alg_cov, method = "hellinger") %>% 
  data.frame()

# Calculate bray curtis distances

alg_dist <- vegdist(alg.hel, method = "bray")

# Run Mantel test

# Does algae community composition co-vary with geographic distance?

mantel(alg_dist, geo_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.03903, p = 0.5936

# Bacterioplankton community vs geography

# Hellinger transform the water OTU table

water_otu <- otu_table(water) %>% t() %>% data.frame()

water_otu.hel <- decostand(water_otu, method = "hellinger")

# Put into phyloseq object

water_hel <- water

otu_table(water_hel) <- otu_table(water_otu.hel, taxa_are_rows = FALSE)

# Calculate unweighted and weighted Unifrac distances

uu_dist_water_hel <- UniFrac(water_hel, weighted = FALSE, normalized = TRUE)
wu_dist_water_hel <- UniFrac(water_hel, weighted = TRUE, normalized = TRUE)

# Run Mantel tests

# Do bacterioplankton communities co-vary with geographic distance?

mantel(uu_dist_water_hel, geo_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.08622, p = 0.2336

mantel(wu_dist_water_hel, geo_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.09264, p = 0.2094

# Do bacterioplankton communities co-vary with environmental conditions?

mantel(uu_dist_water_hel, env_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.1358, p = 0.6965

mantel(wu_dist_water_hel, env_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.0003608, p = 0.4676

#---------Text S6---------

# Which coral species is most abundant, and how much of the community does it make up?

read_csv("data/KARMA_environmental_data.csv") %>%
  select(Site = Name, ACER:SSID) %>% 
  column_to_rownames(var = "Site") %>% 
  mutate(Coral = rowSums(.)) %>% 
  mutate(across(ACER:SSID, .fns = ~./Coral)) %>% 
  select(-Coral) %>% 
  rownames_to_column(var = "Site") %>% 
  pivot_longer(cols = ACER:SSID, names_to = "Species", values_to = "Abundance") %>% 
  select(-Site) %>% 
  group_by(Species) %>% 
  mutate(sd = sd(Abundance),
         Abundance = mean(Abundance)) %>% 
  distinct() %>% 
  arrange(desc(Abundance)) %>% 
  head(1)

# What is the range of values for coral cover?

c(
  read_csv("data/KARMA_environmental_data.csv") %>%
    pull(Coral) %>% 
    min(),
  read_csv("data/KARMA_environmental_data.csv") %>%
    pull(Coral) %>% 
    max())

# Where was coral cover highest?

read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, Coral) %>% 
  arrange(desc(Coral)) %>% 
  head(1)

# Where was coral cover lowest?

read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, Coral) %>% 
  arrange(Coral) %>% 
  head(1)

# Where was coral richness highest?

read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, CoralRich) %>% 
  arrange(desc(CoralRich)) %>% 
  head(1)

# Where was coral richness lowest?

read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, CoralRich) %>% 
  arrange(CoralRich) %>% 
  head(1)

# Where was coral Shannon diversity highest?

read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, CoralShan) %>% 
  arrange(desc(CoralShan)) %>% 
  head(1)

# Where was coral Shannon diversity lowest?

read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, CoralShan) %>% 
  arrange(CoralShan) %>% 
  head(1)

# Where was coral Simpson diversity highest?

read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, CoralSimp) %>% 
  arrange(desc(CoralSimp)) %>% 
  head(1)

# Where was coral Simpson diversity lowest?

read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, CoralSimp) %>% 
  arrange(CoralSimp) %>% 
  head(1)

# Which algae taxon is most abundant, and how much of the community does it make up?

read_csv("data/KARMA_environmental_data.csv") %>%
  select(Site = Name, Macroalgae, Amphiroa:TurfAlgae) %>% 
  column_to_rownames(var = "Site") %>% 
  mutate(across(Amphiroa:TurfAlgae, .fns = ~./Macroalgae)) %>% 
  select(-Macroalgae) %>% 
  rownames_to_column(var = "Site") %>% 
  pivot_longer(cols = Amphiroa:TurfAlgae, names_to = "Taxon", values_to = "Abundance") %>% 
  select(-Site) %>% 
  group_by(Taxon) %>% 
  mutate(sd = sd(Abundance),
         Abundance = mean(Abundance)) %>% 
  distinct() %>% 
  arrange(desc(Abundance)) %>% 
  head(1)

# Where was algae cover highest?

read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, Macroalgae) %>% 
  arrange(desc(Macroalgae)) %>% 
  head(1)

# Where was algae cover lowest?

read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, Macroalgae) %>% 
  arrange(Macroalgae) %>% 
  head(1)

# Where was algae Shannon diversity highest?

read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, MacroShan) %>% 
  arrange(desc(MacroShan)) %>% 
  head(1)

# Where was algae Shannon diversity lowest?

read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, MacroShan) %>% 
  arrange(MacroShan) %>% 
  head(1)

# Where was algae Simpson diversity highest?

read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, MacroSimp) %>% 
  arrange(desc(MacroSimp)) %>% 
  head(1)

# Where was algae Simpson diversity lowest?

read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, MacroSimp) %>% 
  arrange(MacroSimp) %>% 
  head(1)

#---------Bacterial communities among sample types---------

# Remove samples with <2,500 reads from biofilm-only dataset

bio_rare <- subset_samples(bio, sample_sums(bio) > 2500)

# Rarefy to lowest read depth

bio_rare <- 
  rarefy_even_depth(
    bio_rare, 
    sample.size = min(sample_sums(bio_rare)), 
    rngseed = 8, 
    replace = FALSE)

# Calculate diversity indices

bio_rare_div <- adiv_df(bio_rare)

# Determine which statistical tests to use

test_2groups(bio_rare_div, "Shannon", "Substrate")
test_2groups(bio_rare_div, "Simpson", "Substrate")

# Run tests to determine difference in alpha diversity between substrates

t.test(Shannon ~ Substrate, bio_rare_div, var.equal = T) # p = 0.00118
wilcox.test(Simpson ~ Substrate, bio_rare_div) # p = 0.002046

# Which substrate had higher diversity?

bio_rare_div %>% 
  select(Substrate, Shannon, Simpson) %>% 
  group_by(Substrate) %>% 
  mutate(Shannon = mean(Shannon),
         Simpson = mean(Simpson)) %>% 
  distinct()

# Remove samples with <2,500 reads from full dataset

ps_rare <- subset_samples(ps, sample_sums(ps) > 2500)

# Rarefy to lowest read depth
ps_rare <- 
  rarefy_even_depth(
    ps_rare, 
    sample.size = min(sample_sums(ps_rare)), 
    rngseed = 8, 
    replace = FALSE)

# Calculate diversity indices

ps_rare_div <- adiv_df(ps_rare)

# Determine which statistical tests to use

test_3groups(ps_rare_div, "Shannon", "Substrate")
test_3groups(ps_rare_div, "Simpson", "Substrate")

# Run tests to assess differences in alpha diversity by substrate and biofilm vs water

kruskal.test(Shannon ~ Substrate, ps_rare_div) # p = 1.982e-08
kruskal.test(Simpson ~ Substrate, ps_rare_div) # p = 1.148e-07

# Pairwise comparisons

pairwise.wilcox.test(ps_rare_div$Shannon, ps_rare_div$Substrate, p.adjust.method = "BH")
pairwise.wilcox.test(ps_rare_div$Simpson, ps_rare_div$Substrate, p.adjust.method = "BH")

# Which sample type had higher diversity?

ps_rare_div %>% 
  select(Substrate, Shannon, Simpson) %>% 
  group_by(Substrate) %>% 
  mutate(Shannon = mean(Shannon),
         Simpson = mean(Simpson)) %>% 
  distinct()

# Pull sample data for substrate beta diversity comparisons

bio_df <- 
  data.frame(sample_data(bio)) %>% 
  rownames_to_column("SampleID")

# Calculate unweighted and weighted Unifrac distances

uu_dist_bio <- UniFrac(bio, weighted = FALSE, normalized = TRUE)
wu_dist_bio <- UniFrac(bio, weighted = TRUE, normalized = TRUE)

# Tests to assess differences in dispersion between substrates

permutest(betadisper(uu_dist_bio, bio_df$Substrate, bias.adjust = TRUE)) # p = 0.256
permutest(betadisper(wu_dist_bio, bio_df$Substrate, bias.adjust = TRUE)) # p = 0.159

# Pull sample data for sample type beta diversity comparisons

ps_df <- 
  data.frame(sample_data(ps)) %>% 
  rownames_to_column("SampleID")

# Calculate unweighted and weighted Unifrac distances

uu_dist_ps <- UniFrac(ps, weighted = FALSE, normalized = TRUE)
wu_dist_ps <- UniFrac(ps, weighted = TRUE, normalized = TRUE)

# Tests to assess differences in dispersion among sample types

permutest(betadisper(uu_dist_ps, ps_df$Substrate, bias.adjust = TRUE), permutations = 9999) # p = 0.0331
permutest(betadisper(wu_dist_ps, ps_df$Substrate, bias.adjust = TRUE), permutations = 9999) # p = 1e-04

# Tests to assess differences in community composition among sample types

adonis2(uu_dist_ps ~ Substrate, data = ps_df, permutations = 9999) # p = 1e-04
adonis2(wu_dist_ps ~ Substrate, data = ps_df, permutations = 9999) # p = 1e-04

# Pairwise comparisons

pairwiseAdonis::pairwise.adonis2(uu_dist_ps ~ Substrate, data = ps_df, nperm = 9999)
pairwiseAdonis::pairwise.adonis2(wu_dist_ps ~ Substrate, data = ps_df, nperm = 9999)

# Hellinger transform the CCA OTU table

CCA_otu <- otu_table(CCA) %>% t() %>% data.frame()

CCA_otu.hel <- decostand(CCA_otu, method = "hellinger")

# Put into phyloseq object

CCA_hel <- CCA

otu_table(CCA_hel) <- otu_table(CCA_otu.hel, taxa_are_rows = FALSE)

# Calculate unweighted and weighted Unifrac distances

uu_dist_CCA_hel <- UniFrac(CCA_hel, weighted = FALSE, normalized = TRUE)
wu_dist_CCA_hel <- UniFrac(CCA_hel, weighted = TRUE, normalized = TRUE)

# Get sample-level environmental data
env_CCA <- 
  sample_data(CCA) %>% 
  data.frame() %>% 
  select(Depth:MacroSimp)

# Standardize data by centering and scaling
env_stand_CCA <- decostand(env_CCA, method = "standardize")

# Calculate euclidean distance to get environmental distance matrix
env_dist_CCA <- vegdist(env_stand_CCA, method = "euclidian")

# Run Mantel tests

# Do CCA communities co-vary with environmental conditions?

mantel(uu_dist_CCA_hel, env_dist_CCA, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.1681, p = 0.0099

mantel(wu_dist_CCA_hel, env_dist_CCA, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.1265, p = 0.0444

# Hellinger transform the carbonate OTU table

carb_otu <- otu_table(carb) %>% t() %>% data.frame()

carb_otu.hel <- decostand(carb_otu, method = "hellinger")

# Put into phyloseq object

carb_hel <- carb

otu_table(carb_hel) <- otu_table(carb_otu.hel, taxa_are_rows = FALSE)

# Calculate unweighted and weighted Unifrac distances

uu_dist_carb_hel <- UniFrac(carb_hel, weighted = FALSE, normalized = TRUE)
wu_dist_carb_hel <- UniFrac(carb_hel, weighted = TRUE, normalized = TRUE)

# Get sample-level environmental data
env_carb <- 
  sample_data(carb) %>% 
  data.frame() %>% 
  select(Depth:MacroSimp)

# Standardize data by centering and scaling
env_stand_carb <- decostand(env_carb, method = "standardize")

# Calculate euclidean distance to get environmental distance matrix
env_dist_carb <- vegdist(env_stand_carb, method = "euclidian")

# Run Mantel tests

# Do carbonate communities co-vary with environmental conditions?

mantel(uu_dist_carb_hel, env_dist_carb, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.05797, p = 0.7153

mantel(wu_dist_carb_hel, env_dist_carb, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.01428, p = 0.5283

# dbRDA for CCA communities with weighted UniFrac distances

# Get data out of CCA object and scale it

sample_data_CCA <- 
  sample_data(CCA) %>% 
  data.frame() %>% 
  select(Depth:MacroSimp) %>% 
  select(-c(MacroSimp, CoralRich, CoralShan, SDSiteTemp)) %>% 
  mutate(across(c(Depth:MacroShan), scale))

# Perform backwards selection for CCA communities with weighted UniFrac

full_dbrda <- dbrda(wu_dist_CCA_hel ~ ., data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(full_dbrda)$adj.r.squared # 0.1085804

anova.cca(full_dbrda, step = 1000, permutations = 9999, by = "term")

# Remove MacroCover
dbrda1 <- dbrda(wu_dist_CCA_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + NitrateNitrite + Ammonium + Silicate + 
                  CoralCover + CoralSimp + MacroShan, 
                data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda1)$adj.r.squared # 0.1003119

anova.cca(dbrda1, step = 1000, permutations = 9999, by = "term")

# Remove CoralSimp
dbrda2 <- dbrda(wu_dist_CCA_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + NitrateNitrite + Ammonium + Silicate + 
                  CoralCover + MacroShan, 
                data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda2)$adj.r.squared # 0.08759505

anova.cca(dbrda2, step = 1000, permutations = 9999, by = "term")

# Remove Depth
dbrda3 <- dbrda(wu_dist_CCA_hel ~ 
                  AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + NitrateNitrite + Ammonium + Silicate + 
                  CoralCover + MacroShan, 
                data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda3)$adj.r.squared # 0.1044593

anova.cca(dbrda3, step = 1000, permutations = 9999, by = "term")

# Remove NitNit
dbrda4 <- dbrda(wu_dist_CCA_hel ~ 
                  AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + Ammonium + Silicate + 
                  CoralCover + MacroShan, 
                data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda4)$adj.r.squared # 0.1056237

anova.cca(dbrda4, step = 1000, permutations = 9999, by = "term")

# Remove AvgSiteTemp
dbrda5 <- dbrda(wu_dist_CCA_hel ~ 
                  AvgSampTemp + Salinity + 
                  pH + Ammonium + Silicate + 
                  CoralCover + MacroShan, 
                data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda5)$adj.r.squared # 0.1143348

anova.cca(dbrda5, step = 1000, permutations = 9999, by = "term")

# Remove Salinity
dbrda6 <- dbrda(wu_dist_CCA_hel ~ 
                  AvgSampTemp + pH + Ammonium + 
                  Silicate + CoralCover + MacroShan, 
                data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda6)$adj.r.squared # 0.1113296

anova.cca(dbrda6, step = 1000, permutations = 9999, by = "term")

# Remove MacroShan
dbrda7 <- dbrda(wu_dist_CCA_hel ~ 
                  AvgSampTemp + pH + Ammonium + 
                  Silicate + CoralCover, 
                data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda7)$adj.r.squared # 0.1161136

anova.cca(dbrda7, step = 1000, permutations = 9999, by = "term")

# Remove Ammonium (BEST MODEL)
dbrda8 <- dbrda(wu_dist_CCA_hel ~ 
                  AvgSampTemp + pH + Silicate + CoralCover, 
                data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda8)$adj.r.squared # 0.1185445

anova.cca(dbrda8, step = 1000, permutations = 9999, by = "term")

# Remove AvgSampTemp (BEST MODEL)
dbrda9 <- dbrda(wu_dist_CCA_hel ~ pH + Silicate + CoralCover, 
                data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda9)$adj.r.squared # 0.1193836

anova.cca(dbrda9, step = 1000, permutations = 9999, by = "term")
#            Df SumOfSqs      F Pr(>F)    
# pH          1  0.05510 1.5100 0.1079    
# Silicate    1  0.15586 4.2716 0.0004 ***
# CoralCover  1  0.05680 1.5566 0.0943 .  
# Residual   29  1.05817 

anova.cca(dbrda9, step = 1000, permutations = 9999)

# Remove pH
dbrda10 <- dbrda(wu_dist_CCA_hel ~ Silicate + CoralCover, 
                 data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda10)$adj.r.squared # 0.07672386

anova.cca(dbrda10, step = 1000, permutations = 9999, by = "term")

# Remove CoralCover
dbrda11 <- dbrda(wu_dist_CCA_hel ~ Silicate, 
                 data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda11)$adj.r.squared # 0.06347309

anova.cca(dbrda11, step = 1000, permutations = 9999, by = "term")

# Perform backwards selection for CCA communities with unweighted UniFrac

full_dbrda <- dbrda(uu_dist_CCA_hel ~ ., data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(full_dbrda)$adj.r.squared # 0.03472906 (BEST MODEL)

anova.cca(full_dbrda, step = 1000, permutations = 9999, by = "term")
#                Df SumOfSqs      F Pr(>F)  
# Depth           1   0.2406 1.0435 0.2781  
# AvgSampTemp     1   0.2883 1.2505 0.0631 .
# AvgSiteTemp     1   0.2307 1.0007 0.3815  
# Salinity        1   0.2887 1.2522 0.0601 .
# pH              1   0.2277 0.9875 0.4254  
# NitrateNitrite  1   0.2101 0.9112 0.7259  
# Ammonium        1   0.2473 1.0727 0.2229  
# Silicate        1   0.3370 1.4618 0.0133 *
# CoralCover      1   0.2355 1.0216 0.3256  
# MacroCover      1   0.1956 0.8484 0.9476  
# CoralSimp       1   0.2626 1.1391 0.1380  
# MacroShan       1   0.2679 1.1622 0.1234  
# Residual       20   4.6107    

anova.cca(full_dbrda, step = 1000, permutations = 9999) # p = 0.0343

# Remove MacroCover
dbrda1 <- dbrda(uu_dist_CCA_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + NitrateNitrite + Ammonium + Silicate + 
                  CoralCover + CoralSimp + MacroShan, 
                data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda1)$adj.r.squared # 0.0257515

anova.cca(dbrda1, step = 1000, permutations = 9999, by = "term")

# Remove CoralSimp
dbrda2 <- dbrda(uu_dist_CCA_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + NitrateNitrite + Ammonium + Silicate + 
                  CoralCover + MacroShan, 
                data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda2)$adj.r.squared # 0.02653967

anova.cca(dbrda2, step = 1000, permutations = 9999, by = "term")

# Remove Nitnit
dbrda3 <- dbrda(uu_dist_CCA_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + Ammonium + Silicate + 
                  CoralCover + MacroShan, 
                data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda3)$adj.r.squared # 0.02753569

anova.cca(dbrda3, step = 1000, permutations = 9999, by = "term")

# Remove Ammonium
dbrda4 <- dbrda(uu_dist_CCA_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + Silicate + 
                  CoralCover + MacroShan, 
                data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda4)$adj.r.squared # 0.02826122

anova.cca(dbrda4, step = 1000, permutations = 9999, by = "term")

# Remove MacroShan
dbrda5 <- dbrda(uu_dist_CCA_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + Silicate + CoralCover, 
                data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda5)$adj.r.squared # 0.02969245

anova.cca(dbrda5, step = 1000, permutations = 9999, by = "term")
           
# Remove pH
dbrda6 <- dbrda(uu_dist_CCA_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  Silicate + CoralCover, 
                data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda6)$adj.r.squared # 0.02120442

anova.cca(dbrda6, step = 1000, permutations = 9999, by = "term")

# Remove AvgSiteTemp
dbrda7 <- dbrda(uu_dist_CCA_hel ~ 
                  Depth + AvgSampTemp + Salinity + 
                  Silicate + CoralCover, 
                data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda7)$adj.r.squared # 0.02316563

anova.cca(dbrda7, step = 1000, permutations = 9999, by = "term")

# Remove Depth
dbrda8 <- dbrda(uu_dist_CCA_hel ~ 
                  AvgSampTemp + Salinity + Silicate + CoralCover, 
                data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda8)$adj.r.squared # 0.02397269

anova.cca(dbrda8, step = 1000, permutations = 9999, by = "term")

# Remove CoralCov
dbrda9 <- dbrda(uu_dist_CCA_hel ~ 
                  AvgSampTemp + Salinity + Silicate, 
                data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda9)$adj.r.squared # 0.01975365

anova.cca(dbrda9, step = 1000, permutations = 9999, by = "term")

# Remove Salinity
dbrda10 <- dbrda(uu_dist_CCA_hel ~ 
                   AvgSampTemp + Silicate, 
                 data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda10)$adj.r.squared # 0.01484261

anova.cca(dbrda10, step = 1000, permutations = 9999, by = "term")

# Remove AvgSampTemp
dbrda11 <- dbrda(uu_dist_CCA_hel ~ Silicate, 
                 data = sample_data_CCA, na.action = na.exclude)

RsquareAdj(dbrda11)$adj.r.squared # 0.007816596

anova.cca(dbrda11, step = 1000, permutations = 9999, by = "term")

# Perform backwards selection for carbonate communities with weighted UniFrac

# Pull sample data

sample_data_carb <- 
  sample_data(carb) %>% 
  data.frame() %>% 
  select(Depth:MacroSimp) %>% 
  select(-c(MacroSimp, CoralRich, CoralShan, SDSiteTemp)) %>% 
  mutate(across(c(Depth:MacroShan), scale))

# Perform backwards selection

full_dbrda <- dbrda(wu_dist_carb_hel ~ ., data = sample_data_carb, na.action = na.exclude)

RsquareAdj(full_dbrda)$adj.r.squared # 0.04953439

anova.cca(full_dbrda, step = 1000, permutations = 9999, by = "term")

# Remove Ammonium
dbrda1 <- dbrda(wu_dist_carb_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + NitrateNitrite + Silicate + MacroCover +
                  CoralCover + CoralSimp + MacroShan, 
                data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda1)$adj.r.squared # 0.02933197

anova.cca(dbrda1, step = 1000, permutations = 9999, by = "term")

# Remove CoralSimp
dbrda2 <- dbrda(wu_dist_carb_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + NitrateNitrite + Silicate + MacroCover +
                  CoralCover + MacroShan, 
                data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda2)$adj.r.squared # 0.03781429

anova.cca(dbrda2, step = 1000, permutations = 9999, by = "term")

# Remove MacroShan
dbrda3 <- dbrda(wu_dist_carb_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + NitrateNitrite + Silicate + MacroCover + CoralCover, 
                data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda3)$adj.r.squared # 0.05950943

anova.cca(dbrda3, step = 1000, permutations = 9999, by = "term")

# Remove MacroCover
dbrda4 <- dbrda(wu_dist_carb_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + NitrateNitrite + Silicate + CoralCover, 
                data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda4)$adj.r.squared # 0.07025746

anova.cca(dbrda4, step = 1000, permutations = 9999, by = "term")

# Remove CoralCover
dbrda5 <- dbrda(wu_dist_carb_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + NitrateNitrite + Silicate, 
                data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda5)$adj.r.squared # 0.08662797

anova.cca(dbrda5, step = 1000, permutations = 9999, by = "term")

# Remove Depth
dbrda6 <- dbrda(wu_dist_carb_hel ~ 
                  AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + NitrateNitrite + Silicate, 
                data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda6)$adj.r.squared # 0.09743275

anova.cca(dbrda6, step = 1000, permutations = 9999, by = "term")

# Remove AvgSampTemp (BEST MODEL)
dbrda7 <- dbrda(wu_dist_carb_hel ~ 
                  AvgSiteTemp + Salinity + 
                  pH + NitrateNitrite + Silicate, 
                data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda7)$adj.r.squared # 0.1129493

anova.cca(dbrda7, step = 1000, permutations = 9999, by = "term")
#                Df SumOfSqs      F Pr(>F)   
# AvgSiteTemp     1  0.05318 1.3060 0.2090   
# Salinity        1  0.06486 1.5929 0.1344   
# pH              1  0.04459 1.0951 0.2938   
# NitrateNitrite  1  0.14880 3.6541 0.0074 **
# Silicate        1  0.04772 1.1719 0.2635   
# Residual       25  1.01802  

anova.cca(dbrda7, step = 1000, permutations = 9999) # p = 0.0277

# Remove pH
dbrda8 <- dbrda(wu_dist_carb_hel ~ AvgSiteTemp + Salinity + NitrateNitrite + Silicate, 
                data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda8)$adj.r.squared # 0.09642359

anova.cca(dbrda8, step = 1000, permutations = 9999, by = "term")

# Remove AvgSiteTemp
dbrda9 <- dbrda(wu_dist_carb_hel ~ Salinity + NitrateNitrite + Silicate, 
                data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda9)$adj.r.squared # 0.07289754

anova.cca(dbrda9, step = 1000, permutations = 9999, by = "term")

# Remove Salinity
dbrda10 <- dbrda(wu_dist_carb_hel ~ NitrateNitrite + Silicate, 
                 data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda10)$adj.r.squared # 0.04702462

anova.cca(dbrda10, step = 1000, permutations = 9999, by = "term")

# Remove Nitnit
dbrda11 <- dbrda(wu_dist_carb_hel ~ Silicate, 
                 data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda11)$adj.r.squared # 0.04475339

anova.cca(dbrda11, step = 1000, permutations = 9999, by = "term")

# Perform backwards selection for carbonate communities with unweighted UniFrac

full_dbrda <- dbrda(uu_dist_carb_hel ~ ., data = sample_data_carb, na.action = na.exclude)

RsquareAdj(full_dbrda)$adj.r.squared # 0.0003634317

anova.cca(full_dbrda, step = 1000, permutations = 9999, by = "term")

# Remove Ammonium 
dbrda1 <- dbrda(uu_dist_carb_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + NitrateNitrite + Silicate + MacroCover +
                  CoralCover + CoralSimp + MacroShan, 
                data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda1)$adj.r.squared # -0.00120802

anova.cca(dbrda1, step = 1000, permutations = 9999, by = "term")

# Remove  CoralSimp
dbrda2 <- dbrda(uu_dist_carb_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + NitrateNitrite + Silicate + MacroCover +
                  CoralCover + MacroShan, 
                data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda2)$adj.r.squared # 0.001733366

anova.cca(dbrda2, step = 1000, permutations = 9999, by = "term")

# Remove MacroShan
dbrda3 <- dbrda(uu_dist_carb_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + NitrateNitrite + Silicate + MacroCover + CoralCover, 
                data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda3)$adj.r.squared # 0.0100902

anova.cca(dbrda3, step = 1000, permutations = 9999, by = "term")

# Remove silicate
dbrda4 <- dbrda(uu_dist_carb_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + NitrateNitrite + MacroCover + CoralCover, 
                data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda4)$adj.r.squared # 0.009921884

anova.cca(dbrda4, step = 1000, permutations = 9999, by = "term")

# Remove CoralCover (BEST MODEL)
dbrda5 <- dbrda(uu_dist_carb_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  pH + NitrateNitrite + MacroCover, 
                data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda5)$adj.r.squared # 0.01574098

anova.cca(dbrda5, step = 1000, permutations = 9999, by = "term")
#                Df SumOfSqs      F Pr(>F)  
# Depth           1   0.2527 0.9671 0.4920  
# AvgSampTemp     1   0.2404 0.9201 0.5659  
# AvgSiteTemp     1   0.3350 1.2824 0.0507 .
# Salinity        1   0.2714 1.0388 0.3023  
# pH              1   0.2388 0.9139 0.6871  
# NitrateNitrite  1   0.3733 1.4289 0.0163 *
# MacroCover      1   0.2426 0.9285 0.6380  
#Residual       23   6.0089               
anova.cca(dbrda5, step = 1000, permutations = 9999) # p = 0.1555

# Remove pH
dbrda6 <- dbrda(uu_dist_carb_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  NitrateNitrite + MacroCover, 
                data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda6)$adj.r.squared # 0.008308167

anova.cca(dbrda6, step = 1000, permutations = 9999, by = "term")

# Remove MacroCover
dbrda7 <- dbrda(uu_dist_carb_hel ~ 
                  Depth + AvgSampTemp + AvgSiteTemp + Salinity + 
                  NitrateNitrite, 
                data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda7)$adj.r.squared # 0.01155953

anova.cca(dbrda7, step = 1000, permutations = 9999, by = "term")

# Remove AvgSampTemp
dbrda8 <- dbrda(uu_dist_carb_hel ~ 
                  Depth + AvgSiteTemp + Salinity + 
                  NitrateNitrite, 
                data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda8)$adj.r.squared # 0.01436167

anova.cca(dbrda8, step = 1000, permutations = 9999, by = "term")

# Remove Depth
dbrda9 <- dbrda(uu_dist_carb_hel ~ 
                  AvgSiteTemp + Salinity + NitrateNitrite, 
                data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda9)$adj.r.squared # 0.01517665

anova.cca(dbrda9, step = 1000, permutations = 9999, by = "term")

# Remove Salinity
dbrda10 <- dbrda(uu_dist_carb_hel ~ 
                   AvgSiteTemp + NitrateNitrite, 
                 data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda10)$adj.r.squared # 0.005384302

anova.cca(dbrda10, step = 1000, permutations = 9999, by = "term")

# Remove nitnit
dbrda11 <- dbrda(uu_dist_carb_hel ~ AvgSiteTemp, 
                 data = sample_data_carb, na.action = na.exclude)

RsquareAdj(dbrda11)$adj.r.squared # 0.005667692

anova.cca(dbrda11, step = 1000, permutations = 9999, by = "term")

#---------CCA- and carbonate-associated biofilms among sites---------

# Remove samples with <2,500 reads

CCA_rare <- subset_samples(CCA, sample_sums(CCA) > 2500)

# Rarefy CCA communities to lowest read depth
CCA_rare <- 
  rarefy_even_depth(
    CCA_rare, 
    sample.size = min(sample_sums(CCA_rare)), 
    rngseed = 8, 
    replace = FALSE)

# Calculate diversity indices

CCA_rare_div <- adiv_df(CCA_rare)

# Remove sites with <3 samples
CCA_rare_div_filt <-
  CCA_rare_div %>% 
  filter(!(Site %in% 
             (CCA_rare_div %>% 
                group_by(Site) %>% 
                transmute(n = n()) %>% 
                distinct() %>% 
                filter(n < 3) %>%
                pull(Site) %>% 
                as.character())))

# Determine which statistical tests to use

test_3groups(CCA_rare_div_filt, "Shannon", "Site")
test_3groups(CCA_rare_div_filt, "Simpson", "Site")

# Run tests to determine difference in CCA alpha diversity between sites

summary(aov(Shannon ~ Site, CCA_rare_div)) # p = 0.347
summary(aov(Simpson ~ Site, CCA_rare_div)) # p = 0.19

# Pull CCA sample data

CCA_df <- 
  data.frame(sample_data(CCA)) %>% 
  rownames_to_column("SampleID")

# Calculate unweighted and weighted Unifrac distances

uu_dist_CCA <- UniFrac(CCA, weighted = FALSE, normalized = TRUE)
wu_dist_CCA <- UniFrac(CCA, weighted = TRUE, normalized = TRUE)

# Run tests to determine difference in CCA dispersion between sites

permutest(betadisper(uu_dist_CCA, CCA_df$Site, bias.adjust = TRUE)) # p = 0.896
permutest(betadisper(wu_dist_CCA, CCA_df$Site, bias.adjust = TRUE)) # p = 0.837

# Run tests to determine difference in CCA composition between sites

adonis2(uu_dist_CCA ~ Site, data = CCA_df) # p = 0.038
adonis2(wu_dist_CCA ~ Site, data = CCA_df) # p = 0.01

# Pairwise comparisons

pairwiseAdonis::pairwise.adonis2(uu_dist_CCA ~ Site, data = CCA_df) # None significant!
pairwiseAdonis::pairwise.adonis2(wu_dist_CCA ~ Site, data = CCA_df) # None significant!

# Remove samples with <2,500 reads

carb_rare <- subset_samples(carb, sample_sums(carb) > 2500)

# Rarefy carbonate communities to lowest read depth
carb_rare <- 
  rarefy_even_depth(
    carb_rare, 
    sample.size = min(sample_sums(carb_rare)), 
    rngseed = 8, 
    replace = FALSE)

# Calculate diversity indices

carb_rare_div <- adiv_df(carb_rare)

# Remove sites with <3 samples
carb_rare_div_filt <-
  carb_rare_div %>% 
  filter(!(Site %in% 
             (carb_rare_div %>% 
                group_by(Site) %>% 
                transmute(n = n()) %>% 
                distinct() %>% 
                filter(n < 3) %>%
                pull(Site) %>% 
                as.character())))

# Determine which statistical tests to use

test_3groups(carb_rare_div_filt, "Shannon", "Site")
test_3groups(carb_rare_div_filt, "Simpson", "Site")

# Run tests to determine difference in carbonate alpha diversity between sites

summary(aov(Shannon ~ Site, carb_rare_div)) # p = 0.564
summary(aov(Simpson ~ Site, carb_rare_div)) # p = 0.369

# Pull carbonate sample data

carb_df <- 
  data.frame(sample_data(carb)) %>% 
  rownames_to_column("SampleID")

# Calculate unweighted and weighted Unifrac distances

uu_dist_carb <- UniFrac(carb, weighted = FALSE, normalized = TRUE)
wu_dist_carb <- UniFrac(carb, weighted = TRUE, normalized = TRUE)

# Run tests to determine difference in CCA dispersion between sites

permutest(betadisper(uu_dist_carb, carb_df$Site, bias.adjust = TRUE)) # p = 0.902
permutest(betadisper(wu_dist_carb, carb_df$Site, bias.adjust = TRUE)) # p = 0.415

# Run tests to determine difference in CCA composition between sites

adonis2(uu_dist_carb ~ Site, data = carb_df) # p = 0.108
adonis2(wu_dist_carb ~ Site, data = carb_df) # p = 0.158

# Get coordinates per sample
geo_CCA <-
  data.frame(SampleID = uu_dist_CCA %>% labels()) %>% 
  left_join(
    sample_data(CCA) %>% 
      data.frame() %>% 
      rownames_to_column(var = "SampleID"), 
    by = "SampleID") %>% 
  select(SampleID, Site) %>% 
  left_join(geo %>% rownames_to_column(var = "Site"), by = "Site") %>% 
  column_to_rownames(var = "SampleID") %>% 
  select(Longitude, Latitude)

# Calculate geographic distance using Haversine distance
geo_dist_CCA <- 
  distm(geo_CCA, fun = distHaversine) %>% 
  as.dist()

# Run Mantel tests

# Do CCA communities co-vary with geographic distance?

mantel(uu_dist_CCA, geo_dist_CCA, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.0364, p = 0.2226

mantel(wu_dist_CCA, geo_dist_CCA, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.03553, p = 0.2315

# Get coordinates per sample
geo_carb <-
  data.frame(SampleID = uu_dist_carb %>% labels()) %>% 
  left_join(
    sample_data(carb) %>% 
      data.frame() %>% 
      rownames_to_column(var = "SampleID"), 
    by = "SampleID") %>% 
  select(SampleID, Site) %>% 
  left_join(geo %>% rownames_to_column(var = "Site"), by = "Site") %>% 
  column_to_rownames(var = "SampleID") %>% 
  select(Longitude, Latitude)

# Calculate geographic distance using Haversine distance
geo_dist_carb <- 
  distm(geo_carb, fun = distHaversine) %>% 
  as.dist()

# Run Mantel tests

# Do carbonate communities co-vary with geographic distance?

mantel(uu_dist_carb, geo_dist_carb, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.00506, p = 0.5112

mantel(wu_dist_carb, geo_dist_carb, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.01298, p = 0.366

#---------Processes structuring biofilm community and sub-community assembly---------

# Need a phyloseq object grouped at the genus level

# Remove unidentified genera

bio_g <- subset_taxa(bio, Genus != "NA")

# Get taxonomy table

tax_bio_g <- tax_table(bio_g) %>% data.frame()

# Change "uncultured" and "unknown" Genera to be more specific

Genus_fixed <- rep("NA", length(tax_bio_g$Genus))

for(i in 1:length(tax_bio_g$Genus)) {
  
  Genus_fixed[i] <- 
    ifelse(tax_bio_g$Genus[i] == "Unknown_Family", 
           print(paste0("Unknown_", tax_bio_g$Order[i])), 
           print(tax_bio_g$Genus[i]))
  
  Genus_fixed[i] <-
    ifelse(tax_bio_g$Genus[i] == "uncultured",
           print(paste0("Uncultured_", tax_bio_g$Family[i])),
           print(Genus_fixed[i]))
  
  Genus_fixed[i] <-
    ifelse(Genus_fixed[i] == "Uncultured_uncultured",
           print(paste0("Uncultured_", tax_bio_g$Order[i])),
           print(Genus_fixed[i]))
  
  Genus_fixed[i] <-
    ifelse(Genus_fixed[i] == "Uncultured_uncultured",
           print(paste0("Uncultured_", tax_bio_g$Class[i])),
           print(Genus_fixed[i]))
  
  Genus_fixed[i] <-
    ifelse(Genus_fixed[i] == "Uncultured_uncultured",
           print(paste0("Uncultured_", tax_bio_g$Phylum[i])),
           print(Genus_fixed[i]))
  
}

# Replace old Genus names with more specific ones

tax_bio_g <- tax_bio_g %>% mutate(Genus = Genus_fixed)

# Clean up a few more names

tax_bio_g <- 
  tax_bio_g %>% 
  mutate(
    Genus =
      case_when(
        Genus == "Uncultured_Unknown_Family" & Class == "Gammaproteobacteria" ~ "Unknown_Gammaproteobacteria",
        Genus == "Uncultured_Unknown_Family" & Class == "Alphaproteobacteria" ~ "Unknown_Alphaproteobacteria",
        Genus == "Unknown_Gammaproteobacteria_Incertae_Sedis" ~ "Unknown_Gammaproteobacteria",
        Genus != "Uncultured_Unknown_Family" ~ Genus))

# Replace taxonomy table with more specific names

tax_table(bio_g) <- tax_bio_g %>% as.matrix()

# Group biofilm dataset at Genus level

bio_g <- tax_glom(bio_g, taxrank = "Genus")
bio_g <- prune_taxa(taxa_sums(bio_g) > 0, bio_g)

# Make inputs for iCAMP function
# See iCAMP documentation for information about inputs and outputs of these functions
# Takes significant computational power, so necessary outputs to recreate analyses are provided
# Skip to line 1621 for next analysis steps

# Community matrix with sample names as rownames and genera as colnames
comm <- 
  otu_table(bio_g) %>% 
  t() %>% 
  data.frame()

comm_names <- 
  colnames(comm) %>% 
  str_replace("X", "")

colnames(comm) <- comm_names

# Phylogenetic tree in newick format
tree <- phy_tree(bio_g)

# Classification table = taxonomy table with genera as rownames
clas <- 
  tax_table(bio_g) %>% 
  data.frame() %>% 
  select(-Confidence)

# Treatment data table = linking samples with site and substrate
treat <- 
  sample_data(bio_g) %>% 
  data.frame() %>% 
  select(Site, Substrate)

# Environmental data table = linking samples with environmental characteristics
env <- 
  sample_data(bio_g) %>% 
  data.frame() %>% 
  select(Depth:MacroSimp)

# Set parameters
prefix = "KARMA"  # prefix of the output file names. usually use a project ID
rand.time = 1000  # randomization time, 1000 is usually enough
memory.G = 140 # to set the memory size as you need. unit is GB

# Calculate pairwise phylogenetic distance matrix
cat("Calculating pairwise phylogenetic distances\n")

pd.big <- 
  pdist.big(tree = tree, wd = "iCAMP_outputs", 
            nworker = nworker, memory.G = memory.G)

# Assess niche preference difference between species
cat("Assessing niche preferences between species\n")

niche.dif <- 
  dniche(env = env, comm = comm, method = "niche.value",
         nworker = nworker, out.dist=FALSE, bigmemo=TRUE,
         nd.wd="iCAMP_outputs2")

# Phylogenetic binning
cat("Performing phylogenetic binning\n")

# Converts tree to be compatible with bigmemory
tree.rt <- 
  midpoint.root.big(
    tree = tree, 
    pd.desc = pd.big$pd.file,
    pd.spname = pd.big$tip.label, 
    pd.wd = pd.big$pd.wd,
    nworker = nworker)

# Replaces tree with bigmemory tree
tree <- tree.rt$tree

# Set parameters
ds = 0.2
bin.size.limit = 2

# Binning
phylobin <-
  taxa.binphy.big(
    tree = tree, 
    pd.desc = pd.big$pd.file, 
    pd.spname = pd.big$tip.label,
    pd.wd = pd.big$pd.wd, 
    ds = ds, 
    bin.size.limit = bin.size.limit,
    nworker = nworker)

# Test within-bin phylogenetic signal
cat("Testing within-bin phylogenetic signal\n")

# Prep data for test
sp.bin <- phylobin$sp.bin[,3,drop=FALSE]
sp.ra <- colMeans(comm/rowSums(comm))
abcut = 1 # you may remove some species if they are too rare to perform reliable correlation test.
commc = comm[,colSums(comm) >= abcut, drop = FALSE]
dim(commc)
spname.use = colnames(commc)

# Mantel test to evaluate correlation between phylogenetic distance + niche difference
binps <- 
  ps.bin(
    sp.bin = sp.bin,
    sp.ra = sp.ra,
    spname.use = spname.use,
    pd.desc = pd.big$pd.file, 
    pd.spname = pd.big$tip.label, 
    pd.wd = pd.big$pd.wd,
    nd.list = niche.dif$nd,
    nd.spname = niche.dif$names,
    ndbig.wd = niche.dif$nd.wd,
    cor.method = "pearson",
    r.cut = 0.1,
    p.cut = 0.05, 
    min.spn = 5)

# Save these results
write.table(
  data.frame(ds = ds, n.min = bin.size.limit, binps$Index),
  file = paste0("iCAMP_outputs/", prefix,".PhyloSignalSummary.csv"), 
  append = FALSE, 
  quote = FALSE, 
  sep=",", 
  row.names = FALSE,
  col.names = TRUE)

write.table(
  data.frame(ds = ds, n.min = bin.size.limit, binID = rownames(binps$detail), binps$detail),
  file = paste0("iCAMP_outputs/", prefix,".PhyloSignalDetail.csv"),
  append = FALSE, 
  quote = FALSE, 
  sep = ",", 
  row.names = FALSE, 
  col.names = TRUE)

# iCAMP analysis without omitting small bins (merges them with nearby larger bins)
cat("Performing iCAMP analysis\n")

# Set parameters
sig.index = "SES.RC" # This is the "traditional" way. Check normality afterward

# Run iCAMP
icres <-
  icamp.big(
    comm = comm,
    pd.desc = pd.big$pd.file,
    pd.spname=pd.big$tip.label,
    pd.wd = pd.big$pd.wd,
    rand = rand.time,
    tree = tree,
    prefix = prefix,
    ds = ds, 
    pd.cut = NA, 
    sp.check = TRUE,
    phylo.rand.scale = "within.bin", 
    taxa.rand.scale = "across.all",
    phylo.metric = "bMNTD", 
    sig.index = sig.index, 
    bin.size.limit = bin.size.limit, 
    nworker = nworker, 
    memory.G = memory.G, 
    rtree.save = FALSE, 
    detail.save = TRUE, 
    qp.save = TRUE, 
    detail.null = TRUE, 
    ignore.zero = TRUE, 
    output.wd = "iCAMP_outputs", 
    correct.special = TRUE, 
    unit.sum = rowSums(comm), 
    special.method = "depend",
    ses.cut = 1.96, 
    rc.cut = 0.95, 
    conf.cut = 0.975, 
    omit.option = "no",
    meta.ab = NULL)

# Normality test
nntest <- null.norm(icamp.output = icres, p.norm.cut = 0.05, detail.out = FALSE)
# If some ratio values are very high, may need to change to use "Confidence" as sig.index.
saveRDS(nntest, "iCAMP_outputs/nntest.rds")
nntest$summary %>% View()

# Get iCAMP bin-level statistics
cat("Getting iCAMP bin-level statistics\n")

# Calculate bin statistics
icbin <-
  icamp.bins(
    icamp.detail = icres$detail,
    treat = treat,
    clas = clas,
    silent=FALSE, 
    boot = TRUE,
    rand.time = rand.time,
    between.group = TRUE)

# Save these statistics!
save(icbin, file = paste0("iCAMP_outputs/", prefix,".iCAMP.Summary.rda"))
write.csv(icbin$Pt, 
          file = paste0("iCAMP_outputs/", prefix,".ProcessImportance_EachGroup.csv"), 
          row.names = FALSE)
write.csv(icbin$Ptk, 
          file = paste0("iCAMP_outputs/", prefix,".ProcessImportance_EachBin_EachGroup.csv"), 
          row.names = FALSE)
write.csv(icbin$Ptuv, 
          file = paste0("iCAMP_outputs/", prefix,".ProcessImportance_EachTurnover.csv"), 
          row.names = FALSE)
write.csv(icbin$BPtk, 
          file = paste0("iCAMP_outputs/", prefix,".BinContributeToProcess_EachGroup.csv"), 
          row.names = FALSE)
write.csv(
  data.frame(
    ID = rownames(icbin$Class.Bin),
    icbin$Class.Bin,
    stringsAsFactors = FALSE),
  file = paste0("iCAMP_outputs/", prefix,".Taxon_Bin.csv"),
  row.names = FALSE)
write.csv(icbin$Bin.TopClass, 
          file = paste0("iCAMP_outputs/", prefix,".Bin_TopTaxon.csv"), 
          row.names = FALSE)

# Bootstrapping test by substrate
cat("Performing bootstrapping test (Substrate)\n")

# Prepare data
i = 2 # Column with substrate info
treat.use = treat[,i,drop=FALSE]

# Bootstrapping
icboot_sub <- 
  icamp.boot(
    icamp.result = icres$bNTIiRCa,
    treat = treat.use,
    rand.time = rand.time,
    compare = TRUE,
    silent = FALSE,
    between.group = TRUE,
    ST.estimation = TRUE)

# Save bootstrapping results
save(icboot_sub, file = paste0("iCAMP_outputs/", prefix, ".iCAMP.Boot.", colnames(treat)[i], ".rda"))
write.csv(icboot_sub$summary,
          file = paste0("iCAMP_outputs/", prefix, ".iCAMP.BootSummary.", colnames(treat)[i],".csv"),
          row.names = FALSE)
write.csv(icboot_sub$compare,
          file = paste0("iCAMP_outputs/", prefix,".iCAMP.Compare.",colnames(treat)[i],".csv"),
          row.names = FALSE)

# Bootstrapping test by site
cat("Performing bootstrapping test (Site)\n")

# Prepare data
i = 1 # Column with substrate info
treat.use = treat[,i,drop=FALSE]

# Bootstrapping
icboot_site <- 
  icamp.boot(
    icamp.result = icres$bNTIiRCa,
    treat = treat.use,
    rand.time = rand.time,
    compare = TRUE,
    silent = FALSE,
    between.group = TRUE,
    ST.estimation = TRUE)

# Save bootstrapping results
save(icboot_site, file = paste0("iCAMP_outputs/", prefix, ".iCAMP.Boot.", colnames(treat)[i], ".rda"))
write.csv(icboot_site$summary,
          file = paste0("iCAMP_outputs/", prefix, ".iCAMP.BootSummary.", colnames(treat)[i],".csv"),
          row.names = FALSE)
write.csv(icboot_site$compare,
          file = paste0("iCAMP_outputs/", prefix,".iCAMP.Compare.",colnames(treat)[i],".csv"),
          row.names = FALSE)

###END OF iCAMP ANALYSIS###

# What is the relative influence of each assembly process in both substrates?

read_csv("data/KARMA.iCAMP.BootSummary.Substrate.csv") %>% 
  select(Substrate = Group, Process, Mean) %>% 
  filter(Substrate != "CCA_vs_Carbonate") %>% 
  filter(Process != "Stochasticity") %>% 
  mutate(
    Process = 
      case_when(
        Process == "Heterogeneous.Selection" ~ "HeS",
        Process == "Homogeneous.Selection" ~ "HoS",
        Process == "Dispersal.Limitation" ~ "DL",
        Process == "Homogenizing.Dispersal" ~ "HD",
        Process == "Drift.and.Others" ~ "DR")) %>% 
  pivot_wider(names_from = "Process", values_from = "Mean")


# Which assembly processes have significantly different relative influence between substrates?
# Table of p-values

read_csv("data/KARMA.iCAMP.Compare.Substrate.csv") %>% 
  filter(!(str_detect(Group2, "_vs_"))) %>% 
  select(HeS = Heterogeneous.Selection_P.value,
         HoS = Homogeneous.Selection_P.value,
         HD = Homogenizing.Dispersal_P.value,
         DL = Dispersal.Limitation_P.value,
         DR = Drift.and.Others_P.value)

# Identify habitat generalists and habitat specialists

message("There are ", nrow(tax_table(bio_g)), " genera")

# Get otu table

otu_bio_g <- 
  otu_table(bio_g) %>% 
  t() %>% 
  data.frame()

# Edit names

otu_names <- 
  colnames(otu_bio_g) %>% 
  str_replace("X", "")

colnames(otu_bio_g) <- otu_names

# Determine habitat generalists and specialists

genspec_g <- spec.gen(otu_bio_g, niche.width.method = "levins", n = 9999)

# Clean up data frame

genspec_g <- 
  genspec_g %>% 
  rownames_to_column(var = "SeqID") %>% 
  mutate(sign = str_replace(sign, "SPECIALIST", "Specialist"),
         sign = str_replace(sign, "NON SIGNIFICANT", "None"),
         sign = str_replace(sign, "GENERALIST", "Generalist")) %>% 
  rename(Sign = sign)

# Determine rare and abundant taxa

# Calculate relative abundance

bio_g_rel <- transform_sample_counts(bio_g, function(OTU) OTU/sum(OTU))
bio_g_rel <- prune_taxa(taxa_sums(bio_g_rel) > 0, bio_g_rel)

# Transform to data frame

bio_g_rel_df <- psmelt(bio_g_rel)

# Manipulate and assign abundant vs rare labels

abun_data <-
  bio_g_rel_df %>% 
  group_by(OTU) %>% 
  mutate(RelAbun = mean(Abundance)) %>% 
  select(SeqID = OTU, RelAbun) %>% 
  distinct() %>% 
  mutate(
    Category = 
      if_else(
        RelAbun < 0.00001, 
        "Rare", 
        if_else(
          RelAbun >= 0.001, 
          "Abundant", 
          "None"))) %>% 
  select(-RelAbun)

# How many habitat generalists vs specialists are there?

sub_comm_data <- 
  bio_g_rel_df %>% 
  group_by(OTU) %>% 
  transmute(RelAbun = mean(Abundance)) %>% 
  distinct() %>% 
  left_join(genspec_g,
            by = join_by("OTU" == "SeqID")) %>% 
  left_join(abun_data %>% 
              select(OTU = SeqID, Category),
            by = "OTU") %>% 
  rename(Abun_Rare = Category, Gen_Spec = Sign)

data.frame(
  Subcommunity = c("Generalist", "Specialist", "Abundant", "Rare"),
  NumGenera = c(
    genspec_g %>% filter(Sign == "Generalist") %>% nrow(),
    genspec_g %>% filter(Sign == "Specialist") %>% nrow(),
    abun_data %>% filter(Category == "Abundant") %>% nrow(),
    abun_data %>% filter(Category == "Rare") %>% nrow()),
  PercAbun = c(
    sub_comm_data %>% filter(Gen_Spec == "Generalist") %>% pull(RelAbun) %>% sum(),
    sub_comm_data %>% filter(Gen_Spec == "Specialist") %>% pull(RelAbun) %>% sum(),
    sub_comm_data %>% filter(Abun_Rare == "Abundant") %>% pull(RelAbun) %>% sum(),
    sub_comm_data %>% filter(Abun_Rare == "Rare") %>% pull(RelAbun) %>% sum()))

# Connect subcommunity identities to iCAMP bins
# From iCAMP pipeline, this can take a long time, so output files have been provided
# Skip to line 1782 for next analysis steps

# Look at processes shaping subcommunities: habitat generalists vs. specialists

cate <- 
  genspec_g %>% 
  column_to_rownames(var = "SeqID") %>% 
  select(Sign)

iccate_genspec <-
  icamp.cate(
    icamp.bins.result = icbin,
    comm = comm,
    cate = cate,
    treat = treat)

write.csv(iccate_genspec$Ptuvx,
          file = paste0("iCAMP_outputs/", prefix,".iCAMP.GenSpec.Process_EachTurnover_EachCategory.csv"))
write.csv(iccate_genspec$Ptx,
          file = paste0("iCAMP_outputs/", prefix,".iCAMP.GenSpec.Process_EachGroup_EachCategory.csv"))

# Look at processes shaping subcommunities: abundant vs. rare taxa

cate2 <- 
  abun_data %>% 
  column_to_rownames(var = "SeqID")

iccate_abunrare <-
  icamp.cate(
    icamp.bins.result = icbin,
    comm = comm,
    cate = cate2,
    treat = treat)

write.csv(iccate_abunrare$Ptuvx,
          file = paste0("iCAMP_outputs2/", prefix,".iCAMP.AbunRare.Process_EachTurnover_EachCategory.csv"))
write.csv(iccate_abunrare$Ptx,
          file = paste0("iCAMP_outputs2/", prefix,".iCAMP.AbunRare.Process_EachGroup_EachCategory.csv"))

# Read in iCAMP outputs for subcommunity analyses

genspec_ic <-
  read_csv("data/KARMA.iCAMP.GenSpec.Process_EachGroup_EachCategory.csv")

abunrare_ic <-
  read_csv("data/KARMA.iCAMP.AbunRare.Process_EachGroup_EachCategory.csv")

# What is the relative influence of each assembly process within each subcommunity on each substrate?

rbind(
  genspec_ic %>% 
    filter(GroupBasedOn == "Substrate") %>% 
    select(Substrate = Group,
           contains("Generalist") & !(contains("Stochasticity")),
           contains("Specialist") & !(contains("Stochasticity"))) %>% 
    pivot_longer(
      !Substrate, 
      names_to = c("Subcomm", "Process"), 
      names_sep = "\\.", 
      values_to = "RelInf") %>% 
    pivot_wider(names_from = "Process", values_from = "RelInf"),
  abunrare_ic %>% 
    filter(GroupBasedOn == "Substrate") %>% 
    select(Substrate = Group,
           contains("Abundant") & !(contains("Stochasticity")),
           contains("Rare") & !(contains("Stochasticity"))) %>% 
    pivot_longer(
      !Substrate, 
      names_to = c("Subcomm", "Process"), 
      names_sep = "\\.", 
      values_to = "RelInf") %>% 
    pivot_wider(names_from = "Process", values_from = "RelInf")) %>% 
  select(Subcomm, Substrate, HeS:DR) %>% 
  arrange(Subcomm)

#---------Environmental factors driving assembly processes in biofilm communities---------

# Does CCA-associated biofilm composition co-vary with bacterioplankton composition?

CCA_otu_inv <- CCA_otu %>% t() %>% data.frame()
water_otu_inv <- water_otu %>% t() %>% data.frame()

# Arrange CCA OTU table columns in order

CCA_otu_inv <-
  CCA_otu_inv %>% 
  select(all_of(colnames(CCA_otu_inv) %>% sort()))

# Need to triplicate each water sample to align with CCA samples

water_otu_inv <- cbind(water_otu_inv, water_otu_inv, water_otu_inv)

colnames(water_otu_inv) <- make.unique(colnames(water_otu_inv))

# Reorder columns so they align with CCA samples

water_otu_inv <-
  water_otu_inv %>% 
  select(all_of(colnames(water_otu_inv) %>% sort()))

# Rename water columns with CCA sample IDs

colnames(water_otu_inv) <- colnames(CCA_otu_inv)

# Invert again in preparation for Hellinger transformation

CCA_otu_ordered <- CCA_otu_inv %>% t() %>% data.frame()
water_otu_CCA <- water_otu_inv %>% t() %>% data.frame()

# Hellinger transform both datasets

water_otu_CCA.hel <- decostand(water_otu_CCA, method = "hellinger")
CCA_otu_ordered.hel <- decostand(CCA_otu_ordered, method = "hellinger")

# Build new ps objects

water_CCA_mantel <- 
  phyloseq(
    otu_table(water_otu_CCA.hel, taxa_are_rows = FALSE),
    tax_table(water),
    phy_tree(water))

CCA_mantel <- 
  phyloseq(
    otu_table(CCA_otu_ordered.hel, taxa_are_rows = FALSE),
    tax_table(CCA),
    phy_tree(CCA))

# Calculate unweighted and weighted Unifrac distances

uu_dist_wat_CCA_mantel <- UniFrac(water_CCA_mantel, weighted = FALSE, normalized = TRUE)
wu_dist_wat_CCA_mantel <- UniFrac(water_CCA_mantel, weighted = TRUE, normalized = TRUE)

uu_dist_CCA_mantel <- UniFrac(CCA_mantel, weighted = FALSE, normalized = TRUE)
wu_dist_CCA_mantel <- UniFrac(CCA_mantel, weighted = TRUE, normalized = TRUE)

# Does CCA-associated biofilm composition co-vary with bacterioplankton composition?

mantel(uu_dist_CCA_mantel, uu_dist_wat_CCA_mantel, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.05322, p = 0.2748

mantel(wu_dist_CCA_mantel, wu_dist_wat_CCA_mantel, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.02609, p = 0.5888

# Does CCA-associated biofilm composition co-vary with coral composition?

coral_cov_inv <- 
  coral_cov %>% 
  t() %>% 
  data.frame()

# Get triplicate coral community data
coral_cov_inv <- cbind(coral_cov_inv, coral_cov_inv, coral_cov_inv)

# Fix column names to match CCA samples

colnames(coral_cov_inv) <- make.unique(colnames(coral_cov_inv))

# Reorder columns so they align with CCA samples

coral_cov_inv <-
  coral_cov_inv %>% 
  select(Fish.Den, Fish.Den.1, Fish.Den.2, 
         Leo.s.Den, Leo.s.Den.1, Leo.s.Den.2,
         B.S.Deep, B.S.Deep.1, B.S.Deep.2,
         White.Hole, White.Hole.1, White.Hole.2,
         Wicked.Pissa, Wicked.Pissa.1, Wicked.Pissa.2,
         Phlipper.s.Peace, Phlipper.s.Peace.1, Phlipper.s.Peace.2,
         Angus..Aquarium, Angus..Aquarium.1, Angus..Aquarium.2,
         Thalassa.Arg, Thalassa.Arg.1, Thalassa.Arg.2,
         International.Diver, International.Diver.1, International.Diver.2,
         Pillar.Coral, Pillar.Coral.1, Pillar.Coral.2,
         Overheat, Overheat.1, Overheat.2)

# Rename coral columns with CCA sample IDs

colnames(coral_cov_inv) <- colnames(CCA_otu_inv)

# Invert coral cover data frame again to prepare for Hellinger transformation

coral_cov_CCA <- coral_cov_inv %>% t() %>% data.frame()

# Hellinger transform

coral_cov_CCA.hel <- decostand(coral_cov_CCA, method = "hellinger")

# Calculate Bray Curtis distances

coral_cov_CCA_dist <- vegdist(coral_cov_CCA.hel, method = "bray")

# Does CCA-associated biofilm composition co-vary with coral composition?

mantel(uu_dist_CCA_mantel, coral_cov_CCA_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.0004325, p = 0.4688

mantel(wu_dist_CCA_mantel, coral_cov_CCA_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.1112, p = 0.1386

# Does CCA-associated biofilm composition co-vary with algae composition?

algae_cov_inv <- 
  alg_cov %>% 
  t() %>% 
  data.frame()

# Get triplicate algae community data

algae_cov_inv <- cbind(algae_cov_inv, algae_cov_inv, algae_cov_inv)

# Fix column names to match CCA samples

colnames(algae_cov_inv) <- make.unique(colnames(algae_cov_inv))

# Reorder columns so they align with CCA samples

algae_cov_inv <-
  algae_cov_inv %>% 
  select(Fish.Den, Fish.Den.1, Fish.Den.2, 
         Leo.s.Den, Leo.s.Den.1, Leo.s.Den.2,
         B.S.Deep, B.S.Deep.1, B.S.Deep.2,
         White.Hole, White.Hole.1, White.Hole.2,
         Wicked.Pissa, Wicked.Pissa.1, Wicked.Pissa.2,
         Phlipper.s.Peace, Phlipper.s.Peace.1, Phlipper.s.Peace.2,
         Angus..Aquarium, Angus..Aquarium.1, Angus..Aquarium.2,
         Thalassa.Arg, Thalassa.Arg.1, Thalassa.Arg.2,
         International.Diver, International.Diver.1, International.Diver.2,
         Pillar.Coral, Pillar.Coral.1, Pillar.Coral.2,
         Overheat, Overheat.1, Overheat.2)

# Rename algae columns with CCA sample IDs

colnames(algae_cov_inv) <- colnames(CCA_otu_inv)

# Invert algae cover data frame again to prepare for Hellinger transformation

algae_cov_CCA <- algae_cov_inv %>% t() %>% data.frame()

# Hellinger transform

algae_cov_CCA.hel <- decostand(algae_cov_CCA, method = "hellinger")

# Calculate Bray Curtis distances

algae_cov_CCA_dist <- vegdist(algae_cov_CCA.hel, method = "bray")

# Does CCA-associated biofilm composition co-vary with algae composition?

mantel(uu_dist_CCA_mantel, algae_cov_CCA_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.0645, p = 0.2205

mantel(wu_dist_CCA_mantel, algae_cov_CCA_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.01985, p = 0.3886

# Does composition of selection-assembled subcommunity of CCA-associated biofilm co-vary with bacterioplankton composition?

# Get bins that are dominated by selection (HoS or HeS) in CCA

bin_CCA_selection <-
  read_csv("data/KARMA.ProcessImportance_EachBin_EachGroup.csv") %>% 
  filter(Group == "CCA") %>% 
  select(Index, contains("bin")) %>% 
  column_to_rownames(var = "Index") %>% 
  t() %>% 
  data.frame() %>% 
  select(DominantProcess) %>% 
  rownames_to_column(var = "Bin") %>% 
  filter(DominantProcess == "HeS" | DominantProcess == "HoS") %>% 
  mutate(Bin = str_replace(Bin, "b", "B")) %>% 
  pull(Bin)

# Get taxa associated with those bins

taxa_CCA_selection <- 
  read_csv("data/KARMA.Taxon_Bin.csv") %>% 
  select(ID, Bin) %>% 
  filter(Bin %in% bin_CCA_selection) %>% 
  pull(ID)

# Subset phyloseq object to only selection-dominated taxa

CCA_selection <- 
  prune_taxa(taxa_names(CCA_mantel) %in% taxa_CCA_selection, 
             CCA_mantel)

# Get UniFrac distances

uu_dist_CCA_selection <- UniFrac(CCA_selection, weighted = FALSE, normalized = TRUE)
wu_dist_CCA_selection <- UniFrac(CCA_selection, weighted = TRUE, normalized = TRUE)

# Does composition of selection-assembled subcommunity of CCA-associated biofilm co-vary with bacterioplankton composition?

mantel(uu_dist_CCA_selection, uu_dist_wat_CCA_mantel, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.1339, p = 0.0734

mantel(wu_dist_CCA_selection, wu_dist_wat_CCA_mantel, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.1271, p = 0.0766

# Does composition of selection-assembled subcommunity of CCA-associated biofilm co-vary with nearby coral composition?

mantel(uu_dist_CCA_selection, coral_cov_CCA_dist, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.182, p = 0.0337

mantel(wu_dist_CCA_selection, coral_cov_CCA_dist, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.3447, p = 8e-04

# Does composition of selection-assembled subcommunity of CCA-associated biofilm co-vary with nearby macroalgae composition?

mantel(uu_dist_CCA_selection, algae_cov_CCA_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.01133, p = 0.549

mantel(wu_dist_CCA_selection, algae_cov_CCA_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.06313, p = 0.7827

# Does composition of dispersal-assembled subcommunity of CCA-associated biofilm co-vary with bacterioplankton composition?

# Get bins that are dominated by dispersal (HD or DL) in CCA

bin_CCA_dispersal <-
  read_csv("data/KARMA.ProcessImportance_EachBin_EachGroup.csv") %>% 
  filter(Group == "CCA") %>% 
  select(Index, contains("bin")) %>% 
  column_to_rownames(var = "Index") %>% 
  t() %>% 
  data.frame() %>% 
  select(DominantProcess) %>% 
  rownames_to_column(var = "Bin") %>% 
  filter(DominantProcess == "HD" | DominantProcess == "DL") %>% 
  mutate(Bin = str_replace(Bin, "b", "B")) %>% 
  pull(Bin)

# Get taxa associated with those bins

taxa_CCA_dispersal <- 
  read_csv("data/KARMA.Taxon_Bin.csv") %>% 
  select(ID, Bin) %>% 
  filter(Bin %in% bin_CCA_dispersal) %>% 
  pull(ID)

# Subset phyloseq object to only selection-dominated taxa

CCA_dispersal <- 
  prune_taxa(taxa_names(CCA_mantel) %in% taxa_CCA_dispersal, 
             CCA_mantel)

# Get UniFrac distances

uu_dist_CCA_dispersal <- UniFrac(CCA_dispersal, weighted = FALSE, normalized = TRUE)
wu_dist_CCA_dispersal <- UniFrac(CCA_dispersal, weighted = TRUE, normalized = TRUE)

# Does composition of dispersal-assembled subcommunity of CCA-associated biofilm co-vary with bacterioplankton composition?

mantel(uu_dist_CCA_dispersal, uu_dist_wat_CCA_mantel, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.03039, p = 0.6162

mantel(wu_dist_CCA_dispersal, wu_dist_wat_CCA_mantel, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.04883 , p = 0.7084

# Does composition of dispersal-assembled subcommunity of CCA-associated biofilm co-vary with nearby coral composition?

mantel(uu_dist_CCA_dispersal, coral_cov_CCA_dist, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.07286, p = 0.7636

mantel(wu_dist_CCA_dispersal, coral_cov_CCA_dist, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.04106, p = 0.6482

# Does composition of dispersal-assembled subcommunity of CCA-associated biofilm co-vary with nearby macroalgae composition?

mantel(uu_dist_CCA_dispersal, algae_cov_CCA_dist, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.1307, p = 0.0567

mantel(wu_dist_CCA_dispersal, algae_cov_CCA_dist, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.08128, p = 0.1647

# Does carbonate-associated biofilm composition co-vary with bacterioplankton composition?

carb_otu_inv <- carb_otu %>% t() %>% data.frame()
water_otu_inv <- water_otu %>% t() %>% data.frame()

# Arrange carbonate OTU table columns in order

carb_otu_inv <-
  carb_otu_inv %>% 
  select(all_of(colnames(carb_otu_inv) %>% sort()))

# Need to triplicate each water sample to align with carbonate samples

water_otu_inv <- cbind(water_otu_inv, water_otu_inv, water_otu_inv)

colnames(water_otu_inv) <- make.unique(colnames(water_otu_inv))

# Reorder columns so they align with carbonate samples

water_otu_inv <-
  water_otu_inv %>% 
  select(all_of(colnames(water_otu_inv) %>% sort())) %>% 
  select(-K35.2, -K70.2)

# Rename water columns with carbonate sample IDs

colnames(water_otu_inv) <- colnames(carb_otu_inv)

# Invert again in preparation for Hellinger transformation

carb_otu_ordered <- carb_otu_inv %>% t() %>% data.frame()
water_otu_carb <- water_otu_inv %>% t() %>% data.frame()

# Hellinger transform both datasets

water_otu_carb.hel <- decostand(water_otu_carb, method = "hellinger")
carb_otu_ordered.hel <- decostand(carb_otu_ordered, method = "hellinger")

# Build new ps objects

water_carb_mantel <- 
  phyloseq(
    otu_table(water_otu_carb.hel, taxa_are_rows = FALSE),
    tax_table(water),
    phy_tree(water))

carb_mantel <- 
  phyloseq(
    otu_table(carb_otu_ordered.hel, taxa_are_rows = FALSE),
    tax_table(carb),
    phy_tree(carb))

# Calculate unweighted and weighted Unifrac distances

uu_dist_wat_carb_mantel <- UniFrac(water_carb_mantel, weighted = FALSE, normalized = TRUE)
wu_dist_wat_carb_mantel <- UniFrac(water_carb_mantel, weighted = TRUE, normalized = TRUE)

uu_dist_carb_mantel <- UniFrac(carb_mantel, weighted = FALSE, normalized = TRUE)
wu_dist_carb_mantel <- UniFrac(carb_mantel, weighted = TRUE, normalized = TRUE)

# Does carbonate-associated biofilm composition co-vary with bacterioplankton composition?

mantel(uu_dist_carb_mantel, uu_dist_wat_carb_mantel, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.2221, p = 0.0234

mantel(wu_dist_carb_mantel, wu_dist_wat_carb_mantel, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.005961 , p = 0.494

# Does carbonate-associated biofilm composition co-vary with coral composition?

coral_cov_inv <- 
  coral_cov %>% 
  t() %>% 
  data.frame()

# Get triplicate coral community data

coral_cov_inv <- cbind(coral_cov_inv, coral_cov_inv, coral_cov_inv)

# Fix column names to match CCA samples

colnames(coral_cov_inv) <- make.unique(colnames(coral_cov_inv))

# Reorder columns so they align with carbonate samples

coral_cov_inv <-
  coral_cov_inv %>% 
  select(Fish.Den, Fish.Den.1, Fish.Den.2, 
         Leo.s.Den, Leo.s.Den.1, Leo.s.Den.2,
         B.S.Deep, B.S.Deep.1, B.S.Deep.2,
         White.Hole, White.Hole.1, White.Hole.2,
         Wicked.Pissa, Wicked.Pissa.1,
         Phlipper.s.Peace, Phlipper.s.Peace.1, Phlipper.s.Peace.2,
         Angus..Aquarium, Angus..Aquarium.1, Angus..Aquarium.2,
         Thalassa.Arg, Thalassa.Arg.1, Thalassa.Arg.2,
         International.Diver, International.Diver.1, International.Diver.2,
         Pillar.Coral, Pillar.Coral.1,
         Overheat, Overheat.1, Overheat.2)

# Rename coral columns with carbonate sample IDs

colnames(coral_cov_inv) <- colnames(carb_otu_inv)

# Invert coral cover data frame again to prepare for Hellinger transformation

coral_cov_carb <- coral_cov_inv %>% t() %>% data.frame()

# Hellinger transform

coral_cov_carb.hel <- decostand(coral_cov_carb, method = "hellinger")

# Calculate Bray Curtis distances

coral_cov_carb_dist <- vegdist(coral_cov_carb.hel, method = "bray")

# Does CCA-associated biofilm composition co-vary with coral composition?

mantel(uu_dist_carb_mantel, coral_cov_carb_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.04028, p = 0.342 

mantel(wu_dist_carb_mantel, coral_cov_carb_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.05529, p = 0.6865

# Does carbonate-associated biofilm composition co-vary with macroalgae composition?

algae_cov_inv <- 
  alg_cov %>% 
  t() %>% 
  data.frame()

# Get triplicate algae community data

algae_cov_inv <- cbind(algae_cov_inv, algae_cov_inv, algae_cov_inv)

# Fix column names to match carbonate samples

colnames(algae_cov_inv) <- make.unique(colnames(algae_cov_inv))

# Reorder columns so they align with carbonate samples

algae_cov_inv <-
  algae_cov_inv %>% 
  select(Fish.Den, Fish.Den.1, Fish.Den.2, 
         Leo.s.Den, Leo.s.Den.1, Leo.s.Den.2,
         B.S.Deep, B.S.Deep.1, B.S.Deep.2,
         White.Hole, White.Hole.1, White.Hole.2,
         Wicked.Pissa, Wicked.Pissa.1,
         Phlipper.s.Peace, Phlipper.s.Peace.1, Phlipper.s.Peace.2,
         Angus..Aquarium, Angus..Aquarium.1, Angus..Aquarium.2,
         Thalassa.Arg, Thalassa.Arg.1, Thalassa.Arg.2,
         International.Diver, International.Diver.1, International.Diver.2,
         Pillar.Coral, Pillar.Coral.1,
         Overheat, Overheat.1, Overheat.2)

# Rename algae columns with carbonate sample IDs

colnames(algae_cov_inv) <- colnames(carb_otu_inv)

# Invert algae cover data frame again to prepare for Hellinger transformation

algae_cov_carb <- algae_cov_inv %>% t() %>% data.frame()

# Hellinger transform

algae_cov_carb.hel <- decostand(algae_cov_carb, method = "hellinger")

# Calculate Bray Curtis distances

algae_cov_carb_dist <- vegdist(algae_cov_carb.hel, method = "bray")

# Does carbonate-associated biofilm composition co-vary with macroalgae composition?

mantel(uu_dist_carb_mantel, algae_cov_carb_dist, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.1913, p = 0.9858

mantel(wu_dist_carb_mantel, algae_cov_carb_dist,
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.0738, p = 0.7925

# Does composition of selection-assembled subcommunity of carbonate-associated biofilm co-vary with bacterioplankton composition?

# Get bins that are dominated by selection (HoS or HeS) in carbonate

bin_carb_selection <-
  read_csv("data/KARMA.ProcessImportance_EachBin_EachGroup.csv") %>% 
  filter(Group == "Carbonate") %>% 
  select(Index, contains("bin")) %>% 
  column_to_rownames(var = "Index") %>% 
  t() %>% 
  data.frame() %>% 
  select(DominantProcess) %>% 
  rownames_to_column(var = "Bin") %>% 
  filter(DominantProcess == "HeS" | DominantProcess == "HoS") %>% 
  mutate(Bin = str_replace(Bin, "b", "B")) %>% 
  pull(Bin)

# Get taxa associated with those bins

taxa_carb_selection <- 
  read_csv("data/KARMA.Taxon_Bin.csv") %>% 
  select(ID, Bin) %>% 
  filter(Bin %in% bin_carb_selection) %>% 
  pull(ID)

# Subset phyloseq object to only selection-dominated taxa

carb_selection <- 
  prune_taxa(taxa_names(carb_mantel) %in% taxa_carb_selection, 
             carb_mantel)

# Get UniFrac distances

uu_dist_carb_selection <- UniFrac(carb_selection, weighted = FALSE, normalized = TRUE)
wu_dist_carb_selection <- UniFrac(carb_selection, weighted = TRUE, normalized = TRUE)

# Does composition of selection-assembled subcommunity of carbonate-associated biofilm co-vary with bacterioplankton composition?

mantel(uu_dist_carb_selection, uu_dist_wat_carb_mantel, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.03173, p = 0.367

mantel(wu_dist_carb_selection, wu_dist_wat_carb_mantel, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.04578, p = 0.6719

# Does composition of selection-assembled subcommunity of carbonate-associated biofilm co-vary with nearby coral composition?

mantel(uu_dist_carb_selection, coral_cov_carb_dist, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.1345, p = 0.1242

mantel(wu_dist_carb_selection, coral_cov_carb_dist, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.1656, p = 0.0635

# Does composition of selection-assembled subcommunity of carbonate-associated biofilm co-vary with nearby macroalgae composition?

mantel(uu_dist_carb_selection, algae_cov_carb_dist, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.2768, p = 0.9994

mantel(wu_dist_carb_selection, algae_cov_carb_dist, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.2244, p = 0.9979

# Does composition of dispersal-assembled subcommunity of carbonate-associated biofilm co-vary with bacterioplankton composition?

# Get bins that are dominated by dispersal (HD or DL) in carbonate

bin_carb_dispersal <-
  read_csv("data/KARMA.ProcessImportance_EachBin_EachGroup.csv") %>% 
  filter(Group == "Carbonate") %>% 
  select(Index, contains("bin")) %>% 
  column_to_rownames(var = "Index") %>% 
  t() %>% 
  data.frame() %>% 
  select(DominantProcess) %>% 
  rownames_to_column(var = "Bin") %>% 
  filter(DominantProcess == "HD" | DominantProcess == "DL") %>% 
  mutate(Bin = str_replace(Bin, "b", "B")) %>% 
  pull(Bin)

# Get taxa associated with those bins

taxa_carb_dispersal <- 
  read_csv("data/KARMA.Taxon_Bin.csv") %>% 
  select(ID, Bin) %>% 
  filter(Bin %in% bin_carb_dispersal) %>% 
  pull(ID)

# Subset phyloseq object to only selection-dominated taxa

carb_dispersal <- 
  prune_taxa(taxa_names(carb_mantel) %in% taxa_carb_dispersal, 
             carb_mantel)

# Get UniFrac distances

uu_dist_carb_dispersal <- UniFrac(carb_dispersal, weighted = FALSE, normalized = TRUE)
wu_dist_carb_dispersal <- UniFrac(carb_dispersal, weighted = TRUE, normalized = TRUE)

# Does composition of dispersal-assembled subcommunity of carbonate-associated biofilm co-vary with bacterioplankton composition?

mantel(uu_dist_carb_dispersal, uu_dist_wat_carb_mantel, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.2003, p = 0.0241

mantel(wu_dist_carb_dispersal, wu_dist_wat_carb_mantel, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.0402, p = 0.3141

# Does composition of dispersal-assembled subcommunity of carbonate-associated biofilm co-vary with nearby coral composition?

mantel(uu_dist_carb_dispersal, coral_cov_carb_dist, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.08651, p = 0.2007

mantel(wu_dist_carb_dispersal, coral_cov_carb_dist, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = 0.1819, p = 0.0522

# Does composition of dispersal-assembled subcommunity of carbonate-associated biofilm co-vary with nearby macroalgae composition?

mantel(uu_dist_carb_dispersal, algae_cov_carb_dist, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.1228, p = 0.928

mantel(wu_dist_carb_dispersal, algae_cov_carb_dist, 
       method = "spearman", permutations = 9999, na.rm = TRUE)
# r = -0.0937, p = 0.8398 

# Maaslin2 analysis to look for correlation between environmental conditions and iCAMP bins

# Group CCA community at genus level

CCA_g <- subset_samples(bio_g, Substrate == "CCA")
CCA_g <- prune_taxa(taxa_sums(CCA_g) > 0, CCA_g)

# Get bin-taxon key

bin_taxa_key <- 
  read_csv("data/KARMA.Taxon_Bin.csv") %>% 
  select(ID, Bin)

# Add bin to tax table and put back into phyloseq object

CCA_bin <- CCA_g

tax_table(CCA_bin) <-
  tax_table(CCA_g) %>% 
  data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  select(-Confidence) %>% 
  left_join(bin_taxa_key, by = "ID") %>% 
  column_to_rownames(var = "ID") %>% 
  select(Kingdom, Bin, Phylum:Genus) %>% 
  as.matrix()

# Group CCA community by bin

CCA_bin <- tax_glom(CCA_bin, taxrank = "Bin")
CCA_bin <- prune_taxa(taxa_sums(CCA_bin) > 0,  CCA_bin)

#Maaslin2 analysis

# Create Maaslin2 inputs

maaslin_data_CCA <- 
  otu_table(CCA_bin) %>% 
  data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  left_join(bin_taxa_key, by = "ID") %>% 
  select(-ID) %>% 
  column_to_rownames(var = "Bin") %>% 
  t() %>% 
  data.frame()

maaslin_metadata_CCA <- 
  sample_data(CCA_bin) %>% 
  data.frame()

# Run Maaslin2 for CCA

full_model_CCA <-
  Maaslin2(input_data = maaslin_data_CCA, 
           input_metadata = maaslin_metadata_CCA, 
           min_prevalence = 0,
           analysis_method = "NEGBIN",
           transform = "NONE",
           normalization = "TMM",
           output = "data/Maaslin2_full_model_CCA", 
           fixed_effects = c("Depth", "AvgSampTemp", "AvgSiteTemp", "Salinity", "pH",
                             "NitrateNitrite", "Ammonium", "Silicate", "CoralCover", 
                             "MacroCover", "CoralSimp", "MacroShan"),
           max_significance = 0.05)

# Group carbonate community at genus level

carb_g <- subset_samples(bio_g, Substrate == "Carbonate")
carb_g <- prune_taxa(taxa_sums(carb_g) > 0, carb_g)

# Add bin to tax table and put back into phyloseq object

carb_bin <- carb_g

tax_table(carb_bin) <-
  tax_table(carb_g) %>% 
  data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  select(-Confidence) %>% 
  left_join(bin_taxa_key, by = "ID") %>% 
  column_to_rownames(var = "ID") %>% 
  select(Kingdom, Bin, Phylum:Genus) %>% 
  as.matrix()

# Group carbonate community by bin

carb_bin <- tax_glom(carb_bin, taxrank = "Bin")
carb_bin <- prune_taxa(taxa_sums(carb_bin) > 0,  carb_bin)

#Maaslin2 analysis

# Create Maaslin2 inputs

maaslin_data_carb <- 
  otu_table(carb_bin) %>% 
  data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  left_join(bin_taxa_key, by = "ID") %>% 
  select(-ID) %>% 
  column_to_rownames(var = "Bin") %>% 
  t() %>% 
  data.frame()

maaslin_metadata_carb <- 
  sample_data(carb_bin) %>% 
  data.frame()

# Run Maaslin2 for carbonate

full_model_carb <-
  Maaslin2(input_data = maaslin_data_carb, 
           input_metadata = maaslin_metadata_carb, 
           min_prevalence = 0,
           analysis_method = "NEGBIN",
           transform = "NONE",
           normalization = "TMM",
           output = "data/Maaslin2_full_model_carb", 
           fixed_effects = c("Depth", "AvgSampTemp", "AvgSiteTemp", "Salinity", "pH",
                             "NitrateNitrite", "Ammonium", "Silicate", "CoralCover", 
                             "MacroCover", "CoralSimp", "MacroShan"),
           max_significance = 0.05)
