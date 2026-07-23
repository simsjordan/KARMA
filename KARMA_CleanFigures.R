# KARMA Clean Figures Script
# Jordan Sims
# July 2026

# This script contains the code needed to produce all the figures from the KARMA paper.

# INPUTS

# KARMA_V1-filt-tree.rds: a phyloseq object containing the biofilm bacteria and bacterioplankton 16S dataset that has been filtered to remove unassigned taxa, mitochondrial and chloroplast DNA, and taxa found in the negative controls. Negative controls and samples with too few reads have already been removed from the dataset.

# KARMA_environmental_data.csv: CSV file containing all environmental data collected at each site, including metadata, temperature metrics, pH, salinity, nutrient concentrations, benthic community data, coral community data, and macroalgae community data.

# KARMA.xxx.csv: CSV file containing an output produced by the iCAMP function. These outputs can be produced by following the analysis found in KARMA_CleanCode.R. See the iCAMP documentation for details on each file.

# abundant-rare-assignments-g.csv: CSV file identifying each sequence ID (grouped by Genus) as abundant, rare, or neither. This file is an output found in KARMA_CleanCode.R

# generalist-specialist-assignments-g.csv: CSV file containing the niche width of each sequence ID and identifying each sequence ID (grouped by Genus) as a habitat generalist, habitat specialist, or neither. This file is an output found in KARMA_CleanCode.R

# Maaslin2_xxx_significant_results.tsv: TSV file containing an output produced by the Maaslin2 analysis, one for carbonate and one for CCA substrates. These outputs can be produced by following the analysis found in KARMA_CleanCode.R. See the Maaslin2 documentation for details.

# hnd_admbnda_adm2_sinit_20161005.xxx: five map shape files, should be in a directory called "Honduras_Shape"

#---------Read in required packages and data---------

# Packages

library('phyloseq'); packageVersion('phyloseq') # 1.52.0
library('sf'); packageVersion('sf') # 1.0.22
library('tidyverse'); packageVersion('tidyverse') # 2.0.0
library('ggpubr'); packageVersion('ggpubr') # 0.6.2

# Custom functions

# Builds a data frame of Shannon and Simpson diversity values and sample metadata from a phyloseq object

adiv_df <- function(ps) {
  estimate_richness(ps, measures = c("Shannon", "Simpson")) %>% 
    rownames_to_column(var = "Sample") %>% 
    merge(rownames_to_column(data.frame(sample_data(ps)), var = "Sample"), by = "Sample")
}

# Connects points in NMDS plot for shading

find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]

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

# Set color assignments so they are consistent throughout figures

# Reef sites

site_color = 
  c("#f6735f", "#93b7d5", "#165ba9", "#f9ef70", "#ac6c29", "#6c5698", 
    "#faebaa", "#a6be48", "#582a36", "#ab4d56", "#009398")

# Sample type

biofilm_color = c("darkolivegreen3", "rosybrown2")
substrate_color = c("darkolivegreen3", "rosybrown2", "paleturquoise3")

# Subcommunities

gen_spec_color <- c("#96B2F2", "#EDC7DC", "grey80")
abun_rare_color <- c("lightblue", "yellow", "grey80")

# Assembly processes

process_color = c("#173F5F", "#F6D55C", "#3CAEA3", "#ED553B", "#c3d3e0")

# Composition plot colors

coral_color = c("grey80", "thistle3", "coral", "goldenrod1", "darkseagreen4", 
              "darkorchid4", "azure2", "dodgerblue2", "lightgreen", 
              "darkblue", "palevioletred2", "lemonchiffon2", "green4", 
              "darkred", "lightpink", "lightskyblue")

algae_color = c("lavenderblush2", "palevioletred4", "olivedrab4", 
              "chartreuse3", "goldenrod4", "deepskyblue2")

# Assign 16S phyla colors

Other <- c("grey90")
Acido <- c("thistle3")
Actino <- c("coral")
Bactero <- c("goldenrod1")
Bdello <- c("darkblue")
Camp <- c("yellowgreen")
Chloro <- c("mediumpurple1")
Cren <- c("dodgerblue2")
Cyano <- c("lightgreen")
Dada <- c("lightgoldenrod1")
Deino <- c("lightcyan1")
Deferr <- c("lightgoldenrod1")
Depend <- c("lightcyan1")
Desulfo <- c("red")
Ento <- c("palevioletred2")
Firm <- c("lightseagreen")
Fuso <- c("slategray3")
Gemma <- c("darkorchid4")
Lates <- c("aquamarine1")
Marin <- c("darkorchid4")
MBNT <- c("red")
Modul <- c("orchid2")
Myxo <- c("darkred")
NB1 <- c("bisque1")
Nitron <- c("darkblue")
Nitror <- c("dodgerblue2")
PAU <- c("yellowgreen")
Pates <- c("black")
Planct <- c("lightskyblue")
Alpha <- c("lightpink")
Gamma <- c("lightpink3")
Magnet <- c("lightpink4")
SAR <- c("aquamarine1")
Thermo <- c("hotpink2")
Verru <- c("forestgreen")
WPS2 <- c("orchid2")

#---------Fig 1: Map of sample sites---------

# Read in map shape file

map <- 
  st_read("Honduras_Shape/hnd_admbnda_adm2_sinit_20161005.shp") %>% 
  filter(ADM1_ES == "Islas de La Bahia") %>% 
  st_transform(crs = st_crs(4326)) %>% 
  st_crop(xmin = -86.65, xmax = -86.50, ymin = 16.25, ymax = 16.50)

# Get site data points

sites <-
  read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, Latitude, Longitude) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

# Order sites by geography (W --> E)

sites$Site <- 
  factor(
    sites$Site, 
    levels = c("Phlipper's Peace", "Angus' Aquarium", "Fish Den", "Overheat", 
               "Pillar Coral", "White Hole", "International Diver", "B-S Deep",
               "Wicked Pissa", "Thalassa Arg", "Leo's Den"))

# Plot and save map

map %>% 
  ggplot() +
  geom_sf(data = st_union(map), fill = "grey80") +
  geom_sf(data = sites, size = 7, aes(color = Site)) +
  coord_sf(xlim = c(-86.61, -86.51), ylim = c(16.26, 16.37), expand = FALSE) +
  ggspatial::annotation_scale(location = "br", height = unit(0.2, "cm")) +
  ggspatial::annotation_north_arrow(
    location = "br",
    height = unit(1, "cm"), width = unit(1, "cm"),
    pad_x = unit(0.7, "cm"), pad_y = unit(0.8, "cm")) +
  scale_color_manual(values = site_color) +
  guides(alpha = "none") +
  xlab("") + ylab("") +
  theme_bw()

ggsave("Figure1.tiff", width = 150, height = 110, unit = "mm")

#---------Fig 2: Biofilm taxa bar plot---------

# Group biofilm dataset at Class level

bio_c <- tax_glom(bio, taxrank = "Class")

# Edit taxonomy table of bio_c
# Add Magnetococcia, Alpha-, or Gammaproteobacteria to phylum-level ID
# Change Proteobacteria to Pseudomonadota

tax <- 
  tax_table(bio_c) %>% 
  as.data.frame() %>% 
  select(Kingdom, Phylum, Class) %>% 
  mutate(
    Phylum = 
      str_replace(
        Phylum,
        "Proteobacteria", 
        case_when(
          Class == "Gammaproteobacteria" ~ "Pseudomonadota_Gammaproteobacteria",
          Class == "Magnetococcia" ~ "Pseudomonadota_Magnetococcia",
          Class == "Alphaproteobacteria" ~ "Pseudomonadota_Alphaproteobacteria"))) %>% 
  select(-Class)

# Put edited taxonomy table into the phyloseq object

tax_table(bio_c) <- as.matrix(tax)

# Make read counts relative abundance

bio_c_rel <- transform_sample_counts(bio_c, function(OTU) OTU/sum(OTU))
bio_c_rel <- prune_taxa(taxa_sums(bio_c_rel) > 0, bio_c_rel)

# Transform phyloseq object into dataframe

bio_c_rel_df <- psmelt(bio_c_rel)

# Bin all Phyla with <1% abundance into a category called "Other"

bio_c_rel_df <- bio_c_rel_df %>% dplyr::filter(Abundance > 0)
bio_c_rel_df$Phylum <- as.character(bio_c_rel_df$Phylum)
bio_c_rel_df$Phylum[bio_c_rel_df$Abundance < 0.01] <- "Other"

# Order the bars so "Other" is on the top

Phy <- sort(unique(bio_c_rel_df$Phylum), decreasing = FALSE)
Phy <- Phy[-which(Phy == "Other")]
Phy <- c("Other", Phy)

# Order factor

bio_c_rel_df$Phylum <- factor(bio_c_rel_df$Phylum, levels = Phy)
bio_c_rel_df$Phylum <- as.factor(bio_c_rel_df$Phylum)

# Assign colors

fig2_color <- 
  c(Other, Acido, Actino, Bactero, Bdello, Chloro, 
    Cyano, Dada, Deino, Desulfo, Ento, Firm, Gemma, 
    Myxo, NB1, Nitror, PAU, Planct, Alpha, Gamma, Verru)

# Change site names to abbreviations

bio_c_rel_df <-
  bio_c_rel_df %>% 
  mutate(SiteShort = case_when(
    Site == "Angus' Aquarium" ~ "AA",
    Site == "B-S Deep" ~ "BD",
    Site == "Fish Den" ~ "FD",
    Site == "International Diver" ~ "ID",
    Site == "Leo's Den" ~ "LD",
    Site == "Overheat" ~ "OH",
    Site == "Phlipper's Peace" ~ "PP",
    Site == "Pillar Coral" ~ "PC",
    Site == "Thalassa Arg" ~ "TA",
    Site == "White Hole" ~ "WH",
    Site == "Wicked Pissa" ~ "WP"))

# Order abbreviated site names by geography (W --> E)

bio_c_rel_df$SiteShort <- 
  factor(
    bio_c_rel_df$SiteShort, 
    levels = c("PP", "AA", "FD", "OH", "PC", "WH", 
               "ID", "BD", "WP", "TA", "LD"))

# Plot biofilm taxa bar plot

ggplot(data = bio_c_rel_df, 
       aes(x = fct_reorder2(Sample, desc(SiteShort), desc(SiteShort)), 
           y = Abundance, 
           fill = Phylum)) +
  geom_bar(aes(), stat = "identity") +
  facet_grid(cols = vars(Substrate), scales = "free_x", space = "free_x") +
  scale_fill_manual(values = fig2_color) +
  ylab("Relative Abundance of Phyla (>1% Abundance)") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 14)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

ggsave("Figure2_edit.jpg", height = 8, width = 16)

#---------Fig 3: Alpha diversity and composition between sample types---------

# Need to rarefy for alpha diversity analysis

# Remove samples with <2,500 reads
ps_rare <- subset_samples(ps, sample_sums(ps) > 2500)

# Rarefy to lowest read count
ps_rare <- 
  rarefy_even_depth(ps_rare, sample.size = min(sample_sums(ps_rare)), 
                    rngseed = 8, replace = FALSE)
ps_rare <- prune_taxa(taxa_sums(ps_rare) > 0, ps_rare)

# Calculate diversity indices
ps_rare_div <- adiv_df(ps_rare)

# Plot alpha diversity

alpha <- 
  ggarrange(
    ggplot(ps_rare_div, aes(x = Substrate, y = Shannon, fill = Substrate)) +
      geom_boxplot() +
      scale_fill_manual(values = substrate_color) +
      theme_bw() +
      ylab("Shannon Diversity") +
      ggtitle("a)") + 
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(size = 14),
            axis.title.y = element_text(size = 14)) +
      theme(plot.title.position = "plot"), 
    ggplot(ps_rare_div, aes(x = Substrate, y = Simpson, fill = Substrate)) +
      geom_boxplot() +
      scale_fill_manual(values = substrate_color) +
      theme_bw() +
      ylab("Simpson Diversity (1-D)") +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(size = 14),
            axis.title.y = element_text(size = 14)), 
    ncol = 1, 
    common.legend = TRUE, 
    legend = "none", 
    align = "v")

# Calculate unweighted and weighted Unifrac distances

uu_dist_ps <- UniFrac(ps, weighted = FALSE, normalized = TRUE)
wu_dist_ps <- UniFrac(ps, weighted = TRUE, normalized = TRUE)

# Ordinate weighted distances
ord_wu_ps <- ordinate(ps, "NMDS", distance = wu_dist_ps)
nmdsdata_wu_ps <- plot_ordination(ps, ord_wu_ps)$data

# Ordinate unweighted distances
ord_uu_ps <- ordinate(ps, "NMDS", distance = uu_dist_ps)
nmdsdata_uu_ps <- plot_ordination(ps, ord_uu_ps)$data

# Plot community composition
 
beta <-
  ggarrange(
    ggplot(nmdsdata_wu_ps, aes(x = NMDS1, y = NMDS2, col = Substrate)) +
      geom_point(size = 2, alpha = 0.9) + 
      stat_ellipse() +
      scale_color_manual(values = substrate_color) +
      ggtitle("b)      Weighted Unifrac") +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black")) +
      theme(legend.title = element_text(size = 14),
            legend.text = element_text(size = 14)) +
      theme(plot.title.position = "plot"), 
    ggplot(nmdsdata_uu_ps, aes(x = NMDS1, y = NMDS2, col = Substrate)) +
      geom_point(size = 2, alpha = 0.9) + 
      stat_ellipse() +
      scale_color_manual(values = substrate_color) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black")) +
      theme(legend.title = element_text(size = 14),
            legend.text = element_text(size = 14)) +
      ggtitle("Unweighted Unifrac"), 
    ncol = 1, 
    common.legend = TRUE, 
    legend = "right")

ggarrange(alpha, beta, widths = c(1,2))

ggsave("Figure3.tiff", width = 9, height = 8)

#---------Fig 4: Relative influence of assembly processes by substrate---------

# Read in and clean up data

fig4_data <- 
  read_csv("data/KARMA.iCAMP.BootSummary.Substrate.csv") %>% 
  select(Substrate = Group, Process, Lower.whisker:Upper.whisker) %>% 
  rename(Median = "Median...13") %>% 
  filter(Substrate != "CCA_vs_Carbonate") %>% 
  filter(Process != "Stochasticity") %>% 
  mutate(
    Process = 
      case_when(
        Process == "Heterogeneous.Selection" ~ "HeS",
        Process == "Homogeneous.Selection" ~ "HoS",
        Process == "Dispersal.Limitation" ~ "DL",
        Process == "Homogenizing.Dispersal" ~ "HD",
        Process == "Drift.and.Others" ~ "DR"))

# Order assembly process factor

fig4_data$Process <- 
  factor(
    fig4_data$Process, 
    levels = c("HeS", "HoS", "DL", "HD", "DR"))

# Plot boxplots

ggplot(fig4_data, aes(x = Process)) +
  geom_boxplot(
    stat = "identity",
    aes(ymin = Lower.whisker,
        lower = Lower.hinge,
        middle = Median,
        upper = Upper.hinge,
        ymax = Upper.whisker,
        fill = Substrate)) +
  scale_fill_manual(values = biofilm_color) +
  geom_signif(
    stat = "identity",
    data = data.frame(x = 3.78, xend = 4.22, y = 0.05, annotation = "*"),
    aes(x = x, xend = xend, y = y, yend = y, annotation = annotation)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  ylab("Relative Influence") +
  xlab("Assembly Process")

ggsave("Figure4.tiff", height = 6, width = 6)

#---------Fig 5b: Relative influence of assembly processes by substrate---------

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

# Get relative abundance counts

bio_g_rel <- transform_sample_counts(bio_g, function(OTU) OTU/sum(OTU))
bio_g_rel <- prune_taxa(taxa_sums(bio_g_rel) > 0, bio_g_rel)

# Transform phyloseq object to data frame

bio_g_rel_df <- psmelt(bio_g_rel)

# Add sub-community data into data frame

bio_g_rel_df <- 
  bio_g_rel_df %>% 
  group_by(OTU, Substrate, Sample) %>% 
  transmute(RelAbun = mean(Abundance)) %>% 
  distinct() %>% 
  left_join(
    read_csv("data/abundant-rare-assignments-g.csv"),
    by = join_by("OTU" == "SeqID")) %>% 
  left_join(
    read_csv("data/generalist-specialist-assignments-g.csv") %>% 
      select(OTU = SeqID, Sign),
    by = "OTU") %>% 
  rename(Abun_Rare = Category, Gen_Spec = Sign)

# Get relative abundance counts of each subcommunity for each sample

abundance_sub <-
  bio_g_rel_df %>% 
  group_by(Sample, Abun_Rare) %>% 
  transmute(RelAbun = sum(RelAbun),
            Substrate = Substrate) %>% 
  distinct()

habitat_sub <-
  bio_g_rel_df %>% 
  group_by(Sample, Gen_Spec) %>% 
  transmute(RelAbun = sum(RelAbun),
            Substrate = Substrate) %>% 
  distinct()

# Reorder habitat subcommunity levels

habitat_sub$Gen_Spec <- 
  factor(habitat_sub$Gen_Spec, 
         levels = c("Generalist", "Specialist", "None"))

# Reorder abundance subcommunity levels

abundance_sub$Abun_Rare <- 
  factor(abundance_sub$Abun_Rare, 
         levels = c("Abundant", "Rare", "None"))

# Plot boxplots

ggarrange(
  ggplot(habitat_sub, aes(x = Gen_Spec, y = RelAbun, fill = Substrate)) +
    geom_boxplot() +
    scale_fill_manual(values = biofilm_color) +
    theme_bw() +
    theme(axis.title.x = element_blank()) +
    ylim(0, 1) +
    ylab("Relative Abundance"),
  ggplot(abundance_sub, aes(x = Abun_Rare, y = RelAbun, fill = Substrate)) +
    geom_boxplot() +
    scale_fill_manual(values = biofilm_color) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y= element_blank()) +
    ylim(0, 1),
  common.legend = TRUE,
  widths = c(1.1,1))

ggsave("Figure5b.jpg", width = 6, height = 5)

#---------Fig 6: Assembly process by subcommunity + substrate---------

# Get all data in the right format

fig6_data <-
  rbind(
    read_csv("data/KARMA.ProcessImportance_EachGroup.csv") %>% 
      filter(Group == "CCA" | Group == "Carbonate") %>% 
      select(Substrate = Group, HeS:DR) %>% 
      mutate(Community = rep("Whole", 2)) %>% 
      pivot_longer(HeS:DR, names_to = "Process", values_to = "RelInf"),
    
    read_csv("data/KARMA.iCAMP.GenSpec.Process_EachGroup_EachCategory.csv") %>% 
      filter(Group == "CCA" | Group == "Carbonate") %>% 
      select(Substrate = Group, None.HeS:Generalist.Stochasticity) %>% 
      pivot_longer(None.HeS:Generalist.Stochasticity, names_to = "Process", values_to = "RelInf") %>% 
      separate_wider_delim(Process, ".", names = c("Community", "Process")) %>% 
      filter(Community != "None") %>% 
      filter(Process != "Stochasticity"),
    
    read_csv("data/KARMA.iCAMP.AbunRare.Process_EachGroup_EachCategory.csv") %>% 
      filter(Group == "CCA" | Group == "Carbonate") %>% 
      select(Substrate = Group, Abundant.HeS:Rare.Stochasticity) %>% 
      pivot_longer(Abundant.HeS:Rare.Stochasticity, names_to = "Process", values_to = "RelInf") %>% 
      separate_wider_delim(Process, ".", names = c("Community", "Process")) %>% 
      filter(Community != "None") %>% 
      filter(Process != "Stochasticity"))

# Order factors

fig6_data$Process <- 
  factor(fig6_data$Process, 
         levels = c("HeS", "HoS", "DL", "HD", "DR"))

fig6_data$Community <- 
  factor(fig6_data$Community, 
         levels = c("Whole", "Generalist", "Specialist", "Abundant", "Rare"))

ggplot(fig6_data, aes(x = Community, y = RelInf, fill = Process)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = process_color) +
  facet_wrap(~Substrate) + 
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  ylab("Relative Influence")

ggsave("Figure6.tiff", height = 5, width = 7.5)

#---------Fig 7: Correlation heatmap: environmental parameters vs phylogenetic bins by process---------

# Get necessary data

fig7_data <-
  rbind(
    read_tsv("data/Maaslin2_CCA_significant_results.tsv") %>% 
      filter(N.not.0 >= 20) %>% 
      select(Bin = feature, Env = value, Coefficient = coef) %>% 
      mutate(Substrate = rep("CCA", nrow(.))) %>% 
      left_join(
        read_csv("data/KARMA.ProcessImportance_EachBin_EachGroup.csv") %>% 
          filter(Group == "CCA") %>% 
          select(Index, contains("bin")) %>% 
          column_to_rownames(var = "Index") %>% 
          t() %>% 
          data.frame() %>% 
          select(DominantProcess) %>% 
          rownames_to_column(var = "Bin") %>% 
          mutate(Bin = str_replace(Bin, "b", "B")),
        by = "Bin") %>% 
      rename(Process = DominantProcess) %>% 
      mutate(Bin = as.character(str_replace(Bin, "Bin", ""))),
    read_tsv("data/Maaslin2_carbonate_significant_results.tsv") %>% 
      filter(N.not.0 >= 20) %>% 
      select(Bin = feature, Env = value, Coefficient = coef) %>% 
      mutate(Substrate = rep("Carbonate", nrow(.))) %>% 
      left_join(
        read_csv("data/KARMA.ProcessImportance_EachBin_EachGroup.csv") %>% 
          filter(Group == "Carbonate") %>% 
          select(Index, contains("bin")) %>% 
          column_to_rownames(var = "Index") %>% 
          t() %>% 
          data.frame() %>% 
          select(DominantProcess) %>% 
          rownames_to_column(var = "Bin") %>% 
          mutate(Bin = str_replace(Bin, "b", "B")),
        by = "Bin") %>% 
      rename(Process = DominantProcess) %>% 
      mutate(Bin = str_remove(Bin, "Bin"))) %>% 
  mutate(Env = 
           case_when(
             Env == "Ammonium" ~ "Ammonium",
             Env == "AvgSampTemp" ~ "Average Sample Temp",
             Env == "AvgSiteTemp" ~ "Average Site Temp",
             Env == "CoralCover" ~ "Coral Cover",
             Env == "CoralSimp" ~ "Coral Simpson",
             Env == "Depth" ~ "Depth",
             Env == "MacroCover" ~ "Macroalgae Cover",
             Env == "MacroShan" ~ "Macroalgae Shannon",
             Env == "NitrateNitrite" ~ "Nitrate + Nitrite",
             Env == "pH" ~ "pH",
             Env == "Salinity" ~ "Salinity",
             Env == "Silicate" ~ "Silicate"))

bin_taxa <- 
  read_csv("data/KARMA.Taxon_Bin.csv") %>% 
  select(ID, Bin)

bin_key <-
  read_csv("data/KARMA.Bin_TopTaxon.csv") %>% 
  select(Bin,
         Phylum = TopTaxon.Phylum,
         Genus = TopTaxon.Genus) %>% 
  mutate(Bin = str_remove(Bin, "Bin") %>% as.numeric())

# Join BinIDs with taxonomy table

tax_bin <- 
  tax_table(bio_g) %>% 
  data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  left_join(bin_taxa, by = "ID") %>% 
  column_to_rownames(var = "ID")

bio_g_bin <- bio_g

tax_table(bio_g_bin) <- tax_table(tax_bin %>% as.matrix())

# Get relative abundance counts

bio_g_rel_bin <- transform_sample_counts(bio_g_bin, function(OTU) OTU/sum(OTU))
bio_g_rel_bin <- prune_taxa(taxa_sums(bio_g_rel_bin) > 0, bio_g_rel_bin)

# Transform phyloseq object to data frame

bio_g_rel_bin_df <-
  psmelt(bio_g_rel_bin) %>% 
  select(Sample, Substrate, Abundance, Bin) %>% 
  mutate(Bin = str_remove(Bin, "Bin") %>% as.numeric()) %>% 
  group_by(Sample, Substrate, Bin) %>% 
  transmute(Abundance = sum(Abundance)) %>% 
  ungroup() %>% 
  group_by(Substrate, Bin) %>% 
  transmute(Abundance = mean(Abundance)) %>% 
  distinct() %>% 
  left_join(bin_key %>% 
              mutate(Bin = str_remove(Bin, "Bin") %>% as.numeric), 
            by = "Bin")

fig7_data <-
  fig7_data %>% 
  left_join(bio_g_rel_bin_df %>% 
              mutate(Bin = as.character(Bin)), 
            by = c("Bin", "Substrate")) %>% 
  mutate(Abundance = round((Abundance*100), 2),
         Perc = rep("%", nrow(.))) %>% 
  unite("GenusRel", Genus, Abundance, sep = ", ", remove = FALSE) %>% 
  unite("GenusRel", GenusRel, Perc, sep = "")

# Order environmental factors

fig7_data$Env <-
  factor(fig7_data$Env,
         levels = c("Ammonium", "Nitrate + Nitrite", "Silicate",
                    "Depth", "Average Sample Temp", "Average Site Temp",
                    "pH", "Salinity", "Coral Cover", "Coral Simpson",
                    "Macroalgae Cover", "Macroalgae Shannon"))

# Order processes

fig7_data$Process <-
  factor(fig7_data$Process,
         levels = c("HeS", "HoS", "DL", "DR"))

# Order Genus levels by the Phylum they belong to for plot aesthetics

fig7_data$Genus <- 
  factor(
    fig7_data$Genus, 
    levels = 
      fig7_data %>% 
      group_by(Phylum) %>% 
      arrange(Genus, .by_group = TRUE) %>% 
      pull(Genus) %>% 
      unique(), 
    ordered = TRUE)

# Plot heatmap

ggplot(fig7_data, aes(x = Env, y = Genus, fill = Coefficient)) + 
  geom_tile() +
  facet_grid(rows = vars(Process), cols = vars(Substrate), 
             scales = "free", space = "free", switch = "y") + 
  scale_fill_gradientn(colors = c("darkslategray4", "white", "tomato2")) +
  theme_bw() +
  scale_y_discrete(position = "right", limits = rev) +
  theme(strip.background = element_rect("white")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_blank())

ggsave("Figure7.tiff", height = 9.5, width = 7)

#---------Fig S1: Rarefaction curves---------

# Get ASV table

ASV_table_ps <- otu_table(ps) %>% t() %>% data.frame()

# Calculate rarefaction curves

richness_rare <- rarecurve(ASV_table_ps, step = 20, label = FALSE, tidy = TRUE)

# Add sample data

richness_rare <- 
  richness_rare %>% 
  rename(SampleID = Site, Depth = Sample, OTUs = Species) %>% 
  left_join(
    sample_data(ps) %>% 
      data.frame() %>% 
      rownames_to_column(var = "SampleID") %>% 
      select(SampleID, Site, SampleType, Substrate))

# Plot
richness_curve <-
  ggplot(richness_rare, aes(x = Depth, y = OTUs, fill = SampleID)) +
  geom_line(aes(col = Substrate)) +
  scale_color_manual(values = substrate_color) +
  xlab("Read depth") +
  ylab("Number of unique OTUs observed") +
  theme_bw()

# Calculate Shannon diversity rarefaction curve

Shannon_rare <- data.frame()

for (i in seq(from = 100, to = 58579, by = 100)) {
  
  x <- rarefy_even_depth(ps, sample.size = i, rngseed = 88, 
                         replace = FALSE, trimOTUs = TRUE, verbose = FALSE)
  x_div <- estimate_richness(x, measures=c("Shannon"))
  
  Shannon_rare <<-
    rbind(Shannon_rare,
          x_div %>% 
            rownames_to_column(var = "SampleID") %>% 
            mutate(Depth = rep(i, nrow(.))))
  
}

# Clean data

Shannon_rare <- 
  Shannon_rare %>% 
  mutate(SampleID = str_replace(SampleID, "^1$", "K70"))

# Add sample data

Shannon_rare <-
  Shannon_rare %>% 
  left_join(sample_data(ps) %>% 
              data.frame() %>% 
              rownames_to_column(var = "SampleID") %>% 
              select(SampleID, Substrate), 
            by = "SampleID")

# Plot Shannon diversity curve

shannon_curve <-
  ggplot(Shannon_rare, aes(x = Depth, y = Shannon, fill = SampleID)) +
  geom_line(aes(col = Substrate)) +
  xlab("Read depth") +
  ylab("Shannon Diversity") +
  scale_color_manual(values = substrate_color) +
  theme_bw()

# Calculate Simpson diversity rarefaction curve

Simpson_rare <- data.frame()

for (i in seq(from = 100, to = 58579, by = 100)) {
  
  x <- rarefy_even_depth(ps, sample.size = i, rngseed = 88, 
                         replace = FALSE, trimOTUs = TRUE, verbose = FALSE)
  x_div <- estimate_richness(x, measures = c("Simpson"))
  
  Simpson_rare <<-
    rbind(Simpson_rare,
          x_div %>% 
            rownames_to_column(var = "SampleID") %>% 
            mutate(Depth = rep(i, nrow(.))))
}

# Clean data

Simpson_rare <- 
  Simpson_rare %>% 
  mutate(SampleID = str_replace(SampleID, "^1$", "K70"))

# Add sample data

Simpson_rare <-
  Simpson_rare %>% 
  left_join(sample_data(ps) %>% 
              data.frame() %>% 
              rownames_to_column(var = "SampleID") %>% 
              select(SampleID, Substrate), 
            by = "SampleID")

# Plot
simpson_curve <-
  ggplot(Simpson_rare, aes(x = Depth, y = Simpson, fill = SampleID)) +
  geom_line(aes(col = Substrate)) +
  xlab("Read depth") +
  ylab("Simpson Diversity (1-D)") +
  scale_color_manual(values = substrate_color) +
  theme_bw()

# Plot combined figure

ggarrange(richness_curve, shannon_curve, simpson_curve, 
          nrow = 1, common.legend = TRUE, legend = "right")

ggsave("FigureS1.tiff", width = 14, height = 5)

#---------Fig S2: Environmental characteristics by site---------

# Get site-level environmental data

figS2_data <- 
  read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, AvgSiteTemp:Silicate) %>%
  filter(!is.na(Site)) %>% 
  distinct()

# Make short site names

figS2_data <-
  figS2_data %>% 
  mutate(SiteShort = case_when(
    Site == "Angus' Aquarium" ~ "AA",
    Site == "B-S Deep" ~ "BD",
    Site == "Fish Den" ~ "FD",
    Site == "International Diver" ~ "ID",
    Site == "Leo's Den" ~ "LD",
    Site == "Overheat" ~ "OH",
    Site == "Phlipper's Peace" ~ "PP",
    Site == "Pillar Coral" ~ "PC",
    Site == "Thalassa Arg" ~ "TA",
    Site == "White Hole" ~ "WH",
    Site == "Wicked Pissa" ~ "WP"))

# Order short site factor

figS2_data$SiteShort <- 
  factor(figS2_data$SiteShort, 
         levels = c("PP", "AA", "FD", "OH", "PC", "WH", 
                    "ID", "BD", "WP", "TA", "LD"))

figS2_data$Site <- 
  factor(
    figS2_data$Site, 
    levels = c("Phlipper's Peace", "Angus' Aquarium", "Fish Den", "Overheat", 
               "Pillar Coral", "White Hole", "International Diver",
               "B-S Deep", "Wicked Pissa", "Thalassa Arg", "Leo's Den"))

# Make plots

figS2_temp_plot <-
  ggplot(figS2_data, aes(x = SiteShort, y = AvgSiteTemp, color = Site)) +
  geom_linerange(aes(ymin = AvgSiteTemp-SDSiteTemp, ymax = AvgSiteTemp+SDSiteTemp)) +
  geom_point(size = 5) +
  scale_color_manual(values = site_color) +
  theme_bw() +
  ylab("Average site temperature (ºC)") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

figS2_sal_plot <-
  ggplot(figS2_data, aes(x = SiteShort, y = Salinity, color = Site)) +
  geom_point(size = 5) +
  scale_color_manual(values = site_color) +
  theme_bw() +
  ylab("Salinity (ppt)") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

figS2_pH_plot <-
  ggplot(figS2_data, aes(x = SiteShort, y = pH, color = Site)) +
  geom_point(size = 5) +
  scale_color_manual(values = site_color) +
  theme_bw() +
  ylab("pH") +
  theme(axis.title.x = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

figS2_nitnit_plot <-
  ggplot(figS2_data, aes(x = SiteShort, y = NitrateNitrite, color = Site)) +
  geom_point(size = 5) +
  scale_color_manual(values = site_color) +
  theme_bw() +
  ylab("Nitrate + nitrite (µM)") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

figS2_amm_plot <-
  ggplot(figS2_data, aes(x = SiteShort, y = Ammonium, color = Site)) +
  geom_point(size = 5) +
  scale_color_manual(values = site_color) +
  theme_bw() +
  ylab("Ammonium (µM)") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

figS2_sili_plot <-
  ggplot(figS2_data, aes(x = SiteShort, y = Silicate, color = Site)) +
  geom_point(size = 5) +
  scale_color_manual(values = site_color) +
  theme_bw() +
  ylab("Silicate (µM)") +
  theme(axis.title.x = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

ggarrange(figS2_temp_plot, figS2_nitnit_plot, figS2_sal_plot, 
          figS2_amm_plot, figS2_pH_plot, figS2_sili_plot,
          nrow = 3, ncol = 2, common.legend = TRUE, legend = "right", align = "hv")

ggsave("FigureS2.tiff", height = 7, width = 8)

#---------Fig S3: Coral composition, cover, and diversity---------

# Get coral composition data and standardize

coral_data <- 
  read_csv("data/KARMA_environmental_data.csv") %>%
  select(Site = Name, ACER:SSID) %>% 
  column_to_rownames(var = "Site") %>% 
  mutate(Coral = rowSums(.)) %>% 
  mutate(across(ACER:SSID, .fns = ~./Coral)) %>% 
  select(-Coral) %>% 
  rownames_to_column(var = "Site") %>% 
  pivot_longer(cols = ACER:SSID, names_to = "Species", values_to = "Abundance")

# Bin all coral species with <3% abundance into a category called "Other"

coral_data <- coral_data %>% dplyr::filter(Abundance > 0)
coral_data$Species <- as.character(coral_data$Species)
coral_data$Species[coral_data$Abundance < 0.03] <- "Other"

# Order the bars so Other is on the top

coral_vec <- sort(unique(coral_data$Species), decreasing = FALSE)
coral_vec <- coral_vec[-which(coral_vec == "Other")]
coral_vec <- c("Other", coral_vec)

# Order factor

coral_data$Species <- factor(coral_data$Species, levels = coral_vec)
coral_data$Species <- as.factor(coral_data$Species)

# Make short site names

coral_data <-
  coral_data %>% 
  mutate(SiteShort = case_when(
    Site == "Angus' Aquarium" ~ "AA",
    Site == "B-S Deep" ~ "BD",
    Site == "Fish Den" ~ "FD",
    Site == "International Diver" ~ "ID",
    Site == "Leo's Den" ~ "LD",
    Site == "Overheat" ~ "OH",
    Site == "Phlipper's Peace" ~ "PP",
    Site == "Pillar Coral" ~ "PC",
    Site == "Thalassa Arg" ~ "TA",
    Site == "White Hole" ~ "WH",
    Site == "Wicked Pissa" ~ "WP"))

# Order short site factor

coral_data$SiteShort <- 
  factor(coral_data$SiteShort, 
         levels = c("PP", "AA", "FD", "OH", "PC", "WH", 
                    "ID", "BD", "WP", "TA", "LD"))

# Plot coral composition
coral_composition <-
  ggplot(data = coral_data, aes(x = SiteShort, y = Abundance, fill = Species)) +
  geom_bar(aes(), stat = "identity") +
  scale_fill_manual(
    "Coral Species",
    values = coral_color,
    labels = c("Other", 
               expression(italic("Agaricia agaricites")),
               expression(italic("Agaricia tenuifolia")), 
               expression(italic("Diploria labyrinthiformis")),
               expression(italic("Millepora alcicornis")),
               expression(italic("Montastraea cavernosa")),
               expression(italic("Madracis decactis")),
               expression(italic("Orbicella annularis")), 
               expression(italic("Orbicella faveolata")), 
               expression(italic("Orbicella franksii")),
               expression(italic("Porites astreoides")), 
               expression(italic("Porites divaricata")),
               expression(italic("Porites furcata")),
               expression(italic("Porites porites")), 
               expression(italic("Stephanocoenia intersepta")), 
               expression(italic("Siderastrea siderea")))) +
  ylab("Relative Abundance of Species (>2% Abundance)") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank()) +
  theme(legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))

# Get coral diversity data
coral_div_data <- 
  read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, CoralCover, CoralRich, CoralShan, CoralSimp) %>% 
  filter(!is.na(Site)) %>% 
  distinct()

# Make short site names

coral_div_data <-
  coral_div_data %>% 
  mutate(SiteShort = case_when(
    Site == "Angus' Aquarium" ~ "AA",
    Site == "B-S Deep" ~ "BD",
    Site == "Fish Den" ~ "FD",
    Site == "International Diver" ~ "ID",
    Site == "Leo's Den" ~ "LD",
    Site == "Overheat" ~ "OH",
    Site == "Phlipper's Peace" ~ "PP",
    Site == "Pillar Coral" ~ "PC",
    Site == "Thalassa Arg" ~ "TA",
    Site == "White Hole" ~ "WH",
    Site == "Wicked Pissa" ~ "WP"))

# Order short site factor

coral_div_data$SiteShort <- 
  factor(coral_div_data$SiteShort, 
         levels = c("PP", "AA", "FD", "OH", "PC", "WH", 
                    "ID", "BD", "WP", "TA", "LD"))

# Plot coral cover

coral_cover <-
  ggplot(coral_div_data, aes(x = SiteShort, y = CoralCover, color = SiteShort)) +
  geom_point(size = 5) +
  scale_color_manual("Site", values = site_color) +
  theme_bw() +
  ylab("Cover (%)") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

# Plot coral richness

coral_richness <-
  ggplot(coral_div_data, aes(x = SiteShort, y = CoralRich, color = SiteShort)) +
  geom_point(size = 5) +
  scale_color_manual("Site", values = site_color) +
  theme_bw() +
  ylab("Species Richness") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  theme(panel.grid.minor = element_blank())

# Plot coral Shannon diversity

coral_shannon <-
  ggplot(coral_div_data, aes(x = SiteShort, y = CoralShan, color = SiteShort)) +
  geom_point(size = 5) +
  scale_color_manual("Site", values = site_color) +
  theme_bw() +
  ylab("Shannon") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

# Plot coral Simpson diversity

coral_simpson <-
  ggplot(coral_div_data, aes(x = SiteShort, y = CoralSimp, color = SiteShort)) +
  geom_point(size = 5) +
  scale_color_manual("Site", values = site_color) +
  theme_bw() +
  ylab("Simpson") +
  theme(axis.title.x = element_blank()) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) +
  theme(panel.grid.minor = element_blank())

# Put the four diversity plots together

coral_div <- 
  ggarrange(coral_cover, coral_richness, coral_shannon, coral_simpson, 
            ncol = 1, common.legend = TRUE, align = "v", legend = "right")

# Build full plot

ggarrange(coral_composition, coral_div, ncol = 1, heights = c(1.75,2))

ggsave("FigureS3.tiff", width = 7, height = 9.5)

#---------Fig S4: Macroalgae composition, cover, and diversity---------

# Get algae composition data and standardize

algae_data <- 
  read_csv("data/KARMA_environmental_data.csv") %>%
  select(Site = Name, Amphiroa:TurfAlgae) %>% 
  column_to_rownames(var = "Site") %>% 
  mutate(Algae = rowSums(.)) %>% 
  mutate(across(Amphiroa:TurfAlgae, .fns = ~./Algae)) %>% 
  select(-Algae) %>% 
  rownames_to_column(var = "Site") %>% 
  pivot_longer(cols = Amphiroa:TurfAlgae, names_to = "Species", values_to = "Abundance")

# Make short site names

algae_data <-
  algae_data %>% 
  mutate(SiteShort = case_when(
    Site == "Angus' Aquarium" ~ "AA",
    Site == "B-S Deep" ~ "BD",
    Site == "Fish Den" ~ "FD",
    Site == "International Diver" ~ "ID",
    Site == "Leo's Den" ~ "LD",
    Site == "Overheat" ~ "OH",
    Site == "Phlipper's Peace" ~ "PP",
    Site == "Pillar Coral" ~ "PC",
    Site == "Thalassa Arg" ~ "TA",
    Site == "White Hole" ~ "WH",
    Site == "Wicked Pissa" ~ "WP"))

# Order short site factor

algae_data$SiteShort <- 
  factor(algae_data$SiteShort, 
         levels = c("PP", "AA", "FD", "OH", "PC", "WH", 
                    "ID", "BD", "WP", "TA", "LD"))

# Order the taxa so turf is on the top

algae_data$Species <- 
  factor(algae_data$Species, 
         levels = c("TurfAlgae", "Amphiroa", "Dictyota", "Halimeda", "Lobophora", "Cyanobacteria"))
algae_data$Species <- as.factor(algae_data$Species)

# Plot algae composition

algae_composition <-
  ggplot(data = algae_data, aes(x = SiteShort, y = Abundance, fill = Species)) +
  geom_bar(aes(), stat = "identity") +
  scale_fill_manual(
    "Algae Taxa",
    values = algae_color,
    labels = c("Turf algae", 
               expression(italic("Amphiroa")),
               expression(italic("Dictyota")), 
               expression(italic("Halimeda")),
               expression(italic("Lobophora")), 
               "Cyanobacteria")) +
  ylab("Relative Abundance of Taxa") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank()) +
  theme(legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))

# Get algae diversity data

algae_div_data <- 
  read_csv("data/KARMA_environmental_data.csv") %>% 
  select(Site = Name, MacroCover, MacroShan, MacroSimp) %>% 
  filter(!is.na(Site)) %>% 
  distinct()

# Make short site names

algae_div_data <-
  algae_div_data %>% 
  mutate(SiteShort = case_when(
    Site == "Angus' Aquarium" ~ "AA",
    Site == "B-S Deep" ~ "BD",
    Site == "Fish Den" ~ "FD",
    Site == "International Diver" ~ "ID",
    Site == "Leo's Den" ~ "LD",
    Site == "Overheat" ~ "OH",
    Site == "Phlipper's Peace" ~ "PP",
    Site == "Pillar Coral" ~ "PC",
    Site == "Thalassa Arg" ~ "TA",
    Site == "White Hole" ~ "WH",
    Site == "Wicked Pissa" ~ "WP"))

# Order short site factor

algae_div_data$SiteShort <- 
  factor(algae_div_data$SiteShort, 
         levels = c("PP", "AA", "FD", "OH", "PC", "WH", 
                    "ID", "BD", "WP", "TA", "LD"))

# Plot algae cover

algae_cover <-
  ggplot(algae_div_data, aes(x = SiteShort, y = MacroCover, color = SiteShort)) +
  geom_point(size = 5) +
  scale_color_manual("Site", values = site_color) +
  theme_bw() +
  ylab("Cover (%)") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

# Plot algae shannon diversity

algae_shannon <-
  ggplot(algae_div_data, aes(x = SiteShort, y = MacroShan, color = SiteShort)) +
  geom_point(size = 5) +
  scale_color_manual("Site", values = site_color) +
  theme_bw() +
  ylab("Shannon") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

# Plot algae simpson diversity

algae_simpson <-
  ggplot(algae_div_data, aes(x = SiteShort, y = MacroSimp, color = SiteShort)) +
  geom_point(size = 5) +
  scale_color_manual("Site", values = site_color) +
  theme_bw() +
  ylab("Simpson") +
  theme(axis.title.x = element_blank()) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

# Put the three diversity plots together

algae_div <- 
  ggarrange(algae_cover, algae_shannon, algae_simpson, 
            ncol = 1, common.legend = TRUE, align = "v", legend = "right")

# Build full plot

ggarrange(algae_composition, algae_div, ncol = 1, heights = c(1.75,2))

ggsave("FigureS4.tiff", width = 6, height = 8)

#---------Fig S5: Bacterioplankton taxa bar plot---------

# Group bacterioplankton at Class level

water_c <- tax_glom(water, taxrank = "Class")

# Edit taxonomy table of water_c
# Add Magnetococcia, Alpha-, or Gammaproteobacteria to phylum-level ID
# Change Proteobacteria to Pseudomonadota

tax <- 
  tax_table(water_c) %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  select(Kingdom, Phylum, Class) %>% 
  mutate(
    Phylum = 
      str_replace(
        Phylum,
        "Proteobacteria", 
        case_when(
          Class == "Gammaproteobacteria" ~ "Pseudomonadota_Gammaproteobacteria",
          Class == "Magnetococcia" ~ "Pseudomonadota_Magnetococcia",
          Class == "Alphaproteobacteria" ~ "Pseudomonadota_Alphaproteobacteria"))) %>% 
  select(-Class)

# Put edited taxonomy table into the phyloseq object

tax_table(water_c) <- as.matrix(tax)

# Make read counts relative abundance

water_c_rel <- transform_sample_counts(water_c, function(OTU) OTU/sum(OTU))
water_c_rel <- prune_taxa(taxa_sums(water_c_rel) > 0, water_c_rel)

# Transform phyloseq object into dataframe

water_c_rel_df <- psmelt(water_c_rel)

# Bin all Phyla with <1% abundance into a category called "Other"

water_c_rel_df <- water_c_rel_df %>% dplyr::filter(Abundance > 0)
water_c_rel_df$Phylum <- as.character(water_c_rel_df$Phylum)
water_c_rel_df$Phylum[water_c_rel_df$Abundance < 0.01] <- "Other"

# Order the bars so "Other" is on the top

Phy <- sort(unique(water_c_rel_df$Phylum), decreasing = FALSE)
Phy <- Phy[-which(Phy == "Other")]
Phy <- c("Other", Phy)

# Order factor

water_c_rel_df$Phylum <- factor(water_c_rel_df$Phylum, levels = Phy)
water_c_rel_df$Phylum <- as.factor(water_c_rel_df$Phylum)

# Assign colors

color <- c(Other, Actino, Bactero, Bdello, Cyano, Deino, Marin, Alpha, Gamma, SAR, Verru)

# Change site names to abbreviations

water_c_rel_df <-
  water_c_rel_df %>% 
  mutate(SiteShort = case_when(
    Sample == "K49" ~ "AA",
    Sample == "K21" ~ "BD",
    Sample == "K07" ~ "FD",
    Sample == "K63" ~ "ID",
    Sample == "K14" ~ "LD",
    Sample == "K77" ~ "OH",
    Sample == "K42" ~ "PP",
    Sample == "K70" ~ "PC",
    Sample == "K56" ~ "TA",
    Sample == "K28" ~ "WH",
    Sample == "K35" ~ "WP"))

# Order abbreviated site names by geography (W --> E)

water_c_rel_df$SiteShort <- 
  factor(
    water_c_rel_df$SiteShort, 
    levels = c("PP", "AA", "FD", "OH", "PC", "WH", 
               "ID", "BD", "WP", "TA", "LD"))

# Plot bacterioplankton taxa bar plot

ggplot(data = water_c_rel_df, aes(x = SiteShort, y = Abundance, fill = Phylum)) +
  geom_bar(aes(), stat = "identity") +
  scale_fill_manual(values = color) +
  ylab("Relative Abundance of Phyla (>1% Abundance)") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 12),
        axis.title.x = element_blank()) +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "right",
        legend.text = element_text(size = 12))

ggsave("FigureS5.tiff", height = 6, width = 10)

#---------Fig S6: Bacterioplankton alpha + beta diversity by site---------

# Remove samples with <2,500 reads

water_rare <- subset_samples(water, sample_sums(water) > 2500)

# Rarefy to lowest read depth

water_rare <- 
  rarefy_even_depth(
    water_rare, 
    sample.size = min(sample_sums(water_rare)), 
    rngseed = 8, 
    replace = FALSE)

# Calculate diversity indices

water_rare_div <- adiv_df(water_rare)

# Add short site names

water_rare_div <-
  water_rare_div %>% 
  mutate(SiteShort = case_when(
    Site == "Angus' Aquarium" ~ "AA",
    Site == "B-S Deep" ~ "BD",
    Site == "Fish Den" ~ "FD",
    Site == "International Diver" ~ "ID",
    Site == "Leo's Den" ~ "LD",
    Site == "Overheat" ~ "OH",
    Site == "Phlipper's Peace" ~ "PP",
    Site == "Pillar Coral" ~ "PC",
    Site == "Thalassa Arg" ~ "TA",
    Site == "White Hole" ~ "WH",
    Site == "Wicked Pissa" ~ "WP"))

# Order the site names by geography (W --> E)

water_rare_div$SiteShort <- 
  factor(
    water_rare_div$SiteShort, 
    levels = c("PP", "AA", "FD", "OH", "PC", "WH", 
               "ID", "BD", "WP", "TA", "LD"))

water_rare_div$Site <- 
  factor(
    water_rare_div$Site, 
    levels = c("Phlipper's Peace", "Angus' Aquarium", "Fish Den", "Overheat", 
               "Pillar Coral", "White Hole", "International Diver",
               "B-S Deep", "Wicked Pissa", "Thalassa Arg", "Leo's Den"))

# Plot alpha diversity

shannon_water <-
  ggplot(water_rare_div, aes(x = SiteShort, y = Shannon, color = Site)) +
  geom_point(size = 5) +
  scale_color_manual(values = site_color) +
  theme_bw() +
  ylab("Shannon Diversity") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

simpson_water <-
  ggplot(water_rare_div, aes(x = SiteShort, y = Simpson, color = Site)) +
  geom_point(size = 5) +
  scale_color_manual(values = site_color) +
  theme_bw() +
  ylab("Simpson Diversity (1-D)") +
  theme(axis.title.x = element_blank()) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

figS6a <- ggarrange(shannon_water, simpson_water, ncol = 1, 
                   common.legend = TRUE, legend = "none", align = "v")

# Calculate unweighted and weighted Unifrac distances

uu_dist_water <- UniFrac(water, weighted = FALSE, normalized = TRUE)
wu_dist_water <- UniFrac(water, weighted = TRUE, normalized = TRUE)

# Ordinate

ord_wu_water <- ordinate(water, "NMDS", distance = wu_dist_water)
nmdsdata_wu_water <- plot_ordination(water, ord_wu_water)$data

# Order sites geographically (W --> E)

nmdsdata_wu_water$Site <- 
  factor(
    nmdsdata_wu_water$Site, 
    levels = c("Phlipper's Peace", "Angus' Aquarium", "Fish Den", "Overheat", 
               "Pillar Coral", "White Hole", "International Diver",
               "B-S Deep", "Wicked Pissa", "Thalassa Arg", "Leo's Den"))

# Ordinate

ord_uu_water <- ordinate(water, "NMDS", distance = uu_dist_water)
nmdsdata_uu_water <- plot_ordination(water, ord_uu_water)$data

# Order sites geographically (W --> E)

nmdsdata_uu_water$Site <- 
  factor(
    nmdsdata_uu_water$Site, 
    levels = c("Phlipper's Peace", "Angus' Aquarium", "Fish Den", "Overheat", 
               "Pillar Coral", "White Hole", "International Diver",
               "B-S Deep", "Wicked Pissa", "Thalassa Arg", "Leo's Den"))

# Plot weighted UniFrac NMDS

wu_water_nmds <- 
  ggplot(nmdsdata_wu_water, aes(x = NMDS1, y = NMDS2, color = Site)) +
  geom_point(size = 5, alpha = 0.9) + 
  scale_color_manual(values = site_color) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  ggtitle("Weighted Unifrac")

# Plot weighted UniFrac NMDS

uu_water_nmds <- 
  ggplot(nmdsdata_uu_water, aes(x = NMDS1, y = NMDS2, color = Site)) +
  geom_point(size = 5, alpha = 0.9) + 
  scale_color_manual(values = site_color) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  ggtitle("Unweighted Unifrac")

# Put NMDS plots together

figS6b <- 
  ggarrange(wu_water_nmds, uu_water_nmds, nrow = 1, 
            common.legend = TRUE, legend = "right")

# Put alpha and beta diversity plots together

ggarrange(figS6a, figS6b, nrow = 1, common.legend = TRUE, 
          legend = "bottom", widths = c(1,2.75))

ggsave("FigureS6.tiff", width = 12, height = 5)

#---------Fig S7: CCA alpha + beta diversity by site---------

# Remove samples with <2,500 reads

CCA_rare <- subset_samples(CCA, sample_sums(CCA) > 2500)

# Rarefy to lowest read depth

CCA_rare <- 
  rarefy_even_depth(
    CCA_rare, 
    sample.size = min(sample_sums(CCA_rare)), 
    rngseed = 8, 
    replace = FALSE)

# Calculate diversity indices

CCA_rare_div <- adiv_df(CCA_rare)

# Make short site names

CCA_rare_div <-
  CCA_rare_div %>% 
  mutate(SiteShort = case_when(
    Site == "Angus' Aquarium" ~ "AA",
    Site == "B-S Deep" ~ "BD",
    Site == "Fish Den" ~ "FD",
    Site == "International Diver" ~ "ID",
    Site == "Leo's Den" ~ "LD",
    Site == "Overheat" ~ "OH",
    Site == "Phlipper's Peace" ~ "PP",
    Site == "Pillar Coral" ~ "PC",
    Site == "Thalassa Arg" ~ "TA",
    Site == "White Hole" ~ "WH",
    Site == "Wicked Pissa" ~ "WP"))

# Order short site factor

CCA_rare_div$SiteShort <- 
  factor(CCA_rare_div$SiteShort, 
         levels = c("PP", "AA", "FD", "OH", "PC", "WH", 
                    "ID", "BD", "WP", "TA", "LD"))

# Order site factor

CCA_rare_div$Site <- 
  factor(
    CCA_rare_div$Site, 
    levels = c("Phlipper's Peace", "Angus' Aquarium", "Fish Den", "Overheat", 
               "Pillar Coral", "White Hole", "International Diver",
               "B-S Deep", "Wicked Pissa", "Thalassa Arg", "Leo's Den"))

# Plot CCA alpha diversity by site

shannon_CCA <-
  ggplot(CCA_rare_div, aes(x = SiteShort, y = Shannon, fill = Site)) +
  geom_boxplot() +
  scale_fill_manual(values = site_color) +
  theme_bw() +
  ylab("Shannon Diversity") +
  theme(axis.title.x = element_blank()) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

simpson_CCA <-
  ggplot(CCA_rare_div, aes(x = SiteShort, y = Simpson, fill = Site)) +
  geom_boxplot() +
  scale_fill_manual(values = site_color) +
  theme_bw() +
  ylab("Simpson Diversity (1-D)") +
  theme(axis.title.x = element_blank()) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

figS7a <-
  ggarrange(shannon_CCA, simpson_CCA, ncol = 1, common.legend = TRUE, 
            legend = "none", align = "v")

# Calculate unweighted and weighted Unifrac distances

uu_dist_CCA <- UniFrac(CCA, weighted = FALSE, normalized = TRUE)
wu_dist_CCA <- UniFrac(CCA, weighted = TRUE, normalized = TRUE)

# Ordinate weighted distances

ord_wu_CCA <- ordinate(CCA, "NMDS", distance = wu_dist_CCA)
nmdsdata_wu_CCA <- plot_ordination(CCA, ord_wu_CCA)$data
statusplot_wu_CCA <- plot_ordination(CCA, ord_wu_CCA, color = "Site")
hulls_wu_CCA <- plyr::ddply(statusplot_wu_CCA$data, "Site", find_hull)

# Order sites by geography

nmdsdata_wu_CCA$Site <- 
  factor(
    nmdsdata_wu_CCA$Site, 
    levels = c("Phlipper's Peace", "Angus' Aquarium", "Fish Den", "Overheat", 
               "Pillar Coral", "White Hole", "International Diver",
               "B-S Deep", "Wicked Pissa", "Thalassa Arg", "Leo's Den"))

# Ordinate unweighted distances

ord_uu_CCA <- ordinate(CCA, "NMDS", distance = uu_dist_CCA)
nmdsdata_uu_CCA <- plot_ordination(CCA, ord_uu_CCA)$data
statusplot_uu_CCA <- plot_ordination(CCA, ord_uu_CCA, color = "Site")
hulls_uu_CCA <- plyr::ddply(statusplot_uu_CCA$data, "Site", find_hull)

# Order sites by geography
nmdsdata_uu_CCA$Site <- 
  factor(
    nmdsdata_uu_CCA$Site, 
    levels = c("Phlipper's Peace", "Angus' Aquarium", "Fish Den", "Overheat", 
               "Pillar Coral", "White Hole", "International Diver",
               "B-S Deep", "Wicked Pissa", "Thalassa Arg", "Leo's Den"))

# Plot weighted CCA NMDS

wu_CCA_nmds <- 
  ggplot(nmdsdata_wu_CCA, aes(x = NMDS1, y = NMDS2, col = Site, fill = Site)) +
  geom_point(size = 2, alpha = 0.9) + 
  geom_polygon(data = hulls_wu_CCA, alpha = 0.5) + 
  scale_color_manual(values = site_color) +
  scale_fill_manual(values = site_color) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  ggtitle("Weighted Unifrac")

# Plot unweighted CCA NMDS

uu_CCA_nmds <- 
  ggplot(nmdsdata_uu_CCA, aes(x = NMDS1, y = NMDS2, col = Site, fill = Site)) +
  geom_point(size = 2, alpha = 0.9) + 
  geom_polygon(data = hulls_uu_CCA, alpha = 0.5) + 
  scale_color_manual(values = site_color) +
  scale_fill_manual(values = site_color) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  ggtitle("Unweighted Unifrac")

# Put NMDS plots together

figS7b <- 
  ggarrange(wu_CCA_nmds, uu_CCA_nmds, nrow = 1, 
            common.legend = TRUE, legend = "right")

# Put alpha and beta diversity plots together
ggarrange(figS7a, figS7b, widths = c(1,2.4))

ggsave("FigureS7.tiff", width = 12, height = 4)

#---------Fig S8: Carbonate alpha + beta diversity by site---------

# Remove samples with <2,500 reads

carb_rare <- subset_samples(carb, sample_sums(carb) > 2500)

# Rarefy to lowest read depth

carb_rare <- 
  rarefy_even_depth(
    carb_rare, 
    sample.size = min(sample_sums(carb_rare)), 
    rngseed = 8, 
    replace = FALSE)

# Calculate diversity indices

carb_rare_div <- adiv_df(carb_rare)

# Make short site names

carb_rare_div <-
  carb_rare_div %>% 
  mutate(SiteShort = case_when(
    Site == "Angus' Aquarium" ~ "AA",
    Site == "B-S Deep" ~ "BD",
    Site == "Fish Den" ~ "FD",
    Site == "International Diver" ~ "ID",
    Site == "Leo's Den" ~ "LD",
    Site == "Overheat" ~ "OH",
    Site == "Phlipper's Peace" ~ "PP",
    Site == "Pillar Coral" ~ "PC",
    Site == "Thalassa Arg" ~ "TA",
    Site == "White Hole" ~ "WH",
    Site == "Wicked Pissa" ~ "WP"))

# Order short site factor

carb_rare_div$SiteShort <- 
  factor(carb_rare_div$SiteShort, 
         levels = c("PP", "AA", "FD", "OH", "PC", "WH", 
                    "ID", "BD", "WP", "TA", "LD"))

# Order site factor

carb_rare_div$Site <- 
  factor(
    carb_rare_div$Site, 
    levels = c("Phlipper's Peace", "Angus' Aquarium", "Fish Den", "Overheat", 
               "Pillar Coral", "White Hole", "International Diver",
               "B-S Deep", "Wicked Pissa", "Thalassa Arg", "Leo's Den"))

# Plot carbonate alpha diversity by site

shannon_carb <-
  ggplot(carb_rare_div, aes(x = SiteShort, y = Shannon, fill = Site)) +
  geom_boxplot() +
  scale_fill_manual(values = site_color) +
  theme_bw() +
  ylab("Shannon Diversity") +
  theme(axis.title.x = element_blank()) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

simpson_carb <-
  ggplot(carb_rare_div, aes(x = SiteShort, y = Simpson, fill = Site)) +
  geom_boxplot() +
  scale_fill_manual(values = site_color) +
  theme_bw() +
  ylab("Simpson Diversity (1-D)") +
  theme(axis.title.x = element_blank()) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

figS8a <-
  ggarrange(shannon_carb, simpson_carb, ncol = 1, common.legend = TRUE, 
            legend = "none", align = "v")

# Calculate unweighted and weighted Unifrac distances

uu_dist_carb <- UniFrac(carb, weighted = FALSE, normalized = TRUE)
wu_dist_carb <- UniFrac(carb, weighted = TRUE, normalized = TRUE)

# Ordinate weighted distances

ord_wu_carb <- ordinate(carb, "NMDS", distance = wu_dist_carb)
nmdsdata_wu_carb <- plot_ordination(carb, ord_wu_carb)$data
statusplot_wu_carb <- plot_ordination(carb, ord_wu_carb, color = "Site")
hulls_wu_carb <- plyr::ddply(statusplot_wu_carb$data, "Site", find_hull)

# Order sites by geography

nmdsdata_wu_carb$Site <- 
  factor(
    nmdsdata_wu_carb$Site, 
    levels = c("Phlipper's Peace", "Angus' Aquarium", "Fish Den", "Overheat", 
               "Pillar Coral", "White Hole", "International Diver",
               "B-S Deep", "Wicked Pissa", "Thalassa Arg", "Leo's Den"))

# Ordinate unweighted distances

ord_uu_carb <- ordinate(carb, "NMDS", distance = uu_dist_carb)
nmdsdata_uu_carb <- plot_ordination(carb, ord_uu_carb)$data
statusplot_uu_carb <- plot_ordination(carb, ord_uu_carb, color = "Site")
hulls_uu_carb <- plyr::ddply(statusplot_uu_carb$data, "Site", find_hull)

# Order sites by geography
nmdsdata_uu_carb$Site <- 
  factor(
    nmdsdata_uu_carb$Site, 
    levels = c("Phlipper's Peace", "Angus' Aquarium", "Fish Den", "Overheat", 
               "Pillar Coral", "White Hole", "International Diver",
               "B-S Deep", "Wicked Pissa", "Thalassa Arg", "Leo's Den"))

# Plot weighted carbonate NMDS

wu_carb_nmds <- 
  ggplot(nmdsdata_wu_carb, aes(x = NMDS1, y = NMDS2, col = Site, fill = Site)) +
  geom_point(size = 2, alpha = 0.9) + 
  geom_polygon(data = hulls_wu_carb, alpha = 0.5) + 
  scale_color_manual(values = site_color) +
  scale_fill_manual(values = site_color) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  ggtitle("Weighted Unifrac")

# Plot unweighted carbonate NMDS

uu_carb_nmds <- 
  ggplot(nmdsdata_uu_carb, aes(x = NMDS1, y = NMDS2, col = Site, fill = Site)) +
  geom_point(size = 2, alpha = 0.9) + 
  geom_polygon(data = hulls_uu_carb, alpha = 0.5) + 
  scale_color_manual(values = site_color) +
  scale_fill_manual(values = site_color) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  ggtitle("Unweighted Unifrac")

# Put NMDS plots together

figS8b <- 
  ggarrange(wu_carb_nmds, uu_carb_nmds, nrow = 1, 
            common.legend = TRUE, legend = "right")

# Put alpha and beta diversity plots together
ggarrange(figS8a, figS8b, widths = c(1,2.4))

ggsave("FigureS8.tiff", width = 12, height = 4)
