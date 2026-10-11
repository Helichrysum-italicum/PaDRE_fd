# ==============================================================================
# 0. SETUP & LIBRARIES
# ==============================================================================
library(tidyverse)
library(ggplot2)
library(sf)
library(rWCVP)
library(rWCVPdata)
library(patchwork)
library(viridis)

# Set base paths
setwd("xxx")


# ==============================================================================
# 1. TRY DATABASE COVERAGE (BARPLOT)
# ==============================================================================
# Load TRY data
try_data <- read.csv(paste0("0_trait_dataset/try_wide_output.csv"))

# Rename columns to manageable names
try_data <- try_data %>% rename(
  "Mycorrhiza" = "Mycorrhiza.type",
  "N_fixation" = "Plant.nitrogen.N..fixation.capacity",
  "Growth_form" = "Plant.growth.form",
  "Dispersal" = "Dispersal.syndrome",
  "Flowering_time" = "Plant.reproductive.phenology.timing..flowering.time.",
  "Pollination" = "Pollination.syndrome",
  "Growth_form_det" = "Plant.growth.form.detailed.consolidated",
  "Height" = "Plant.height.vegetative_m",
  "Seed_mass" = "Seed.dry.mass_mg",
  "Wood_density" = "Stem.specific.density..SSD..stem.dry.mass.per.stem.fresh.volume..or.wood.density_g.cm3",
  "Leaf_thickness" = "Leaf.thickness_mm",
  "Leaf_N_mass" = "Leaf.nitrogen..N..content.per.leaf.dry.mass_mg.g",
  "LDMC" = "Leaf.dry.mass.per.leaf.fresh.mass..leaf.dry.matter.content..LDMC._g.g",
  "Germination_rate" = "Seed.germination.rate..germination.efficiency._.",
  "Leaf_N_area" = "Leaf.nitrogen..N..content.per.leaf.area_g.m.2",
  "Chlorophyll_area" = "Leaf.chlorophyll.content.per.leaf.area_g.m.2",
  "Seed_length" = "Seed.length_mm",
  "Root_depth" = "Root.rooting.depth_m",
  "Chlorophyll_SPAD" = "Leaf.chlorophyll.content.per.leaf.area_SPAD.Chlorophyll.Index",
  "Fruit_length" = "Fruit.length_mm",
  "Fruit_mass" = "Fruit.dry.mass_mg",
  "Chlorophyll_CCM" = "Leaf.chlorophyll.content.per.leaf.area_CCM.Chlorophyll.Content.Index",
  "Chlorophyll_FCI" = "Leaf.chlorophyll.content.per.leaf.area_Falker.Chlorofilog.Chlorophyll.Content.Index.FCI",
  "Chlorophyll_CCI" = "Leaf.chlorophyll.content.per.leaf.area_CCI",
  "Leaf_area" = "Leaf.area..in.case.of.compound.leaves..leaflet..undefined.if.petiole.is.in..or.excluded._mm2",
  "SLA" = "Leaf.area.per.leaf.dry.mass..specific.leaf.area..SLA.or.1.LMA...undefined.if.petiole.is.in..or.excluded_mm2.mg.1"
)

# Calculate coverage percentage for selected columns
plot_data_try <- try_data %>%
  select(Height:SLA) %>%
  pivot_longer(everything(), names_to = "Trait", values_to = "Value") %>%
  group_by(Trait) %>%
  summarise(
    Present = sum(!is.na(Value)) / n(),
    Missing = sum(is.na(Value)) / n()
  ) %>%
  pivot_longer(c(Present, Missing), names_to = "Status", values_to = "Percentage")

# Generate horizontal bar plot
try_numerical_coverage <- ggplot(plot_data_try, aes(x = Trait, y = Percentage, fill = Status)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_fill_manual(values = c("Present" = "#008080", "Missing" = "#BDC3C7")) +
  labs(title = "Coverage by numerical trait (TRY)", x = "Trait", y = "Share (0-1)") +
  theme_minimal()

# Save PDF
ggsave(filename = paste0("accessory_visual/1_try_numerical_coverage.pdf"), 
       plot = try_numerical_coverage, width = 8, height = 6, device = "pdf")

# ==============================================================================
# 2. IMPUTATION SPECIES POOL COVERAGE (BARPLOT)
# ==============================================================================
# Define target columns for traits
cols_to_plot <- c("Seed_mass_mg", "SLA_mm2_mg", "LDMC_g_g", "Max_plant_height_m")

# Load imputation species pool
species_pool <- readRDS(paste0("species_pool_for_imputation.rds"))

# Filter to keep only target columns and drop empty ones
species_pool <- species_pool %>%
  select(any_of(cols_to_plot)) %>%
  select(where(~ !all(is.na(.)) && !all(. == 0, na.rm = TRUE)))

# Calculate coverage percentage
plot_data_pool <- species_pool %>%
  pivot_longer(everything(), names_to = "Trait", values_to = "Value") %>%
  group_by(Trait) %>%
  summarise(
    Present = sum(!is.na(Value)) / n(),
    Missing = sum(is.na(Value)) / n()
  ) %>%
  pivot_longer(c(Present, Missing), names_to = "Status", values_to = "Percentage")

# Generate vertical bar plot
species_pool_missing <- ggplot(plot_data_pool, aes(x = Trait, y = Percentage, fill = Status)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Present" = "#008080", "Missing" = "#BDC3C7")) +
  labs(title = "Imputation species pool data coverage", x = "Trait", y = "Proportion (0-1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save PDF
ggsave(filename = paste0("accessory_visual/2_species_pool_coverage.pdf"), 
       plot = species_pool_missing, width = 6, height = 6, device = "pdf")

# ==============================================================================
# 3. MERGED DATASET COVERAGE (BARPLOT)
# ==============================================================================
# Load merged database
all_db <- readRDS(paste0("0_trait_dataset/GIFT_try_bien_merged.RDS"))

# Calculate coverage percentage
plot_data_all <- all_db %>%
  select(all_of(cols_to_plot)) %>%
  pivot_longer(everything(), names_to = "Trait", values_to = "Value") %>%
  group_by(Trait) %>%
  summarise(
    Present = sum(!is.na(Value)) / n(),
    Missing = sum(is.na(Value)) / n()
  ) %>%
  pivot_longer(c(Present, Missing), names_to = "Status", values_to = "Percentage")

# Generate vertical bar plot
all_datasets_missing <- ggplot(plot_data_all, aes(x = Trait, y = Percentage, fill = Status)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Present" = "#008080", "Missing" = "#BDC3C7")) +
  labs(title = "Trait Data Coverage (Merged Database)", x = "Trait", y = "Proportion (0-1)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save PDF
ggsave(filename = paste0("accessory_visual/3_merged_database_coverage.pdf"), 
       plot = all_datasets_missing, width = 6, height = 6, device = "pdf")

# ==============================================================================
# 4. SPATIAL COVERAGE MAPPING (rWCVP)
# ==============================================================================
# Load global species distributions (native, non-extinct only)
# We use distributions instead of checklist to have the specific area_code_l3
traits_data <- readRDS("C:/Users/nicosia/Desktop/PaDRE_FD/data/0_trait_dataset/GIFT_try_bien_merged.RDS")

checklist <- wcvp_checklist(taxon_rank = "species",
                            synonyms = FALSE,
                            native = TRUE, 
                            introduced = FALSE,
                            extinct = FALSE, 
                            hybrids = FALSE,
                            infraspecies = FALSE
)

data_with_traits <- checklist %>%
  left_join(traits_data, by = c( "taxon_name" = "species_name"))

# Calculate percentage of non-NA values per trait per Botanical Country (Level 3)
coverage_data <- data_with_traits %>%
  group_by(area_code_l3) %>%
  summarise(
    Seed_mass_cov = sum(!is.na(Seed_mass_mg)) / n() * 100,
    SLA_cov = sum(!is.na(SLA_mm2_mg)) / n() * 100,
    LDMC_cov = sum(!is.na(LDMC_g_g)) / n() * 100,
    Height_cov = sum(!is.na(Max_plant_height_m)) / n() * 100,
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(Mean_coverage = mean(c(Seed_mass_cov, SLA_cov, LDMC_cov, Height_cov), na.rm = TRUE)) %>%
  ungroup()

# Load WGSRPD Level 3 spatial polygons and join with our coverage metrics
wgsrpd3_geo <- rWCVPdata::wgsrpd3
map_data <- wgsrpd3_geo %>%
  left_join(coverage_data, by = c("LEVEL3_COD" = "area_code_l3"))

# Helper function to generate a consistent map plot
create_trait_map <- function(data, column, title) {
  ggplot(data) +
    geom_sf(aes(fill = .data[[column]]), color = "white", linewidth = 0.1) +
    scale_fill_viridis_c(
      option = "viridis", 
      na.value = "grey90", 
      name = "% Cov", 
      limits = c(0, 100) 
    ) +
    theme_void() +
    labs(title = title) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(b = 5)),
      legend.position = "bottom",       
      legend.key.width = unit(1.2, "cm"),
      legend.key.height = unit(0.3, "cm"),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
}

# Generate individual trait maps
p_seed <- create_trait_map(map_data, "Seed_mass_cov", "Seed Mass")
p_sla <- create_trait_map(map_data, "SLA_cov", "Specific Leaf Area")
p_ldmc <- create_trait_map(map_data, "LDMC_cov", "Leaf Dry Matter Content")
p_height <- create_trait_map(map_data, "Height_cov", "Plant Height")

# Combine the 4 traits into a single layout using Patchwork
combined_traits_plot <- (p_seed + p_sla) / (p_ldmc + p_height) +
  plot_annotation(
    title = "Trait Coverage Percentage by Botanical Country",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

# Create standalone map for the mean coverage
mean_coverage_plot <- create_trait_map(map_data, "Mean_coverage", "Mean Coverage (4 Traits)") +
  theme(plot.title = element_text(size = 16))

# Save Patchwork PDF (needs larger dimensions for 2x2 grid)
ggsave(filename = paste0("accessory_visual/4_spatial_traits_patchwork.pdf"), 
       plot = combined_traits_plot, width = 14, height = 10, device = "pdf")

# Save Mean Coverage PDF
ggsave(filename = paste0("accessory_visual/5_spatial_mean_coverage.pdf"), 
       plot = mean_coverage_plot, width = 10, height = 6, device = "pdf")
