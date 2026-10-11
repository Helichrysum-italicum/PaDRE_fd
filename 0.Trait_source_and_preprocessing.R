# ==============================================================================
# PROJECT:     PaDRE_FD (Functional Diversity)
# AUTHOR:      Samuele Nicosia
# DATE:        2026-04-20
# DESCRIPTION: 
# load functional traits data from TRY-db and pre-process them
# following steps suggested by RTY paper from Ying Lam
# DOI: 10.1002/ece3.11292 
# and biological limits defined by the handbook of traits Pérez-Harguindeguy
# http://dx.doi.org/10.1071/BT12225_CO
# ==============================================================================


# ==============================================================================
                        # IMPORT AND PRE_PROCESS TRY DATA
# ==============================================================================
# ---- 0. LIBRARIES & SETUP ----

library(rtry)
library(tidyverse)
library(dplyr)
library(GIFT)
library(BIEN)


#set working directory

setwd("xxxx")


# ---- 1. LOAD TRY TRAITS ----

#load all your raw traits and bind in a single table for easier management
set1 <- rtry_import("0_trait_dataset/44430.txt")
set2 <- rtry_import("0_trait_dataset/47839.txt ")
traits <- rbind(set1, set2)

#remove imported datasets and perform garbage collection to free up memory
rm(set1)
rm(set2)
gc()


# ---- 1.5 CHECK TRAITS AND FILTER OUTLIERS ----


#check traits inside your dataset
traitid <- rtry_explore(traits, TraitID, TraitName)


# filter outliers based on error risk (standard deviation)
traits_filtered <- rtry_exclude(traits,
                                ErrorRisk > 3,
                                baseOn = ObsDataID)

# remove juvenile observation and filtering all possible misspelling
traits_filtered <- rtry_exclude(
  traits_filtered,
  (DataID == 413) & grepl("juvenil|junv|juve|seedl|sapl|young|immature|weeks|month", 
                          OrigValueStr, 
                          ignore.case = TRUE),
  baseOn = ObservationID
)


# ---- 2. SPLIT TRAITS BY DATA TYPE ----

numerical  <- c(4, 6, 14, 26, 27, 46, 47, 50, 95, 413, 918, 919, 3113, 3117, 3568)
# 4: Stem specific density # 6: Root rooting depth
# 14: Leaf N (mass) # 26: Seed dry mass
# 27: Seed length # 46: Leaf thickness
# 47: LDMC # 50: Leaf N (area)
# 95: Seed germination rate # 413: Chlorophyll content
# 918: Fruit length # 919: Fruit dry mass
# 3113: Leaf area # 3117: SLA # 3568: Flower size

height <- c(3106)
# 3106: Plant height

categorical  <- c(7, 8, 28, 29, 42, 335, 3401)
# 7: Mycorrhiza type # 8: Plant nitrogen(N) fixation capacity
# 28: Dispersal syndrome # 29: Pollination syndrome
# 42: Plant growth form # 335: Plant reproductive phenology timing
# 3401: Plant growth form detailed consolidated

#filter traits using IDs
data_num <- rtry_select_row(traits_filtered, TraitID %in% numerical)
data_cat <- rtry_select_row(traits_filtered, TraitID %in% categorical)
data_height <- rtry_select_row(traits_filtered, TraitID %in% height)


#remove NA and invalid data 
data_cat <- rtry_select_col(data_cat, AccSpeciesName,TraitName,OrigValueStr)%>%
  filter(!is.na(OrigValueStr))

data_height <- rtry_select_col(data_height, AccSpeciesName,TraitName,StdValue,UnitName)%>%
  filter(!is.na(StdValue))

data_num <- rtry_select_col(data_num, AccSpeciesName,TraitName,StdValue,UnitName)%>%
  filter(!is.na(StdValue))

# ---- 3. FILTER USING BIOLOGICAL LIMITS ----

#keep only data within biologically possible values

data_height <- data_height%>%
  filter(
    case_when(
      TraitName == "Plant height vegetative" ~ (StdValue > 0 & StdValue <= 140),
      TRUE ~ TRUE ))

data_num <- data_num %>%
  mutate(StdValue = as.numeric(StdValue)) %>%
  filter(
    case_when(
      TraitName == "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded" ~ (StdValue > 0 & StdValue <= 300),
      TraitName == "Leaf dry mass per leaf fresh mass (leaf dry matter content, LDMC)" ~ (StdValue >= 0.03 & StdValue <= 0.80),
      TraitName == "Seed dry mass" ~ (StdValue > 0 & StdValue <= 10000000),
      TRUE ~ TRUE 
    )
  )


# ---- 4. CALCULATE SPECIES LEVEL MEANS  ----

#function get_p95 obtain the 95 quantile of the values at species level 
#avoid including non-labeled juvenile records for plant height
#and outliers within global biological boundaries

get_p95 <- function(x) {
  x <- x[!is.na(x) & x > 0] 
  if (length(x) == 0) return(NA) 
  quantile(x, probs = 0.95, names = FALSE)
}
#obtain the mode for categorical values where mean cannot be calculated

get_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#use geometric mean for numerical values as adviced in 
#Diaz 2016 global spectrum of plants form and functions

get_geomean <- function(x) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) == 0) return(NA)
  exp(mean(log(x)))
}

#compute the operation using try_trans wider functions above are in values_fn
data_num_wide <- rtry_trans_wider(data_num,
                                  names_from = c(TraitName,UnitName),
                                  values_from = c(StdValue),
                                  values_fn = list(StdValue = get_geomean))

data_height_wide <- rtry_trans_wider(data_height,
                                     names_from = c(TraitName,UnitName),
                                     values_from = c(StdValue),
                                     values_fn = list(StdValue = get_p95))

data_cat_wide <- rtry_trans_wider(data_cat,
                                  names_from = c(TraitName),
                                  values_from = OrigValueStr,
                                  values_fn = list(OrigValueStr = get_mode))



#join back all data type to complete the species-by-trait matrix
data_final <- rtry_join_outer(data_cat_wide, data_height_wide, baseOn = AccSpeciesName)
data_final <-  rtry_join_outer(data_final, data_num_wide, baseOn = AccSpeciesName)

#save the raw species-mean level matrix for all traits
write_csv(data_final, "0_trait_dataset/try_wide_output.csv")

#select traits with high coverage for downstream operations
data_selected <- data_final %>%
  select(
    species_name      = `AccSpeciesName`,
    growth_form       = `Plant growth form detailed consolidated`,
    Height_m          = `Plant height vegetative_m`,
    SLA_mm2_mg        = `Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded_mm2 mg-1`,
    LDMC_g_g          = `Leaf dry mass per leaf fresh mass (leaf dry matter content, LDMC)_g/g`,
    Seed_mass_mg      = `Seed dry mass_mg`,
    Leaf_area_mm2     = `Leaf area (in case of compound leaves: leaflet, undefined if petiole is in- or excluded)_mm2`,
    Stem_density_g_cm3 = `Stem specific density (SSD, stem dry mass per stem fresh volume) or wood density_g/cm3`,
    Leaf_N_mg_g       = `Leaf nitrogen (N) content per leaf dry mass_mg/g`
  )


write_csv(data_selected, "0_trait_dataset/try_wide_output_selected.csv")


# ==============================================================================
                      # DOWNLOAD, JOIN AND MERGE GIFT DATA
# ==============================================================================
# ---- 5. IMPORT GIFT DATA ----

meta <- GIFT_traits_meta()
#	Woodiness 1.1.1
# Climber 1.4.1
# Climber 1.4.2
# Parasite 1.5.1
# Epiphyte 1.3.1
# Aquaticity 1.7.1
# Growth form "simple" 1.2.1
# Growth form "expanded" 1.2.2
# raunkiaer life forms Life form "expanded" 2.3.2
# raunkiaer life forms Life form "simple" 2.3.1
# 6.1.2 Elevational_range_max
# 6.1.1 Elevational_range_min
# Habitat 6.3.1
# max plant height 1.6.2
# Seed mass 3.2.3
# Specific Leaf Area (SLA) 4.1.3
# Leaf_length_max 4.6.2
# Leaf_width_max 4.7.2

GIFT_traits <- GIFT_traits(
  trait_IDs = c(
    "1.1.1", "1.4.1","1.4.2", "1.5.1", "1.3.1", "1.7.1",
    "1.2.1", "1.2.2", "2.3.2", "2.3.1",         
    "6.1.2", "6.1.1", "6.3.1", "1.6.2", "3.2.3", "4.1.3",
    "4.6.2", "4.7.2"
  ),
  agreement = 0.66,
  bias_ref = TRUE,
  bias_deriv = TRUE,
  api = "https://gift.uni-goettingen.de/api/extended/",
  GIFT_version = "latest"
)

GIFT_traits <- GIFT_traits %>%
  rename(
    species_name = work_species,
    Woodiness = trait_value_1.1.1,
    Climber = trait_value_1.4.1,
    Climber_type = trait_value_1.4.2,
    Parasite = trait_value_1.5.1,
    Epiphyte = trait_value_1.3.1,
    Aquaticity = trait_value_1.7.1,
    Growth_form_simple = trait_value_1.2.1,
    Growth_form_expanded = trait_value_1.2.2,
    Raunkiaer_life_forms_expanded = trait_value_2.3.2,
    Raunkiaer_life_forms_simple = trait_value_2.3.1,
    Elevational_range_max = trait_value_6.1.2,
    Elevational_range_min = trait_value_6.1.1,
    Habitat = trait_value_6.3.1,
    Max_plant_height = trait_value_1.6.2,
    Seed_mass = trait_value_3.2.3,
    SLA = trait_value_4.1.3,
    Leaf_width_max = trait_value_4.7.2,
    Leaf_length_max = trait_value_4.6.2
  )


GIFT_traits <- GIFT_traits %>%
  select(
    species_name,
    Woodiness,
    Climber,
    Climber_type,
    Parasite,
    Epiphyte,
    Aquaticity,
    Growth_form_simple,
    Growth_form_expanded,
    Raunkiaer_life_forms_expanded,
    Raunkiaer_life_forms_simple,
    Elevational_range_max,
    Elevational_range_min,
    Habitat,
    Max_plant_height,
    Seed_mass,
    SLA)%>% 
  mutate(
    Seed_mass = Seed_mass * 1000, #mutate unit in mg
    SLA       = SLA /10 # mutate unit in mm2_mg
  )

#saveRDS(GIFT_traits, file="trait_datasets/GIFT_raw_traits.RDS")


# ---- 6. MERGE GIFT AND TRY----

try_and_GIFT <- GIFT_traits %>%
  full_join(
    data_selected %>% select(species_name, Height_m, Seed_mass_mg, SLA_mm2_mg, LDMC_g_g), 
    by = "species_name"
  ) %>%
  mutate(
    Max_plant_height_m = coalesce(Max_plant_height, Height_m),
    Seed_mass_mg       = coalesce(Seed_mass, Seed_mass_mg),
    SLA_mm2_mg         = coalesce(SLA, SLA_mm2_mg)
  ) %>%
  select(-Height_m, -Seed_mass, -SLA, -Max_plant_height)


#saveRDS(try_and_GIFT, file="trait_datasets/GIFT_try_merged.RDS")
#try_and_GIFT <- readRDS("trait_datasets/GIFT_try_merged.RDS")

# ==============================================================================
#                         DOWNLOAD, JOIN AND MERGE BIEN DATA
# ==============================================================================

# ---- 8. DOWNLOAD BIEN DATA ----


# Download specific trait data for all species

height_data       <- BIEN_trait_trait(trait = "whole plant height")
SLA_data          <- BIEN_trait_trait(trait = "leaf area per leaf dry mass")
LDMC_data         <- BIEN_trait_trait(trait = "leaf dry mass per leaf fresh mass")


save(height_data,
     SLA_data,
     LDMC_data,
     file = "0_trait_dataset/BIEN_traits_data.RData")

#load("trait_datasets/BIEN_traits_data.RData")


# ---- 9. PREPROCESS-SUMMARIZE-UNIFORM UNITS FOR BIEN DATA ----
#define as numeric
height_data$trait_value <- as.numeric(height_data$trait_value)
#filter extreme outliers
height_data <- subset(height_data, trait_value > 0 & trait_value <= 140)

#group by name and summarise
max_height_data <- height_data %>%
  filter(!is.na(scrubbed_species_binomial)) %>%
  group_by(scrubbed_species_binomial) %>%
  summarise(max_height_m = quantile(trait_value, probs = 0.95, na.rm = TRUE))

#free up memory
rm(height_data)

#define as numeirc
SLA_data$trait_value <- as.numeric(SLA_data$trait_value)

#filter extreme outliers
SLA_data <- subset(SLA_data, trait_value >= 0 & trait_value <= 300)

#summarise
mean_SLA <- SLA_data %>%
  filter(!is.na(scrubbed_species_binomial)) %>%
  group_by(scrubbed_species_binomial) %>%
  summarise(mean_SLA_m2_kg = exp(mean(log(trait_value), na.rm = TRUE)))

#free up memory
rm(SLA_data)

#define as numeric
LDMC_data$trait_value <- as.numeric(LDMC_data$trait_value)

#convert unit to match TRY and GIFT
LDMC_data$trait_value <- LDMC_data$trait_value / 1000

#group by specie and summarise
mean_LDMC <- LDMC_data %>%
  filter(!is.na(scrubbed_species_binomial)) %>%
  group_by(scrubbed_species_binomial) %>%
  summarise(mean_LDMC_g_g = exp(mean(log(trait_value), na.rm = TRUE)))
#free up memory
rm(LDMC_data)

# ---- 10. MERGE WITH GIFT AND TRY ----

try_GIFT_BIEN <- try_and_GIFT %>%
  full_join(max_height_data, by = c("species_name" = "scrubbed_species_binomial")) %>%
  full_join(mean_SLA, by = c("species_name" = "scrubbed_species_binomial")) %>%
  full_join(mean_LDMC, by = c("species_name" = "scrubbed_species_binomial")) %>%
  #fill voids
  mutate(
    Max_plant_height_m = coalesce(Max_plant_height_m, max_height_m),
    SLA_mm2_mg         = coalesce(SLA_mm2_mg, mean_SLA_m2_kg),
    LDMC_g_g           = coalesce(LDMC_g_g, mean_LDMC_g_g)
  ) %>%
  select(-max_height_m, -mean_SLA_m2_kg, -mean_LDMC_g_g)

#remove zeros 
try_GIFT_BIEN <- try_GIFT_BIEN %>%
  mutate(Seed_mass_mg = na_if(Seed_mass_mg, 0))

#save file
saveRDS(try_GIFT_BIEN, file="0_trait_dataset/GIFT_try_bien_merged.RDS")

