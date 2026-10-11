# ==============================================================================
# PROJECT:     PaDRE_FD (Functional Diversity)
# AUTHOR:      Samuele Nicosia
# DATE:        2026-04-20
# DESCRIPTION: 
# Downloads native vascular plant checklists from GIFT database
# Filter native flora and subset it by nation creating a list object
# join to the list of each country the species extirpated within the county
# ==============================================================================

# ---- 0. LIBRARIES & SETUP ----
library(GIFT)
library(dplyr)
library(purrr)
library(tidyr)
library(purrr)
library(dplyr)
library(sf)
library(rnaturalearth)
library(countrycode)
library(GIFT)
library(giscoR)


setwd("xxx")


# ---- 1. SELECT TARGET COUNTRIES AND PREPARE SHAPFILES ----

# Target list 
target_countries <- c(
  "Serbia", "Macedonia",  "Austria", "Belgium",
  "Bulgaria","Moldova", "Cyprus","France", "Czech Republic", "United Kingdom", "Iceland",
   "Germany", "Luxembourg", "Netherlands", "Switzerland", "Denmark", "Estonia",
  "Finland", "Ireland", "Latvia", "Lithuania", "Norway", "Sweden", "Albania",
  "Croatia", "Greece", "Italy", "Malta", "Portugal", "Spain", "Poland",
  "Slovakia", "Ukraine", "Hungary", "Romania", "Belarus","Slovenia" )


#name each element of the vector with the name contained
target_countries <- set_names(target_countries)

#download GIFT taxonomy
taxonomy <- GIFT_taxonomy()%>%
  mutate(genus = taxon_name)

#download GIFT shapefile
shapes <- GIFT_shapes()
# saveRDS(shapes, file="data/shapesGIFT.RDS")
# shapes <- readRDS("data/shapesGIFT.RDS")

# ---- 2. COUNTRY EXTRACTION & LIST CREATION ----
country_checklists <- imap(target_countries, function(country, list_name) {
  
  
  country_shape <- shapes %>% 
    filter( geo_entity == country)%>%
    select(geometry)
  
  #safety check to visualize if the correct geomtry has been exported
  #NB: for france and netherlands dependencies might be included in the shape
  #check it by plotting the shape and checking with your eyes
  
  #plot(st_geometry(country_shape), main = paste("Shapefile of:", country), col = "turquoise4")
  
  # Download the checklist for the country shape from GIFT
  
  #request native flora
  gift_data <- GIFT_checklists(
    geo_type = "All",
    remove_overlap = FALSE,
    shp = country_shape,
    floristic_group = "native",
    overlap = "centroid_inside",
    api = "https://gift.uni-goettingen.de/api/extended/",
    GIFT_version = "latest"
  )
  
  
  #manage croatia differently due to the centroid of the country
  #falling out of the country borders
  
  # ---- CROATIA: filter by manually validated list_IDs ----
  if (country == "Croatia") {
    
    gift_data <- GIFT_checklists(
      geo_type = "All",
      remove_overlap = FALSE,
      floristic_group = "native",
      shp = country_shape,
      overlap = "shape_intersect",
      api = "https://gift.uni-goettingen.de/api/extended/",
      GIFT_version = "latest"
    )
    
    croatia_list_ids <- c(14533, 14115, 10353, 13141, 13167, 13206, 12585)
    gift_data[[2]] <- gift_data[[2]] %>%
      filter(list_ID %in% croatia_list_ids)}
  
  #----- END OF CROATIA MANUAL MANAGEMENT ------
  
  
  # Join genus to the species list
  gift_data[[2]] <- gift_data[[2]] %>%
    left_join(
      taxonomy %>% select(genus, taxon_ID),
      by = c("genus_ID" = "taxon_ID")
    )
  
  # merge all checklists produced by gift together and remove duplicates
  clean_checklist <- gift_data[[2]] %>%
    distinct(work_species, .keep_all = TRUE) %>%
    filter(!is.na(work_species)) %>%
    filter(native == 1) %>%
    #select columns to retain
    select(family, genus, species_name = work_species, author = work_author)
  
  
  return(clean_checklist)
})


#save as rds and re load if necessary
# saveRDS(country_checklists, file = "data/countries_checklists.rds")
# country_checklists <- readRDS("data/countries_checklists.rds")


# ---- 3. JOIN EXTIRPATED SPECIES  ----

#load extirpated species list and red list categorized species 
ext_species <- read.csv("0_extirpation_list/extinct_species_RL_WCVP_updated.csv", fileEncoding = "latin1")
all_species <- read.csv("0_extirpation_list/rl_combined_REDLIST_updated.csv", fileEncoding = "latin1")


# rename columns and add extirpated status
ext_clean <- ext_species %>%
  select(Country,family, genus, species_name = accepted_taxon_name, author = accepted_taxon_authors) %>%
  mutate(status = "extirpated")


#select only necessary columns
rl_clean <- all_species %>%
  select(Country, family, genus, species_name = accepted_taxon_name, author = accepted_taxon_authors, Redlist_cat)

rm(ext_species)
rm(all)

# join extirpated species to each country subset and add redlist cat for all categorized species
complete_checklists <- imap(country_checklists, function(.x, .y) {
  
  current_country <- .y
  
  # Subset and de-duplicate Red List data (keep only species_name and category)
  country_rl_data <- rl_clean %>% 
    filter(Country == current_country) %>%
    distinct(species_name, .keep_all = TRUE) %>%
    select(species_name, rl_cat = Redlist_cat)
  
  # Subset and de-duplicate Extirpated data 
  #remove country to avoid conflicts
  country_ex_data <- ext_clean %>% 
    filter(Country == current_country) %>%
    distinct(species_name, .keep_all = TRUE) %>%
    select(-Country)
  
  # join data to the national flora checklist
  joined_df <- .x %>%
    full_join(country_ex_data, by = c("family", "genus", "species_name", "author")) %>%
    left_join(country_rl_data, by = "species_name")  
  
  
  final_checklist <- joined_df %>%
    mutate(
      status = replace_na(status, "extant"),
      
      # Apply and join back status 
      Redlist_cat = case_when(
        status == "extirpated" ~ "extirpated",
        !is.na(rl_cat)          ~ rl_cat,
        TRUE                    ~ "not assessed"
      )
    ) %>%
    
    select(-any_of(c("rl_cat")))
  
  return(final_checklist)
})


# saveRDS(complete_checklists, file = "countries_checklists_and_extinct.rds")
# complete_checklists <- readRDS("countries_checklists_and_extinct.rds")

# ---- 6. LOAD PRE-PROCESSED FUNCTIONAL TRAITS ----

traits_wide <- readRDS("0_trait_dataset/GIFT_try_bien_merged.RDS")

#map traits to each country's checklist
checklists_with_traits <- map(complete_checklists, function(.x) {
  
  .x %>%
    left_join(traits_wide, by = "species_name")
})



# ---- 9. CALCULATE PERCENTAGE OF MISSING TRAITS BY STATUS ----

#calculate missing summary via gruop by status

missing_traits_summary <- map(checklists_with_traits, function(.x) {
  
  .x %>%
    # group and calculate missing % by 
    group_by(status) %>%
    summarise(
      across(
        .cols = c(Max_plant_height_m, Seed_mass_mg, SLA_mm2_mg, LDMC_g_g), 
        .fns = ~ round(mean(is.na(.)) * 100, 2)
      ),
      .groups = "drop"
    ) %>%
    
    # calculate the overral mean missing
    rowwise() %>%
    mutate(
      average = round(mean(c_across(-status), na.rm = TRUE), 2)
    ) %>%
    ungroup() %>%
    
    # Pivot to have one row for each country
    pivot_wider(
      names_from = status,
      values_from = !status, 
      names_glue = "missing_{.value}_{status}"
    )
})


# ---- 10. CONVERT LIST TO A SINGLE GLOBAL DATAFRAME ----
global_missing_traits_df <- missing_traits_summary %>%
  list_rbind(names_to = "Country")

#  filter countries with mean missing < 50%

selected_countries_df <- global_missing_traits_df %>%
 filter(missing_average_extant < 50 & missing_average_extirpated < 50 )

selected_countries_list <- checklists_with_traits[unique(selected_countries_df$Country)]
saveRDS(selected_countries_list, file = "selected_countries_for_imputation.rds")



# ---- 10. BUILD SPECIES POOL FOR GAP_FILLING AT CONTINENT LEVEL ----

# 1. base pool with countries whith mean missing <50%
base_pool <- selected_countries_list %>%
  list_rbind() %>%
  distinct(species_name, .keep_all = TRUE) 

# 2. excluded countries
dropped_countries <- setdiff(names(checklists_with_traits), unique(selected_countries_df$Country))

# 3. Excluded species with traits to strentghen phylogenetic signals
extra_species_pool <- checklists_with_traits[dropped_countries] %>%
  list_rbind() %>%
  filter(rowSums(!is.na(select(., Max_plant_height_m, Seed_mass_mg, SLA_mm2_mg, LDMC_g_g))) >= 1) %>%
  distinct(species_name, .keep_all = TRUE)

continent_pool <- bind_rows(base_pool, extra_species_pool)%>%
  distinct(species_name, .keep_all = TRUE)

# Export the final nested list containing Taxonomy, Status, Red List and Functional Traits
saveRDS(continent_pool, file = "species_pool_for_imputation.rds")
