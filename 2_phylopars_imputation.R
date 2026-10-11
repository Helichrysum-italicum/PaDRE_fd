# ==============================================================================
# PROJECT:     PaDRE_FD (Functional Diversity)
# AUTHOR:      Samuele Nicosia
# DATE:        2026-04-20
# DESCRIPTION: 
# use GBOT and phylomaker to produce a phylogenetic tree for the species pool
# use phylopars to impute (gap fill) missing data and re-subset into national lists
# ==============================================================================

library(Rphylopars)
library(tidyverse)
library(V.PhyloMaker2)
library(stringr)
library(adespatial)
library(adephylo)
library(ape)
library(reshape2)
library(gt)


setwd("xxx")
continent_pool <- readRDS(file = "species_pool_for_imputation.rds")



# ---- 1. PREPARE TAXONOMY  ----

taxonomy <- continent_pool %>%
  select( family, genus, species_name)%>%
  rename( species = species_name)

taxonomy <- taxonomy %>%
  relocate(species, genus, family)

taxonomy <- taxonomy %>%
  mutate(family = case_when(
    genus == "Comandra"    ~ "Comandraceae",
    genus == "Corrigiola"  ~ "Caryophyllaceae",
    genus == "Cyrtomium"   ~ "Dryopteridaceae",
    genus == "Cystopteris" ~ "Cystopteridaceae",
    genus == "Diplazium"   ~ "Athyriaceae",
    genus == "Metaxya"     ~ "Metaxyaceae",
    genus == "Woodsia"     ~ "Woodsiaceae",
    genus == "Lapiedra" ~ "Amaryllidaceae",
    TRUE                   ~ family 
  ))


# ---- 2. BUILD PHYLOGENETIC TREE  ----

#phylogenetic tree using all default setting
# 
#  result <- phylo.maker(
#    sp.list = taxonomy,
#    tree = GBOTB.extended.TPL,
#    nodes = nodes.info.1.TPL,
#    output.sp.list = TRUE,
#    output.tree = FALSE,
#    scenarios = "S3",
#    r = 1)
# 

# saveRDS(result, file="phylomakertree.RDS")

result <- readRDS("phylomakertree.RDS")

#export result for the tree
tree <- result$scenario.3


# ---- 3. PREPARE TRAITS FOR IMPUTATION PROCESS  ----

traits_for_imputation <- continent_pool %>%
  rename(species = species_name) %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
    select(species ,Max_plant_height_m,Seed_mass_mg, SLA_mm2_mg, LDMC_g_g) %>%
    relocate(species)


trait_cols <- setdiff(names(traits_for_imputation), "species")

traits_log <- traits_for_imputation
traits_log[trait_cols] <- log(traits_log[trait_cols])

traits_scaled <- traits_log

scaled_matrix <- scale(traits_log[trait_cols])

trait_means <- attr(scaled_matrix, "scaled:center")
trait_sds <- attr(scaled_matrix, "scaled:scale")


traits_scaled[trait_cols] <- scaled_matrix


# ---- 4. RUN PHYLOPARS  ----
set.seed(42)

phylo_imputation <- phylopars(
  model = "BM",
  trait_data = traits_scaled,
  tree = tree,
  pheno_error = FALSE,
  phylo_correlated = TRUE,
  npd = TRUE,
  usezscores = FALSE,
  full_alpha = FALSE,
)

# ---- 5. UNSCALE DATA  ----

imputed_data <- as.data.frame(phylo_imputation$anc_recon)%>%
  rownames_to_column( var = "species")

imputed_data <- subset(imputed_data, species %in% traits_scaled$species)

unscaled_matrix <- sweep(imputed_data[trait_cols], 2, trait_sds, FUN = "*")
unscaled_matrix <- sweep(unscaled_matrix, 2, trait_means, FUN = "+")

imputed_final <- imputed_data


imputed_final[trait_cols] <- exp(unscaled_matrix)

imputed_final <- imputed_final %>% 
  mutate(species = str_replace_all(species, "_", " "))


# ---- 6. LONG TRASFORM FOR QUALITY CHECK  ----

long_raw <- traits_for_imputation %>%
  select(species, all_of(trait_cols)) %>%
  pivot_longer(cols = -species, names_to = "trait", values_to = "value") %>%
  mutate(source = "raw") %>%
  filter(!is.na(value))

long_imp <- imputed_final %>%
  select(species, all_of(trait_cols)) %>%
  pivot_longer(cols = -species, names_to = "trait", values_to = "value") %>%
  mutate(source = "Imputed") %>%
  filter(!is.na(value))

long_both <- bind_rows(long_raw, long_imp)

means <- long_both %>%
  group_by(trait, source) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    mean_log_val = mean(log(value), na.rm = TRUE), 
    .groups = "drop"
  )

# ---- 7. PLOT TRAITS  ----

phylo_dist <- ggplot(long_both, aes(x = value, fill = source, colour = source)) +
  geom_density(alpha = 0.5, linewidth = 0.8) +
  geom_vline(data = means, aes(xintercept = mean_val, colour = source),
             linetype = "dashed", linewidth = 1) +
  facet_wrap(~ trait, scales = "free") +
  scale_fill_manual(values = c("raw" = "lightblue", "Imputed" = "orange")) +
  scale_colour_manual(values = c("raw" = "blue", "Imputed" = "darkorange")) +
  theme_minimal() +
  labs(
    title = "Trait data with NAs vs Imputed trait densities - phylopars",
    subtitle = "Dashed lines = means (Original Scale)",
    x = "Trait value", 
    y = "Density",
    fill = "Data", 
    colour = "Data"
  ) +
  theme(legend.position = "bottom")

phylo_dist_log <- ggplot(long_both, aes(x = log(value), fill = source, colour = source)) +
  geom_density(alpha = 0.5, linewidth = 0.8) +
  geom_vline(data = means, aes(xintercept = mean_log_val, colour = source),
             linetype = "dashed", linewidth = 1) +
  facet_wrap(~ trait, scales = "free") +
  scale_fill_manual(values = c("raw" = "lightblue", "Imputed" = "orange")) +
  scale_colour_manual(values = c("raw" = "blue", "Imputed" = "darkorange")) +
  theme_minimal() +
  labs(
    title = "Trait data with NAs vs Imputed trait densities - Phylopars",
    subtitle = "Dashed lines = means (Log Transformed)",
    x = "Trait log(value)", 
    y = "Density",
    fill = "Data", 
    colour = "Data"
  ) +
  theme(legend.position = "bottom")

# ---- 7.5 COVARIANCE / CORRELATION CHECK (PRE vs POST IMPUTATION) ----

raw_complete <- traits_for_imputation %>%
  select(all_of(trait_cols)) %>%
  drop_na()

imp_values <- imputed_final %>%
  select(all_of(trait_cols))


cor_raw <- cor(raw_complete,  use = "pairwise.complete.obs", method = "pearson")
cor_imp <- cor(imp_values,    use = "complete.obs",          method = "pearson")
cor_diff <- cor_imp - cor_raw

melt_raw  <- melt(cor_raw)  %>% mutate(matrix = "Raw (complete cases)")
melt_imp  <- melt(cor_imp)  %>% mutate(matrix = "Imputed")
melt_diff <- melt(cor_diff) %>% mutate(matrix = "Difference (Imputed - Raw)")

melt_all  <- bind_rows(melt_raw, melt_imp, melt_diff)

phylo_cov <- ggplot(melt_all, aes(Var1, Var2, fill = value)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(value, 2)), size = 3) +
  facet_wrap(~ matrix) +
  scale_fill_gradient2(
    low  = "#2166AC", mid = "white", high = "#D6604D",
    midpoint = 0, limits = c(-1, 1), name = "r"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    title    = "Trait correlation structure: Raw vs Imputed (Phylopars)",
    x = NULL, y = NULL
  )


phylopars_plot <- list(phylo_dist, phylo_dist_log, phylo_cov)
saveRDS(phylopars_plot, file="1_phylopars/phylopars_validation_plots.RDS")

# ---- 8. GAP FILL NATIONAL CHECKLISTS  ----

selected_countries_list <- readRDS("selected_countries_for_imputation.rds")


gap_filled_countries <- map(selected_countries_list, function(.x){
 .x <- .x %>%
   select(- Max_plant_height_m, -SLA_mm2_mg, -LDMC_g_g, -Seed_mass_mg)
 
 .x <- .x %>%
   left_join(imputed_final%>%select(species, Max_plant_height_m, Seed_mass_mg, SLA_mm2_mg, LDMC_g_g), by = c("species_name"= "species"))
})


saveRDS(gap_filled_countries, file = "1_phylopars/gap_filled_selected_countries_phylopars.rds")


# ==============================================================================
# ---- 9. CROSS VALIDATION  ----
# ==============================================================================


cv_iterations <- 100   
mask_fraction <- 0.01

calculate_nrmse <- function(observed, predicted) {
  rmse <- sqrt(mean((observed - predicted)^2, na.rm = TRUE))
  sd_obs <- sd(observed, na.rm = TRUE)
  return(rmse / sd_obs)
}

# ------------------------------------------------------------------------------
# ---- 9A. LEAVE 1% OF OBSERVATIONS OUT  ----
# ------------------------------------------------------------------------------
cat("Leave 1% of OBSERVATIONS out...")
cv_obs_results <- list()

for (i in 1:cv_iterations) {
  cat(sprintf("iteration N %d/%d\n", i, cv_iterations))
  
  masked_data_obs <- traits_scaled
  masked_indices <- list()
  
  # mask
  for (trait in trait_cols) {
    observed_rows <- which(!is.na(masked_data_obs[[trait]]))
    n_mask <- max(1, round(length(observed_rows) * mask_fraction))
    to_mask <- sample(observed_rows, n_mask)
    
    masked_indices[[trait]] <- to_mask
    masked_data_obs[to_mask, trait] <- NA
  }
  
  # gap filling
  p_imp <- phylopars(
    model = "BM", trait_data = masked_data_obs, tree = tree, 
    pheno_error = FALSE, phylo_correlated = TRUE, npd = TRUE, 
    usezscores = FALSE, full_alpha = FALSE
  )
  
  # results
  imp_df <- as.data.frame(p_imp$anc_recon)
  
  for (trait in trait_cols) {
    idx <- masked_indices[[trait]]
    sp_names <- masked_data_obs$species[idx] # masked names
    
    true_vals <- traits_scaled[idx, trait]
    pred_vals <- imp_df[sp_names, trait]    
    
    cv_obs_results[[length(cv_obs_results) + 1]] <- data.frame(
      iteration = i, trait = trait, true_value = true_vals, imputed_value = pred_vals
    )
  }
}
df_obs_results <- bind_rows(cv_obs_results)


# ------------------------------------------------------------------------------
# ---- 9B. LEAVE 1% OF SPECIES OUT  ----
# ------------------------------------------------------------------------------
cat("Leave 1% of SPECIES out...")
cv_sp_results <- list()
all_species <- traits_scaled$species

for (i in 1:cv_iterations) {
  cat(sprintf("Iteration N %d/%d\n", i, cv_iterations))
  
  masked_data_sp <- traits_scaled
  n_sp_mask <- max(1, round(length(all_species) * mask_fraction))
  sp_to_mask <- sample(all_species, n_sp_mask)
  
  # index numbers
  idx_mask <- which(masked_data_sp$species %in% sp_to_mask)
  
  # save true values
  true_vals_list <- list()
  for (trait in trait_cols) {
    true_vals_list[[trait]] <- traits_scaled[idx_mask, trait]
  }
  
  # delete all traits
  masked_data_sp[idx_mask, trait_cols] <- NA
  
  # gap filling
  p_imp_sp <- phylopars(
    model = "BM", trait_data = masked_data_sp, tree = tree, 
    pheno_error = FALSE, phylo_correlated = TRUE, npd = TRUE, 
    usezscores = FALSE, full_alpha = FALSE
  )
  
  imp_df_sp <- as.data.frame(p_imp_sp$anc_recon)
  
  # results
  for (trait in trait_cols) {
    true_vals <- true_vals_list[[trait]]
    
    # evaluate error
    valid_idx <- !is.na(true_vals) 
    
    if (any(valid_idx)) {
      sp_valid <- sp_to_mask[valid_idx]
      pred_vals <- imp_df_sp[sp_valid, trait]
      
      cv_sp_results[[length(cv_sp_results) + 1]] <- data.frame(
        iteration = i, trait = trait, species = sp_valid,
        true_value = true_vals[valid_idx], imputed_value = pred_vals
      )
    }
  }
}
df_sp_results <- bind_rows(cv_sp_results)


# ------------------------------------------------------------------------------
# ---- 9C. SUMMARY----
# ------------------------------------------------------------------------------

# stats for validation at obs level
stats_obs <- df_obs_results %>%
  group_by(trait) %>%
  summarise(
    R_Observation_out = cor(true_value, imputed_value, use = "complete.obs", method = "pearson"),
    NRMSE_Observation_out = calculate_nrmse(true_value, imputed_value),
    .groups = "drop"
  )

# stats for validation at species level
stats_sp <- df_sp_results %>%
  group_by(trait) %>%
  summarise(
    R_Species_out = cor(true_value, imputed_value, use = "complete.obs", method = "pearson"),
    NRMSE_Species_out = calculate_nrmse(true_value, imputed_value),
    .groups = "drop"
  )

#merge results
final_validation_stats_phylo <- left_join(stats_obs, stats_sp, by = "trait")

cat("=======================================================")
cat("CROSS VALIDATION RESULTS (PHYLOPARS)")
cat("=======================================================")
print(final_validation_stats_phylo %>% gt() %>% fmt_number(columns = -trait, decimals = 3))

write_csv(final_validation_stats_phylo, "1_phylopars/phylopars_validation_statistics.csv")



