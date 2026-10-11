# ==============================================================================
# PROJECT:     PaDRE_FD (Functional Diversity)
# AUTHOR:      Samuele Nicosia
# DATE:        2026-05-28
# DESCRIPTION: 
# Prepare taxonomy and build phylogenetic tree using V.PhyloMaker2.
# Extract phylogenetic eigenvectors and optimize missForest imputation error.
# Gap fill missing trait data using missForest and re-subset into national lists.
# ==============================================================================

library(tidyverse)
library(stringr)
library(visdat)
library(randomForest)
library(V.PhyloMaker2)
library(missForest)
library(adespatial)
library(adephylo)
library(ape)
library(gt)
library(reshape2)
library(RSpectra)
library(PVR)
library(doParallel)


num_cores <- min(4, max(1, detectCores() - 1))
cl <- makeCluster(num_cores)
registerDoParallel(cl)
setwd("xxx")

# Load data
sp_list <- readRDS("species_pool_for_imputation.rds")
continent_pool <- sp_list


# ---- 1. PREPARE TAXONOMY  ----

taxonomy <- sp_list %>%
  select(family, genus, species_name) %>%
  rename(species = species_name) %>%
  relocate(species, genus, family) %>%
  mutate(family = case_when(
    genus == "Comandra"    ~ "Comandraceae",
    genus == "Corrigiola"  ~ "Caryophyllaceae",
    genus == "Cyrtomium"   ~ "Dryopteridaceae",
    genus == "Cystopteris" ~ "Cystopteridaceae",
    genus == "Diplazium"   ~ "Athyriaceae",
    genus == "Metaxya"     ~ "Metaxyaceae",
    genus == "Woodsia"     ~ "Woodsiaceae",
    TRUE                   ~ family 
  ))


# ---- 2. BUILD PHYLOGENETIC TREE & EIGENVECTORS ----

# # Build phylogenetic tree using V.PhyloMaker2
result <- phylo.maker(
  sp.list = taxonomy,
  tree = GBOTB.extended.TPL,
  nodes = nodes.info.1.TPL,
  output.sp.list = TRUE,
  output.tree = FALSE,
  scenarios = "S3",
  r = 1
)
#result <- readRDS("phylomakertree.RDS")

phyloeigens <- PVRdecomp(result$scenario.3, scale = T)

# saveRDS(phyloeigens, file="phyloeigens.RDS")
# phyloeigens <- readRDS("phyloeigens.RDS")

# ---- 3. PREPARE TRAITS FOR IMPUTATION PROCESS  ----

trait_names <- c("LDMC_g_g", "Seed_mass_mg", "SLA_mm2_mg", "Max_plant_height_m")  
n_traits <- length(trait_names)

species_ip_ready <- sp_list %>%
  select(species_name, all_of(trait_names)) 

# log1p transformation
species_log <- species_ip_ready %>%
  mutate(across(where(is.numeric), log1p))

# Extract means and standard deviations for back-transformation
trait_means <- sapply(species_log %>% select(where(is.numeric)), mean, na.rm = TRUE)
trait_sds   <- sapply(species_log %>% select(where(is.numeric)), sd, na.rm = TRUE)

# Scaling
species_scaled <- species_log %>%
  mutate(across(where(is.numeric), ~ scale(.)[,1]))


vecnum <- 50

eigenvector <- as.data.frame(phyloeigens@Eigen$vectors[, 1:vecnum]) %>%
  mutate(species_name = gsub("_", " ", phyloeigens@phylo$tip.label))

trait_imp <- left_join(species_scaled, eigenvector, by = "species_name") %>%
  column_to_rownames("species_name")


# ---- 4. ERROR OPTIMISATION LOOP (MISSFOREST) ----

set.seed(42) 

error_df <- data.frame(matrix(NA, nrow = vecnum, ncol = n_traits + 1))
colnames(error_df) <- c("k", paste0("OOB_", trait_names))

for (n in 1:vecnum) {
  df_iter <- trait_imp[, 1:(n_traits + n)]
  
  fit <- missForest(df_iter, maxiter = 5, ntree = 300, variablewise = TRUE, parallelize = "forests")
  
  error_df[n, 1] <- n
  for (t in seq_along(trait_names)) {
    error_df[n, t + 1] <- fit$OOBerror[t]
  }
}

# Visual error table 
error_df %>%
  gt() %>%
  data_color(
    columns = -k,
    fn = function(x) {
      ifelse(x == min(x, na.rm = TRUE), "lightgreen", "white")
    }
  )

# Select the best number of eigenvectors (k) for each trait
best_k <- error_df %>%
  pivot_longer(-k, names_to = "trait", values_to = "error") %>%
  group_by(trait) %>%
  slice_min(error, n = 1) %>%
  slice_min(k, n = 1) %>% 
  ungroup()

k_min_Height <- best_k %>% filter(trait == "OOB_Max_plant_height_m") %>% pull(k)
k_min_LDMC   <- best_k %>% filter(trait == "OOB_LDMC_g_g") %>% pull(k)
k_min_SLA    <- best_k %>% filter(trait == "OOB_SLA_mm2_mg") %>% pull(k)
k_min_Seed   <- best_k %>% filter(trait == "OOB_Seed_mass_mg") %>% pull(k)


# ---- 5. RUN MISSFOREST IMPUTATION ----

set.seed(42) 
Height_ideal <- missForest(trait_imp[, 1:(n_traits + k_min_Height)], maxiter = 20, ntree = 500, variablewise = TRUE, parallelize = "forests") 
LDMC_ideal   <- missForest(trait_imp[, 1:(n_traits + k_min_LDMC)],   maxiter = 20, ntree = 500, variablewise = TRUE, parallelize = "forests") 
SLA_ideal    <- missForest(trait_imp[, 1:(n_traits + k_min_SLA)],    maxiter = 20, ntree = 500, variablewise = TRUE, parallelize = "forests") 
Seed_ideal   <- missForest(trait_imp[, 1:(n_traits + k_min_Seed)],   maxiter = 20, ntree = 500, variablewise = TRUE, parallelize = "forests") 

best_Height  <- tibble(Max_plant_height_m     = Height_ideal$ximp$Max_plant_height_m, species = rownames(Height_ideal$ximp))
best_LDMC    <- tibble(LDMC_g_g     = LDMC_ideal$ximp$LDMC_g_g,       species = rownames(LDMC_ideal$ximp))
best_SLA     <- tibble(SLA_mm2_mg   = SLA_ideal$ximp$SLA_mm2_mg,      species = rownames(SLA_ideal$ximp))
best_Seed    <- tibble(Seed_mass_mg = Seed_ideal$ximp$Seed_mass_mg,   species = rownames(Seed_ideal$ximp))

# Merge imputed traits
impute_out <- best_Height %>%
  left_join(best_LDMC, by = "species") %>%
  left_join(best_SLA,  by = "species") %>%
  left_join(best_Seed, by = "species") %>%
  column_to_rownames(var = "species")


# ---- 6. UNSCALE DATA ----

imputed_final_4traits <- impute_out
current_traits <- colnames(imputed_final_4traits)

for(trait in current_traits) {
  m <- trait_means[[trait]]
  s <- trait_sds[[trait]]
  
  # 1. Reverse standardization
  log_imputed <- (imputed_final_4traits[[trait]] * s) + m
  
  # 2. Reverse log1p (expm1)
  real_value <- as.numeric(expm1(log_imputed))
  
  # 3. Overwrite
  imputed_final_4traits[[trait]] <- real_value
}

imputed_final <- imputed_final_4traits %>%
  rownames_to_column("species_name")

# Patching the original dataset
df_clean <- sp_list %>%
  rows_patch(imputed_final, by = "species_name")


# ---- 7. LONG TRANSFORM FOR QUALITY CHECK  ----

long_raw <- species_ip_ready %>%
  select(species_name, all_of(trait_names)) %>%
  pivot_longer(cols = -species_name, names_to = "trait", values_to = "value") %>%
  mutate(source = "raw") %>%
  filter(!is.na(value))

long_imp <- imputed_final %>%
  select(species_name, all_of(trait_names)) %>%
  pivot_longer(cols = -species_name, names_to = "trait", values_to = "value") %>%
  mutate(source = "Imputed") %>%
  filter(!is.na(value))

long_both <- bind_rows(long_raw, long_imp)

means <- long_both %>%
  group_by(trait, source) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    mean_log_val = mean(log(value), na.rm = TRUE),
    .groups = "drop")
    

# ---- 7A. PLOT TRAITS  ----

trait_distribution <- ggplot(long_both, aes(x = value, fill = source, colour = source)) +
  geom_density(alpha = 0.5, linewidth = 0.8) +
  geom_vline(data = means, aes(xintercept = mean_val, colour = source),
             linetype = "dashed", linewidth = 1) +
  facet_wrap(~ trait, scales = "free") +
  scale_fill_manual(values = c("raw" = "lightblue", "Imputed" = "orange")) +
  scale_colour_manual(values = c("raw" = "blue", "Imputed" = "darkorange")) +
  theme_minimal() +
  labs(
    title = "Trait data with NAs vs Imputed trait densities - missForest",
    subtitle = "Dashed lines = means (Original Scale)",
    x = "Trait value", 
    y = "Density",
    fill = "Data", 
    colour = "Data"
  ) +
  theme(legend.position = "bottom")

trait_distribution_log <- ggplot(long_both, aes(x = log(value), fill = source, colour = source)) +
  geom_density(alpha = 0.5, linewidth = 0.8) +
  geom_vline(data = means, aes(xintercept = mean_log_val, colour = source),
             linetype = "dashed", linewidth = 1) +
  facet_wrap(~ trait, scales = "free") +
   + scale_fill_manual(values = c("raw" = "lightblue", "Imputed" = "orange")) +
  scale_colour_manual(values = c("raw" = "blue", "Imputed" = "darkorange")) +
  theme_minimal() +
  labs(
    title = "Trait data with NAs vs Imputed trait densities - missForest",
    subtitle = "Dashed lines = means (Log Transformed)",
    x = "Trait log(value)", 
    y = "Density",
    fill = "Data", 
    colour = "Data"
  ) +
  theme(legend.position = "bottom")

# ---- 7B. COVARIANCE / CORRELATION CHECK (PRE vs POST IMPUTATION) ----

# Extract raw data (complete cases only)
raw_complete <- species_ip_ready %>%
  select(all_of(trait_names)) %>%
  drop_na()

# Extract final imputed values
imp_values <- imputed_final %>%
  select(all_of(trait_names))

# Calculate correlation matrices
cor_raw <- cor(raw_complete, use = "pairwise.complete.obs", method = "pearson")
cor_imp <- cor(imp_values,   use = "complete.obs",          method = "pearson")
cor_diff <- cor_imp - cor_raw



# Prepare data for ggplot
melt_raw  <- melt(cor_raw)  %>% mutate(matrix = "Raw (complete cases)")
melt_imp  <- melt(cor_imp)  %>% mutate(matrix = "Imputed (missForest)")
melt_diff <- melt(cor_diff) %>% mutate(matrix = "Difference (Imputed - Raw)")

melt_all  <- bind_rows(melt_raw, melt_imp, melt_diff)

p_cov <- ggplot(melt_all, aes(Var1, Var2, fill = value)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(value, 2)), size = 3) +
  facet_wrap(~ matrix) +
  scale_fill_gradient2(
    low  = "#2166AC", mid = "white", high = "#D6604D",
    midpoint = 0, limits = c(-1, 1), name = "r"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    title    = "Trait correlation structure: Raw vs Imputed (missForest)",
    x = NULL, y = NULL
  )

missforest_plot <- list(trait_distribution, trait_distribution_log, p_cov)
saveRDS(missforest_plot, file= "1_missforest/missoforest_validation_plots.RDS")

# ---- 8. GAP FILL NATIONAL CHECKLISTS ----

selected_countries_list <- readRDS("selected_countries_for_imputation.rds")

gap_filled_countries <- map(selected_countries_list, function(.x){
  
  # Remove trait columns potentially already present with NAs
  .x <- .x %>%
    select(-any_of(c("Max_plant_height_m", "SLA_mm2_mg", "LDMC_g_g", "Seed_mass_mg")))
  
  # Merge the new imputed traits from the 'imputed_final' dataframe
  .x <- .x %>%
    left_join(imputed_final %>% select(species_name, all_of(trait_names)), 
              by = "species_name")
})

# Save the gap-filled national checklists in RDS
saveRDS(gap_filled_countries, file = "1_missforest/gap_filled_selected_countries_mf.rds")

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
  
  masked_data_obs <- trait_imp
  masked_indices <- list()
  
  # mask
  for (trait in trait_names) {
    observed_rows <- which(!is.na(masked_data_obs[[trait]]))
    n_mask <- max(1, round(length(observed_rows) * mask_fraction))
    to_mask <- sample(observed_rows, n_mask)
    
    masked_indices[[trait]] <- to_mask
    masked_data_obs[to_mask, trait] <- NA
  }
  
  # gap filling
  Height_cv <- missForest(masked_data_obs[, 1:(n_traits + k_min_Height)], maxiter = 5, ntree = 300, variablewise = TRUE, parallelize = "forests")
  LDMC_cv   <- missForest(masked_data_obs[, 1:(n_traits + k_min_LDMC)],   maxiter = 5, ntree = 300, variablewise = TRUE, parallelize = "forests")
  SLA_cv    <- missForest(masked_data_obs[, 1:(n_traits + k_min_SLA)],    maxiter = 5, ntree = 300, variablewise = TRUE, parallelize = "forests")
  Seed_cv   <- missForest(masked_data_obs[, 1:(n_traits + k_min_Seed)],   maxiter = 5, ntree = 300, variablewise = TRUE, parallelize = "forests")
  
  imp_list <- list(
    Max_plant_height_m     = Height_cv$ximp$Max_plant_height_m,
    LDMC_g_g     = LDMC_cv$ximp$LDMC_g_g,
    SLA_mm2_mg   = SLA_cv$ximp$SLA_mm2_mg,
    Seed_mass_mg = Seed_cv$ximp$Seed_mass_mg
  )
  
  # results
  for (trait in trait_names) {
    idx <- masked_indices[[trait]]
    
    true_vals <- trait_imp[idx, trait]
    pred_vals <- imp_list[[trait]][idx]
    
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
all_species <- rownames(trait_imp)

for (i in 1:cv_iterations) {
  cat(sprintf("iteration N %d/%d\n", i, cv_iterations))
  
  masked_data_sp <- trait_imp
  n_sp_mask <- max(1, round(length(all_species) * mask_fraction))
  sp_to_mask <- sample(all_species, n_sp_mask)
  
  # index numbers
  idx_mask <- which(rownames(masked_data_sp) %in% sp_to_mask)
  
  # save true values
  true_vals_list <- list()
  for (trait in trait_names) {
    true_vals_list[[trait]] <- trait_imp[idx_mask, trait]
  }
  
  # delete all traits
  masked_data_sp[idx_mask, trait_names] <- NA
  
  # gap filling
  Height_cv <- missForest(masked_data_sp[, 1:(n_traits + k_min_Height)], maxiter = 5, ntree = 300, variablewise = TRUE, parallelize = "forests")
  LDMC_cv   <- missForest(masked_data_sp[, 1:(n_traits + k_min_LDMC)],   maxiter = 5, ntree = 300, variablewise = TRUE, parallelize = "forests")
  SLA_cv    <- missForest(masked_data_sp[, 1:(n_traits + k_min_SLA)],    maxiter = 5, ntree = 300, variablewise = TRUE, parallelize = "forests")
  Seed_cv   <- missForest(masked_data_sp[, 1:(n_traits + k_min_Seed)],   maxiter = 5, ntree = 300, variablewise = TRUE, parallelize = "forests")
  
  imp_list <- list(
    Max_plant_height_m     = Height_cv$ximp$Max_plant_height_m,
    LDMC_g_g     = LDMC_cv$ximp$LDMC_g_g,
    SLA_mm2_mg   = SLA_cv$ximp$SLA_mm2_mg,
    Seed_mass_mg = Seed_cv$ximp$Seed_mass_mg
  )
  
  # results
  for (trait in trait_names) {
    true_vals <- true_vals_list[[trait]]
    
    # evaluate error
    valid_idx <- !is.na(true_vals)
    
    if (any(valid_idx)) {
      valid_row_idx <- idx_mask[valid_idx]
      pred_vals <- imp_list[[trait]][valid_row_idx]
      
      cv_sp_results[[length(cv_sp_results) + 1]] <- data.frame(
        iteration = i, trait = trait, species = rownames(trait_imp)[valid_row_idx],
        true_value = true_vals[valid_idx], imputed_value = pred_vals
      )
    }
  }
}
df_sp_results <- bind_rows(cv_sp_results)


# ------------------------------------------------------------------------------
# ---- 9C. SUMMARY ----
# ------------------------------------------------------------------------------

# stats for validation at obs level
stats_obs <- df_obs_results %>%
  group_by(trait) %>%
  summarise(
    R_Observation_out    = cor(true_value, imputed_value, use = "complete.obs", method = "pearson"),
    NRMSE_Observation_out = calculate_nrmse(true_value, imputed_value),
    .groups = "drop"
  )

# stats for validation at species level
stats_sp <- df_sp_results %>%
  group_by(trait) %>%
  summarise(
    R_Species_out    = cor(true_value, imputed_value, use = "complete.obs", method = "pearson"),
    NRMSE_Species_out = calculate_nrmse(true_value, imputed_value),
    .groups = "drop"
  )

# merge results
final_validation_stats_mf <- left_join(stats_obs, stats_sp, by = "trait")

cat("=======================================================")
cat("CROSS VALIDATION RESULTS (MISSFOREST)")
cat("=======================================================")
print(final_validation_stats_mf %>% gt() %>% fmt_number(columns = -trait, decimals = 3))


write_csv(final_validation_stats_mf, "1_missforest/missforest_validation_statistics.csv")

