# ==============================================================================
# PROJECT:     PaDRE_FD (Functional Diversity)
# AUTHOR:      Samuele Nicosia
# DESCRIPTION:  Null Model for Extirpated vs Random Species
#The null model tests whether extirpated species occupy a significantly different 
#functional space than expected by chance, given a random draw of the same 
# size from the full flora.
# ==============================================================================

# ---- 0. SETUP ----
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)
library(hypervolume)
library(agricolae)
library(BAT)
library(purrr)
library(forcats)
library(ggtext)
library(tools)

setwd("xxx")

gap_filled_countries <- readRDS("1_missforest/gap_filled_selected_countries_mf.rds")

gap_filled_countries <- map(gap_filled_countries, function(df) {
  filter(df, SLA_mm2_mg >= 1 & LDMC_g_g > 0.04)
})

# ---- 1. PCA & TRAIT SPACE PREPARATION ----

# Drop countries with fewer than 10 extirpated species (insufficient for hypervolume)
gap_filled_countries <- discard(gap_filled_countries, function(.x, country_name) {
  sum(.x$status == "extirpated") < 10
})

pca_list <- imap(gap_filled_countries, function(.x, country_name) {
  message("Generating PCA for: ", country_name)
  
  scaled_species <- .x %>%
    mutate(species_name = make.unique(as.character(species_name))) %>%
    select(species_name, status, Max_plant_height_m, SLA_mm2_mg, LDMC_g_g, Seed_mass_mg) %>%
    mutate(
      Seed_mass_resid = residuals(lm(log10(Seed_mass_mg) ~ log10(Max_plant_height_m), data = .))
    ) %>%
    mutate(
      Height             = as.numeric(scale(log10(Max_plant_height_m))),
      SLA                = as.numeric(scale(log10(SLA_mm2_mg))),
      LDMC               = as.numeric(scale(log10(LDMC_g_g))),
      `seed mass`        = as.numeric(scale(log10(Seed_mass_mg))),
      `seed mass resid`  = as.numeric(scale(Seed_mass_resid))
    ) %>%
    mutate(status = if_else(status == "extirpated", "extirpated", "extant"))
  
  PCA <- prcomp(scaled_species[, c("Height", "SLA", "LDMC", "seed mass")])
  
  PCAvalues <- data.frame(
    Species = scaled_species$species_name,
    Status  = scaled_species$status,
    PCA$x
  )
  
  pc_cols <- c("PC1", "PC2", "PC3")
  traits_mat <- as.matrix(PCAvalues[, pc_cols])
  rownames(traits_mat) <- PCAvalues$Species
  
  overall_bandwidth <- estimate_bandwidth(PCAvalues[, pc_cols])
  
  list(
    PCAvalues         = PCAvalues,
    Traits            = traits_mat,
    overall_bandwidth = overall_bandwidth
  )
})


# ==============================================================================
# ---- 2. COMPUTE OBSERVED METRICS ----
#
# Observed metrics are computed ONCE per country on ALL extirpated species.
#   SES = (X_obs - mean(X_null)) / sd(X_null)
# ==============================================================================

observed_metrics <- imap(pca_list, function(pca_data, country_name) {
  message("Computing observed metrics for: ", country_name)
  
  PCAvalues        <- pca_data$PCAvalues
  overall_bandwidth <- pca_data$overall_bandwidth
  extirpated_sp    <- PCAvalues$Species[PCAvalues$Status == "extirpated"]
  
  # Build trait matrix for ALL extirpated species (observed community)
  traits_obs <- PCAvalues %>%
    filter(Species %in% extirpated_sp) %>%
    select(Species, PC1, PC2, PC3) %>%
    column_to_rownames("Species") %>%
    as.matrix()
  
  comm_obs <- matrix(
    1, nrow = 1, ncol = nrow(traits_obs),
    dimnames = list("Observed", rownames(traits_obs))
  )
  
  hv_obs <- BAT::kernel.build(
    comm = comm_obs, trait = traits_obs,
    method.hv = "gaussian", kde.bandwidth = overall_bandwidth,
    abund = FALSE, cores = 0
  )
  
  list(
    obs_alpha    = as.numeric(kernel.alpha(hv_obs)),
    obs_evenness = as.numeric(kernel.evenness(hv_obs)),
    n_extirpated = length(extirpated_sp)
  )
})


# ==============================================================================
# ---- 3. GENERATE NULL DISTRIBUTIONS ----
# ==============================================================================

nperm <- 499  # number of null permutations

null_distributions <- imap(pca_list, function(pca_data, country_name) {
  message("Generating null distribution for: ", country_name)
  
  PCAvalues         <- pca_data$PCAvalues
  overall_bandwidth  <- pca_data$overall_bandwidth
  n_extirpated      <- observed_metrics[[country_name]]$n_extirpated
  
  map_dfr(seq_len(nperm), function(i) {
    
    # Random draw of n_extirpated species from the entire flora
    samp_null <- sample(PCAvalues$Species, n_extirpated, replace = FALSE)
    
    traits_null <- PCAvalues %>%
      filter(Species %in% samp_null) %>%
      select(Species, PC1, PC2, PC3) %>%
      column_to_rownames("Species") %>%
      as.matrix()
    
    comm_null <- matrix(
      1, nrow = 1, ncol = nrow(traits_null),
      dimnames = list("Null", rownames(traits_null))
    )
    
    hv_null <- BAT::kernel.build(
      comm = comm_null, trait = traits_null,
      method.hv = "gaussian", kde.bandwidth = overall_bandwidth,
      abund = FALSE, cores = 0
    )
    
    tibble(
      Iteration     = i,
      Null_Richness = as.numeric(kernel.alpha(hv_null)),
      Null_Evenness = as.numeric(kernel.evenness(hv_null))
    )
    
  }, .progress = TRUE)
})

saveRDS(null_distributions, file = "null_data.RDS")
 null_distributions <- readRDS("null_data.RDS")


# ==============================================================================
# ---- 4. SIGNIFICANCE TESTING ----
#
# SES = (X_obs - mean(X_null)) / sd(X_null)
# ==============================================================================

compute_ses_stats <- function(obs_value, null_vector) {
  mean_null <- mean(null_vector)
  sd_null   <- sd(null_vector)
  
  ses <- (obs_value - mean_null) / sd_null
  
  
  p_low  <- (sum(null_vector <= obs_value) + 1) / (length(null_vector) + 1)
  p_high <- (sum(null_vector >= obs_value) + 1) / (length(null_vector) + 1)
  
  p_value <- 2 * min(p_low, p_high)
  
  p_value <- min(p_value, 1)
  list(ses = ses, p_value = p_value, mean_null = mean_null, sd_null = sd_null)
}

# Apply the test to every country and metric
ses_results <- imap_dfr(null_distributions, function(null_dist, country_name) {
  
  obs  <- observed_metrics[[country_name]]
  
  rich <- compute_ses_stats(obs$obs_alpha,    null_dist$Null_Richness)
  even <- compute_ses_stats(obs$obs_evenness, null_dist$Null_Evenness)
  
  tibble(
    Country       = country_name,
    N_Extirpated  = obs$n_extirpated,
    
    Obs_Richness  = obs$obs_alpha,
    Mean_Null_Richness = rich$mean_null,
    SES_Richness  = rich$ses,
    Pval_Richness = rich$p_value,
    
    Obs_Evenness  = obs$obs_evenness,
    Mean_Null_Evenness = even$mean_null,
    SES_Evenness  = even$ses,
    Pval_Evenness = even$p_value
  )
}) %>%
  arrange(desc(N_Extirpated)) %>%
  mutate(Country = fct_reorder(Country, N_Extirpated))

write.csv(ses_results, "2_result_charts/SES_Summary_Statistics.csv", row.names = FALSE)


# ==============================================================================
# ---- 5. PREPARE DATA FOR PLOTTING ----
# ==============================================================================

# Long-format null distributions, factor levels matched to ses_results ordering
null_data_long <- imap_dfr(null_distributions, function(null_dist, country_name) {
  null_dist %>% mutate(Country = country_name)
}) %>%
  mutate(Country = factor(Country, levels = levels(ses_results$Country)))


# ==============================================================================
# ---- 6. PLOT FUNCTION: OBSERVED VALUE vs NULL DISTRIBUTION ----
# ==============================================================================

plot_null_vs_obs <- function(null_df, obs_df, metric_name, title) {
  
  null_col <- paste0("Null_", metric_name)
  obs_col  <- paste0("Obs_",  metric_name)
  pval_col <- paste0("Pval_", metric_name)
  
  current_countries <- levels(droplevels(obs_df$Country))
  
  # Build y-axis labels with significance annotation
  labels <- setNames(map_chr(current_countries, function(cn) {
    row   <- obs_df %>% filter(Country == cn)
    p     <- row[[pval_col]]
    n_ext <- row$N_Extirpated
    
    # Two-tailed stars: significant if p < 0.025 (either tail at alpha = 0.05)
    p_stars <- case_when(
      p <= 0.001 ~ "***",
      p <= 0.010 ~ "**",
      p <= 0.050 ~ "*",
      TRUE ~ ""
    )
    
    if (p_stars != "") {
      paste0(cn, " (*n* = ", n_ext, ")<sup>", p_stars, "</sup>")
    } else {
      paste0(cn, " (*n* = ", n_ext, ")")
    }
  }), current_countries)
  
  ggplot() +
    # Null distribution (grey violin)
    geom_violin(
      data = null_df,
      aes(x = Country, y = .data[[null_col]], fill = "Null distribution"),
      color = "grey50", scale = "width", alpha = 0.7
    ) +
    # Observed value for extirpated species (purple diamond)
    geom_point(
      data = obs_df,
      aes(x = Country, y = .data[[obs_col]], color = "Observed (extirpated)"),
      size = 3.5, shape = 18
    ) +
    scale_fill_manual(
      name   = NULL,
      values = c("Null distribution" = "grey80")
    ) +
    scale_color_manual(
      name   = NULL,
      values = c("Observed (extirpated)" = "#9B4F96")
    ) +
    guides(
      fill  = guide_legend(override.aes = list(alpha = 0.7, size = 5, shape = 22, color = "grey50")),
      color = guide_legend(override.aes = list(size = 4, shape = 18))
    ) +
    scale_x_discrete(labels = labels) +
    coord_flip() +
    labs(title = title, x = NULL, y = paste("Functional", metric_name)) +
    theme_bw(base_size = 9) +
    theme(
      plot.title         = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.y        = element_markdown(size = 9, lineheight = 1.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      legend.position    = "bottom",
      legend.text        = element_text(size = 9)
    )
}


# ==============================================================================
# ---- 7. GENERATE AND SAVE PLOTS ----
# ==============================================================================

layout_theme <- theme(plot.margin = margin(8, 16, 8, 8))

# --- A. All countries ---
full_plot <- (
  plot_null_vs_obs(null_data_long, ses_results, "Richness", "Functional Richness\n(Observed vs Null)") +
    plot_null_vs_obs(null_data_long, ses_results, "Evenness", "Functional Evenness\n(Observed vs Null)")
) + 
  plot_layout(guides = "collect") & 
  layout_theme & 
  theme(legend.position = "bottom")

ggsave(
  plot     = full_plot,
  filename = "2_result_charts/FD_SES_null_model_full.pdf",
  width = 12, height = 14, units = "in"
)

# --- B. Selected countries only ---
selected_countries <- c("Malta", "Croatia", "Austria", "Norway",
                        "Portugal", "Germany", "Czech Republic", "Poland")

null_sel  <- null_data_long %>% filter(Country %in% selected_countries) %>% mutate(Country = droplevels(Country))
stats_sel <- ses_results    %>% filter(Country %in% selected_countries) %>% mutate(Country = droplevels(Country))

selected_plot <- (
  plot_null_vs_obs(null_sel, stats_sel, "Richness", "Functional Richness\n(Observed vs Null)") +
    plot_null_vs_obs(null_sel, stats_sel, "Evenness", "Functional Evenness\n(Observed vs Null)")
) + 
  plot_layout(guides = "collect") & 
  layout_theme & 
  theme(legend.position = "bottom")

ggsave(
  plot     = selected_plot,
  filename = "2_result_charts/FD_SES_null_model_selected.pdf",
  width = 10, height = 6, units = "in"
)