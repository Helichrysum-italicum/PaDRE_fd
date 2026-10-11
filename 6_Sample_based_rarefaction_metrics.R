# ==============================================================================
# PROJECT:     PaDRE_FD (Functional Diversity)
# AUTHOR:      Samuele Nicosia
# DESCRIPTION: Rarefaction-based null model — Extirpated vs Extant comparison
#At each iteration, min_rar species are randomly drawn from the
#extirpated pool AND from the extant pool separately. FD metrics
#are computed for both groups, building two parallel distributions
#of 499 values each. The test then compares whether the two
#distributions differ systematically.
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
library(tidyverse)

setwd("xxx")

gap_filled_countries <- readRDS("1_missforest/gap_filled_selected_countries_mf.rds")

gap_filled_countries <- map(gap_filled_countries, function(df) {
  filter(df, SLA_mm2_mg >= 1 & LDMC_g_g > 0.04)
})


# ==============================================================================
# ---- 1. PCA & TRAIT SPACE PREPARATION ----
# ==============================================================================

# Drop countries with fewer than 10 extirpated species
gap_filled_countries <- discard(gap_filled_countries, function(.x, country_name) {
  sum(.x$status == "extirpated") < 10
})

pca_list <- imap(gap_filled_countries, function(.x, country_name) {
  message("Generating PCA for: ", country_name)
  
  scaled_species <- .x %>%
    mutate(species_name = make.unique(as.character(species_name))) %>%
    select(species_name, status, Growth_form_simple,
           Max_plant_height_m, SLA_mm2_mg, LDMC_g_g, Seed_mass_mg) %>%
    mutate(
      Seed_mass_resid = residuals(lm(log10(Seed_mass_mg) ~ log10(Max_plant_height_m), data = .))
    ) %>%
    mutate(
      Height             = as.numeric(scale(log10(Max_plant_height_m))),
      SLA                = as.numeric(scale(log10(SLA_mm2_mg))),
      LDMC               = as.numeric(scale(log10(LDMC_g_g))),
      `seed mass`        = as.numeric(scale(log10(Seed_mass_mg))),
      `seed mass residual` = as.numeric(scale(Seed_mass_resid))
    ) %>%
    mutate(status = if_else(status == "extirpated", "extirpated", "extant"))
  
  PCA <- prcomp(scaled_species[, c("Height", "SLA", "LDMC", "seed mass")])
  
  PCAvalues <- data.frame(
    Species = scaled_species$species_name,
    Status  = scaled_species$status,
    PCA$x
  )
  
  pc_cols   <- c("PC1", "PC2", "PC3")
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
# ---- 2. RAREFACTION NULL MODEL — METRIC COMPUTATION ----
# ==============================================================================

nperm   <- 499
min_rar <- 6

null_model_list <- imap(pca_list, function(pca_data, country_name) {
  message("Running rarefaction null model for: ", country_name)
  
  PCAvalues         <- pca_data$PCAvalues
  overall_bandwidth <- pca_data$overall_bandwidth
  
  extirpated_sp <- PCAvalues$Species[PCAvalues$Status == "extirpated"]
  extant_sp     <- PCAvalues$Species[PCAvalues$Status == "extant"]
  
  map_dfr(seq_len(nperm), function(i) {
    
    # Random draw of min_rar species from each pool independently
    samp_ext  <- sample(extirpated_sp, min_rar, replace = FALSE)
    samp_rest <- sample(extant_sp,     min_rar, replace = FALSE)
    
    # Build trait matrices
    traits_extirpated <- PCAvalues %>%
      filter(Species %in% samp_ext) %>%
      select(Species, PC1, PC2, PC3) %>%
      column_to_rownames("Species") %>%
      as.matrix()
    
    traits_extant <- PCAvalues %>%
      filter(Species %in% samp_rest) %>%
      select(Species, PC1, PC2, PC3) %>%
      column_to_rownames("Species") %>%
      as.matrix()
    
    # Build community matrices (presence-only, single row)
    comm_extirpated <- matrix(
      1, nrow = 1, ncol = nrow(traits_extirpated),
      dimnames = list("Extirpated_Community", rownames(traits_extirpated))
    )
    comm_extant <- matrix(
      1, nrow = 1, ncol = nrow(traits_extant),
      dimnames = list("Extant_Community", rownames(traits_extant))
    )
    
    # Build hypervolumes
    hv_extirpated <- BAT::kernel.build(
      comm = comm_extirpated, trait = traits_extirpated,
      method.hv = "gaussian", kde.bandwidth = overall_bandwidth,
      abund = FALSE, cores = 0
    )
    hv_extant <- BAT::kernel.build(
      comm = comm_extant, trait = traits_extant,
      method.hv = "gaussian", kde.bandwidth = overall_bandwidth,
      abund = FALSE, cores = 0
    )
    
    tibble(
      Iteration              = i,
      FD_Richness_Extirpated = as.numeric(kernel.alpha(hv_extirpated)),
      FD_Evenness_Extirpated = as.numeric(kernel.evenness(hv_extirpated)),
      FD_Richness_Extant     = as.numeric(kernel.alpha(hv_extant)),
      FD_Evenness_Extant     = as.numeric(kernel.evenness(hv_extant))
    )
    
  }, .progress = TRUE)
})

saveRDS(null_model_list, file = "rarefaction_data.RDS")
null_model_list <- readRDS("rarefaction_data.RDS")


n_extirpated_country <- map_dbl(pca_list, ~ sum(.x$PCAvalues$Status == "extirpated"))


# ==============================================================================
# ---- 4. SIGNIFICANCE TESTING ----
# ==============================================================================

null_model_list <- imap(null_model_list, function(null_dist, country_name) {
  message("Running significance test for: ", country_name)
  
  kw_test <- function(col_extirpated, col_extant) {
    kruskal.test(
      x = c(null_dist[[col_extirpated]], null_dist[[col_extant]]),
      g = factor(rep(c("Extirpated", "Extant"), each = nrow(null_dist)))
    )
  }
  
  kw_rich <- kw_test("FD_Richness_Extirpated", "FD_Richness_Extant")
  kw_even <- kw_test("FD_Evenness_Extirpated", "FD_Evenness_Extant")
  
  # Attach p-values as constant columns for downstream use in plots
  null_dist %>% mutate(
    kruskal_alpha    = kw_rich$p.value,
    kruskal_evenness = kw_even$p.value
  )
})


# ==============================================================================
# ---- 5. DATA PREPARATION FOR PLOTTING ----
# ==============================================================================

all_data <- imap_dfr(null_model_list, function(df, country_name) {
  df %>%
    pivot_longer(
      cols          = c(FD_Richness_Extirpated, FD_Evenness_Extirpated,
                        FD_Richness_Extant,     FD_Evenness_Extant),
      names_to      = c(".value", "Plot_Group"),
      names_pattern = "(FD_Richness|FD_Evenness)_(Extirpated|Extant)"
    ) %>%
    mutate(Country = country_name)
}) %>%
  mutate(
    n_ext_country = n_extirpated_country[Country],
    Country       = fct_reorder(Country, n_ext_country)
  )

selected_countries <- c("Malta", "Croatia", "Austria", "Norway",
                        "Portugal", "Germany", "Czech Republic", "Poland")

data_selected <- all_data %>% filter(Country %in% selected_countries)


# ==============================================================================
# ---- 6. PLOT FUNCTION ----
# ==============================================================================

plot_violin <- function(df, y_var, pval_col, title) {
  
  current_countries <- levels(droplevels(df$Country))
  
  labels <- setNames(map_chr(current_countries, function(cn) {
    p     <- df %>% filter(Country == cn) %>% pull(.data[[pval_col]]) %>% first()
    n_ext <- n_extirpated_country[[cn]]
    
    p_stars <- case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    )
    
    if (p_stars != "") {
      paste0(cn, " (*n* = ", n_ext, ")<sup>", p_stars, "</sup>")
    } else {
      paste0(cn, " (*n* = ", n_ext, ")")
    }
  }), current_countries)
  
  ggplot(df, aes(x = Country, y = .data[[y_var]], fill = Plot_Group)) +
    geom_violin(
      width = 0.8, trim = FALSE, scale = "width", alpha = 0.65,
      position = position_dodge(0.8), adjust = 1.1,
      linejoin = "round", lineend = "round"
    ) +
    geom_boxplot(
      aes(group = interaction(Country, Plot_Group)),
      width = 0.2, fill = "white", outlier.shape = NA,
      position = position_dodge(0.8)
    ) +
    scale_x_discrete(labels = labels) +
    scale_fill_manual(
      values = c("Extant" = "grey65", "Extirpated" = "#9B4F96"),
      breaks = c("Extant", "Extirpated"),
      name   = NULL
    ) +
    coord_flip() +
    labs(title = title, x = NULL, y = toTitleCase(gsub("_", " ", y_var))) +
    theme_bw(base_size = 9) +
    theme(
      plot.title         = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.y        = element_markdown(size = 9, lineheight = 1.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    )
}


# ==============================================================================
# ---- 7. GENERATE AND SAVE PLOTS ----
# ==============================================================================

layout_theme <- theme(legend.position = "bottom", plot.margin = margin(8, 16, 8, 8))

# --- A. All countries ---
full_plot <- (
  plot_violin(all_data, "FD_Richness", "kruskal_alpha",    "Functional Richness across Countries") +
    plot_violin(all_data, "FD_Evenness", "kruskal_evenness", "Functional Evenness across Countries")
) + plot_layout(guides = "collect") & layout_theme

ggsave(
  plot     = full_plot,
  filename = "2_result_charts/FD_rarefaction_violins_full.pdf",
  width = 15, height = 18, units = "in"
)

# --- B. Selected countries only ---
selected_plot <- (
  plot_violin(data_selected, "FD_Richness", "kruskal_alpha",    "Functional Richness across Countries") +
    plot_violin(data_selected, "FD_Evenness", "kruskal_evenness", "Functional Evenness across Countries")
) + plot_layout(guides = "collect") & layout_theme

ggsave(
  plot     = selected_plot,
  filename = "2_result_charts/FD_rarefaction_violins_selected.pdf",
  width = 10, height = 8, units = "in"
)