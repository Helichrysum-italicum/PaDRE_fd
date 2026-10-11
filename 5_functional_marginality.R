# ==============================================================================
# PROJECT:     PaDRE_FD (Functional Diversity)
# AUTHOR:      Samuele Nicosia
# DATE:        2026-04-20
# DESCRIPTION: 
#use BAT library to compute hypervolume and origiality metrics at country level
#export violin plots with results.
# ==============================================================================

# ---- 0. SETUP   ----
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggExtra)
library(ggside)
library(hypervolume)
library(agricolae)
library(patchwork) 
library(ggnewscale)
library(stringr)
library(BAT)
library(purrr)
library(ggtext)  
library(tidyverse)
library(tools)

setwd("xxx")
set.seed(42)

gap_filled_countries <- readRDS("1_missforest/gap_filled_selected_countries_mf.rds")

gap_filled_countries <- map(gap_filled_countries, function(df) {
  filter(df, SLA_mm2_mg >= 1 & LDMC_g_g > 0.04)
})

# ---- 1. COMPUTE METRICS   ----

master_metrics_list <- imap(gap_filled_countries, function(.x, .y) {
  
  message("working on: ", .y)
  
  # A. prepare and scale data
  scaled_species <- .x %>%
    mutate(species_name = make.unique(as.character(species_name))) %>%
    select(species_name, status, Growth_form_simple, Max_plant_height_m, SLA_mm2_mg, LDMC_g_g, Seed_mass_mg) %>%
    mutate(
      Seed_mass_resid = residuals(lm(log10(Seed_mass_mg) ~ log10(Max_plant_height_m), data = .))
    ) %>%
    mutate(
      Height      = as.numeric(scale(log10(Max_plant_height_m))),
      SLA         = as.numeric(scale(log10(SLA_mm2_mg))),
      LDMC        = as.numeric(scale(log10(LDMC_g_g))),
      `seed mass` = as.numeric(scale(log10(Seed_mass_mg))),
      `seed mass residual` = as.numeric(scale(Seed_mass_resid))     
    ) %>%
    mutate(
      status = case_when(
        status == "extirpated" ~ "extirpated",
        TRUE ~ "extant"
      )
    )
  
  # B. PCA
  PCA <- prcomp(scaled_species[, c("Height","SLA", "LDMC" 
                                   ,"seed mass"
                                   #,"seed mass residual"
  )])
  PCAvalues <- data.frame(Species = scaled_species$species_name, Status = scaled_species$status, PCA$x)
  PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)
  
  # C. Hypervolume & BAT Metrics
  pc_cols <- c("PC1", "PC2", "PC3")
  overall_bandwidth <- estimate_bandwidth(PCAvalues[, pc_cols])
  
  traits <- as.matrix(PCAvalues[, pc_cols])
  mode(traits) <- "numeric"
  rownames(traits) <- PCAvalues$Species
  
  comm_matrix_single <- matrix(rep(1, nrow(PCAvalues)), nrow = 1, ncol = nrow(PCAvalues), byrow = TRUE)
  rownames(comm_matrix_single) <- "country_flora"
  colnames(comm_matrix_single) <- PCAvalues$Species
  
  #depending on your versions of the BAT package and R you might encounter 
  #an error related to the dimensionality of the community matrix
  #add the empty community to fix, BAT will automatically drop the empty community
  
  # comm_matrix_single <- rbind(comm_matrix_single,
  #                        matrix(0, nrow = 1, ncol = ncol(comm_matrix_single),
  #                               dimnames = list("Dummy", colnames(comm_matrix_single))))
  
  hv_country_bat <- kernel.build(comm = comm_matrix_single, trait = traits, 
                                 method.hv = "gaussian", kde.bandwidth = overall_bandwidth)
  
  sp_originality <- as.data.frame(kernel.originality(hv_country_bat, frac = 0.1)) %>% rename(sp_originality = 1)  
  sp_contribution <- as.data.frame(kernel.contribution(hv_country_bat)) %>% rename(sp_contribution = 1)
  
  # Combine metrics
  combined_metrics <- cbind(sp_originality, sp_contribution) %>%
    rownames_to_column(var = "Species") %>%
    left_join(PCAvalues %>% select(Species, Status) %>% distinct(), by = "Species") %>%
    mutate(Status = as.factor(Status))
  
  return(list(
    PCA_Data = list(Values = PCAvalues, Loadings = PCAloadings),
    Metrics = combined_metrics
  ))
})

saveRDS(master_metrics_list, file = "originality_metrics.RDS")

master_metrics_list <- readRDS("originality_metrics.RDS")
# ---- 2. CREATE LISTS   ----


#  METRICS LIST
metrics_list <- map(master_metrics_list, function(.x) {
  .x$Metrics %>% select(Species, Status, sp_originality, sp_contribution)
})

#  KRUSKAL WALLIS TEST

kw_tests_list <- map(metrics_list, function(df) {

  
  kw_orig <- kruskal(df$sp_originality, df$Status, console = FALSE)
  kw_contr <- kruskal(df$sp_contribution, df$Status, console = FALSE)
    list(
    Originality = kw_orig,
    Contribution = kw_contr
  )
})


# ==============================================================================
# ---- 3. COMPOSE PLOTS FOR SELECTED COUNTRIES ----
# ==============================================================================



# 1. COMBINE DATA & CLEANUP
all_data <- imap_dfr(metrics_list, ~ mutate(.x, Country = .y)) %>%
  mutate(
    #rename groups
    Plot_Group = if_else(Status == "extirpated", "Extirpated", "Extant"),
    #count extirpated
    n_ext_country = sum(Plot_Group == "Extirpated", na.rm = TRUE),
    .by = Country
  ) %>%
#order by extirpated
    mutate(Country = fct_reorder(Country, n_ext_country)) 


# 2. selection for smaller visualization
selected_countries <- c("Malta", "Croatia", "Austria", "Norway", "Portugal", "Germany", "Czech Republic", "Poland")

data_selected <- all_data %>% filter(Country %in% selected_countries)


# 3.  PLOTTING FUNCTION
plot_violin <- function(df, y_var, test_name, title) {
  
  current_countries <- levels(droplevels(df$Country))
  #label by p-value
  labels <- setNames(map_chr(current_countries, function(cn) {
    p <- kw_tests_list[[cn]][[test_name]]$statistics$p.chisq
    
    p_stars <- case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
    n_ext <- sum(df$Country == cn & df$Plot_Group == "Extirpated", na.rm = TRUE)
    
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
      name = NULL
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

# 4. PLOT GENERATION AND ASSEMBLY
layout_theme <- theme(legend.position = "bottom", plot.margin = margin(8, 16, 8, 8))

# Plot A all countries
full_plot <- (plot_violin(all_data, "sp_originality", "Originality", "Originality across Countries") + 
                plot_violin(all_data, "sp_contribution", "Contribution", "Contribution across Countries")) +
  plot_layout(guides = "collect")+ plot_annotation(title = "Marginality Metrics", 
                                                 theme = theme(plot.title = element_text(size = 20, hjust = 0.5))) & 
  layout_theme 

ggsave(
  plot = full_plot, 
  filename = "2_result_charts/marginality_full.pdf", 
  width = 13, height = 18, units = "in"
)

#plot B selection of countrs
selected_plot <- (plot_violin(data_selected, "sp_originality", "Originality", "Originality across Countries") + 
                    plot_violin(data_selected, "sp_contribution", "Contribution", "Contribution across Countries")) +
  plot_layout(guides = "collect") + plot_annotation(title = "Marginality Metrics", 
                                                  theme = theme(plot.title = element_text(size = 20, hjust = 0.5))) & 
  layout_theme

ggsave(
  plot = selected_plot,
  filename = "2_result_charts/marginality_selected.pdf",
  width = 10, height = 8, units = "in"
)


