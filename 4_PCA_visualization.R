# ==============================================================================
# PROJECT:     PaDRE_FD (Functional Diversity)
# AUTHOR:      Samuele Nicosia
# DATE:        2026-04-20
# DESCRIPTION: 
# Compute principal component analysis on log-z transformed data and produce 
# graphic visualization for selected countries
# ==============================================================================

# ---- 0. SETUP   ----
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggside)
library(ggrepel)
library(patchwork)
library(ggnewscale)
library(showtext)
library(purrr) 
library(vegan)
library(purrr)
library(cowplot)
font_add(family = "Georgia", regular = "georgia.ttf")
showtext_auto()

setwd("xxx")


#load gpa filled countries list
gap_filled_countries <- readRDS("1_missforest/gap_filled_selected_countries_mf.rds")

gap_filled_countries <- map(gap_filled_countries, function(df) {
  filter(df, SLA_mm2_mg >= 1 & LDMC_g_g > 0.04)
})

# ---- 1. COMPUTE PCA FOR EACH COUNTRY   ----
pca_data_list <- imap(gap_filled_countries, function(.x, .y) {
  
  # A. prepare and scale data
  scaled_species <- .x %>%
    mutate(species_name = make.unique(as.character(species_name))) %>%
    select(species_name, status, Growth_form_simple, Max_plant_height_m,SLA_mm2_mg, LDMC_g_g,Seed_mass_mg) %>%
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
  
  # B. compute PCA
  PCA <- prcomp(scaled_species[, c("Height","SLA", "LDMC" 
                                   ,"seed mass"
                                   #,"seed mass residual"
                                   )])
  
  if (PCA$rotation["Height", "PC1"] < 0) 
    {PCA$rotation[, "PC1"] <- -PCA$rotation[, "PC1"]
    PCA$x[, "PC1"]        <- -PCA$x[, "PC1"]}
  
  if (PCA$rotation["Height", "PC2"] < 0) 
    {PCA$rotation[, "PC2"] <- -PCA$rotation[, "PC2"]
    PCA$x[, "PC2"]        <- -PCA$x[, "PC2"]}
  
    #explained variance by axis
    pca_var <- PCA$sdev^2
  pca_var_perc <- round(pca_var / sum(pca_var) * 100, 1) 
  
    #extract PCA values
    PCAvalues <- data.frame(
    Species = scaled_species$species_name, 
    status = scaled_species$status, 
    Growth_form_simple = scaled_species$Growth_form_simple,
    PCA$x)
    
    #extract PCA loadings
  PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)
  
  return(
    list(
    PCAvalues = PCAvalues,
    PCAloadings = PCAloadings,
    PCAvariance = pca_var_perc,
    country_name = .y)
    )
})



# ---- 2. PLOT CHARTS   ----


pca_plots_list <- map(pca_data_list, function(.x) {
  
  PCAvalues <- .x$PCAvalues
  PCAloadings <- .x$PCAloadings
  PCAvariance <- .x$PCAvariance 
  country_name <- .x$country_name
  
  message("Generating plot for: ", country_name)
  
    #enlarge loading arrows
    arrow_scale <- 3.5 
  
plot_pca <- ggplot(data = PCAvalues, mapping = aes(x = PC1, y = PC2, color= status)) +
  coord_fixed() +  
  #coord_fixed(xlim = c(-5, 5), ylim = c(-5, 5)) +
    
    # A. Bins
    stat_density_2d(geom = "polygon", contour = TRUE, 
                    aes(x = PC1, y = PC2, fill = after_stat(level)),
                    inherit.aes = FALSE, bins = 10, alpha = 0.85, 
                    #colour = "grey90", linewidth = 0.01,
                    show.legend = FALSE) +
    scale_fill_distiller(name = "level", palette = "Greys", direction = 1, guide = "none") +
    
    new_scale_fill() + 
    
    # B. Points
  geom_jitter(data = subset(PCAvalues, status == "extant"),
              aes(color = Growth_form_simple), alpha = 0.3, size = 1, show.legend = FALSE) +
  geom_jitter(data = subset(PCAvalues, status == "extirpated"),
              shape = 17, alpha = 1, size = 2, show.legend = FALSE) +
  
    scale_color_manual(values = c(
      "extirpated" = "#9B4F96",  
      "tree"       = "#66C2A5",  
      "shrub"      = "#E5C494",  
      "herb"       = "#A6D854"   
    )) + 
  stat_ellipse(
    data = subset(PCAvalues, 
                  status == "extant" & 
                    !is.na(Growth_form_simple) & 
                    Growth_form_simple %in% names(which(table(Growth_form_simple) >= 4))
    ), 
    aes(x = PC1, y = PC2, color = Growth_form_simple),
    inherit.aes = FALSE,
    type = "t",
    level = 0.90,        
    linewidth = 0.8,   
    show.legend = FALSE
  ) +
    
    
    scale_fill_manual(values = c("extirpated" = "#9B4F96", "extant" = "gray65")) +
    # C. Loadings
    geom_segment(data = PCAloadings, 
                 aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
                 inherit.aes = FALSE, arrow = arrow(length = unit(0.15, "cm")),
                 colour = "black", linewidth = 0.4) +
    geom_text_repel(data = PCAloadings,
                    aes(x = PC1 * arrow_scale, y = PC2 * arrow_scale, label = Variables),
                    inherit.aes = FALSE, size = 2.5, fontface = "bold",
                    box.padding = 0.5, point.padding = 0.1, min.segment.length = Inf) +
    
    # D. Marginal density
    geom_xsidedensity(aes(y = after_stat(density), fill = status), 
                      alpha = 0.3, linewidth = 0.7, adjust = 0.9, show.legend = FALSE) + 
    geom_ysidedensity(aes(x = after_stat(density), fill = status), 
                      alpha = 0.3, linewidth = 0.7, adjust = 0.9, show.legend = FALSE) +
    
    labs(
      x = paste0("PC1 (", PCAvariance[1], "%)"),
      y = paste0("PC2 (", PCAvariance[2], "%)")
    ) +
    
    # E. Theme adjustments
    theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.line = element_blank(),
      legend.position = c(0.02, 0.02),
      legend.justification = c(0,0),
      legend.background = element_blank(), 
      legend.key.size = unit(0.3, "cm"), 
      legend.text = element_text(size = 8), 
      legend.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(family = "Georgia"),
      ggside.panel.scale.x = 0.2, ggside.panel.scale.y = 0.2, 
      ggside.panel.border = element_blank(), ggside.panel.background = element_blank(),
      ggside.axis.text = element_blank(), ggside.axis.ticks = element_blank()
    ) +
    ggtitle(country_name, )
  
  return(plot_pca)
})


# ---- 3. MERGE CHARTS WITH PATCHWORK AND SAVE   ----

selected<- pca_plots_list[c("Malta", "Croatia","Austria","Norway","Portugal", "Germany","Czech Republic","Poland")]


legend_data <- data.frame(
  Group = factor(c("Tree", "Shrub", "Herb", "Extirpated"),
                 levels = c("Tree", "Shrub", "Herb", "Extirpated")),
  x = 1, y = 1
)


legend_plot <- ggplot(legend_data, aes(x, y, color = Group, shape = Group, size = Group)) +
  geom_point() +
  scale_color_manual(values = c(
    "Tree"       = "#66C2A5",
    "Shrub"      = "#E5C494",
    "Herb"       = "#A6D854",
    "Extirpated" = "#9B4F96"
  )) +
  scale_shape_manual(values = c(
    "Tree" = 16, "Shrub" = 16, "Herb" = 16, "Extirpated" = 17
  )) +
  scale_size_manual(values = c(
    "Tree" = 3.5, "Shrub" = 3.5, "Herb" = 3.5, "Extirpated" = 4.5
  )) +
  labs(color = "Plant group & status", shape = "Plant group & status",
       size = "Plant group & status") +
  guides(color = guide_legend(nrow = 1)) +
  theme_void(base_size = 10) +
  theme(
    legend.direction = "horizontal",
    legend.title      = element_text(size = 10, face = "bold", family = "sans", colour = "black"),
    legend.text       = element_text(size = 9,  family = "sans", colour = "black"),
    legend.key        = element_blank(),
    legend.spacing.x  = unit(0.35, "cm")
  )

standalone_legend <- ggdraw(get_legend(legend_plot))

ggsave("2_result_charts/Standalone_Legend.pdf",
       plot = standalone_legend, device = "pdf",
       width = 140, height = 18, units = "mm")



selected_dashboard <- ((wrap_plots(selected, ncol = 4)) / standalone_legend) + 
  plot_layout(heights = c(0.9, 0.10))


all_dashboard <- ((wrap_plots(pca_plots_list, ncol = 4)) / standalone_legend) + 
  plot_layout(heights = c(0.95, 0.05))


ggsave("2_result_charts/PCA_Dashboard_selected_Countries_ellipse.pdf",
       plot = selected_dashboard, device = "pdf",
       width =440 , height =  220, units = "mm")

ggsave("2_result_charts/PCA_Dashboard_all_Countries_ellipse.pdf",
       plot = all_dashboard, device = "pdf",
       width =440 , height =  660, units = "mm")

pdf("2_result_charts/PCA_Dashboard_all_single_Countries_ellipse.pdf", 
    width = 8.27, height = 8.27)
walk(pca_plots_list, function(p) {
  
  page <- (wrap_elements(p) / standalone_legend)+ 
    plot_layout(heights = c(0.9, 0.1))
  print(page)
})
dev.off()
