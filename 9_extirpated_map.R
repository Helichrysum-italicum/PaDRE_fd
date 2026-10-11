 library(ggplot2)
library(tidyverse)
library(sf)
library(giscoR)
library(visdat)
library(dplyr)
library(patchwork)

setwd("xxx")

all_sp <- read.csv("data/0_extirpation_list/extinct_species_RL_WCVP_updated.csv", encoding = "Latin1") %>%
  select(accepted_taxon_name, family, ISO3_Code) %>%
  distinct(accepted_taxon_name, .keep_all = TRUE)

families <- all_sp%>%
  count(family, sort=TRUE)%>%
  mutate( family = fct_lump_n(family, n = 20, w = n, other_level = "Other"))%>%
  group_by(family)%>%
  summarise(n = sum(n))%>%
  mutate( pct = n / sum(n)*100)%>%
  filter(family != "Other")

families <- ggplot(families, aes(x = reorder(family, n), y = n, fill = n)) +
  geom_col(show.legend = FALSE) + 
  coord_flip() +
  scale_fill_gradient(
    low = "#e0e0e0",
    high = "#003366",
    breaks = c(0, 8,6)) +
  geom_text(aes(label = paste0(round(pct, 1), "%")), 
            hjust = -0.2, size = 4, family = "Georgia") +
  labs(
    title = "20 most abundant families",
    subtitle = "Based on extirpated species dataset",
    x = NULL, 
    y = "Species Count"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    text = element_text(family = "Georgia"),
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank() 
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))


families






all_sp <- read.csv("data/0_extirpation_list/extinct_species_RL_WCVP_updated.csv", encoding = "Latin1")
  

geodat <- all_sp %>%
  select(accepted_taxon_name, ISO3_Code)%>%
  group_by(ISO3_Code) %>%
  summarise(n_species = n_distinct(accepted_taxon_name))

world_map <- gisco_get_countries(resolution = "60", year = "2020")

world_robinson <- world_map %>%
  left_join(geodat, by = c("ISO3_CODE" = "ISO3_Code")) %>%
  st_transform(crs = "ESRI:54030")

extirpation_map <- ggplot(world_robinson) +
  geom_sf(aes(fill = n_species), color = "white", size = 0.05) +
  scale_fill_gradient(
    low = "#e0e0e0",
    high = "#003366",    
    trans = "log10",
    na.value = "#f5f5f5", 
    breaks = c(1, 725),
    labels = c("Minimum Extirpation", "Maximum Extirpation"),
    guide = guide_colorbar(
      title = "Extirpation Density",
      title.position = "bottom",
      barwidth = 15,
      ticks = FALSE,  
      draw.ulim = FALSE,   
      draw.llim = FALSE,
      frame.colour = "black",
      label.theme = element_text(face = "italic", size = 9)
    )
  )+
  theme_void() + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12, family = "Georgia"), 
    legend.title = element_text(size = 12, family = "Georgia", face = "bold"),
    plot.title = element_text(family = "Georgia", face = "bold", hjust = 0.5, size = 15),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(title = "Distribution of extirpated species")
extirpation_map


# Save the family bar chart
ggsave(
  filename = paste0("data/accessory_visual/family_distribution_ext.pdf"),
  plot = families, device = cairo_pdf, width = 20, height = 10, units = "in", limitsize = FALSE
)

# Save the global extirpation map
ggsave(
  filename = paste0("data/accessory_visual/global_map_ext.pdf"),
  plot = extirpation_map, device = cairo_pdf, width = 20, height = 10, units = "in", limitsize = FALSE
)
