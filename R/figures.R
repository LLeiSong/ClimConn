# ============================================================
# Script: figures.R
#
# Purpose:
#   Generate all manuscript figures from pre-computed results:
#     - Fig. 3: Block size vs. dispersal radius correlation curves
#     - Fig. 4 & S: Connectivity ensembles and typology maps
#     - Fig. S1: Model evaluation (AUC, reconstruction correlations)
#     - Fig. S2: Species loss over time and current richness map
#     - Fig. S3: Dispersal rate distribution by mammal order.
#
# Inputs:
#   data/                                 core data and shapefiles
#     colombia.geojson
#     Env/M_variable/vars.tif
#     smr_species.csv
#     suit_cors.csv
#     dispersal/species_dispersal_rate.csv
#   results/omni/                         main Omniscape connectivity
#     conn_lyrs.tif
#   results/omni_single/                  single-run connectivity
#     conn_lyrs.tif
#   results/block_radius/                 block–radius sensitivity
#     correlation_block_size.csv
#     fitted_curves.csv
#   results/typology/                     connectivity typology
#     conn_types_*.tif
#     connectivity_typology.csv
#
# Outputs (in results/figures/):
#   Figure3_block_radius.png
#   Figure4_conn_all.png
#   Figure_s_conn.png
#   Figure5_conn_types.png
#   Figure_s_conn_types_simple.png
#   Figure_s_model_eval.png
#   Figure_s_species.png
#   Figure_s_dispersal.png
#
# Author: Lei Song <lei.song@rutgers.edu>
# Last updated: 2025-11-21
# ============================================================

# Load libraries
library(ggplot2)
library(tidyterra)
library(ggsci)
library(here)
library(ggpubr)
library(smoothr)
library(stars)
library(dplyr)
library(ini)
library(terra)
library(stringr)
library(patchwork)
library(ggsankey)
library(ggpattern)
library(elevatr)
library(rnaturalearth)
library(ggnewscale)
sf_use_s2(FALSE)

# Setting
root_dir <- here()
data_dir <- file.path(root_dir, "data")
result_dir <- file.path(root_dir, "results/omni")
result_cmp <- file.path(root_dir, "results/omni_single")
fig_dir <- file.path(root_dir, "results/figures")

# Figures
# Fig.2 uses the following results as example
sp <- "Aotus.griseimembra"
scenario <- "ssp126"

#### Figure S4 ####
bs_sens_dir <- file.path(root_dir, "results/block_radius")
sens_cors <- read.csv(file.path(bs_sens_dir, "correlation_block_size.csv"))
fitted_curves <- read.csv(file.path(bs_sens_dir, "fitted_curves.csv"))

ggplot(sens_cors, 
       aes(x = radius_size, y = cor, color = as.factor(block_size))) +
  geom_point(size = 0.7) + 
  geom_line(
    data = fitted_curves, 
    aes(radius_size, pred, color = as.factor(block_size)),
    linewidth = 0.9) +
  scale_color_npg(name = "Block size (pixel)") +
  xlab("Dispersal distance (pixel)") +
  ylab("Correlation with no block") +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  theme_pubclean()

ggsave(file.path(fig_dir, "Figure3_block_radius.png"),
       width = 5, height = 4, dpi = 300, bg = "white")

#### Figure 3 ####
# Load data
conns <- rast(file.path(result_dir, "conn_lyrs_a_mean.tif"))
col <- st_read(file.path(data_dir, "colombia.geojson")) %>% 
  st_transform(st_crs(conns))
conns <- conns %>% crop(col) %>% mask(col)
aoi <- st_read(file.path(result_dir, "fig3_aoi.geojson")) %>% st_buffer(0.3)

conns_to_cmp <- rast(file.path(result_cmp, "conn_lyrs_a_mean.tif")) %>% 
  trim() %>% crop(col) %>% mask(col)

# Plot
nms <- names(conns)
nms <- nms[str_detect(nms, "ssp126")]

lyrs <- c(trim(conns_to_cmp), trim(subset(conns, nms)))
names(lyrs) <- letters[1:nlyr(lyrs)]

ggplot() +
  geom_spatraster(data = lyrs) +
  scale_fill_grass_c(
    "Connectivity\n(Mean current flow)", palette = "inferno", 
    na.value = "transparent",
    labels = scales::label_number(accuracy = 0.01)) +
  facet_wrap(~lyr, ncol = 4, nrow = 2, strip.position = "left") +
  geom_sf(data = aoi, color = "white", fill = "transparent", 
          linewidth = 0.3) +
  geom_sf(data = col, color = "black", fill = "transparent", 
          linewidth = 0.4) +
  theme_void() +
  theme(legend.position = "right",
        legend.key.width = unit(0.2, "cm"),
        legend.key.height = unit(0.6, "cm"),
        strip.text = element_text(
          size = 12, face = "bold", vjust = 0.9, hjust = 1),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        panel.spacing = unit(0, "lines")) + 
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))

ggsave(file.path(fig_dir, "Figure3_compare.png"),
       width = 6.5, height = 3.5, dpi = 300, bg = "white")

#### Figure 4 ####
# Load data
a_means <- rast(file.path(result_dir, "conn_lyrs_a_mean.tif"))
col <- st_read(file.path(data_dir, "colombia.geojson")) %>% 
  st_transform(st_crs(a_means))
a_means <- a_means %>% crop(col) %>% mask(col)

g_means <- rast(file.path(result_dir, "conn_lyrs_g_mean.tif")) %>% 
  crop(col) %>% mask(col)

sums <- rast(file.path(result_dir, "conn_lyrs_sum.tif")) %>% 
  crop(col) %>% mask(col)

# Single species
sp <- "Puma.concolor"
conn_sp <- rast(
  file.path(result_dir, "ssp126", sp, "2081-2100/cum_currmap.tif"))
conn_sp <- trim(stretch(conn_sp, minv = 0, maxv = 1))

# Plot
nms <- names(a_means)
nms <- nms[str_detect(nms, "2081-2100_ssp126")]

lyrs <- c(trim(subset(a_means, nms)), trim(subset(g_means, nms)),
          trim(subset(sums, nms)))
names(lyrs) <- letters[1:nlyr(lyrs)]

g1 <- ggplot() +
  geom_spatraster(data = conn_sp) +
  scale_fill_grass_c(
    expression(atop(Connectivity, italic("(Puma concolor)"))),
    palette = "inferno", 
    na.value = "transparent",
    labels = scales::label_number(accuracy = 0.1)) +
  geom_sf(data = col, color = "black", fill = "transparent", 
          linewidth = 0.4) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.key.width = unit(0.7, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.margin = margin(0, 0, 0, 0, "mm")) + 
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))

g2 <- ggplot() +
  geom_spatraster(data = trim(subset(a_means, nms)) * 10) +
  scale_fill_grass_c(
    expression(atop(Connectivity, "(Arithmetic mean 10"^-1*")")), 
    palette = "inferno", 
    na.value = "transparent",
    labels = scales::label_number(accuracy = 1)) +
  geom_sf(data = col, color = "black", fill = "transparent", 
          linewidth = 0.4) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.key.width = unit(0.7, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.margin = margin(0, 0, 0, 0, "mm")) + 
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))

g3 <- ggplot() +
  geom_spatraster(data = trim(subset(g_means, nms)) * 100) +
  scale_fill_grass_c(
    expression(atop(Connectivity, "(Geometric mean 10"^-2*")")), 
    palette = "inferno", 
    na.value = "transparent",
    labels = scales::label_number(accuracy = 1)) +
  geom_sf(data = col, color = "black", fill = "transparent", 
          linewidth = 0.4) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.key.width = unit(0.7, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.margin = margin(0, 0, 0, 0, "mm")) + 
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))

g4 <- ggplot() +
  geom_spatraster(data = trim(subset(sums, nms)) / 10) +
  scale_fill_grass_c(
    expression(atop(Connectivity, "(Sum 10"^1*")")), 
    palette = "inferno", 
    na.value = "transparent",
    labels = scales::label_number(accuracy = 1)) +
  geom_sf(data = col, color = "black", fill = "transparent", 
          linewidth = 0.4) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.key.width = unit(0.7, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.margin = margin(0, 0, 0, 0, "mm")) + 
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))

ggarrange(g1, g2, g3, g4, nrow = 1, labels = letters[1:4],
          font.label = list(size = 12, face = "bold"))

ggsave(file.path(fig_dir, "Figure4_ensemble.png"),
       width = 6.5, height = 3.5, dpi = 300, bg = "white")

#### Figure 5 and Ss ####
## Figure Ss
### Continuous maps
# Re-order the names
nms <- names(a_means)
nms <- c(nms[str_detect(nms, "ssp126")], 
         nms[str_detect(nms, "ssp370")], 
         nms[str_detect(nms, "ssp585")])

items <- list("a_mean" = a_means, "g_mean" = g_means, "sum" = sums)

for (i in names(items)){
  lyrs <- trim(subset(items[[i]], nms))
  names(lyrs) <- letters[1:nlyr(lyrs)]
  
  if (i == "a_mean"){
    label <- "Connectivity (Arithmetic mean)"
    acc <- 0.01
  } else if (i == "g_mean"){
    label <- "Connectivity (Geometric mean)"
    acc <- 0.01
  } else{
    label <- "Connectivity (Sum)"
    acc <- 1
  }
  
  ggplot() +
    geom_spatraster(data = lyrs) +
    scale_fill_grass_c(
      label, palette = "inferno", 
      na.value = "transparent",
      labels = scales::label_number(accuracy = acc)) +
    facet_wrap(~lyr, ncol = 4, nrow = 3, strip.position = "left") +
    geom_sf(data = col, color = "black", fill = "transparent", 
            linewidth = 0.4) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.key.width = unit(2, "cm"),
          legend.key.height = unit(0.2, "cm"),
          strip.text = element_text(
            size = 12, face = "bold", vjust = 0.9, hjust = 1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          panel.spacing = unit(0, "lines")) + 
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
  
  ggsave(file.path(fig_dir, sprintf("Figure_s_%s.png", i)),
         width = 6.5, height = 7.5, dpi = 300, bg = "white")
}

### Categorical maps
typology_dir <- file.path(root_dir, "results/classification")
classes <- read.csv(file.path(typology_dir, "connectivity_typology.csv"))

colors <- c("#9bc1bc", "#127475",
            "#f6ae2d", "#f26419", "#f5dfbb")

items <- c("a_mean", "g_mean", "sum")
for (item in items){
  fnames <- list.files(
    typology_dir, pattern = sprintf("%s.tif$", item), full.names = TRUE)
  lyrs <- rast(fnames)
  names(lyrs) <- str_extract(basename(fnames), "[0-9]{4}-[0-9]{4}_ssp[0-9]{3}")
  lyrs <- trim(subset(lyrs, nms))
  
  levels(lyrs) <- lapply(1:nlyr(lyrs), function(x) classes[, c("ID", "name")]) 
  names(lyrs) <- letters[1:nlyr(lyrs)]
  
  ggplot() +
    geom_spatraster(data = lyrs, na.rm = TRUE) +
    scale_fill_manual(
      name = "Connectivity type", values = colors, 
      breaks = c(levels(lyrs[[1]])[[1]]$name[-1], "Background"),
      na.translate = FALSE) +
    facet_wrap(~lyr, ncol = 4, nrow = 3, strip.position = "left") +
    geom_sf(data = col, color = "black", fill = "transparent", 
            linewidth = 0.4) +
    theme_void() +
    theme(legend.position = "bottom",
          strip.text = element_text(
            size = 12, face = "bold", vjust = 0.9, hjust = 1),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12, face = "bold"),
          panel.spacing = unit(0, "lines")) + 
    guides(fill = guide_legend(
      title.position = "top", title.hjust = 0.5, nrow = 2))
  
  ggsave(file.path(fig_dir, sprintf("Figure_s_types_%s.png", item)),
         width = 6.5, height = 7.5, dpi = 300, bg = "white")
}

## Table Ss
# Check area numbers
items <- c("a_mean", "g_mean", "sum")
areas <- lapply(items, function(item){
  fnames <- list.files(
    typology_dir, pattern = sprintf("%s.tif$", item), full.names = TRUE)
  lyrs <- rast(fnames)
  names(lyrs) <- str_extract(basename(fnames), "[0-9]{4}-[0-9]{4}_ssp[0-9]{3}")
  lyrs <- trim(subset(lyrs, nms))
  
  lapply(names(lyrs), function(nm){
    zonal(cellSize(lyrs[[nm]], unit = "km"), 
          lyrs[[nm]], fun = "sum") %>% data.frame() %>% 
      arrange(names(.)[1]) %>% 
      select(area) %>% mutate(area = area / sum(area) * 100) %>% 
      rename({{nm}} := "area")
  }) %>% bind_cols() %>% 
    mutate(type = levels(lyrs[[1]])[[1]] %>% 
             arrange(value) %>% pull(class), .before = 1) %>% 
    mutate(method = item, .before = 1)
}) %>% bind_rows()

areas <- areas %>% 
  mutate(type = factor(
    type, levels = c("F1", "F2", "I1", "I2", "Background"),
    labels = c("Facilitative area", "Highly facilitative area",
               "Impeded area", "Highly impeded area", "Background"))) %>% 
  arrange(method, type)
  
write.csv(areas, file.path(fig_dir, "areas_classes.csv"), row.names = FALSE)

## Figure 5.
# Plot
cont_maps <- list(
  "a_mean" = subset(a_means, names(a_means)[str_detect(names(a_means), "ssp585")][c(1, 4)]), 
  "g_mean" = subset(g_means, names(g_means)[str_detect(names(g_means), "ssp585")][c(1, 4)]), 
  "sum" = subset(sums, names(sums)[str_detect(names(sums), "ssp585")][c(1, 4)]))

items <- c("a_mean", "g_mean", "sum")

fig_list <- lapply(items, function(item){
  fnames <- list.files(
    typology_dir, pattern = sprintf("ssp585_%s.tif$", item), full.names = TRUE)
  lyrs <- rast(fnames[c(1, 4)])
  levels(lyrs) <- lapply(1:nlyr(lyrs), function(x) classes[, c("ID", "name")])
  names(lyrs) <- str_extract(
    basename(fnames[c(1, 4)]), "[0-9]{4}-[0-9]{4}")
  
  # Get the stable classes
  msk <- mask(lyrs[[1]], lyrs[[1]] == lyrs[[2]], 
              maskvalues = 0, updatevalue = NA)
  msk[msk == 0] <- NA; msk <- patch_clean(msk); msk[msk == 0] <- NA
  msk <- as.polygons(msk) %>% st_as_sf()
  
  cont_lyrs <- cont_maps[[item]]
  
  if (item == "a_mean"){
    acc <- 10
    n <- c(0, 1, 2)
  } else if (item == "g_mean"){
    acc <- 100
    n <- c(0.001, 1, 2)
  } else{
    acc <- 0.1
    n <- c(0, 2, 4)
  }
  
  g10 <- ggplot() +
    geom_spatraster(data = cont_lyrs[[1]] * acc) +
    scale_fill_grass_c(
      "", palette = "inferno", 
      na.value = "transparent", 
      breaks = n,
      labels = scales::label_number(accuracy = 1)) +
    geom_sf(data = col, color = "black", fill = "transparent", 
            linewidth = 0.4) +
    ggtitle(sprintf("Conn.\u00d7%s", acc)) + theme_void() +
    theme(legend.position = "right",
          legend.key.width = unit(0.1, "cm"),
          legend.key.height = unit(0.25, "cm"),
          legend.spacing = unit(0, "lines"),
          legend.text = element_text(size = 9),
          plot.title = element_text(size = 9, hjust = 1),
          plot.title.position = "plot") +
    guides(fill = guide_colourbar(title.position = "right"))
  
  g1 <- ggplot() +
    geom_spatraster(data = hill, show.legend = FALSE) +
    scale_fill_gradient(low = "black", high = "white", na.value = NA) +
    new_scale_fill() + # From 'ggnewscale' package
    geom_spatraster(data = lyrs[[1]], na.rm = TRUE, alpha = 0.6) +
    scale_fill_manual(
      name = "Connectivity type", values = colors, 
      breaks = c(levels(lyrs[[1]])[[1]]$name[-1], "Background"),
      na.translate = FALSE) +
    labs(caption = names(lyrs[[1]])) +
    geom_sf_pattern(
      data = msk,
      fill = NA,
      color = NA, 
      pattern = "crosshatch",
      pattern_colour = "black",
      pattern_fill = "black", 
      pattern_spacing = 0.01, 
      pattern_angle = 45,
      pattern_size = 0.01) +
    geom_sf(data = col, color = "black", fill = "transparent", 
            linewidth = 0.4) +
    theme_void() +
    theme(legend.position = "none",
          panel.spacing = unit(0, "lines"),
          plot.caption = element_text(
            hjust = 0.15,
            vjust = 20,
            size = 10),
          plot.title = element_text(color = "transparent"))
  
  # Combine the plots
  g1 <- g1 + 
    inset_element(g10, left = 0.55, bottom = 0.61, right = 1, top = 1,
                  align_to = 'panel')
  
  g20 <- ggplot() +
    geom_spatraster(data = cont_lyrs[[2]] * acc) +
    scale_fill_grass_c(
      "", palette = "inferno", 
      na.value = "transparent", 
      breaks = n,
      labels = scales::label_number(accuracy = 1)) +
    geom_sf(data = col, color = "black", fill = "transparent", 
            linewidth = 0.4) +
    ggtitle(sprintf("Conn.\u00d7%s", acc)) + theme_void() +
    theme(legend.position = "right",
          legend.key.width = unit(0.1, "cm"),
          legend.key.height = unit(0.25, "cm"),
          legend.spacing = unit(0, "lines"),
          legend.text = element_text(size = 9),
          plot.title = element_text(size = 9, hjust = 1),
          plot.title.position = "plot") +
    guides(fill = guide_colourbar(title.position = "right"))
  
  g2 <- ggplot() +
    geom_spatraster(data = hill, show.legend = FALSE) +
    scale_fill_gradient(low = "black", high = "white", na.value = NA) +
    new_scale_fill() + # From 'ggnewscale' package
    geom_spatraster(data = lyrs[[2]], na.rm = TRUE, alpha = 0.6) +
    scale_fill_manual(
      name = "Connectivity type", values = colors, 
      breaks = c(levels(lyrs[[1]])[[1]]$name[-1], "Background"),
      na.translate = FALSE) +
    labs(caption = names(lyrs[[2]])) +
    geom_sf_pattern(
      data = msk,
      fill = NA,
      color = NA, 
      pattern = "crosshatch",
      pattern_colour = "black",
      pattern_fill = "black", 
      pattern_spacing = 0.01, 
      pattern_angle = 45,
      pattern_size = 0.01) +
    geom_sf(data = col, color = "black", fill = "transparent", 
            linewidth = 0.4) +
    theme_void() +
    theme(legend.position = "none",
          panel.spacing = unit(0, "lines"),
          plot.caption = element_text(
            hjust = 0.15,
            vjust = 20,
            size = 10),
          plot.title = element_text(color = "transparent"))
  
  # Combine the plots
  g2 <- g2 + 
    inset_element(g20, left = 0.55, bottom = 0.61, right = 1, top = 1,
                  align_to = 'panel')
  
  ggarrange(g1, NULL, g2, nrow = 3, heights = c(1, -0.1, 1))
})


# Make the legend
legend <- ggplot() +
  geom_spatraster(data = lyrs, na.rm = TRUE) +
  scale_fill_manual(
    name = "Connectivity type",
    values = colors,
    breaks = c(levels(lyrs[[1]])[[1]]$name[-1], "Background"),
    na.translate = FALSE
  ) +
  facet_wrap(~ lyr, ncol = 4, nrow = 3) +
  geom_sf_pattern(
    data = col,
    aes(pattern = "Overlay"),
    fill = NA,
    color = NA,
    pattern_colour = "black",
    pattern_fill = "black",
    pattern_spacing = 0.01,
    pattern_angle = 45,
    pattern_size = 0.01,
    inherit.aes = FALSE,
    show.legend = TRUE
  ) +
  scale_pattern_manual(name = "Stable area",
                       values = c("Overlay" = "crosshatch")) +
  theme_void() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    panel.spacing = unit(0, "lines")
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 2,
      override.aes = list(
        pattern = "none",
        pattern_colour = NA,
        pattern_fill = NA
      )
    ),
    pattern = guide_legend(title.position = "top", title.hjust = 0.5)
  )

legend <- get_legend(legend)

g <- ggarrange(plotlist = fig_list, ncol = 3, 
          labels = c("a (Arithmetic mean)", "b (Geometric mean)", "c (Sum)"), 
          label.x = 0.5, label.y = 1, hjust = 0.5, vjust = 1,
          font.label = list(size = 12))
ggarrange(g, NULL, legend, nrow = 3, heights = c(6.2, -0.2, 0.7))

ggsave(file.path(fig_dir, "Figure5_conn_types_simple.png"),
       width = 6.5, height = 7, dpi = 300, bg = "white")

#### Figure S1 ####
# Colombia geography map
# 1. Get Colombia boundary
colombia <- ne_countries(scale = "medium", country = "colombia", returnclass = "sf")
rivers_global <- ne_download(
  scale = 10, type = "rivers_lake_centerlines", 
  category = "physical", returnclass = "sf")
rivers <- st_intersection(rivers_global, colombia)

# 2. Fetch Elevation (DEM)
# z = 6 is a good balance between detail and speed; increase for higher resolution
elev <- get_elev_raster(colombia, z = 6, clip = "locations")
elev <- rast(elev)
names(elev) <- "elevation"

# Define a rough line for the Andes through Colombia
ranges <- data.frame(
  name = c("Cordillera Occidental", "Cordillera Central", "Cordillera Oriental"),
  x = c(-77.2, -75.8, -74),
  y = c(3.0, 4.5, 5),
  angle = c(68, 68, 55))

basins <- data.frame(
  name = c("Orinoco Basin", "Amazon Basin"),
  x = c(-70, -72),
  y = c(5, 0),
  angle = c(25, 15))

# 2. Calculate Hillshade
slope <- terrain(elev, "slope", unit = "radians")
aspect <- terrain(elev, "aspect", unit = "radians")
hillshade <- shade(slope, aspect, angle = 45, direction = 315)

# 3. Plot with tidyterra
ggplot() +
  geom_spatraster(data = hillshade) +
  scale_fill_gradient(
    low = "black", high = "white", guide = "none", na.value = NA) +
  new_scale_fill() +
  geom_spatraster(data = elev, aes(fill = elevation), alpha = 0.8) +
  scale_fill_hypso_tint_c(palette = "dem_poster", na.value = NA) +
  labs(fill = "Elevation (m)") +
  geom_sf(data = rivers, color = "#5dade2", linewidth = 0.3) +
  geom_sf(data = colombia, fill = "transparent", 
          color = "black", linewidth = 0.6) +
  geom_text(data = ranges, aes(x = x, y = y, label = name, angle = angle),
            fontface = "bold", size = 4, color = "black") +
  geom_text(data = basins, aes(x = x, y = y, label = name, angle = angle),
            size = 6, color = "white") +
  theme_void() +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11)) +
  coord_sf(expand = FALSE)

ggsave(file.path(fig_dir, "Figure_s_col_geo.png"),
       width = 6, height = 6, dpi = 500)

#### Figure S2 ####
sps <- list.files(file.path(data_dir, "mamiferos"))
evals <- lapply(sps, function(sp){
  sp_dir <- file.path(data_dir, "mamiferos", sp)
  eval <- read.csv(file.path(sp_dir, "eval_results_enmeval", "best_models.csv"))
  data.frame(sp = sp, auc = c(mean(eval$auc.train), mean(eval$auc.val.avg)), 
             type = c("Training", "Test"))
}) %>% bind_rows() %>% 
  mutate(type = factor(type, levels = c("Training", "Test")))

g1 <- ggplot(data = evals, 
            aes(x = auc, after_stat(density), fill = type)) +
  geom_density(alpha = 0.8) +
  xlab("AUC") +
  ylab("Density") +
  scale_fill_bmj(name = "Type") +
  theme_pubclean(base_size = 12) +
  theme(axis.text = element_text(color = "black"))

suit_cor <- read.csv(file.path(data_dir, "suit_cors.csv"))

g2 <- ggplot(data = suit_cor, aes(x = cor)) +
  geom_boxplot(outliers = FALSE, fill = "#C62584") +
  xlab("Correlation") + ylab("") +
  theme_pubclean(base_size = 12) +
  theme(axis.text = element_text(color = "black"),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_line(color = "white"),
        axis.title.y = element_text(angle = 0),
        axis.ticks.y = element_blank())

ggarrange(g1, g2, nrow = 2, heights = c(6, 2), labels = c("a", "b"))

ggsave(file.path(fig_dir, "Figure_s_model_eval.png"),
       width = 4, height = 4, dpi = 500, bg = "white")

#### Figure S3 ####
# Species loss
years <- c("1970-2000", "2021-2040", 
           "2041-2060", "2061-2080", "2081-2100")
ssps <- c("ssp126", "ssp370", "ssp585")

smr_species <- read.csv(file.path(data_dir, "smr_species.csv"))
dat <- smr_species %>% filter(big_habitat == 1) %>% 
  mutate(time_periods = ifelse(
    is.na(time_periods), "1970-2000", time_periods)) %>% 
  group_by(time_periods, ssp) %>% summarise(n = n()) %>% 
  mutate(time_periods = factor(time_periods, levels = years, labels = years),
         ssp = factor(ssp, levels = ssps, labels = ssps))

g1 <- ggplot(data = dat) +
  geom_point(aes(x = time_periods, y = n / max(n) * 100, color = ssp)) + 
  geom_line(aes(x = as.integer(time_periods), 
                y = n / max(n) * 100, color = ssp)) +
  xlab("Time period") + ylab("Percentage of species") +
  scale_color_bmj(name = "Climate scenario") + 
  theme_pubclean(base_size = 10) + 
  theme(axis.text = element_text(color = "black"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

# Species richness
cur_path <- "reconstruct/current/MAXENT/%s_10_MAXENT.tif"
cur_suit_path <- "reconstruct/current/MAXENT/%s_MAXENT.tif"

col <- st_read(file.path(data_dir, "colombia.geojson"))
template <- rast(file.path(data_dir, "Env/M_variable/vars.tif"), lyrs = 1)
values(template) <- 0
template <- mask(crop(template, col), col)

richness <- template
suitability <- template

sps <- smr_species %>% filter(big_habitat == 1) %>% 
  mutate(time_periods = ifelse(
    is.na(time_periods), "1970-2000", time_periods)) %>% 
  filter(time_periods == "1970-2000") %>% pull(sp) %>% unique()

for (sp in sps){
  sp_dir <- file.path(data_dir, "mamiferos", sp)
  
  cur <- rast(file.path(
    sp_dir, sprintf(cur_path, gsub("\\.", "_", sp)))) %>% 
    crop(richness) %>% mask(richness) %>% 
    terra::extend(richness, fill = NA)
  
  suit <- rast(file.path(
    sp_dir, sprintf(cur_suit_path, gsub("\\.", "_", sp)))) %>% 
    crop(suitability) %>% mask(suitability) %>% 
    terra::extend(suitability, fill = NA)
  
  richness <- sum(c(richness, cur), na.rm = TRUE)
  suitability <- sum(c(suitability, suit), na.rm = TRUE)
}

suitability <- suitability / length(sps)

g2 <- ggplot() +
  geom_spatraster(data = richness) +
  scale_fill_grass_c(
    name = "No of species", 
    palette = "plasma", na.value = "transparent") +
  geom_sf(data = col, fill = "transparent") +
  theme_void() +
  guides(fill = guide_colorbar(
    title.hjust = 0.5, title.position = "top",
    barheight = unit(3, "mm"))) +
  theme(legend.position = "bottom",
        plot.margin = unit(c(0, -1, 0, -1), "cm"))

ggarrange(g2, g1, ncol = 2, labels = c("a", "b"), widths = c(1, 2))

ggsave(file.path(fig_dir, "Figure_s_species.png"),
       width = 6, height = 3, dpi = 500, bg = "white")

#### Figure S4 ####
dispersal_rate <- read.csv(file.path(
  data_dir, "dispersal/species_dispersal_rate.csv"))

dispersal_rate <- rbind(dispersal_rate, 
                        dispersal_rate %>% mutate(order = "All"))

dispersal_rate <- dispersal_rate %>% group_by(order) %>% 
  summarise(q1 = quantile(dispersal_rate, 0.25),
            q3 = quantile(dispersal_rate, 0.75),
            median = median(dispersal_rate))

order_order <- dispersal_rate %>% arrange(-median) %>% pull(order)
order_order <- c("All", setdiff(order_order, "All"))

dispersal_rate <- dispersal_rate %>% 
  mutate(order = factor(order, levels = order_order, labels = order_order))

ggplot(data = dispersal_rate) +
  geom_errorbar(
    data = dispersal_rate,
    aes(x = order, ymin = q1, ymax = q3), color = "#EABE51",
    width = 0, linewidth = 1.5) +
  geom_point(
    data = dispersal_rate,
    aes(y = median, x = order),
    color ="white", size = 2) +
  geom_point(
    data = dispersal_rate,
    aes(y = median, x = order), size = 1.6, color = "black") +
  coord_flip() + 
  labs(y = "Dispersal rate (km/year)",
       x = element_blank()) +
  theme_pubclean(base_size = 12) +
  theme(axis.text = element_text(
    color = "black"),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(
      linetype = "dotted", color = "lightgrey"),
    plot.margin = unit(c(0, 0.5, 0.2, 0.5), "cm"))

ggsave(file.path(fig_dir, "Figure_s_dispersal.png"),
       width = 5, height = 3.5, dpi = 500, bg = "white")
