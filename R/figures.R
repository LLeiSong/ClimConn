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
#   Figure5_conn_types.png (optional, commented in code)
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

#### Figure 3 ####
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

#### Figure 4 and Ss ####
# Load data
conns <- rast(file.path(result_dir, "conn_lyrs.tif"))
col <- st_read(file.path(data_dir, "colombia.geojson")) %>% 
  st_transform(st_crs(conns))
conns <- conns %>% crop(col) %>% mask(col)

conns_to_cmp <- rast(file.path(result_cmp, "conn_lyrs.tif")) %>% trim() %>% 
  crop(col) %>% mask(col)

# Plot
nms <- names(conns)
nms <- nms[str_detect(nms, "ssp126")]

lyrs <- c(trim(conns_to_cmp), trim(subset(conns, nms)))
names(lyrs) <- letters[1:nlyr(lyrs)]

g1 <- ggplot() +
  geom_spatraster(data = lyrs) +
  scale_fill_grass_c(
    "Connectivity (Cummulated current flow)", palette = "inferno", 
    na.value = "transparent",
    labels = scales::label_number(accuracy = 0.01)) +
  facet_wrap(~lyr, ncol = 4, nrow = 2, strip.position = "left") +
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

# ggsave(file.path(fig_dir, "Figure4_conns.png"),
#        width = 6.5, height = 5.3, dpi = 300, bg = "white")

nms <- names(conns)
nms <- c(nms[str_detect(nms, "ssp370")], nms[str_detect(nms, "ssp585")])

lyrs <- trim(subset(conns, nms))
names(lyrs) <- letters[1:nlyr(lyrs)]

ggplot() +
  geom_spatraster(data = lyrs) +
  scale_fill_grass_c(
    "Connectivity (Cummulated current flow)", palette = "inferno", 
    na.value = "transparent",
    labels = scales::label_number(accuracy = 0.01)) +
  facet_wrap(~lyr, ncol = 4, nrow = 2, strip.position = "left") +
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

ggsave(file.path(fig_dir, "Figure_s_conn.png"),
       width = 6.5, height = 5.3, dpi = 300, bg = "white")

## Another part of Figure 4.
# Load data
typology_dir <- file.path(root_dir, "results/typology")
fnames <- list.files(typology_dir, pattern = ".tif$", full.names = TRUE)
lyrs <- rast(fnames)
names(lyrs) <- str_extract(basename(fnames), "[0-9]{4}-[0-9]{4}_ssp[0-9]{3}")

col <- st_read(file.path(data_dir, "colombia.geojson")) %>% 
  st_transform(st_crs(lyrs))

lyrs <- mask(crop(lyrs, col), col)
classes <- read.csv(file.path(typology_dir, "connectivity_typology.csv"))

colors <- c("#9bc1bc", "#127475",
            "#f6ae2d", "#f26419", "#f5dfbb")

# Check area numbers
areas <- lapply(names(lyrs), function(nm){
  zonal(cellSize(lyrs[[nm]], unit = "km"), 
        lyrs[[nm]], fun = "sum") %>% data.frame() %>% 
    arrange(names(.)[1]) %>% 
    select(area) %>% mutate(area = area / sum(area) * 100) %>% 
    rename({{nm}} := "area")
}) %>% bind_cols() %>% 
  mutate(type = levels(lyrs[[1]])[[1]] %>% 
           arrange(value) %>% pull(class), .before = 1)

# Plot
nms <- names(lyrs)
dat <- subset(lyrs, nms[str_detect(nms, "ssp126")])
levels(dat) <- lapply(1:nlyr(dat), function(x) classes[, c("ID", "name")]) 
names(dat) <- letters[9:(8 + nlyr(dat))]

g2 <- ggplot() +
  geom_spatraster(data = dat, na.rm = TRUE) +
  scale_fill_manual(
    name = "Connectivity type", values = colors, 
    breaks = c(levels(dat[[1]])[[1]]$name[-1], "Background"),
    na.translate = FALSE) +
  facet_wrap(~lyr, ncol = 4, nrow = 1, strip.position = "left") +
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

# ggsave(file.path(fig_dir, "Figure5_conn_types.png"),
#        width = 6.5, height = 3.2, dpi = 300, bg = "white")

ggarrange(g1, g2, ncol = 1, nrow = 2, heights = c(5.3, 3.2))

ggsave(file.path(fig_dir, "Figure4_conn_all.png"),
       width = 6.5, height = 8.5, dpi = 300, bg = "white")

dat <- subset(lyrs, c(nms[str_detect(nms, "ssp370")], 
                      nms[str_detect(nms, "ssp585")]))
levels(dat) <- lapply(1:nlyr(dat), function(x) classes[, c("ID", "name")])
names(dat) <- letters[1:nlyr(dat)]

ggplot() +
  geom_spatraster(data = dat, na.rm = TRUE) +
  scale_fill_manual(
    name = "Connectivity type", values = colors, 
    breaks = c(levels(dat[[1]])[[1]]$name[-1], "Background"),
    na.translate = FALSE) +
  facet_wrap(~lyr, ncol = 4, nrow = 2, strip.position = "left") +
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

ggsave(file.path(fig_dir, "Figure_s_conn_types_simple.png"),
       width = 6.5, height = 5.5, dpi = 300, bg = "white")

#### Figure S1 ####
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

#### Figure S2 ####
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

#### Figure S3 ####
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
