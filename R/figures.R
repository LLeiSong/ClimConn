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

#### Figure 4, S4 and S5 ####
# Figure 4
typology_dir <- file.path(root_dir, "results/typology")
fnames <- list.files(typology_dir, pattern = ".tif$", full.names = TRUE)
lyrs <- rast(fnames)
names(lyrs) <- str_extract(basename(fnames), "[0-9]{4}-[0-9]{4}_ssp[0-9]{3}")

col <- st_read(file.path(data_dir, "colombia.geojson")) %>% 
  st_transform(st_crs(lyrs))

lyrs <- mask(crop(lyrs, col), col)

conns <- rast(file.path(result_dir, "conn_lyrs.tif")) %>% 
  crop(col) %>% mask(col)

colors <- c("#f5dfbb", "#9bc1bc", "#127475", "#333d29",
            "#f6ae2d", "#f26419", "#582f0e", "#4361ee")

# Make data for transition
years <- c("2021-2040", "2041-2060", "2061-2080", "2081-2100")
ssps <- c("ssp126", "ssp370", "ssp585")

trans_df <- lapply(ssps, function(ssp){
  lyrs_ssp <- subset(lyrs, str_detect(names(lyrs), ssp))
  
  vals <- values(c(cellSize(lyrs_ssp[[1]], unit = "km"), lyrs_ssp))
  vals <- vals[complete.cases(vals), , drop = FALSE] %>% 
    as_tibble()
  names(vals) <- c("area", years)
  
  vals %>% 
    count(`2021-2040`, `2041-2060`, `2061-2080`, `2081-2100`, wt = area) %>% 
    make_long(`2021-2040`, `2041-2060`, `2061-2080`, `2081-2100`, value = n) %>% 
    mutate(ssp = ssp)
}) %>% bind_rows() %>% data.frame() %>% 
  mutate(node = factor(
    node, levels = levels(lyrs[[1]])[[1]][["ID"]],
    labels = levels(lyrs[[1]])[[1]][["class"]]),
    next_node = factor(
      next_node, levels = levels(lyrs[[1]])[[1]][["ID"]],
      labels = levels(lyrs[[1]])[[1]][["class"]]))

# Make plot
g1 <- ggplot() +
  geom_spatraster(data = conns$`2021-2040_ssp126`, na.rm = TRUE) +
  scale_fill_grass_c("Conn.          ", palette = "inferno", 
                     na.value = "transparent") +
  geom_sf(data = col, color = "black", fill = "transparent", 
          linewidth = 0.6) +
  theme_void() +
  theme(legend.position = c(1.4, 0.55),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.width = unit(0.3, "cm"))

g2 <- ggplot() +
  geom_spatraster(data = lyrs$`2021-2040_ssp126`, na.rm = TRUE) +
  scale_fill_manual(
    name = "Connectivity type", values = colors, na.translate = FALSE) +
  geom_sf(data = col, color = "black", 
          fill = "transparent", linewidth = 0.6) +
  theme_void() +
  theme(legend.position = "none",
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11))

g1 <- g2 + inset_element(g1, left = 0.55, bottom = 0.62, right = 1, top = 1)

dat <- trans_df %>% filter(ssp == "ssp126")

g2 <- ggplot(dat, aes(
  x = x,
  next_x = next_x,
  node = node,
  next_node = next_node,
  fill = node,
  value = value
)) + xlab("") +
  geom_sankey(flow.alpha = 0.7) +
  scale_fill_manual(name = "Type", values = colors) +
  theme_sankey() +
  theme(text = element_text(size = 11),
        legend.position = "bottom",
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11),
        axis.text.x = element_text(size = 11, angle = 90))

ggarrange(ggarrange(g1, NULL, g2 + theme(legend.position = "none"), ncol = 3,
          widths = c(1, 0.1, 1),
          labels = c("a", "", "b")),
          as_ggplot(get_legend(g2)), nrow = 2, heights = c(6, 1))

ggsave(file.path(fig_dir, "Figure4_conn_types.png"),
       width = 6.5, height = 5, dpi = 300, bg = "white")

# Fig. S4
figs <- lapply(names(lyrs)[-1], function(nm){
  g1 <- ggplot() +
    geom_spatraster(data = conns[[nm]], na.rm = TRUE) +
    scale_fill_grass_c("Conn.", palette = "inferno", 
                       na.value = "transparent") +
    geom_sf(data = col, color = "black", fill = "transparent", 
            linewidth = 0.4) +
    theme_void() +
    theme(legend.position = c(1.35, 0.6),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.key.width = unit(0.2, "cm"))
  
  ggplot() +
    geom_spatraster(data = lyrs[[nm]], na.rm = TRUE) +
    scale_fill_manual(
      name = "Connectivity type", values = colors, na.translate = FALSE) +
    geom_sf(data = col, color = "black", 
            fill = "transparent", linewidth = 0.4) +
    ggtitle(nm) +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(
            hjust = -1.3, vjust = -68, size = 16, face = "bold")) +
    inset_element(g1, left = 0.55, bottom = 0.62, right = 1, top = 1)
})

ggarrange(
  plotlist = c(figs, list(as_ggplot(
    get_legend(g2 + theme(
      legend.position = "right",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 16)))))), ncol = 4, nrow = 3)

ggsave(file.path(fig_dir, "Figure_s_conn_types.png"),
       width = 13, height = 12, dpi = 300, bg = "white")

# Fig. S5
figs <- lapply(ssps[-1], function(scenario){
  dat <- trans_df %>% filter(ssp == scenario)
  
  ggplot(dat, aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = node,
    value = value
  )) + xlab("") + ggtitle(scenario) +
    geom_sankey(flow.alpha = 0.7) +
    scale_fill_manual(name = "Type", values = colors) +
    theme_sankey() +
    theme(text = element_text(size = 11),
          legend.position = "bottom",
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 11),
          axis.text.x = element_text(size = 11, angle = 90),
          plot.title = element_text(
            hjust = 0.5, size = 16, face = "bold"))
})

ggarrange(plotlist = figs, ncol = 2, nrow = 1,
          common.legend = TRUE, legend = "bottom")

ggsave(file.path(fig_dir, "Figure_s_transition.png"),
       width = 6.5, height = 5, dpi = 300, bg = "white")

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
