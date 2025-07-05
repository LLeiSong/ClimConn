library(ggplot2)
library(tidyterra)
library(ggsci)
library(here)
library(ggpubr)
library(smoothr)
library(stars)
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

# Fig. 3
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
  theme_pubclean(base_size = 12) + 
  theme(axis.text = element_text(color = "black"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

# Fig. S1
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

# Fig. S2
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

g1 <- ggplot() +
  geom_spatraster(data = richness) +
  scale_fill_grass_c(
    name = "No of species", 
    palette = "plasma", na.value = "transparent") +
  geom_sf(data = col, fill = "transparent") +
  theme_void() +
  guides(fill = guide_colorbar(title.vjust = 0.8)) +
  theme(legend.position = "bottom",
        plot.margin = unit(c(0, -1, 0, -1), "cm"))

g2 <- ggplot() +
  geom_spatraster(data = suitability) +
  scale_fill_grass_c(
    name = "Suitability", 
    palette = "viridis", na.value = "transparent") +
  geom_sf(data = col, fill = "transparent") +
  theme_void() +
  guides(fill = guide_colorbar(title.vjust = 0.8)) +
  theme(legend.position = "bottom",
        plot.margin = unit(c(0, -1, 0, -1), "cm"))

ggarrange(g1, g2, ncol = 2, labels = c("a", "b"))

ggsave(file.path(fig_dir, "Figure_s_species.png"),
       width = 6, height = 4.4, dpi = 500, bg = "white")


# Fig. S3
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
