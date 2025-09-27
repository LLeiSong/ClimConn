library(readxl)
root_dir <- "/home/lsong/ClimConn"

# Get species list
species_list <- list.files(file.path(root_dir, "data/mamiferos"))
species_list <- gsub("\\.", " ", species_list) %>% 
  data.frame(species = .)

# Load species trait
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.gd0m3
traits_mammals <- read_excel(
  file.path(root_dir, "data/dispersal", 
            "Generation\ Lenght\ for\ Mammals.xlsx"), sheet = 1)

# Calculate the natal dispersal distance
## https://www.pnas.org/doi/epdf/10.1073/pnas.1116791109
traits_mammals <- traits_mammals %>% 
  mutate(AdultBodyMass_kg = AdultBodyMass_g / 1000) %>% 
  select(Order, Scientific_name, AdultBodyMass_kg, GenerationLength_d) %>% 
  mutate(dispersal_distance = ifelse(Order == "Carnivora",
                                     3.45 * AdultBodyMass_kg^0.89,
                                     1.45 * AdultBodyMass_kg^0.54)) %>% 
  mutate(dispersal_rate = dispersal_distance / (GenerationLength_d / 365))

# Attach the values to species
species_dispersal <- left_join(
  species_list, traits_mammals, 
  by = c('species' = 'Scientific_name'))

# Deal with missing species for dispersal distance
## Using the median of the order
order_dispersal <- traits_mammals %>% 
  filter(!is.na(Order)) %>% 
  group_by(Order) %>% 
  summarise(dispersal_rate = median(dispersal_rate, na.rm = TRUE))

# Concatenate everything together
class_search <- classification(species_left$species, db = 'gbif')
class_search <- lapply(species_left$species, function(sp){
  ord <- class_search[[sp]]
  ord <- ord[ord$rank == "order", "name"]
  data.frame(species = sp, order = ord)
}) %>% bind_rows()

# Replace two missing orders with similar order
class_search <- class_search %>% 
  mutate(order = ifelse(order == "Soricomorpha", "Rodentia", order)) %>% 
  mutate(order = ifelse(order == "Artiodactyla", "Perissodactyla", order))

species_left <- species_dispersal %>% filter(is.na(Order))
species_left <- left_join(species_left, class_search, by = "species") %>% 
  select(species, order) %>% 
  left_join(order_dispersal, by = c("order" = "Order"))

species_dispersal <- species_dispersal %>% 
  filter(!is.na(Order)) %>% 
  select(species, Order, dispersal_rate) %>% 
  rename(order = Order) %>% 
  rbind(species_left)

write.csv(species_dispersal, 
          file.path(root_dir, "data/dispersal/species_dispersal_rate.csv"),
          row.names = FALSE)
