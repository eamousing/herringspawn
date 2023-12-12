#!/usr/bin/env Rscript

# Analysis script: Herring spawning migration and bioenergetics 2016 year-class
#
# Author: Erik Askov Mousing, IMR, erik.askov.mousing@hi.no
# Date: 2023-10-20
#
# Data dependencies:
#
# Output:
#   The script returns the results from the migration/bioenergetics model
#   as csv-files per year and a single plot.
#
# To run the script from the terminal on UNIX/Linux:
#
#   > Rscript herringspawn_main.R

# Suppress warnings loading packages
suppress_warning <- TRUE
if (suppress_warning) {
  options(warn=-1)
  options(dplyr.summarise.inform = FALSE)
  suppressPackageStartupMessages(library("tidyverse"))
  suppressPackageStartupMessages(library("ggplot2"))
  suppressPackageStartupMessages(library("patchwork"))
  suppressPackageStartupMessages(library("maps"))
  suppressPackageStartupMessages(library("geosphere"))
  suppressPackageStartupMessages(library("xlsx"))
  suppressPackageStartupMessages(library("R.matlab"))  
}

# Load additional libraries and functions
library(tidyverse)
library(ggplot2)
library(patchwork)
library(maps)
library(geosphere)
library(xlsx)
library(R.matlab)
source("herringspawn_functions.R")

# Write out some setup information
print("Analysis script: Herring spawning migration and bioenergetics 2016 year-class")
print("Author: Erik Askov Mousing, Institute of Marine Research, 2023")
print("")
print("Analysis setup:")
print(paste('Git branch:', system("git rev-parse --short HEAD", intern = TRUE)))
print("")
print("Software versions:")
print(paste(R.version.string))
print(paste('tidyverse version', packageVersion('tidyverse')))
print(paste('ggplot2 version', packageVersion('ggplot2')))
print(paste('patchwork version', packageVersion('patchwork')))
print(paste('maps version', packageVersion('maps')))
print(paste('geosphere version', packageVersion('geosphere')))
print(paste('xlsx version', packageVersion('xlsx')))
print(paste('R.matlab version', packageVersion('R.matlab')))
print("")

#### Setup environment and prepare simulations ####

print("1/5: Reading data and setting up simulations...")

# Setup directory paths
data_dir <- "./data/"
results_dir <- "./results/"
# Make results directory if non-existent
system("mkdir -p results")

# Read migration route and currents data
her_cog <- get_route_positions(data_dir)
her_cog <- her_cog[2:13,]
her_cur <- get_currents(data_dir, her_cog)

print(her_cog)

# Reduce by 25% to account for deeper position of fish
her_cur$u <- her_cur$u*0.75
her_cur$v <- her_cur$v*0.75

# Total distance
tot_dist <- 0
trans_dist <- data.frame(trans = vector(length = 11), dist = vector(length = 11))
for (i in 1:11) {
  t_dist <- distm(c(her_cog$lon[i], her_cog$lat[i]), c(her_cog$lon[i+1], her_cog$lat[i+1]))
  tot_dist <- tot_dist + t_dist
  trans_dist$trans[i] = i
  trans_dist$dist[i] = tot_dist / 1000
}
tot_dist <- tot_dist/1000
trans_dist$dist_rel <- trans_dist$dist / max(trans_dist$dist) * 100

# Remove 1994 due to missing data
her_cur <- her_cur[her_cur$year != 1994, ]

#### Define simulation setup ####
temp <- 8 # Water temperature
years <- unique(her_cur$year) # Years to simulate
transects <- unique(her_cur$trans) # Transects to swim

# Bioenergic parameters 
bio_params <- list(
  "ra" = 0.0033, # intercept for max standard resp (g o2 g-1 d-1)
  "rb" = -0.227, # slope for maximum standard respiration
  "rq" = 0.0548, # Temperature dependence of temp on resting metabolism
  "rt0" = 0.03, # Temperature impact on activity
  "oxy_coef" = 20083.0,
  "sspeed" = 1.0, # relative swimming speed (bl sec-1)
  "fed" = 11000, # Fish energy density (J g-1)
  "wc" = 0.68, # Water content (%)
  "ed_fat" = 41, # Energy density fat (kJ g-1)
  "ed_gonad" = 8, # Energy density gonads (kJ g-1)
  "ed_solid" = 21, # Energy density solids (kJ g-1)
  "fat_frac" = 0.20, # Initial fat fraction
  "gonad_frac" = 0.10, # Initial gonad fraction
  "solid_frac" = 0.70 # Initial solid fraction
)

#### Simulation ####
print("2/5: Running the simulations...")

# Simulation 1: Observed fish size in 2020
fish_len <- 27.5 # Fish length (cm)
fish_wei <- 168 # 180 # Fish weight (grams)
fish_dwei <- fish_wei * (1 - bio_params$wc)
fish_energy <- calc_tot_energy(fish_dwei) # Fish total energy content
sim1 <- simulate_migration(fish_len, fish_wei, 2020, transects, bio_params)
sim1$energy_rel <- sim1$energy_accum / fish_energy
sim1_end <- sim1[sim1$trans == max(transects), ]

# Simulation 2: Observed fish size in 2021
fish_len <- 29 # Fish length (cm)
fish_wei <- 210 # 205 # Fish weight (grams)
fish_dwei <- fish_wei * (1 - bio_params$wc) # Fish dry weight
fish_energy <- calc_tot_energy(fish_dwei) # Fish total energy content
sim2 <- simulate_migration(fish_len, fish_wei, 2021, transects, bio_params)
sim2$energy_rel <- sim2$energy_accum / fish_energy
sim2_end <- sim2[sim2$trans == max(transects), ]

# Simulation 3: Observed fish size in 2022
fish_len <- 30 # Fish length (cm)
fish_wei <- 242 # Fish weight (grams)
fish_dwei <- fish_wei * (1 - bio_params$wc) # Fish dry weight
fish_energy <- calc_tot_energy(fish_dwei) # Fish total energy content
sim3 <- simulate_migration(fish_len, fish_wei, 2022, transects, bio_params)
sim3$energy_rel <- sim3$energy_accum / fish_energy
sim3_end <- sim3[sim3$trans == max(transects), ]

# Simulation 4: Observed fish size in 2023
fish_len <- 31.5 # Fish length (cm)
fish_wei <- 284 # Fish weight (grams)
fish_dwei <- fish_wei * (1 - bio_params$wc) # Fish dry weight
fish_energy <- calc_tot_energy(fish_dwei) # Fish total energy content
sim4 <- simulate_migration(fish_len, fish_wei, 2023, transects, bio_params)
sim4$energy_rel <- sim4$energy_accum / fish_energy
sim4_end <- sim4[sim4$trans == max(transects), ]

# Write model output to files
print("3/5: Writing output to fil...")
write.csv2(sim1, paste0(results_dir, "mod_results_observed_2020.csv"), row.names = F)
write.csv2(sim2, paste0(results_dir, "mod_results_observed_2021.csv"), row.names = F)
write.csv2(sim3, paste0(results_dir, "mod_results_observed_2022.csv"), row.names = F)
write.csv2(sim4, paste0(results_dir, "mod_results_observed_2023.csv"), row.names = F)  



# Concatenate results from observed fish sizes
sim_obs <- rbind(sim1, sim2, sim3, sim4)
sim_obs$energy_rel <- sim_obs$energy_rel*100 # Convert to %

#### Plots ####
print("4/5: Making plot...")

# Prepare data for plotting
world <- map_data("world")

# Read spawning location data
her_df <- read.xlsx2(paste0(data_dir, "LUF11_2018_AllVessels_Herring_NASC.xlsx"), 1)
names(her_df) <- c("date", "start_log", "stop_log", "lat", "lon", "NASC")
her_df$date <- as.Date(her_df$date, format = "%Y%m%d")
her_df$lon <- as.numeric(her_df$lon)
her_df$lat <- as.numeric(her_df$lat)
her_df$NASC <- as.numeric(her_df$NASC)

# Read currents data for 2020
jan_feb <- readMat(paste0(data_dir, "currents1993_2023jan_feb.mat"))
yr_index <- which(jan_feb$year == 2020)
uv_2020 <- data.frame(
  lon = as.vector(jan_feb$xlon),
  lat = as.vector(jan_feb$ylat),
  u = as.vector(jan_feb$ueast[,,yr_index]),
  v = as.vector(jan_feb$vnorth[,,yr_index])
)
uv_2020$w <- sqrt(uv_2020$u^2 + uv_2020$v^2)

# Calculate the maximum number of transects before the K threshold
k1 <- round(get_k_trans(sim1, 0.65), 2)
k2 <- round(get_k_trans(sim2, 0.65), 2)
k3 <- round(get_k_trans(sim3, 0.65), 2)
k4 <- round(get_k_trans(sim4, 0.65), 2)

# Estimate position along transect
loc1 <- estimate_loc(her_cog, k1)
loc2 <- estimate_loc(her_cog, k2)
loc3 <- estimate_loc(her_cog, k3)
loc4 <- estimate_loc(her_cog, k4)



# Create data frames for plotting transects swimmed before the K threshold
dist_obs <- data.frame(
  year = c("2020", "2021", "2022", "2023"),
  lon = c(loc1[1], loc2[1], loc3[1], loc4[1]),
  lat = c(loc1[2], loc2[2], loc3[2], loc4[2])
)

# Create plots

# Map
p_dist <- ggplot() +
  geom_map(data = world, map = world, aes(long, lat, map_id=region), fill = "lightgreen", color="darkgreen") +
  geom_point(data = her_df, aes(x = lon, y = lat), col = "wheat", alpha = 0.1) +
  geom_segment(
    aes(x = her_cog$lon[1:14], 
        y = her_cog$lat[1:14],
        xend = her_cog$lon[2:15],
        yend = her_cog$lat[2:15]),
    arrow = arrow(length=unit(.3, 'cm'), type = "open"),
    color = "darkblue") +
  geom_point(aes(x = lon, y = lat), her_cog, color = "darkblue") +
  geom_point(aes(x = lon, y = lat, col = as.factor(year)), dist_obs, shape = 16, size = 4, show.legend = FALSE) +
  geom_point(aes(x = lon, y = lat), dist_obs, shape = 21, size = 4, show.legend = FALSE) +
  geom_text(aes(x = lon, y = lat, label = as.factor(year)), dist_obs, hjust = 1.5, size = 3) +
  scale_y_continuous(limits = c(62,71), expand = c(0, 0)) +
  scale_x_continuous(limits = c(4,18), expand = c(0, 0)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

# Currents
p_cur2020 <- ggplot() +
  geom_map(data = world, map = world, aes(long, lat, map_id=region), fill = "lightgreen", color="darkgreen") +
  geom_segment(data = uv_2020, aes(x = lon, xend = lon + u, y = lat, yend = lat + v, colour = w), 
               arrow = arrow(angle = 15, length = unit(0.03, "inches"), type = "closed"), alpha = 0.3, show.legend = F) +
  scale_color_gradient(low = "black", high = "red") +
  scale_y_continuous(limits = c(62,71), expand = c(0, 0)) +
  scale_x_continuous(limits = c(4,18), expand = c(0, 0)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Condition after spawning
p_cond_k <- ggplot(sim_obs, aes(x = trans)) +
  geom_hline(yintercept = 0.65, linetype="dashed") +
  scale_x_continuous(breaks=c(1:11), labels=c(1:11),limits=c(1,11)) +
  geom_line(aes(y = cond_k, col = as.factor(year)), show.legend = F) +
  xlab("Transect") +
  ylab("CF after spawning") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Swimming speed
p_sspeed <- ggplot(sim_obs, aes(x = trans)) +
  geom_line(aes(y = sspeed, col = as.factor(year)), show.legend = F) +
  scale_x_continuous(breaks=c(1:11), labels=c(1:11),limits=c(1,11)) +
  xlab("Transect") +
  ylab("Swimming speed (bl/s)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Weight
p_wei <- ggplot(sim_obs, aes(x = trans)) +
  geom_line(aes(y = wei, col = as.factor(year)), show.legend = F) +
  scale_x_continuous(breaks=c(1:11), labels=c(1:11),limits=c(1,11)) +
  xlab("Transect") +
  ylab("Weight (g)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Energy used
p_rel_energy <- ggplot(sim_obs, aes(x =trans)) +
  geom_line(aes(y = energy_rel, col = as.factor(year)), show.legend = F) +
  scale_x_continuous(breaks=c(1:11), labels=c(1:11),limits=c(1,11)) +
  xlab("Transect") +
  ylab("Energy use (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Create custom layout
my_layout <- "
ABCD
ABEF
"
# Combine plots
p_dist + p_cur2020 + p_wei + p_rel_energy + p_cond_k + p_sspeed + 
  plot_layout(guides = "collect", design = my_layout) & 
  plot_annotation(tag_levels = "a") & theme(legend.position = 'bottom', legend.title = element_blank())
ggsave(paste0(results_dir, "migration.tiff"), scale = 1.1, dpi = 300, width = 10, height = 4.5, compression = "lzw")

print("5/5: Analysis completed!")
