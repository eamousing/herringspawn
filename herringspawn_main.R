

#### Setup environment and prepare simulations ####
# Setup directory paths
data_dir <- "./data/"
results_dir <- "./results/"

# Load additional libraries and functions
library(tidyverse)
library(ggplot2)
library(patchwork)
library(maps)
library(geosphere)
require(xlsx)
source("herringspawn_functions.R")

# Read migration route and currents data
her_cog <- get_route_positions(data_dir)
her_cog <- her_cog[2:13,]
her_cur <- get_currents(data_dir, her_cog)

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
transects <- unique(her_cur$trans)

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

# Simulation 1: Generic fish (k=0.9, 27 cm)
fish_len <- 27 # Fish length (cm)
fish_wei <- 177 # Fish weight (grams)
fish_dwei <- fish_wei * (1 - bio_params$wc) # Fish dry weight
fish_energy <- calc_tot_energy(fish_dwei) # Fish total energy content
sim1r <- simulate_migration(fish_len, fish_wei, years, transects, bio_params)
sim1r$energy_rel <- sim1r$energy_accum / fish_energy
sim1_2020 <- sim1r[sim1r$year == 2020, ]
sim1_2021 <- sim1r[sim1r$year == 2021, ]
sim1_2022 <- sim1r[sim1r$year == 2022, ]
sim1_2023 <- sim1r[sim1r$year == 2023, ]
sim1_end <- sim1r[sim1r$trans == max(transects), ]

# Simulation 2: Generic fish (k=0.9, 30 cm)
fish_len <- 30 # Fish length (cm)
fish_wei <- 243 # Fish weight (grams)
fish_dwei <- fish_wei * (1 - bio_params$wc) # Fish dry weight
fish_energy <- calc_tot_energy(fish_dwei) # Fish total energy content
sim2r <- simulate_migration(fish_len, fish_wei, years, transects, bio_params)
sim2r$energy_rel <- sim2r$energy_accum / fish_energy
sim2_2020 <- sim2r[sim2r$year == 2020, ]
sim2_2021 <- sim2r[sim2r$year == 2021, ]
sim2_2022 <- sim2r[sim2r$year == 2022, ]
sim2_2023 <- sim2r[sim2r$year == 2023, ]
sim2_end <- sim2r[sim2r$trans == max(transects), ]

# Simulation 3: Generic fish (k=0.9, 33 cm)
fish_len <- 33 # Fish length (cm)
fish_wei <- 323 # Fish weight (grams)
fish_dwei <- fish_wei * (1 - bio_params$wc) # Fish dry weight
fish_energy <- calc_tot_energy(fish_dwei) # Fish total energy content
sim3r <- simulate_migration(fish_len, fish_wei, years, transects, bio_params)
sim3r$energy_rel <- sim3r$energy_accum / fish_energy
sim3_2020 <- sim3r[sim3r$year == 2020, ]
sim3_2021 <- sim3r[sim3r$year == 2021, ]
sim3_2022 <- sim3r[sim3r$year == 2022, ]
sim3_2023 <- sim3r[sim3r$year == 2023, ]
sim3_end <- sim3r[sim3r$trans == max(transects), ]

# Simulation 4: Generic fish (k=0.9, 36 cm)
fish_len <- 36 # Fish length (cm)
fish_wei <- 419 # Fish weight (grams)
fish_dwei <- fish_wei * (1 - bio_params$wc) # Fish dry weight
fish_energy <- calc_tot_energy(fish_dwei) # Fish total energy content
sim4r <- simulate_migration(fish_len, fish_wei, years, transects, bio_params)
sim4r$energy_rel <- sim4r$energy_accum / fish_energy
sim4_2020 <- sim4r[sim4r$year == 2020, ]
sim4_2021 <- sim4r[sim4r$year == 2021, ]
sim4_2022 <- sim4r[sim4r$year == 2022, ]
sim4_2023 <- sim4r[sim4r$year == 2023, ]
sim4_end <- sim4r[sim4r$trans == max(transects), ]

# Simulation 5: Observed fish size in 2020
fish_len <- 27.5 # Fish length (cm)
fish_wei <- 168 # 180 # Fish weight (grams)
fish_dwei <- fish_wei * (1 - bio_params$wc)
fish_energy <- calc_tot_energy(fish_dwei) # Fish total energy content
sim5r <- simulate_migration(fish_len, fish_wei, 2020, transects, bio_params)
sim5r$energy_rel <- sim5r$energy_accum / fish_energy
sim5_end <- sim5r[sim5r$trans == max(transects), ]

# Simulation 6: Observed fish size in 2021
fish_len <- 29 # Fish length (cm)
fish_wei <- 210 # 205 # Fish weight (grams)
fish_dwei <- fish_wei * (1 - bio_params$wc) # Fish dry weight
fish_energy <- calc_tot_energy(fish_dwei) # Fish total energy content
sim6r <- simulate_migration(fish_len, fish_wei, 2021, transects, bio_params)
sim6r$energy_rel <- sim6r$energy_accum / fish_energy
sim6_end <- sim6r[sim6r$trans == max(transects), ]

# Simulation 7: Observed fish size in 2022
fish_len <- 30 # Fish length (cm)
fish_wei <- 242 # Fish weight (grams)
fish_dwei <- fish_wei * (1 - bio_params$wc) # Fish dry weight
fish_energy <- calc_tot_energy(fish_dwei) # Fish total energy content
sim7r <- simulate_migration(fish_len, fish_wei, 2022, transects, bio_params)
sim7r$energy_rel <- sim7r$energy_accum / fish_energy
sim7_end <- sim7r[sim7r$trans == max(transects), ]

# Simulation 8: Observed fish size in 2023
fish_len <- 31.5 # Fish length (cm)
fish_wei <- 284 # Fish weight (grams)
fish_dwei <- fish_wei * (1 - bio_params$wc) # Fish dry weight
fish_energy <- calc_tot_energy(fish_dwei) # Fish total energy content
sim8r <- simulate_migration(fish_len, fish_wei, 2023, transects, bio_params)
sim8r$energy_rel <- sim8r$energy_accum / fish_energy
sim8_end <- sim8r[sim8r$trans == max(transects), ]

# Concatenate results from observed fish sizes
sim_obs_end <- rbind(sim5_end, sim6_end, sim7_end, sim8_end)

# Concatenate results from all simulations
sim1_end$sim <- "sim1"
sim2_end$sim <- "sim2"
sim3_end$sim <- "sim3" 
sim4_end$sim <- "sim4"
sim_obs_end$sim <- "sim999"

sims <- bind_rows(sim1_end, sim2_end, sim3_end, sim4_end)


#### Plots ####
world <- map_data("world")

# Read spawning location data
her_df <- read.xlsx2(paste0(data_dir, "LUF11_2018_AllVessels_Herring_NASC.xlsx"), 1)
names(her_df) <- c("date", "start_log", "stop_log", "lat", "lon", "NASC")
her_df$date <- as.Date(her_df$date, format = "%Y%m%d")
her_df$lon <- as.numeric(her_df$lon)
her_df$lat <- as.numeric(her_df$lat)
her_df$NASC <- as.numeric(her_df$NASC)

# Calculate the maximum number of transects before the K threshold
dist1_2020 <- round(get_k_trans(sim1_2020, 0.65), 0)
dist2_2020 <- round(get_k_trans(sim2_2020, 0.65), 0)
dist3_2020 <- round(get_k_trans(sim3_2020, 0.65), 0)
dist4_2020 <- round(get_k_trans(sim4_2020, 0.65), 0)
dist5_2020 <- round(get_k_trans(sim5r, 0.65), 0)
dist1_2021 <- round(get_k_trans(sim1_2021, 0.65), 0)
dist2_2021 <- round(get_k_trans(sim2_2021, 0.65), 0)
dist3_2021 <- round(get_k_trans(sim3_2021, 0.65), 0)
dist4_2021 <- round(get_k_trans(sim4_2021, 0.65), 0)
dist5_2021 <- round(get_k_trans(sim6r, 0.65), 0)
dist1_2022 <- round(get_k_trans(sim1_2022, 0.65), 0)
dist2_2022 <- round(get_k_trans(sim2_2022, 0.65), 0)
dist3_2022 <- round(get_k_trans(sim3_2022, 0.65), 0)
dist4_2022 <- round(get_k_trans(sim4_2022, 0.65), 0)
dist5_2022 <- round(get_k_trans(sim7r, 0.65), 0)
dist1_2023 <- round(get_k_trans(sim1_2023, 0.65), 0)
dist2_2023 <- round(get_k_trans(sim2_2023, 0.65), 0)
dist3_2023 <- round(get_k_trans(sim3_2023, 0.65), 0)
dist4_2023 <- round(get_k_trans(sim4_2023, 0.65), 0)
dist5_2023 <- round(get_k_trans(sim8r, 0.65), 0)

# Create data frames for plotting
dist_2020 <- data.frame(
  len = c("27 cm", "30 cm", "33 cm", "36 cm", "2016 year-class"),
  lon = c(her_cog$lon[dist1_2020+1], her_cog$lon[dist2_2020+1]-0.3, her_cog$lon[dist3_2020+1]+0.3, her_cog$lon[dist4_2020+1], her_cog$lon[dist5_2020+1]),
  lat = c(her_cog$lat[dist1_2020+1], her_cog$lat[dist2_2020+1], her_cog$lat[dist3_2020+1], her_cog$lat[dist4_2020+1], her_cog$lat[dist5_2020+1])
)

dist_2021 <- data.frame(
  len = c("27 cm", "30 cm", "33 cm", "36 cm", "2016 year-class"),
  lon = c(her_cog$lon[dist1_2021+1], her_cog$lon[dist2_2021+1]-0.3, her_cog$lon[dist3_2021+1]+0.3, her_cog$lon[dist4_2021+1], her_cog$lon[dist5_2021+1]),
  lat = c(her_cog$lat[dist1_2021+1], her_cog$lat[dist2_2021+1], her_cog$lat[dist3_2021+1], her_cog$lat[dist4_2021+1], her_cog$lat[dist5_2021+1])
)

dist_2022 <- data.frame(
  len = c("27 cm", "30 cm", "33 cm", "36 cm", "2016 year-class"),
  lon = c(her_cog$lon[dist1_2022+1], her_cog$lon[dist2_2022+1]-0.3, her_cog$lon[dist3_2022+1]-0.3, her_cog$lon[dist4_2022+1]+0.3, her_cog$lon[dist5_2022+1]+0.3),
  lat = c(her_cog$lat[dist1_2022+1], her_cog$lat[dist2_2022+1], her_cog$lat[dist3_2022+1], her_cog$lat[dist4_2022+1], her_cog$lat[dist5_2022+1])
)

dist_2023 <- data.frame(
  len = c("27 cm", "30 cm", "33 cm", "36 cm", "2016 year-class"),
  lon = c(her_cog$lon[dist1_2023+1], her_cog$lon[dist2_2023+1], her_cog$lon[dist3_2023+1]-0.3, her_cog$lon[dist4_2023+1], her_cog$lon[dist5_2023+1]+0.3),
  lat = c(her_cog$lat[dist1_2023+1], her_cog$lat[dist2_2023+1], her_cog$lat[dist3_2023+1], her_cog$lat[dist4_2023+1], her_cog$lat[dist5_2023+1])
)

# Create maps

custom_theme <- theme(axis.text.x=element_blank(), #remove x axis labels
                      axis.ticks.x=element_blank(), #remove x axis ticks
                      axis.text.y=element_blank(),  #remove y axis labels
                      axis.ticks.y=element_blank(),  #remove y axis ticks
                      axis.title.x=element_blank(), 
                      axis.title.y=element_blank()
)

p0_2020 <- ggplot() +
  geom_map(data = world, map = world, aes(long, lat, map_id=region), fill = "grey", color="darkgrey") +
  geom_point(data = her_df, aes(x = lon, y = lat), col = "lightblue", alpha = 0.1) +
  geom_segment(
    aes(x = her_cog$lon[1:14], 
        y = her_cog$lat[1:14],
        xend = her_cog$lon[2:15],
        yend = her_cog$lat[2:15]),
    arrow = arrow(length=unit(.3, 'cm'), type = "open"),
    color = "darkblue") +
  geom_point(aes(x = lon, y = lat), her_cog, color = "darkblue") +
  geom_point(aes(x = lon, y = lat), dist_2020, color = c(2,3,4,5,1), shape = 16, size = 3) +
  scale_y_continuous(limits = c(61,71), expand = c(0, 0)) +
  scale_x_continuous(limits = c(4,18), expand = c(0, 0)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()
# p0_2020

p0_2021 <- ggplot() +
  geom_map(data = world, map = world, aes(long, lat, map_id=region), fill = "grey", color="darkgrey") +
  geom_point(data = her_df, aes(x = lon, y = lat), col = "lightblue", alpha = 0.1) +
  geom_segment(
    aes(x = her_cog$lon[1:14], 
        y = her_cog$lat[1:14],
        xend = her_cog$lon[2:15],
        yend = her_cog$lat[2:15]),
    arrow = arrow(length=unit(.3, 'cm'), type = "open"),
    color = "darkblue") +
  geom_point(aes(x = lon, y = lat), her_cog, color = "darkblue") +
  geom_point(aes(x = lon, y = lat), dist_2021, color = c(2,3,4,5,1), shape = 16, size = 3) +
  scale_y_continuous(limits = c(61,71), expand = c(0, 0)) +
  scale_x_continuous(limits = c(4,18), expand = c(0, 0)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()
# p0_2021

p0_2022 <- ggplot() +
  geom_map(data = world, map = world, aes(long, lat, map_id=region), fill = "grey", color="darkgrey") +
  geom_point(data = her_df, aes(x = lon, y = lat), col = "lightblue", alpha = 0.1) +
  geom_segment(
    aes(x = her_cog$lon[1:14], 
        y = her_cog$lat[1:14],
        xend = her_cog$lon[2:15],
        yend = her_cog$lat[2:15]),
    arrow = arrow(length=unit(.3, 'cm'), type = "open"),
    color = "darkblue") +
  geom_point(aes(x = lon, y = lat), her_cog, color = "darkblue") +
  geom_point(aes(x = lon, y = lat), dist_2022, color = c(2,3,4,5,1), shape = 16, size = 3) +
  scale_y_continuous(limits = c(61,71), expand = c(0, 0)) +
  scale_x_continuous(limits = c(4,18), expand = c(0, 0)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()
# p0_2022

p0_2023 <- ggplot() +
  geom_map(data = world, map = world, aes(long, lat, map_id=region), fill = "grey", color="darkgrey") +
  geom_point(data = her_df, aes(x = lon, y = lat), col = "lightblue", alpha = 0.1) +
  geom_segment(
    aes(x = her_cog$lon[1:14], 
        y = her_cog$lat[1:14],
        xend = her_cog$lon[2:15],
        yend = her_cog$lat[2:15]),
    arrow = arrow(length=unit(.3, 'cm'), type = "open"),
    color = "darkblue") +
  geom_point(aes(x = lon, y = lat), her_cog, color = "darkblue") +
  geom_point(aes(x = lon, y = lat), dist_2023, color = c(2,3,4,5,1), shape = 16, size = 3) +
  scale_y_continuous(limits = c(61,71), expand = c(0, 0)) +
  scale_x_continuous(limits = c(4,18), expand = c(0, 0)) +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()
# p0_2023

# Calculate v anomaly
current <- data.frame(
  year = years,
  v_mean = tapply(her_cur$v, her_cur$year, median)
)
current$v_anom <- current$v_mean - mean(her_cur$v)

p0 <- ggplot(current, aes(x = year, y = v_anom)) +
  geom_line(col = "darkgreen") +
  geom_point(size = 2, col = "darkgreen") +
  # geom_hline(yintercept = 0, linetype="dashed") +
  xlab("Year") +
  ylab("Mean northwards current speed (m/s)") +
  theme_bw()
p0

p1 <- ggplot(sims, aes(x = year)) +
  geom_line(aes(y = sdist_accum, col = sim)) +
  geom_point(aes(y = sdist_accum, col = sim), sim_obs_end, size = 2) +
  scale_color_manual(
    name = "", 
    labels = c("27 cm", "30 cm", "33 cm", "36 cm", "2016 year-class"),
    values = c(2,3,4,5,1)) +
  xlab("Year") +
  ylab("Actual swimming distance (km)") +
  theme_bw()
p1

p2 <- ggplot(sims, aes(x = year)) +
  geom_line(aes(y = energy_rel, col = sim)) +
  geom_point(aes(y = energy_rel, col = sim), sim_obs_end, size = 2) +
  scale_color_manual(
    name = "", 
    labels = c("27 cm", "30 cm", "33 cm", "36 cm", "2016 year-class"),
    values = c(2,3,4,5,1)) +
  xlab("Year") +
  ylab("Energy loss (%)") +
  theme_bw()
p2

p3 <- ggplot(sims, aes(x = year)) +
  geom_line(aes(y = cond_k, col = sim)) +
  geom_point(aes(y = cond_k, col = sim), sim_obs_end, size = 2) +
  scale_color_manual(
    name = "", 
    labels = c("27 cm", "30 cm", "33 cm", "36 cm", "2016 year-class"),
    values = c(2,3,4,5,1)) +
  xlab("Year") +
  ylab("CF after spawning") +
  theme_bw()
p3

# Combine plots
p0_2020 + p0_2021 + p0_2022 + p0_2023 + p1 + p2 + p3 + guide_area() + plot_layout(guides = "collect", ncol = 4) & plot_annotation(tag_levels = "a")
ggsave("migration.tiff", scale = 1.1, dpi = 300, width = 9.5, height = 6, compression = "lzw")
