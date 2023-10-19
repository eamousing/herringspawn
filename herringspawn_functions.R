get_k_trans <- function(sim_df, thold) {
  mod <- lm(sim_df$trans ~ sim_df$cond_k)
  return(mod$coefficients[2] * thold + mod$coefficients[1])
}

get_k_dist <- function(sim_df, thold) {
  mod <- lm(sim_df$sdist_accum ~ sim_df$cond_k)
  return(mod$coefficients[2] * thold + mod$coefficients[1])
}

calc_tot_energy <- function(dwei, params) {
  tot_energy <- (dwei * bio_params$fat_frac * bio_params$ed_fat) +
    (dwei * bio_params$gonad_frac * bio_params$ed_gonad) +
    (dwei * bio_params$solid_frac * bio_params$ed_solid)
  return(tot_energy)
}

get_route_positions <- function(data_dir) {
  require(xlsx)
  
  # Read spawning location data
  her_df <- read.xlsx2(paste0(data_dir, "LUF11_2018_AllVessels_Herring_NASC.xlsx"), 1)
  names(her_df) <- c("date", "start_log", "stop_log", "lat", "lon", "NASC")
  her_df$date <- as.Date(her_df$date, format = "%Y%m%d")
  her_df$lon <- as.numeric(her_df$lon)
  her_df$lat <- as.numeric(her_df$lat)
  her_df$NASC <- as.numeric(her_df$NASC)
  
  # Calculate center of gravity
  lat_slices <- seq(min(her_df$lat), max(her_df$lat), 0.5)
  her_cog <- data.frame(
    lat = vector(length = length(lat_slices)-1),
    lon = vector(length = length(lat_slices)-1))
  for (i in 1:(length(lat_slices)-1)) {
    her_subset <- her_df[her_df$lat >= lat_slices[i] &
                           her_df$lat < lat_slices[i+1],]
    for (j in 1:length(her_subset$lat)) {
      her_cog$lat[i] <- her_cog$lat[i] +
        her_subset$lat[j] * her_subset$NASC[j]
      her_cog$lon[i] <- her_cog$lon[i] +
        her_subset$lon[j] * her_subset$NASC[j]
    }
    her_cog$lat[i] <- her_cog$lat[i] / sum(her_subset$NASC)
    her_cog$lon[i] <- her_cog$lon[i] / sum(her_subset$NASC)
  }
  her_cog$lat <- rev(her_cog$lat)
  her_cog$lon <- rev(her_cog$lon)
  her_cog$x <- 0
  her_cog$y <- 0
  her_cog <- as_tibble(her_cog)
  return(her_cog)
}

get_currents <- function(data_dir, positions) {
  source("bilin_inv.R")
  require(R.matlab)
  
  # Read data from file
  jan_feb <- readMat(paste0(data_dir, "currents1993_2023jan_feb.mat"))
  
  # Find domain coordinates for spawning positions
  for (i in 1:dim(positions)[1]) {
    xy <- bilin_inv(positions$lon[i], positions$lat[i], jan_feb$xlon, jan_feb$ylat, 105, 53)
    positions$x[i] <- xy[1]
    positions$y[i] <- xy[2]
  }
  
  # Extract currents speed along swimming vectors
  tr_out <- data.frame(
    year = NULL, x = NULL, y = NULL, trans = NULL,
    vnorth = NULL, ueast = NULL
  )
  
  for (j in 1:dim(jan_feb$year)[2]) {
    for (i in 1:(dim(positions)[1]-1)) {
      x <- positions$x[i:(i+1)]
      y <- positions$y[i:(i+1)]
      f <- approxfun(x, y = y)
      df <- data.frame(
        x = seq(as.integer(positions$x[i]), as.integer(positions$x[i+1]))
      )
      df$y <- as.integer(f(seq(as.integer(positions$x[i]),
                               as.integer(positions$x[i+1]))))
      if (i == 1) {
        tr_df <- data.frame(year = 1992+j, x = df$x, y = df$y, trans = i)
      } else {
        tr_df <- rbind(tr_df, 
                       data.frame(year = 1992+j, x = df$x, y = df$y, trans = i))
      }
    }
    tr_df <- na.omit(tr_df)
    for (i in 1:dim(tr_df)[1]) {
      tr_df$vnorth[i] <- jan_feb$vnorth[tr_df$x[i], tr_df$y[i], j]
      tr_df$ueast[i] <- jan_feb$ueast[tr_df$x[i], tr_df$y[i], j]
    }
    tr_out <- rbind(tr_out, tr_df)
  }
  
  # Calculate mean currents speed for each transect and year
  cur_temp <- tr_out %>% select(year, trans, vnorth, ueast) %>%
    group_by(year, trans) %>% 
    summarise(v = mean(vnorth), u = mean(ueast))
  
  return(cur_temp)
}

energy_consumption_kJ <- function(len, wei, temp, sspeed, param) {
  rft <- exp(param$rq * temp) # Temperature impact
  res <- param$ra * wei^param$rb # Size impact
  vel <- sspeed * len # Swimming speed (cm s-1)
  activ <- exp(param$rt0 * vel) # Activity impact
  res <- res * rft * activ # specific Respiration (o2)
  res <- res * param$oxy_coef * 0.75 # specific respiration (J g-1)
  res <- res * wei # Respiration (J)
  res <- res / 1000 # J -> kJ
  res <- res / 86400 # kJ days-1 -> kJ seconds-1
  return(res)
}

#' Calculates distance and energy consumption for swimming a transect
#' while under the influence of advection
#'
#' @param lon0 Start position longitude
#' @param lat0 Start position latitude
#' @param lon1 End position longitude
#' @param lat1 End position latitude
#' @param len Fish length (cm)
#' @param wei Fish weight (grams)
#' @param te Water temperature (degC)
#' @param u Eastern current speed (m sec-1)
#' @param v Northern current speed (m sec-1)
#'
#' @return Vector with distance covered (km), time used (days) and energy used (kJ)
swim_transect_fixed <- function(lon0, lat0, lon1, lat1, len, u, v, param) {
  require(geosphere)
  r_earth <- 6375 # Earth radius (km)
  
  # Calculate unmodified swimming distance and time
  dist0 <- distm(c(lon0, lat0), c(lon1, lat1))
  t0 <- dist0 / (len * param$sspeed * 0.01)
  
  vel_east <- u * t0 / 1000
  vel_north <- v * t0 / 1000
  
  # Iterative estimation of distance and time taking into account currents
  t_diff <- 1e6 # Time accuracy
  while (t_diff > 60) { # Repeat until accuracy is below 60 seconds
    # Add advection to starting point
    lon_a <- lon0 + (vel_east / r_earth) * (180 / pi) * cos(lat0 * pi / 180)
    lat_a <- lat0 + (vel_north / r_earth) * (180 / pi)
    
    # Get new distance and time
    dist_a <- distm(c(lon_a, lat_a), c(lon1, lat1))
    t_a <- dist_a / (len * param$sspeed * 0.01)
    
    # Calculate time difference and adjust starting point
    t_diff <- t_a - t0
    vel_east <- vel_east + u * t_diff / 1000
    vel_north <- vel_north + v * t_diff / 1000
    t0 <- t_a
  }
  
  # Calculate energy consumption
  dist_out <- dist_a / 1000
  t_out <- t_a
  # t_out <- t_a / 86400
  # energy_out <- energy_consumption_kJ(len, wei, te, param) * t_a
  # return(c(dist_out, t_out, energy_out))
  return(c(dist_out, t_out))
}

swim_transect_variable <- function(lon0, lat0, lon1, lat1, len, u, v, param) {
  require(geosphere)
  r_earth <- 6375.0 # Earth radius (km)
  
  # Calculate unmodified swimming distance and time
  dist0 <- distm(c(lon0, lat0), c(lon1, lat1))
  t0 <- dist0/(len*param$sspeed*0.01)
  
  # Calculate advection displacement vectors
  vel_east <- u*t0/1000.0
  vel_north <- v*t0/1000.0
  
  # Find new starting point assuming advection
  lon_a <- lon0 + (vel_east/r_earth)*(180.0/pi)*cos(lat0*pi/180.0)
  lat_a <- lat0 + (vel_north/r_earth)*(180.0/pi)
  
  # Get swimming distance and time with advection
  dist_a <- distm(c(lon_a, lat_a), c(lon1, lat1))
  
  # Get new swimming speed
  sspeed <- dist_a/(t0*len*0.01)
  
  # Check that new swimming time == t0
  t_a <- dist_a/(len*sspeed)*100.0
  if (diff(c(t0, t_a)) > 1e-3) {
    print("WARNING: t0 and t_a are not equal")
  }

  return(c(dist_a/1000.0, t0, sspeed))
}

simulate_migration <- function(fish_len, fish_wei, yrs, trans, params) {
  # Create return data frame
  sim_out <- data.frame(
    year = NULL, # Year
    trans = NULL, # Transect
    len = NULL, # Length
    wei = NULL, # Wet weight
    fat = NULL, # Fat weight (g DW)
    gonad = NULL, # Gonad weight (g DW)
    solid = NULL, # Somatic weight (g DW)
    sspeed = NULL, # Swimming speed (m s-1)
    energy_tot = NULL, # Total energy content (kJ)
    cond = NULL, # Condition factor
    cond_k = NULL, # Condition factor after spawning
    sdist = NULL, # Swimming distance (km)
    stime = NULL, # Swimming time (days)
    energy = NULL, # Energy consumption (kJ)
    sdist_accum = NULL, # Accumulated swimming distance (km)
    stime_accum = NULL, # Accumulated swimming time (days)
    energy_accum = NULL # Accumulated energy consumption (kJ)
  )
  
  # Year loop
  inc <- 0
  for (yr in 1:length(yrs)) {
    
    # Create fish traits
    sind <- list()
    sind$len <- fish_len
    sind$wei <- fish_wei
    sind$dwei <- sind$wei * (1 - params$wc)
    sind$fat_mass <- sind$dwei * params$fat_frac
    sind$gonad_mass <- sind$dwei * params$gonad_frac
    sind$solid_mass <- sind$dwei * params$solid_frac
    # Get currents for transect
    cur_yr <- her_cur[her_cur$year == yrs[yr], ]
    
    # Transect loop
    day <- 0
    for (i in trans) {
      inc <- inc + 1
      x0 <- her_cog$lon[i]
      y0 <- her_cog$lat[i]
      x1 <- her_cog$lon[i+1]
      y1 <- her_cog$lat[i+1]
      u0 <- cur_yr$u[i]
      v0 <- cur_yr$v[i]
      
      # Get swimming distance and time used
      swim <- swim_transect_variable(x0, y0, x1, y1, sind$len, u0, v0, params)
      
      # Calculate energy consumption
      energy_use <- energy_consumption_kJ(sind$len, sind$wei, temp, swim[3], params)*swim[2]
      
      # Get fat burned
      fat_use_dwei <- energy_use / params$ed_fat
      solid_use_dwei <- 0
      
      # If fat is is depleted, use somatic tissue to fuel migration
      if (sind$fat_mass - fat_use_dwei < 0) {
        fat_use_dwei <- 0
        solid_use_dwei <- energy_use / params$ed_solid
      }
      
      # Calculate somatic tissue converted to gonads
      ig_old <- sind$gonad_mass / sind$dwei
      day <- day + swim[2] / 86400
      ig_target <- (0.2561 * day + 9.7436) / 100
      ig_inc <- 0
      gonad_inc <- 0
      # print(c(day, sind$dwei, ig_old, ig_target))
      if (ig_old < ig_target) {
        ig_inc <- ig_target - ig_old
        gonad_inc <- ig_inc * sind$dwei
      }
      
      # Update state
      sind$fat_mass <- sind$fat_mass - fat_use_dwei
      sind$gonad_mass <- sind$gonad_mass + gonad_inc
      sind$solid_mass <- sind$solid_mass - (gonad_inc / params$ed_solid) * params$ed_gonad - solid_use_dwei
      sind$dwei <- sind$fat_mass + sind$gonad_mass + sind$solid_mass
      sind$wei <- sind$dwei / (1 - params$wc)

      # Concatenate results
      sim_out <- 
        rbind(sim_out,
              data.frame(
                year = yrs[yr],
                trans = i,
                len = sind$len,
                wei = sind$wei,
                fat = sind$fat_mass,
                gonad = sind$gonad_mass,
                solid = sind$solid_mass,
                sspeed = swim[3],
                energy_tot = sind$fat_mass * params$ed_fat + sind$gonad_mass * params$ed_gonad + sind$solid_mass * params$ed_solid,
                cond = 100 * (sind$wei / sind$len^3),
                cond_k = 100 * (((sind$fat_mass + sind$solid_mass) / (1 - params$wc)) / sind$len^3),
                sdist = swim[1],
                stime = swim[2] / 86400,
                energy = energy_use
              ))
    }
  }
  
  # Calculate accumulated distance, time and energy use
  sim_out$sdist_accum <- 0
  sim_out$stime_accum <- 0
  sim_out$energy_accum <- 0
  inc <- 0
  for (j in 1:length(yrs)) {
    nb_trans <- max(sim_out$trans[sim_out$year == yrs[j]])
    for (i in 1:nb_trans) {
      inc <- inc + 1
      if (i == 1) {
        sim_out$sdist_accum[inc] <- sim_out$sdist[inc]
        sim_out$stime_accum[inc] <- sim_out$stime[inc]
        sim_out$energy_accum[inc] <- sim_out$energy[inc]
      } else {
        sim_out$sdist_accum[inc] <- sim_out$sdist_accum[inc-1] + sim_out$sdist[inc]
        sim_out$stime_accum[inc] <- sim_out$stime_accum[inc-1] + sim_out$stime[inc]
        sim_out$energy_accum[inc] <- sim_out$energy_accum[inc-1] + sim_out$energy[inc]
      }
    }
  }
  
  # Return results
  return(sim_out)
}