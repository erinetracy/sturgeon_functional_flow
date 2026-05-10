# =============================================================================
# GREEN STURGEON TELEMETRY ANALYSIS — UPDATED 2026
# Author: E. Tracy
# Last updated: April 2026
#
# Description: Full pipeline for green sturgeon acoustic telemetry data.
# Builds on cleaned detection events from detection_cleaning_gsflow.R and
# migration status classification. Calculates migration routes, flow and
# temperature covariates, runs GLMs and Boosted Regression Trees.
#
# Study periods:
#   Route selection:    2007-2017 (requires Delta/Sac discriminating receivers)
#   Migration success:  2007-2022 (requires Benicia/Carquinez coverage only)
# =============================================================================

# =============================================================================
# REQUIRED DATASETS
# =============================================================================
#
# INPUT FILES:
#   cleaned_data/events_with_receivergroups_032026.csv  — cleaned detection
#     events with receiver group labels, water year, and migration status
#   flow/HMC_variables_final_FULL.csv                  — annual HMC predictors
#   flow/VON_annual_flow_result_updated_FULL.csv        — annual VON predictors
#   cleaned_data/daily_gauge_data_final.csv             — daily flow (cfs)
#   cleaned_data/Rio_Vista_daily_temp.csv               — daily temperature
#
# GENERATED OUTPUTS:
#   cleaned_data/up_migration_routes.csv               — route strings per fish
#   cleaned_data/glm_data.csv                          — final analysis dataset
#   cleaned_data/brt_seed_results_full.csv             — BRT stability results
#   cleaned_data/brt_summary_table.csv                 — BRT summary table
#   gs_analysis_workspace.RData                        — full workspace save

# =============================================================================
# 0. LIBRARIES
# =============================================================================

library(glatos)        # acoustic telemetry processing
library(data.table)    # fast data manipulation
library(lubridate)     # date/time handling
library(tidyverse)     # core data wrangling + ggplot2
library(dplyr)         # explicit namespace for joins/filters
library(ggplot2)       # plotting
library(stringdist)    # Levenshtein distance (kept for reference)
library(corrplot)      # correlation matrix plots
library(reshape2)      # melt() for correlation heatmap
library(gbm)           # gradient boosted models
library(dismo)         # gbm.step wrapper for BRT
library(patchwork)     # combine ggplot panels
library(ggspatial)     # scale bar and north arrow
library(sf)            # spatial data
library(tigris)        # Census TIGER shapefiles

# =============================================================================
# 1. SET WORKING DIRECTORY AND LOAD DATA
# =============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Subfolder paths
cleaned_data <- "cleaned_data"
flow         <- "flow"
figures      <- "figures"
raw_data     <- "raw_data"

# Load cleaned detection events from 01_detection_cleaning.R
# This file already has: group (1-24), receiver_group, water_year,
# migration status (up_complete, up_incomplete, down_complete, etc.)
events_2026 <- read.csv(
  file.path(cleaned_data, "events_with_receivergroups_032026.csv")
) %>%
  mutate(
    first_detection = as.POSIXct(first_detection, tz = "UTC"),
    last_detection  = as.POSIXct(last_detection,  tz = "UTC")
  )
#had to update new receivers and used previous receiver group numbers to do this 
#so its approximately right, maybe a few errors but works for analysis of sac vs. delta

# =============================================================================
# 2. FILTER TO UPSTREAM MIGRANTS ONLY
# =============================================================================
# Use migration status from cleaning script:
#   up_complete   = fish detected at ocean side → benicia/carquinez →
#                   spawning grounds (confirmed full upstream migration)
#   up_incomplete = fish detected at ocean side → benicia/carquinez but
#                   never reached spawning grounds (failed/incomplete migration)
#
# NOTE: Entry point is benicia/carquinez — the Carquinez Strait — which
# confirms fish have committed to the freshwater migration corridor.
# Golden Gate/Bay detections must precede Benicia to exclude downstream
# migrants. This is enforced in the migration status classification in
# 01_detection_cleaning.R via the ocean_before_bc flag.
#
# NOTE: 5 fish-year combinations within 2007-2017 missed Benicia/Carquinez
# but reached upstream receivers (max group 22-24), attributed to receiver
# gaps rather than true bypass. These fish received NA status in the
# classification and are excluded from analysis. Their exclusion has
# negligible effect on results (<1% of sample).

upstream_migrants <- events_2026 %>%
  filter(status %in% c("up_complete", "up_incomplete"))

cat("Total upstream migrant detections:", nrow(upstream_migrants), "\n")
cat("Unique fish:", n_distinct(upstream_migrants$animal_id), "\n")
cat("Water years:", paste(range(upstream_migrants$water_year), collapse = "-"), "\n")

# =============================================================================
# 3. BUILD MIGRATION SUMMARY (one row per fish-year)
# =============================================================================
# Captures key timestamps for each fish-year combination:
#   first_detection_startup_1  = first detection at group 1 (SFE entry)
#   last_detection_startup_1   = last detection at group 1
#   last_detection_startup_2   = last detection at group 2 (Benicia Bridge)
#   first_detection_endup      = first detection at group 24 (spawning grounds)
#   migration_binary           = 1 (up_complete) or 0 (up_incomplete)

migration_summary <- upstream_migrants %>%
  group_by(animal_id, water_year, status) %>%
  summarise(
    first_detection_startup_1 = suppressWarnings(
      min(first_detection[group == 1], na.rm = TRUE)),
    last_detection_startup_1  = suppressWarnings(
      max(last_detection[group == 1],  na.rm = TRUE)),
    last_detection_startup_2  = suppressWarnings(
      max(last_detection[group == 2],  na.rm = TRUE)),
    first_detection_endup     = suppressWarnings(
      min(first_detection[group == 24], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    # Convert Inf/-Inf to NA for datetime columns only
    # (occurs when fish never detected at that group)
    first_detection_startup_1 = if_else(is.infinite(first_detection_startup_1),
                                        as.POSIXct(NA), first_detection_startup_1),
    last_detection_startup_1  = if_else(is.infinite(last_detection_startup_1),
                                        as.POSIXct(NA), last_detection_startup_1),
    last_detection_startup_2  = if_else(is.infinite(last_detection_startup_2),
                                        as.POSIXct(NA), last_detection_startup_2),
    first_detection_endup     = if_else(is.infinite(first_detection_endup),
                                        as.POSIXct(NA), first_detection_endup),
    migration_binary = if_else(status == "up_complete", 1L, 0L),
    migration_type   = "up"
  )

cat("Migration summary rows:", nrow(migration_summary), "\n")
cat("up_complete:",   sum(migration_summary$migration_binary == 1), "\n")
cat("up_incomplete:", sum(migration_summary$migration_binary == 0), "\n")


# Fix status for Benicia-gap fish in events_2026
events_2026 <- events_2026 %>%
  mutate(status = case_when(
    # GS0303 2010: valid upstream migration entering via group 4 after
    # failed first attempt Nov 2009, reached spawning grounds March 2010
    animal_id == "UCDHIST-GS0303-2006-06-28" & water_year == 2010 ~ "up_complete",
    # GS0704, GS0723, GS0746 2012: clean upstream migration sequences
    # ocean → golden gate → bay → Sacramento mainstem → spawning grounds
    # missed Benicia due to detection gap — confirmed legitimate migrants
    animal_id == "UCDHIST-GS0704-2011-08-11" & water_year == 2012 ~ "up_complete",
    animal_id == "UCDHIST-GS0723-2011-08-19" & water_year == 2012 ~ "up_complete",
    animal_id == "UCDHIST-GS0746-2011-09-02" & water_year == 2012 ~ "up_complete",
    TRUE ~ status
  ))

# Verify
events_2026 %>%
  filter(animal_id %in% c("UCDHIST-GS0303-2006-06-28",
                          "UCDHIST-GS0704-2011-08-11",
                          "UCDHIST-GS0723-2011-08-19",
                          "UCDHIST-GS0746-2011-09-02")) %>%
  distinct(animal_id, water_year, status) %>%
  arrange(animal_id, water_year)

# Save updated events
write.csv(events_2026,
          file.path(cleaned_data, "events_with_receivergroups_032026.csv"),
          row.names = FALSE)


# =============================================================================
# 4. ROUTE CLASSIFICATION — 2007-2017 ONLY
# =============================================================================

# --- 4a. Filter to up-migration window for route analysis (2007-2017) ---
# Window: group 1 entry → group 24 arrival
# Using detections_up prevents downstream passage at group 2 from
# contaminating migration start time — mirrors original analysis approach
detections_up <- upstream_migrants %>%
  filter(water_year >= 2007 & water_year <= 2017) %>%
  group_by(animal_id, water_year) %>%
  filter(any(group == 1)) %>%
  filter(
    first_detection >= suppressWarnings(
      min(first_detection[group == 1], na.rm = TRUE)) &
      (is.na(suppressWarnings(
        min(first_detection[group == 24], na.rm = TRUE))) |
         first_detection <= suppressWarnings(
           min(first_detection[group == 24], na.rm = TRUE)))
  ) %>%
  ungroup() %>%
  arrange(animal_id, water_year, first_detection)

cat("Detections in up-migration window (2007-2017):", nrow(detections_up), "\n")
cat("Unique fish-years:", n_distinct(paste(detections_up$animal_id,
                                           detections_up$water_year)), "\n")

# --- 4b. Filter to up-migration window for success analysis (2007-2022) ---
# Window: group 2 entry (Benicia) → group 24 arrival
# All fish in new dataset must have hit Benicia so start from group 2
detections_up_success <- upstream_migrants %>%
  filter(water_year >= 2007 & water_year <= 2022) %>%
  group_by(animal_id, water_year) %>%
  filter(any(group == 2)) %>%
  filter(
    first_detection >= suppressWarnings(
      min(first_detection[group == 2], na.rm = TRUE)) &
      (is.na(suppressWarnings(
        min(first_detection[group == 24], na.rm = TRUE))) |
         first_detection <= suppressWarnings(
           min(first_detection[group == 24], na.rm = TRUE)))
  ) %>%
  ungroup() %>%
  arrange(animal_id, water_year, first_detection)

cat("Detections in up-migration window (2007-2022):", nrow(detections_up_success), "\n")

# --- 4c. Build unique-receiver route strings ---
fish_routes <- detections_up %>%
  arrange(animal_id, water_year, first_detection) %>%
  group_by(animal_id, water_year) %>%
  summarize(
    route = paste(unique(group[!is.na(group)]), collapse = " -> "),
    .groups = "drop"
  )

cat("Total routes built:", nrow(fish_routes), "\n")
cat("Routes containing NA:", sum(str_detect(fish_routes$route, "NA")), "\n")

# --- 4d. Remove sparse routes ---
fish_routes_filtered <- fish_routes %>%
  filter(!route %in% c(
    "1 -> 4 -> 5 -> 24",
    "1 -> 5 -> 6 -> 24",
    "1 -> 24",
    "1"
  ))

cat("Routes after sparse filter:", nrow(fish_routes_filtered), "\n")

# --- 4e. Route classification ---
# Preserve validated labels from original glm_data for existing fish
# Manually classify only new fish added in updated dataset

# New fish not in original glm_data
new_fish_routes <- fish_routes_filtered %>%
  anti_join(glm_data, by = c("animal_id", "water_year"))

cat("New fish-years to classify:", nrow(new_fish_routes), "\n")

# Classify new fish — only 1 reached group 24 (Delta route)
new_fish_labelled <- new_fish_routes %>%
  mutate(
    route_cluster = case_when(
      route == "1 -> 2 -> 4 -> 5 -> 6 -> 13 -> 17 -> 20 -> 24" ~ 1,
      TRUE ~ NA_real_
    ),
    migration_binary = case_when(
      route == "1 -> 2 -> 4 -> 5 -> 6 -> 13 -> 17 -> 20 -> 24" ~ 1L,
      TRUE ~ 0L
    )
  ) %>%
  dplyr::select(animal_id, water_year, route, route_cluster, migration_binary)

# Get validated labels from original glm_data for existing fish
original_fish_labelled <- fish_routes_filtered %>%
  inner_join(
    glm_data %>%
      dplyr::select(animal_id, water_year, route_cluster, migration_binary) %>%
      mutate(route_cluster    = as.double(route_cluster),
             migration_binary = as.integer(migration_binary)),
    by = c("animal_id", "water_year")
  )

# Combine
fish_routes_labelled <- bind_rows(original_fish_labelled, new_fish_labelled)

cat("\nRoute classification summary:\n")
fish_routes_labelled %>% count(route_cluster, migration_binary)

write.csv(fish_routes_labelled,
          file.path(cleaned_data, "up_migration_routes.csv"),
          row.names = FALSE)

# =============================================================================
# 5. BUILD MIGRATION SUMMARIES
# =============================================================================

# --- 5a. Route analysis migration summary (2007-2017) ---
# Built from detections_up — already filtered to up-migration window
# so max(last_detection[group==2]) correctly captures upstream passage only
route_data <- detections_up %>%
  group_by(animal_id, water_year) %>%
  summarise(
    first_detection_startup_1 = suppressWarnings(
      min(first_detection[group == 1], na.rm = TRUE)),
    last_detection_startup_1  = suppressWarnings(
      max(last_detection[group == 1], na.rm = TRUE)),
    last_detection_startup_2  = suppressWarnings(
      max(last_detection[group == 2], na.rm = TRUE)),
    first_detection_endup     = suppressWarnings(
      min(first_detection[group == 24], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    first_detection_startup_1 = if_else(is.infinite(first_detection_startup_1),
                                        as.POSIXct(NA), first_detection_startup_1),
    last_detection_startup_1  = if_else(is.infinite(last_detection_startup_1),
                                        as.POSIXct(NA), last_detection_startup_1),
    last_detection_startup_2  = if_else(is.infinite(last_detection_startup_2),
                                        as.POSIXct(NA), last_detection_startup_2),
    first_detection_endup     = if_else(is.infinite(first_detection_endup),
                                        as.POSIXct(NA), first_detection_endup)
  ) %>%
  left_join(
    fish_routes_labelled %>%
      dplyr::select(animal_id, water_year, route_cluster, migration_binary),
    by = c("animal_id", "water_year")
  ) %>%
  mutate(migration_type = "up")

cat("Route data rows:", nrow(route_data), "\n")
cat("Delta:", sum(route_data$route_cluster == 1, na.rm = TRUE), "\n")
cat("Mainstem Sac:", sum(route_data$route_cluster == 0, na.rm = TRUE), "\n")
cat("Incomplete:", sum(is.na(route_data$route_cluster)), "\n")

# --- 5b. Migration success summary (2007-2022) ---
# Built from detections_up_success — filtered to up-migration window
# Waiting for 2018-2022 flow metrics before running BRT
success_data <- detections_up_success %>%
  group_by(animal_id, water_year) %>%
  summarise(
    last_detection_startup_2 = suppressWarnings(
      max(last_detection[group == 2], na.rm = TRUE)),
    first_detection_endup    = suppressWarnings(
      min(first_detection[group == 24], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    last_detection_startup_2 = if_else(is.infinite(last_detection_startup_2),
                                       as.POSIXct(NA), last_detection_startup_2),
    first_detection_endup    = if_else(is.infinite(first_detection_endup),
                                       as.POSIXct(NA), first_detection_endup)
  ) %>%
  left_join(
    upstream_migrants %>%
      distinct(animal_id, water_year, status) %>%
      mutate(migration_binary = if_else(status == "up_complete", 1L, 0L)),
    by = c("animal_id", "water_year")
  ) %>%
  filter(!is.na(status), status != "incomplete_dead") %>%
  mutate(migration_type = "up")

cat("\nSuccess data rows:", nrow(success_data), "\n")
cat("Complete:", sum(success_data$migration_binary == 1), "\n")
cat("Incomplete:", sum(success_data$migration_binary == 0), "\n")

write.csv(success_data,
          file.path(cleaned_data, "success_data_2007_2022.csv"),
          row.names = FALSE)



# =============================================================================
# 6. REPEAT SPAWNER ANALYSIS
# =============================================================================

repeat_spawners <- route_data %>%
  group_by(animal_id) %>%
  summarize(
    unique_routes     = n_distinct(route_cluster, na.rm = TRUE),
    routes_selected   = paste(sort(unique(route_cluster)), collapse = ", "),
    routes_per_year   = paste(unique(paste(water_year, route_cluster,
                                           sep = "-")), collapse = ", ")
  )

cat("\nFish detected in multiple years:", nrow(repeat_spawners), "\n")
cat("Fish that switched routes:\n")
repeat_spawners %>% filter(unique_routes > 1) %>% print()

# =============================================================================
# 7. MIGRATION TIMING
# =============================================================================

# --- 7a. Estimate mean travel time from group 1 to group 2 ---
avg_travel_time_hours <- route_data %>%
  filter(
    !is.na(last_detection_startup_1),
    !is.na(last_detection_startup_2),
    !is.na(route_cluster)
  ) %>%
  mutate(time_diff_hr = as.numeric(
    difftime(last_detection_startup_2, last_detection_startup_1,
             units = "hours")
  )) %>%
  filter(time_diff_hr >= 0, time_diff_hr <= 48) %>%
  summarise(mean_time_diff_hr = mean(time_diff_hr, na.rm = TRUE)) %>%
  pull(mean_time_diff_hr)

cat("Mean group 1 to 2 travel time (hours):", round(avg_travel_time_hours, 2), "\n")

# --- 7b. Define migration start time per fish ---
# Route data: use group 2 if available, impute from group 1 + travel time
# if fish missed Benicia. These fish were validated in original analysis.
# Success data: all fish hit Benicia so use group 2 directly.
route_data <- route_data %>%
  mutate(
    last_detection_startup = case_when(
      !is.na(last_detection_startup_2) ~ last_detection_startup_2,
      !is.na(last_detection_startup_1) ~ last_detection_startup_1 +
        lubridate::dseconds(avg_travel_time_hours * 3600),
      TRUE ~ as.POSIXct(NA)
    )
  )

success_data <- success_data %>%
  mutate(last_detection_startup = last_detection_startup_2)

cat("Route data with startup:", sum(!is.na(route_data$last_detection_startup)), "\n")
cat("Success data with startup:", sum(!is.na(success_data$last_detection_startup)), "\n")


# =============================================================================
# 8. FLOW COVARIATE CALCULATION — SUCCESSFUL MIGRATIONS
# =============================================================================
# Mean daily discharge at each gauge during each fish's migration window
# (group 2 exit → group 24 arrival).

gauge_data <- read.csv(file.path(cleaned_data, "daily_gauge_data_final.csv"))

# Helper function to calculate mean discharge per fish-year-gauge
calculate_discharge <- function(data, gauge_data) {
  
  # Expand to fish × year × gauge grid
  fish_gauge <- expand.grid(
    animal_id  = unique(data$animal_id),
    water_year = unique(data$water_year),
    gauge      = unique(gauge_data$gauge),
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      data %>%
        dplyr::select(animal_id, water_year,
                      last_detection_startup, first_detection_endup),
      by = c("animal_id", "water_year")
    ) %>%
    filter(!is.na(last_detection_startup) & !is.na(first_detection_endup))
  
  # Join flow data and flag dates within migration window
  fish_gauge %>%
    left_join(gauge_data, by = "gauge") %>%
    mutate(
      in_migration          = date >= last_detection_startup &
        date <= first_detection_endup,
      flow_cfs_in_migration = ifelse(in_migration, flow_cfs, NA_real_)
    ) %>%
    group_by(animal_id, water_year, gauge) %>%
    summarise(
      mean_discharge = mean(flow_cfs_in_migration, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(mean_discharge = ifelse(is.nan(mean_discharge),
                                   NA_real_, mean_discharge)) %>%
    pivot_wider(
      names_from  = gauge,
      values_from = mean_discharge,
      names_glue  = "{gauge}_mean_discharge"
    )
}

# Calculate discharge for route and success datasets
route_discharge   <- calculate_discharge(route_data,   gauge_data)
success_discharge <- calculate_discharge(success_data, gauge_data)

# =============================================================================
# 9. FLOW COVARIATE CALCULATION — FAILED MIGRATIONS
# =============================================================================
# For incomplete migrations, window = group 2 start → last receiver visited.

calculate_failed_discharge <- function(data, events, gauge_data,
                                       avg_travel_time_hours) {
  # Filter failed fish
  failed_ids <- data %>%
    filter(migration_binary == 0)
  
  # Build startup time — handle datasets with or without group 1 column
  if ("last_detection_startup_1" %in% names(failed_ids)) {
    failed_ids <- failed_ids %>%
      dplyr::select(animal_id, water_year,
                    last_detection_startup_1, last_detection_startup_2) %>%
      mutate(
        last_detection_startup_final = coalesce(
          last_detection_startup_2,
          last_detection_startup_1 + lubridate::dhours(avg_travel_time_hours)
        )
      )
  } else {
    failed_ids <- failed_ids %>%
      dplyr::select(animal_id, water_year, last_detection_startup_2) %>%
      mutate(last_detection_startup_final = last_detection_startup_2)
  }
  
  # Get last detection at highest group visited
  failed_last <- events %>%
    mutate(last_detection = as.POSIXct(last_detection)) %>%
    semi_join(failed_ids, by = c("animal_id", "water_year")) %>%
    group_by(animal_id, water_year) %>%
    filter(group == max(group, na.rm = TRUE)) %>%
    summarise(last_detection_final = max(last_detection, na.rm = TRUE),
              .groups = "drop")
  
  # Merge start/end timestamps
  failed_info <- failed_ids %>%
    left_join(failed_last, by = c("animal_id", "water_year")) %>%
    filter(!is.na(last_detection_startup_final) &
             !is.na(last_detection_final))
  
  # Calculate discharge
  failed_info %>%
    dplyr::select(animal_id, water_year,
                  last_detection_startup_final, last_detection_final) %>%
    crossing(gauge = unique(gauge_data$gauge)) %>%
    left_join(gauge_data, by = "gauge",
              relationship = "many-to-many") %>%
    mutate(
      in_migration = date >= as.Date(last_detection_startup_final) &
        date <= as.Date(last_detection_final),
      flow_cfs_in_migration = ifelse(in_migration, flow_cfs, NA_real_)
    ) %>%
    group_by(animal_id, water_year, gauge) %>%
    summarise(mean_discharge = mean(flow_cfs_in_migration, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(mean_discharge = ifelse(is.nan(mean_discharge),
                                   NA_real_, mean_discharge)) %>%
    pivot_wider(
      names_from  = gauge,
      values_from = mean_discharge,
      names_glue  = "{gauge}_mean_discharge"
    )
}

# Rerun with detections_up as events source
route_failed_discharge <- calculate_failed_discharge(
  route_data, detections_up, gauge_data, avg_travel_time_hours)

success_failed_discharge <- calculate_failed_discharge(
  success_data, detections_up_success, gauge_data, avg_travel_time_hours)

cat("Route failed discharge rows:", nrow(route_failed_discharge), "\n")
cat("Route failed has HMC:", sum(!is.na(route_failed_discharge$HMC_mean_discharge)), "\n")
cat("Success failed discharge rows:", nrow(success_failed_discharge), "\n")
cat("Success failed has HMC:", sum(!is.na(success_failed_discharge$HMC_mean_discharge)), "\n")


# =============================================================================
# 10. COMBINE SUCCESSFUL + FAILED DISCHARGE → FINAL GLM DATASETS
# =============================================================================

combine_discharge <- function(data, successful_discharge, failed_discharge) {
  data %>%
    left_join(successful_discharge, by = c("animal_id", "water_year")) %>%
    left_join(
      failed_discharge %>%
        rename(HMC_failed = HMC_mean_discharge,
               VON_failed = VON_mean_discharge),
      by = c("animal_id", "water_year")
    ) %>%
    mutate(
      HMC_mean_discharge = coalesce(HMC_mean_discharge, HMC_failed),
      VON_mean_discharge = coalesce(VON_mean_discharge, VON_failed)
    ) %>%
    dplyr::select(-HMC_failed, -VON_failed) %>%
    mutate(
      migration_duration = as.numeric(
        difftime(first_detection_endup, last_detection_startup,
                 units = "days")
      )
    )
}

glm_data_route   <- combine_discharge(route_data,
                                      route_discharge,
                                      route_failed_discharge)
glm_data_success <- combine_discharge(success_data,
                                      success_discharge,
                                      success_failed_discharge)

# Check
cat("GLM route rows:", nrow(glm_data_route), "\n")
cat("Has HMC — route:", sum(!is.na(glm_data_route$HMC_mean_discharge)), "\n")
cat("NA HMC — route:", sum(is.na(glm_data_route$HMC_mean_discharge)), "\n")

cat("\nGLM success rows:", nrow(glm_data_success), "\n")
cat("Has HMC — success:", sum(!is.na(glm_data_success$HMC_mean_discharge)), "\n")
cat("NA HMC — success:", sum(is.na(glm_data_success$HMC_mean_discharge)), "\n")

write.csv(glm_data_route,
          file.path(cleaned_data, "glm_data_route.csv"), row.names = FALSE)
write.csv(glm_data_success,
          file.path(cleaned_data, "glm_data_success.csv"), row.names = FALSE)

# =============================================================================
# 11. TEMPERATURE COVARIATE CALCULATION
# =============================================================================
# Mean water temperature ±1 week around each fish's group 2 departure time.

Rio_Vista_daily <- read.csv(
  file.path(flow, "Rio_Vista_daily_temp.csv")
) %>%
  mutate(date = as.Date(date))

calculate_temperature <- function(data, temp_data) {
  data %>%
    dplyr::select(animal_id, water_year, last_detection_startup) %>%
    mutate(
      last_detection_startup = as.Date(last_detection_startup),
      date_start = last_detection_startup - 7,
      date_end   = last_detection_startup + 7
    ) %>%
    rowwise() %>%
    mutate(dates = list(seq(date_start, date_end, by = "day"))) %>%
    unnest(cols = c(dates)) %>%
    ungroup() %>%
    left_join(temp_data, by = c("dates" = "date")) %>%
    mutate(
      period = case_when(
        dates < last_detection_startup ~ "before",
        dates > last_detection_startup ~ "after",
        TRUE                           ~ "exclude"
      )
    ) %>%
    filter(period != "exclude") %>%
    group_by(animal_id, water_year, period) %>%
    summarise(mean_temp = mean(daily_temp, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from  = period,
      values_from = mean_temp,
      names_glue  = "mean_temp_{period}_1week"
    )
}

temp_route   <- calculate_temperature(glm_data_route,   Rio_Vista_daily)
temp_success <- calculate_temperature(glm_data_success, Rio_Vista_daily)

# Join temperature to GLM datasets
glm_data_route   <- glm_data_route   %>%
  left_join(temp_route,   by = c("animal_id", "water_year"))
glm_data_success <- glm_data_success %>%
  left_join(temp_success, by = c("animal_id", "water_year"))

# =============================================================================
# 12. GLMs
# =============================================================================
# Route selection:     route_cluster ~ predictors  (binomial, 2007-2017)
# Migration success:   migration_binary ~ predictors (binomial, 2007-2022)
# Migration duration:  migration_duration ~ predictors (gaussian log, 2007-2017)

# --- Route selection ---
glm_route_flow <- glm(route_cluster ~ HMC_mean_discharge,
                      data = glm_data_route, family = binomial)
summary(glm_route_flow)

glm_route_verona <- glm(route_cluster ~ VON_mean_discharge,
                        data = glm_data_route, family = binomial)
summary(glm_route_verona)

glm_route_year <- glm(route_cluster ~ water_year,
                      data = glm_data_route, family = binomial)
summary(glm_route_year)

glm_route_temp <- glm(route_cluster ~ mean_temp_before_1week,
                      data = glm_data_route, family = binomial)
summary(glm_route_temp)

# --- Migration success ---
glm_success_flow <- glm(migration_binary ~ HMC_mean_discharge,
                        data = glm_data_success, family = binomial)
summary(glm_success_flow)

glm_success_temp <- glm(migration_binary ~ mean_temp_before_1week,
                        data = glm_data_success, family = binomial)
summary(glm_success_temp)

# --- Migration duration ---
glm_duration_route <- glm(migration_duration ~ route_cluster,
                          data = glm_data_route,
                          family = gaussian(link = "log"))
summary(glm_duration_route)

glm_duration_flow <- glm(migration_duration ~ VON_mean_discharge,
                         data = glm_data_route,
                         family = gaussian(link = "log"))
summary(glm_duration_flow)

glm_duration_temp <- glm(migration_duration ~ mean_temp_before_1week,
                         data = glm_data_route,
                         family = gaussian(link = "log"))
summary(glm_duration_temp)
hist(glm_data_route$migration_duration,
     main = "Distribution of migration duration (days)")

save.image(file = "gs_analysis_workspace_2026.RData")

# Save final GLM datasets
write.csv(glm_data_route,
          file.path(cleaned_data, "glm_data_route.csv"),
          row.names = FALSE)

write.csv(glm_data_success,
          file.path(cleaned_data, "glm_data_success.csv"),
          row.names = FALSE)

write.csv(success_data,
          file.path(cleaned_data, "success_data_2007_2022.csv"),
          row.names = FALSE)

cat("Saved glm_data_route:", nrow(glm_data_route), "rows\n")
cat("Saved glm_data_success:", nrow(glm_data_success), "rows\n")
cat("Saved success_data:", nrow(success_data), "rows\n")

# Join age estimates to success data
glm_data_success_age <- glm_data_success %>%
  left_join(
    fish_lengths_fixed %>% 
      dplyr::select(animal_id, age_est, length_cm),
    by = "animal_id"
  )

cat("Fish with age in success dataset:", 
    sum(!is.na(glm_data_success_age$age_est)), "\n")

# Quick GLM to test age effect on migration success
glm_age <- glm(migration_binary ~ age_est,
               data = glm_data_success_age,
               family = binomial)
summary(glm_age)
#older fish are statistically more likely to be successful but age calc was messy (had to estimate)
# GLM with fork length as predictor
# Rebuild with only plausible lengths
fish_lengths_plausible <- events_2026 %>%
  filter(animal_id %in% upstream_migrants$animal_id) %>%
  distinct(animal_id, length, lengthtype) %>%
  filter(!is.na(length)) %>%
  group_by(animal_id) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    length_cm = length * 100
  ) %>%
  filter(length_cm >= 60)

# Join to success data
glm_data_success_length <- glm_data_success %>%
  left_join(
    fish_lengths_plausible %>%
      dplyr::select(animal_id, length_cm, lengthtype),
    by = "animal_id"
  )

cat("Fish with plausible length in success dataset:",
    sum(!is.na(glm_data_success_length$length_cm)), "\n")

# GLM with plausible lengths only
glm_length_clean <- glm(migration_binary ~ length_cm,
                        data = glm_data_success_length,
                        family = binomial)
summary(glm_length_clean)

# =============================================================================
# 13. SUMMARY BAR CHART — MIGRATION ROUTES BY WATER YEAR (2007-2017)
# =============================================================================

route_graph <- glm_data_route %>%
  mutate(
    route_label = case_when(
      route_cluster == 1   ~ "Delta",
      route_cluster == 0   ~ "Mainstem Sacramento",
      is.na(route_cluster) ~ "Incomplete migration"
    ),
    route_label = factor(
      route_label,
      levels = c("Mainstem Sacramento", "Delta", "Incomplete migration")
    )
  )

route_counts <- route_graph %>%
  group_by(water_year, route_label) %>%
  summarise(count = n(), .groups = "drop")

route_props <- route_counts %>%
  group_by(water_year) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Shared theme
shared_theme <- theme_minimal(base_family = "Times New Roman") +
  theme(
    text             = element_text(family = "Times New Roman", color = "black"),
    axis.title       = element_text(color = "black", size = 12),
    axis.text        = element_text(color = "black", size = 10),
    legend.title     = element_text(color = "black", size = 11),
    legend.text      = element_text(color = "black", size = 10),
    panel.grid.minor = element_blank()
  )

# Panel A — raw counts
plot_counts <- ggplot(route_counts,
                      aes(x = factor(water_year), y = count,
                          fill = route_label)) +
  geom_col(width = 0.7, position = "stack") +
  scale_fill_brewer(palette = "Set2", name = "Migration route") +
  labs(x = NULL, y = "Number of fish", tag = "A") +
  shared_theme

# Panel B — proportions
plot_props <- ggplot(route_props,
                     aes(x = factor(water_year), y = proportion,
                         fill = route_label)) +
  geom_col(width = 0.7, position = "stack") +
  scale_fill_brewer(palette = "Set2", name = "Migration route") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Water year", y = "Proportion of fish", tag = "B") +
  shared_theme

# Combined panel
plot_counts / plot_props +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(figures, "gs_migration_routes.pdf"),
       width = 8, height = 8, units = "in", dpi = 300)

# =============================================================================
# 14. MIGRATION TIMING — ARRIVAL AT BENICIA BRIDGE (Group 2)
# =============================================================================

first_det_group2 <- upstream_migrants %>%
  filter(group == 2,
         water_year >= 2007 & water_year <= 2022) %>%
  arrange(animal_id, water_year, first_detection) %>%
  group_by(animal_id, water_year) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    first_date  = as.Date(first_detection),
    first_doy   = yday(first_date),
    dummy_date  = as.Date(first_doy - 1, origin = "2000-01-01"),
    outcome     = if_else(status == "up_complete", "Complete", "Incomplete")
  )

ggplot(first_det_group2, aes(x = factor(water_year), y = dummy_date)) +
  geom_boxplot(fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(aes(color = outcome), width = 0.15, size = 1.5, alpha = 0.8) +
  scale_color_manual(
    values = c("Complete" = "#2166ac", "Incomplete" = "#d73027"),
    name   = "Migration outcome"
  ) +
  scale_y_date(date_labels = "%b") +
  theme_bw(base_family = "Times New Roman") +
  labs(
    x     = "Water Year",
    y     = "Migration Start (Month of Water Year)",
    title = "Timing of First Benicia Bridge Detection by Water Year"
  ) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1,
                                    family = "Times New Roman"),
    axis.text.y      = element_text(family = "Times New Roman"),
    axis.title       = element_text(family = "Times New Roman"),
    legend.text      = element_text(family = "Times New Roman"),
    legend.title     = element_text(family = "Times New Roman"),
    panel.grid.minor = element_blank()
  )

# =============================================================================
# 15. PREDICTOR VARIABLE CORRELATION SCREENING
# =============================================================================
# Spearman correlation — appropriate for ordinal + continuous predictors.
# Variables with |r| > 0.8 not entered in same BRT/GLM.

predictors_HMC_full <- read.csv(
  file.path(flow, "HMC_variables_final_FULL.csv")
) 

predictors_VON_full <- read.csv(
  file.path(flow, "VON_annual_flow_result_updated_FULL.csv")
) %>% rename(water_year = Year)

# Correlation matrix — HMC
cor_matrix_HMC <- cor(
  predictors_HMC_full %>% dplyr::select(-water_year),
  method = "spearman", use = "pairwise.complete.obs"
)

# Subset to just the 7 final variables
pred_vars <- c("FA_pulse_dif", "Peak_Fre_2", "SP_ROC", "SP_Dur", 
               "SP_Mag", "Wet_BFL_Mag_50", "Wet_Tim")

# HMC correlation matrix
cor_matrix_HMC <- cor(
  predictors_HMC %>% dplyr::select(all_of(pred_vars)),
  method = "spearman",
  use = "complete.obs"
)

# VON correlation matrix
cor_matrix_VON <- cor(
  predictors_VON %>% dplyr::select(all_of(pred_vars)),
  method = "spearman",
  use = "complete.obs"
)

# Plot HMC
corrplot(cor_matrix_HMC, method = "color", type = "lower",
         col         = colorRampPalette(c("blue", "white", "red"))(200),
         tl.col      = "black", tl.srt = 45, addCoef.col = "black",
         title       = "Spearman Correlation — HMC Predictors",
         mar         = c(0, 0, 2, 0))

# Plot VON
corrplot(cor_matrix_VON, method = "color", type = "lower",
         col         = colorRampPalette(c("blue", "white", "red"))(200),
         tl.col      = "black", tl.srt = 45, addCoef.col = "black",
         title       = "Spearman Correlation — VON Predictors",
         mar         = c(0,0,2,0))

# Heatmap version
melt(cor_matrix_HMC) %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0) +
  theme_minimal(base_family = "Times New Roman") +
  labs(title = "Spearman Correlation Heatmap — HMC", fill = "Correlation") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# =============================================================================
# 16. BOOSTED REGRESSION TREES (BRT)
# =============================================================================
# Workflow:
#   Step 1 — Load and select final 7 predictor variables
#   Step 2 — Build BRT input datasets (one row per fish-year)
#   Step 3 — Fit full models with set.seed(33)
#   Step 4 — Inspect relative influence, note top predictors
#   Step 5 — Multi-seed stability analysis (10 seeds)
#   Step 6 — Partial dependence plots

library(dismo)
library(patchwork)

# --- 16a. Select final 7 predictor variables ---
pred_vars <- c("FA_pulse_dif", "Wet_BFL_Mag_50", "Wet_Tim",
               "SP_Mag", "SP_ROC", "SP_Dur", "Peak_Fre_2")

predictors_HMC <- read.csv(file.path(flow, "HMC_variables_final_FULL.csv")) %>%
  dplyr::select(water_year, all_of(pred_vars))

predictors_VON <- read.csv(file.path(flow, "VON_annual_flow_result_updated_FULL.csv")) %>%
  rename(water_year = Year) %>%
  dplyr::select(water_year, all_of(pred_vars))

# --- 16b. Build BRT input datasets ---

# Route selection
route_base <- glm_data_route %>%
  filter(!is.na(route_cluster)) %>%
  dplyr::select(animal_id, water_year, route_cluster)

final_data_route_HMC <- route_base %>%
  left_join(predictors_HMC, by = "water_year") %>%
  filter(complete.cases(.))

final_data_route_VON <- route_base %>%
  left_join(predictors_VON, by = "water_year") %>%
  filter(complete.cases(.))

# Migration duration
duration_base <- glm_data_route %>%
  filter(!is.na(migration_duration)) %>%
  dplyr::select(animal_id, water_year, migration_duration)

final_data_duration_HMC <- duration_base %>%
  left_join(predictors_HMC, by = "water_year") %>%
  filter(complete.cases(.))

final_data_duration_VON <- duration_base %>%
  left_join(predictors_VON, by = "water_year") %>%
  filter(complete.cases(.))

# Migration success — note: waiting for 2018-2022 flow metrics
success_base <- glm_data_success %>%
  filter(!is.na(migration_binary)) %>%
  dplyr::select(animal_id, water_year, migration_binary)

final_data_success_HMC <- success_base %>%
  left_join(predictors_HMC, by = "water_year") %>%
  filter(complete.cases(.))

final_data_success_VON <- success_base %>%
  left_join(predictors_VON, by = "water_year") %>%
  filter(complete.cases(.))

cat("Route HMC n:", nrow(final_data_route_HMC), "\n")
cat("Route VON n:", nrow(final_data_route_VON), "\n")
cat("Duration HMC n:", nrow(final_data_duration_HMC), "\n")
cat("Duration VON n:", nrow(final_data_duration_VON), "\n")
cat("Success HMC n:", nrow(final_data_success_HMC), "\n")
cat("Success VON n:", nrow(final_data_success_VON), "\n")

# --- 16c. Helper functions ---

summarize_brt <- function(model, name) {
  null_dev      <- model$self.statistics$mean.null
  resid_dev     <- model$self.statistics$mean.resid
  dev_explained <- (null_dev - resid_dev) / null_dev * 100
  cv_deviance   <- mean(model$cv.statistics$deviance.mean, na.rm = TRUE)
  
  cat("=====", name, "=====\n")
  cat("N trees:            ", model$n.trees, "\n")
  cat("Deviance explained: ", round(dev_explained, 2), "%\n")
  cat("Mean CV deviance:   ", round(cv_deviance, 4), "\n\n")
}

plot_partial_dependence <- function(model, vars_to_plot,
                                    y_label, family = "bernoulli") {
  pd_list <- list()
  for (var in vars_to_plot) {
    pd <- gbm::plot.gbm(model, i.var = var, return.grid = TRUE)
    pd <- pd %>%
      mutate(
        fitted    = if (family == "bernoulli") exp(y) / (1 + exp(y)) else y,
        predictor = var
      )
    pd_list[[var]] <- pd
  }
  pd_all <- bind_rows(pd_list)
  
  plots <- lapply(vars_to_plot, function(var) {
    pd_sub <- pd_all %>% filter(predictor == var)
    ggplot(pd_sub, aes(x = .data[[var]], y = fitted)) +
      geom_line(color = "steelblue", linewidth = 0.8) +
      labs(x = var, y = y_label) +
      theme_minimal(base_family = "Times New Roman") +
      theme(axis.title = element_text(size = 11),
            axis.text  = element_text(size = 9))
  })
  
  wrap_plots(plots) +
    plot_layout(guides = "collect") &
    theme(axis.title.y = element_blank(),
          text = element_text(family = "Times New Roman"))
}

# =============================================================================
# 16d. FIT FULL MODELS — set.seed(33)
# =============================================================================
set.seed(33)

# Route selection
brt_route_HMC <- gbm.step(
  data            = as.data.frame(final_data_route_HMC),
  gbm.x           = which(names(final_data_route_HMC) %in% pred_vars),
  gbm.y           = which(names(final_data_route_HMC) == "route_cluster"),
  family          = "bernoulli",
  tree.complexity = 2, learning.rate = 0.001, bag.fraction = 0.75,
  n.folds = 5, n.minobsinnode = 3, verbose = TRUE, plot.main = FALSE
)

brt_route_VON <- gbm.step(
  data            = as.data.frame(final_data_route_VON),
  gbm.x           = which(names(final_data_route_VON) %in% pred_vars),
  gbm.y           = which(names(final_data_route_VON) == "route_cluster"),
  family          = "bernoulli",
  tree.complexity = 2, learning.rate = 0.001, bag.fraction = 0.75,
  n.folds = 5, n.minobsinnode = 3, verbose = TRUE, plot.main = FALSE
)

# Migration duration
brt_duration_HMC <- gbm.step(
  data            = as.data.frame(final_data_duration_HMC),
  gbm.x           = which(names(final_data_duration_HMC) %in% pred_vars),
  gbm.y           = which(names(final_data_duration_HMC) == "migration_duration"),
  family          = "gaussian",
  tree.complexity = 2, learning.rate = 0.001, bag.fraction = 0.75,
  n.folds = 5, n.minobsinnode = 3, verbose = TRUE, plot.main = FALSE
)

brt_duration_VON <- gbm.step(
  data            = as.data.frame(final_data_duration_VON),
  gbm.x           = which(names(final_data_duration_VON) %in% pred_vars),
  gbm.y           = which(names(final_data_duration_VON) == "migration_duration"),
  family          = "gaussian",
  tree.complexity = 2, learning.rate = 0.001, bag.fraction = 0.75,
  n.folds = 5, n.minobsinnode = 3, verbose = TRUE, plot.main = FALSE
)

# Migration success
brt_success_HMC <- gbm.step(
  data            = as.data.frame(final_data_success_HMC),
  gbm.x           = which(names(final_data_success_HMC) %in% pred_vars),
  gbm.y           = which(names(final_data_success_HMC) == "migration_binary"),
  family          = "bernoulli",
  tree.complexity = 2, learning.rate = 0.001, bag.fraction = 0.75,
  n.folds = 5, n.minobsinnode = 3, verbose = TRUE, plot.main = FALSE
)

brt_success_VON <- gbm.step(
  data            = as.data.frame(final_data_success_VON),
  gbm.x           = which(names(final_data_success_VON) %in% pred_vars),
  gbm.y           = which(names(final_data_success_VON) == "migration_binary"),
  family          = "bernoulli",
  tree.complexity = 2, learning.rate = 0.001, bag.fraction = 0.75,
  n.folds = 5, n.minobsinnode = 3, verbose = TRUE, plot.main = FALSE
)

# Relative influence tables
cat("\n--- Route HMC ---\n");    summary(brt_route_HMC)
cat("\n--- Route VON ---\n");    summary(brt_route_VON)
cat("\n--- Success HMC ---\n");  summary(brt_success_HMC)
cat("\n--- Success VON ---\n");  summary(brt_success_VON)
cat("\n--- Duration HMC ---\n"); summary(brt_duration_HMC)
cat("\n--- Duration VON ---\n"); summary(brt_duration_VON)

# --- 16e. Multi-seed stability analysis ---
seeds <- c(33, 42, 123, 456, 789, 1000, 1234, 5678, 9999, 2024)

model_specs <- list(
  list(name     = "Route — HMC",
       data     = final_data_route_HMC,
       response = "route_cluster",
       family   = "bernoulli"),
  list(name     = "Route — VON",
       data     = final_data_route_VON,
       response = "route_cluster",
       family   = "bernoulli"),
  list(name     = "Success — HMC",
       data     = final_data_success_HMC,
       response = "migration_binary",
       family   = "bernoulli"),
  list(name     = "Success — VON",
       data     = final_data_success_VON,
       response = "migration_binary",
       family   = "bernoulli"),
  list(name     = "Duration — HMC",
       data     = final_data_duration_HMC,
       response = "migration_duration",
       family   = "gaussian"),
  list(name     = "Duration — VON",
       data     = final_data_duration_VON,
       response = "migration_duration",
       family   = "gaussian")
)

all_results <- list()

for (spec in model_specs) {
  cat("Running:", spec$name, "\n")
  seed_results <- list()
  
  for (s in seeds) {
    set.seed(s)
    
    model <- tryCatch(
      gbm.step(
        data            = as.data.frame(spec$data),
        gbm.x           = which(names(spec$data) %in% pred_vars),
        gbm.y           = which(names(spec$data) == spec$response),
        family          = spec$family,
        tree.complexity = 2, learning.rate = 0.001, bag.fraction = 0.75,
        n.folds         = 5, n.minobsinnode = 3,
        verbose         = FALSE, plot.main = FALSE
      ),
      error = function(e) {
        cat("  Error at seed", s, ":", conditionMessage(e), "\n")
        NULL
      }
    )
    
    if (is.null(model)) next
    
    null_dev      <- model$self.statistics$mean.null
    resid_dev     <- model$self.statistics$mean.resid
    dev_explained <- (null_dev - resid_dev) / null_dev * 100
    cv_dev        <- mean(model$cv.statistics$deviance.mean, na.rm = TRUE)
    auc <- if (spec$family == "bernoulli") {
      mean(model$cv.statistics$discrimination.mean, na.rm = TRUE)
    } else {
      NA
    }
    
    seed_results[[as.character(s)]] <- data.frame(
      model         = spec$name,
      seed          = s,
      dev_explained = dev_explained,
      cv_deviance   = cv_dev,
      auc           = auc,
      n_trees       = model$n.trees
    )
  }
  
  all_results[[spec$name]] <- bind_rows(seed_results)
  cat("  Done\n")
}

all_results_df <- bind_rows(all_results)

summary_table <- all_results_df %>%
  group_by(model) %>%
  summarise(
    mean_dev_explained  = round(mean(dev_explained, na.rm = TRUE), 2),
    range_dev_explained = paste0(round(min(dev_explained, na.rm = TRUE), 2),
                                 " - ",
                                 round(max(dev_explained, na.rm = TRUE), 2)),
    mean_cv_deviance    = round(mean(cv_deviance, na.rm = TRUE), 4),
    range_cv_deviance   = paste0(round(min(cv_deviance, na.rm = TRUE), 4),
                                 " - ",
                                 round(max(cv_deviance, na.rm = TRUE), 4)),
    mean_auc            = round(mean(auc, na.rm = TRUE), 3),
    range_auc           = paste0(round(min(auc, na.rm = TRUE), 3),
                                 " - ",
                                 round(max(auc, na.rm = TRUE), 3)),
    .groups = "drop"
  )

cat("\nBRT Performance Summary:\n")
print(summary_table)

write.csv(all_results_df,
          file.path(cleaned_data, "brt_seed_results_full.csv"),
          row.names = FALSE)
write.csv(summary_table,
          file.path(cleaned_data, "brt_summary_table.csv"),
          row.names = FALSE)

save.image("gs_analysis_workspace_2026.RData")

# --- 16f. Partial dependence plots ---
# --- 16f. Partial dependence plots ---
plot_partial_dependence_both <- function(model_HMC, model_VON, 
                                         y_label, family = "bernoulli",
                                         n_HMC = 3, n_VON = 3,
                                         layout = "vertical") {
  
  var_labels <- c(
    "FA_pulse_dif"   = "Fall Pulse Difference (cfs)",
    "Wet_BFL_Mag_50" = "Wet Season Median Magnitude (cfs)",
    "Wet_Tim"        = "Wet Season Timing (day of year)",
    "SP_Mag"         = "Spring Recession Magnitude (cfs)",
    "SP_ROC"         = "Spring Rate of Change (proportion/day)",
    "SP_Dur"         = "Spring Duration (days)",
    "Peak_Fre_2"     = "Peak Flow Frequency (count)"
  )
  
  total_panels <- max(n_HMC, n_VON)
  
  make_row <- function(model, n_top) {
    top_n <- summary(model, plotit = FALSE) %>%
      arrange(desc(rel.inf)) %>%
      slice(1:n_top) %>%
      pull(var)
    
    pd_list <- list()
    for (var in top_n) {
      pd <- gbm::plot.gbm(model, i.var = var, return.grid = TRUE)
      pd <- pd %>%
        mutate(
          fitted    = if (family == "bernoulli") exp(y) / (1 + exp(y)) else y,
          predictor = var
        )
      pd_list[[var]] <- pd
    }
    pd_all <- bind_rows(pd_list)
    
    plots <- lapply(top_n, function(var) {
      pd_sub <- pd_all %>% filter(predictor == var)
      ggplot(pd_sub, aes(x = .data[[var]], y = fitted)) +
        geom_line(color = "steelblue", linewidth = 0.8) +
        labs(x = var_labels[var], y = y_label) +
        theme_minimal(base_family = "Times New Roman") +
        theme(
          axis.title   = element_text(size = 10, family = "Times New Roman"),
          axis.text    = element_text(size = 8,  family = "Times New Roman"),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank()
        )
    })
    
    n_blank <- total_panels - n_top
    if (n_blank > 0) {
      blank <- replicate(n_blank, plot_spacer(), simplify = FALSE)
      plots <- c(plots, blank)
    }
    
    wrap_plots(plots, nrow = 1)
  }
  
  row_HMC <- make_row(model_HMC, n_HMC)
  row_VON <- make_row(model_VON, n_VON)
  
  # Layout switch
  if (layout == "horizontal") {
    combined <- (row_HMC | row_VON)
  } else {
    combined <- (row_HMC / row_VON)
  }
  
  combined +
    plot_annotation(
      tag_levels = "a",
      tag_suffix = ")",
      theme = theme(
        plot.tag = element_text(size = 12, face = "bold",
                                family = "Times New Roman")
      )
    )
}

# Route — top 3 for both, vertical (stacked)
png(file.path(figures, "route_brt_plots.png"),
    width = 3600, height = 2400, res = 300)
plot_partial_dependence_both(brt_route_HMC, brt_route_VON,
                             "Probability of Delta Route",
                             family = "bernoulli",
                             n_HMC = 3, n_VON = 3,
                             layout = "vertical")
dev.off()

# Success — top 1 each, horizontal (side by side)
png(file.path(figures, "success_brt_plots.png"),
    width = 2400, height = 1200, res = 300)
plot_partial_dependence_both(brt_success_HMC, brt_success_VON,
                             "Probability of Migration Success",
                             family = "bernoulli",
                             n_HMC = 1, n_VON = 1,
                             layout = "horizontal")
dev.off()

# Duration — top 1 each, horizontal (side by side)
png(file.path(figures, "duration_brt_plots.png"),
    width = 2400, height = 1200, res = 300)
plot_partial_dependence_both(brt_duration_HMC, brt_duration_VON,
                             "Migration Duration (days)",
                             family = "gaussian",
                             n_HMC = 1, n_VON = 1,
                             layout = "horizontal")
dev.off()


# =============================================================================
# 17. SAVE WORKSPACE
# =============================================================================


