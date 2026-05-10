#==============================================================================
# GREEN STURGEON DETECTION CLEANING
# Script: 01_detection_cleaning.R
# Author: Erin Tracy
# Last updated: April 2026
#
# PURPOSE:
# Download raw green sturgeon detections from PATH database, run through
# GLATOS false detection filter, reduce to detection events, and assign
# receiver group labels from ArcGIS metadata.
#
# INPUTS:
#   - PATH database (requires UCD credentials and campus VPN)
#   - arc_receivers_update.csv: receiver metadata with group labels from ArcGIS
#
# OUTPUTS:
#   - events_with_receivergroups_032026.csv: cleaned detection events with
#     receiver group labels, water year, and migration status
#
# NOTE: Database query requires connection to UCD VPN if running remotely.
#   Raw detection file (green_sturgeon_detection_OTN_012926.csv) can be used
#   to skip the database query if already downloaded.
#==============================================================================
library(data.table)
library(lubridate)
library(tidyverse)
library(RPostgres)
library(DBI)
library(dbplyr)
library(glatos)

#install.packages('glatos', repos = c('https://ocean-tracking-network.r-universe.dev', 'https://cloud.r-project.org'))

#==============================================================================
# SECTION 1: DOWNLOAD DETECTIONS FROM PATH DATABASE
# Skip this section if green_sturgeon_detection_OTN_012926.csv already exists
#==============================================================================

# Connect to PATH database - requires UCD credentials in .Renviron file

readRenviron("./.Renviron")
con <- DBI::dbConnect(
  RPostgres::Postgres(),
  host = Sys.getenv("dbhost"),
  user = Sys.getenv("dbuser"),
  password = Sys.getenv("dbpasswd"),
  dbname = "pathnode"
)

#THIS WORKS PERFECT FOR RUNNING DETECTIONS THROUGH GLATOS
query <- "
select 
det.collectioncode AS detectedby,
anm.collectioncode AS collectioncode,
anm.catalognumber AS catalognumber,
anm.scientificname AS scientificname,
anm.commonname AS commonname,
det.datelastmodified AS datelastmodified,
det.collectioncode :: text AS detectedby,
det.collectioncode :: text AS receiver_group,
CASE WHEN det.collectornumber ::text ~~ '%(%' :: text THEN (det.station || '(' :: text) || split_part(det.collectornumber :: text, '(' :: text, 2) ELSE det.station END AS station,
det.collectornumber :: text AS receiver,
det.bottom_depth,
det.receiver_depth,
det.fieldnumber :: text AS tagname,
CASE WHEN det.fieldnumber :: text ~~ 'A%' :: text THEN (
  split_part(det.fieldnumber :: text, '-' :: text, 1) || '-' :: text
) || split_part(det.fieldnumber :: text, '-' :: text, 2) ELSE NULL :: text END AS codespace,
det.sensorname :: text AS sensorname,
det.sensorraw :: text AS sensorraw,
CASE WHEN det.sensorraw IS NULL THEN 'pinger' :: text ELSE det.sensortype :: text END AS sensortype,
det.sensorvalue :: numeric AS sensorvalue,
det.sensorunit :: text AS sensorunit,
det.datecollected,
'UTC' :: text AS timezone,
round(det.longitude :: numeric, 5) AS longitude,
round(det.latitude :: numeric, 5) AS latitude,
det.the_geom,
det.yearcollected :: text AS yearcollected,
det.monthcollected :: text AS monthcollected,
det.daycollected :: text AS daycollected,
det.julianday,
det.timeofday,
local_area,
det.notes,
citation,
(det.collectioncode :: text || '-' :: text) || det.catalognumber :: text AS unqdetecid 
from  (SELECT * FROM obis.otn_detections_early
       UNION ALL
       SELECT * FROM obis.otn_detections_2007
       UNION ALL
       SELECT * FROM obis.otn_detections_2008
       UNION ALL
       SELECT * FROM obis.otn_detections_2009
       UNION ALL
       SELECT * FROM obis.otn_detections_2010
       UNION ALL
       SELECT * FROM obis.otn_detections_2011
       UNION ALL
       SELECT * FROM obis.otn_detections_2012
       UNION ALL
       SELECT * FROM obis.otn_detections_2013
       UNION ALL
       SELECT * FROM obis.otn_detections_2014
       UNION ALL
       SELECT * FROM obis.otn_detections_2015
       UNION ALL
       SELECT * FROM obis.otn_detections_2016
       UNION ALL
       SELECT * FROM obis.otn_detections_2017
       UNION ALL
       SELECT * FROM obis.otn_detections_2018
       UNION ALL
       SELECT * FROM obis.otn_detections_2019
       UNION ALL
       SELECT * FROM obis.otn_detections_2020
       UNION ALL
       SELECT * FROM obis.otn_detections_2021
       UNION ALL
       SELECT * FROM obis.otn_detections_2022
       UNION ALL
       SELECT * FROM obis.otn_detections_2023
       UNION ALL
       SELECT * FROM obis.otn_detections_2024
       UNION ALL
       SELECT * FROM obis.otn_detections_2025
) det
left JOIN obis.otn_animals anm ON det.relatedcatalogitem = anm.catalognumber 
left join obis.otn_resources res on det.collectioncode = res.collectioncode
where det.relationshiptype = 'ANIMAL'
and LOWER(anm.commonname) = 'green sturgeon'
order by det.datecollected"

# Execute the query and fetch the results
gs_detections_original <- dbGetQuery(con, query)
write.csv(gs_detections_original, "green_sturgeon_detection_OTN_012926.csv", row.names = FALSE)

#only run if need tag metadata
#obis otn_animals metadata to join with detections/receiver data 
query <- "
SELECT *
FROM obis.otn_animals
WHERE LOWER(commonname) = 'green sturgeon'
"
otn_animals <- dbGetQuery(con, query)
write.csv(otn_animals, "green_sturgeon_animalmetadata_OTN_020726.csv", row.names = FALSE)

animal_lookup <- otn_animals %>%
  dplyr::select(catalognumber, age, lifestage, length, lengthtype, sex) %>%
  distinct(catalognumber, .keep_all = TRUE)

events <- events %>%
  left_join(
    animal_lookup,
    by = c("animal_id" = "catalognumber")
  )

#==============================================================================
# SECTION 2: LOAD RAW DETECTIONS
# Start here if raw detection file already exists
#==============================================================================
gs_detections_original <- read.csv("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/raw_data/green_sturgeon_detection_OTN_012926/green_sturgeon_detection_OTN_012926.csv") 

#==============================================================================
# SECTION 3: GLATOS FALSE DETECTION FILTER
#==============================================================================

#Need to rename columns to run through GLATOS
colnames(gs_detections_original)[c(3, 8, 10, 13, 14, 18, 19, 20, 22, 23)] <- c("animal_id", "glatos_array", "receiver_sn", "transmitter_id", "transmitter_codespace", "sensor_value", "sensor_unit", "detection_timestamp_utc", "deploy_long", "deploy_lat")
#making sure timezone write class
class(gs_detections_original$detection_timestamp_utc)
gs_detections_original$detection_timestamp_utc <-
  as.POSIXct(gs_detections_original$detection_timestamp_utc,format = "%Y-%m-%d %H:%M:%S",
             tz = "UTC")

#False detection filter
#Write the filtered data to a new det_filtered object
#Doesn't delete rows, adds new column if detection was filtered out
detections_filtered <- false_detections(gs_detections_original, tf=3600, show_plot=TRUE)
#chapter 3 results: The filter identified 53505 (0.78%) of 6858330 detections as potentially false
#re-run 020126 The filter identified 55540 (0.8%) of 6965102 detections as potentially false

# Filter based on the column if you're happy with it.
detections_filtered <- detections_filtered[detections_filtered$passed_filter == 1,]
nrow(detections_filtered) # Check that its smaller than before


#==============================================================================
  # SECTION 4: REDUCE TO DETECTION EVENTS
  #==============================================================================
# Reduce Detections to Detection Events ####
events <- detection_events(detections_filtered,
                           location_col = 'station', # combines events across different receivers in a single array
                           time_sep=Inf)

#==============================================================================
# SECTION 5: REMOVE DOUBLE-TAGGED FISH
#==============================================================================
#now exclude 1 of the tags for double tagged fish
double_tagged_remove <- c(
  # Original double-tagged fish
  "UCDHIST-GS0576-2011-05-06",
  "UCDHIST-GS0775-2012-04-23",
  "UCDHIST-GS0787-2012-05-07",
  "UCDHIST-GS0937-2011-06-08",
  "UCDHIST-GS0939-2011-06-08",
  # GS9 prefix duplicates - same fish as GS0 counterparts
  "UCDHIST-GS9455-2010-05-04",
  "UCDHIST-GS9457-2010-05-06",
  "UCDHIST-GS9458-2010-05-08",
  "UCDHIST-GS9539-2010-07-20",
  "UCDHIST-GS9540-2010-07-23",
  "UCDHIST-GS9573-2011-05-06"
)

# ADD this line after the double_tagged_remove vector
events <- events %>%
  filter(!animal_id %in% double_tagged_remove)

cat("Fish remaining after removing double-tagged:", n_distinct(events$animal_id), "\n")

#==============================================================================
# SECTION 6: ADD WATER YEAR
#==============================================================================

events <- events %>%
  mutate(
    year       = year(first_detection),
    month      = month(first_detection),
    water_year = ifelse(month > 9, year + 1, year)
  )

###############################################################################################
#Only do first time when looking for fish to exclude
#Make Abacus plots to look for dropped tags, deaths, and duplicates 
# Ensure the date column is converted to Date format
events_plot <- events$last_detection <- as.Date(events$last_detection)
# Define the output PDF file
output_file <- "C:/Users/etracy1/Desktop/Backup/R_directory/ST_telemetry/abacus_plots_location_020226.pdf"

# Open a PDF device for multi-page output
pdf(output_file, width = 8, height = 6)  # Adjust width and height as needed

# Get unique animal IDs
unique_animals <- unique(events_plot$animal_id)

# Loop through each animal_id
for (animal in unique_animals) {
  # Filter data for the current animal
  animal_data <- events_plot %>%
    filter(animal_id == animal)
  
  # Create the abacus plot for the current animal
  p <- ggplot(animal_data, aes(x = last_detection, y = mean_latitude)) +
    geom_point(size = 3, color = "blue") +
    labs(
      title = paste("Abacus Plot for GS:", animal),
      x = "Last Detection",
      y = "Mean Latitude"
    ) +
    scale_x_date(
      date_breaks = "1 week",          # Display labels at weekly intervals
      date_labels = "%b %d '%y"       # Format as "Month Day 'Year" (e.g., Jan 01 '25)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate labels for readability
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
  # Print the plot to the current page in the PDF
  print(p)
}
# Close the PDF device
dev.off()
# Message to indicate the PDF is ready
cat("Abacus plots have been saved to", output_file, "\n")

#Filter
# Filtering animals based on abacus plots
events$animal_id <- as.character(events$animal_id)
# List of animal_ids to exclude

# Remove spaces
events$animal_id <- trimws(events$animal_id)
any(grepl("^\\s|\\s$", events$animal_id))

#animals excluded based on dropped tag potential (didn't analyze this time)
animal_ids_to_exclude <- c()
events <- events[!(events$animal_id %in% animal_ids_to_exclude), ]



###########################################################################################################
#convert from standard to fork length then get length at age

#have to assume fork length for NA length type for 48 fish

# Rebuild with CA-specific Von Bertalanffy
# Lt = 155.27 × (1 - e^(-0.125(t + 1.318)))
# Solve for t: t = -ln(1 - FL/155.27) / 0.125 - 1.318
#Ulaski and Quist et al. 2021
# updated ulaski to the fish base estimates because our fish were bigger

#i think some of these were entered wrong they are too small 
fish_lengths_final <- events_2026 %>%
  filter(animal_id %in% upstream_migrants$animal_id) %>%
  distinct(animal_id, length, lengthtype) %>%
  filter(!is.na(length)) %>%
  group_by(animal_id) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    length_cm = length * 100,
    lengthtype_assumed = if_else(is.na(lengthtype), "FORK", lengthtype),
    FL_cm = case_when(
      lengthtype_assumed == "FORK"  ~ length_cm,
      lengthtype_assumed == "TOTAL" ~ (length_cm + 4.6131) / 1.1374,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(FL_cm >= 60)

cat("Fish with plausible FL:", nrow(fish_lengths_final), "\n")
summary(fish_lengths_final$FL_cm)
#if i just exclude the 69 fish that have a questionable length we assumed was decimeters
metadata <- read.csv("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_flow_phd/raw_data/greensturgeon_metadata_path.csv")
# Extract Fish_ID from animal_id to join to metadata
small_fish <- events_2026 %>%
  filter(animal_id %in% upstream_migrants$animal_id) %>%
  distinct(animal_id, length, lengthtype) %>%
  filter(!is.na(length)) %>%
  group_by(animal_id) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    length_cm = length * 100,
    Fish_ID = str_extract(animal_id, "GS\\d+")
  ) %>%
  filter(length_cm < 60) %>%
  left_join(metadata, by = "Fish_ID")

print(small_fish, n = 70)

write.csv(small_fish, "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_flow_phd/raw_data/greensturgeon_smallfish_path.csv")


# Then I did a bunch of exploring the receivers and abacus plots for time to make the groups
# Receivers were grouped to ensure no time gaps and confluences were monitored
# Can be found in gs_analysis.R code the file is receiver_group_final.csv

###################################################
#sort into up and down migrants 
#Battalie 2024 used Approximate Coordinates (RM 200 / 322 km) which is approximately hamilton city 
# slightly south at lat -122 and long 39.7
# our lowest point for group 24 is -121.9745, 39.73131

#==============================================================================
# SECTION 7: ASSIGN RECEIVER GROUP LABELS FROM ARCGIS METADATA
#==============================================================================

# Load receiver metadata - authoritative source for receiver group labels
# Generated from ArcGIS spatial join of receiver locations to route polygons
receiver_metadata <- read.csv(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/cleaned_data/arc_receivers_update.csv"
)

# Build location to receiver_group lookup from metadata
# DECISION: Use receiver_metadata as authoritative source for group labels
# rather than the collectioncode from the database query, which is often
# the tagging project code rather than a meaningful spatial group
correct_groups <- receiver_metadata %>%
  distinct(relatedcatalogitem, receiver_group) %>%
  filter(!is.na(relatedcatalogitem), relatedcatalogitem != "",
         !is.na(receiver_group), receiver_group != "") %>%
  rename(location = relatedcatalogitem,
         receiver_group_correct = receiver_group) %>%
  # Fix known conflicts identified during review
  mutate(receiver_group_correct = case_when(
    location == "SR_FEATHER2_RT" ~ "sacramento",
    location == "RICHBR_22_2015" ~ "bay",
    TRUE ~ receiver_group_correct
  )) %>%
  distinct(location, .keep_all = TRUE)

group <- read.csv("C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/cleaned_data/events_with_receivergroups_032026.csv")
#group was in a different receiver metadata file, it can be joined here

# Remove unused event column from GLATOS output
events <- events %>% dplyr::select(-event)

cat("Unique receiver groups assigned:", length(unique(events$receiver_group)), "\n")
cat("Unique fish:", n_distinct(events$animal_id), "\n")
cat("Water years covered:", range(events$water_year, na.rm = TRUE), "\n")

#==============================================================================
# SECTION 8: SAVE CLEANED EVENTS
#==============================================================================

write.csv(
  events,
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/cleaned_data/events_with_receivergroups_032026.csv",
  row.names = FALSE
)

cat("Cleaned events saved to events_with_receivergroups_032026.csv\n")
#==============================================================================
# GREEN STURGEON MIGRATION STATUS CLASSIFICATION

# PURPOSE:
# Classify each unique animal_id x water_year combination into migration
# status categories based on detection history relative to key receiver arrays.
# Also joins status back to events for use in downstream analysis.
#
# INPUTS:
#   - events_with_receivergroups_032026.csv: cleaned detection events
#     (output from 01_detection_cleaning.R)
#
# OUTPUTS:
#   - migration_status: data frame with animal_id, water_year, status
#   - events: events data frame with status column added
#   - events_with_receivergroups_032026.csv: updated with status column
# MIGRATION STATUS CATEGORIES:
#   up_complete:     fish detected at ocean side receivers (golden_gate or bay)
#                    THEN benicia/carquinez THEN spawning ground — confirmed
#                    full upstream spawning migration
#   up_incomplete:   fish detected at ocean side receivers THEN benicia/carquinez
#                    but never reached spawning ground — failed or incomplete
#                    upstream migration
#   down_complete:   fish detected moving downstream from spawning ground to
#                    golden gate — confirmed full downstream migration
#   down_incomplete: fish started downstream migration but did not reach ocean
#   incomplete_dead: fish confirmed or likely dead based on detection history
#                    review — labeled for specific water years only
#   bad:             ambiguous movement, delta resident, or data issues
#==============================================================================

#==============================================================================
# SECTION 1: LOAD DATA
#==============================================================================

events <- read.csv(
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/cleaned_data/events_with_receivergroups_032026.csv"
)

events <- events %>%
  mutate(first_detection = as.POSIXct(first_detection, tz = "UTC"),
         last_detection  = as.POSIXct(last_detection,  tz = "UTC"))

#==============================================================================
# SECTION 2: CLASSIFY MIGRATION STATUS
#==============================================================================

# DECISION: Use benicia OR carquinez as migration entry point
# RATIONALE: Both receiver arrays are at the Carquinez Strait — the narrow
# passage between San Pablo Bay and Suisun Bay. Any fish detected at either
# array has entered the freshwater migration corridor. Using both ensures
# maximum detection coverage given occasional receiver gaps at either array.

migration_status <- events %>%
  # Step 1: Get all fish that passed benicia or carquinez
  filter(receiver_group %in% c("benicia", "carquinez")) %>%
  distinct(animal_id, water_year) %>%
  # Step 2: Join summary of detection history for each fish x water_year
  left_join(
    events %>%
      arrange(animal_id, water_year, first_detection) %>%
      group_by(animal_id, water_year) %>%
      dplyr::summarise(
        lifestage  = first(lifestage),
        first_group = first(receiver_group),
        last_group  = last(receiver_group),
        # Earliest detection at each key receiver group
        gg  = suppressWarnings(min(first_detection[receiver_group == "golden_gate"],  na.rm = TRUE)),
        bay = suppressWarnings(min(first_detection[receiver_group == "bay"],           na.rm = TRUE)),
        bc  = suppressWarnings(min(first_detection[receiver_group %in% c("benicia", "carquinez")], na.rm = TRUE)),
        sac = suppressWarnings(min(first_detection[receiver_group == "sacramento"],    na.rm = TRUE)),
        sg  = suppressWarnings(min(first_detection[receiver_group == "spawning_ground"], na.rm = TRUE)),
        .groups = "drop"
      ),
    by = c("animal_id", "water_year")
  ) %>%
  mutate(
    # Boolean flags for detection at each key receiver group
    has_gg  = !is.na(gg)  & !is.infinite(gg),
    has_bay = !is.na(bay) & !is.infinite(bay),
    has_bc  = !is.na(bc)  & !is.infinite(bc),
    has_sac = !is.na(sac) & !is.infinite(sac),
    has_sg  = !is.na(sg)  & !is.infinite(sg),
    # Ocean side receiver detected before benicia/carquinez
    ocean_before_bc = (has_gg & gg < bc) | (has_bay & bay < bc),
    # Spawning ground detected before benicia/carquinez (downstream migrant)
    sg_before_bc = has_sg & sg < bc,
    # Classify migration status
    status = case_when(
      # DOWNSTREAM COMPLETE: sg then bc then gg - full downstream migration
      sg_before_bc & has_bc & has_gg & bc < gg              ~ "down_complete",
      # DOWNSTREAM COMPLETE: sac then bc then gg - fish may have spawned
      # downstream of spawning_ground receivers
      has_sac & sac < bc & has_gg & bc < gg                 ~ "down_complete",
      # DOWNSTREAM INCOMPLETE: sg then bc but never reached gg
      sg_before_bc & has_bc & !has_gg                       ~ "down_incomplete",
      # UPSTREAM COMPLETE: ocean side then bc then sg
      ocean_before_bc & has_sg & bc < sg                    ~ "up_complete",
      # UPSTREAM INCOMPLETE: ocean side then bc but never reached sg
      ocean_before_bc & !has_sg                             ~ "up_incomplete",
      # BAD: ambiguous, delta resident, or data issues
      TRUE                                                   ~ "bad"
    )
  )

cat("Migration status counts:\n")
print(migration_status %>% count(status))
#==============================================================================
# SECTION 3: MANUAL CORRECTIONS
#==============================================================================

# UCDHIST-GS0823: anomalous detection pattern suggesting tag shed or error
# CDFWA15-1219838: detection sequence inconsistent with upstream migration
# Manual corrections based on individual detection history review
migration_status <- migration_status %>%
  mutate(status = case_when(
    animal_id == "CDFWA15-1219838-2016-09-08" & water_year == 2017 ~ "bad",
    animal_id == "UCDHIST-GS0823-2012-07-03"  & water_year == 2021 ~ "bad",
    # GS0512 WY2015: interior_delta -> benicia -> bay movement is downstream
    # migration not upstream - fish never came from ocean side before Benicia
    animal_id == "UCDHIST-GS0512-2012-04-20"  & water_year == 2015 ~ "bad",
    TRUE ~ status
  ))
#==============================================================================
# SECTION 4: LABEL INCOMPLETE_DEAD FISH
#==============================================================================

# These fish are labeled incomplete_dead for specific water years only
# based on individual detection history review indicating likely mortality.
# Multi-year tagged fish may have been alive in other years so we only
# label the specific water year of likely mortality, not all years.
dead_fish_years <- data.frame(
  animal_id = c(
    "UCDHIST-GS0276-2005-08-20",
    "UCDHIST-GS0488-2011-08-10",
    "UCDHIST-GS0516-2012-04-26",
    "UCDHIST-GS0821-2012-07-02",
    "UCDHIST-GS0814-2012-07-01",
    "CDFWA15-1306970-2018-12-18"
  ),
  water_year = as.integer(c(2008, 2016, 2012, 2014, 2020, 2020))
)

migration_status <- migration_status_clean %>%
  left_join(dead_fish_years %>% mutate(is_dead = TRUE),
            by = c("animal_id", "water_year")) %>%
  mutate(status = ifelse(!is.na(is_dead), "incomplete_dead", status)) %>%
  dplyr::select(-is_dead)

# Verify GS0488
migration_status %>%
  filter(animal_id == "UCDHIST-GS0488-2011-08-10") %>%
  dplyr::select(animal_id, water_year, status)

migration_status <- migration_status %>%
  mutate(status = case_when(
    # GS0303 2010: valid upstream migration entering via group 4 after
    # failed first attempt in Nov 2009. Reached spawning grounds March 2010.
    animal_id == "UCDHIST-GS0303-2006-06-28" & water_year == 2010 ~ "up_complete",
    TRUE ~ status
  ))

# Verify overall counts
migration_status %>% count(status)
write.csv(
  migration_status,
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/cleaned_data/migratory_status_03162026.csv",
  row.names = FALSE
)

#==============================================================================
# SECTION 5: JOIN STATUS BACK TO EVENTS
#==============================================================================

# Join status back to events - ensure one status per fish x water_year
events <- events %>%
  dplyr::select(-any_of("status")) %>%
  left_join(
    migration_status %>% 
      dplyr::select(animal_id, water_year, status),
    by = c("animal_id", "water_year")
  )

# Verify no duplicates
duplicate_check <- events %>%
  distinct(animal_id, water_year, status) %>%
  group_by(animal_id, water_year) %>%
  filter(n() > 1)

cat("Fish x water_year with multiple statuses:", nrow(duplicate_check), "\n")



#==============================================================================
# SECTION 6: SAVE
#==============================================================================

# Save updated events with status
write.csv(
  events,
  "C:/Users/eetracy/Desktop/R_directory/ST_telemetry/gs_multistate/cleaned_data/events_with_receivergroups_032026.csv",
  row.names = FALSE
)

cat("\nEvents with status saved to events_with_receivergroups_032026.csv\n")
cat("migration_status object ready for use in 04_multistate_data_prep.R\n")


#exploring age as a covariate
#perry 2018 included length as a covariate 
# Get one length per fish (some fish have multiple rows)
fish_lengths <- events %>%
  distinct(animal_id, length, age) %>%
  filter(!is.na(length)) %>%
  # If a fish has multiple length records take the first/most common
  group_by(animal_id) %>%
  dplyr::summarise(
    length_m = first(length),
    age_at_tag = first(age),
    .groups = "drop"
  )

# Join to model fish
model_fish_size <- detection_history %>%
  dplyr::select(animal_id, water_year) %>%
  left_join(fish_lengths, by = "animal_id") %>%
  mutate(
    # Approximate age at migration
    tag_year = as.integer(substr(animal_id, 
                                 nchar(animal_id)-9, 
                                 nchar(animal_id)-6)),
    years_since_tag = water_year - tag_year
  )

# Summary
summary(model_fish_size$length_m)
summary(model_fish_size$years_since_tag)
hist(model_fish_size$years_since_tag, 
     main = "Years between tagging and migration",
     xlab = "Years since tagging")
hist(model_fish_size$length_m,
     main = "Fork length at tagging (m)",
     xlab = "Fork length (m)")
#lengths are all over the place (0.4 to 2) because of juveniles and time since taggin 
#to migration is long 4-9 years old
