# sturgeon_functional_flow
Code repository for PhD work with green sturgeon and functional flow movement analysis
## Data
Raw detection data and intermediate cleaned datasets are available upon request from the corresponding author. The analysis-ready datasets provided here are sufficient to reproduce all GLM and BRT analyses presented in the manuscript."
Raw acoustic telemetry detections were queried from the Pacific Animal Telemetry Hub (PATH) database. The following receiver arrays were used:

- **Golden Gate / Bay** (group 1) — San Francisco Bay estuary entry
- **Benicia / Carquinez** (group 2) — tidal freshwater entry
- **Groups 3-23** — Sacramento River and Delta intermediate receivers
- **Group 24** — Spawning ground receivers (upper Sacramento River)

Fish were classified as `up_complete` if detected at group 24 (spawning grounds) and `up_incomplete` if detected at Benicia but not reaching group 24. Route selection analysis was restricted to 2007-2017 when Delta-discriminating receivers (Georgiana Slough, Moki Delta Cross) were active. Migration success analysis extended through 2022.

Discharge data were downloaded from the California Data Exchange Center (CDEC):
- Hamilton City gauge (HMC, sensor 41)
- Verona gauge (VON, sensor 41)
- Rio Vista temperature (RVB, sensor 25)

Functional flow metrics were calculated following Yarnell et al. (2020) for water years 2007-2022.

## Analysis Pipeline

### detection_cleaning_gsflow.R
- Queries and filters detections from PATH database
- Assigns receiver groups (1-24) based on spatial location
- Classifies migration status (up_complete, up_incomplete, down_complete, etc.)
- Removes double-tagged fish and applies manual corrections for validated migrants with detection gaps
- Outputs: `events_with_receivergroups_032026.csv`

### gs_analysis_clean_final.R
- Builds route strings and classifies Delta vs. mainstem Sacramento routes
- Calculates mean discharge during each fish's migration window (Hamilton City and Verona gauges)
- Calculates mean water temperature one week prior to Benicia Bridge detection
- Fits binomial GLMs for route selection and migration success
- Fits Gaussian GLMs for migration duration
- Fits boosted regression trees (BRTs) for all three response variables
- Generates partial dependence plots and summary figures

## Dependencies

```r
# Core packages
library(tidyverse)
library(lubridate)
library(glatos)      # GLATOS detection filtering

# Spatial
library(sf)
library(tigris)
library(RANN)        # Nearest-neighbor group assignment

# Modeling
library(dismo)       # BRT via gbm.step
library(patchwork)   # Figure composition

# Data download
# CDEC gauge data downloaded via custom cdec_query() function in 02_analysis.R
```

Install all packages with:
```r
install.packages(c("tidyverse", "lubridate", "sf", "tigris", 
                   "RANN", "dismo", "patchwork"))
```

## Citation

Tracy, E.E., Walter, J., Yarnell, S., Rypel, A.L., and Fangue, N.A. (*in prep*). Green sturgeon migratory movements are related to flow variations in the Sacramento-San Joaquin Delta.

## Contact

Erin E. Tracy — eetracy@ucdavis.edu

## License

This repository is made available for reproducibility purposes. Please contact the corresponding author before using data or code for purposes beyond reproducing the analyses presented in the associated manuscript.
