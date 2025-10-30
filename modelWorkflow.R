library(tidyverse)
library(sf)
sf_use_s2(FALSE)
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("OptimSPD/src/SdpMigC.cpp")
source("OptimSPD/R/OptimSPD.R")
source("functions/fnc.R")

### Grid ###
hex <- st_read("inputData/hexagonalGrid_100km_laea_EuropeNoWater.gpkg")

europe      <- rnaturalearth::ne_countries(scale = "medium", continent = "europe", returnclass = "sf")
europe_laea <- st_transform(europe, st_crs(hex)) %>% st_intersection(st_bbox(hex) %>% st_as_sfc())

ggplot() +
  geom_sf(data = europe_laea) +
  geom_sf(data = hex, fill = NA)


### Environment (merge hexagon properties with species wintering nd breeding ranges) ###
siteEnv  <- readRDS("./inputData/EnvTab_forModel.rds") %>% st_as_sf()

# plot(siteEnv %>% filter(yr == 2016) %>% group_by(yr, hex_id) %>% summarise(open = max(openingDay_wijk)) %>% dplyr::select(open))
# plot(siteEnv %>% filter(yr == 2016, yday==1) %>% dplyr::select(hex_type))


### Extract body mass ###
# from Dunning 2007. A transcript of the body masses for waterfowl (Anatidae family) is provided in the scripts folder
BMdf <- read.csv("inputData/Anatidae_bodymass_fromHandbookbyDunning.csv") %>%
  filter(grepl("Anser albifrons", Species.Name)) %>% 
  group_by(Common.Name) %>% summarise(mean = mean(Mean),
                                      min  = min(Mean),
                                      max  = max(Mean))

### User-defined parameters ###
parms <- list(
  
  ### Project/species name
  projectName = "Anser albifrons",
  
  ### Model run parameters
  startDay = 61,  ## Start date given as DOY (day of the year), see above why 61. This will be treated as "day 1" in the rest of the parametrisation.
  MaxRun   = 150,  ## Maximum time run in days, starting from startDay. This means the mode will run from March 1st till July.
  arrWindow = 14,  ## Arrival window in days, needed for calculating terminal reward
  NSites = length(unique(siteEnv$hex_id)), ## N. of all potential stopover (breeding sites are excluded one by one at a later stage)
  
  ### Species-specific Parameters
  order = "Anseriformes", # species order, used to calculate BMR. It can be NULL, in which case a generic equation calculating BMR for flying birds is used
  meanBodyMass_g = BMdf %>% pull(mean),
  minBodyMass_g  = BMdf %>% pull(min),
  maxBodyMass_g  = BMdf %>% pull(max),
  MaxX = 100,  ## Levels of body condition (N. steps dividing max body mass - lean min body mass )
  B0 = 3,       ## Future reproductive success (N. juveniles survived in the future)
  
  # Flight capacities:
  flightStrategy = NULL, # one of "landSoarer" or NULL (others). For estimation of flight cost based on Scacco et al. 2026
  speed_ms = 12.6, # Ground speed: either one numeric value in m/s or NULL. If provided, ground speed is assumed to be constant. If NULL, ground speed will be calculated based on wind support. For that, the arguments windTab_u/windTab_v must be provided within the makeSiteTab() function below.
  #speed_from_windSup = NULL, # Not yet implemented. Potentially an alternative equation to express ground speed as a function of wind support. If this is NULL and windSpeed_ms is also NULL, the default equation from Scacco et al. will be used to compute daily ground speeds from daily wind support values.
  
  # Decision Error
  decError = 4000 ## the higher (i.e. the closer 1/decErr to 0) the more suboptimal decisions are allowed so that similarly optimal routes are kept
  
)


### All combinations of start (wintering) and end (breeding) hexagons ###
StartEndTab <- siteEnv %>% filter(!is.na(hex_type), hex_type %in% c('wintering', 'breeding'), !duplicated(hex_id)) %>%
  dplyr::select(hex_id, hex_type) %>% st_drop_geometry() %>% 
  expand_grid(wintering = hex_id[hex_type=="wintering"], breeding = hex_id[hex_type!="wintering"])

  
### Formatting of remaining parameters for the model ### 
sdp_sites <- makeSiteTab(siteEnv    = siteEnv,
                         year       = 2016,
                         wintSites  = StartEndTab %>% slice(1) %>% pull(wintering), 
                         breedSites = StartEndTab %>% slice(1) %>% pull(breeding), 
                         listParms  = parms)

### Model ### 
bwd <- bwdIteration(sdp_sites, pbar = FALSE)

fwd <- fwdSimulation(bwd, NrInd = 100, start_t = 1, start_site = 1, start_x = c(25, 35))

save(bwd,fwd, file="output/bwd_fwd_w1130_b2644.rdata")

### Plotting ### 
condProfile(simu = fwd[50:51,,], model = sdp_sites)
pltNetwork(simu = fwd, model = sdp_sites, map = europe_laea)
