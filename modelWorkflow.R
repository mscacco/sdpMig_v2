library(tidyverse)
library(sf)
sf_use_s2(FALSE)
library(Rcpp)
library(RcppArmadillo)

# Load R and C++ model functions
Rcpp::sourceCpp("functions/src/SdpMigC.cpp")
source("functions/R/OptimSPD_fnc.R")

#___________
### Grid ###
hex <- st_read("inputData/hexagonalGrid_100km_laea_EuropeNoWater.gpkg")

europe      <- rnaturalearth::ne_countries(scale = "medium", continent = "europe", returnclass = "sf")
europe_laea <- st_transform(europe, st_crs(hex)) %>% st_intersection(st_bbox(hex) %>% st_as_sfc())

ggplot() + geom_sf(data = europe_laea) + geom_sf(data = hex, fill = NA)

#__________________
### Environment ###
siteEnv  <- readRDS("./inputData/EnvTab_forModel.rds") %>% st_as_sf()

ggplot() + ggtitle("Opening Day") + geom_sf(data = hex, aes(fill=NULL)) +
  geom_sf(data = filter(siteEnv, yday==1), aes(fill = openingDay_0)) + #opening day (one value per hex)
  scale_fill_viridis_c(name = "DOY", na.value = "transparent") + theme_minimal()

ggplot() + ggtitle("Geese ranges") + geom_sf(data = hex, aes(fill=NULL)) +
  geom_sf(data = filter(siteEnv, yday==1), aes(fill = hex_type)) +
  scale_fill_manual(values=c("darkorange","grey","blue"), name="Hex type", na.value = "transparent") +  theme_minimal()

#______________________
### Start/End sites ###
# All combinations of start (wintering) and end (breeding) hexagons
allWintSites <- unique(siteEnv$hex_id[siteEnv$hex_type == "wintering"])
allBreedSites <- unique(siteEnv$hex_id[siteEnv$hex_type == "breeding"])
StartEndTab <- expand.grid(wintering = allWintSites, breeding = allBreedSites)

#________________
### Sites Tab ###
# Input tab for model below, with one row per site and with columns in the correct order
siteTab <- siteEnv %>% filter(!duplicated(hex_id)) %>%
  select(hex_id, openingDay_0, x, y) %>%
  st_drop_geometry()

#________________
### Body mass ###
# from Dunning 2007. A transcript of the body masses for waterfowl (Anatidae family) is provided in inputData
BMdf <- read.csv("inputData/Anatidae_bodymass_fromHandbookbyDunning.csv") %>%
  filter(grepl("Anser albifrons", Species.Name)) %>% 
  group_by(Common.Name) %>% summarise(mean = mean(Mean),
                                      min  = min(Mean),
                                      max  = max(Mean))

#______________________________
### User-defined parameters ###
parms <- list(
  
  # Project/species name
  projectName = "Anser albifrons",
  
  # Model run parameters
  startDOY = 61,  ## Start date given as DOY (day of the year), see above why 61. This will be treated as "day 1" in the rest of the parametrisation.
  endDOY   = 180,  ## End date as DOY. This means the mode will run for nDays = endDOY-startDOY.
  arrWindow = 14,  ## Arrival window in days, needed for calculating terminal reward
  NSites = length(unique(siteEnv$hex_id)), ## N. of all potential stopover (breeding sites are excluded one by one at a later stage)
  
  # Species-specific Parameters
  order = "Anseriformes", # species order, used to calculate BMR. It can be NULL, in which case a generic equation calculating BMR for flying birds is used
  meanBodyMass_g = BMdf$mean,
  minBodyMass_g  = BMdf$min,
  maxBodyMass_g  = BMdf$max,
  MaxX = 100,  ## Levels of body condition (N. steps dividing max body mass - lean min body mass )
  B0 = 3,       ## Future reproductive success (N. juveniles survived in the future)
  
  # Flight capacities:
  flightStrategy = NULL, # one of "landSoarer" or NULL (others). For estimation of flight cost based on Scacco et al. 2026
  speed_ms = 12.6, # Ground speed: either one numeric value in m/s or NULL. If provided, ground speed is assumed to be constant. If NULL, ground speed will be calculated based on wind support. For that, the arguments windTab_u/windTab_v must be provided within the makeSiteTab() function below.
  #speed_from_windSup = NULL, # Not yet implemented. Potentially an alternative equation to express ground speed as a function of wind support. If this is NULL and windSpeed_ms is also NULL, the default equation from Scacco et al. will be used to compute daily ground speeds from daily wind support values.
  
  # Decision Error
  decError = 4000 ## the higher (i.e. the closer 1/decErr to 0) the more suboptimal decisions are allowed so that similarly optimal routes are kept
  
)


#_____________
### Model ###

# (one model run per combination of wintering-breeding sites)
i=1
#lapply(1:nrow(StartEndTab), function(i){

  sdp_sites <- makeSiteTab(siteTab    = siteTab,
                           wintSite  = StartEndTab[i, "wintering"], 
                           breedSite = StartEndTab[i, "breeding"], 
                           listParms  = parms)
  bwd <- bwdIteration(sdp_sites, pbar = FALSE)
  
  fwd <- fwdSimulation(bwd, NrInd = 100, start_t = 1, start_site = 1, start_x = c(25, 35))
  
  save(bwd,fwd, file="output/bwd_fwd_w1130_b2644.rdata")
#})

# Plot simulation output
condProfile(simu = fwd, model = sdp_sites)
pltNetwork(simu = fwd, model = sdp_sites, map = europe_laea)

