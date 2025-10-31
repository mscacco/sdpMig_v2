library(tidyverse)
library(sf)
sf_use_s2(FALSE)
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("simeon/functions/src/SdpMigC.cpp")
source("./simeon/functions/R/fnc_simeon.R")

### Grid
hex <- st_read("simeon/inputData/hexagonalGrid_100km_laea_EuropeNoWater.gpkg")

europe      <- rnaturalearth::ne_countries(scale = "medium", continent = "europe", returnclass = "sf")
europe_laea <- st_transform(europe, st_crs(hex)) %>% st_intersection(st_bbox(hex) %>% st_as_sfc())

ggplot() +
  geom_sf(data = europe_laea) +
  geom_sf(data = hex, fill = NA)


### Environment
siteEnv  <- readRDS("simeon/inputData/allENV_AllYears.rds") %>%
  left_join(readRDS("simeon/inputData/HexagonTypes_AnserAlbifrons.rds") %>% dplyr::select(hex_id, hex_type) , by = "hex_id") %>%
  filter(hex_id %in% hex$hex_id) %>% st_as_sf()

# plot(siteEnv %>% filter(yr == 2016) %>% group_by(yr, hex_id) %>% summarise(open = max(openingDay_wijk)) %>% dplyr::select(open))
# plot(siteEnv %>% filter(yr == 2016, yday==1) %>% dplyr::select(hex_type))


### Species parameters
BMdf <- read.csv("simeon/inputData/Anatidae_bodymass_fromHandbookbyDunning.csv") %>%
  filter(grepl("Anser albifrons", Species.Name)) %>% 
  group_by(Common.Name) %>% summarise(mean = mean(Mean),
                                      min  = min(Mean),
                                      max  = max(Mean))


### Parameters ###

parms <- list(
  
  ### Project/species name
  projectName = "Anser albifrons",
  
  ### Model run parameters ###
  startDay = 61,  
  MaxRun   = 180,  
  arrWindow = 14, 
  NSites = length(unique(siteEnv$hex_id)),
  
  ### Species-specific Parameters ###
  order = "Anseriformes",
  meanBodyMass_g = BMdf %>% pull(mean),
  minBodyMass_g  = BMdf %>% pull(min),
  maxBodyMass_g  = BMdf %>% pull(max),
  MaxX = 100,  
  B0 = 3,
  
  # Flight capacities:
  flightStrategy = NULL, 
  speed_ms = 12.6, 
  
  # Decision Error
  decError = 4000
  
)


### start /end
StartEndTab <- siteEnv %>% filter(!is.na(hex_type), hex_type %in% c('wintering', 'breeding'), !duplicated(hex_id)) %>%
  dplyr::select(hex_id, hex_type) %>% st_drop_geometry() %>% 
  expand_grid(wintering = hex_id[hex_type=="wintering"], breeding = hex_id[hex_type!="wintering"])

  
sdp_sites <- makeSiteTab(siteEnv    = siteEnv,
                         year       = 2016,
                         wintSites  = StartEndTab %>% slice(1) %>% pull(wintering), 
                         breedSites = StartEndTab %>% slice(1) %>% pull(breeding), 
                         listParms  = parms)

bckw <- bwdIteration(sdp_sites, FALSE)


fwd <- fwdSimulation(bckw, NrInd = 100, start_t = 1, start_site = 1, start_x = c(10, 20))

condProfile(simu = fwd, model = sdp_sites)
pltNetwork(simu = fwd, model = sdp_sites, map = europe_laea)
