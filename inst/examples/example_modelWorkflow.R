
library(ggplot2)
library(sf)
sf_use_s2(FALSE)

### If not ran as package, import model functions as:
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("src/SdpMigC.cpp") # C++ functions
source("R/OptimSPD_fnc_v2.R") # R functions

# Once all is running smoothly we delete the above and run once:
# Rcpp::compileAttributes() and then:
# library(sdpMig_v2)

#__________________
### Environment ###

hex <- st_read("./inst/data/hexagonalGrid_100km_laea_EuropeNoWater.gpkg")
siteTab  <- readRDS("./inst/data/sample_siteTab.rds") |> st_as_sf()

ggplot() + ggtitle("Opening Day") + geom_sf(data = hex, aes(fill=NULL)) +
  geom_sf(data = siteTab, aes(fill = openingDay)) +
  scale_fill_viridis_c(name = "DOY", na.value = "transparent") + theme_minimal()

#______________________
### Start/End sites ###
# All combinations of start (wintering) and end (breeding) hexagons
allWintSites <- unique(siteTab$hex_id[siteTab$hex_type == "wintering"])
allBreedSites <- unique(siteTab$hex_id[siteTab$hex_type == "breeding"])
StartEndTab <- expand.grid(wintering = allWintSites, breeding = allBreedSites)


#_______________________________________
### Exclude low suitability hexagons ###

# This is just an example filtering that is not biologically supported in this case 
# Keep hexagons with suitability >= 50% but still including wintering and breeding grounds
siteTab <- filter(siteTab, suitability>=0.5 | hex_id %in% c(allWintSites,allBreedSites))

ggplot() + ggtitle("Suitability > 50%") + 
  geom_sf(data = hex, aes(fill=NULL)) +
  geom_sf(data = filter(st_as_sf(siteTab)), 
          aes(fill = suitability)) +
  scale_fill_viridis_c(name = "% Suitability", na.value = "transparent") + theme_minimal()


#________________
### Body mass ###
# from Dunning 2007. A transcript of the body masses for waterfowl (Anatidae family) is provided in inputData
BMdf <- read.csv("./inst/data/Anatidae_bodymass_fromHandbookbyDunning.csv") %>%
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
  startDOY = 61,  ## Start date given as DOY (day of the year). This will be treated as "day 1" in the rest of the parametrisation.
  endDOY   = 180,  ## End date as DOY. This means the mode will run for nDays = endDOY-startDOY.
  arrWindow = 14,  ## Arrival window in days, needed for calculating terminal reward
  NSites = nrow(siteTab), ## N. of all potential stopover (breeding sites are excluded one by one at a later stage)
  
  # Species-specific Parameters
  order = "Anseriformes", # species order, used to calculate BMR. It can be NULL, in which case a generic equation calculating BMR for flying birds is used
  meanBodyMass_g = BMdf$mean,
  minBodyMass_g  = BMdf$min,
  maxBodyMass_g  = BMdf$max,
  MaxX = 100,  ## Levels of body condition (N. steps dividing max body mass - lean min body mass )
  B0 = 3,       ## Future reproductive success (N. juveniles survived in the future)
  
  # Flight capacities:
  f = 0.7,  ## proportion of body reserve (100) that can be devoted to flying
  flightStrategy = NULL, # one of "landSoarer" or NULL (others). For estimation of flight cost based on Scacco et al. 2026
  speed_ms = 12.6, # Ground speed: either one numeric value in m/s or NULL. If provided, ground speed is assumed to be constant. If NULL, ground speed will be calculated based on wind support. For that, the arguments windTab_u/windTab_v must be provided within the makeSiteTab() function below.
  #speed_from_windSup = NULL, # Not yet implemented. Potentially an alternative equation to express ground speed as a function of wind support. If this is NULL and windSpeed_ms is also NULL, the default equation from Scacco et al. will be used to compute daily ground speeds from daily wind support values.
  
  # Decision Error
  decError = 4000 ## the higher (i.e. the closer 1/decErr to 0) the more suboptimal decisions are allowed so that similarly optimal routes are kept
  
)


#_____________
### Model ###

# Test one of many combinations of wintering-breeding sites
# Note that computation time grows fast with increasing number of hexagons
i=1
wId <- StartEndTab[i, "wintering"]
bId <- StartEndTab[i, "breeding"]

ggplot() + geom_sf(data = hex, aes(fill=NULL)) +
  geom_sf(data = filter(siteTab, hex_id==StartEndTab[i, "wintering"]), aes(fill = "Wintering site")) +
  geom_sf(data = filter(siteTab, hex_id==StartEndTab[i, "breeding"]), aes(fill = "Breeding site")) +
  scale_fill_manual(name = "", values = c("Wintering site" = "darkblue","Breeding site" = "darkred")) +
  theme_minimal() + theme(legend.position = c(0.8, 0.05), legend.justification = c(0.5, 0))  # anchor from bottom-center

# Format parameters for model run
sdp_sites <- makeSiteTab(siteTab   = st_drop_geometry(siteTab),
                         wintSite  = wId, 
                         breedSite = bId, 
                         listParms = parms)
# Backward iteration
bwd <- bwdIteration(sdp_sites, pbar = FALSE)

# Forward simulation
fwd <- fwdSimulation(bwd, NrInd = 100, start_t = 1, start_site = 1, start_x = c(25, 35))


# Plot simulation output
condProfile(simu = fwd, model = sdp_sites)
pltNetwork(simu = fwd, model = sdp_sites, map = europe_laea)

