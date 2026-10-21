
#library(tidyverse)
library(dplyr)
library(sf)
library(geosphere)
#library(spdep)

# for plotting
library(rnaturalearth)
library(ggplot2)
library(patchwork)
sf_use_s2(FALSE)

### If not ran as package, import model functions as:
library(Rcpp)
library(RcppProgress)
library(RcppArmadillo)
Rcpp::sourceCpp("src/SdpMigC.cpp") # C++ functions
source("R/OptimSPD_fnc_v4_dynamicGain.R") # R functions

#_______________________
# Once all is running smoothly we delete the above and run once:
# Rcpp::compileAttributes() and then:
# Once the sdpMig2 code will be turned into a package the above can be deleted and it will be enough to load:
# library(remotes)
# remotes::install_github("mscacco/sdpMig_v2")
# library(sdpMig_v2)
#_______________________

#__________________
### Environment ###

# This is a spatial object with one observation per hexagon-doy combination + geometry
siteTab  <- readRDS("./inst/data/sample_siteTab_yday-GDDjerk.rds") |> st_as_sf()
nrow(siteTab) == length(unique(siteTab$hex_id)) * length(unique(siteTab$yday))

ggplot() + ggtitle("GDD jerk - day 61") +
  geom_sf(data = siteTab[siteTab$yday == 61,], aes(fill = GDDjerk_wijk)) +
  scale_fill_viridis_c(name = "GDD jerk", na.value = "transparent") + theme_minimal()


#______________________
### Start/End sites ###

# All combinations of start (wintering) and end (breeding) hexagons
allWintSites <- unique(siteTab$hex_id[which(siteTab$hex_type == "wintering")])
allBreedSites <- unique(siteTab$hex_id[which(siteTab$hex_type == "breeding")])
StartEndTab <- expand.grid(wintering = allWintSites, breeding = allBreedSites)

#________________
### Species Body mass ###

# from Dunning 2007. A transcript of the body masses for waterfowl (Anatidae family) is provided in inputData
BMdf <- read.csv("./inst/data/Anatidae_bodymass_fromHandbookbyDunning.csv") |>
  filter(grepl("Anser albifrons", Species.Name)) |>
  group_by(Common.Name) |> summarise(mean = mean(Mean),
                                      min  = min(Mean),
                                      max  = max(Mean))

# After making a package this should be available with
# data("....")


#_________________________________________
### User-defined parameters for model ####

listParms <- list(
  
  # Species name or name of specific parameters combination being tested
  projectName = "Anser albifrons",
  
  # Model run parameters
  startDOY = 61,   ## Start date given as DOY (day of the year). This will be treated as "day 1" in the rest of the parametrisation.
  endDOY   = 180,  ## End date as DOY. This means the mode will run for nDays = endDOY-startDOY.
  arrWindowBef = 10,  ## Arrival window in days before the opening day of the breeding site, needed for calculating terminal reward
  arrWindowAft = 10,  ## Arrival window in days after the opening day of the breeding site, needed for calculating terminal reward
  
  # Species-specific Parameters
  order = "Anseriformes", # species order, used to calculate BMR. It can be NULL, in which case a generic equation calculating BMR for flying birds is used
  meanBodyMass_g = BMdf$mean,
  minBodyMass_g  = BMdf$min,
  maxBodyMass_g  = BMdf$max,
  MaxX = 100,     ## Levels of body condition (N. steps dividing max body mass - lean min body mass )
  B0 = 3,         ## Future reproductive success (N. juveniles survived in the future)
  
  # Flight capacities:
  f = 0.5,  ## proportion of body reserve (100) that can be devoted to flying
  speed_ms = 12.6, # Ground speed: either one numeric value in m/s or NULL. If provided, ground speed is assumed to be constant. If NULL, ground speed will be calculated based on wind support. For that, the arguments windTab_u/windTab_v must be provided within the makeSiteTab() function below.
  #speed_from_windSup = NULL, # Not yet implemented. Potentially an alternative equation to express ground speed as a function of wind support. If this is NULL and windSpeed_ms is also NULL, the default equation from Scacco et al. will be used to compute daily ground speeds from daily wind support values.
  
  # Decision Error
  decError = 4000 ## the higher, the more optimal decisions are allowed so that similarly optimal routes are kept (the closer 1/decErr is to 0)
  
)


#_____________
### Model ###

# Test one of many combinations of wintering-breeding sites
# Note that computation time grows fast with increasing number of hexagons
i=1
wId <- StartEndTab[i, "wintering"]
bId <- StartEndTab[i, "breeding"]

ggplot() +
  geom_sf(data = filter(siteTab, yday==1), aes(fill = NULL)) +
  geom_sf(data = filter(siteTab, yday==1, hex_id==StartEndTab[i, "wintering"]), aes(fill = "Wintering site")) +
  geom_sf(data = filter(siteTab, yday==1, hex_id==StartEndTab[i, "breeding"]), aes(fill = "Breeding site")) +
  scale_fill_manual(name = "", values = c("Wintering site" = "darkblue","Breeding site" = "darkred")) +
  theme_minimal() + theme(legend.position = c(0.8, 0.05), legend.justification = c(0.5, 0))  # anchor from bottom-center

# Format site tab
siteTab_format <- siteTab %>% st_drop_geometry() %>% 
  dplyr::select(hex_id, yday, x, y, GDDjerk_wijk) %>%
  as.data.frame()

# Format parameters for model run
sdp_sites <- makeSDPMig(siteTab   = siteTab_format,
                        wintSite  = wId, 
                        breedSite = bId, 
                        listParms = listParms)
# Backward iteration
bwd <- bwdIteration(sdp_sites, pbar = FALSE)

# Forward simulation
simuObj <- fwdSimulation(bwd, NrInd = 1000, start_t = 1, start_site = 1, 
                         start_x = c(10, 20)) # low starting conditions, as geese start fattening just before starting migration


#_____________
### Output ###

# Quick simple plots for this specific combination
# Europe boundaries for plotting
europe <- ne_countries(scale = "medium", continent = "europe", returnclass = "sf")
europe_laea <- st_crop(st_transform(europe, st_crs(siteTab)), st_bbox(siteTab))

pltNetworkGG(simuObj, map = europe_laea)

condProfile(simuObj)

# Compute and plot simple output metrics for this specific combination
out <- computeOutputs(simuObj)

str(out)

hist(out$migPerInd$migrStartDOY)
hist(out$migPerInd$migrEndDOY)
hist(out$migPerInd$migrDurationDays)
hist(out$migPerInd$Nstops)
hist(out$migPerInd$meanStag)


