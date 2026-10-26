# Example workflow for the v5 R functions (R/OptimSPD_fnc_v5_dynamicGain.R). listParms$order is not needed in v5 (no BMR calculation).
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

### Load the package (installed with remotes::install_github("mscacco/sdpMig_v2"), or from the repo folder with devtools::load_all())
library(sdpMig2)

# Without installing, the functions can also be loaded from the repo folder:
# Rcpp::sourceCpp("src/SdpMigC.cpp") # C++ functions. The original C++ is in inst/archive/SdpMigC_original.cpp
# source("R/OptimSPD_fnc_v5_dynamicGain.R") # R functions (v5)

#__________________
### Environment ###

# This is a spatial object with one observation per hexagon-doy combination + geometry
siteTab  <- readRDS(system.file("extdata", "sample_siteTab_yday-GDDjerk.rds", package = "sdpMig2")) |> st_as_sf()
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
BMdf <- read.csv(system.file("extdata", "Anatidae_bodymass_fromHandbookbyDunning.csv", package = "sdpMig2")) |>
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
  meanBodyMass_g = BMdf$mean,
  minBodyMass_g  = BMdf$min,
  maxBodyMass_g  = BMdf$max,
  MaxX = 100,     ## Levels of body condition (N. steps dividing max body mass - lean min body mass )
  B0 = 3,         ## Future reproductive success (N. juveniles survived in the future)
  
  # Flight capacities:
  f = 1,  ## proportion of MaxX usable as flight fuel; (1 - f) * MaxX is a reserve that is never used for flight (a bird with x <= (1 - f) * MaxX cannot fly)
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
hist(out$migPerInd$migrDurationDays)
# Number of stopovers (wintering and breeding sites excluded) per outcome:
# arrived alive within the arrival window, arrived alive outside the window, unsuccessful (died or did not reach the breeding site).
# For unsuccessful individuals Nstops counts the stops made until death or the end of the run.
nsGroups <- list("Arrived within window"  = out$migPerInd$inWindow,
                 "Arrived outside window" = out$migPerInd$success & !out$migPerInd$inWindow,
                 "Unsuccessful"           = !out$migPerInd$success)
nsBreaks <- seq(min(out$migPerInd$Nstops) - 0.5, max(out$migPerInd$Nstops) + 0.5, 1)
for (g in names(nsGroups)) {
  hist(out$migPerInd$Nstops[nsGroups[[g]]], breaks = nsBreaks,
       main = paste0(g, " (n = ", sum(nsGroups[[g]]), ")"), xlab = "Number of stopovers")
}
hist(out$migPerInd$meanStag)

# Arrival at the breeding site (migrEndDOY is NA for individuals that died or did not reach the breeding site)
out$migStats # SuccessProp (breeding site reached alive at any time) and SuccessWindowProp (reached within the arrival window)
# nRoutes and topRouteShare: number of distinct routes and share of the most common route (individuals reaching the breeding site alive)
arrWindowDOY <- simuObj@Sites$breedOpeningDOY + c(-listParms$arrWindowBef, listParms$arrWindowAft) # arrival window in DOY
inW  <- out$migPerInd$inWindow
outW <- !out$migPerInd$inWindow & !is.na(out$migPerInd$migrEndDOY)
xlim <- range(c(out$migPerInd$migrEndDOY, arrWindowDOY), na.rm = TRUE)
# Arrival DOY of individuals arriving within the arrival window
hist(out$migPerInd$migrEndDOY[inW], breaks = seq(xlim[1]-0.5, xlim[2]+0.5, 1), xlim = xlim,
     main = paste0("Arrival within window (n = ", sum(inW), ")"), xlab = "Arrival DOY at breeding site")
abline(v = arrWindowDOY, lty = 2) # window limits
# Arrival DOY of individuals arriving at the breeding site outside the arrival window
hist(out$migPerInd$migrEndDOY[outW], breaks = seq(xlim[1]-0.5, xlim[2]+0.5, 1), xlim = xlim,
     main = paste0("Arrival outside window (n = ", sum(outW), ")"), xlab = "Arrival DOY at breeding site")
abline(v = arrWindowDOY, lty = 2)
# Arrival relative to the opening day of the breeding site (all individuals that reached the breeding site)
hist(out$migPerInd$arrVsOpeningDays)
# Terminal reward at arrival (arrival time and body condition combined), individuals that reached the breeding site alive.
# Arrivals outside the arrival window all have arrivalTR = B0 (dashed line), so they appear as a peak at B0.
ggplot(out$migPerInd[!is.na(out$migPerInd$arrivalTR), ], aes(x = arrivalTR)) +
  geom_density(fill = "orange", colour = "orange", alpha = 0.35, linewidth = 0.5) +
  geom_vline(xintercept = simuObj@Species$B0, linetype = "dashed") +
  labs(x = "Terminal reward at arrival (dashed: B0, arrival outside window)", y = "Density",
       title = paste0("n = ", sum(!is.na(out$migPerInd$arrivalTR)), " individuals reaching the breeding site alive")) +
  theme_minimal()
