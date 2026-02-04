
library(data.table)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)

#source("./Scripts/WaterfowlMigration/PLOT_functions.R")

### Grid and borders for maps ###
hex <- st_read("inputData/hexagonalGrid_100km_laea_EuropeNoWater.gpkg")
europe <- ne_countries(scale = "medium", continent = "europe", returnclass = "sf")
europe_laea <- st_transform(europe, st_crs(hex))

### Import model output ###
outFld <- "/home/mscacco/ownCloud - mscacco@ab.mpg.de@owncloud.gwdg.de/Martina/ProgettiVari/WSL/KappaFlu/Outputs/Simulations/"
load(paste0(outFld,"bwd_fwd_w1130_b2644_Nov4.rdata")) #bwd, fwd

### Count number of days spent in each site ###
# in dim 2 of fwd we have: day, site, x, fly/not, death. 
# Fly/not takes either the ID of the next site, or 0 for alive individuals not using the site or NA for dead individuals.
# Death takes 0 for alive and 1 for dead.
# daysPerSite has dim1=NrInds, dim2=NrSites, dim3=nDays. 
# sites are ordered according to bwd@Sites$SitesID
daysPerSite <- matrix(NA, nrow = dim(fwd)[1], ncol = bwd@Init$NSites + 1)
for (i in 1:nrow(daysPerSite)) {
  tmp <- fwd[i, , ] # for each individual i, extract dim2 and dim3
  if (sum(tmp[5, ]) > 0) # extract column death
    tmp <- tmp[, 1:(min(which(tmp[5, ] == 1)))] # exclude days after death
  daysPerSite[i, ] <- apply(cbind(0:bwd@Init$NSites), 1, function(x) { #for each site
    ifelse(all(tmp[2, ] < x), NA, sum(tmp[2, ] == x, # count how many days it spent in each site
                                      na.rm = T)) 
  })
}


### Plot migratory route per individual ###
dim(daysPerSite)
daysPerSite[,1:10] #staging in first 10 sites for example
# ID of wintering site is first, breeding site is last
wintID <- bwd@Sites$SitesID[1]
breedID <- bwd@Sites$SitesID[bwd@Init$NSites]

# Store individual migratory routes
# each row of this matrix represents an individual route, indicating where the individual stopped and for how long
migRoutes <- rbindlist(lapply(1:nrow(daysPerSite), function(i){ # for each individual
  stopSites <- which(daysPerSite[i,]>0) # index of stopover sites
  if(length(stopSites)>0){ # exclude dead individuals
    staging <- daysPerSite[i,stopSites]
    #stopSitesID <- rownames(bwd@Sites$crds)[stopSites] # hexID corresponding to the stopover sites
    stopSitesID <- bwd@Sites$SitesID[stopSites] # hexID corresponding to the stopover sites, including wintering and breeding
    indRoute <- data.frame(ind = i,
                           stopSitesID=as.character(stopSitesID), # corresponding to hexID
                           stagDays=staging)
    return(indRoute)}
}))
# Before plotting, extract coordinates of hexagons centroids in laea
hex$hex_id <- as.character(hex$hex_id)
hexCentroids_laea <- st_centroid(hex)
stopoverHex <- merge(migRoutes, hex, by.x="stopSitesID", by.y="hex_id", all.x=T)
migRoutes_pts <- st_as_sf(merge(migRoutes, hexCentroids_laea, by.x="stopSitesID", by.y="hex_id", all.x=T))
migRoutes_lines <- migRoutes_pts %>% 
  arrange(ind, stopSitesID) %>%  # make sure points are in order
  group_by(ind) %>% summarise(do_union = FALSE) %>%      # keep all points per ind
  st_cast("LINESTRING") # transform the points in lines

# Plot migratory routes per individual + stopover duration (staging time)
ggplot() +
  geom_sf(data = europe_laea, fill = NA, col = "black", size = 1.2) +
  geom_sf(data = hex, fill = NA, col = "darkgrey") +
  #geom_sf(data = migRoutes_lines, aes(color = factor(ind)), linewidth = 1.5, alpha = 0.1) + # one color per ind
  geom_sf(data = migRoutes_lines, col = "black", linewidth = 1.5, alpha = 0.2) +
  geom_sf(data = migRoutes_pts, aes(size = stagDays), shape = 21, col = "black", fill = "darkorange") +
  scale_size(range = c(1, 10)) + # scale all point sizes to show bigger
  geom_sf(data = hex[hex$hex_id==wintID,], fill = "blue") +
  geom_sf(data = hex[hex$hex_id==breedID,], fill = "red") +
  coord_sf(
    xlim = st_bbox(hex)[c("xmin", "xmax")],
    ylim = st_bbox(hex)[c("ymin", "ymax")],
    expand = FALSE
  ) +
  theme_minimal()  

### Important output measures ###
# individuals per site, mortality, staging times per site, body condition per site, fuelling rates per site 
# in daysPerSite, sites are ordered according to bwd@Sites$SitesID

# 1. N. of individuals per site used
ind_site <- apply(daysPerSite, 2, function(x){sum(x!=0, na.rm=T)})
names(ind_site) <- bwd@Sites$SitesID
(ind_site[ind_site != 0])

# 2. N. of stopovers per individual excluding start/end
stop_ind <- apply(daysPerSite, 1, function(x){sum(x!=0, na.rm=T)})
stop_ind <- stop_ind - 2 # excluding start and arrival sites
hist(stop_ind, main = "N. of stopovers per individual")

# 3. Mean staging times per site across individuals 
#daysPerSite[apply(daysPerSite, 1, function(x) any(is.na(x))),] <- cbind(rep(NA, ncol(daysPerSite)))
mean.staging <- apply(daysPerSite, 2, mean, na.rm = T)
summary(mean.staging[mean.staging!=0])

# 4. Timing of migration
dim(fwd)

# 5. Mortality proportion in population
nDeaths <- sum(apply(fwd[,5,], 1, function(x) any(x==1)))
nDeaths/dim(fwd)[1] #deaths/nInds, proportion of deaths in the population

# 4. Mean body condition per site over time across individuals ????
meanX <- matrix(ncol = bwd@Init$MaxT, nrow = bwd@Init$NSites+1) # nSites per nDays
meanX <- meanX[,1:dim(fwd)[3]]  # this is to be removed after the MaxT problem is solved
for(i in 1:ncol(meanX)) { # for each day, average the body condition of all individuals per site
  t1 <- aggregate(fwd[,3,i], by = list(Site=fwd[,2,i]), FUN = function(x) mean(x, na.rm = T))
  meanX[t1$Site+1, i]  <- t1$x # and add x to corresponding day and site in the matrix
}
t(meanX)
# Mean body condition per site overall
(mean.bodycond <- apply(t(meanX), 2, weighted.mean, w= , na.rm=T))
(mean.bodycond <- apply(t(meanX), 2, mean, na.rm=T))


