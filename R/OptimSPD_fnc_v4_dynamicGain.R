#!! Major change relative to previous version is the implementation of a dynamic suitability layer that changes with DOY
# This replaced the previous suitability layer and changed the interaction with gain

setClass(
  "SDPMig",
  slots = c(Name    = 'character',
            Init    = "list",
            Species = "list",
            Sites   = "list",
            Results = "list")
)

Mat3Array      <- function(x) {
  
  rows = length(x)
  cols = length(x[[1]])
  d3 = length(x[[1]][[1]])
  
  out <- array(dim = c(rows, cols, d3))
  
  for(r in 1:rows) {
    for(c in 1:cols) {
      for(d.1 in 1:d3) {
        out[r, c, d.1] <- x[[r]][[c]][[d.1]]
      }
    }
  }
  out 
}

CalcTR         <- function(x, t) {
  TR <- approx(parms$xFTReward, parms$yFTReward, t, rule = 3)$y
  s  <- parms$w * (x - parms$xc)
  
  if(x == 0) { SR <- 0 } else {
    SR <- (((exp(s) - exp(-s)) / (exp(s) + exp(-s))) + 1)/2
  }
  TR * SR + parms$B0
}


### Sampling individuals per hexagon based on probabilities introducing real stochasticity
# probs is a vector of relative probabilities per hexagon/pixel (their sum should equal 1)
# returns a vector of the same length as prob, with number of individuals
find_nInds <- function(probs){
  as.integer(table(
    factor(sample(seq_len(length(probs)), size = 1000, replace = TRUE, prob = probs),
           levels = seq_len(length(probs)))
  ))}


#############################
### Create Parameters Tab ###
#############################

# Inputs are: 
# 1. a data.frame or matrix with nrow being nSites * 366 (year days) with 5 columns in this order: 
# - 1st col = sites ID (hexagons id), 
# - 2nd col = day of the year (Jan 1st being day 1), 
# - 3rd and 4th columns are coordinates (long and lat in decimal degrees),
# - 5th column is GDD jerk in its original units
#
# 2. two integers indicating the site ID (hexagon ID) of the wintering and breeding site. These should be among the IDs in the first column of data.
# 3. the list of parameters from previous steps,

makeSDPMig <- function(siteTab,
                       wintSite, 
                       breedSite,
                       listParms,
                       BMR_file = "./inst/data/BMRcalculation_perOrder_McNab2008.csv",
                       pred_a1  = 2, pred_a2  = 2,
                       b0 = 1e-04, b1 = 0.001, b2 = 1e-04,
                       stoch = 1
){
  require(geosphere)
  require(terra)
  require(sf)
  require(data.table)
  
  #___Species parameters to leave constant________
  w         <- 0.028  ## Parameters for sigmoidal TR function: slope
  xc        <- 55     ## Parameters for sigmoidal TR function: body condition (x) at inflection point of sigmoid
  max_u     <- 1.0    ## maximum foraging intensity (max prop. of day spent foraging)
  MaxReward <- 2      ## value within arrival window, for terminal reward calculation 

  #___Sort hex-day and subset Site Tab for Model Days___
  # !! matrix fills below (deiM) assume rows ordered by site then day, with all days present for every site
  siteTab <- siteTab[order(siteTab[,1], siteTab[,2]), ]
  modelRun <- length(listParms$startDOY:listParms$endDOY)
  siteTabDay <- siteTab[siteTab[,2] >= listParms$startDOY & # one row per site-day combination
                          siteTab[,2] <= listParms$endDOY,]
  stopifnot(all(table(siteTabDay[,1]) == modelRun)) # every site has every model day
  
  #___Sorting sites (breeding site last, wintering site first)___
  hexIDs        <- unique(siteTab[,1])
  hex_idx       <- 1:length(hexIDs) 
  hex_idx_new   <- c(which(hexIDs==wintSite), hex_idx[hex_idx!=which(hexIDs==wintSite) & hex_idx!=which(hexIDs==breedSite)], which(hexIDs==breedSite))
  hexIDs_sort   <- hexIDs[hex_idx_new]
  siteTab_sort  <- siteTab[!duplicated(siteTab[,1]),][hex_idx_new,] # one row per site
  
  #___Basal Metabolic Rate__________ 
  {
    # Calculate Basal Metabolic Rate based on body mass (1 Watt = 1 Joule/second = 3.6 kJ/hour)
    # Calculation order-specific, based on McNab 2008 "Ecological factors affect the level and scaling of avian BMR"
    # File provided together with the code in the scripts folder
    BMRtab <- read.csv(BMR_file)
    
    if(listParms$order %in% BMRtab$Order){
      idx <- which(BMRtab$Order == listParms$order)
      BMR_kJhour <- BMRtab$coefficient[idx] * (listParms$meanBodyMass_g ^ BMRtab$power[idx])
    } else if(!listParms$order %in% BMRtab$Order | is.null(listParms$order)){
      warning("Order not found, using generic BMR calculation for flying birds: BMR kJ/h = 0.105 * (mass ^ 0.708)")
      BMR_kJhour <- 0.105 * (listParms$meanBodyMass_g ^ 0.708) # general calculation if order is missing, parallel to reg line of mammals when excluding non flying birds (also McNab 2008)
    }
  }
  
  BMR_kJday <- BMR_kJhour * 24
  
  #___Energy per body condition level_______ 
  {
    # 1 g fat = 39.6 kJ from Jenni & Jenni-Eiermann 1998 "Fuel Supply and Metabolic Constraints in Migrating Birds"
    e <- 39.6 # e is the energetic value per gram of fat in kJ/g from the above source
    maxFat_g <- (listParms$maxBodyMass_g - listParms$minBodyMass_g)
    maxEnergy_kJ <- maxFat_g * e # this is the maximum energy available when body condition are max
    energyPerX <- maxEnergy_kJ/listParms$MaxX # energetic equivalent (in kJ) per body condition unit (X level)
    
    #___Daily Energy Expenditure DEE / Field Metabolic Rate FMR per site__________
    # In Lisovski et al. PNAS, DEE is assumed to be 1.5 * BMR
    # But according to Nagy 2005: "there is no trend in present data suggesting that FMR slopes are identical to BMR slopes." (therefore that FMR is a constant multiple of BMR)
    # FMR in kJ/day = 10.5 * (mass in g ^0.681), For birds based on 95 species, from Nagy 2005 "Field metabolic rate and body size"
    FMR_kjday <- 10.5 * (listParms$meanBodyMass_g^0.681) #This gives us 2181.735 kJ/day
    
    # This must be converted in units of X (body condition). How many X does FMR corresponds to? 
    # We obtain it by dividing FMR for the energy equivalent of a unit of X:
    dee_inX <- FMR_kjday / energyPerX # equivalent to 5.96 units of X
    
    # Alternatively we could use 1.5 * BMR as approximation of DEE
    #dee_inX <- (BMR_kJday * 1.5) / energyPerX # equivalent to only 3.28 units of X, repeated for all sites and all days. 2*BMR is 4.38 units of X and 2.5*BMR 5.47 (almost as FMR)
    
    deeMat <- matrix(dee_inX, nrow = length(hexIDs), ncol = modelRun)
  }
  
  #___ Gain per site __________
  
  # Calculate the maximum Daily Metabolised Energy Intake (DMEmax) according to Kirkwood 1983 “A limit to metabolisable energy intake in mammals and birds”.
  DMEmax <- 1713 * ((listParms$minBodyMass_g/1000) ^ 0.72)
  DMEmax_inX <- DMEmax / energyPerX # gain = DME in units of body condition X, equal to 7.9

  # DMEmax already excludes the BMR (or FMR), so in order to obtain the maximum gain, BMR (or FMR) have to be summed to DMEmax
  gainMax <- DMEmax_inX + dee_inX

  ### Gain per day is 0 before opening day, then it is modulated by suitability between 0 and DMEmax_inX
  gain.x <- seq(1:modelRun)
  gain.p  <- rep(stoch, modelRun) # stoch is the same all across
  gain.y <- matrix(gainMax, nrow = length(hexIDs), ncol = modelRun) # we could use 0 or a minimum dee_inX
  
  #### !!!! We changed gain to being dynamic and changing daily with daily suitability values
  deiM <- matrix(siteTabDay[,5], nrow = nrow(gain.y), 
                 ncol = modelRun, byrow = T)[hex_idx_new, ] # sort
  
  deiMnew <- t(apply(deiM, 1, function(x) {
    mx <- which(x==max(x, na.rm = T))
    x[mx:length(x)] <- max(x, na.rm = T)
    x
  }))
  deiM_scale <- scales::rescale(deiMnew, c(dee_inX/gainMax, 1), range(deiMnew))
  
  gain.y <- gain.y * deiM_scale
  # matplot(t(deiMnew), type = "l", lty = 1, lwd = 1) # modified gdd jerk time series
  # matplot(t(deiM_scale), type = "l", lty = 1, lwd = 1) # modified gdd jerk time series
  # matplot(t(gain.y), type = "l", lty = 1, lwd = 1) # gain time series
  # hist(gain.y[,35]) # gain at migration peak (around doy 80)
  # hist(gain.y[,modelRun]) # gain at last model day
  # ggplot() + geom_sf(data=st_as_sf(siteTab_sort, coords = c(3,4), crs = 4326), aes(col = gain.y[,1]))
  
  # Fix gain of wintering site always to their DEE + the mean gain of all sites, assuming the bird knows the wintering ground and can always find resources
  gain.y[which(hexIDs_sort==wintSite),] <- gainMax #dee_inX * 1.1
  
  # Accuracy (foraging success, sampled from a normal distribution)
  ZStdNorm = c(-2.5, -2.0, -1.5, -1.0, -0.5,  0.0,  0.5,  1.0,  1.5,  2.0,  2.5)
  PStdNorm = c(0.0092, 0.0279, 0.0655, 0.1210, 0.1747, 0.2034, 0.1747, 0.1210, 0.0655, 0.0279, 0.0092)
  
  #___Exclude unsuitable sites______
  
  # To speed up computation, we exclude sites that are consistently unsuitable (gain < FMR) throughout the year.
  # But we make sure that wint and breed sites (first and last) are kept
  
  unsuitSites <- apply(gain.y, 1, function(x) all(x <= dee_inX))
  unsuitSites[c(1, length(unsuitSites))] <- FALSE
  # ggplot() + geom_sf(data=st_as_sf(siteTab_sort, coords = c(3,4), crs = 4326), aes(col = unsuitSites))
  
  # Remove the sites consistently from all matrices and vectors
  deeMat <- deeMat[! unsuitSites,]
  gain.y <- gain.y[! unsuitSites,]
  hexIDs_sort <- hexIDs_sort[! unsuitSites]
  siteTab_sort <- siteTab_sort[! unsuitSites,]
  
  #___ Terminal Reward__________
  # Calculate terminal reward for breeding sites (nrow = N sites, ncol defining the "curve shape" of the terminal reward):
  # x = Model day at which each site opens, y = reward. Breeding site is the last in the vector!
  
  jerkBreed <- gain.y[nrow(gain.y),]
  openingDayBreed <- min(which(jerkBreed == max(jerkBreed)))
  
  xFTReward <- c(1, 
                 openingDayBreed - listParms$arrWindowBef -1, openingDayBreed - listParms$arrWindowBef,
                 openingDayBreed + listParms$arrWindowAft, openingDayBreed + listParms$arrWindowAft + 1, 
                 length(gain.x))
  yFTReward <- c(0, 0, MaxReward, MaxReward, 0, 0)
  # plot(xFTReward, yFTReward, type="l")
  
  # ___ Predation Rates per site______
  b0 <- rep(b0, nrow(siteTab_sort)) # Background predation risk (e.g. number of predators in the area)
  b1 <- rep(b1, nrow(siteTab_sort)) # Body reserve dependent risk (e.g. the fatter the less mobile to escape)
  b2 <- rep(b2, nrow(siteTab_sort)) # Foraging intensity dependent risk (e.g. the more they feed the less they see predators)

  # !!! Increase predation risk outside wintering grounds (1st site) to "simulate" the risk of migration relative to staying at wintering ground
  b0[1] <- b0[1] * 0.001
  b1[1] <- b1[1] * 0.001
  b2[1] <- b2[1] * 0.001
  
  #___ Distances between all sites__________
  # Matrix of distances between all sites, coords are assumed to be in decimal degrees
  distMat <- as.matrix(terra::distance(x = siteTab_sort[,3:4], lonlat = TRUE, method = "geo", unit="km"))
  # plot distances from wintering ground (first site)
  # ggplot(siteTab_sort, aes(x = x, y = y, color = distMat[,1])) +
  #   geom_point(shape=16, size=3) +
  #   scale_color_viridis_b() + coord_equal()
  
  
  #___Flight parameters_______
  # 1. Flight cost
  cot <- exp(1.587) * ((listParms$meanBodyMass_g/1000) ^ -0.328) * (listParms$meanBodyMass_g/1000) # 9.124 J/m
  
  # 2. Max distance Dmax that can be covered with max fat storage (Xmax) and flight range parameter c
  Dmax_km <- ((maxEnergy_kJ * 1000) / cot)/1000 # 4014.469 km (against 4328 km in silke's paper and 2500 km in Ulysses' thesis)
  
  # 3. The c parameter is just a calibration step that ensures that range maximum = Dmax (range max happening when x = MaxX and f = 1)
  # This allows to reduce distance between stopovers according to f: when f > 1, range will be < Dmax_km (range is defined in the bwdIteration in the cpp file)
  c_new <- Dmax_km / (1 - (1 + 1)^(-0.5)) 
  
  # 4. Flight speed (ground speed) provided by the user in m/s and here converted in km covered per day
  speed_kmDay <- listParms$speed_ms * 86.4
  
  
  # Fill and return a SDPMig object for the next steps
  new(
    "SDPMig", 
    Name = listParms$projectName,
    Init = list(
      MinT   = listParms$startDOY,
      MaxT   = listParms$endDOY,
      NSites = nrow(siteTab_sort),
      MaxX   = listParms$MaxX
    ),
    Species = list(
      B0    = listParms$B0,
      w     = w,
      xc    = xc,
      max_u = max_u,
      f     = listParms$f,
      c     = c_new,
      # Wind and flight speed can vary per site so might go into site parameters but for now we leave them here for consistency with previous models
      WindAssist = 0, #can be either one value or an array [nhex, nhex, ndays] 
      WindProb   = 1, #can be either one value or NULL
      speed = speed_kmDay,
      decError   = listParms$decError,
      ZStdNorm   = ZStdNorm,
      PStdNorm   = PStdNorm,
      xFTReward  = xFTReward,
      yFTReward  = yFTReward
    ),
    Sites = list(
      SitesID = hexIDs_sort,
      breedSites = breedSite,
      wintSites  = wintSite,
      crds = siteTab_sort[,3:4],
      dist = distMat,
      
      pred_a1 = pred_a1,
      pred_a2 = pred_a2,
      
      b0 = b0,
      b1 = b1,
      b2 = b2,
      
      expend  =  deeMat,
      gain    =  list(gain.x = gain.x, 
                      gain.y = gain.y, 
                      gain.p = gain.p)
    ),
    Results = list(
      FitnessMatrix   = NA,
      DecisionMatrix  = NA,
      ProbMatrix      = NA
    )
  )
}


##########################
### Backward iteration ###
##########################

bwdIteration <- function(obj, pbar = FALSE) {
  
  Init(obj@Init$MaxT - obj@Init$MinT, # int
       obj@Init$NSites, # int
       obj@Init$MaxX, # int
       obj@Species$w, # double
       obj@Species$xc, # double
       obj@Species$B0, # double
       obj@Sites$b0, # arma::vec
       obj@Sites$b1, # arma::vec
       obj@Sites$b2, # arma::vec
       obj@Sites$pred_a1, # double
       obj@Sites$pred_a2, # double
       obj@Species$c, # double
       obj@Species$speed, # double
       obj@Species$f, # double
       obj@Species$WindAssist, # arma::vec
       obj@Species$WindProb, # arma::vec
       obj@Species$ZStdNorm, # arma::vec
       obj@Species$PStdNorm, # arma::vec
       obj@Species$xFTReward, # arma::vec
       obj@Species$yFTReward, # arma::vec
       obj@Species$decError, # double
       as.matrix(obj@Sites$dist), # arma::mat
       obj@Sites$gain$gain.x, # arma::vec
       obj@Sites$gain$gain.y, # arma::mat
       obj@Sites$gain$gain.p, # arma::vec
       as.matrix(obj@Sites$expend) # arma::mat
  ) 
  
  out <- BackwardIteration(pbar = pbar)
  
  obj@Results$FitnessMatrix <- out[[1]]
  DM <- array(dim = c(dim(out[[2]]),2))
  DM[,,,1] <- out[[2]]
  DM[,,,2] <- out[[3]]
  obj@Results$DecisionMatrix <- DM
  PM <- array(dim = c(dim(out[[4]]),2))
  PM[,,,1] <- out[[4]]
  PM[,,,2] <- out[[5]]
  obj@Results$ProbMatrix <- PM
  
  obj
}


############################
#### Forward Simulation#####
############################

fwdSimulation <- function(model, NrInd, start_t, start_site, start_x) {
  
  InitSim(model@Init$MaxT, 
          model@Init$NSites, 
          model@Init$MaxX,
          model@Species$w,
          model@Species$xc,
          model@Species$B0,
          model@Sites$b0,
          model@Sites$b1,
          model@Sites$b2,
          model@Sites$pred_a1,
          model@Sites$pred_a2,
          model@Species$c,
          model@Species$speed,
          model@Species$f,
          model@Species$WindAssist,
          model@Species$WindProb,
          model@Species$ZStdNorm,
          model@Species$PStdNorm,
          model@Species$xFTReward,
          model@Species$yFTReward,
          model@Species$decError,
          model@Sites$dist,
          model@Sites$gain$gain.x,
          model@Sites$gain$gain.y,
          model@Sites$gain$gain.p,
          model@Sites$expend)
  
  
  x <- round(runif(NrInd, start_x[1], start_x[2]),0)
  
  if(length(start_site)>1 & length(start_site)<start_x[1]) {
    stop("start_site must have same length as numbers of individuals or a single site.")
  }
  if(length(start_site)==1) start_site <- rep(start_site, NrInd)
  
  SimOut = array(dim = c(length(x), 6, dim(model@Results$FitnessMatrix)[1]))
  
  ### First entry
  for(i in 1:dim(SimOut)[1]) {
    SimOut[i, ,start_t] <- c(start_t, start_site[i], x[i], 0, 0, 0)
  }
  
  
  ## SimOut: 1 = time, 2 = site, 3 = x, 4 = decision, 5 = flying, 6 = dead {
  for(time in 1:(dim(SimOut)[3]-1)) {
    
    for(ind in 1:dim(SimOut)[1]) {
      
      ## Not dead, not arrived, not flying
      if(!SimOut[ind, 6, time] & 
         sum(SimOut[ind, 2, ] >= nrow(model@Sites$crds), na.rm = T)<1 & !SimOut[ind, 5, time]) {
        
        ## Decision (when negative feeding, when positive it's the site to fly)
        if(runif(1) <  model@Results$ProbMatrix[SimOut[ind, 2, time], time, SimOut[ind, 3, time], 1]) {
          decision  <- model@Results$DecisionMatrix[SimOut[ind, 2, time], time, SimOut[ind, 3, time], 1]
        } else {
          decision  <- model@Results$DecisionMatrix[SimOut[ind, 2, time], time, SimOut[ind, 3, time], 2]
        }
        
        ## Action
        if(decision>=0) { ## Flying
          
          fl_help = simFlying(decision, time-1, SimOut[ind, 2, time]-1, SimOut[ind, 3, time])
          
          nextt = fl_help[1] + 1
          if(nextt<=time) nextt <- time+1
          if(nextt>dim(SimOut)[3]) time <- dim(SimOut)[3]
          
          nextx = fl_help[2] + 1
          if(nextx <  0) {
            nextx = 0
            dead  = 1 } else  dead = 0
          if(nextx > model@Init$MaxX) nextx = model@Init$MaxX
          
          SimOut[ind,,nextt] = c(nextt, decision+1, nextx, NA, 0, dead)
          if(nextt>(time+1)) SimOut[ind,5:6,(time+1):(nextt-1)] = cbind(1,0)
          if(SimOut[ind, 6, nextt])  SimOut[ind, 6, nextt:dim(SimOut)[3]] = 1 ## if dead make dead till the end
          
          if(SimOut[ind, 2, nextt]==nrow(model@Sites$crds)) {
            SimOut[ind,2:6, nextt:dim(SimOut)[3]] <- SimOut[ind, 2:6, nextt]
            SimOut[ind,1,   nextt:dim(SimOut)[3]] <- seq(nextt, dim(SimOut)[3])
          }
          
        } else { ## Feeding
          
          fo_help = simForaging(abs(decision+1.0), time-1, SimOut[ind, 2, time]-1, SimOut[ind, 3, time])
          
          newx = fo_help[1]+1
          dead = fo_help[2]
          if(newx<=0) dead = 1
          
          if(newx > model@Init$MaxX) newx = model@Init$MaxX
          
          SimOut[ind,,time+1] = c(time+1, SimOut[ind, 2, time], newx, abs(decision+1.0), 0, dead)
          
          if(SimOut[ind, 6, time+1])  SimOut[ind, 6, (time+1):dim(SimOut)[3]] = 1 ## if dead make dead till the end
          
        }
        
      }
      
    } ## Ind loop
  } ## time loop
  
  SimOut       <- SimOut[,-5,] 
  SimOut[,2,] <- SimOut[,2,]-1
  
  # Add names to the columns in the second dimension for interpretability
  # The final SimOut obj has dimensions: n.Ind, 5, n.Days. 
  # The 5 info in dim 2 are: model day, site, x, fly/not, death. 
  # Fly/not takes either the ID of the next site, or 0 for alive individuals not using the site or NA for dead individuals.
  # Death takes 0 for alive and 1 for dead.
  dimnames(SimOut)[[2]] <- c("moday","siteIdx","x","flyProb","death")
  
  # Replace the probability matrices with the output simulation. This way all parms information are kept in the output (which will still be a sdpMig object)
  model@Results <- list(Fwd = SimOut) 
  
  return(model)
}

#############################
#### Plotting functions #####
#############################

### Plot body condition of single simulation over time
condProfile <- function(simuObj) {
  
  trkS  <- simuObj@Results$Fwd[,2,]+1
  trkX  <- simuObj@Results$Fwd[,3,]
  trkD  <- simuObj@Results$Fwd[,5,]
  
  cls <- rainbow(max(trkS, na.rm = T))[sample(1:max(trkS, na.rm = T))]
  
  opar <- par(mfrow = c(1,1))
  matplot(simuObj@Init$MinT:simuObj@Init$MaxT, t(trkX), type = "n", 
          xaxt = "n", xlab  = "", ylab = "Body condition", las = 1)
  for(i in 1:dim(simuObj@Results$Fwd)[1]) {
    if(all(trkD[i,]==0)){
      tmp <- cbind(simuObj@Init$MinT:simuObj@Init$MaxT, trkS[i,], trkX[i,])[!is.na(trkS[i,]),]
      points(tmp[,1], tmp[,3], type = "l", col = "grey90")
      points(tmp[,1], tmp[,3], pch = 16, cex = 0.5, col = cls[tmp[,2]])
    }
  }
  axis(1, at = seq(simuObj@Init$MinT,simuObj@Init$MaxT, by = 15), labels = format(as.POSIXct("2012-01-01")+seq(simuObj@Init$MinT,simuObj@Init$MaxT, by = 15)*24*60*60, "%b-%d"))
  par(opar)
}

### Plot migratory routes of single simulation as network

pltNetworkGG <- function(simuObj, map = europe_laea) {
  
  require(data.table)
  require(ggplot2)
  
  sitesCrds <- simuObj@Sites$crds
  dead      <- apply(simuObj@Results$Fwd[,5,], 1, sum, na.rm = TRUE)
  trkS      <- simuObj@Results$Fwd[dead < 1, 2, ] + 1
  
  diagT <- apply(trkS, 1, function(x) {
    tt <- as.data.frame(table(x[!is.na(x)]))
    m  <- merge(data.frame(Var1 = 1:nrow(sitesCrds)),
                tt, all.x = TRUE)
    m[nrow(m),2] <- NA
    ifelse(is.na(m[,2]), 0, m[,2])
  })
  
  tmp.out <- matrix(0,
                    ncol = nrow(sitesCrds),
                    nrow = nrow(sitesCrds))
  
  diag(tmp.out) <- apply(diagT, 1, sum)
  
  transT <- as.data.frame(
    table(
      apply(
        do.call(
          rbind, apply(trkS, 1, function(x) {
            tmp01 <- cbind(
              x[!is.na(x)][-sum(!is.na(x))],
              x[!is.na(x)][-1]
            )
            
            tmp01[
              tmp01[,1] != tmp01[,2] &
                !is.na(tmp01[,1]) &
                !is.na(tmp01[,2]),
            ]
          })
        ),
        1,
        function(y) paste(y, collapse = "_")
      )
    )
  )
  
  trans <- cbind(
    t(apply(transT, 1,
            function(z)
              as.numeric(strsplit(z[1], "_")[[1]]))),
    transT[,2])
  
  ## Project coordinates once
  pts_sf <- st_as_sf(as.data.frame(sitesCrds),
                     coords = c("x", "y"), crs = 4326) |>
    st_transform(st_crs(map))
  
  coords <- cbind(sitesCrds, st_coordinates(pts_sf))
  
  names(coords)[(ncol(coords)-1):ncol(coords)] <- c("X", "Y")
  
  ## Segment dataframe
  seg_df <- data.frame(
    x    = coords$X[trans[,1]],
    y    = coords$Y[trans[,1]],
    xend = coords$X[trans[,2]],
    yend = coords$Y[trans[,2]],
    freq = trans[,3])
  
  ## Point dataframe
  site_df <- data.frame(X = coords$X, Y = coords$Y, d = diag(tmp.out))
  
  site_df$start <- FALSE
  site_df$start[1] <- TRUE
  site_df$end <- FALSE
  site_df$end[nrow(site_df)] <- TRUE
  
  ggplot() +
    geom_sf(data = map, fill = "black", color = "lightgrey") +
    geom_segment(data = seg_df, 
                 aes(x = x, y = y, xend = xend, yend = yend, linewidth = freq), 
                 colour = "orange") +
    scale_linewidth(range = c(0.1, 3), guide = "none") +
    geom_point(data = subset(site_df, !start & !end),
               aes(X, Y, size = d), shape = 21,
               fill = adjustcolor("orange", alpha.f = 0.55),
               colour = "white", stroke = 0.3) +
    geom_point(data = subset(site_df, start), aes(X, Y),
               shape = 23, size = 5, fill = adjustcolor("turquoise", alpha.f = 1),
               colour = "white", stroke = 0.3) +
    geom_point(data = subset(site_df, end), aes(X, Y),
               shape = 23, size = 5, fill = adjustcolor("orchid", alpha.f = 1),
               colour = "white", stroke = 0.3) +
    scale_size(range = c(0.1, 8), guide = "none") +
    coord_sf() +
    theme_void()
}

### Comparison/validation of parameters between sets of simulations (plots and summary table)

make_density <- function(data, var, xlab, group, cols) {
  require(ggplot2)
  ggplot(st_drop_geometry(data), aes(x = .data[[var]], fill = .data[[group]], colour = .data[[group]])) +
    geom_density(alpha = 0.35, linewidth = 0.5) +
    scale_fill_manual(values = cols)   +
    scale_colour_manual(values = cols) +
    labs(x = xlab, y = "Density") +
    theme_minimal() +
    theme(legend.position = "none")
}

make_bar <- function(data, var, xlab, group, cols) {
  require(ggplot2)
  data %>%
    st_drop_geometry() %>%
    count(.data[[group]], val = .data[[var]]) %>%
    group_by(.data[[group]]) %>%
    mutate(prop = n / sum(n)) %>%
    ggplot(aes(x = val, y = prop, fill = .data[[group]])) +
    geom_col(position = "dodge", alpha = 0.8, width = 0.7) +
    scale_fill_manual(values = cols) +
    labs(x = xlab, y = "Proportion") +
    theme_minimal() +
    theme(legend.position = "none")
}

make_table <- function(data, vars, group) {
  require(ggplot2)
  require(tidyr)
  data %>%
    group_by(.data[[group]]) %>%
    summarise(across(all_of(vars), list(
      min    = \(x) min(x,    na.rm = TRUE),
      max    = \(x) max(x,    na.rm = TRUE),
      mean   = \(x) mean(x,   na.rm = TRUE),
      median = \(x) median(x, na.rm = TRUE),
      sd     = \(x) sd(x,     na.rm = TRUE)
    ), .names = "{.fn}_{.col}"), .groups = "drop") %>%
    pivot_longer(-all_of(group), 
                 names_to    = c("stat", "metric"),
                 names_pattern = "^(min|max|mean|median|sd)_(.*)",
                 values_to = "value") %>%
    pivot_wider(names_from = c(all_of(group), "stat"),  values_from = "value")
}



#############################
#### Output functions #####
#############################

# This function calculates the following measures: 
# individuals per site, mortality, staging times per site, body condition per site, fuelling rates per site, migratory path, successful migr prop and death prop 
# !Reminder for the calculations! simuObj@Results$Fwd has dimensions:
# n.Ind, 5, n.Days. The 5 info in dim 2 are: day, site, x, fly/not, death. 
# Fly/not takes either the ID of the next site, or 0 for alive individuals not using the site or NA for dead individuals.
# Death takes 0 for alive and 1 for dead.

computeOutputs <- function(simuObj){
  
  require(data.table)
  require(dplyr)
  
  combID <- paste0("w",simuObj@Sites$wintSites,"_b",simuObj@Sites$breedSites)
  Ninds <- dim(simuObj@Results$Fwd)[1]
  Ndays <- dim(simuObj@Results$Fwd)[3]
  
  # 1. Migratory routes and staging time per individual. Each row of this matrix represents an individual route, indicating where the individual stopped (which hex ID) and for how long
  migRoutes <- rbindlist(lapply(1:Ninds, function(i){ # for each individual
    
    tmp <- as.data.frame(t(simuObj@Results$Fwd[i,,]))
    usedSites <- unique(tmp$siteIdx[!is.na(tmp$siteIdx)]) # index of stopover sites (C++)
    
    if(length(usedSites)>0){ # exclude individuals who died before starting migration
      
      indRoute <- rbindlist(lapply(usedSites, function(s){
        
        entranceDay <- min(tmp$moday[which(tmp$siteIdx==s)])
        exitDay <- max(tmp$moday[which(tmp$siteIdx==s)])
        entranceX <- tmp$x[which(tmp$moday == entranceDay)]
        exitX <- tmp$x[which(tmp$moday == exitDay)]
        
        indSite <- data.frame(combID=combID,
                              ind = i,
                              indID = paste(combID, i, sep = "_"),
                              used_siteIdx = s+1,
                              used_hexID = as.character(simuObj@Sites$SitesID[s+1]), # extract corresponding to hexID
                              arrivalDOY = entranceDay + simuObj@Init$MinT, # convert to day of year
                              departureDOY = exitDay + simuObj@Init$MinT,
                              stopDays = exitDay - entranceDay,
                              arrivalCond = entranceX,
                              departureCond = exitX,
                              stopGain = exitX - entranceX)
        # add death information
        if(sum(tmp$death[which(tmp$siteIdx==s)]) > 0){ # if the individual die in this site death=TRUE
          indSite$death <- T
        }else{
          indSite$death <- F
        }
        
        return(indSite)
      }))
      # add whether migration was successful (if last stopover is the breeding ground)
      if(indRoute$used_hexID[nrow(indRoute)] == simuObj@Sites$breedSites){
        indRoute$success <- TRUE
      }else{
        indRoute$success <- FALSE
      }
      
      return(indRoute)}
  }))
  # add staging days of breed and wint to different columns
  migRoutes[which(migRoutes$used_hexID == simuObj@Sites$wintSites),
            `:=`(
              wintDays = stopDays,
              wintGain = stopGain)] # what they stored at wintering ground before starting migration
  migRoutes[which(migRoutes$used_hexID == simuObj@Sites$breedSites),
            breedDays := stopDays]
  # remove entrance information (day and gain) from wintering ground because meaningless (they reflect our decision on model initial parameter rather then biology)
  migRoutes[which(migRoutes$used_hexID %in% simuObj@Sites$wintSites), 
            c("arrivalDOY","arrivalCond","stopGain","stopDays")] <- NA
  # remove exit information (day and gain) from breeding ground because meaningless (they reflect our decision on model end rather than any biology)
  migRoutes[which(migRoutes$used_hexID %in% simuObj@Sites$breedSites), 
            c("departureDOY","departureCond","stopGain","stopDays")] <- NA
  
  # 2. Migration summary per individual. 
  # First the timing, start end and duration. Site 0 is the wintering ground (first site index in C++) and site Nsites - 1 is the breeding ground (again because C++ starts counting from 0)
  MigrationPerInd <- data.frame(combID = combID,
                                ind = 1:Ninds, # one row per ind
                                indID = paste(combID, 1:Ninds, sep = "_"),
                                migrStartDOY = NA, migrEndDOY = NA, migrDurationDays = NA) 
  
  for (i in 1:Ninds) { # for each individual i
    tmp <- t(simuObj@Results$Fwd[i, , ])
    migrStart <- max(tmp[,"moday"][tmp[,"siteIdx"] == 0], na.rm=T) # last model day at which the individual is at wintering ground
    if (any(tmp[,"death"] == 1) || # if individual is dead or
        all(tmp[,"siteIdx"] != (simuObj@Init$NSites - 1), na.rm = TRUE) # unsuccessful (did not reach breeding ground)
    ){
      migrEnd <- NA # migration did not end (NA)
    }else{
      migrEnd <- min(tmp[,"moday"][tmp[,"siteIdx"] == simuObj@Init$NSites - 1], na.rm=T) # first model day at which breeding ground is reached
    }
    migrDuration <- migrEnd - migrStart
    # store info in day of year instead of model days
    MigrationPerInd[i, 4:6] <- cbind(migrStart+simuObj@Init$MinT, 
                                     migrEnd+simuObj@Init$MinT,
                                     migrDuration)  
  }
  # Then add N. of stopovers, avg staging, death and success per individual excluding start/end (one value per ind)
  NstopsPerInd <- migRoutes %>%
    group_by(combID, ind) %>%
    summarise(success = first(success),
              death = any(death == 1),
              Nstops = n_distinct(used_siteIdx) - if_else(success, 2L, 1L), # drop 2 stops (breed and wint) for successful individuals but only 1 (wint) for unsuccessful ones
              meanStag = mean(stopDays[-1], na.rm = TRUE),
              .groups = "drop")
  # Merge both
  migPerInd <- merge(MigrationPerInd, NstopsPerInd, by=c("combID","ind"))
  
  # 3. Migration summary stats. Prop of mortality and of successful migration (all that arrived within the expected window)
  SuccessProp <- sum(migPerInd$success)/Ninds #success/nInds
  MortalityProp <- sum(migPerInd$death)/Ninds #deaths/nInds, proportion of deaths in the population
  # same as: sum(apply(simuObj@Results$Fwd[,5,], 1, function(x) any(x==1)))/Ninds
  migrStats <- data.frame(combID = combID,
                          Ninds = Ninds,
                          SuccessProp = SuccessProp,
                          MortalityProp = MortalityProp)
  
  # 4. Migration summary per site (hex ID). Unused sites are excluded.
  migPerSite <- group_by(migRoutes, used_siteIdx, used_hexID) %>% 
    summarise(Ninds = n(),
              mean_wintDays = mean(wintDays, na.rm = T),
              mean_wintGain = mean(wintGain, na.rm = T),
              mean_stopDays = mean(stopDays, na.rm = T),
              mean_breedDays = mean(breedDays, na.rm = T),
              mean_arrivalDOY = mean(arrivalDOY, na.rm = T),
              mean_departureDOY = mean(departureDOY, na.rm = T),
              mean_arrivalCond = mean(arrivalCond, na.rm = T),
              mean_departureCond = mean(departureCond, na.rm = T),
              mean_stopGain = mean(stopGain, na.rm = T),
              mean_mortality = mean(as.numeric(death)),
              .groups = "drop")
  
  # 5. Body condition data.frame per individual day and site (hex ID) 
  condDf <- rbindlist(lapply(1:Ninds, function(i){
    tmp <- simuObj@Results$Fwd[i,,]
    day_x <- as.data.frame(t(tmp[c("moday","x","death"),]))
    day_x <- day_x[cumsum(day_x$death) <= 1, ]   # keep rows up to and including death (removing days after death)
    data.frame(combID = combID,
               ind = i,
               indID = paste(combID, i, sep = "_"),
               modDay = day_x$moday,
               doy = day_x$moday + simuObj@Init$MinT,
               bodyCond = day_x$x,
               delta_BodyCond = day_x$x - lag(day_x$x, 1),
               death = day_x$death)
  }))
  
  # Return all measurements per scenario
  return(list(combID = combID, 
              migStats = migrStats,
              migRoutes = migRoutes, 
              migPerInd = migPerInd,
              migPerSite = migPerSite,
              BodyCondProf = condDf))
}










# simNetwork <- function(simu, model, crds_ind, map = eaafMap, plot = FALSE) {
#   
#   dead_orNA <- apply(cbind(apply(simu[,5,], 1, sum, na.rm = T)==0, 
#                            apply(simu[,2,], 1, function(x) any(x[!is.na(x)]==max(simu[,2,], na.rm = T)))), 1, 
#                      function(x) all(x))
#   trkS      <- simu[dead_orNA,2,]+1
#   
#   sites <- tibble(index = model@Sites$index, ts = apply(trkS, 1, function(y) {
#     seqTmp <- y[!is.na(y)][max(which(y[!is.na(y)]==1)):min(which(y[!is.na(y)]==max(y[!is.na(y)])))]
#     tibble(sites = seqTmp) %>% group_by(sites) %>% summarise(count = n()) %>%
#       full_join(tibble(sites = 1:length(model@Sites$index)), by = "sites") %>% arrange(sites) %>%
#       mutate(count = ifelse(is.na(count), 0, count)) %>% pull(count)
#   }) %>% rowSums()) %>% group_by(index) %>% summarise(ts = sum(ts)) %>%
#     full_join(tibble(index = 1:nrow(crds_ind)), by = "index") %>%
#     arrange(index) %>% mutate(ts = ifelse(is.na(ts), 0, ts)) %>% pull(ts)
#   
#   if(plot) {
#     opar <- par(mar = c(0,0,0,0))
#     plot(eaafMap$grid, col = "grey80", lty = 3)
#     plot(eaafMap$map, col = "grey90", border = "grey60", add = T)
#     plot(eaafMap$bbox, add = T, border = "grey60")
#   }
#   
#   ### Relative Site Use
#   relSite <- do.call("rbind", lapply(1:nrow(trkS), function(y) {
#     tmp <- trkS[y,]
#     tibble(ts = model@Sites$index[tmp[!is.na(tmp)][(max(which(tmp[!is.na(tmp)]==1))+1):(min(which(tmp[!is.na(tmp)]==max(tmp[!is.na(tmp)])))-1)]]) %>%
#       group_by(ts) %>% summarise(days = n()) %>% arrange(desc(days)) %>% mutate(site = 1:nrow(.)) %>%
#       select(site, ts, days)
#   }))
#   migrPhen <- do.call("rbind", lapply(1:nrow(trkS), function(y) {
#     tmp <- trkS[y,]
#     tibble(dep = (max(which(!is.na(tmp) & tmp==min(tmp, na.rm=T))+1))+model@Init$MinT,
#            arr = (min(which(!is.na(tmp) & tmp==max(tmp, na.rm=T))))+model@Init$MinT)
#   })) %>% mutate(dur = arr-dep)
#   
#   
#   
#   transTab <- apply(trkS, 1, function(y) {
#     tmp00 <- model@Sites$index[y[!is.na(y)][max(which(y[!is.na(y)]==1)):min(which(y[!is.na(y)]==max(y[!is.na(y)])))]]
#     tmp01 <- cbind(tmp00[!is.na(tmp00)][-sum(!is.na(tmp00))], tmp00[!is.na(tmp00)][-1])
#     as.data.frame(tmp01[tmp01[,1]!=tmp01[,2] & !is.na(tmp01[,1]) & !is.na(tmp01[,2]),])
#   }) %>% do.call("rbind",.) %>% setNames(c("a", "b")) %>% group_by(a, b) %>% summarise(count = n())
#   
#   trans <- matrix(0, ncol = length(sites), nrow = length(sites))
#   trans[transTab %>% select(a, b) %>% as.matrix()] <- transTab %>% pull(count)
# 
#   
#   if(plot) {
#     transT   <- tibble(a = rep(1:nrow(trans), nrow(trans)), b = rep(1:nrow(trans), each = nrow(trans)), count = c(trans)) %>% filter(count > 0)
#     trans_sf <- lapply(1:nrow(transT), function(t) st_linestring(crds_ind[unlist(transT[t,1:2]),])) %>% st_sfc() %>% st_set_crs(st_crs(eaafMap$map)) %>%
#       st_sf() %>% mutate(trans = transT[,3]) %>% select(trans, geometry)
#     plot(trans_sf %>% st_geometry(), add = T, lwd = approx(c(1, 100), c(0, 5), transT$count)$y, col = adjustcolor("grey30", alpha.f = 0.7))
#     points(crds_ind, cex = approx(range(sites), c(0, 8), sites)$y,
#            pch = 21, lwd = 1.5, col = "grey40", bg = adjustcolor("orange", alpha.f = 0.6))
#   }
#   
#   diag(trans) <- sites
#   
#   ### fuelling
#   trkS  <- simu[,2,]+1
#   trkX  <- simu[,3,]
#   trkD  <- simu[,5,]
#   
#   fueling <- lapply(1:dim(simu)[1], function(i) {
#     if(all(trkD[i,]==0)){
#       cbind(i, model@Init$MinT:model@Init$MaxT, trkS[i,], trkX[i,])[!is.na(trkS[i,]),]
#     } else NULL
#   }) %>% Reduce("rbind", .)
#   
#   
#   fitness <- lapply(1:dim(simu)[1], function(x) {
#     if(all(simu[x,5,]==0)) {
#       t <- min(which(simu[x,2,]==max(simu[x,2,], na.rm = T)))
#       matrix(c(t, simu[x,3,t]), ncol = 2)
#     } else NULL
#     }) %>% do.call("rbind",.) %>% apply(., 1, function(f) model@Results$FitnessMatrix[f[1], dim(model@Results$FitnessMatrix)[2],f[2]])
#   
#   
#   list(trans, relSite, migrPhen, fueling, fitness)
#   
# } 

