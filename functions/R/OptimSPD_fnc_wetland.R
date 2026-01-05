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

### Michaelis–Menten / hyperbolic saturation function
# Y is the maximum value of the variable of interest (e.g. gain)
# x is the variable Y varies upon (e.g. daylength, suitability)
# K is a curvature parameter that governs saturation rate when x is between 0 and 1
# otherwise K is the value of x at which y reaches the half-max 
mm_fun <- function(Ymax, x, K) Ymax * (x / (x + K))
# K is based on the fraction yf of Y that we expect to reach at the x0 value of x
# e.g. we expect to reach 0.95% of the maximum gain at suitability x0=0.5
findK <- function(x0, yf) x0 * (1 - yf) / yf  

#############################
### Create Parameters Tab ###
#############################

# Major changes are marked with !!!

# Inputs are: 
# 1. a data.frame or matrix of 5 columns in this order: 1st col = sites ID (hexagons id), 2nd col = opening day as day of the year (Jan 1st being day 1), 3rd and 4th columns are coordinates (long and lat in decimal degrees), 5th column is a suitability layer between 0 and 1. Potential additional columns are ignored.
# 2. two integers indicating the site ID (hexagon ID) of the wintering and breeding site. These should be among the IDs in the first column of data.
# 3. the list of parameters from previous steps,
# 4. optionally we can have matrices of environmental values (e.g. windTab) with dimensions nrow = Nsites and ncol = doy to make distance dependent on wind

makeSiteTab <- function(siteTab,
                        wintSite, 
                        breedSite,
                        listParms,
                        BMR_file = "inputData/BMRcalculation_perOrder_McNab2008.csv",
                        pred_a1  = 2, pred_a2  = 2,
                        b0 = 1e-04, b1 = 0.001, b2 = 1e-04,
                        stoch = 1
){
  require(geosphere)
  require(terra)
  require(sf)
  require(data.table)
  
  #___ Species parameters to leave constant________
  w         <- 0.028  ## Parameters for sigmoidal TR function: slope
  xc        <- 55     ## Parameters for sigmoidal TR function: body condition (x) at inflection point of sigmoid
  max_u     <- 1.0    ## maximum foraging intensity (max prop. of day spent foraging)
  MaxReward <- 2      ## value within arrival window, for terminal reward calculation 
  
  #___Sorting sites (breeding site last, wintering site first)___
  hexIDs    <- siteTab[,1]
  hex_idx <- 1:length(hexIDs) 
  hex_idx_new   <- c(which(hexIDs==wintSite), hex_idx[hex_idx!=which(hexIDs==wintSite) & hex_idx!=which(hexIDs==breedSite)], which(hexIDs==breedSite))
  hexIDs_sort <- hexIDs[hex_idx_new]
  siteTab_sort <- siteTab[hex_idx_new,]
  
  #___Model Days___
  modelRun <- length(listParms$startDOY:listParms$endDOY)
  
  #___ Basal Metabolic Rate__________ 
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
    
    #___Daily Energy Expenditure DEE / Field Metabolic Rate FMR__________
    # In Lisovski et al. PNAS, DEE is assumed to be 1.5 * BMR
    # But according to Nagy 2005: "there is no trend in present data suggesting that FMR slopes are identical to BMR slopes." (therefore that FMR is a constant multiple of BMR)
    # FMR in kJ/day = 10.5 * (mass in g ^0.681), For birds based on 95 species, from Nagy 2005 "Field metabolic rate and body size"
    FMR_kjday <- 10.5 * (listParms$meanBodyMass_g^0.681) #This gives us 2181.735 kJ/day
    
    # This must be converted in units of X (body condition). How many X does FMR corresponds to? 
    # We obtain it by dividing FMR for the energy equivalent of a unit of X:
    # !!! dee_inX <- FMR_kjday / energyPerX # equivalent to 5.96 units of X
    # Alternatively we could use 1.5 * BMR as approximation of DEE
    dee_inX <- (BMR_kJday * 1.5) / energyPerX # equivalent to only 3.28, repeated for all sites and all days
    deeMat <- matrix(dee_inX, nrow = nrow(siteTab_sort), ncol = modelRun)
  }
  
  #___ Gain / maximum Daily Metabolised Energy Intake (DMEmax)__________
  
  # Calculate the maximum Daily Metabolised Energy Intake (DMEmax) according to Kirkwood 1983 “A limit to metabolisable energy intake in mammals and birds”
  DMEmax <- 1713 * ((listParms$minBodyMass_g/1000) ^ 0.72) 
  DMEmax_inX <- DMEmax / energyPerX # gain = DME in units of body condition X, equal to 7.9
  
  # Opening days in model days. It cannot be negative, so if the startDay is after the opening day, we declare the site as open from the 1st day of model run (day 1)
  openDOY <- siteTab_sort[,2]
  openMod <- (openDOY - listParms$startDOY)
  
  ### Gain per day is dee_inX * 1.1 before opening day
  # !!! Martina changed gain.y function
  gain.x <- seq(1:modelRun)
  gain.p  <- rep(stoch, modelRun) # stoch is the same all across
  gain.y <- matrix(0, nrow = nrow(siteTab_sort), ncol = modelRun) #we could use 0 or a minimum dee_inX * 1.1
  ### option 1: Gain is max after opening day, one value per hex
  # for (i in 1:nrow(gain.y)) { # calculate gain per hexagon
  #   if (openMod[i] <= 0) {
  #     gain.y[i, ] <- DMEmax_inX
  #   }else if (openMod[i] < modelRun) {
  #     gain.y[i, (openMod[i]:modelRun)] <- DMEmax_inX
  #   }}
  ### option 2: Gain depends on suitability (wetland prop) after opening day (suitability is static, so one value per hexagon)
  for (i in 1:nrow(gain.y)) { # calculate gain per hexagon
    gainVal <- mm_fun(Ymax=DMEmax_inX, x=siteTab_sort[i,5], K=findK(x0=0.6, yf=0.95)) # at suit 0.6 gain is 0.95% of max 
    if (openMod[i] <= 0) {
        gain.y[i, ] <- gainVal
      }else if (openMod[i] < modelRun) {
        gain.y[i, (openMod[i]:modelRun)] <- gainVal
      }}
  ### option 3: Gain dependent on daylength after opening day (gain per hexagon changing daily after opening)
  # Calculate daylength to derive energy intake per hour of light and adjust gain.y relative to number of foraging hours
  # daylength <- t(apply(siteTab_sort, 1, function(x) {
  #   round(geosphere::daylength(as.numeric(x["y"]), 1:366), 1)
  # }))
  #ggplot(siteTab_sort, aes(x = x, y = y, fill = daylength[,172])) + geom_point(shape=21, size=3) + coord_equal()
  #for (i in 1:nrow(gain.y)) { # calculate gain per hexagon
    #if (openMod[i] <= 0) {
      #h <- daylength[i, (listParms$startDOY:listParms$endDOY)]
      #gain.y[i, ] <- mm_fun(Ymax=DMEmax_inX, x=h, K=round(mean(daylength))) # Michaelis–Menten / hyperbolic saturation. where K is the hour of the day at which we expect gain to reach half maximum (here mean of the daylength across all sites)
      #gain.y[i, ] <- DMEmax_inX * (h / max(daylength)) #linear increase of gain with daylength
      # if site opens before the END of model run, calculate new gain from respective opening day (in mod day and doy)
    #} else if (openMod[i] < modelRun) {
      #h <- daylength[i, (openDOY[i]:(listParms$endDOY + 1))]
      #gain.y[i, (openMod[i]:modelRun)] <- mm_fun(Ymax=DMEmax_inX, x=h, K=round(mean(daylength))) 
      #gain.y[i, (openMod[i]:modelRun)] <- DMEmax_inX * (h / max(daylength)) 
    #}}
  
  # Fix gain of wintering site always to maximum, assuming the bird knows the wintering ground and can always find resources
  gain.y[which(hexIDs_sort==wintSite),] <- DMEmax_inX
  # range(gain.y)
  # plot gain at the beginning middle and end of model run
  # ggplot(siteTab_sort, aes(x = x, y = y, fill = gain.y[,1])) + geom_point(shape=21, size=3) + coord_equal()
  # ggplot(siteTab_sort, aes(x = x, y = y, fill = gain.y[,70])) + geom_point(shape=21, size=3) + coord_equal()
  # ggplot(siteTab_sort, aes(x = x, y = y, fill = gain.y[,120])) + geom_point(shape=21, size=3) + coord_equal()
  
  # Accuracy (foraging success, sampled from a normal distribution)
  ZStdNorm = c(-2.5, -2.0, -1.5, -1.0, -0.5,  0.0,  0.5,  1.0,  1.5,  2.0,  2.5)
  PStdNorm = c(0.0092, 0.0279, 0.0655, 0.1210, 0.1747, 0.2034, 0.1747, 0.1210, 0.0655, 0.0279, 0.0092)
  
  #___ Terminal Reward__________
  # Calculate terminal reward for breeding sites (nrow = N sites, ncol defining the "curve shape" of the terminal reward):
  # x = Model day at which each site opens, y = reward. Breeding site is the last in the vector!
  xFTReward <- c(1, 
                 openMod[length(openMod)]-1, openMod[length(openMod)],
                 openMod[length(openMod)] + listParms$arrWindow, openMod[length(openMod)] + listParms$arrWindow + 1, 
                 length(gain.x))
  yFTReward <- c(0, 0, MaxReward, MaxReward, 0, 0)
  # plot(xFTReward, yFTReward, type="l")
 
  # ___ Predation Rates______
  b0 <- rep(b0, listParms$NSites) # Background predation risk (e.g. number of predators in the area)
  b1 <- rep(b1, listParms$NSites) # Body reserve dependent risk (e.g. the fatter the less mobile to escape)
  b2 <- rep(b2, listParms$NSites) # Foraging intensity dependent risk (e.g. the more they feed the less they see predators)
  # !!! Tested but excluded because did not make a big difference, reduced predation at wintering ground:
  # For wintering site (first value), set predation to three orders of magnitude less
  # b0[1] <- b0[1] * 0.001
  # b1[1] <- b1[1] * 0.001
  # b2[1] <- b2[1] * 0.001
  
  #___ Distances between all sites__________
  # Matrix of distances between all sites, coords are assumed to be in decimal degrees
  distMat <- as.matrix(terra::distance(x = siteTab_sort[,3:4], lonlat = TRUE, method = "geo", unit="km"))
  # plot distances from wintering ground (first site)
  # ggplot(siteTab_sort, aes(x = x, y = y, fill = distMat[,1])) + geom_point(shape=21, size=3) + coord_equal()
  
  
  #___Flight parameters_______
  # 1. Flight cost
  cot <- exp(1.587) * ((listParms$meanBodyMass_g/1000) ^ -0.328) * (listParms$meanBodyMass_g/1000) # 9.124 J/m
  
  # 2. Max distance Dmax that can be covered with max fat storage (Xmax) and flight range parameter c
  Dmax_km <- ((maxEnergy_kJ * 1000) / cot)/1000 # 4014.469 km (against 4328 km in silke's paper and 2500 km in Ulysses' thesis)
  # !!! Test to make less distance possible
  #Dmax_km <- Dmax_km * 0.3
  
  xf      <- listParms$f * listParms$MaxX # fraction xf of the maximum fuel load xmax dedicated to flight, change based on f
  c       <- Dmax_km/(1 - (1 + xf/listParms$MaxX)^(-0.5)) # flight range parameter 13706.25 (against 14776 in Silke's paper and 8535 in Ulysses' thesis)
  
  # 3. Flight speed is provided by the user in m/s and here converted in km covered per day
  speed_kmDay <- listParms$speed_ms * 86.4 # We obtain 1088.64 km/day compared to 1440 km/day in Silke's paper but they used airspeed. I think we should use ground speed.
  
  # Fill and return a SDPMig object for the next steps
  new(
    "SDPMig", 
    Name = listParms$projectName,
    Init = list(
      MinT   = listParms$startDOY,
      MaxT   = listParms$endDOY,
      NSites = listParms$NSites,
      MaxX   = listParms$MaxX
    ),
    Species = list(
      B0    = listParms$B0,
      w     = w,
      xc    = xc,
      max_u = max_u,
      f     = listParms$f,
      c     = c,
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
  
  SimOut
}

#############################
#### Plotting functions #####
#############################

### Plot body condition over time
condProfile <- function(simu, model) {
  trkS  <- simu[,2,]+1
  trkX  <- simu[,3,]
  trkD  <- simu[,5,]
  
  cls <- rainbow(max(trkS, na.rm = T))[sample(1:max(trkS, na.rm = T))]
  
  opar <- par(mfrow = c(1,1))
  matplot(model@Init$MinT:model@Init$MaxT, t(trkX), type = "n", 
          xaxt = "n", xlab  = "", ylab = "Body condition", las = 1)
  for(i in 1:dim(simu)[1]) {
    if(all(trkD[i,]==0)){
      tmp <- cbind(model@Init$MinT:model@Init$MaxT, trkS[i,], trkX[i,])[!is.na(trkS[i,]),]
      points(tmp[,1], tmp[,3], type = "l", col = "grey90")
      points(tmp[,1], tmp[,3], pch = 16, cex = 0.5, col = cls[tmp[,2]])
    }
  }
  axis(1, at = seq(model@Init$MinT,model@Init$MaxT, by = 15), labels = format(as.POSIXct("2012-01-01")+seq(model@Init$MinT,model@Init$MaxT, by = 15)*24*60*60, "%b-%d"))
  par(opar)
}

### Plot migratory routes as network
pltNetwork <- function(simu, model, map = europe_laea) {
  
  sitesCrds <- model@Sites$crds
  dead      <- apply(simu[,5,], 1, sum, na.rm = T)
  trkS      <- simu[dead<1,2,]+1
  
  diagT <-   apply(trkS, 1, function(x) {
    tt <- as.data.frame(table(x[!is.na(x)]))
    m  <- merge(data.frame(Var1 = 1:nrow(sitesCrds)), tt, all.x = T)
    m[nrow(m),2] <- NA
    ifelse(is.na(m[,2]), 0, m[,2])
  })
  tmp.out   <- matrix(0, ncol = nrow(sitesCrds), nrow = nrow(sitesCrds)) 
  diag(tmp.out) <- apply(diagT, 1, sum)
  
  transT <- as.data.frame(table(apply(do.call("rbind", apply(trkS, 1, function(x) {
    tmp01 <- cbind(x[!is.na(x)][-sum(!is.na(x))], x[!is.na(x)][-1])
    as.data.frame(tmp01[tmp01[,1]!=tmp01[,2] & !is.na(tmp01[,1]) & !is.na(tmp01[,2]),])
  })), 1, function(y) paste(y, collapse = "_"))))
  trans <- cbind(t(apply(transT, 1, function(z) as.numeric(strsplit(z[1], "_")[[1]]))), transT[,2])
  
  
  lwd.f <- approxfun(c(0, max(trans[,3])), c(0.1, 5), rule = 3)
  
  opar <- par(mfrow = c(1,1))
  plot(map %>% st_geometry(), col = "grey90", border = "grey50")
  pr <- st_as_sf(model@Sites$crds %>% as_data_frame(), coords = c("x", "y")) %>% st_set_crs(4326) %>% st_transform(st_crs(map)) %>% st_coordinates()
  segments(pr[trans[,1],1], pr[trans[,1],2], pr[trans[,2],1], pr[trans[,2],2], lwd = lwd.f(trans[,3]), col = "grey20")
  
  site01 <- cbind(sitesCrds, d = diag(tmp.out))
  cex.f <- approxfun(x = c(1, max(site01[-1,3])), y = c(0.1, 15), rule = 3)
  
  cex <- cex.f(site01[,3])
  cex[model@Sites$start] <- 2
  
  pch = rep(21, nrow(site01))
  pch[model@Sites$start] <- 23
  
  pts <- st_as_sf(site01[,1:2] %>% as_data_frame(), coords = c("x", "y")) %>% st_set_crs(4326) %>% st_transform(st_crs(map)) %>% st_coordinates()
  
  points(pts, cex = cex, pch = pch, col  = "grey10",
         bg = c("white", rep(adjustcolor("orange", alpha.f = 0.55), nrow(site01)-1)), lwd = 1.4)
  points(pr[1,1], pr[1,2], pch = 21, bg = adjustcolor("white", alpha.f = 0.3))
  par(opar)
  
  cat(sum(dead<1)/length(dead))
  
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

