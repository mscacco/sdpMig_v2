// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include <iostream>

// Namespace for Parameters
namespace sdp {

// Init
int MaxT;
int NSites;
int MaxX;

// Individual
double B0;
double w;
double  xc;
arma::vec b0;
arma::vec b1;
arma::vec b2;
double pred_a1;
double pred_a2;
double c;
double speed;
double max_u;
double f;
double flightMort; // mortality risk per day of flight
arma::vec WindAssist;
arma::vec WindProb;
arma::vec ZStdNorm;
arma::vec PStdNorm;
double decError;

// Sites  
arma::mat dist;
arma::mat x_gain;
arma::mat y_gain;
arma::mat p_gain;
arma::vec nTR_x;
arma::vec nTR_y;
arma::mat expend;

// Output
arma::cube FMatrix;
arma::cube DMatrix1;
arma::cube DMatrix2;
arma::cube PMatrix1;
arma::cube PMatrix2;

// Internal paramters
double Fintensity = 0.0;

}

// Eval (not exported):
// Function use to truncate value to range between min and max
int Eval(
    double x, 
    double min, 
    double max
) 
{
  return (x <= min) ? 0 : ((x >= (max)) ? 2 : 1);
}


// Decision Error (not exported)
// Function to calculate decision error based on the importance of the decision
double DecError(
    double var,
    double decError
)
{
  double help1 = 1.0 + (decError * var);
  double help  = 1.0 / help1;
  return help;
}

// Round (not exported)
// Function equivalent to round(value, 0) in R
int Round(double value)
{  
  // Round half up, correct also for negative values ((int) truncates towards 0, so the old version gave e.g. -0.3 -> 1, -1.7 -> -1)
  return (int)(floor(value + 0.5));
}

// Sigmoidal TR function (not exported)
inline double FNtanh(double g)
{
  return ((( (exp(g) - exp(-g)) / (exp(g) + exp(-g)) ) + 1.0)/ 2.0);
}

// InterpolateTR
// Function to derive Terminal Reward for t
double InterpolateTR (int c_tr,
                      arma::vec  c_nTR_x,
                      arma::vec  c_nTR_y) 
{
  int i = 0;
  double help = 0.0;
  double help_t = (double)(c_tr);
  double min = 0.0; double max = 2.0;
  
  while((c_nTR_x(i) < help_t) && (i < (c_nTR_x.size()))) ++i;
  if((i==0) || (c_nTR_x(i) - c_nTR_x[i-1]) == 0)
    help = c_nTR_y(i);
  else
    help = c_nTR_y(i-1) + 
      (c_nTR_y(i) - c_nTR_y(i-1))/
        (c_nTR_x(i) - c_nTR_x(i-1)) * (help_t - c_nTR_x(i-1));
  
  if (help < min) help = min;
  if (help > max) help = max;
  return help;
}

// CalculateTR
// Function to calculate the Terminal Reward
double CalculateTR(int c_x, 
                   int c_t, 
                   double     c_w, 
                   double     c_xc, 
                   const int  c_B0,
                   arma::vec  c_nTR_x,
                   arma::vec  c_nTR_y) 
{
  double TimeReward, StateReward, TR_tmp = 0.0;
  
  TimeReward = InterpolateTR(c_t, c_nTR_x, c_nTR_y);
  if (c_x == 0) StateReward = 0.0;
  else StateReward = FNtanh(c_w * ((double)(c_x) - c_xc));
  
  TR_tmp = TimeReward * StateReward + c_B0;
  return TR_tmp;
}


////////////////////////////////////////////////////////////////
////// Functions used in Backward and Forward Simulation ///////
////////////////////////////////////////////////////////////////

// Interpolate
// Get gain values for time and site
// Interpolation required as originally gain function was provided as a list of
// (time, gain) values that needed to be interpolated

double Interpolate(
    const int& time,
    const int& s,
    double min, 
    double max,
    bool p
)
{
  int i = 1;
  double help;
  double nGains = sdp::x_gain.n_cols;
  
  if (!p) 
  {
    while((sdp::x_gain(s, i) < time) && (i < (nGains-1))) ++i;
    if ((i == 1) || ((sdp::x_gain(s, i) - sdp::x_gain(s, i-1)) == 0))
      help = sdp::y_gain(s, i);
    else
      help = sdp::y_gain(s, i-1) + (sdp::y_gain(s, i) - sdp::y_gain(s, i-1))/
        (sdp::x_gain(s, i) - sdp::x_gain(s, i-1))*(time - sdp::x_gain(s, i-1));
    
    if (help < min) help = min;
    if (help > max) help = max;
    
    return help;
  } 
  else
  {
    while((sdp::x_gain(s, i) < time) && (i < (nGains-1))) ++i;
    if ((i == 1) || ((sdp::x_gain(s, i) - sdp::x_gain(s, i-1)) == 0))
      help = sdp::p_gain(s, i);
    else
      help =  sdp::p_gain(s, i-1) + (sdp::p_gain(s, i) - sdp::p_gain(s, i-1))/
        (sdp::x_gain(s, i) - sdp::x_gain(s, i-1))*(time - sdp::x_gain(s, i-1));
    
    if (help < min) help = min;
    if (help > max) help = max;
    
    return help;
  }
} // end Interpolate


// --------------------      Predation.    ------------------------------------
// Calculates the predation risk based on background predation risk, foraging-
// intensity-dependent predation risk and body-reserves-dependent predation risk

double Predation (
    const int& time, 
    const int& site, 
    const int& x,
    double u,
    double forage
)
{
  const double rew_tol = 0.0000005;
  double netgain, cor_u, cor_x;
  
  cor_u = u;
  cor_x = (double)(x);
  netgain = ((u * forage) - sdp::expend(site, time));

  if (netgain <= 0.0) netgain = rew_tol;
  if (cor_u   <= 0.0) cor_u   = rew_tol;
  if (cor_x   <= 0.0) cor_x   = rew_tol;
  
  double help1 = exp((sdp::pred_a1 + 1) * log(cor_x + netgain));
  double help2 = exp((sdp::pred_a1 + 1) * log(cor_x));
  double help3 = (sdp::pred_a1 + 1.0) * netgain;
  double help4 = exp(sdp::pred_a2 * log(cor_u));
  double help5 = sdp::b0(site) + sdp::b1(site) * ((help1 - help2)/help3) * sdp::b2(site) * help4;
  
  // double help1 = sdp::b1(site) * (pow(cor_u, sdp::pred_a1));
  // double help2 = sdp::b2(site) * (pow(cor_x/((double)(sdp::MaxX)), sdp::pred_a2));
  // double help5 = sdp::b0(site) + help1 + help2;
  
  return help5;
} // end Predation


// FindFitness value
// Interpolates FitnessValue for time, site, x and feeding intensity
double FindF (
    const int& time, 
    const int& site, 
    const int& x,
    int accuracy, 
    double gain, 
    double fx, 
    double u
)
{
  double res1, res2, part1, part2, interpolReward;
  res1 = 0.0; res2 = 0.0; part1 = 0.0; part2 = 0.0;
  
  double expenditure = sdp::expend(site, time);
  double nextx = (double)(x) + gain * u - expenditure;
  switch (Eval(nextx, 0.0, (double)(sdp::MaxX)))
  {
  case 0: interpolReward = sdp::FMatrix(time+1, site, 0); break;
  case 1: {    res1 = nextx - (int)(nextx);
    res2 = (int)(nextx) + 1.0 - nextx;
    
    part1 = res1 * sdp::FMatrix(time+1, site, (int)(nextx)+1);
    part2 = res2 * sdp::FMatrix(time+1, site, (int)(nextx));
    interpolReward = part1 + part2;
    break;
  }    
  case 2: interpolReward = sdp::FMatrix(time+1, site, sdp::MaxX); break;
  }
  double fx_new = fx + (sdp::PStdNorm(accuracy) * interpolReward);
  return fx_new;
} // end FindF


// -----------------         Foraging.       -----------------------------------
// Major Foraging function: 
// Calculates the optimal foraging intensity, given current body condition, expected 
// new body condition and the (predation and starvation) risk involved 

double Foraging(const int& time, 
                const int& site, 
                const int& x
)
{
  const double r = 0.61803399;
  const double tol = 0.0000005;
  double mean, gain, SD, Freward;
  double c, u0, u1, u2, u3, f1, f2, hold1, hold2;
  double hold1_old, hold2_old, f1_old, f2_old;
  int NStdNorm = sdp::ZStdNorm.size();
  
  c = 1.0 - r;
  u0 = 0.0;
  u1 = r;
  u3 = 1.0;
  u2 = u1 + c * (u3 - u1);
  
  f1 = 0.0; f2 = 0.0;
  
  mean = sdp::y_gain(site,time);
  SD   = sdp::p_gain(time);
  
  for (int accuracy = 0; accuracy < NStdNorm; ++accuracy)
  {
    gain = mean + sdp::ZStdNorm(accuracy) * SD;
    f1_old = f1;
    f2_old = f2;
    f1 = FindF(time, site, x, accuracy, gain, f1_old, u1);
    f2 = FindF(time, site, x, accuracy, gain, f2_old, u2);
  }
  
  f1_old = f1;
  f2_old = f2;
  // Predation() takes the food gain per day (mean), not the fitness of the option (f1, f2, hold1, hold2 were passed before)
  f1 = (1.0 - Predation(time, site, x, u1, mean)) * f1;
  f2 = (1.0 - Predation(time, site, x, u2, mean)) * f2;
  
  while ((fabs(u3 - u0)) > tol)
  {
    if (f2 > f1)
    {
      u0 = u1;
      u1 = u2;
      u2 = (r * u1) + (c * u3);
      hold2 = 0.0;
      for (int accuracy = 0; accuracy < NStdNorm; ++accuracy)
      {
        gain = mean + sdp::ZStdNorm(accuracy) * SD;
        hold2_old = hold2;
        hold2 = FindF(time, site, x, accuracy, gain, hold2_old, u2);
      }
      f1 = f2;
      f2 = ((1.0 - Predation(time, site, x, u2, mean)) * hold2);
    }
    if (f2 <= f1)
    {
      u3 = u2;
      u2 = u1;
      u1 = (r * u2) + (c * u0);
      hold1 = 0.0;
      for (int accuracy = 0; accuracy < NStdNorm; ++accuracy)
      {
        gain = mean + sdp::ZStdNorm(accuracy) * SD;
        hold1_old = hold1;
        hold1 = FindF(time, site, x, accuracy, gain, hold1_old, u1);
      }
      f2 = f1;
      f1 = ((1.0 - Predation(time, site, x, u1, mean)) * hold1);
    }
  }  //end of while loop
  
  if (f1 <= f2)
  {
    sdp::Fintensity = u2;
    Freward = f2;
  }
  else //(f1 > f2)
  {
    sdp::Fintensity = u1;
    Freward = f1;
  }
  return Freward;
} // end Foraging --------------------------------------------------------------


// -----------------------   Flying    ----------------------------------------

// Major Flying Function - 
// calculates consequences of flying from current site "site" to a destination 
// site "dep_site", based on body resources required to cover the distance and
// time required to reach destination

// [[Rcpp::export]]

double Flying(
    const int& time, 
    const int& site, 
    const int& x, 
    int dep_site)
{
  double totalD, Whold, range, Sqr_c, Sqr_ca;
  double distance, t, f, interpolReward, nextx;
  double res1, res2, tim1, tim2, part1, part2, part3, part4;
  
  interpolReward = 0.0;
  totalD = 0.0;
  Whold = 0.0;
  f = sdp::f; //1.0;
  
  // Flight reserve: (1 - f) * MaxX is body storage that can never be used for flight; only the storage above it is fuel.
  // A bird with x <= xReserve cannot fly. After the flight the reserve is still there (nextx = xReserve + remaining fuel).
  // Same rule as in simFlying().
  double xReserve = (1.0 - f) * sdp::MaxX;
  double usableF  = std::max(0.0, (double)(x) - xReserve) / sdp::MaxX;
  range = sdp::c * (1.0 - 1.0 / sqrt(1.0 + usableF));

  totalD = sdp::dist(site, dep_site); // distance between starting site and destination
  for (int h = 0; h < sdp::WindProb.size(); ++h) // wind support
  {
    distance = totalD * (1.0 + sdp::WindAssist(h));

    Sqr_c = (sdp::c * sdp::c);
    Sqr_ca = (sdp::c - (range - distance))*(sdp::c - (range - distance));
    nextx = xReserve + ((Sqr_c/Sqr_ca) - 1.0) * sdp::MaxX; // reserve + remaining usable fuel
    // Flights last whole days, rounded up (e.g. 0.3 -> 1 day, 1.3 -> 2 days), so a flight never reads the time slice
    // that is being computed (same-day arrivals read sites not yet updated at this time step). Same rule as in simFlying().
    t = (double)(time) + std::ceil(distance / sdp::speed);
    // The part of the last flight day not spent flying is spent on the ground: it costs the daily ground expenditure
    // (proportional to the time left), otherwise a short flight is cheaper than a day on the ground. Same rule as in simFlying().
    nextx -= sdp::expend(dep_site, time) * (std::ceil(distance / sdp::speed) - distance / sdp::speed);

    if (t >= (sdp::MaxT)) interpolReward = 0.0;
    else  if (distance > range) interpolReward = 0.0; // destination out of range of the usable fuel
    else  if (nextx <= 0.0) interpolReward = 0.0;
    else
    {
      res1 = nextx - (int)(nextx);
      res2 = (int)(nextx) + 1.0 - nextx;
      tim1 = t - (int)(t);
      tim2 = (int)(t) + 1.0 - t;
      
      if ((nextx+1.0) > (double)(sdp::MaxX))
      {
        part1 = tim1*res1*sdp::FMatrix((int)(t)+1, dep_site, sdp::MaxX);
        part3 = tim2*res1*sdp::FMatrix((int)(t),   dep_site, sdp::MaxX);
      }
      else
      {
        part1 = tim1*res1*sdp::FMatrix((int)(t)+1, dep_site, (int)(nextx)+1);
        part3 = tim2*res1*sdp::FMatrix((int)(t),   dep_site, (int)(nextx)+1);
      }
      part2 = tim1*res2*sdp::FMatrix((int)(t)+1, dep_site, (int)(nextx));
      part4 = tim2*res2*sdp::FMatrix((int)(t),   dep_site, (int)(nextx));
     
      interpolReward = part1+part2+part3+part4;
      interpolReward *= (1.0 - sdp::flightMort * distance / sdp::speed); // survival during the flight (risk per day of flight time)
    }
    double help_old = Whold;
    Whold = help_old + (sdp::WindProb(h) * interpolReward);
  } //End h
  return Whold;
}  // end Flying ---------------------------------------------------------------


////////////////////////////////////////////////////////////////
////// Export Functions: Backward Iteration ////////////////////
////////////////////////////////////////////////////////////////


// [[Rcpp::export]]
void Init( int MaxT,
           int NSites,
           int MaxX,
           double w,
           double xc,
           int B0,
           Rcpp::NumericVector b0,
           Rcpp::NumericVector b1,
           Rcpp::NumericVector b2,
           double pred_a1,
           double pred_a2,
           double c,
           double speed,
           double f,
           double flightMort,
           Rcpp::NumericVector WindAssist,
           Rcpp::NumericVector WindProb,
           Rcpp::NumericVector ZStdNorm,
           Rcpp::NumericVector PStdNorm,
           Rcpp::NumericVector nTR_x,
           Rcpp::NumericVector nTR_y,
           double decError,
           arma::mat dist,
           Rcpp::NumericVector x_gain,
           arma::mat y_gain,
           Rcpp::NumericVector p_gain,
           arma::mat expend
           )
{
  // Basic Parms
  sdp::MaxT   = MaxT; //Maximum time the model runs
  sdp::NSites = NSites; // Maximum number of sites
  sdp::MaxX   = MaxX;  // Maximum level of body reserves
  
  // Individual
  sdp::B0 = B0;
  sdp::w = w;
  sdp::xc = xc;
  sdp::b0 = b0;
  sdp::b1 = b1;
  sdp::b2 = b2;
  sdp::pred_a1 = pred_a1;
  sdp::pred_a2 = pred_a2;
  sdp::c = c;
  sdp::speed = speed;
  sdp::f = f;
  sdp::flightMort = flightMort;
  sdp::WindAssist = WindAssist;
  sdp::WindProb = WindProb;
  sdp::ZStdNorm = ZStdNorm;
  sdp::PStdNorm = PStdNorm;
  sdp::nTR_x = nTR_x;
  sdp::nTR_y = nTR_y;
  sdp::decError = decError;
  
  // Sites  
  sdp::dist   = dist;
  sdp::x_gain = x_gain;
  sdp::y_gain = y_gain;
  sdp::p_gain = p_gain;
  sdp::expend = expend;
 
  sdp::Fintensity = 0.0; 
}


// [[Rcpp::export]]
Rcpp::List BackwardIteration(bool pbar) {
  
  /// Terminal Reward
  arma::cube FM_TR = arma::zeros<arma::cube>(sdp::MaxT+1, sdp::NSites, sdp::MaxX+1);
  
  for (int r = 0; r <= sdp::MaxT; ++r) //loop through all time-steps
  {
    for (int c = 0; c <= sdp::MaxX; ++c) //loop through all body condition levels
    {
      // Calculate terminal reward
      FM_TR(r, (sdp::NSites-1), c) = CalculateTR(c, r, sdp::w, sdp::xc, sdp::B0, sdp::nTR_x, sdp::nTR_y);
    }
  }
  
  for (int r = 0; r < (sdp::NSites-1); ++r)
  {
    for (int c = 0; c <= sdp::MaxX; ++c)
    {
      FM_TR(sdp::MaxT, r, c) = sdp::B0;
    }
  }
  
  sdp::FMatrix  = FM_TR;
  sdp::DMatrix1 = arma::zeros<arma::cube>(sdp::NSites, sdp::MaxT+1, sdp::MaxX+1);
  sdp::DMatrix2 = arma::zeros<arma::cube>(sdp::NSites, sdp::MaxT+1, sdp::MaxX+1);
  sdp::PMatrix1 = arma::zeros<arma::cube>(sdp::NSites, sdp::MaxT+1, sdp::MaxX+1);
  sdp::PMatrix2 = arma::zeros<arma::cube>(sdp::NSites, sdp::MaxT+1, sdp::MaxX+1);
  

  arma::vec  Mrew(sdp::NSites+1);
  double     max_Reward, decision;
  
  Rcpp::List out(5);

  // Progress bar
  // Progress p(sdp::MaxT*(sdp::NSites), pbar);

  for (int time = (sdp::MaxT-1); time >= 0; --time)
  {
    // The breeding site (NSites-1) is excluded: it is absorbing, its value is the terminal reward set in FM_TR above.
    // Including it overwrote the terminal reward with the value of foraging/flying on from the breeding site.
    for (int site = 0; site < (sdp::NSites - 1); ++site)
    {
      sdp::FMatrix(time, site, 0) = 0.0;
      // p.increment(); // update progress

      for(int x = 1; x <= sdp::MaxX; ++x) {

        sdp::Fintensity = 0.0;
        max_Reward = 0.0;
        decision = 0.0;
        
        double f_help = Foraging(time, site, x); // calculate reward from foraging
        sdp::FMatrix(time, site, x) = f_help;
        
       
        decision   = - sdp::Fintensity - 1.0; // optimal foraging intensity
        sdp::DMatrix1(site, time, x) = decision;

//        max_Reward = Foraging(time, site, x); // why call Foraging twice??
        max_Reward = f_help;                  // shouldn't it be just like that?

                
// Caclulate reward from flying to another site
        for (int bb = 0; bb < sdp::NSites; ++bb) Mrew(bb) = 0.0;
        for (int dest = 0; dest < sdp::NSites; ++dest) // loops thru all sites - why not exclude current?
        {
          if (dest == site) Mrew[dest] = 0.0;
          else {
          float help_z = Flying(time, site, x, dest);
          Mrew[dest] = help_z;
          }
        }

        double help_flying = 0.0;
        int best_bb = 0;
        for (int bb = 0; bb < sdp::NSites; ++bb)
        {
          if (Mrew[bb] > help_flying)
          {
            help_flying = Mrew(bb);
            best_bb = bb;
            sdp::DMatrix2(site, time, x) = (double)(bb);
            if (help_flying > max_Reward) max_Reward = help_flying;
          }
        }

        if (help_flying <= 0.0) {
          // No destination can be reached: the bird forages. Before, the flight decision stayed at its default
          // (site 0, the wintering site) and was taken with a small probability through decision error.
          sdp::DMatrix2(site, time, x) = arma::datum::nan;
          sdp::PMatrix1(site, time, x) = 1.0;
          sdp::PMatrix2(site, time, x) = 0.0;
        } else {
        sdp::PMatrix1(site, time, x) = DecError((max_Reward - sdp::FMatrix(time, site, x))/max_Reward, sdp::decError);
        sdp::PMatrix2(site, time, x) = DecError((max_Reward - Mrew(best_bb))/max_Reward, sdp::decError);

        double sum_action = sdp::PMatrix1(site, time, x) + sdp::PMatrix2(site, time, x);

        if (sum_action > 0.0) {
        sdp::PMatrix1(site, time, x) /= sum_action;
        sdp::PMatrix2(site, time, x) /= sum_action;
        }
        else {
          sdp::PMatrix1(site, time, x) = 0.0;
          sdp::PMatrix2(site, time, x) = 0.0;
        }
        } // end else (a destination can be reached)
        sdp::FMatrix(time, site, x)   = max_Reward;
        
      } // end x
    } // end site
  } // end time

  
  out[0] = sdp::FMatrix;
  out[1] = sdp::DMatrix1;
  out[2] = sdp::DMatrix2;
  out[3] = sdp::PMatrix1;
  out[4] = sdp::PMatrix2;
  
  return out;
}


////////////////////////////////////////////////////////////////
////// Export Functions: Forward Simulation ////////////////////
////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
void InitSim (int MaxT,
              int NSites,
              int MaxX,
              double w,
              double xc,
              int B0,
              Rcpp::NumericVector b0,
              Rcpp::NumericVector b1,
              Rcpp::NumericVector b2,
              double pred_a1,
              double pred_a2,
              double c,
              double speed,
              double f,
              double flightMort,
              Rcpp::NumericVector WindAssist,
              Rcpp::NumericVector WindProb,
              Rcpp::NumericVector ZStdNorm,
              Rcpp::NumericVector PStdNorm,
              Rcpp::NumericVector nTR_x,
              Rcpp::NumericVector nTR_y,
              double decError,
              arma::mat dist,
              Rcpp::NumericVector x_gain,
              arma::mat y_gain,
              Rcpp::NumericVector p_gain,
              arma::mat expend)
{
  // Basic Parms
  sdp::MaxT   = MaxT;
  sdp::NSites = NSites;
  sdp::MaxX   = MaxX;
  
  // Individual
  sdp::B0 = B0;
  sdp::w = w;
  sdp::xc = xc;
  sdp::b0 = b0;
  sdp::b1 = b1;
  sdp::b2 = b2;
  sdp::pred_a1 = pred_a1;
  sdp::pred_a2 = pred_a2;
  sdp::c = c;
  sdp::speed = speed;
  sdp::f = f;
  sdp::flightMort = flightMort;
  sdp::WindAssist = WindAssist;
  sdp::WindProb = WindProb;
  sdp::ZStdNorm = ZStdNorm;
  sdp::PStdNorm = PStdNorm;
  sdp::nTR_x = nTR_x;
  sdp::nTR_y = nTR_y;
  sdp::decError = decError;
  
  // Sites  
  sdp::dist   = dist;
  sdp::x_gain = x_gain;
  sdp::y_gain = y_gain;
  sdp::p_gain = p_gain;
  sdp::expend = expend;

}


// [[Rcpp::export]]
arma::vec simForaging(double f_intensity, int time, int site, int x)
{
  double gain, new_x = 0.00, mean, SD, Predat, Runif, pre_x;
  double LB, UB;
  int hit = 0, dead = false;
  int NStdNorm = sdp::ZStdNorm.size();
  arma::vec out = arma::zeros<arma::vec>(2);
 
  Runif = R::runif(0,1);
  LB = 0.0;
  UB = 0.0;
  for (int h = 0; h < NStdNorm; ++h)
  {
    UB += sdp::PStdNorm(h);
    if ((Runif > LB) && (Runif <= UB))
      hit = h;
    LB += sdp::PStdNorm(h);
  }
  
  mean = sdp::y_gain(site, time);
  SD =   sdp::p_gain(time);
  
  double temp_gain = mean + sdp::ZStdNorm(hit) * SD;
  gain = temp_gain;
  if (gain < 0.0) gain = 0.0;
  
  pre_x = double(x);
  double expenditure = sdp::expend(site, time);
  
  new_x = pre_x + f_intensity * gain - expenditure;
  
  
  switch(Eval(new_x, 0.0, sdp::MaxX))
  {
  case 0: x = 0; break;
  case 1: x = Round(new_x); break;
  case 2: x = sdp::MaxX; break;
  };
  
  // Same net gain and same rule for netgain <= 0 as in the backward iteration (Predation())
  double netgain = f_intensity * gain - expenditure;
  if (netgain <= 0.0) netgain = 0.0000005;
  if ((f_intensity > 0.0) && (pre_x > 0.0))
  {
    double help1 = exp((sdp::pred_a1 + 1.0) * log(pre_x + netgain));
    double help2 = exp((sdp::pred_a1 + 1.0) * log(pre_x));
    double help3 = (sdp::pred_a1 + 1.0) * netgain;
    double help4 = exp(sdp::pred_a2 * log(f_intensity));
    Predat = sdp::b0(site) + sdp::b1(site) *
      ((help1 - help2)/help3) * sdp::b2(site) * help4;
  }
  else Predat = sdp::b0(site);
  
  if (Predat > 0.0)
  {
    if (R::runif(0, 1) < Predat) dead = true;
  }

  
  out(0) = Round(new_x);
  out(1) = dead;
  return out;
}

// [[Rcpp::export]]
arma::vec simFlying(int decision, int time, int site, int x)
{
  double nextx, total_D, distance, range;
  double Sqr_c, Sqr_ca, t, f, UB, LB, Runif;
  int hit = 0;
  int NWind = sdp::WindAssist.size();
  arma::vec out = arma::zeros<arma::vec>(2);
  f = sdp::f;
  
  int dest_site = decision;
  total_D = 0.0;
  total_D = sdp::dist(site, dest_site);

  Runif = R::runif(0,1);
  LB = 0.0;
  UB = 0.0;
  for (int h = 0; h < NWind; ++h)
  {
    UB += sdp::WindProb(h);
    if ((Runif > LB ) && (Runif <= UB))
      hit = h;
    LB += sdp::WindProb(h);
  }
  distance = total_D * (1.0 + sdp::WindAssist(hit));
  
  // Flight reserve, same rule as in Flying(): only the storage above (1 - f) * MaxX is fuel, the reserve is kept after the flight.
  double xReserve = (1.0 - f) * sdp::MaxX;
  double usableF  = std::max(0.0, (double)(x) - xReserve) / (double)(sdp::MaxX);
  range = sdp::c * (1.0 - 1.0 / sqrt(1.0 + usableF));

  Sqr_c  = (sdp::c * sdp::c);
  Sqr_ca = (sdp::c - (range - distance))*(sdp::c - (range - distance));
  nextx  = xReserve + (((Sqr_c/Sqr_ca) - 1.0) * (double)(sdp::MaxX)); // reserve + remaining usable fuel
  // Ground expenditure for the part of the last flight day not spent flying, same rule as in Flying()
  nextx -= sdp::expend(dest_site, time) * (std::ceil(distance / sdp::speed) - distance / sdp::speed);
  if (R::runif(0, 1) < sdp::flightMort * distance / sdp::speed) nextx = 0.0; // death during the flight (risk per day of flight time): x = 0, dead
  if (distance > range) nextx = 0.0; // destination out of range (only through decision error): x = 0, the R wrapper marks the bird as dead

  
  t =  time + std::ceil(distance / sdp::speed); // whole days, rounded up, same rule as in Flying()
  
  out(0) = Round(t);
  out(1) = Round(nextx);
  return out;
}