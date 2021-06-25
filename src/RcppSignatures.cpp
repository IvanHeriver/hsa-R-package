#include <Rcpp.h>
// #include <RcppArmadillo.h>
using namespace Rcpp;


//------------------------------------------------------------------------------
//' compute a baseflow time series using the Lyne & Hollick algorithm.
//' 
//' This function compute baseflow time series using the algorithm of Lyne and Hollick
//' (1979)
//' 
//' @param Q numeric vector. Streamflow vector.
//' @param k float. parameter.
//' @return Returns a baseflow time series.
//' @references 
//' V. Lyne and M. Hollick,
//' “Stochastic time-variable rainfall-runoff modelling,”
//' in Hydrology and Water Ressources Symposium, Perth, Australia, 1979, pp. 89–92.
//' @export
// [[Rcpp::export]]
NumericVector baseflow_LyneHollick(NumericVector Q, double k) {
  int n = Q.size();
  NumericVector Qqf(n);
  Qqf[0] = 0;
  for (int i = 1; i < n; ++i) {
    if (NumericVector::is_na(Q[i])) {
      Qqf[i] = NA_REAL;
    } else {
      if (NumericVector::is_na(Q[i - 1])) {
        Qqf[i] = 0;
      } else {
        Qqf[i] = k * Qqf[i-1] + (1 + k) / 2 * (Q[i] - Q[i-1]);
        if (Qqf[i] < 0){
          Qqf[i] = 0 ;
        }
      }
    }
  }
  return Q - Qqf;
}

//------------------------------------------------------------------------------
//' compute a baseflow time series using the Chapman & Maxwell algorithm.
//' 
//' This function compute baseflow time series using the algorithm of Chapman and Maxwell
//' (1996)
//' 
//' @param Q numeric vector. Streamflow vector.
//' @param k float. parameter.
//' @return Returns a baseflow time series.
//' @references 
//' T. G. Chapman and A. I. Maxwell, 
//' “Baseflow separation - comparison of numerical methods with tracer experiments,” Hobart, Tasmania, 1996.
//' @export
// [[Rcpp::export]]
NumericVector baseflow_ChapmanMaxwell(NumericVector Q, double k) {
  int n = Q.size();
  NumericVector Qbf(n);
  Qbf[0] = Q[0];
  for (int i = 1; i < n; ++i) {
    if (NumericVector::is_na(Q[i])) {
      Qbf[i] = NA_REAL;
    } else {
      if (NumericVector::is_na(Q[i - 1])) {
        Qbf[i] = Q[i];
      } else {
        Qbf[i] = k / (2 - k) * Qbf[i - 1] + (1 - k) / (2 - k) * Q[i];
        if (Qbf[i] > Q[i]){
          Qbf[i] = Q[i] ;
        }
      }
    }
  }
  return Qbf;
}

//------------------------------------------------------------------------------
//' compute a baseflow time series using the Boughton algorithm.
//' 
//' This function compute baseflow time series using the algorithm of Boughton
//' described in Chapman (1999)
//' 
//' @param Q numeric vector. Streamflow vector.
//' @param k float. parameter.
//' @param C float. parameter.
//' @return Returns a baseflow time series.
//' @references 
//' T. Chapman,
//' “A comparison of algorithms for streamflow recession and basefow separation,”
//' Hydrological Processes, vol. 13, pp. 701–714, 1999.
//' @export
// [[Rcpp::export]]
NumericVector baseflow_Boughton(NumericVector Q, double k, double C) {
  int n = Q.size();
  NumericVector Qbf(n);
  Qbf[0] = Q[0];
  for (int i = 1; i < n; ++i) {
    if (NumericVector::is_na(Q[i])) {
      Qbf[i] = NA_REAL;
    } else {
      if (NumericVector::is_na(Q[i - 1])) {
        Qbf[i] = Q[i];
      } else {
        Qbf[i] = k / (1 + C) * Qbf[i - 1] + C / (1 + C) * Q[i];
        if (Qbf[i] > Q[i]){
          Qbf[i] = Q[i];
        }
      }
    }
  }
  return Qbf;
}


//------------------------------------------------------------------------------
//' compute a baseflow time series using the Eckhardt algorithm.
//' 
//' This function compute baseflow time series using the algorithm of Eckhardt (2005, 2008)
//' 
//' @param Q numeric vector. Streamflow vector.
//' @param a float. parameter.
//' @param BFImax float. parameter.
//' @return Returns a baseflow time series.
//' @references 
//'b K. Eckhardt,
//' “How to construct recursive digital filters for baseflow separation,” 
//' Hydrological Processes, vol. 19, pp. 507–515, 2005.
//' K. Eckhardt, 
//' “A comparison of baseflow indices, which were calculated with seven different baseflow separation methods,”
//' Journal of Hydrology, vol. 352, pp. 168–173, 2008, doi: 10.1016/j.jhydrol.2008.01.005.
//' @export
// [[Rcpp::export]]
NumericVector baseflow_Eckhardt(NumericVector Q, double a, double BFImax) {
  int n = Q.size();
  NumericVector Qbf(n);
  Qbf[0] = Q[0];
  for (int i = 1; i < n; ++i) {
    if (NumericVector::is_na(Q[i])) {
      Qbf[i] = NA_REAL;
    } else {
      if (NumericVector::is_na(Q[i - 1])) {
        Qbf[i] = Q[i];
      } else {
        Qbf[i] = ((1 - BFImax) * a * Qbf[i - 1] + (1 - a) * BFImax * Q[i]) / (1 - a * BFImax);
        if (Qbf[i] > Q[i]){
          Qbf[i] = Q[i];
        }
      }
    }
  }
  return Qbf;
}

// [[Rcpp::export]]
IntegerVector ioh_min_pivots(NumericVector Q, int d, double k) {
  int n = Q.size();
  
  // get non-overlapping d-day windows minima
  std::vector<int> imin;
  IntegerVector j = seq_len(d);
  NumericVector tmp(d);
  for (int i=0; (i + d - 1) < n; i += d) {
    tmp = Q[j + i - 1];
    imin.push_back(which_min(tmp) + i);
  }

  // look for pivot points in minima (i.e. when both previous and next values is above)
  int nmin = imin.size();
  LogicalVector valid_min = rep(false, nmin);
  int curr_i = 0, curr_ib = 0, curr_ia = 0;
  for (int i = 1; i < (nmin - 1); ++i) {
    curr_i = imin[i];
    curr_ib = imin[i - 1];
    curr_ia = imin[i + 1];
    valid_min[i] = (k * Q[curr_i] < Q[curr_ib]) & (k * Q[curr_i] < Q[curr_ia]);
  }
  
  // return valid pivot point
  IntegerVector out(nmin);
  out = imin;
  return out[valid_min];
}

// [[Rcpp::export]]
IntegerVector rec_events(IntegerVector imax, IntegerVector imin) {
  int n = imax.size(), m = imin.size();
  int j = 0, i = 0, i_done=0;
  IntegerVector new_imin(n);
  std::fill(new_imin.begin(), new_imin.end(), IntegerVector::get_na());
  for (int k = 0; k < (n - 1); ++k) { // loop over all index in imax (but the last one)
    j = i_done; // take the first index to consider to look in imin (previously saved)
    i = imin[j]; 
    while (i < imax[k + 1]) { // as long as investigated imin is below the next imax
      if (i > imax[k]) { // if investigated imin is above the current imax ==> job's done
        new_imin[k] = i;
        i_done = j + 1;
        break;
      } // else: we keep looking
      j += 1;
      if (j >= m) { // except if the end of imin is reached
        break;
      }
      i = imin[j]; // next investigated imin
    }
    if (j >= m) { // no point continuing if end of imin is reached
      break;
    }
  }
  // now the special case of the last value of imax
  if (j < m) { // only if imin wasn't completly searched
    j = i_done;
    i = imin[j];
    while (j < m) {
      if (i > imax[n - 1]) {
        new_imin[n - 1] = i;
        break;
      }
      j += 1;
      if (j >= m) { // except if the end of imin is reached
        break;
      }
      i = imin[j];
    }
  }
  return new_imin;
}


/****R

*/

