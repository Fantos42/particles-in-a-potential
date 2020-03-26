// [[Rcpp::plugins("cpp11")]]
#include <Rcpp.h>
#include <vector>
using namespace Rcpp;


// =======================================================
// ========================================= IAP_CPP
// =======================================================
// [[Rcpp::export]]
double IAT_CPP(NumericVector x){
  const int n = x.size(); //number of particles
  std::vector<double> X (n, 0);
  for (int i = 0; i < n; i++) X [i] = x(i);
  
  const int maxlag = std::max(3., n/2.);
  
  double mu = 0;
  for (int i = 0; i < n; i++) mu += X[i];
  mu /= n;
  
  double Ga0 = 0;
  double Ga1 = 0;
  double Ga2 = 0;
  double Ga3 = 0;
  
  double lg;
  
  lg = 0;
  for (int i = lg; i < n; i++) Ga0 += (X[i-lg] - mu) * (X[i] - mu) / (n-lg-1.);
  lg = 1;
  for (int i = lg; i < n; i++) Ga1 += (X[i-lg] - mu) * (X[i] - mu) / (n-lg-1.);
  lg = 2;
  for (int i = lg; i < n; i++) Ga2 += (X[i-lg] - mu) * (X[i] - mu) / (n-lg-1.);
  lg = 3;
  for (int i = lg; i < n; i++) Ga3 += (X[i-lg] - mu) * (X[i] - mu) / (n-lg-1.);
  
  double Ga [] = {Ga0+Ga1, Ga2+Ga3};
  
  double m = 1;
  double IAT = Ga[0];
  while ( (Ga[1]>0.0) && (Ga[1] < Ga[0] || Ga[1]>0.1) ) {
    m += 1;
    if (2*m+1 > maxlag) {
      Rcout << "Not enough data, maxlag= " << maxlag << std::endl;
      return maxlag;
    }
    Ga[0] = Ga[1];
    lg = 2*m;
    Ga[1] = 0;
    for (int i = lg; i < n; i++) Ga[1] += (X[i-lg] - mu) * (X[i] - mu) / (n-lg-1.);
    lg = 2*m + 1;
    for (int i = lg; i < n; i++) Ga[1] += (X[i-lg] - mu) * (X[i] - mu) / (n-lg-1.);
    IAT += Ga[0];
  }
  Rcout << "Interrupted IAT calculation at dT = " << m*2 << ". Ga[dT] = " << Ga[1] << std::endl;
  //Rcout << "Breakpoint: " << 2*m << std::endl;
  IAT = (-1. + 2. * IAT / Ga0) / 2.;
  return IAT;
}

// =======================================================
// ========================================= Calculate Energy
// =======================================================
// [[Rcpp::export]]
double E_CPP(NumericMatrix x, NumericVector q){
  int n = q.size(); //number of particles
  int d = x.ncol(); //dimensions
  
  double a = 0;
  double b = 0;
  double sum = 0;
  for(int i = 0; i < n-1; i++){
    for(int j = i+1; j < n; j++){
      sum = 0;
      for(int k = 0; k < d; k++){
        sum += (x(i,k) - x(j,k))*(x(i,k) - x(j,k));
      }
      a += q[i] * q[j] / sqrt(sum);
      b += 1 / pow(sum, 4); 
    }
  }
  return (a+b);
}

// =======================================================
// ========================================= estimate_clstSize
// =======================================================

int getFreeIndex(std::vector<int> idx, int k) {
  const int nPart = idx.size();
  
  for (int i = 0; i < nPart; i++) {
    if (idx[i]==1) k--;
    
    if (k == 0) {
      return i;
    }
  }
  
  return -1;
}

double getDist(std::vector<double> a, std::vector<double> b) {
  const int dim = a.size();
  
  double dist = 0;
  for (int i = 0; i < dim; i++) {
    dist += (a[i] - b[i])*(a[i] - b[i]);
  }
  return std::sqrt(dist);
}

std::vector<int> estimate_clstSize(std::vector<std::vector<double>> X) {
  const int nPart = X.size();
  
  const double min_dist = 2;
  
  std::vector<std::vector<int>> clusters;
  std::vector<int>              cluster_sizes;
  
  std::vector<int> idx(nPart, +1);
  int nFree = nPart;
  
  while (nFree > 0) 
  {
    // Choose first free particle for new cluster
    int k = -1;
    for (int i = 0; i < nPart; i++) 
    {
      if (idx[i] == 1) 
      {
        k = i;
        idx[i] = 0;
        nFree--;
        break;
      }
    }
    
    // Create new cluster vector
    std::vector<int> temp (nPart, -1);
    clusters.push_back(temp);
    clusters.back()[0] = k;
    cluster_sizes.push_back(1);
    
    
    // Iterate through all cluster particles and look for neighbours.
    for (int j = 0; j < cluster_sizes.back(); j++) 
    {
      const int clst_idx =clusters.back()[j];
      // Add all neighbours of particle j.
      for (int k = 0; k < nFree; k++) 
      {
        int free_idx = getFreeIndex(idx, k+1);
        
        if (getDist(X[clst_idx], X[free_idx]) < min_dist) {
          clusters.back()[cluster_sizes.back()]  = free_idx;
          cluster_sizes.back()  += 1; // Current cluster has grown by 1
          idx[free_idx]          = 0; // Particle has been assigned to a cluster
          nFree                 -= 1; // One further particle is assigned to a cluster
          //k                     -= 1; // Neccessary for correct procedure
        }
        
      }
    }
  }
  
  return cluster_sizes;
}

// =======================================================
// ========================================= other_observables
// =======================================================

std::vector<double> estimate_particle_dists(std::vector<std::vector<double>> X){
  const unsigned int nPart = X.size();
  const unsigned int d = X[0].size();
  // double dist = 0;
  // const double normalization = 2/nPart/(nPart-1);
  std::vector<double> particle_dist_vec;
  
  for(unsigned int i = 0; i < nPart; i++){
    for(unsigned int j = i+1; j < nPart; j++){
      double sum = 0;
      for(unsigned int dim = 0; dim < d; dim++){
        sum += (X[i][dim] - X[j][dim])*(X[i][dim] - X[j][dim]);
      }
      // dist += sqrt(sum);
      particle_dist_vec.push_back(sqrt(sum));
    }
  }
  
  return(particle_dist_vec);
}

std::vector<double> estimate_nearest_neighbour_dists(std::vector<std::vector<double>> X){
  const unsigned int nPart = X.size();
  double min_dist = 1e7;
  // double sum_min_dists = 0;
  std::vector<double> min_dist_vec;
  
  for(unsigned int i = 0; i < nPart; i++){
    for(unsigned int j = 0; j < nPart; j++){
      if(i == j) continue;
      if(getDist(X[i], X[j]) < min_dist) min_dist = getDist(X[i], X[j]);
    }
    min_dist_vec.push_back(min_dist);
    // sum_min_dists += min_dist;
  }
  
  return(min_dist_vec);
}

std::vector<int> estimate_Nparticles_neighbourhood(std::vector<std::vector<double>> X){
  const unsigned int nPart = X.size();
  const double min_dist = 2;
  std::vector<int> N_neigh; 
  
  for(unsigned int i = 0; i < nPart; i++){
    unsigned int neigh = 0;
    for(unsigned int j = 0; j < nPart; j++){
      if(i == j) continue;
      if(getDist(X[i], X[j]) < min_dist) neigh++; 
    }
    N_neigh.push_back(neigh);
  }
  
  return(N_neigh);
}

// =======================================================
// ========================================= MCMC_CPP
// =======================================================
// [[Rcpp::export]]
Rcpp::List MCMC_CPP(int nIt, int nPart, double vol, double t, double sigma, const double N_burnIn, NumericMatrix X0, NumericVector q){
  int d = X0(1,Rcpp::_).size();
  const double l = std::pow(vol, 1./d);
  const double t0= t;
  
  const double T_burnIn_1 = N_burnIn * 3./4.;
  const double T_burnIn_2 = N_burnIn;
  // ----------------------- Observables
  NumericVector dE(nIt, 0.0);
  dE[0] = E_CPP(X0, q);
  // -----------------------
  
  std::vector<std::vector<double>> X(nPart, std::vector<double>(d,0));
  for (int i = 0; i < nPart; i++) {
    for (int j = 0; j < d; j++) {
      X[i][j] = X0(i,j);
    }
  }
  double accept_rate = 0;
  int k  = -1;
  for (int n = 1; n < nIt; n++) {
    k = std::floor(R::runif(0,nPart));
    // ---------------------------------------------------------------------------
    for(int j = 0; j < nPart; j++){
      if(j == k) continue;
      
      double sum = 0;
      for(int dim = 0; dim < d; dim++){
        sum += (X[k][dim] - X[j][dim])*(X[k][dim] - X[j][dim]);
      }
      
      dE[n] -= q[k] * q[j] / sqrt(sum) + 1/ pow(sum, 4);
    }
    // ---------------------------------------------------------------------------
    // Make new candidates (Save old coordinates in case of rejection)
    std::vector<double> old_x(d,0);
    for (int i = 0; i < d; i++) old_x[i] = X[k][i];
    
    std::vector<double> shift_vec(d,0);
    if (d == 1) {
      const double r   = R::rnorm(0, sigma);
      shift_vec = {r};
    } else if (d == 2) {
      const double r   = R::rnorm(0, sigma);
      const double phi = R::runif(0, 1);
      shift_vec = {cos(phi*M_PI)*r, sin(phi*M_PI)*r};
    } else if (d == 3) {
      const double r   = R::rnorm(0, sigma);
      const double phi = R::runif(0, 1);
      const double theta=R::runif(0, 1);
      shift_vec = {cos(phi*M_PI)*sin(theta*M_PI)*r, sin(phi*M_PI)*sin(theta*M_PI)*r, cos(theta*M_PI)*r};
    }
    for(int i=0; i < d; i++){
      X[k][i] += shift_vec[i]; 
      // Periodic boundary wrap
      if      (X[k][i] < -l/2) X[k][i] += l;
      else if (X[k][i] > +l/2) X[k][i] -= l;
    }
    // ---------------------------------------------------------------------------
    for(int j = 0; j < nPart; j++){
      if(j == k) continue;
      
      double sum = 0;
      for(int l = 0; l < d; l++){
        sum += (X[k][l] - X[j][l])*(X[k][l] - X[j][l]);
      }
      
      dE[n] += q[k] * q[j] / sqrt(sum) + 1/ pow(sum, 4);
    }
    // ---------------------------------------------------------------------------
    if (n < T_burnIn_1) {
      t = std::exp(std::log(10) * (+2. + (double)n / T_burnIn_1 * (std::log10(t0)-2.)));
    } else {
      t = t0;
    }
    //if (n == T_burnIn_1-2) Rcout << "Temperature: " << t << ", t0:" << t0 << std::endl;
    if (R::runif(0,1) < std::exp(-1./t * dE[n])) {
      if (n > T_burnIn_2) accept_rate += 1;
      dE[n] = dE[n-1] + dE[n];
    } else {
      for (int i = 0; i < d; i++) X[k][i] = old_x[i];
      dE[n] = dE[n-1];
    }
  }
  //Rcout << "Acceptance: " << accept_rate/nIt*100. << "%" << std::endl;
  
  for (int i = 0; i < nPart; i++) {
    for (int j = 0; j < d; j++) {
      X0(i,j) = X[i][j];
    }
  }
  
  //Estimate Observables
  // ---------------------------------------------------------------------------
  //Estimate cluster sizes
  std::vector<int> cluster_sizes = estimate_clstSize(X); 
  
  //Estimate distances between particles of final config
  std::vector<double> particle_dists = estimate_particle_dists(X);
  
  //Estimate nearest neighbour distance of each particle of final config
  std::vector<double> nearest_neigh_dists = estimate_nearest_neighbour_dists(X);
  
  //Estimate number of particles in direct neighbourhood of final config
  std::vector<int> N_particles_neighbourhood = estimate_Nparticles_neighbourhood(X);
  // ---------------------------------------------------------------------------
  
  Rcpp::List rList = Rcpp::List::create(
    Named("E")          = dE, 
    Named("X")          = X0, 
    Named("Accpetance") = accept_rate/(nIt-T_burnIn_2), 
    Named("q")          = cluster_sizes, 
    Named("lambda")     = particle_dists, 
    Named("omega")      = nearest_neigh_dists, 
    Named("kappa")      = N_particles_neighbourhood
    );
  return rList;
}


/*** R
print("Loaded MCMC_functions.cpp.")
*/
