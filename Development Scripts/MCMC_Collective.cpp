// [[Rcpp::plugins("cpp11")]]

// ---------------------------------------------------------------------------
// MCMC_Collective.cpp is an alternative, though deprecated by now, version
// MCMC_functions.cpp in which we implemented a collective sampling
// variant of the Random-Walk Metropolis algorithm that we developed
// ourselves. It's efficiency prooved to be equal to the simple approach while
// requiring a lot more computation power. Because of this the ansatz was
// discarded
// ---------------------------------------------------------------------------

#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// =======================================================
// ========================================= Calculate Energy
// =======================================================
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
// [[Rcpp::export]]
double getDist(std::vector<double> a, std::vector<double> b, const double l) {
  const int dim = a.size();
  double dist = 0;
  for (int i = 0; i < dim; i++) {
    if        ((a[i] - b[i]) > +l/2){
      b[i] += l;
    } else if ((a[i] - b[i]) < -l/2){
      b[i] -= l;
    }
    dist += (a[i] - b[i])*(a[i] - b[i]);
  }
  return std::sqrt(dist);
}

double estimate_clstSize(std::vector<std::vector<double>> X, const double L) {
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

        if (getDist(X[clst_idx], X[free_idx], L) < min_dist) {
          clusters.back()[cluster_sizes.back()]  = free_idx;
          cluster_sizes.back()  += 1; // Current cluster has grown by 1
          idx[free_idx]          = 0; // Particle has been assigned to a cluster
          nFree                 -= 1; // One further particle is assigned to a cluster
          k                     -= 1; // Neccessary for correct procedure
        }

      }
    }
  }

  // Calc mean cluster size and return
  double mean_size = 0;
  for (unsigned int i = 0; i < cluster_sizes.size(); i++) {
    mean_size += cluster_sizes[i];
  }

  return mean_size/cluster_sizes.size();
}

// =======================================================
// ========================================= MCMC_CPP
// =======================================================
// [[Rcpp::export]]
Rcpp::List MCMC_CPP_Collective(int nIt, int nPart, double vol, double t, double sigma, NumericMatrix X0, NumericVector q){
  int d = X0(1,Rcpp::_).size();
  const double l = std::pow(vol, 1./d);
  const double t0= t;
  const double p_bond = 0.5;

  const double T_burnIn_1 = 1e5;
  const double T_burnIn_2 = 1e5;
  const double min_dist = 2;
  // ----------------------- Observables
  NumericVector dE(nIt, 0.0);
  dE[0] = E_CPP(X0, q);

  const int N_clstSize_samples = 10;
  int k_clstSize               = 0;
  const int n_min_clstSize           = nIt-1e4;
  std::vector<double> clstSize(N_clstSize_samples, -1);
  const double stp_clstSize    = (nIt - n_min_clstSize)/N_clstSize_samples;
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
    // ---------------------------------------------------------------------------
    // Sample bonds to create a cluster
    k = std::floor(R::runif(0,nPart));
    int nFree = nPart-1;
    //Rcout << "Markov Chain step: " << n << ", seedparticle: " << k << std::endl;

    std::vector<int> bonded_particles;
    bonded_particles.push_back(k);
    std::vector<int> id(nPart, 1);
    id[k] = 0;
    std::vector<int> not_bonded_particles;
    std::vector<std::vector<int>> failed_bonds;

    for(unsigned int i = 0; i < bonded_particles.size(); i++){
      for(int j = 0; j < nFree; j++){
        int free_id = getFreeIndex(id, j+1);
        if(getDist(X[bonded_particles[i]], X[free_id], l) < min_dist){
          if(R::runif(0,1) < p_bond){
            bonded_particles.push_back(free_id);
            id[free_id] = 0;
            nFree--;
            j--;
          }else{
            std::vector<int> tmp = {bonded_particles[i], free_id};
            failed_bonds.push_back(tmp);
          }
        }
      }
    }

    // Determine how many failed outgoing bonds exist
    int n_out_bef = 0;
    for(unsigned int i = 0; i < failed_bonds.size(); i++){
      for(unsigned int j = 0; j < bonded_particles.size(); j++){
        if(failed_bonds[i][1] == bonded_particles[j]) continue;
        n_out_bef++;
      }
    }

    // Save the set of particles not assigned to the cluster for easier access
    for(int u = 0; u < nPart; u++){
      for(unsigned int j = 0; j < bonded_particles.size(); j++){
        if(u == bonded_particles[j]) continue;
        not_bonded_particles.push_back(u);
      }
    }

    // ---------------------------------------------------------------------------
    // Calculate the necessary energy terms before shifting
    for  (unsigned int i = 0; i < bonded_particles.size()    ; i++){
      for(unsigned int j = 0; j < not_bonded_particles.size(); j++){
        const int a = bonded_particles[i];
        const int b = not_bonded_particles[j];
        double sum = 0;
        for(int dim = 0; dim < d; dim++){
          sum += (X[a][dim] - X[b][dim])*(X[a][dim] - X[b][dim]);
        }

        dE[n] -= q[a] * q[b] / sqrt(sum) + 1/ pow(sum, 4);
      }
    }

    // ---------------------------------------------------------------------------
    // Make new candidate by shifting (Save old coordinates in case of rejection)
    std::vector<std::vector<double>> old_vals;
    std::vector<double> old_x(d,0);
    for(unsigned int w = 0; w < bonded_particles.size(); w++){
      for (int i = 0; i < d; i++) old_x[i] = X[bonded_particles[w]][i];
      old_vals.push_back(old_x);
    }


    // Make shift vector
    std::vector<double> shift_vec(d,0);
    if        (d == 1) {
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

    // Apply shift vector to all cluster particles
    for(unsigned int p = 0; p < bonded_particles.size(); p++){
      const int a = bonded_particles[p];
      for(int i=0; i < d; i++){
        X[a][i] += shift_vec[i];
        // Periodic boundary wrap
        if      (X[a][i] < -l/2) X[a][i] += l;
        else if (X[a][i] > +l/2) X[a][i] -= l;
      }
    }

    // ---------------------------------------------------------------------------
    // Calculate the necessary energy terms after shifting
    for(  unsigned int i = 0; i < bonded_particles.size()    ; i++){
      for(unsigned int j = 0; j < not_bonded_particles.size(); j++){
        double sum = 0;
        const int a = bonded_particles[i];
        const int b = not_bonded_particles[j];
        for(int dim = 0; dim < d; dim++){
          sum += (X[a][dim] - X[b][dim])*(X[a][dim] - X[b][dim]);
        }

        dE[n] += q[a] * q[b] / sqrt(sum) + 1/ pow(sum, 4);
      }
    }

    // ---------------------------------------------------------------------------
    // Determine outgoing failed bonds after shift
    int n_out_aft = 0;
    for(  unsigned int i = 0; i < bonded_particles.size()    ; i++){
      for(unsigned int j = 0; j < not_bonded_particles.size(); j++){
       if(getDist(X[bonded_particles[i]], X[not_bonded_particles[j]], l) < min_dist){
         n_out_aft++;
       }
      }
    }

    //if (bonded_particles.size() > 2) {
    //  Rcout << n << ": Determined cluster " << bonded_particles.size() << ", " << n_out_bef << ", " << n_out_aft << std::endl;
    //}

    // ---------------------------------------------------------------------------
    // To accelerate the burn-in, while in the first T_burnIn_1 steps the temperature is reduced from 10^2 to t0 exponentially
    if (n < T_burnIn_1) {
      t = std::exp(std::log(10) * (+2. + (double)n / T_burnIn_1 * (std::log10(t0)-2.)));
    } else {
      t = t0;
    }

    // Accept-Reject step of new candidate
    if (R::runif(0,1) < std::exp(-1./t * dE[n])*pow(1-p_bond, (n_out_aft - n_out_bef))) {
    // if (R::runif(0,1) < std::exp(-1./t * dE[n])) {
      if (n > T_burnIn_2) accept_rate += 1;
    } else {
      dE[n] = 0;
      for(unsigned int w = 0; w < bonded_particles.size(); w++){
        const int a = bonded_particles[w];
        for (int i = 0; i < d; i++) X[a][i] = old_vals[w][i];
      }
    }
    // ---------------------------------------------------------------------------
    if ((n-n_min_clstSize) >= k_clstSize * stp_clstSize) {
      clstSize[k_clstSize] = estimate_clstSize(X, l);
      k_clstSize++;
    }
    // ---------------------------------------------------------------------------
  }


  for (int i = 0; i < nPart; i++) {
    for (int j = 0; j < d; j++) {
      X0(i,j) = X[i][j];
    }
  }

  // Calc mean cluster size
  double mean_size = 0;
  for (unsigned int i = 0; i < clstSize.size(); i++) {
    mean_size += clstSize[i];
  }
  mean_size = mean_size/clstSize.size();

  Rcpp::List rList = Rcpp::List::create(Named("dE")=dE, Named("X")=X0, Named("Accpetance")=accept_rate/(nIt-T_burnIn_2), Named("MeanClusterSize")=mean_size);
  return rList;
}

/*** R
print("Loaded MCMC_Collective.cpp.")
*/
