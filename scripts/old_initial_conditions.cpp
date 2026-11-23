// seir_model.cpp
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix get_initial_conditions(NumericVector pop,
                                     double I0_total,
                                     NumericVector vacc_cov,
                                     NumericVector ve_inf,
                                     double trans,
                                     NumericVector susc,
                                     double lat_per,
                                     double inf_per,
                                     NumericMatrix cij,
                                     double dt = 0.01) {
  
  int ng = pop.size();
  double lat_rate = 2.0 / lat_per;
  double inf_rate = 2.0 / inf_per;
  
  // Total unvaccinated susceptibles for proportional seeding
  double pop_unvacc_susc = 0.0;
  for (int i = 0; i < ng; i++)
    pop_unvacc_susc += (pop[i] - vacc_cov[i]) * susc[i];
  
  // Seed 1 exposed person proportionally across groups
  NumericVector S(ng), E1(ng), E2(ng), I1(ng), I2(ng), R(ng), V(ng);
  NumericVector Sv(ng), E1v(ng), E2v(ng), I1v(ng), I2v(ng), Rv(ng);
  for (int i = 0; i < ng; i++) {
    E1[i] = ((pop[i] - vacc_cov[i]) * susc[i]) / pop_unvacc_susc;
    S[i]  = ((pop[i] - vacc_cov[i]) * susc[i]) - E1[i];  // reduce S by E1 share, not by 1
    R[i]  = (pop[i] - vacc_cov[i]) * (1.0 - susc[i]);
    Sv[i] = (1.0 - ve_inf[i]) * vacc_cov[i] * susc[i];
    Rv[i] = (1.0 - ve_inf[i]) * vacc_cov[i] * (1.0 - susc[i]);
    V[i]  = ve_inf[i] * vacc_cov[i];
  }
  
  NumericVector I_tot(ng), lambda(ng), newInf(ng), newInfv(ng);
  NumericVector progE1(ng), progE2(ng), progI1(ng), progI2(ng);
  NumericVector progE1v(ng), progE2v(ng), progI1v(ng), progI2v(ng);
  
  bool keep_going = true;
  int max_steps = 10000000;  // safety limit: ~100,000 days at dt=0.01
  for (int step = 0; step < max_steps && keep_going; step++) {
    
    for (int i = 0; i < ng; i++) I_tot[i] = I1[i] + I2[i] + I1v[i] + I2v[i];
    
    for (int i = 0; i < ng; i++) {
      double s = 0.0;
      for (int j = 0; j < ng; j++) s += cij(i, j) * I_tot[j];
      lambda[i] = trans * s;
    }
    
    for (int i = 0; i < ng; i++) {
      newInf[i]  = lambda[i] * S[i];   progE1[i] = lat_rate * E1[i];
      progE2[i]  = lat_rate * E2[i];   progI1[i] = inf_rate * I1[i];
      progI2[i]  = inf_rate * I2[i];
      newInfv[i] = lambda[i] * Sv[i];  progE1v[i] = lat_rate * E1v[i];
      progE2v[i] = lat_rate * E2v[i];  progI1v[i] = inf_rate * I1v[i];
      progI2v[i] = inf_rate * I2v[i];
    }
    
    for (int i = 0; i < ng; i++) {
      double dn = dt * newInf[i], dnv = dt * newInfv[i];
      S[i]   -= dn;
      E1[i]  += dn - dt * progE1[i];
      E2[i]  += dt * (progE1[i]  - progE2[i]);
      I1[i]  += dt * (progE2[i]  - progI1[i]);
      I2[i]  += dt * (progI1[i]  - progI2[i]);
      R[i]   += dt * progI2[i];
      Sv[i]  -= dnv;
      E1v[i] += dnv - dt * progE1v[i];
      E2v[i] += dt * (progE1v[i] - progE2v[i]);
      I1v[i] += dt * (progE2v[i] - progI1v[i]);
      I2v[i] += dt * (progI1v[i] - progI2v[i]);
      Rv[i]  += dt * progI2v[i];
    }
    
    // Stop when total infectious reaches target
    double I_total = 0.0;
    for (int i = 0; i < ng; i++) I_total += I1[i] + I2[i] + I1v[i] + I2v[i];
    if (I_total >= I0_total) keep_going = false;
  }
  
  // Return ng x 12 matrix: columns = S, E1, E2, I1, I2, R, Sv, E1v, E2v, I1v, I2v, Rv
  NumericMatrix out(ng, 12);
  for (int i = 0; i < ng; i++) {
    out(i, 0) = S[i];
    out(i, 1) = E1[i];
    out(i, 2) = E2[i];
    out(i, 3) = I1[i];
    out(i, 4) = I2[i];
    out(i, 5) = R[i];
    out(i, 6) = Sv[i];
    out(i, 7) = E1v[i];
    out(i, 8) = E2v[i];
    out(i, 9) = I1v[i];
    out(i, 10) = I2v[i];
    out(i, 11) = Rv[i];
  }
  
  return out;
  
}


