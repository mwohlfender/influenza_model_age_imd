// seir_model.cpp
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// ── Helper: R0 from NGM dominant eigenvalue ────────────────────────────────
// Mirrors R0_func(): ngm[j,k] = trans * susc[k] * cij[j,k] * pop[k] * inf_per
// (column k = source group; population-weights the contact matrix)
// Returns R0; fills dom_eigvec with the normalised dominant right eigenvector
static double compute_R0(const NumericVector& susc,
                         double               inf_per,
                         double               trans,
                         const NumericMatrix& cij,
                         const NumericVector& pop,
                         const NumericVector& vacc_cov,
                         const NumericVector& ve_inf,
                         arma::vec&           dom_eigvec) {   // <-- new output arg
  int ng = pop.size();
  arma::mat ngm(ng, ng);
  for (int k = 0; k < ng; k++) {
    double eff_pop_k = pop[k] - ve_inf[k] * vacc_cov[k];
    for (int j = 0; j < ng; j++)
      ngm(j, k) = trans * susc[k] * cij(j, k) * eff_pop_k * inf_per;
  }
  
  arma::cx_vec eigvals;
  arma::cx_mat eigvec;
  arma::eig_gen(eigvals, eigvec, ngm);   // right eigenvectors in columns
  
  double R0 = 0.0;
  int    dom_idx = 0;
  for (int i = 0; i < (int)eigvals.n_elem; i++) {
    if (std::abs(eigvals(i).imag()) < 1e-6 && eigvals(i).real() > R0) {
      R0      = eigvals(i).real();
      dom_idx = i;
    }
  }
  
  // Extract dominant right eigenvector; take abs() to guarantee positive entries
  dom_eigvec.set_size(ng);
  for (int i = 0; i < ng; i++)
    dom_eigvec(i) = std::abs(eigvec(i, dom_idx).real());
  dom_eigvec /= arma::sum(dom_eigvec);   // normalise to sum = 1
  
  return R0;
}

// ── Helper: exact r for 2E-2I Erlang model ─────────────────────────────────
// Euler-Lotka for 2 Erlang stages each:
//   R0 = [(1 + r/lat_rate)(1 + r/inf_rate)]^2
// Square-root → quadratic in r:
//   r²/(lat_rate·inf_rate) + r·(1/lat_rate + 1/inf_rate) + (1 − √R0) = 0
static double compute_r(double R0, double lat_rate, double inf_rate) {
  if (R0 <= 1.0) {
    Rcpp::warning("R0 <= 1: epidemic cannot grow; r set to 0.");
    return 0.0;
  }
  double a =  1.0 / (lat_rate * inf_rate);
  double b =  1.0 / lat_rate + 1.0 / inf_rate;
  double c =  1.0 - std::sqrt(R0);           // c < 0 when R0 > 1
  // Positive root of the quadratic
  return (-b + std::sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
}

// ── Exported: standalone R0 calculator (matches R0_func) ───────────────────
// [[Rcpp::export]]
double compute_R0_cpp(NumericVector susc, double inf_per, double trans,
                      NumericMatrix cij,  NumericVector pop,
                      NumericVector vacc_cov, NumericVector ve_inf) {
  arma::vec tmp;
  return compute_R0(susc, inf_per, trans, cij, pop, vacc_cov, ve_inf, tmp);
}


// [[Rcpp::export]]
NumericMatrix run_seir_cpp_v2(NumericVector pop,
                              double        I0,      // target total infectious (I1+I2+I1v+I2v) at t = 0
                              NumericVector vacc_cov,
                              NumericVector ve_inf,
                              double        trans,
                              NumericVector susc,
                              double        lat_per,
                              double        inf_per,
                              NumericMatrix cij,
                              double        t_end = 250.0,
                              double        dt    = 0.001) {
  
  int    ng       = pop.size();
  double lat_rate = 2.0 / lat_per;
  double inf_rate = 2.0 / inf_per;
  
  // ── Initialise compartments at disease-free equilibrium ──────────────────
  NumericVector S(ng),  E1(ng),  E2(ng),  I1(ng),  I2(ng),  R(ng),  V(ng),  cumI(ng);
  NumericVector Sv(ng), E1v(ng), E2v(ng), I1v(ng), I2v(ng), Rv(ng), cumIv(ng);
  
  for (int i = 0; i < ng; i++) {
    S[i]  = (pop[i] - vacc_cov[i]) * susc[i];
    R[i]  = (pop[i] - vacc_cov[i]) * (1.0 - susc[i]);
    Sv[i] = (1.0 - ve_inf[i]) * vacc_cov[i] * susc[i];
    Rv[i] = (1.0 - ve_inf[i]) * vacc_cov[i] * (1.0 - susc[i]);
    V[i]  = ve_inf[i] * vacc_cov[i];
    // All E and I compartments remain at zero; set by burn-in below
  }
  
  // ── Seed: 1 exposed person in group 0, unvaccinated stream ───────────────
  // The burn-in ODE will naturally propagate this seed across all groups and
  // both streams, converging to the dominant eigenvector distribution.
  E1[0] = 1.0;
  S[0] -= 1.0;
  
  // ── Derivative vectors (reused in both loops) ─────────────────────────────
  NumericVector I_tot(ng), lambda(ng), newInf(ng), newInfv(ng);
  NumericVector progE1(ng), progE2(ng), progI1(ng), progI2(ng);
  NumericVector progE1v(ng), progE2v(ng), progI1v(ng), progI2v(ng);
  
  // ── Burn-in: run ODE until total infectious reaches I0 ───────────────────
  // After several serial intervals the system settles into the dominant
  // eigenvector distribution regardless of where the seed was placed.
  int  max_burnin_steps = (int)(180.0 / dt);  // cap at 180 days
  bool target_reached   = false;
  double I_total = 0.0;
  
  for (int step = 0; step < max_burnin_steps; step++) {
    
    for (int i = 0; i < ng; i++) I_tot[i] = I1[i]+I2[i]+I1v[i]+I2v[i];
    for (int i = 0; i < ng; i++) {
      double s = 0.0;
      for (int j = 0; j < ng; j++) s += cij(i, j) * I_tot[j];
      lambda[i] = trans * s;
    }
    for (int i = 0; i < ng; i++) {
      newInf[i]  = lambda[i]*S[i];
      progE1[i]  = lat_rate*E1[i];
      progE2[i]  = lat_rate*E2[i];
      progI1[i]  = inf_rate*I1[i];
      progI2[i]  = inf_rate*I2[i];
      
      newInfv[i] = lambda[i]*Sv[i];
      progE1v[i] = lat_rate*E1v[i];
      progE2v[i] = lat_rate*E2v[i];
      progI1v[i] = inf_rate*I1v[i];
      progI2v[i] = inf_rate*I2v[i];
    }
    for (int i = 0; i < ng; i++) {
      double dn = dt*newInf[i], dnv = dt*newInfv[i];
      S[i]   -= dn;
      E1[i]  += dn             - dt*progE1[i];
      E2[i]  += dt*(progE1[i]  - progE2[i]);
      I1[i]  += dt*(progE2[i]  - progI1[i]);
      I2[i]  += dt*(progI1[i]  - progI2[i]);
      R[i]   += dt* progI2[i];
      Sv[i]  -= dnv;
      E1v[i] += dnv            - dt*progE1v[i];
      E2v[i] += dt*(progE1v[i] - progE2v[i]);
      I1v[i] += dt*(progE2v[i] - progI1v[i]);
      I2v[i] += dt*(progI1v[i] - progI2v[i]);
      Rv[i]  += dt* progI2v[i];
    }
    
    // Check stopping condition after each step
    I_total = 0.0;
    for (int i = 0; i < ng; i++) I_total += I1[i]+I2[i]+I1v[i]+I2v[i];
    if (I_total >= I0) { target_reached = true; break; }
  }
  
  if (!target_reached) {
    
    Rcpp::warning("Burn-in: total infectious did not reach I0 within 180 days. Check R0 > 1.");
    
    arma::vec dom_eigvec;
    double R0 = compute_R0(susc, inf_per, trans, cij, pop, vacc_cov, ve_inf, dom_eigvec);
    
    Rcpp::Rcout << "R0: " << R0 << " I_total: " << I_total << " I0: " << I0 << " trans: " << trans << std::endl;
    
  }
  
  // ── Scale all E and I compartments so that total infectious = I0 ──────────
  // S and Sv are adjusted to absorb the extra seeded and preserve the exact
  // population balance (S + E1+E2+I1+I2 + R = pop - vacc_cov for each group).
  // R, Rv and V are left unchanged.
  double I_burnin = 0.0;
  for (int i = 0; i < ng; i++) I_burnin += I1[i]+I2[i]+I1v[i]+I2v[i];
  
  double scale = (I_burnin > 0.0) ? I0 / I_burnin : 1.0;
  
  for (int i = 0; i < ng; i++) {
    double EI_pre  = E1[i]  + E2[i]  + I1[i]  + I2[i];
    double EIv_pre = E1v[i] + E2v[i] + I1v[i] + I2v[i];
    
    E1[i]  *= scale;  E2[i]  *= scale;  I1[i]  *= scale;  I2[i]  *= scale;
    E1v[i] *= scale;  E2v[i] *= scale;  I1v[i] *= scale;  I2v[i] *= scale;
    
    // Absorb the additional seeded individuals from S and Sv
    S[i]  -= (E1[i]  + E2[i]  + I1[i]  + I2[i])  - EI_pre;
    Sv[i] -= (E1v[i] + E2v[i] + I1v[i] + I2v[i]) - EIv_pre;
  }
  
  // ── Set cumI at t = 0 to those currently infectious ──────────────────────
  for (int i = 0; i < ng; i++) {
    cumI[i]  = I1[i] + I2[i];
    cumIv[i] = I1v[i] + I2v[i];
  }
  
  // ── Output setup ─────────────────────────────────────────────────────────
  int n_steps        = (int)(t_end / dt) + 1;
  int n_out          = (int)(t_end)      + 1;
  int steps_per_unit = (int)std::round(1.0 / dt);
  NumericMatrix out(n_out, 2*ng + 1);
  int out_row = 0;
  
  out(out_row, 0) = 0.0;
  for (int i = 0; i < ng; i++) {
    out(out_row, i+1)    = cumI[i];
    out(out_row, i+1+ng) = cumIv[i];
  }
  out_row++;
  
  // ── Main ODE loop ─────────────────────────────────────────────────────────
  double t = 0.0;
  for (int step = 0; step < n_steps - 1; step++) {
    
    for (int i = 0; i < ng; i++) I_tot[i] = I1[i]+I2[i]+I1v[i]+I2v[i];
    for (int i = 0; i < ng; i++) {
      double s = 0.0;
      for (int j = 0; j < ng; j++) s += cij(i, j) * I_tot[j];
      lambda[i] = trans * s;
    }
    for (int i = 0; i < ng; i++) {
      newInf[i]  = lambda[i]*S[i];    progE1[i]  = lat_rate*E1[i];
      progE2[i]  = lat_rate*E2[i];    progI1[i]  = inf_rate*I1[i];
      progI2[i]  = inf_rate*I2[i];
      newInfv[i] = lambda[i]*Sv[i];   progE1v[i] = lat_rate*E1v[i];
      progE2v[i] = lat_rate*E2v[i];   progI1v[i] = inf_rate*I1v[i];
      progI2v[i] = inf_rate*I2v[i];
    }
    for (int i = 0; i < ng; i++) {
      double dn = dt*newInf[i], dnv = dt*newInfv[i];
      S[i]   -= dn;
      E1[i]  += dn             - dt*progE1[i];
      E2[i]  += dt*(progE1[i]  - progE2[i]);
      I1[i]  += dt*(progE2[i]  - progI1[i]);
      I2[i]  += dt*(progI1[i]  - progI2[i]);
      R[i]   += dt* progI2[i];
      Sv[i]  -= dnv;
      E1v[i] += dnv            - dt*progE1v[i];
      E2v[i] += dt*(progE1v[i] - progE2v[i]);
      I1v[i] += dt*(progE2v[i] - progI1v[i]);
      I2v[i] += dt*(progI1v[i] - progI2v[i]);
      Rv[i]  += dt* progI2v[i];
      cumI[i]  += dt * progE2[i];    // E2 → I1 flow (onset of infectiousness)
      cumIv[i] += dt * progE2v[i];
    }
    
    t += dt;
    if ((step+1) % steps_per_unit == 0 && out_row < n_out) {
      out(out_row, 0) = std::round(t);
      for (int i = 0; i < ng; i++) {
        out(out_row, i+1)    = cumI[i];
        out(out_row, i+1+ng) = cumIv[i];
      }
      out_row++;
    }
  }
  
  return out;
}


