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

// ── Main model ─────────────────────────────────────────────────────────────
// [[Rcpp::export]]
NumericMatrix run_seir_cpp_v1(NumericVector pop,
                           double        I0,     // total infectious (I1+I2) at t=0
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
  
  // ── Compute R0 and exact growth rate r ──────────────────────────────────
  arma::vec dom_eigvec;
  double R0 = compute_R0(susc, inf_per, trans, cij, pop, vacc_cov, ve_inf, dom_eigvec);
  double r  = compute_r(R0, lat_rate, inf_rate);
  
  // ── Eigenvalue-consistent compartment ratios at rate r ──────────────────
  //   I2 = x,  I1 = x·(inf_rate+r)/inf_rate
  //   E2 = I1·(inf_rate+r)/lat_rate
  //   E1 = E2·(lat_rate+r)/lat_rate
  // ai_u = total seeded (E1+E2+I1+I2) for the unvaccinated stream of group i
  // => i2 = ai_u / denom_I,  where denom_I is the sum of all 4 compartment ratios
  double ratio_I1_I2 = (inf_rate + r) / inf_rate;
  double ratio_E2_I1 = (inf_rate + r) / lat_rate;
  double ratio_E1_E2 = (lat_rate + r) / lat_rate;
  double denom_I = 1.0 + ratio_I1_I2 + ratio_E2_I1 * ratio_I1_I2 + ratio_E1_E2 * ratio_E2_I1 * ratio_I1_I2;
  
  // I0 = total infectious (I1+I2+I1v+I2v); scale to total seeded
  double total_seeded = I0 * denom_I / (1.0 + ratio_I1_I2);
  
  // ── Initialise compartments ──────────────────────────────────────────────
  NumericVector S(ng),   E1(ng),  E2(ng),  I1(ng),  I2(ng),  R(ng),  V(ng),  cumI(ng);
  NumericVector Sv(ng),  E1v(ng), E2v(ng), I1v(ng), I2v(ng), Rv(ng), cumIv(ng);
  
  for (int i = 0; i < ng; i++) {
    // --- split seed between streams proportional to susceptible pool sizes ---
    double s_i  = (pop[i] - vacc_cov[i]) * susc[i];
    double sv_i = (1.0 - ve_inf[i]) * vacc_cov[i] * susc[i];
    double f_v  = sv_i / (s_i + sv_i);          // fraction in vaccinated stream
    
    double ai_u = total_seeded * dom_eigvec(i) * (1.0 - f_v);  // unvaccinated
    double ai_v = total_seeded * dom_eigvec(i) * f_v;           // vaccinated
    
    // --- unvaccinated stream ---
    double i2  = ai_u / denom_I;
    double i1  = ratio_I1_I2 * i2;
    double e2  = ratio_E2_I1 * i1;
    double e1  = ratio_E1_E2 * e2;
    
    // --- vaccinated stream ---
    double i2v = ai_v / denom_I;
    double i1v = ratio_I1_I2 * i2v;
    double e2v = ratio_E2_I1 * i1v;
    double e1v = ratio_E1_E2 * e2v;
    
    S[i]   = s_i  - (e1  + e2  + i1  + i2 );
    E1[i] = e1;
    E2[i] = e2;
    I1[i] = i1;
    I2[i] = i2;
    
    Sv[i]  = sv_i - (e1v + e2v + i1v + i2v);
    E1v[i] = e1v;
    E2v[i] = e2v;
    I1v[i] = i1v;
    I2v[i] = i2v;
    
    V[i]  = ve_inf[i] * vacc_cov[i];
    R[i]  = (pop[i] - vacc_cov[i]) * (1.0 - susc[i]);
    Rv[i] = (1.0 - ve_inf[i]) * vacc_cov[i] * (1.0 - susc[i]);
    
    cumI[i] = 0;
    cumIv[i] = 0;
    
  }
  
  // ── Output setup ─────────────────────────────────────────────────────────
  int n_steps        = (int)(t_end / dt) + 1;
  int n_out          = (int)(t_end)      + 1;
  int steps_per_unit = (int)std::round(1.0 / dt);
  NumericMatrix out(n_out, 2 * ng + 1);
  int out_row = 0;
  
  out(out_row, 0) = 0.0;
  for (int i = 0; i < ng; i++) {
    out(out_row, i + 1)      = cumI[i];
    out(out_row, i + 1 + ng) = cumIv[i];
  }
  out_row++;
  
  // ── ODE loop (unchanged) ─────────────────────────────────────────────────
  NumericVector I_tot(ng), lambda(ng), newInf(ng), newInfv(ng);
  NumericVector progE1(ng), progE2(ng), progI1(ng), progI2(ng);
  NumericVector progE1v(ng), progE2v(ng), progI1v(ng), progI2v(ng);
  
  double t = 0.0;
  for (int step = 0; step < n_steps - 1; step++) {
    
    for (int i = 0; i < ng; i++) I_tot[i] = I1[i] + I2[i] + I1v[i] + I2v[i];
    
    for (int i = 0; i < ng; i++) {
      double s = 0.0;
      for (int j = 0; j < ng; j++) s += cij(i, j) * I_tot[j];
      lambda[i] = trans * s;
    }
    
    for (int i = 0; i < ng; i++) {
      newInf[i]  = lambda[i] * S[i];
      progE1[i] = lat_rate * E1[i];
      progE2[i]  = lat_rate * E2[i];
      progI1[i] = inf_rate * I1[i];
      progI2[i]  = inf_rate * I2[i];
      newInfv[i] = lambda[i] * Sv[i];
      progE1v[i] = lat_rate * E1v[i];
      progE2v[i] = lat_rate * E2v[i];
      progI1v[i] = inf_rate * I1v[i];
      progI2v[i] = inf_rate * I2v[i];
    }
    
    for (int i = 0; i < ng; i++) {
      double dn = dt * newInf[i], dnv = dt * newInfv[i];
      S[i]   -= dn;
      E1[i]  += dn          - dt * progE1[i];
      E2[i]  += dt * (progE1[i]  - progE2[i]);
      I1[i]  += dt * (progE2[i]  - progI1[i]);
      I2[i]  += dt * (progI1[i]  - progI2[i]);
      R[i]   += dt *  progI2[i];
      Sv[i]  -= dnv;
      E1v[i] += dnv         - dt * progE1v[i];
      E2v[i] += dt * (progE1v[i] - progE2v[i]);
      I1v[i] += dt * (progE2v[i] - progI1v[i]);
      I2v[i] += dt * (progI1v[i] - progI2v[i]);
      Rv[i]  += dt *  progI2v[i];
      cumI[i]  += dt * progE2[i];
      cumIv[i] += dt * progE2v[i];
    }
    
    t += dt;
    if ((step + 1) % steps_per_unit == 0 && out_row < n_out) {
      out(out_row, 0) = std::round(t);
      for (int i = 0; i < ng; i++) {
        out(out_row, i + 1)      = cumI[i];
        out(out_row, i + 1 + ng) = cumIv[i];
      }
      out_row++;
    }
  }
  
  return out;
}


