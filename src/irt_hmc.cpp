// 2PL IRT item-parameter sampler (alpha, beta | theta) wired to WALNUTS.
//
// This is the compiled engine behind bairrtt. It draws the IRT item parameters
// (discrimination alpha, difficulty beta) conditional on the person abilities
// theta, using the gradient-based WALNUTS sampler (a NUTS variant). theta is
// held FIXED for this update: in the surrounding model it feeds non-
// differentiable BART surfaces and is sampled elsewhere by Metropolis, but
// conditional on theta the (alpha, beta) log-posterior is smooth, so a gradient
// sampler applies.
//
// The target (IrtLogpGrad) and the persistent sampler state (IrtSampler) live
// in bairrtt_types.h so that they are complete types in both this file and the
// generated RcppExports.cpp; this file adds the lifecycle entry points:
//   create_sampler / warmup_draw / freeze / draw_sample.
// warmup_draw drives AdaptiveWalnuts (adapting step size + mass); freeze() hands
// off to a non-adapting WalnutsSampler; draw_sample drives that.
//
// The exported symbols are dot-prefixed and treated as internal; bairrtt wraps
// them in validated, documented R functions (see R/engine.R), and the C side
// ASSUMES those preconditions hold -- e.g. create_sampler does log(alpha_init)
// without checking positivity, and to_natural assumes an even-length position.
// Call the `.irt_*` symbols directly at your own risk. WALNUTS needs
// C++20 (concepts) and Eigen -- both are arranged by the package build:
// CXX_STD = CXX20 in src/Makevars, Eigen via LinkingTo: RcppEigen, and the
// vendored WALNUTS headers via -I../inst/include.

#include "../inst/include/bairrtt_types.h"

// Exposed so the R side can finite-difference check the gradient and inspect
// the target. Not on the sampling hot path.
// [[Rcpp::export(.irt_logdensity_grad)]]
Rcpp::List irt_logdensity_grad(const Eigen::Map<VectorXd> par,
                               const Eigen::Map<MatrixXd> Y,
                               const Eigen::Map<VectorXd> theta,
                               double beta_sd = 10.0) {
  IrtLogpGrad f(Y, theta, beta_sd);
  double lp = 0.0;
  VectorXd grad;
  f(par, lp, grad);
  return Rcpp::List::create(Rcpp::_["value"] = lp, Rcpp::_["gradient"] = grad);
}

// Lifecycle (interleaved warmup, then freeze):
//   create_sampler()                     -- once, before the outer Gibbs loop
//   for each OUTER BURN-IN scan:          -- adaptation sees the moving theta
//     warmup_draw(s, theta)               -- one AdaptiveWalnuts step; tuning adapts
//   freeze(s)                             -- once; AdaptiveWalnuts::sampler() handoff
//   for each OUTER SAMPLING scan:
//     draw_sample(s, theta)               -- one fixed-tuning WalnutsSampler step
//
// Both warmup_draw and draw_sample condition on the current theta (the functor
// is held by const reference inside WALNUTS, so updating theta in place is seen
// by the next transition) and return the (alpha, beta) draw on the natural
// scale. Warmup draws are NOT valid for inference but are still needed to
// condition the rest of each outer burn-in scan.

// (alpha_raw, beta) -> (alpha, beta) on the natural scale, for returning to R.
static VectorXd to_natural(const VectorXd& position) {
  const int n_items = static_cast<int>(position.size() / 2);
  VectorXd out(2 * n_items);
  out.head(n_items) = position.head(n_items).array().exp().matrix();  // alpha = exp(alpha_raw)
  out.tail(n_items) = position.tail(n_items);                         // beta
  return out;
}

// Build the sub-sampler once, before the outer Gibbs loop, and start it in the
// adapting phase.
// [[Rcpp::export(.irt_create_sampler)]]
Rcpp::XPtr<IrtSampler> create_sampler(const Eigen::Map<MatrixXd> Y,
                                      const Eigen::Map<VectorXd> alpha_init,
                                      const Eigen::Map<VectorXd> beta_init,
                                      double beta_sd = 10.0,
                                      double step_size = 0.1,
                                      int seed = 1) {
  const int n_items = static_cast<int>(Y.cols());
  Rcpp::XPtr<IrtSampler> s(new IrtSampler(), true);

  s->logp_grad = IrtLogpGrad(Y, VectorXd::Zero(Y.rows()), beta_sd);  // theta set per draw
  s->position.resize(2 * n_items);
  s->position.head(n_items) = alpha_init.array().log().matrix();     // to unconstrained scale
  s->position.tail(n_items) = beta_init;
  s->step_size = step_size;
  s->rng.seed(static_cast<std::mt19937_64::result_type>(seed));

  // Start adaptation: identity initial mass, user-supplied initial step size.
  walnuts::InitChainConfig init(s->step_size, s->position,
                                VectorXd::Ones(2 * n_items));
  s->adapter.emplace(s->rng, s->handler, s->logp_grad, init, s->warmup_cfg,
                     s->sampling_cfg);
  return s;
}

// One ADAPTING step, conditioning on the current theta. Call once per outer
// burn-in scan. Returns the (alpha, beta) draw -- invalid for inference, but
// needed to condition the rest of the scan. Adaptation averages over the
// sequence of theta it is shown.
// [[Rcpp::export(.irt_warmup_draw)]]
VectorXd warmup_draw(Rcpp::XPtr<IrtSampler> s, const Eigen::Map<VectorXd> theta) {
  if (s->frozen)   Rcpp::stop("warmup_draw() called after freeze(); use draw_sample()");
  if (!s->adapter) Rcpp::stop("sampler not initialized");
  s->logp_grad.theta = theta;   // condition on this scan's theta

  (*s->adapter)();              // one adaptive transition; updates tuning
  s->position = s->handler.position;
  ++s->warmup_iter;
  return to_natural(s->position);
}

// Freeze adaptation: hand the adapted tuning to a non-adapting sampler. Call
// once, between the outer burn-in and sampling phases.
// [[Rcpp::export(.irt_freeze)]]
void freeze(Rcpp::XPtr<IrtSampler> s) {
  if (s->frozen)   Rcpp::stop("freeze() called twice");
  if (!s->adapter) Rcpp::stop("sampler not initialized");
  s->sampler.emplace(s->adapter->sampler());  // AdaptiveWalnuts -> WalnutsSampler
  s->adapter.reset();                          // adaptation done; safe (rng/handler/functor owned here)
  s->frozen = true;
}

// One FIXED-TUNING step, conditioning on the current theta. Call once per outer
// sampling scan, after freeze(). Returns the (alpha, beta) draw to keep.
// [[Rcpp::export(.irt_draw_sample)]]
VectorXd draw_sample(Rcpp::XPtr<IrtSampler> s, const Eigen::Map<VectorXd> theta) {
  if (!s->frozen)  Rcpp::stop("draw_sample() called before freeze(); use warmup_draw()");
  s->logp_grad.theta = theta;   // condition on this scan's theta

  (*s->sampler)();              // one fixed-tuning transition
  s->position = s->handler.position;
  return to_natural(s->position);
}

// Inspect the current tuning -- useful for watching adaptation during burn-in.
// [[Rcpp::export(.irt_sampler_tuning)]]
Rcpp::List sampler_tuning(Rcpp::XPtr<IrtSampler> s) {
  return Rcpp::List::create(
      Rcpp::_["frozen"]      = s->frozen,
      Rcpp::_["warmup_iter"] = static_cast<double>(s->warmup_iter),
      Rcpp::_["step_size"]   = s->handler.step_size,
      Rcpp::_["inv_mass"]    = s->handler.inv_mass);
}
