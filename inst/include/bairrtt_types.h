// Shared C++ definitions for the WALNUTS item-parameter engine.
//
// This header is auto-included into the Rcpp-generated src/RcppExports.cpp
// (Rcpp::compileAttributes pulls in <pkg>_types.h) as well as into
// src/irt_hmc.cpp, so anything needed by the exported function *signatures*
// lives here.
//
// IrtSampler is defined here IN FULL (not merely forward-declared) on purpose:
// it is handled as Rcpp::XPtr<IrtSampler>, whose default finalizer runs
// `delete` on the pointer. If RcppExports.cpp saw only a forward declaration,
// that `delete` would act on an incomplete type (-Wdelete-incomplete) and
// -- because the finalizer template has external linkage and is also
// instantiated in irt_hmc.cpp -- linking could pick the incomplete-type
// instantiation, skipping ~IrtSampler and leaking the WALNUTS/Eigen state on
// every garbage collection. A complete type in every translation unit avoids
// that.
#pragma once

#include <cmath>     // std::exp, std::log1p
#include <cstddef>   // std::size_t
#include <optional>
#include <random>

// RcppEigen's bundled Eigen and the WALNUTS headers trip a handful of warnings
// under -Wall -Wextra that we neither own nor can fix. Silence them just around
// these includes (pattern from stan4bart/src/interruptable_sampler.hpp). Once
// RcppEigen.h is included here (behind its include guard), RcppExports.cpp's own
// later `#include <RcppEigen.h>` is a no-op, so this covers both translation
// units.
#if (defined(__clang__) && (__clang_major__ > 3 || (__clang_major__ == 3 && __clang_minor__ >= 7))) || \
    (defined(__GNUC__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)))
#  define BAIRRTT_SUPPRESS_DIAGNOSTIC 1
#endif

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS 1
#ifdef BAIRRTT_SUPPRESS_DIAGNOSTIC
#  ifdef __clang__
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wunknown-pragmas"
#    pragma clang diagnostic ignored "-Wunused-variable"
#    pragma clang diagnostic ignored "-Wunused-but-set-variable"
#    pragma clang diagnostic ignored "-Wunused-parameter"
#    pragma clang diagnostic ignored "-Wunused-local-typedef"
#    pragma clang diagnostic ignored "-Wunused-function"
#    pragma clang diagnostic ignored "-Wsign-compare"
#    pragma clang diagnostic ignored "-Wignored-qualifiers"
#    pragma clang diagnostic ignored "-Wshorten-64-to-32"
#    pragma clang diagnostic ignored "-Wmismatched-tags"
#    if __has_warning("-Wdeprecated-copy")
#      pragma clang diagnostic ignored "-Wdeprecated-copy"
#    endif
#  else
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wunknown-pragmas"
#    pragma GCC diagnostic ignored "-Wunused-variable"
#    pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#    pragma GCC diagnostic ignored "-Wunused-parameter"
#    pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#    pragma GCC diagnostic ignored "-Wunused-function"
#    pragma GCC diagnostic ignored "-Wsign-compare"
#    pragma GCC diagnostic ignored "-Wignored-qualifiers"
#    if __GNUC__ >= 6
#      pragma GCC diagnostic ignored "-Wignored-attributes"
#    endif
#    if __GNUC__ >= 9
#      pragma GCC diagnostic ignored "-Wdeprecated-copy"
#    endif
#  endif
#endif

#include <RcppEigen.h>

#include <walnuts/adaptive_walnuts.hpp>
#include <walnuts/config.hpp>
#include <walnuts/walnuts.hpp>

#ifdef BAIRRTT_SUPPRESS_DIAGNOSTIC
#  ifdef __clang__
#    pragma clang diagnostic pop
#  else
#    pragma GCC diagnostic pop
#  endif
#endif

using Eigen::MatrixXd;
using Eigen::VectorXd;

// numerically stable log(1 / (1 + exp(-x)))
inline double log_sigmoid(double x) {
  if (x >= 0.0)
    return -std::log1p(std::exp(-x));
  else
    return x - std::log1p(std::exp(x));
}

// ---- TARGET DENSITY ---------------------------------------------------------
// Model:
//   eta_ij  = alpha_i * (theta_j - beta_i)
//   y_ij    ~ Bernoulli(plogis(eta_ij))
//   alpha_i ~ Exponential(1)   (positive; sampled as alpha_i = exp(alpha_raw_i))
//   beta_i  ~ Normal(0, beta_sd)
//
// Position vector `par` (length 2 * n_items) is c(alpha_raw[1..I], beta[1..I]).
struct IrtLogpGrad {
  MatrixXd Y;        // n_persons x n_items, 0/1 (owned copy)
  VectorXd theta;    // length n_persons; refreshed each Gibbs scan
  double beta_sd;

  IrtLogpGrad() : beta_sd(10.0) {}
  IrtLogpGrad(const MatrixXd& Y_, const VectorXd& theta_, double beta_sd_)
      : Y(Y_), theta(theta_), beta_sd(beta_sd_) {}

  // matches walnuts::LogpGrad
  void operator()(const VectorXd& par, double& lp, VectorXd& grad) const {
    const int n_persons = static_cast<int>(Y.rows());
    const int n_items   = static_cast<int>(Y.cols());
    const double inv_beta_var = 1.0 / (beta_sd * beta_sd);

    lp = 0.0;
    if (grad.size() != 2 * n_items) grad.resize(2 * n_items);

    for (int i = 0; i < n_items; ++i) {
      const double alpha_raw = par(i);
      const double alpha     = std::exp(alpha_raw);
      const double beta      = par(n_items + i);

      double g_alpha   = 0.0;  // d logLik / d alpha_i = sum_j (y - p)(theta_j - beta)
      double sum_resid = 0.0;  // sum_j (y - p), reused for the beta gradient

      for (int j = 0; j < n_persons; ++j) {
        const double d   = theta(j) - beta;
        const double eta = alpha * d;
        const double y   = Y(j, i);

        lp += log_sigmoid((2.0 * y - 1.0) * eta);   // y log p + (1-y) log(1-p)

        const double p     = 1.0 / (1.0 + std::exp(-eta));
        const double resid = y - p;                 // d logLik / d eta
        g_alpha   += resid * d;                     // d eta / d alpha = (theta - beta)
        sum_resid += resid;
      }

      const double g_beta = -alpha * sum_resid;     // d eta / d beta = -alpha

      // priors (constants dropped)
      lp += -alpha + alpha_raw;                      // alpha ~ Exp(1), log scale (+Jacobian)
      lp += -0.5 * beta * beta * inv_beta_var;       // beta ~ Normal(0, beta_sd)

      // gradient
      grad(i)           = g_alpha * alpha + (1.0 - alpha);  // chain rule + d/d alpha_raw of prior
      grad(n_items + i) = g_beta - beta * inv_beta_var;
    }
  }
};

// ---- SAMPLER STATE ----------------------------------------------------------
// Minimal handler satisfying walnuts::ChainHandler (hence SampleHandler too):
// records the most recent draw and the latest tuning so the R layer can read
// them back. WALNUTS pushes draws out through these callbacks; there is no
// public position getter on the samplers.
struct LatestDraw {
  VectorXd position;
  double   lp = 0.0;
  double   step_size = 0.0;
  VectorXd inv_mass;
  void on_sample(const VectorXd& p, double l) { position = p; lp = l; }
  void on_warmup(const VectorXd& p, double l, double s, const VectorXd& m) {
    position = p; lp = l; step_size = s; inv_mass = m;
  }
  void on_warmup_complete(double s, const VectorXd& m) { step_size = s; inv_mass = m; }
};

using IrtAdapter = walnuts::AdaptiveWalnuts<IrtLogpGrad, std::mt19937_64, LatestDraw>;
using IrtWalnuts = walnuts::WalnutsSampler<IrtLogpGrad, std::mt19937_64, LatestDraw>;

// Persistent state for the (alpha, beta) sub-sampler across the whole Gibbs run.
// WALNUTS holds logp_grad, handler, rng, and the configs by reference, so they
// must live here and the IrtSampler must not move -- it lives on the heap behind
// the XPtr.
struct IrtSampler {
  IrtLogpGrad logp_grad;   // data + current theta + beta_sd
  LatestDraw  handler;     // captures each draw / latest tuning out of WALNUTS
  VectorXd    position;    // (alpha_raw, beta): init, then a mirror of handler.position
  double      step_size;   // initial step size (Adam starts here)
  std::mt19937_64 rng;
  bool        frozen = false;
  std::size_t warmup_iter = 0;

  // WALNUTS holds these two by const reference, so they must outlive the adapter.
  walnuts::WarmupConfig   warmup_cfg   = walnuts::WarmupConfigBuilder().build();
  walnuts::SamplingConfig sampling_cfg = walnuts::SamplingConfigBuilder().build();

  std::optional<IrtAdapter> adapter;   // live during warmup
  std::optional<IrtWalnuts> sampler;   // live after freeze()
};
