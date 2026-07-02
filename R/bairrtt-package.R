#' bairrtt: causal inference with a latent IRT confounder
#'
#' Fits a causal model in which a latent person trait, measured by a
#' two-parameter logistic (2PL) item-response model, confounds a binary
#' treatment and a continuous outcome. Two BART surfaces (from \pkg{dbarts})
#' model the outcome and the treatment assignment as functions of the latent
#' trait; the IRT item parameters are drawn conditional on the trait by the
#' gradient-based WALNUTS sampler, and the trait by per-person Metropolis.
#'
#' The headline entry point is [irt_causal_bart()], which runs the whole Gibbs
#' sampler on an arbitrary item-response matrix. [simulate_irt_causal()]
#' generates data from the model for testing and examples. The WALNUTS
#' item-parameter sampler is exported on its own via [irt_item_sampler()] and
#' friends.
#'
#' @useDynLib bairrtt, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @keywords internal
"_PACKAGE"
