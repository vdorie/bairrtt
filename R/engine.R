# WALNUTS item-parameter sampler: exported R surface over the compiled engine.
# The `.irt_*` functions are the Rcpp-generated .Call shims (internal); these
# wrappers add validation, coercion, and documentation.

#' 2PL IRT item-parameter sampler (WALNUTS)
#'
#' Creates a persistent sampler for the item parameters of a two-parameter
#' logistic (2PL) item-response model, conditional on the person abilities
#' `theta`. The item parameters --- discrimination `alpha` (positive, with an
#' `Exponential(1)` prior, sampled on the log scale) and difficulty `beta`
#' (`Normal(0, beta_sd)`) --- are drawn by the gradient-based WALNUTS sampler.
#'
#' The intended lifecycle, holding `theta` fixed within each draw but letting it
#' change between draws, is: create the sampler once, call [irt_warmup()] once
#' per outer warm-up scan (adaptation of step size and mass), [irt_freeze()]
#' once to hand off to fixed tuning, then [irt_draw()] once per kept scan. See
#' [irt_causal_bart()] for the full model that drives it.
#'
#' @param responses Integer/numeric matrix of 0/1 item responses, persons in
#'   rows and items in columns (`n_persons` x `n_items`). Coerced to double.
#'   `NA` marks an unanswered item and is dropped from the item likelihood; an
#'   item with nothing observed is then drawn from its prior.
#' @param alpha Numeric vector of initial discriminations, length `n_items`,
#'   all positive.
#' @param beta Numeric vector of initial difficulties, length `n_items`.
#' @param beta_sd Prior standard deviation for the item difficulties.
#' @param step_size Initial leapfrog step size for adaptation.
#' @param seed Integer seed for the sampler's own random-number generator (it is
#'   independent of R's RNG).
#'
#' @return An object of class `"irt_item_sampler"` wrapping an external pointer
#'   to the compiled sampler, plus the `n_persons` and `n_items` dimensions.
#'
#' @seealso [irt_warmup()], [irt_freeze()], [irt_draw()], [irt_tuning()]
#' @export
irt_item_sampler <- function(
  responses,
  alpha,
  beta,
  beta_sd = 10,
  step_size = 0.1,
  seed = 1L
) {
  responses <- as_response_matrix(responses)
  n_items <- ncol(responses)
  alpha <- as.double(alpha)
  beta <- as.double(beta)
  if (length(alpha) != n_items) {
    stop("'alpha' must have length ncol(responses) = ", n_items)
  }
  if (length(beta) != n_items) {
    stop("'beta' must have length ncol(responses) = ", n_items)
  }
  if (anyNA(alpha) || any(alpha <= 0)) {
    stop("'alpha' must be positive and non-missing")
  }
  if (anyNA(beta)) {
    stop("'beta' must be non-missing")
  }
  if (!is.numeric(beta_sd) || length(beta_sd) != 1L || beta_sd <= 0) {
    stop("'beta_sd' must be a single positive number")
  }
  if (!is.numeric(step_size) || length(step_size) != 1L || step_size <= 0) {
    stop("'step_size' must be a single positive number")
  }

  ptr <- .irt_create_sampler(
    responses,
    alpha,
    beta,
    as.double(beta_sd),
    as.double(step_size),
    as.integer(seed)
  )
  structure(
    list(ptr = ptr, n_persons = nrow(responses), n_items = n_items),
    class = "irt_item_sampler"
  )
}

#' Advance an IRT item-parameter sampler
#'
#' `irt_warmup()` takes one adapting step (step size and mass are updated) and
#' `irt_draw()` one fixed-tuning step; both condition on the supplied `theta`
#' and return the current `c(alpha, beta)` draw on the natural scale (`alpha`
#' first, then `beta`, each of length `n_items`). Call [irt_freeze()] once
#' between the warm-up and sampling phases.
#'
#' Warm-up draws are not valid for inference --- they exist to condition the
#' rest of an outer warm-up scan while the tuning adapts.
#'
#' @param sampler An `"irt_item_sampler"` from [irt_item_sampler()].
#' @param theta Numeric vector of person abilities, length `n_persons`.
#'
#' @return A numeric vector of length `2 * n_items`: `c(alpha, beta)`.
#' @name irt_step
#' @seealso [irt_item_sampler()]
#' @export
irt_warmup <- function(sampler, theta) {
  .irt_draw_dispatch(sampler, theta, .irt_warmup_draw)
}

#' @rdname irt_step
#' @export
irt_draw <- function(sampler, theta) {
  .irt_draw_dispatch(sampler, theta, .irt_draw_sample)
}

#' Freeze an IRT item-parameter sampler
#'
#' Ends adaptation and hands the adapted tuning to a fixed-tuning sampler. Call
#' once, after the last [irt_warmup()] and before the first [irt_draw()]. It is
#' an error to freeze twice. Freezing before any [irt_warmup()] call (i.e. with
#' `irt_tuning(sampler)$warmup_iter == 0`) is allowed but skips adaptation
#' entirely, leaving the sampler at its initial `step_size` and identity mass.
#'
#' @param sampler An `"irt_item_sampler"` from [irt_item_sampler()].
#' @return The `sampler`, invisibly.
#' @seealso [irt_item_sampler()]
#' @export
irt_freeze <- function(sampler) {
  check_sampler(sampler)
  .irt_freeze(sampler$ptr)
  invisible(sampler)
}

#' Inspect an IRT item-parameter sampler's tuning
#'
#' Reports the sampler's current adaptation state. During burn-in the
#' `step_size` and `inv_mass` are still adapting, so this is most useful for
#' watching adaptation progress; after [irt_freeze()] they are fixed at the
#' values the kept draws use.
#'
#' @param sampler An `"irt_item_sampler"` from [irt_item_sampler()].
#' @return A list with `frozen` (logical), `warmup_iter` (count of warm-up
#'   steps taken), `step_size`, and `inv_mass` (the diagonal inverse mass).
#' @seealso [irt_item_sampler()]
#' @export
irt_tuning <- function(sampler) {
  check_sampler(sampler)
  .irt_sampler_tuning(sampler$ptr)
}

#' 2PL IRT log-posterior and gradient for the item parameters
#'
#' Evaluates the (alpha, beta) conditional log-posterior and its exact gradient
#' used by [irt_item_sampler()]. Exposed mainly to finite-difference check the
#' gradient. The position vector is on the sampler's unconstrained scale,
#' `c(log(alpha), beta)`.
#'
#' @param par Numeric vector `c(log(alpha), beta)`, length `2 * n_items`.
#' @param responses 0/1 response matrix (`n_persons` x `n_items`); `NA` cells
#'   are dropped from the likelihood and its gradient.
#' @param theta Numeric vector of person abilities, length `n_persons`.
#' @param beta_sd Prior standard deviation for the item difficulties.
#' @return A list with `value` (the log-posterior) and `gradient`.
#' @export
irt_item_logdensity <- function(par, responses, theta, beta_sd = 10) {
  responses <- as_response_matrix(responses)
  par <- as.double(par)
  theta <- as.double(theta)
  if (length(par) != 2L * ncol(responses)) {
    stop("'par' must have length 2 * ncol(responses)")
  }
  if (length(theta) != nrow(responses)) {
    stop("'theta' must have length nrow(responses)")
  }
  .irt_logdensity_grad(par, responses, theta, as.double(beta_sd))
}

#' @export
print.irt_item_sampler <- function(x, ...) {
  frozen <- tryCatch(
    isTRUE(.irt_sampler_tuning(x$ptr)$frozen),
    error = function(e) NA
  )
  cat(sprintf(
    "<irt_item_sampler: %d persons x %d items, %s>\n",
    x$n_persons,
    x$n_items,
    if (isTRUE(frozen)) "frozen (sampling)" else "adapting (warm-up)"
  ))
  invisible(x)
}

# --- internal helpers --------------------------------------------------------

# Coerce to a plain double 0/1 matrix (Eigen::Map<MatrixXd> needs double).
# NA is kept: a missing cell is dropped from the item likelihood, in the C++
# gradient block and in the theta Metropolis ratio alike. See
# docs/design/missing-responses.md. NaN counts as missing too, so this predicate
# and the C++ std::isnan() agree on exactly one set of cells.
as_response_matrix <- function(responses) {
  if (is.data.frame(responses)) {
    responses <- as.matrix(responses)
  }
  if (!is.matrix(responses)) {
    stop("'responses' must be a matrix or data frame")
  }
  if (!is.numeric(responses)) {
    stop("'responses' must be numeric (0/1/NA)")
  }
  storage.mode(responses) <- "double"
  # The observed mask leads the chain so that a missing cell short-circuits to
  # FALSE rather than NA -- any(NA) is NA, and if (NA) is an error.
  observed <- !is.na(responses)
  if (any(observed & responses != 0 & responses != 1)) {
    stop("'responses' must contain only 0, 1, and NA")
  }
  responses
}

# One item bank per latent trait. A bare matrix (or data frame) is the
# single-trait case; a list of them gives one bank per trait. Item counts may
# differ between banks, but every bank must cover the same persons in the same
# order -- row j is person j throughout.
as_response_banks <- function(responses) {
  banks <- if (!is.list(responses) || is.data.frame(responses)) {
    list(as_response_matrix(responses))
  } else {
    if (length(responses) == 0L) {
      stop("'responses' must contain at least one item bank")
    }
    lapply(responses, as_response_matrix)
  }
  n_persons <- vapply(banks, nrow, integer(1L))
  if (any(n_persons != n_persons[1L])) {
    stop(
      "all item banks in 'responses' must have the same number of rows ",
      "(persons); got ",
      paste(n_persons, collapse = ", ")
    )
  }
  for (k in seq_along(banks)) {
    warn_all_missing(banks[[k]], if (length(banks) == 1L) NULL else k)
  }
  banks
}

# Neither case is an error: both are legitimate in a matrix-sampled design and
# neither breaks the sampler. An item nobody answered has a flat likelihood, so
# its (alpha, beta) are prior draws; a person with nothing observed in a bank
# gets a theta informed only by the prior and by y and z. Both are also a good
# sign the data were assembled wrong. See "Degenerate data" in
# docs/design/missing-responses.md.
warn_all_missing <- function(bank, k = NULL) {
  where <- if (is.null(k)) "" else sprintf(" in item bank %d", k)
  observed <- !is.na(bank)
  empty_items <- which(colSums(observed) == 0L)
  if (length(empty_items) > 0L) {
    warning(sprintf(
      paste0(
        "%d item(s)%s have no observed response (column %s); their 'alpha' ",
        "and 'beta' are draws from the prior."
      ),
      length(empty_items),
      where,
      abbreviate_indices(empty_items)
    ))
  }
  empty_persons <- which(rowSums(observed) == 0L)
  if (length(empty_persons) > 0L) {
    warning(sprintf(
      paste0(
        "%d person(s)%s have no observed response (row %s); their 'theta' is ",
        "informed only by the prior and by 'y' and 'z'."
      ),
      length(empty_persons),
      where,
      abbreviate_indices(empty_persons)
    ))
  }
  invisible(NULL)
}

# The first few indices, so a warning about 900 empty rows stays readable.
abbreviate_indices <- function(i, n_show = 5L) {
  if (length(i) <= n_show) {
    return(paste(i, collapse = ", "))
  }
  paste0(paste(i[seq_len(n_show)], collapse = ", "), ", ...")
}

# Observed covariates as a data frame, or NULL. Kept as a frame rather than a
# matrix so dbarts sees factors as factors; the latent regression takes
# model.matrix() of the same frame.
as_covariate_frame <- function(x, n_persons, n_traits) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.matrix(x) || is.vector(x)) {
    x <- as.data.frame(x)
  }
  if (!is.data.frame(x)) {
    stop("'x' must be a matrix, data frame, or NULL")
  }
  if (nrow(x) != n_persons) {
    stop("'x' must have n_persons = ", n_persons, " rows")
  }
  if (ncol(x) == 0L) {
    return(NULL)
  }
  if (anyNA(x)) {
    stop("'x' must not contain NA")
  }
  # A name collision would silently shadow a trait column in the BART frames.
  reserved <- c("y", "z", "theta", paste0("theta", seq_len(n_traits)))
  clash <- intersect(names(x), reserved)
  if (length(clash) > 0L) {
    stop(
      "'x' must not use the reserved column name(s) ",
      paste(sQuote(clash), collapse = ", ")
    )
  }
  x
}

# The trait columns as a data frame, for dbarts()'s and predict()'s interfaces.
trait_frame <- function(theta, theta_names) {
  out <- as.data.frame(theta)
  names(out) <- theta_names
  out
}

check_sampler <- function(sampler) {
  if (!inherits(sampler, "irt_item_sampler")) {
    stop("'sampler' must be an 'irt_item_sampler' (see irt_item_sampler())")
  }
}

.irt_draw_dispatch <- function(sampler, theta, fn) {
  check_sampler(sampler)
  theta <- as.double(theta)
  if (length(theta) != sampler$n_persons) {
    stop("'theta' must have length n_persons = ", sampler$n_persons)
  }
  fn(sampler$ptr, theta)
}
