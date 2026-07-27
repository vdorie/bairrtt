#' Simulate data from the IRT-confounded causal model
#'
#' Draws a data set from the model that [irt_causal_bart()] fits: person
#' abilities `theta ~ N(0, 1)`, 2PL item responses with difficulties
#' `b ~ N(0, 1)` and discriminations `a ~ N(1, 0.2)`, treatment assignment
#' `z ~ Bernoulli(plogis(sum_k theta_k))` (so the traits confound), and a
#' continuous outcome with a constant treatment effect,
#' `y = sum_k prognostic_k * theta_k + ate * z + N(0, 1)`.
#'
#' With `n_traits > 1` each trait gets its own item bank; every item loads on
#' exactly one trait (the disjoint-bank case, which is the identified one --- see
#' "Identification" in `docs/design/multi-trait.md`). The banks are drawn
#' independently and the traits are independent `N(0, 1)`.
#'
#' Two deliberate mismatches with the model [irt_causal_bart()] fits (they are
#' benign, not bugs): the assignment here uses a **logit** link while the fitted
#' assignment BART uses a probit link --- BART is flexible enough to absorb the
#' shape difference, and the assignment model is only a balancing nuisance term;
#' and discriminations are drawn `N(1, 0.2)`, which can technically be negative
#' (negligibly so at the default `sd`), whereas the fitted model places a
#' strictly positive `Exp(1)` prior on discrimination. Keep the discrimination
#' `sd` small so reverse-keyed items do not appear.
#'
#' @param n_persons Number of persons (rows of each response matrix).
#' @param n_items Number of items per trait (columns). Length 1 (the same count
#'   for every trait) or `n_traits`.
#' @param n_traits Number of latent traits `K`. At `K = 1` (the default) the
#'   return value is the single-trait one described below.
#' @param ate True (constant) treatment effect.
#' @param prognostic Coefficient of each trait in the outcome model. Length 1
#'   (recycled) or `n_traits`.
#' @param seed Optional integer seed; if supplied, set before drawing.
#'
#' @return A list with the observed data --- `responses`, `y`, `z` --- and the
#'   ground truth: `theta`, item discriminations `alpha`, item difficulties
#'   `beta`, `prognostic`, and `ate`. At `n_traits = 1`, `responses` is an
#'   `n_persons` x `n_items` 0/1 matrix, `theta` a length-`n_persons` vector,
#'   and `alpha`/`beta` length-`n_items` vectors. At `n_traits > 1`, `responses`
#'   is a length-`K` list of response matrices, `theta` an `n_persons` x `K`
#'   matrix, and `alpha`/`beta` length-`K` lists of vectors --- exactly the
#'   shapes [irt_causal_bart()] takes and returns.
#'
#' @examples
#' sim <- simulate_irt_causal(n_persons = 200, n_items = 30, ate = -0.2, seed = 1)
#' dim(sim$responses)
#' sim$ate
#'
#' # two traits, two item banks
#' sim2 <- simulate_irt_causal(
#'   n_persons = 200, n_items = c(20, 30), n_traits = 2,
#'   ate = -0.2, prognostic = c(1.2, 0.8), seed = 1
#' )
#' sapply(sim2$responses, dim)
#' dim(sim2$theta)
#'
#' @importFrom stats rnorm rbinom plogis
#' @export
simulate_irt_causal <- function(
  n_persons = 1000L,
  n_items = 100L,
  n_traits = 1L,
  ate = -0.2,
  prognostic = 1,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n_persons <- as.integer(n_persons)
  n_traits <- as.integer(n_traits)
  if (length(n_traits) != 1L || is.na(n_traits) || n_traits < 1L) {
    stop("'n_traits' must be a single integer >= 1")
  }
  n_items <- as.integer(n_items)
  if (length(n_items) == 1L) {
    n_items <- rep(n_items, n_traits)
  }
  if (length(n_items) != n_traits) {
    stop("'n_items' must have length 1 or 'n_traits' = ", n_traits)
  }
  prognostic <- as.double(prognostic)
  if (length(prognostic) == 1L) {
    prognostic <- rep(prognostic, n_traits)
  }
  if (length(prognostic) != n_traits) {
    stop("'prognostic' must have length 1 or 'n_traits' = ", n_traits)
  }

  # standard-normal population abilities, one column per trait
  theta <- matrix(rnorm(n_persons * n_traits), n_persons, n_traits)

  beta <- vector("list", n_traits) # item difficulty, per bank
  alpha <- vector("list", n_traits) # item discrimination, per bank
  responses <- vector("list", n_traits)
  for (k in seq_len(n_traits)) {
    p_k <- n_items[k]
    beta[[k]] <- rnorm(p_k)
    alpha[[k]] <- rnorm(p_k, mean = 1, sd = 0.2)
    th_mat <- matrix(theta[, k], n_persons, p_k)
    b_mat <- matrix(beta[[k]], n_persons, p_k, byrow = TRUE)
    a_mat <- matrix(alpha[[k]], n_persons, p_k, byrow = TRUE)
    p_mat <- plogis(a_mat * (th_mat - b_mat))
    responses[[k]] <- matrix(
      rbinom(n_persons * p_k, 1, p_mat),
      n_persons,
      p_k
    )
  }

  # Potential outcomes with a constant treatment effect; treatment depends on
  # every trait, so every trait is a confounder.
  z <- rbinom(n_persons, 1, plogis(as.vector(theta %*% rep(1, n_traits))))
  mu <- as.vector(theta %*% prognostic)
  y1 <- mu + ate + rnorm(n_persons)
  y0 <- mu + rnorm(n_persons)
  y <- ifelse(z == 1, y1, y0)

  # K = 1 unwraps to the single-trait shapes irt_causal_bart() takes directly.
  if (n_traits == 1L) {
    responses <- responses[[1L]]
    theta <- theta[, 1L]
    alpha <- alpha[[1L]]
    beta <- beta[[1L]]
  }

  list(
    responses = responses,
    y = y,
    z = z,
    theta = theta,
    alpha = alpha,
    beta = beta,
    prognostic = prognostic,
    ate = ate
  )
}
