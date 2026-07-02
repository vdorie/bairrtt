#' Simulate data from the IRT-confounded causal model
#'
#' Draws a data set from the model that [irt_causal_bart()] fits: person
#' abilities `theta ~ N(0, 1)`, 2PL item responses with difficulties
#' `b ~ N(0, 1)` and discriminations `a ~ N(1, 0.2)`, treatment assignment
#' `z ~ Bernoulli(plogis(theta))` (so `theta` confounds), and a continuous
#' outcome with a constant treatment effect,
#' `y = prognostic * theta + ate * z + N(0, 1)`.
#'
#' Two deliberate mismatches with the model [irt_causal_bart()] fits (they are
#' benign, not bugs): the assignment here uses a **logit** link
#' (`plogis(theta)`) while the fitted assignment BART uses a probit link --- BART
#' is flexible enough to absorb the shape difference, and the assignment model is
#' only a balancing nuisance term; and discriminations are drawn `N(1, 0.2)`,
#' which can technically be negative (negligibly so at the default `sd`), whereas
#' the fitted model places a strictly positive `Exp(1)` prior on discrimination.
#' Keep the discrimination `sd` small so reverse-keyed items do not appear.
#'
#' @param n_persons Number of persons (rows of the response matrix).
#' @param n_items Number of items (columns).
#' @param ate True (constant) treatment effect.
#' @param prognostic Coefficient of `theta` in the outcome model.
#' @param seed Optional integer seed; if supplied, set before drawing.
#'
#' @return A list with the observed data --- `responses` (`n_persons` x
#'   `n_items` 0/1 matrix), `y`, `z` --- and the ground truth: `theta`, item
#'   discriminations `alpha`, item difficulties `beta`, and `ate`.
#'
#' @examples
#' sim <- simulate_irt_causal(n_persons = 200, n_items = 30, ate = -0.2, seed = 1)
#' dim(sim$responses)
#' sim$ate
#'
#' @importFrom stats rnorm rbinom plogis
#' @export
simulate_irt_causal <- function(n_persons = 1000L, n_items = 100L,
                                ate = -0.2, prognostic = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_persons <- as.integer(n_persons)
  n_items <- as.integer(n_items)

  theta <- rnorm(n_persons)               # standard-normal population ability
  beta <- rnorm(n_items)                  # item difficulty
  alpha <- rnorm(n_items, mean = 1, sd = 0.2)  # item discrimination

  th_mat <- matrix(theta, n_persons, n_items)
  b_mat <- matrix(beta, n_persons, n_items, byrow = TRUE)
  a_mat <- matrix(alpha, n_persons, n_items, byrow = TRUE)
  p_mat <- plogis(a_mat * (th_mat - b_mat))
  responses <- matrix(rbinom(n_persons * n_items, 1, p_mat), n_persons, n_items)

  # Potential outcomes with a constant treatment effect; treatment depends on
  # theta, so theta is a confounder.
  z <- rbinom(n_persons, 1, plogis(theta))
  y1 <- prognostic * theta + ate + rnorm(n_persons)
  y0 <- prognostic * theta + rnorm(n_persons)
  y <- ifelse(z == 1, y1, y0)

  list(responses = responses, y = y, z = z,
       theta = theta, alpha = alpha, beta = beta, ate = ate)
}
