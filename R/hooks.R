# as_draws_array/as_draws_df are posterior's generics; posterior is
# Suggests-only, so the methods register dynamically whenever posterior's
# namespace loads - before or after bairrtt, and without forcing the load here.
# envir locates the generic (posterior's namespace, not bairrtt's):
# registerS3method stores the method in the generic's own S3 table. Pattern
# taken from dbarts's R/hooks.R.
.onLoad <- function(libname, pkgname) {
  register_posterior_methods <- function(...) {
    ns <- asNamespace("posterior")
    registerS3method(
      "as_draws_array",
      "irt_causal_fit",
      as_draws_array.irt_causal_fit,
      envir = ns
    )
    registerS3method(
      "as_draws_df",
      "irt_causal_fit",
      as_draws_df.irt_causal_fit,
      envir = ns
    )
  }
  setHook(packageEvent("posterior", "onLoad"), register_posterior_methods)
  if (isNamespaceLoaded("posterior")) register_posterior_methods()
}
