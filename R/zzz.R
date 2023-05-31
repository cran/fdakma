.onAttach <- function(libname, pkgname) {
  msg <- "The {.pkg {pkgname}} package is only a wrapper around the faster and
  more complete {.pkg fdacluster} package, which exposes the k-means alignment
  algorithm to the user. It is mainly kept alive because it was advertised under
  this name when the seminal paper was published. We strongly encourage R users
  to switch to using the {.pkg fdacluster} package from now on."
  rlang::inform(cli::format_inline(msg), class = "packageStartupMessage")
}
