.onLoad <- function(libname, pkgname) {
  cli::cli_alert_warning("The {.pkg fdakma} package is only a wrapper around the
  faster and more complete {.pkg fdacluster} package, which exposes the k-means
  alignment algorithm to the user. It is mainly kept alive because it was
  advertised under this name when the seminal paper was published. We strongly
  encourage R users to switch to using the {.pkg fdacluster} from now on.")
}
