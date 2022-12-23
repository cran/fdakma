#' K-Means Alignment Algorithm
#'
#' This is an implementation of the k-means alignment algorithm originally
#' described in Sangalli et al. (2010), with improvements as proposed in Vantini
#' (2012).
#'
#' @inherit fdacluster::kma
#' @export
#' @importFrom fdacluster kma
#'
#' @examples
#' res <- kma(
#'   fdacluster::simulated30$x,
#'   fdacluster::simulated30$y,
#'   seeds = c(1, 21),
#'   n_clust = 2,
#'   center_method = "medoid",
#'   warping_method = "affine",
#'   dissimilarity_method = "pearson"
#' )
kma <- function(x, y, seeds = NULL, warping_options = c(0.15, 0.15),
                n_clust = 1, maximum_number_of_iterations = 100,
                number_of_threads = 1, parallel_method = 0,
                distance_relative_tolerance = 0.001, use_fence = FALSE,
                check_total_dissimilarity = TRUE, use_verbose = TRUE,
                compute_overall_center = FALSE, warping_method = "affine",
                center_method = "mean", dissimilarity_method = "l2",
                optimizer_method = "bobyqa") {
  fdacluster::kma(x, y, seeds, warping_options, n_clust,
                  maximum_number_of_iterations, number_of_threads,
                  parallel_method, distance_relative_tolerance, use_fence,
                  check_total_dissimilarity, use_verbose,
                  compute_overall_center, warping_method, center_method,
                  dissimilarity_method, optimizer_method)
}
