#' K-Means Alignment Algorithm
#'
#' This is an implementation of the k-means alignment algorithm originally
#' described in Sangalli et al. (2010), with improvements as proposed in Vantini
#' (2012).
#'
#' @inherit fdacluster::fdakmeans
#' @export
#'
#' @examples
#' #----------------------------------
#' # Extracts 15 out of the 30 simulated curves in `simulated30_sub` data set
#' idx <- c(1:5, 11:15, 21:25)
#' x <- fdacluster::simulated30$x[idx, ]
#' y <- fdacluster::simulated30$y[idx, , ]
#'
#' #----------------------------------
#' # Runs a k-means clustering with affine alignment, searching for 2 clusters
#' res <- kma(
#'   x,
#'   y,
#'   n_clusters = 2,
#'   seeds = c(1, 11),
#'   warping_class = "affine",
#'   centroid_type = "medoid",
#'   metric = "pearson"
#' )
kma <- function(x, y = NULL,
                n_clusters = 1L,
                seeds = NULL,
                seeding_strategy = c("kmeans++", "exhaustive-kmeans++", "exhaustive", "hclust"),
                warping_class = c("affine", "dilation", "none", "shift", "srsf"),
                centroid_type = "mean",
                metric = c("l2", "pearson"),
                cluster_on_phase = FALSE,
                use_verbose = TRUE,
                warping_options = c(0.15, 0.15),
                maximum_number_of_iterations = 100L,
                number_of_threads = 1L,
                parallel_method = 0L,
                distance_relative_tolerance = 0.001,
                use_fence = FALSE,
                check_total_dissimilarity = TRUE,
                compute_overall_center = FALSE,
                add_silhouettes = TRUE,
                expand_domain = TRUE) {
  fdacluster::fdakmeans(
    x = x, y = y,
    n_clusters = n_clusters,
    seeds = seeds,
    seeding_strategy = seeding_strategy,
    warping_class = warping_class,
    centroid_type = centroid_type,
    metric = metric,
    cluster_on_phase = cluster_on_phase,
    use_verbose = use_verbose,
    warping_options = warping_options,
    maximum_number_of_iterations = maximum_number_of_iterations,
    number_of_threads = number_of_threads,
    parallel_method = parallel_method,
    distance_relative_tolerance = distance_relative_tolerance,
    use_fence = use_fence,
    check_total_dissimilarity = check_total_dissimilarity,
    compute_overall_center = compute_overall_center,
    add_silhouettes = add_silhouettes,
    expand_domain = expand_domain
  )
}
