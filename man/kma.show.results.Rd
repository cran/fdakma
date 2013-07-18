\name{kma.show.results}
\alias{kma.show.results}
\title{
Auxiliary function plotting results of \link{kma} function.
}

\description{
kma.show.results graphically shows the output results of \link{kma} function or one of the outuput results of \link{kma.compare} function.% 

Four or six plots are generated (depending on the presence of \code{y0} and/or \code{y1} among the inputs):%

- Plot of original functions with center (if \code{y0} is given as input in \link{kma}/\link{kma.compare} function).%

- Plot of aligned functions with centers (if \code{y0} is given as input in \link{kma}/\link{kma.compare} function).% 

- Plot of original function first derivatives with center (if \code{y1} is given as input and the chosen \code{similarity.method} concerns function first derivatives).% 

- Plot of aligned function first derivatives with centers (if \code{y1} is given as input and the chosen \code{similarity.method} concerns function first derivatives).%

- Plot of warping functions.%

- Boxplot of similarity/dissimilarity indexes of original and aligned functions.%

}

\usage{
kma.show.results(Result)
}

\arguments{
  \item{Result}{list:
  output of \link{kma} function or one of the outputs of \link{kma.compare} functions.
  }
}

\references{
Sangalli, L.M., Secchi, P., Vantini, S., Vitelli, V., 2010. \emph{"K-mean alignment for curve clustering"}. Computational Statistics and Data Analysis, 54, 1219-1233.
}

\author{
Mirco Patriarca, Laura Sangalli, Piercesare Secchi, Simone Vantini, Valeria Vitelli.
}

\seealso{
  \code{\link{kma}, \link{kma.compare}, \link{kma.similarity}, \link{fdakma}, \link{kma.data}}
}

\examples{
data(kma.data)

x <- kma.data$x # abscissas
y0 <- kma.data$y0 # evaluations of original functions
y1 <- kma.data$y1 # evaluations of original function first derivatives

\dontrun{
# kma function with 2 clusters, allowing affine 
# transformation for the abscissas and considering 
# 'd1.pearson' as similarity.method.
kma.show.results_example1 <- kma (
  x=x, y0=y0, y1=y1, n.clust = 2, 
  warping.method = 'affine', 
  similarity.method = 'd1.pearson',
  center.method = 'k-means', 
  seeds = c(1,21)
)

# Example: kma.show.results shows the results of kma function
kma.show.results(kma.show.results_example1)


# Example using outputs of kma.compare function

# Results of kma function with 3 different 
# numbers of clusters (1,2,3) combined with four alignment
# methods ('NOalignment' by default, 'shift', 'dilation',
# 'affine') and considering 'd1.pearson' as similarity.method.
kma.show.results_example2 <- kma.compare (
  x=x, y0=y0, y1=y1, n.clust = 1:3, 
  warping.method = c('affine'), 
  similarity.method = 'd1.pearson',
  center.method = 'k-means', 
  seeds = c(1,21,30),
  plot.graph=1)

names (kma.show.results_example2)

# To see results for kma function with n.clust=2 
# and warping.method='affine'.
kma.show.results (kma.show.results_example2$Result.affine[[2]])

# Labels assigned to each function for the 
# kma function with n.clust=2 and warping.method='affine'.
kma.show.results_example2$Result.affine[[2]]$labels
}


}
