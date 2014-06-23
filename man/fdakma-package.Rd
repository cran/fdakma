\name{fdakma-package}
\alias{fdakma-package}
\alias{fdakma}
\docType{package}
\title{Functional Data Analysis: K-Mean Alignment}

\description{
fdakma jointly performs clustering and alignment of a functional dataset (multidimensional or unidimensional functions).
}

\details{
\tabular{ll}{
Package: \tab fdakma\cr
Type: \tab Package\cr
Version: \tab 1.1.1\cr
Date: \tab 2014-06-23\cr
License: \tab GPL-3\cr
}
}

\author{
Mirco Patriarca, Laura Sangalli, Piercesare Secchi, Simone Vantini, Valeria Vitelli, Alice Parodi.
}

\references{
Sangalli, L.M., Secchi, P., Vantini, S., Vitelli, V., 2010. \emph{"K-mean alignment for curve clustering"}. Computational Statistics and Data Analysis, 54, 1219-1233.}

\seealso{
  \code{\link{kma.compare}, \link{kma.similarity}, \link{kma.data}, \link{kma}, \link{kma.show.results} }
}

\examples{
data(kma.data)

x <- kma.data$x # abscissas
y0 <- kma.data$y0 # evaluations of original functions
y1 <- kma.data$y1 # evaluations of original function first derivatives

\dontrun{
# Plot of original functions
matplot(t(x),t(y0), type='l', xlab='x', ylab='orig.func')
title ('Original functions')

# Plot of original function first derivatives
matplot(t(x),t(y1), type='l', xlab='x', ylab='orig.deriv')
title ('Original function first derivatives')


# Example: result of kma function with 2 clusters, 
# allowing affine transformation for the abscissas
# and considering 'd1.pearson' as similarity.method.
fdakma_example <- kma (
  x=x, y0=y0, y1=y1, n.clust = 2, 
  warping.method = 'affine', 
  similarity.method = 'd1.pearson',
  center.method = 'k-means', 
  seeds = c(1,21)
)

kma.show.results(fdakma_example)

names(fdakma_example)

# Labels assigned to each function
fdakma_example$labels

# Total shifts and dilations applied to the original 
# abscissa to obtain the aligned abscissa
fdakma_example$shift
fdakma_example$dilation
}

}

\keyword{Functional Data Analysis}
\keyword{Alignment}
\keyword{Registration}
\keyword{K-Mean Clustering}
\keyword{Time Warping}
