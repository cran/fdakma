\name{kma.data}
\alias{kma.data}
\docType{data}
\title{Simulated Data}

\description{
kma.data is a functional dataset displaying both amplitude and phase variability.
}

\usage{data(kma.data)}

\format{
List of 3 elements:%

 $x : abscissa values where each function is evaluated%
 
 $y0: evaluations of the original functions on the abscissa grid \code{kma.data$x}%
 
 $y1: evaluations of the original function first derivatives on the abscissa grid \code{kma.data$x}.
 
}

\references{
Sangalli, L.M., Secchi, P., Vantini, S., Vitelli, V., 2010. \emph{"K-mean alignment for curve clustering"}. Computational Statistics and Data Analysis, 54, 1219-1233.}

\seealso{
  \code{\link{kma.compare}, \link{kma.similarity}, \link{fdakma}, \link{kma}, \link{kma.show.results} }
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
}
}
