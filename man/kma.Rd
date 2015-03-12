\name{kma}
\alias{kma}
\title{Clustering and alignment of functional data}

\description{
kma jointly performs clustering and alignment of a functional dataset (multidimensional or unidimensional functions). 
To run kma function with different numbers of clusters and/or different alignment methods see \link{kma.compare}. }

\usage{
kma(x, y0 = NULL, y1 = NULL, n.clust = 1, warping.method = "affine",%
similarity.method = "d1.pearson", center.method = "k-means", seeds = NULL,%
optim.method = "L-BFGS-B", span = 0.15, t.max = 0.1, m.max = 0.1, n.out = NULL,% 
tol = 0.01, fence = TRUE, iter.max = 100, show.iter = 0, nstart=2, return.all=FALSE,%
check.total.similarity=FALSE)
}

\arguments{
  \item{x}{matrix \emph{n.func} X \emph{grid.size} or vector \emph{grid.size}: 
  the abscissa values where each function is evaluated. \emph{n.func}: number of functions in the dataset. \emph{grid.size}: maximal number of abscissa values where each function is evaluated. The abscissa points may be unevenly spaced and they may differ from function to function. \code{x} can also be a vector of length \emph{grid.size}. In this case, \code{x} will be used as abscissa grid for all functions.}

  \item{y0}{matrix \emph{n.func} X \emph{grid.size} or array \emph{n.func} X \emph{grid.size} X \emph{d}: evaluations of the set of original functions on the abscissa grid \code{x}. \emph{n.func}: number of functions in the dataset. \emph{grid.size}: maximal number of abscissa values where each function is evaluated. \emph{d}: (only if the sample is multidimensional) number of function components, i.e. each function is a \emph{d}-dimensional curve. Default value of \code{y0} is \code{NULL}. The parameter \code{y0} must be provided if the chosen \code{similarity.method} concerns original functions.}

  \item{y1}{matrix \emph{n.func} X \emph{grid.size} or array \emph{n.func} X \emph{grid.size} X \emph{d}: evaluations of the set of original functions first derivatives on the abscissa grid \code{x}. Default value of \code{y1} is \code{NULL}. The parameter \code{y1} must be provided if the chosen \code{similarity.method} concerns original function first derivatives.}

  \item{n.clust}{scalar: required number of clusters. Default value is \code{1}. Note that if \code{n.clust=1} kma performs only alignment without clustering.}

  \item{warping.method}{character: type of alignment required. If \code{warping.method='NOalignment'} kma performs only k-mean clustering (without alignment). If \code{warping.method='affine'} kma performs alignment (and possibly clustering) of functions using linear affine transformation as warping functions, i.e., \emph{x.final = dilation*x + shift}. If \code{warping.method='shift'} kma allows only shift, i.e., \emph{x.final = x + shift}. If \code{warping.method='dilation'} kma allows only dilation, i.e., \emph{x.final = dilation*x}. Default value is \code{'affine'}.}

  \item{similarity.method}{character: required similarity measure. Possible choices are: \code{'d0.pearson'}, \code{'d1.pearson'}, \code{'d0.L2'}, \code{'d1.L2'}, \code{'d0.L2.centered'}, \code{'d1.L2.centered'}. Default value is \code{'d1.pearson'}. See \link{kma.similarity} for details.}

  \item{center.method}{character: type of clustering method to be used. Possible choices are: \code{'k-means'} and \code{'k-medoids'}. Default value is \code{'k-means'}.}

 \item{seeds}{vector \emph{max(n.clust)} or matrix \emph{nstart} X \emph{n.clust}: indexes of the functions to be used as initial centers. If it is a matrix, each row contains the indexes of the initial centers of one of the \code{nstart} initializations. In the case where not all the values of \code{seeds} are provided, those not provided are randomly chosen among the \code{n.func} original functions. If \code{seeds=NULL} all the centers are randomly chosen. Default value of \code{seeds} is \code{NULL}}.


  \item{optim.method}{character: optimization method chosen to find the best warping functions at each iteration. Possible choices are: \code{'L-BFGS-B'} and \code{'SANN'}. See \link{optim} function for details. Default method is \code{'L-BFGS-B'}.
}

  \item{span}{scalar: the span to be used for the \link{loess} procedure in the center estimation step when \code{center.method='k-means'}. Default value is 0.15. If \code{center.method='k-medoids'} value of \code{span} is ignored.}

  \item{t.max}{scalar: \code{t.max} controls the maximal allowed shift, at each iteration, in the alignment procedure with respect to the range of curve domains. \code{t.max} must be such that \emph{0<t.max<1} (e.g., \code{t.max=0.1} means that shift is bounded, at each iteration, between \emph{-0.1*range(x)} and \emph{+0.1*range(x)}). Default value is \code{0.1}. If \code{warping.method='dilation'} value of \code{t.max} is ignored.}

  \item{m.max}{scalar: \code{m.max} controls the maximal allowed dilation, at each iteration, in the alignment procedure. \code{m.max} must be such that \emph{0<m.max<1} (e.g., \code{m.max=0.1} means that dilation is bounded, at each iteration, between \emph{1-0.1} and \emph{1+0.1} ). Default value is \code{0.1}. If \code{warping.method='shift'} value of \code{m.max} is ignored.}
  
  \item{n.out}{scalar: the desired length of the abscissa for computation of the similarity indexes and the centers. Default value is \code{round(1.1*grid.size)}.}
  
  \item{tol}{scalar: the algorithm stops when the increment of similarity of each function with respect to the corrispondent center is lower than \code{tol}. Default value is \code{0.01}.}

  \item{fence}{boolean: if \code{fence=TRUE} a control is activated at the end of each iteration. The aim of the control is to avoid shift/dilation outlighers with respect to their computed distributions. If \code{fence=TRUE} the running time can increase considerably. Default value of \code{fence} is TRUE.}

  \item{iter.max}{scalar: maximum number of iterations in the k-mean alignment cycle. Default value is \code{100}.}
  
  \item{show.iter}{boolean: if \code{show.iter=TRUE} kma shows the current iteration of the algorithm. Default value is \code{FALSE}.}
 
 \item{nstart}{scalar: number of initializations with different seeds. Default value is \code{2}. This parameter is used only if \code{center.method} is \code{'k-medoids'}. When \code{center.method = 'k-means'} one initialization is performed.}
	
	\item{return.all}{boolean: if \code{return.all=TRUE} the results of all the \code{nstart} initializations are returned; the output is a list of length \code{nstart}. If \code{return.all=FALSE} only the best result is provided (the one with higher mean similarity if \code{similarity.method} is \code{'d0.pearson'} or\code{'d1.pearson'}, or the one with lower distance if \code{similarity.method} is \code{'d0.L2'}, \code{'d1.L2'}, \code{'d0.L2.centered'} or \code{'d1.L2.centered'}). Default value is \code{FALSE}.}

\item{check.total.similarity}{boolean: if \code{check.total.similarity=TRUE} at each iteration the algorithm checks if there is a decrease of the total similarity and stops. In the affermative case the result obtained in the penultimate iteration is returned. Default value is \code{FALSE}}
}


\value{
The function output is a list containing the following elements:

\item{iterations}{scalar: total number of iterations performed by kma function.}

\item{x}{as input.}

\item{y0}{as input.}

\item{y1}{as input.}

\item{n.clust}{as input.}

\item{warping.method}{as input.}

\item{similarity.method}{as input.}

\item{center.method}{as input.}

\item{x.center.orig}{vector \emph{n.out}: abscissa of the original center.}

\item{y0.center.orig}{matrix \emph{1} X \emph{n.out}: the unique row contains the evaluations of the original function center. If \code{warping.method='k-means'} there are two scenarios: if \code{similarity.method='d0.pearson'} or \code{'d0.L2'} or \code{d0.L2.centered} the original function center is computed via \link{loess} procedure applied to original data; if \code{similarity.method='d1.pearson'} or \code{'d1.L2'} or \code{d1.L2.centered} it is computed by integration of first derivatives center \code{y1.center.orig} (the integration constant is computed minimizing the sum of the weighed L2 distances between the center and the original functions). If \code{warping.method='k-medoids'} the original function center is the medoid of original functions.}

\item{y1.center.orig}{matrix \emph{1} X \emph{n.out}: the unique row contains the evaluations of the original function first derivatives center. If \code{warping.method='k-means'} the original center is computed via \link{loess} procedure applied to original function first derivatives. If \code{warping.method='k-medoids'} the original center is the medoid of original functions.}

\item{similarity.orig}{vector: original similarities between the original functions and the original center.}

\item{x.final}{matrix \emph{n.func} X \emph{grid.size}: aligned abscissas.}

\item{n.clust.final}{scalar: final number of clusters. Note that, when \code{center.method='k.means'}, the parameter \code{n.clust.final} may differ from initial number of clusters (i.e., from \code{n.clust}) if some clusters are found to be empty. In this case a warning message is issued.}

\item{x.centers.final}{vector \emph{n.out}: abscissas of the final function centers and/or of the final function first derivatives centers.}

\item{y0.centers.final}{matrix \emph{n.clust.final} X \emph{n.out}: rows contain the evaluations of the final functions centers. \code{y0.centers.final} is \code{NULL} if \code{y0} is not given as input.}

\item{y1.centers.final}{matrix \emph{n.clust.final} X \emph{n.out}: rows contains the evaluations of the final derivatives centers. \code{y1.centers.final} is \code{NULL} if the chosen similarity measure does not concern function first derivatives.}

\item{labels}{vector: cluster assignments.}

\item{similarity.final}{vector: similarities between each function and the center of the cluster the function is assigned to.}

\item{dilation.list }{list: dilations obtained at each iteration of kma function.}

\item{shift.list }{list: shifts obtained at each iteration of kma function.}

\item{dilation}{vector: dilation applied to the original abscissas \code{x} to obtain the aligned abscissas \code{x.final}.}

\item{shift}{vector: shift applied to the original abscissas \code{x} to obtain the aligned abscissas \code{x.final}.}

}

\author{
Alice Parodi, Mirco Patriarca, Laura Sangalli, Piercesare Secchi, Simone Vantini, Valeria Vitelli.
}

\references{
Sangalli, L.M., Secchi, P., Vantini, S., Vitelli, V., 2010. \emph{"K-mean alignment for curve clustering"}. Computational Statistics and Data Analysis, 54, 1219-1233.

Sangalli, L.M., Secchi, P., Vantini, S., 2014. \emph{"Analysis of AneuRisk65 data: K-mean Alignment"}. Electronic Journal of Statistics, Special Section on "Statistics of Time Warpings and Phase Variations", Vol. 8, No. 2, 1891-1904.
}


\seealso{
  \code{\link{kma.compare}, \link{kma.similarity}, \link{fdakma}, \link{kma.data}, \link{kma.show.results} }
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
kma_example <- kma (
  x=x, y0=y0, y1=y1, n.clust = 2, 
  warping.method = 'affine', 
  similarity.method = 'd1.pearson',
  center.method = 'k-means', 
  seeds = c(1,21)
  )

kma.show.results(kma_example)

names(kma_example)

# Labels assigned to each function
kma_example$labels

# Total shifts and dilations applied to the original 
# abscissa to obtain the aligned abscissa
kma_example$shift
kma_example$dilation
}

}
