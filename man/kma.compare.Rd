\name{kma.compare}
\alias{kma.compare}
\title{kma.compare runs \link{kma} with different numbers of clusters and different warping methods.}

\description{
In kma.compare the user can specify multiple values for n.clust and warping.method.
kma.compare runs the K-Mean Alignment algorithm (\link{kma} function) for all couples of specified values of \code{n.clust} and \code{warping.method}.}

\usage{
kma.compare(x, y0 = NULL, y1 = NULL, n.clust = c(1, 2),%
warping.method = c("NOalignment", "shift", "dilation", "affine"),%
similarity.method = "d1.pearson", center.method = "k-means", seeds = NULL,%
optim.method = "L-BFGS-B", span = 0.15, t.max = 0.1, m.max = 0.1, n.out = NULL,%
tol = 0.01, fence = TRUE, iter.max = 100, show.iter = 0, plot.graph = 0, %
nstart = 2, return.all = FALSE)
}

\arguments{
  \item{x}{matrix \emph{n.func} X \emph{grid.size} or vector \emph{grid.size}: 
  the abscissa values where each function is evaluated. \emph{n.func}: number of functions in the dataset. \emph{grid.size}: maximal number of abscissa values where each function is evaluated. The abscissa points may be unevenly spaced and they may differ from function to function. \code{x} can also be a vector of length \emph{grid.size}. In this case, \code{x} will be used as abscissa grid for all functions.}

  \item{y0}{matrix \emph{n.func} X \emph{grid.size} or array \emph{n.func} X \emph{grid.size} X \emph{d}: evaluations of the set of original functions on the abscissa grid \code{x}. \emph{n.func}: number of functions in the dataset. \emph{grid.size}: maximal number of abscissa values where each function is evaluated. \emph{d}: (only if the sample is multidimensional) number of function components, i.e. each function is a \emph{d}-dimensional curve. Default value of \code{y0} is \code{NULL}. The parameter \code{y0} must be provided if the chosen \code{similarity.method} concerns original functions.}

  \item{y1}{matrix \emph{n.func} X \emph{grid.size} or array \emph{n.func} X \emph{grid.size} X \emph{d}: evaluations of the set of original functions first derivatives on the abscissa grid \code{x}. Default value of \code{y1} is \code{NULL}. The parameter \code{y1} must be provided if the chosen \code{similarity.method} concerns original function first derivatives.}

   \item{n.clust}{vector: \code{n.clust} contains the numbers of clusters with which kma.compare runs \link{kma} function. Default value is \code{c(1,2)}. See details.
}
  
  \item{warping.method}{vector: \code{warping.method} contains the types of alignment with which kma.compare runs \link{kma} function. See details.
}

  \item{similarity.method}{character: required similarity measure. Possible choices are: \code{'d0.pearson'}, \code{'d1.pearson'}, \code{'d0.L2'}, \code{'d1.L2'}, \code{'d0.L2.centered'}, \code{'d1.L2.centered'}. Default value is \code{'d1.pearson'}. See \link{kma.similarity} for details.}

  \item{center.method}{character: type of clustering method to be used. Possible choices are: \code{'k-means'} and \code{'k-medoids'}. Default value is \code{'k-means'}.}

  \item{seeds}{vector \emph{max(n.clust) or matrix \emph{nstart} X \emph{n.clust}}: indexes of the functions to be used as initial centers. If it is a matrix, each row contains the indexes of the initial centers of one of the \code{nstart} initializations; if not all the values of \code{seeds} are provided, the ones not introduced are randomly chosen among the \code{n.func} original functions. If \code{seeds=NULL} all the centers are randomly chosen. Default value of \code{seeds} is \code{NULL}.
}

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
  
  \item{plot.graph}{boolean:
  if \code{plot.graph=TRUE}, kma.compare plots a graphic with the means of similarity indexes as ordinate and the number of clusters as abscissa. Default value is \code{FALSE}.
}

	\item{nstart}{scalar: number of initializations with different seeds. Default value is \code{1}.}
	
	\item{return.all}{boolean: if \code{return.all=TRUE} the results of all the \code{nstart} initializations are return; the output is a list of length \code{nstart}. If \code{return.all=FALSE} only the best result is provided (the one with higher mean similarity if \code{similarity.method} is \code{'d0.pearson'} or\code{'d1.pearson'}, or the one with lower similarity if \code{similarity.method} is \code{'d0.L2'}, \code{'d1.L2'}, \code{'d0.L2.centered'} or \code{'d1.L2.centered'}), Default value is \code{FALSE}.}

}

\details{
Example of use: if \code{n.clust=c(1,2,3)} and \code{warping.method=c('shift','affine')}, kma.compare runs \link{kma} function with number of clusters equal to \code{1}, \code{2} and \code{3} using \code{warping.method='shift'} and \code{warping.method='affine'}.
}
  
\value{
The function output is a list containing the following elements:     

  \item{Result.NOalignment}{list of outputs of \link{kma} function with \code{warping.type='NOalignment'}. The sublist \code{Result.NOalignment[[k]]} corresponds to the results when number of clusters is \code{n.clust[k]}. Note that if \code{'NOalignment'} is not chosen as \code{warping.type}, then \code{Result.NOalignment} will be \code{NULL}.}
  
  \item{Result.shift}{list of outputs of \link{kma} function with \code{warping.type='shift'}. The sublist \code{Result.shift[[k]]} corresponds to the results when number of clusters is \code{n.clust[k]}. Note that if \code{'shift'} is not chosen as \code{warping.type}, then \code{Result.shift} will be \code{NULL}.}
  
  \item{Result.dilation}{list of outputs of \link{kma} function with \code{warping.type='dilation'}. The sublist \code{Result.dilation[[k]]} corresponds to the results when number of clusters is \code{n.clust[k]}. Note that if \code{'dilation'} is not chosen as \code{warping.type}, then \code{Result.dilation} will be \code{NULL}.}
  
  \item{Result.affine}{list of outputs of \link{kma} function with \code{warping.type='affine'}. The sublist \code{Result.affine[[k]]} corresponds to the results when number of clusters is \code{n.clust[k]}. Note that if \code{'affine'} is not chosen as \code{warping.type}, then \code{Result.affine} will be \code{NULL}.}
  
  \item{n.clust}{as input.}
  
  \item{mean.similarity.NOalignment}{vector: mean similarity indexes of functions after running \link{kma} function with all elements of \code{n.clust} and \code{warping.type='NOalignment'}. \code{mean.similarity.NOalignment} contains the ordinates of the black curve ("without alignment" in the legend) of the output graphic of the kma.compare function (if \code{plot.graph=1}). }
  
  \item{mean.similarity.shift}{vector: mean similarity indexes of curves after running \link{kma} function with all elements of \code{n.clust} and \code{warping.type='shift'}. \code{mean.similarity.shift} contains the ordinates of the blue curve ("shift" in the legend) of the output graphic of the kma.compare function (if \code{plot.graph=1}). }
  
  \item{mean.similarity.dilation}{vector: mean similarity indexes of curves after running \link{kma} function with all elements of \code{n.clust} and \code{warping.type='dilation'}. \code{mean.similarity.dilation} contains the ordinates of the green curve ("dilation" in the legend) of the output graphic of the kma.compare function (if \code{plot.graph=1}). }
  
  \item{mean.similarity.affine}{vector: mean similarity indexes of curves after running \link{kma} function with all elements of \code{n.clust} and \code{warping.type='affine'}. \code{mean.similarity.affine} contains the ordinates of the orange curve ("affine" in the legend) of the output graphic of the kma.compare function (if \code{plot.graph=1}). }
              
}

\references{
Sangalli, L.M., Secchi, P., Vantini, S., Vitelli, V., 2010. \emph{"K-mean alignment for curve clustering"}. Computational Statistics and Data Analysis, 54, 1219-1233.

Sangalli, L.M., Secchi, P., Vantini, S., 2014. \emph{"Analysis of AneuRisk65 data: K-mean Alignment"}. Electronic Journal of Statistics, Special Section on "Statistics of Time Warpings and Phase Variations", Vol. 8, No. 2, 1891-1904.
}

\author{
Alice Parodi, Mirco Patriarca, Laura Sangalli, Piercesare Secchi, Simone Vantini, Valeria Vitelli.
}


\seealso{
  \code{\link{kma}, \link{kma.similarity}, \link{fdakma}, \link{kma.data}, \link{kma.show.results} }
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


# Example: results of kma function with 3 different 
# numbers of clusters (1,2,3) combined with four alignment
# methods ('NOalignment' by default, 'shift', 'dilation',
# 'affine') and considering 'd1.pearson' as similarity.method.
kma.compare_example <- kma.compare (
  x=x, y0=y0, y1=y1, n.clust = 1:3, 
  warping.method = c('affine'), 
  similarity.method = 'd1.pearson',
  center.method = 'k-means', 
  seeds = c(1,21,30),
  plot.graph=1)

names (kma.compare_example)

# To see results for kma function with n.clust=2 
# and warping.method='affine'.
kma.show.results (kma.compare_example$Result.affine[[2]])


# Labels assigned to each function for the 
# kma function with n.clust=2 and warping.method='affine'.
kma.compare_example$Result.affine[[2]]$labels
}

}
