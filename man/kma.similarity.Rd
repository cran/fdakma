\name{kma.similarity}
\alias{kma.similarity}
\title{Similarity/dissimilarity index between two functions}

\description{
kma.similarity computes a similarity/dissimilarity measure between two functions \eqn{f} and \eqn{g}. Users can choose among different types of measures.
}

\usage{
kma.similarity(x.f = NULL, y0.f = NULL, y1.f = NULL,%
x.g = NULL, y0.g = NULL, y1.g = NULL, similarity.method, unif.grid = TRUE)
}

\arguments{

  \item{x.f}{vector \emph{length.f}: abscissa grid where function \eqn{f}  and his first derivatives \eqn{f'} is evaluated. \emph{length.f}: numbrt of abscissa values where \eqn{f} is evaluated. \code{x.f} must always be provided.}
  
  \item{y0.f}{vector \emph{length.f} or matrix \emph{length.f} X \emph{d}: evaluations of function \eqn{f} on the abscissa grid \code{x.f}. \emph{length.f}: number of abscissa values where \eqn{f} is evaluated. \emph{d} (only if \eqn{f} and \eqn{g} are multidimensional) number of function's components, i.e. \eqn{f} is \eqn{d}-dimensional curve.  Default value of \code{y0.f} is \code{NULL}. The vector\code{y0.f} must be provided if the chosen \code{similarity.method} concerns original functions.}
  
  \item{y1.f}{vector \emph{length.f} or matrix \emph{length.f} X \emph{d}: evaluations of \eqn{f} first derivative, i.e., \eqn{f'}, on the abscissa grid \code{x.f}. Default value of \code{y1.f} is \code{NULL}. The vector \code{y1.f} must be provided if the chosen \code{similarity.method} concerns function first derivatives.}
  
  \item{x.g}{vector \emph{length.g}: abscissa grid where function \eqn{g}  and his first derivatives \eqn{g'} is evaluated. \emph{length.g}: numbrt of abscissa values where \eqn{g} is evaluated. \code{x.g} must always be provided.}
  
  \item{y0.g}{vector  \emph{length.g} or matrix \emph{length.g} X \emph{d}: evaluations of function \eqn{g} on the abscissa grid \code{x.g}.  \emph{length.g}: number of abscissa values where \eqn{g} is evaluated. \emph{d} (only if \eqn{f} and \eqn{g} are multidimensional) number of function's components, i.e. \eqn{g} is \eqn{d}-dimensional curve. Default value of \code{y0.g} is \code{NULL}. The vector \code{y0.g} must be provided if the chosen \code{similarity.method} concerns original functions.}
  
  \item{y1.g}{vector \emph{length.g} or matrix \emph{length.g} X \emph{d}: evaluations of \eqn{g} first derivative, i.e., \eqn{g'}, on the abscissa grid \code{x.g}. Default value is of \code{y1.g} \code{NULL}. The vector \code{y1.g} must be provided if the chosen \code{similarity.method} concerns function first derivatives.}
  
  \item{similarity.method}{character: similarity/dissimilarity between \eqn{f} and \eqn{g}. Possible choices are: \code{'d0.pearson'}, \code{'d1.pearson'}, \code{'d0.L2'}, \code{'d1.L2'},%
\code{'d0.L2.centered'}, \code{'d1.L2.centered'}. Default value is \code{'d1.pearson'}. See details.}
  
  \item{unif.grid}{boolean: if equal to \code{TRUE} the similarity measure is computed over an uniform grid built in the intersection domain of the two functions, that is an additional discretization is performed. If equal to \code{FALSE} the additional discretization is not performed, so the functions are supposed to be already defined on the same abscissa grid and the grid is supposed to be fine enough to well compute similarity.
}
}


\details{
We report the list of the currently available similarities/dissimilarities. Note that all norms and inner products are computed over \eqn{D}, that is the intersection of the domains of \eqn{f} and \eqn{g}. \eqn{\overline{f}} and \eqn{\overline{g}} denote the mean value, respectively, of functions \eqn{f} and \eqn{g}.%


1. \code{'d0.pearson'}: this similarity measure is the cosine of the angle between the two functions \eqn{f} and \eqn{g}.

\deqn{
\frac{<f,g>_{L^2}}{\|{f}\|_{L^2} \|{g}\|_{L^2}}}{%
             <f-g>_{L^2}/\|{f}\|_{L^2} \|{g}\|_{L^2}}
 
2. \code{'d1.pearson'}: this similarity measure is the cosine of the angle between the two function derivatives \eqn{f'} and \eqn{g'}.

\deqn{
\frac{<f',g'>_{L^2}}{\|{f'}\|_{L^2} \|{g'}\|_{L^2}}}{%
             <f'-g'>_{L^2}/\|{f'}\|_{L^2} \|{g'}\|_{L^2}}

3. \code{'d0.L2'}: this dissimilarity measure is the L2 distance of the two functions \eqn{f} and \eqn{g} normalized by the length of the common domain \eqn{D}.

\deqn{
\frac{\|{f-g}\|_{L^2}}{|D|}}{%
             \|{f-g}\|_{L^2}/|D|}


4. \code{'d1.L2'}: this dissimilarity measure is the L2 distance of the two function first derivatives \eqn{f'} and \eqn{g'} normalized by the length of the common domain \eqn{D}.

\deqn{
\frac{\|{f'-g'}\|_{L^2}}{|D|}}{%
             \|{f'-g'}\|_{L^2}/|D|}

5. \code{'d0.L2.centered'}: this dissimilarity measure is the L2 distance of \eqn{f-\overline{f}} and \eqn{g-\overline{g}} normalized by the length of the common domain \eqn{D}.

\deqn{
\frac{\|{(f-\overline{f})-(g-\overline{g})}\|_{L^2}}{|D|}}{%
             \|{(f-\overline{f})-(g-\overline{g})}\|_{L^2}/|D|}
             

6. \code{'d1.L2.centered'}: this dissimilarity measure is the L2 distance of \eqn{f'-\overline{f'}} and \eqn{g'-\overline{g'}} normalized by the length of the common domain \eqn{D}.

\deqn{
\frac{\|{(f'-\overline{f'})-(g'-\overline{g'})}\|_{L^2}}{|D|}}{%
             \|{(f'-\overline{f'})-(g'-\overline{g'})}\|_{L^2}/|D|}

For multidimensional functions, if \code{similarity.method='d0.pearson'} or \code{'d1.pearson'} the similarity/dissimilarity measure is computed via the average of the indexes in all directions.%

The coherence properties specified in Sangalli et al. (2010) implies that if \code{similarity.method} is set to \code{'d0.L2'}, \code{'d1.L2'}, \code{'d0.L2.centered'} or \code{'d1.L2.centered'}, value of \code{warping.method} must be \code{'shift'} or \code{'NOalignment'}. If \code{similarity.method} is set to \code{'d0.pearson'} or \code{'d1.pearson'} all values for \code{warping.method} are allowed.%
}

\value{
scalar: similarity/dissimilarity measure between the two functions \eqn{f} and \eqn{g} computed via the similarity/dissimilarity measure specified.
}

\references{
Sangalli, L.M., Secchi, P., Vantini, S., Vitelli, V., 2010. \emph{"K-mean alignment for curve clustering"}, Computational Statistics and Data Analysis, 54, 1219-1233, 2010.
}

\author{
Mirco Patriarca, Laura Sangalli, Piercesare Secchi, Simone Vantini, Valeria Vitelli, Alice Parodi.
}

\seealso{
  \code{\link{kma}, \link{kma.compare}, \link{kma.show.results}, \link{fdakma}, \link{kma.data}}
}

\examples{
data(kma.data)

x.f <- kma.data$x # abscissas of f and f'
x.g <- kma.data$x # abscissas of g and g'

y0.f <- kma.data$y0[1,] # evaluations of f on the abscissa grid x.f
y1.f <- kma.data$y1[1,] # evaluations of f' on the abscissa grid x.f
y0.g <- kma.data$y0[3,] # evaluations of g on the abscissa grid x.g
y1.g <- kma.data$y1[3,] # evaluations of g' on the abscissa grid x.g


\dontrun{
  # Plot of the two functions f and g
  plot(t(x.g),t(y0.g), type='l', xlab='x', ylab='y',, col='red')
  points(t(x.f),t(y0.f), type='l')
  title ('f and g')
  legend('bottomleft', legend=c('f','g'), 
         col=c('black','red'), lty=c(1,1), cex = 0.5)
}

# Example: 'd0.pearson'
kma.similarity (x.f=x.f, y0.f=y0.f, x.g=x.g,  y0.g=y0.g, similarity.method='d0.pearson')

# Example: 'd0.L2'
kma.similarity (x.f=x.f, y0.f=y0.f, x.g=x.g,  y0.g=y0.g, similarity.method='d0.L2')

\dontrun{
  # Plot of the two function first derivatives f' and g'
  plot(t(x.g),t(y1.g), type='l', xlab='x', ylab='y', col='red')
  points(t(x.f),t(y1.f), type='l')
  title ("f' and  g'")
  legend('bottomleft', legend=c("f'","g'"), 
         col=c('black','red'), lty=c(1,1), cex = 0.5)
}

# Example: 'd1.pearson'
kma.similarity (x.f=x.f, y1.f=y1.f, x.g=x.g,  y1.g=y1.g, similarity.method='d1.pearson')

# Example: 'd1.L2'
kma.similarity (x.f=x.f, y1.f=y1.f, x.g=x.g,  y1.g=y1.g, similarity.method='d1.L2')

}

\keyword{Similarity}
