#'Multiply imputed data set (\code{mids})
#'
#'The \code{mids} object contains a multiply imputed data set. The \code{mids} object is
#'generated by the \code{mice()} and \code{mice.mids()} functions. The \code{mids}
#'class of objects has methods for the following generic functions:
#'\code{print}, \code{summary}, \code{plot}.
#'
#'@section Slots: 
#'  \describe{
#'    \item{\code{.Data}:}{Object of class \code{"list"} containing the 
#'    following slots:}
#'    \item{\code{call}:}{The call that created the object.}
#'    \item{\code{data}:}{A copy of the incomplete data set.}
#'    \item{\code{where}:}{The \code{where} argument of the \code{mice()} function.}
#'    \item{\code{m}:}{The number of imputations.}
#'    \item{\code{nmis}:}{An array containing the number of missing observations per column.}
#'    \item{\code{imp}:}{A list of \code{ncol(data)} components with the generated multiple
#'imputations. Each part of the list is a \code{nmis[j]} by \code{m} matrix of
#'imputed values for variable \code{j}.}
#'    \item{\code{method}:}{A vector of strings of \code{length(ncol(data))} specifying the
#'elementary imputation method per column.}
#'    \item{\code{predictorMatrix}:}{A square matrix of size \code{ncol(data)}
#'containing integers specifying the predictor set.}
#'    \item{\code{visitSequence}:}{The sequence in which columns are visited.}
#'    \item{\code{post}:}{A vector of strings of length \code{ncol(data)} with
#'commands for post-processing}
#'    \item{\code{seed}:}{The seed value of the solution.}
#'    \item{\code{iteration}:}{Last Gibbs sampling iteration number.}
#'    \item{\code{lastSeedValue}:}{The most recent seed value.}
#'    \item{\code{chainMean}:}{A list of \code{m} components. Each component is a
#'\code{length(visitSequence)} by \code{maxit} matrix containing the mean of
#'the generated multiple imputations. The array can be used for monitoring
#'convergence.  Note that observed data are not present in this mean.}
#'    \item{\code{chainVar}:}{A list with similar structure of \code{chainMean},
#'containing the covariances of the imputed values.}
#'    \item{\code{loggedEvents}:}{A \code{data.frame} with six columns containing warnings, corrective actions, and other inside info.}
#'    \item{\code{pad}:}{A list containing various settings of the padded imputation
#'model, i.e. the imputation model after creating dummy variables. Normally,
#'this list is only useful for error checking. List members are \code{pad$data}
#'(data padded with columns for factors), \code{pad$predictorMatrix} (predictor
#'matrix for the padded data), \code{pad$method} (imputation methods applied to
#'the padded data), the vector \code{pad$visitSequence} (the visit sequence
#'applied to the padded data), \code{pad$post} (post-processing commands for
#'padded data) and \code{categories} (a matrix containing descriptive
#'information about the padding operation).}
#'\item{\code{loggedEvents}:}{A matrix with six columns containing a record of
#'automatic removal actions. It is \code{NULL} is no action was made.  At
#'initialization the program does the following three actions: 
#'1. A variable that contains missing values, that is not imputed and that is used as a
#'predictor is removed, 2. a constant variable is removed, and 3. a collinear
#'variable is removed. During iteration, the program does the following
#'actions: 1. one or more variables that are linearly dependent are removed
#'(for categorical data, a 'variable' corresponds to a dummy variable), and 2.
#'proportional odds regression imputation that does not converge and is
#'replaced by \code{polyreg}. Column \code{it} is the iteration number at which
#'the record was added, \code{im} is the imputation number, \code{co} is the
#'column number in the data, \code{dep} is the name of the name of the
#'dependent variable, \code{meth} is the imputation method used, and \code{out}
#'is a (possibly long) character vector with the names of the altered or
#'removed predictors.}
#'}
#'
#' @note Many of the functions of the \code{mice} package do not use the S4 class definitions, 
#' and instead rely on the S3 list equivalent \code{oldClass(obj) <- "mids"}.
#' 
#'@name mids-class
#'@rdname mids-class
#'@aliases mids-class mids
#'@author Stef van Buuren, Karin Groothuis-Oudshoorn, 2000
#'@seealso \code{\link{mice}}, \code{\link[=mira-class]{mira}}, \code{\link[=mipo-class]{mipo}}
#'@references van Buuren S and Groothuis-Oudshoorn K (2011). \code{mice}:
#'Multivariate Imputation by Chained Equations in \code{R}. \emph{Journal of
#'Statistical Software}, \bold{45}(3), 1-67.
#'\url{http://www.jstatsoft.org/v45/i03/}
#'@keywords classes
#'@export
setClass("mids",
         representation(
             call      = "call",
             data      = "data.frame" ,
             where     = "matrix",
             m         = "numeric",
             nmis      = "integer",
             imp       = "list",
             method    = "character",
             predictorMatrix = "matrix",
             visitSequence = "numeric",
             post      = "character",
             seed      = "numeric",
             iteration = "integer",
             lastSeedValue = "integer",
             chainMean = "array",
             chainVar  = "array",
             loggedEvents = "data.frame",
             pad       = "list"),
         contains  = "list"
)
