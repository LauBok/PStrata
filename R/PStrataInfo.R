#' Create an object that defines the principal strata
#' 
#' \code{PStrataInfo} is a class of object that defines all principal strata
#' to be considered, by specifying the potential value of each post-randomization
#' confounding variable under each treatment arm.
#' 
#' Since definition of the principal strata appears fundamental and essential in 
#' principal stratification analyses, the creation of such an object is designed
#' to be user-friendly - various ways are accommodated to create a \code{PStrataInfo}
#' object, some possibly preferable over others under different settings.
#' 
#' 
#' @examples
#' PStrataInfo(strata = c(n = "00*", c = "01", a = "11"))
#' PStrataInfo(
#'   strata = list(n = c(0, 0), c = c(0, 1), a = c(1, 1)), 
#'   ER = c(TRUE, FALSE, FALSE)
#' )
#' PStrataInfo(
#'   strata = list(n = c(0, 0), c = c(0, 1), a = c(1, 1)), 
#'   ER = c("n")
#' )
#' 
#' @param strata a list or a vector defining all principal strata. Details of
#' the syntax are given in 'Details' below.
#' @param ER a vector indicating on which strata exclusion restriction is assumed. Details are
#' given in 'Details' below.
#' 
#' @returns an object of class \code{PSStrataInfo}, which is a list of the following components.
#' \describe{
#'   \item{num_strata}{number of principal strata defined}
#'   \item{num_treatment}{number of treatment arms}
#'   \item{num_postrand_var}{number of post-randomization variables}
#'   \item{max_postrand_level}{integer vector, the biggest number used by each post-randomization variable}
#'   \item{strata_matrix}{integer matrix, each row corresponding to one stratum and
#'   each column corresponding to one treatment arm.
#'   The matrix is designed only for internal use.
#'   }
#'   \item{ER_list}{logical vector, each component corresponding to one stratum, indicating whether 
#'   ER is assumed for the specific stratum}
#'   \item{strata_names}{character vector, the names of all strata}
#' }
#' 
#' @details 
#' There are mainly two ways to easily create a \code{PStrataInfo} object.
#' 
#' \subsection{By string}{
#' To define the principal strata by strings, the \code{strata} argument
#' should receive a named vector, each component being the description of one strata
#' with the name of that strata. The naming does not affect the actual inference,
#' but informative names can be helpful for users to distinguish among strata.
#' 
#' Each stratum is defined by the potential values of the post-randomization confounding
#' variable \eqn{D} under each treatment arm. By convention, assume that the K treatment arms are
#' numbered from 0 to K-1. Then, each stratum is defined by the tuple \eqn{(D(0), \ldots, D(K-1))},
#' which can be written compactly as a string. For example, under binary treatment,
#' the never-takers (i.e. \eqn{D(0) = D(1) = 0}) can be represented by string \code{"00"} and 
#' the compliers (i.e. \eqn{D(0) = 0, D(1) = 1}) can be represented by string \code{"01"}.
#' Note that the value that the post-randomization confounding variable can take is limited between 0 to 9
#' for the string to be parsed correctly. This should be more than enough in most of the applications, and
#' in cases where a number above 10 is needed, please create the \code{PStrataInfo} object by matrix (see below).
#' 
#' When multiple post-randomization confounding variables exist, the string for each confounding variable
#' is concatenated with the symbol "\code{|}". For example, if \eqn{D_0} and \eqn{D_1} are both binary 
#' post-randomization confounding variables, the stratum defined by \eqn{D_0(0) = D_0(1) = 0, D_1(0) = 0, D_1(1) = 1}
#' can be represented by string \code{"00|11"}. The order of these confounding variables should be the same
#' as they appear in the \code{S.formula} parameter in \code{\link{PSObject}}.
#' 
#' A common assumption in practice is the exclusion restriction (ER) assumption, which assumes that 
#' the causal effect of the treatment on the outcome is totally realized through the post-randomization 
#' confounding variables. For example, the ER assumption on the stratum of never-takers can be interpreted
#' as the outcome is identically distributed across the treated and control group, because all causal effect
#' of the treatment is realized through the post-randomization variable, which is the same (0) under both 
#' treatment arms. To assume ER for some stratum, simply put an asterisk "\code{*}" at the end of the string,
#' such as "00*" for the never-taker stratum. \emph{Note that under the context of multiple post-randomization
#' variables, the package treats all such variables as a unity. The outcome is assumed to be identical under
#' different treatment arms only when all post-randomization variables remain the same under these treatment arms.}
#' 
#' Another way to specify the stratum where ER is assumed is to use the \code{ER} argument. It either takes
#' a logical vector of the same length of \code{strata} with \code{TRUE} indicating ER is assumed and \code{FALSE}
#' otherwise, or takes a character vector with the names of all strata where ER is to be assumed upon.
#' When names to the strata are not provided in \code{strata}, the strata can be referred to by their 
#' canonical name, which is the string used to define the stratum with asterisks removed. For example, 
#' the strata \code{"00|11*"} can be referred to with name "00|11".
#' }
#' 
#' \subsection{By matrix}{
#' To define the principal strata by matrices, the \code{strata} argument
#' should receive a named list, each component being a matrix. The number of rows matches the number
#' of post-randomization variables, and the number of columns matches that of possible treatment arms.
#' For any fixed row \eqn{i}, column \eqn{j} stores the potential value of the \eqn{i}-th post-randomization
#' variable under treatment arm \eqn{j}.
#' 
#' When this approach is used, there is no shorthand to specify ER assumption. The \code{ER} argument is 
#' required to do this.
#' }
#' 
#' \bold{Warning:} When ER assumption is specified in both \code{strata} and \code{ER} argument, the shorthand
#' notation for ER in \code{strata} is ignored, and a warning is given regardless of whether the specification
#' given by \code{strata} and \code{ER} actually match.
#' 
#' @export
PStrataInfo <- function (strata, ER = NULL) {
  # shorthand
  strata_without_asterisk <- NULL
  ER_list <- NULL
  if (is.character(strata)) {
    strata_names <- names(strata)
    if (is.null(strata_names)) strata_names <- rep("", length(strata))
    ER_list <- stringr::str_ends(strata, stringr::fixed('*'))
    strata_without_asterisk <- stringr::str_remove(strata, stringr::fixed('*'))
    strata_names <- ifelse(strata_names == "", strata_without_asterisk, strata_names)
    strata_split <- stringr::str_split(strata_without_asterisk, stringr::fixed("|"))
    strata_matrix <- lapply(
                       lapply(
                         lapply(
                           lapply(strata_split, stringr::str_split, ''), 
                           function(x) do.call(cbind, x)),
                         apply, c(1,2), as.integer
                       ), t
                     )
  }
  if (is.list(strata)) {
    strata_matrix <- lapply(strata, function(x) {
      if (!is.matrix(x))
        matrix(x, nrow = 1)
      else
        x
    })
    strata_names <- names(strata)
    stopifnot("All strata must be named if the matrix
              initialization is used." = !is.null(strata_names) || all(strata_names != ""))
  }
  num_strata <- length(strata_matrix)
  num_treatment <- ncol(strata_matrix[[1]])
  num_postrand_var <- nrow(strata_matrix[[1]])
  stacked_strata_matrix <- abind::abind(strata_matrix, along = 3)
  max_postrand_level <- apply(stacked_strata_matrix, 1, max)
  strata_matrix <- apply(stacked_strata_matrix, c(3,2), 
        function(x) t(c(x, 0)) %*% c(1, cumprod(max_postrand_level + 1)))

  if (!is.null(ER)) {
    if (!is.null(ER_list) && any(ER_list)) 
      warning("Shorthand specification of ER is ignored and overwritten, since the ER argument is provided with value.")
    if (is.logical(ER))
      ER_list <- ER
    else {
      ER_list <- rep(FALSE, length(strata))
      for (word in ER) {
        if (word %in% strata_names) {
          ER_list[which(strata_names == word)] <- TRUE
        }
        else if (!is.null(strata_without_asterisk) && word %in% strata_without_asterisk){
          ER_list[which(strata_without_asterisk == word)] <- TRUE
        }
        else {
          warning(paste0("The strata named ", word, " is undefined, thus ignored."))
        }
      }
    }
  }
  
  return (structure(
    list(
      num_strata = num_strata,
      num_treatment = num_treatment,
      num_postrand_var = num_postrand_var,
      max_postrand_level = max_postrand_level,
      strata_matrix = strata_matrix,
      ER_list = ER_list,
      strata_names = strata_names
    ),
    class = "PStrataInfo"
  ))
}
