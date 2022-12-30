#' The family function for survival data
#' 
#' Construct a family object for survival data
#' 
#' @param method the parametric method used for survival data. Can be Cox or AFT.
#' @param link a link function, currently only identity is implemented and used
#' @return A \code{family} object
#' @export
survival <- function(method = "Cox", link = "identity") {
  if (method == 'Cox')
    family = "survival_Cox"
  if (method == "AFT")
    family = "survival_AFT"
  structure(list(family = family, link = link),
            class = "family")
}
