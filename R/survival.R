#' The family function for survival data
#'
#' Construct a family object for survival data
#'
#' @param method the parametric method used for survival data. Can be Cox or AFT.
#' @param link a link function, currently only identity is implemented and used
#' @return A \code{family} object
#' @export
survival <- function(method = c("Cox", "AFT"), link = "identity") {
  method <- match.arg(method)
  family <- paste0("survival_", method)
  structure(list(family = family, link = link), class = "family")
}
