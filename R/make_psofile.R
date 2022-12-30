#' Write a pso (\pkg{PSObject}) file
#' 
#' Generate a text file that includes essential information, 
#' which is read by \code{to_stan()} and converted to a stan file for sampling.
#' 
#' @param PSobject an object of class \code{PSobject}
#' @param filename (optional) string. If not \code{NULL}, the pso file will be saved via
#' \code{\link{writeLines}} in a text file named after the string supplied.
#'
#' @return A list of strings, representing a file.  
#' 
make_psofile <- function(PSobject, filename = NULL) {
  
  prior_names <- c("intercept", "coefficient", "sigma", 
                   "alpha", "lambda", "theta")
  
  if (!is.null(filename))
    fileConn <- file(paste0(filename, ".pso"))
    
  lines <- c()
  
  SZDG <- matrix(nrow = 0, ncol = 4)
  G_id <- 0
  
  for (strata_id in 1:PSobject$strata_info$num_strata) {
    S_id <- strata_id - 1
    for (treatment_id in 1:PSobject$strata_info$num_treatment) {
      Z_id <- treatment_id - 1
      D_id <- PSobject$strata_info$strata_matrix[strata_id, treatment_id]
      tmp_G_id <- NULL
      if (PSobject$strata_info$ER_list[strata_id] && nrow(SZDG) > 0) {
        for (i in 1:nrow(SZDG)) {
          if (SZDG[i, 1] == S_id && SZDG[i, 3] == D_id) {
            tmp_G_id <- SZDG[i, 4]
          }
        }
      }
      if (is.null(tmp_G_id)) {
        G_id <- G_id + 1
        tmp_G_id <- G_id
      }
      line <- paste("SZDG", S_id, Z_id, D_id, tmp_G_id, sep = " ")
      SZDG <- rbind(SZDG, c(S_id, Z_id, D_id, tmp_G_id))
      lines <- c(lines, line)
    }
  }
  lines <- c(lines, "")
  
  family_line <- paste("Y", PSobject$Y.family$family, PSobject$Y.family$link)
  lines <- c(lines, family_line, "")
  
  if (length(PSobject$survival.time.points) == 1)
    survtime_line <- paste("T", PSobject$survival.time.points)
  else
    survtime_line <- paste("T", length(PSobject$survival.time.points))
  lines <- c(lines, survtime_line, "")
  
  random_line_S <- paste("random S", length(PSobject$S.formula$random_eff_list))
  random_line_Y <- paste("random Y", length(PSobject$Y.formula$random_eff_list))
  lines <- c(lines, random_line_S, random_line_Y, "")
  
  for(name in prior_names){
    prior_curr <- eval(parse(text = paste0("PSobject$prior_", name)))
    prior_args <- prior_curr$args
    prior_line <- paste("prior", name, prior_curr$name,
                         length(prior_args), 
                         paste0(unlist(prior_args), collapse = " "))

    lines <- c(lines, prior_line)
  }
  
  if (!is.null(filename)) {
    writeLines(lines, fileConn)
    close(fileConn)
  }
  
  meta_data <- list(
    SZDG = SZDG
  )
  return (list(pso_lines = lines, meta_data = meta_data))
}
