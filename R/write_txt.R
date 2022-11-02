#' Write a PSObject(.pso) file
#' 
#' Generate a text file that includes essential information, which is read by \code{to_stan()} and converted to a stan file for sampling.
#' 
#' @param S.formula a formula for the principal strata model
#' @param Y.formula a formula for the outcome model
#' @param Y.family a family object, indicating the family of the model
#' @param strata a vector of integers, each integer referring to a stratum that is included in the model
#' @param ER a vector of integers, indicating the strata on which the exclusion restriction assumption is assumed.
#' @param prior_... the prior distribution for the corresponding parameters.
#' @param survival.time.points integer, the number of points to evaluate.
#' @return A text file, readable by \code{to_stan()}.
#' 
#' @seealso \code{PStrata}
write.txt <- function(S.formula, Y.formula, Y.family, strata, ER,
                      prior_intercept = prior_flat(),
                      prior_coefficient = prior_normal(),
                      prior_sigma = prior_inv_gamma(),
                      prior_alpha = prior_inv_gamma(),
                      prior_lambda = prior_inv_gamma(),
                      prior_theta = prior_normal(),
                      survival.time.points = 50,
                      filename = NULL){
  
  count.res.var <- function(formula) {
    count.ast <- function(AST) {
      if (is.name(AST))
        return (1)
      else if (is.call(AST))
        return (sum(sapply(AST[-1], count.ast)))
      else 
        return (0)
    }
    
    LHS <- if(length(formula) == 3) formula[[2]] else NULL
    return (count.ast(LHS))
  }
  
  tostr.strata <- function(stratum, num_of_var) {
    if (is.character(stratum)) {
      stopifnot(length(stratum) == 2 * num_of_var)
      return (stratum)
    }
    bits <- c()
    for (i in 1:(2 * num_of_var)){
      bits <- c(stratum %% 2, bits)
      stratum <- stratum %/% 2
    }
    return (paste0(bits, collapse = ""))
  }
  
  if (is.character(strata))
    strata <- sapply(strata, strtoi, base = 2)
  
  prior_names <- c("intercept", "coefficient", "sigma", 
                   "alpha", "lambda", "theta")
  
  if (!is.null(filename))
    fileConn <- file(paste0(filename, ".pso"))
  
  P <- count.res.var(S.formula) - 1
  lines <- c()
  S_id <- 0
  G_id <- 0
  
  SZDG <- matrix(nrow = 0, ncol = 4)
  
  for(stratum in strata){
    S_bin <- tostr.strata(stratum, P)
    S_bin0 <- substr(S_bin, start=1, stop=P)
    D0 <- strtoi(S_bin0, base = 2)
    S_bin1 <- substr(S_bin, start=P+1, stop=P*2)
    D1 <- strtoi(S_bin1, base = 2)
    
    line0 <- paste0("SZDG ", S_id, " 0 ", D0, " ", G_id)
    SZDG <- rbind(SZDG, c(S_id, 0, D0, G_id))

    if(!stratum %in% ER) 
      G_id <- G_id + 1
    
    line1 <- paste0("SZDG ", S_id, " 1 ", D1, " ", G_id)
    SZDG <- rbind(SZDG, c(S_id, 1, D1, G_id))
    
    lines <- c(lines, line0, line1)
    
    S_id <- S_id + 1
    G_id <- G_id + 1
  }
  lines <- c(lines, "")
  
  family_line <- paste0("Y ", Y.family$family, " ", Y.family$link)
  lines <- c(lines, family_line, "")
  
  survtime_line <- paste0("T ", survival.time.points)
  lines <- c(lines, survtime_line, "")
  
  random_line_S <- paste0("random S ", length(lme4::findbars(S.formula)))
  random_line_Y <- paste0("random Y ", length(lme4::findbars(Y.formula)))
  lines <- c(lines, random_line_S, random_line_Y, "")
  
  for(name in prior_names){
    prior_curr <- eval(parse(text = paste0("prior_", name)))
    prior_args <- prior_curr$args
    prior_line <- paste0("prior ", name, " ", prior_curr$name, " ",
                         length(prior_args), " ",
                         paste0(unlist(prior_args), collapse = " "))

    lines <- c(lines, prior_line)
  }
  
  if (!is.null(filename)) {
    writeLines(lines, fileConn)
    close(fileConn)
  }
  
  meta_data <- list(
    SZDG = SZDG,
    strata = strata
  )
  return (list(pso_lines = lines, meta_data = meta_data))
}
