#' Sample with \code{Stan}
#' 
#' Sample from the posterior distribution by calling \code{\link[rstan:stan]{stan}}. 
#' Check \code{\link[rstan:stan]{stan}} for details of the arguments.
#' 
#' @export
PSSample <- function(file, model_name = "anon_model", model_code = "", fit = NA, 
                     data = list(), pars = NA,
                     chains = 4, iter = 2000, warmup = floor(iter/2), thin = 1, 
                     init = "random", seed = sample.int(.Machine$integer.max, 1), 
                     algorithm = c("NUTS", "HMC", "Fixed_param"), 
                     control = NULL, sample_file = NULL, diagnostic_file = NULL, 
                     save_dso = TRUE, verbose = FALSE, include = TRUE,
                     cores = getOption("mc.cores", 1L),
                     open_progress = interactive() && !isatty(stdout()) &&
                       !identical(Sys.getenv("RSTUDIO"), "1"),
                     ...,
                     boost_lib = NULL, eigen_lib = NULL) {
  rstan::stan(file, model_name = model_name, model_code = model_code, fit = fit,
              data = data, pars = pars, chains = chains, iter = iter, warmup = warmup,
              thin = thin, init = init, seed = seed, algorithm = algorithm, control = control,
              sample_file = sample_file, diagnostic_file = diagnostic_file, save_dso = save_dso,
              verbose = verbose, include = include, cores = cores, open_progress = open_progress,
              ..., boost_lib = boost_lib, eigen_lib = eigen_lib)
}