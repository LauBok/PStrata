get_strata <- function(monotonicity) {
  if (monotonicity == "default")
    return (c("00", "01", "11"))
  else if (monotonicity == "none")
    return (c("00", "01", "10", "11"))
  else if (monotonicity == "strong")
    return (c("00", "01"))
}