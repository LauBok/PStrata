extend <- function(arr, axis = c()) {
  n <- length(dim(arr)) + length(axis)
  new_dim <- numeric(n)
  new_dimname <- list()
  new_dim[axis] <- 1
  new_dim[-axis] <- dim(arr)
  new_dimname[axis] <- NULL
  new_dimname[(1:n)[-axis]] <- dimnames(arr)
  return (array(arr, new_dim, new_dimname))
}

broadcast <- function(arr, to_dim) {
  `%[%` <- function(x, indexList) do.call("[", c(list(x), indexList, drop = F))
  from_dim <- dim(arr)
  multiplier <- to_dim / from_dim
  arr %[% lapply(1:length(multiplier), function(i) rep(1:from_dim[i], times = multiplier[i]))
}
