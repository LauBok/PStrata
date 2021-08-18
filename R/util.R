nice.print <- function(str, width, sep = 2, pad = "-", newline = T){
  n_pad_left <- (width - stringr::str_length(str)) %/% 2 - sep
  n_pad_right <- width - 2 * sep - stringr::str_length(str) - n_pad_left
  return (paste0(
    paste0(rep(pad, n_pad_left), collapse = ''),
    paste0(rep(' ', sep), collapse = ''),
    str,
    paste0(rep(' ', sep), collapse = ''),
    paste0(rep(pad, n_pad_right), collapse = ''),
    ifelse(newline, '\n', '')
  ))
}

get_mean <- function(str){
  eval(parse(text = substitute(str)))
}






