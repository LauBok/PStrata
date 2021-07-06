byRef <- function(..., envir=parent.frame(), inherits=TRUE) {
  cl <- match.call(expand.dots = TRUE)
  cl[c(1, match(c("envir", "inherits"), names(cl), 0L))] <- NULL
  for (x in as.list(cl)) {
    s <- substitute(x)
    sx <- do.call(substitute, list(s), envir=envir)
    dx <- deparse(sx)
    expr <- substitute(assign(dx, s, envir=parent.frame(), inherits=inherits))
    do.call(on.exit, list(expr, add=TRUE), envir=envir)
  }
}

get_symbol_ID <- function(symbol_name, symbol_list){
  a <- which(symbol_list == as.character(symbol_name))
  if (length(a) > 0) return(list(id = a, symbol_list = symbol_list))
  # add to symbol
  symbol_list <- c(symbol_list, as.character(symbol_name))
  return (list(id = length(symbol_list), symbol_list = symbol_list))
}

.convert_call <- function(call_obj, symbol_list){
  for (i in 2:length(call_obj)){
    if (class(call_obj[[i]]) == "name"){
      tmp <- get_symbol_ID(call_obj[[i]], symbol_list)
      symbol_list <- tmp$symbol_list
      call_obj[[i]] <- as.name(paste0("$", tmp$id))
    }
    else if (class(call_obj[[i]]) == "call"){
      tmp <- .convert_call(call_obj[[i]], symbol_list)
      call_obj[[i]] <- tmp$obj
      symbol_list <- tmp$symbol_list
    }
  }
  return (list(obj = call_obj, symbol_list = symbol_list))
}

.format_call <- function(call_obj){
  gsub("I", "", deparse(call_obj, backtick = F))
}

gen_name <- function(class, stratum, desc, name = NA){
  .str <- paste(stratum, collapse = '_')
  .base <- paste0(class, .str, '_', desc)
  .res <- ifelse(is.na(name), .base, paste0(.base, '_', name))
  .res <- gsub("I\\(|\\)", "", .res)
  .res <- gsub('\\^', "pow", .res)
  .res <- gsub('\\*', "mult", .res)
  .res <- gsub('\\+', "add", .res)
  .res <- gsub('\\-', "minus", .res)
  return (gsub("[^0-9a-zA-z_]|\\^", '_', .res))
}

expand_formula <- function(formula, symbol_list, start_number = 1){
  .ft <- terms(formula)
  .tmp <- .convert_call(attr(.ft, "variables"), symbol_list)
  .vars <- sapply(as.list(.tmp$obj)[-1], .format_call)
  symbol_list <- .tmp$symbol_list
  .factors <- attr(.ft, "factors")
  .factor_col <- function(col){ return(paste(.vars[col == 1], collapse = "*"))}
  if (length(.factors) == 0)
    .terms <- c()
  else
    .terms <- apply(.factors, 2, .factor_col)
  has_intercept <- attr(.ft, "intercept") == 1
  if (has_intercept){
    res_base <- paste0("@", start_number)
    start_number <- start_number + 1
  }
  else{
    res_base <- c()
  }
  .num <- paste("@", start_number:(start_number + length(.terms) - 1), sep = "")
  res <- paste0(c(res_base, purrr::map2(.num, .terms, function(x, y) paste0(x, '*', y))), 
                collapse = " + ")
  
  attr(res, "count") <- length(.terms)
  attr(res, "intercept") <- has_intercept
  attr(res, "term.names") <- attr(.ft, "term.labels")
  return (list(formula = res, symbol_list = symbol_list))
}

nice.print <- function(str, width, sep = 2, pad = "-", newline = T){
  n_pad_left <- (width - str_length(str)) %/% 2 - sep
  n_pad_right <- width - 2 * sep - str_length(str) - n_pad_left
  return (paste0(
    paste0(rep(pad, n_pad_left), collapse = ''),
    paste0(rep(' ', sep), collapse = ''),
    str,
    paste0(rep(' ', sep), collapse = ''),
    paste0(rep(pad, n_pad_right), collapse = ''),
    ifelse(newline, '\n', '')
  ))
}

substitute <- function(str){
  paste(
    "function(x, p)",
    gsub("\\$([0-9]*)", "x[\\1]", gsub("@([0-9]*)", "p[\\1]", str)),
    collapse = " "
  )
}

get_mean <- function(str){
  eval(parse(text = substitute(str)))
}






