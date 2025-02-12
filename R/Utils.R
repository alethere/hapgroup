
#Utils ------


#' Utility function for R6 classes
#'
#' @param priv private object to return
#' @param warn.msg warning message if user tries to modify this field
active_access <- function(priv,warn.msg = "Field cannot be modified"){
  eval(substitute(function(value){
    if(missing(value)){
      eval(priv)
    }else{
      warning(warn.msg,call. = FALSE)
    }
  }))
}

#' Utility function to produce alphanumeric codes
#'
#' It takes any set of numbers and returns a 0-padded
#' alphanumeric code with a specificed prefix
#'
#' @param x vector of objects to generate labels for
#' @param prefix character to be used as prefix for the labels
numerate <- function(x, prefix = "V"){
  fmt <- paste0(prefix,"%0",ceiling(log10(length(x) + 1)),"d")
  sprintf(fmt,1:length(x))
}
