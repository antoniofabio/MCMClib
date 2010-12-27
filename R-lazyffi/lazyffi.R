library(Rffi)

make.fun <- function(symbol) {
  force(symbol)
  dictionary <- list(double = doubleType, integer = sint32Type,
                     character = stringType)
  function(..., RETURN=voidType) {
    args <- list(...)
    modes <- sapply(args, storage.mode)
    ffiTypes <- vector(length(modes), mode="list")
    if(length(ffiTypes) > 0) {
      inDict <- modes %in% names(dictionary)
      ffiTypes[inDict] <- dictionary[modes[inDict]]
      ffiTypes[!inDict] <- pointerType
    }
    cif <- prepCIF(RETURN, ffiTypes)
    ans <- callCIF(cif, symbol, ...)
    if(identical(RETURN, voidType)) {
      return(invisible(NULL))
    } else {
      if(is.list(ans) && ("value" %in% names(ans))) {
        return(ans$value)
      } else {
        return(ans)
      }
    }
  }
}

dyn.import <- function(symbs, envir=parent.frame()) {
  for(symbol in symbs) {
    assign(symbol, make.fun(symbol), envir=envir)
  }
}
