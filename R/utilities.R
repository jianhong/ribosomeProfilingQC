#' @importFrom txdbmaker makeTxDbFromGFF
prepareTxDb <- function(gtf, param='gtf', ...){
  txdb <- NULL
  if(is(gtf, "TxDb")){
    txdb <- gtf
  }else{
    if(is.character(gtf)){
      gtf <- gtf[1]
      suppressWarnings(suppressMessages(txdb <- makeTxDbFromGFF(gtf, ...)))
    }
  }
  if(!is(txdb, "TxDb")){
    stop("Can not determine annotations from ", param, " parameter.")
  }
  return(txdb)
}
