#' Class \code{"cvgd"}
#' @description An object of class \code{"cvgd"}
#' represents output of \code{coverageDepth}.
#' @aliases cvgd
#' @rdname cvgd-class
#' @slot coverage \code{"list"}, list of \link[IRanges:AtomicList]{CompressedRleList},
#' specify the coverage of features of each sample.
#' @slot granges \link[GenomicRanges:GRangesList-class]{CompressedGRangesList},
#' specify the features.
#' @import methods
#' @import IRanges
#' @import GenomicRanges
#' @exportClass cvgd
#' @examples
#' cvgd()
setClass(Class = "cvgd",
         representation = representation(
           coverage = "list",
           granges  = "CompressedGRangesList"
         ),
         prototype = prototype(
           coverage = list(RleList(compress = TRUE)),
           granges  = GRangesList(compress = TRUE)
         ),
         validity=function(object){
           sapply(object@coverage, function(.ele){
             if(!is(.ele, "CompressedRleList")){
               return("coverage should be a list of CompressedRleList.")
             }
           })
           n <- lengths(object@coverage)
           n <- unique(n)
           if(length(n)!=1){
             return("Lengths of coverage should be identical for each element.")
           }
           if(n!=length(object@granges)){
             return("Lengths of coverage should be identical to the length of granges.")
           }
           return(TRUE)
         })

#' @rdname cvgd-class
#' @param \dots Each argument in \dots becomes an slot in the new \code{"cvgd"}-class.
#' @export
cvgd <- function(...){
  new("cvgd", ...)
}

#' @rdname cvgd-class
#' @exportMethod `$`
#' @param x cvgd object.
setMethod("$", "cvgd", function(x, name) slot(x, name))
#' @rdname cvgd-class
#' @param name A literal character string or a name (possibly backtick quoted).
#' @param value value to replace.
#' @exportMethod `$<-`
setReplaceMethod("$", "cvgd",
                 function(x, name, value){
                   slot(x, name, check = TRUE) <- value
                   x
                 })


#' @rdname cvgd-class
#' @exportMethod `[[`
#' @param i,j indices specifying elements to extract or replace.
#' @param exact see \link[base]{Extract}
setMethod("[[", "cvgd", function(x, i, j, ..., exact=TRUE) slot(x, i))
#' @rdname cvgd-class
#' @exportMethod `[[<-`
setReplaceMethod("[[", "cvgd",
                 function(x, i, ..., value){
                   slot(x, i, check = TRUE) <- value
                   x
                 })

#' @rdname cvgd-class
#' @exportMethod show
#' @param object cvgd object.
#' @importFrom utils head
#'
setMethod("show", "cvgd", function(object){
  samples <- names(object@coverage)
  cat("cvgd object")
  cat(paste(length(samples), "samples:",
              paste(samples, collapse = ", ")))
  if(length(samples)>0){
    cat(paste("\nhead of coverage for sample\n", samples[1]))
    head(object@coverage[[samples[1]]])
  }
  cat("\nhead of granges\n")
  head(object@granges)
})
