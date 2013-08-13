#
# Defines the studyDesign class using
# the S4 class system
#
# Author: Sarah Kreidler
# Created: 8/11/2013
#
#

setClass ("studyDesign",
          representation ( id = "character",
                           description = "character",
                           X = "matrix",
                           B = "matrix",
                           Sigma = "matrix",
                           C = "matrix",
                           U = "matrix",
                           ThetaNull = "matrix"
                           ),
          prototype ( id ="generic Two-sample T-test",
                      description ="A two group study design",
                      X = diag(2),
                      B = as.matrix(c(1,0),nrow=2),
                      Sigma = diag(1),
                      C = as.matrix(c(1,-1), nrow=1),
                      U = diag(1),
                      ThetaNull = as.matrix(c(0)))
          )

setMethod("id", "studyDesign",function)




function (generic, signature, file, external = FALSE, where = topenv(parent.frame())) 
{
  fdef <- getGeneric(generic, where = where)
  if (is.null(fdef)) {
    fdef <- implicitGeneric(generic, where = where)
    if (is.null(fdef)) 
      stop(gettextf("No function definition found for %s", 
                    sQuote(generic)), domain = NA)
  }
  else {
    generic <- fdef@generic
  }
  signature <- matchSignature(signature, fdef)
  if (length(signature) == 0) 
    signature <- "ANY"
  sigNames <- fdef@signature
  length(sigNames) <- length(signature)
  method <- function() {
  }
  formals(method) <- formals(fdef)
  body(method) <- quote({
    stop("Need a definition for the method here")
  })
  methodName <- paste(c(generic, signature), collapse = "_")
  if (missing(file)) 
    file <- paste0(methodName, ".R")
  output <- c(paste0("setMethod(\"", generic, "\","), paste0("    signature(", 
                                                             paste0(sigNames, " = \"", signature, "\"", collapse = ", "), 
                                                             "),"))
  method <- deparse(method)
  if (identical(external, FALSE)) 
    output <- c(output, paste0("    ", method), ")")
  else {
    if (is(external, "character")) 
      methodName <- toString(external)
    method[[1L]] <- paste0("`", methodName, "` <- ", method[[1L]])
    output <- c(method, "", output, paste0("  `", methodName, 
                                           "`)"))
  }
  writeLines(output, file)
  message("Skeleton of method written to ", if (is.character(file)) 
    file
          else "connection")
  invisible(file)
}


              
test = new("studyDesign")