#
# R interface to the SAS mixed model simulator
# 
# Author: Sarah Kreidler
# Created: 8/14/2013
#
#

#
# Generate an SAS/IML string representation of an R matrix
#
matrixToIML <- function(name, m) {
  rowDataStr = paste(sapply(1:nrow(m),
                      function(i){paste("\t", paste(m[i,], collapse=" "), 
                                        ifelse(i<nrow(m),",",""))}),collapse="\n")
  imlStr = paste(name, "= {\n", rowDataStr, "\n};\n")
  return(imlStr)
}

#
# Generate SAS code for mixed model simulation
#
generateSASCode.mixedSimulation = function(path) {
  imlCode = paste(
    
    
    )
  write(imlCode, file=path)
  
  
#   /*
#     * Generated SAS code for simulating a linear mixed model
#   *
#     * Author: Sarah Kreidler
#   * Date: 8/12/2013
#   */
#     ```{r setup, echo=FALSE}
#   
#   ```
#   
#   PROC IML SYMSIZE=1000 WORKSIZE=2000;
#   
#   %INCLUDE "calculatePowerKenwardRoger.sxs"/NOSOURCE2;
#   
#   ```{r, echo=FALSE, comment=""}
#   X=diag(4)
#   matrixToIML("X",X)
#   ```
#   
#   X = Xessence@J(10,1,1);
#   
#   C = {1 -1};
#   
#   SigmaS = I(NROW(X));
#   
#   Beta = {
#     1,
#     0
#   };
#   
#   thetaNull = {0};
#   alpha = {0.05};
#   
#   print X;
#   print sigmaS;
#   print Beta;
#   print C;
#   
#   do i = 0 to 3;
#   power = power // calculatePowerKenwardRoger(X, i*Beta, C, SigmaS, thetaNull, alpha);
#   end;
#   print power;
#   
#   quit;
#   
#   
#   /*
#     *
#     * POwerlib equivalent
#   *
#     */
#     PROC IML SYMSIZE=1000 WORKSIZE=2000;
#   %INCLUDE "C:\KeithMullerSoftware\power\Iml\POWERLIB21.IML" /NOSOURCE2;
#   
#   * Define inputs to power program;
#   ALPHA = 0.05;
#   SIGMA = {1};
#   SIGSCAL = {1};
#   
#   ESSENCEX = I(2);
#   REPN = { 10 };
#   
#   BETA = {1 0}`;
#   BETASCAL = DO(0,3,1);
#   C = {1 -1};
#   
#   OPT_OFF= {C U};
#   ROUND = 15;
#   RUN POWER;
#   
#   QUIT;
  
}

#
#
#
# calculateEmpiricalPower
#
# Calculates empirical power for the given study design object.
# Executes a SAS macro from the command line to do so.
# Results are written to the specified data set
#
# Arguments:
# studyDesign - an object describing the study design matrices
# outputDatasetName - name of the output dataset
#
# Returns:
# None, but writes output dataset
#
# throws:
# simpleError on failure
#
calculateEmpiricalPower <- function(studyDesign, 
                                    outputDatasetName="empiricalPower") {
  sascode = generateSASCode(studyDesign, outputDatasetName)
  cat(sascode,file=paste(OUTPUT_DATA_DIR,"gen-",studyDesig$id, "",sep=""),sep="\n")
}

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