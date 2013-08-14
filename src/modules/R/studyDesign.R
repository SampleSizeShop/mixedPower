#
# Defines the studyDesign class using
# the S4 class system
#
# Author: Sarah Kreidler
# Created: 8/11/2013
#
#

setClass ("mixedStudyDesign",
          representation ( name = "character",
                           description = "character",
                           X = "matrix",
                           beta = "matrix",
                           Sigma = "matrix",
                           C = "matrix",
                           thetaNull = "matrix"
                           sascall = "character"
                           ),
          prototype ( name ="",
                      description ="",
                      X = matrix(),
                      beta = matrix(),
                      Sigma = matrix(),
                      C = matrix(),
                      thetaNull = matrix(),
                      sascall = ""
          )
)

## validation routines ##

setGeneric("validObject", function(object) standardGeneric("validObject"))
setMethod("validObject", "mixedStudyDesign", function(object){
  # TODO
  return(TRUE);
})

## Getters / Setters ##
setGeneric("name", function(object) standardGeneric("name"))
setMethod("name", "mixedStudyDesign", function(object){return(object@name)})
setGeneric("name<-", function(object, value) standardGeneric("name<-"))
setReplaceMethod("name", "mixedStudyDesign", function(object, value){
  object@name <- value
  validObject(object)
  object
  })

setGeneric("description", function(object) standardGeneric("description"))
setMethod("description", "mixedStudyDesign", function(object){return(object@description)})
setGeneric("description<-", function(object, value) standardGeneric("description<-"))
setReplaceMethod("description", "mixedStudyDesign", function(object, value){
  object@description <- value
  validObject(object)
  object
})

setGeneric("XMatrix", function(object) standardGeneric("XMatrix"))
setMethod("XMatrix", "mixedStudyDesign", function(object){return(object@X)})
setGeneric("XMatrix<-", function(object, value) standardGeneric("XMatrix<-"))
setReplaceMethod("XMatrix", "mixedStudyDesign", function(object, value){
  object@X <- value
  validObject(object)
  object
})

setGeneric("BetaMatrix", function(object) standardGeneric("BetaMatrix"))
setMethod("BetaMatrix", "mixedStudyDesign", function(object){return(object@Beta)})
setGeneric("BetaMatrix<-", function(object, value) standardGeneric("BetaMatrix<-"))
setReplaceMethod("BetaMatrix", "mixedStudyDesign", function(object, value){
  object@Beta <- value
  validObject(object)
  object
})

setGeneric("SigmaMatrix", function(object) standardGeneric("SigmaMatrix"))
setMethod("SigmaMatrix", "mixedStudyDesign", function(object){return(object@Sigma)})
setGeneric("SigmaMatrix<-", function(object, value) standardGeneric("SigmaMatrix<-"))
setReplaceMethod("SigmaMatrix", "mixedStudyDesign", function(object, value){
  object@Sigma <- value
  validObject(object)
  object
})

setGeneric("contrast", function(object) standardGeneric("contrast"))
setMethod("contrast", "mixedStudyDesign", function(object){return(object@C)})
setGeneric("contrast<-", function(object, value) standardGeneric("contrast<-"))
setReplaceMethod("contrast", "mixedStudyDesign", function(object, value){
  object@C <- value
  validObject(object)
  object
})


setGeneric("thetaNullMatrix", function(object) standardGeneric("thetaNullMatrix"))
setMethod("thetaNullMatrix", "mixedStudyDesign", function(object){return(object@ThetaNull)})
setGeneric("thetaNullMatrix<-", function(object, value) standardGeneric("thetaNullMatrix<-"))
setReplaceMethod("thetaNullMatrix", "mixedStudyDesign", function(object, value){
  object@ThetaNull <- value
  validObject(object)
  object
})


setGeneric("sascall", function(object) standardGeneric("sascall"))
setMethod("sascall", "mixedStudyDesign", function(object){return(object@sascall)})
setGeneric("sascall<-", function(object, value) standardGeneric("sascall<-"))
setReplaceMethod("sascall", "mixedStudyDesign", function(object, value){
  object@ThetaNull <- value
  validObject(object)
  object
})
              
test = new("mixedStudyDesign")