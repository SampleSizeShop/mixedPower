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
                           thetaNull = "matrix",
                           sascall = "character"
                           ),
          prototype ( name ="two sample t-test",
                      description ="A design which tests for the difference of 
                      means between two independent groups",
                      X = diag(2),
                      beta = matrix(c(1,0),nrow=2),
                      Sigma = matrix(c(1)),
                      C = matrix(c(1,-1), nrow=1),
                      thetaNull = matrix(c(0)),
                      sascall = "* define the mixed model fitting macro; 
* this must contain a 'by setID' statement, but can otherwise;
* be defined as needed by the model;
%macro fitMixedModel(datasetName);
  proc mixed data=&datasetName;
		model y = A B / noint solution ddfm=KR;
		by setID;
		contrast \"treatment effect\" A 1 B -1;
	run;
%mend;"
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

setGeneric("betaMatrix", function(object) standardGeneric("betaMatrix"))
setMethod("betaMatrix", "mixedStudyDesign", function(object){return(object@beta)})
setGeneric("betaMatrix<-", function(object, value) standardGeneric("betaMatrix<-"))
setReplaceMethod("betaMatrix", "mixedStudyDesign", function(object, value){
  object@beta <- value
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
setMethod("thetaNullMatrix", "mixedStudyDesign", function(object){return(object@thetaNull)})
setGeneric("thetaNullMatrix<-", function(object, value) standardGeneric("thetaNullMatrix<-"))
setReplaceMethod("thetaNullMatrix", "mixedStudyDesign", function(object, value){
  object@thetaNull <- value
  validObject(object)
  object
})


setGeneric("sascall", function(object) standardGeneric("sascall"))
setMethod("sascall", "mixedStudyDesign", function(object){return(object@sascall)})
setGeneric("sascall<-", function(object, value) standardGeneric("sascall<-"))
setReplaceMethod("sascall", "mixedStudyDesign", function(object, value){
  object@sascall <- value
  validObject(object)
  object
})
              
test = new("mixedStudyDesign")