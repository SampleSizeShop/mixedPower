#
# Generate data sets
#
#
# Generate SAS code for the R matrix
#

#
# Build the stacked mixed model X matrix from the given design
#
getXStacked = function(design) {
  # determine the number of unique treatment groups
  numGroups = length(unique(sapply(design@xPatternList, function(x) {
    return(x@group)
  })))
  
  # get the max number of planned observations 
  maxObs = nrow(design@Sigma)
  
  # form the data frame and column names
  X = data.frame(id=numeric())
  for(grp in 1:numGroups) {
    for(obs in 1:maxObs) {
      X[,paste(c("trt", grp, "_rep", obs),collapse="")] <- numeric(0)
    }
  }
  
  isu = 1
  for(i in 1:length(design@xPatternList)) {
    pattern = design@xPatternList[[i]]
    # get the between/within portions of the design
    Xi = (matrix(diag(numGroups)[pattern@group,], nrow=1) %x%
            diag(maxObs)[pattern@observations,])
    for(j in 1:pattern@size) {
      row = data.frame(cbind(rep(isu, nrow(Xi)), Xi))
      names(row) = names(X)
      X = rbind(X, row)
      isu = isu + 1
    } 
  }
  
  return(X)
}

#
# Build the stacked covariance from the given design
#
getSigmaStacked = function(design) {
  # determine the number of unique treatment groups
  numGroups = length(unique(sapply(design@xPatternList, function(x) {
    return(x@group)
  })))
  
  # get the max number of planned observations 
  maxObs = nrow(design@Sigma)
  
  # form the data frame and column names
  sigmaList = list()
  isu = 1
  for(pattern in design@xPatternList) {
    designWithin = diag(maxObs)[pattern@observations,]
    SigmaD = designWithin %*% design@Sigma %*% t(designWithin)
    for(i in 1:pattern@size) {
      sigmaList[[isu]] = SigmaD
      isu = isu + 1
    }
  }
  
  return(as.matrix(bdiag(sigmaList)))
}

#
# Generate random data sets for the given design
#
simulateData = function(design, replicates=10000, blockSize=1000, 
                        outputDir=".", filePrefix="simData") {
  if (class(design) != "design.mixed") {
    stop("design must have class 'design.mixed'")
  }
  
  X = getXStacked(design)
  SigmaS = getSigmaStacked(design)
  
  # calculate mean
  XB = as.matrix(X[,2:ncol(X)]) %*% design@beta
  # total sample size
  totalN = length(unique(X$id)) 
  
  ### determine the number of sets and the size of each set ###
  
  # when the blocksize does not divide evenly into the
  # total replicates, the size of the last set is the remainder
  lastSetSize = replicates %% blockSize;
  numSets = floor(replicates / blockSize)
  if (lastSetSize > 0) {
    numSets = numSets + 1
  }
  
  ### buil the data set blocks ###
  sapply(1:numSets, function(setNumber) {
    # set start and end iteration numbers for this block
    startIter = ((setNumber-1)*blockSize) + 1
    if (setNumber == numSets && lastSetSize > 0 && lastSetSize != blockSize) {
      numDataSets = lastSetSize
    } else {
      numDataSets = blockSize
    }
    endIter = startIter + numDataSets - 1;
    cat("Generating data sets ", startIter, " to ", endIter, "\n")
    
    # generate numDataSets number of data sets
    dataSet = do.call("rbind", 
                      lapply(1:numDataSets, function(blockNumber) {
                        errorMatrix = mvrnorm(n = 1, 
                                  mu=rep(0, nrow(SigmaS)), 
                                  SigmaS)
                        yData = XB + matrix(errorMatrix)
                        dataBlock = data.frame(setID=((setNumber-1)*blockSize+blockNumber),Y=yData, X=X)
                        return(dataBlock)
                      }))
    # write to the output directory
    filename = paste(outputDir,"/", filePrefix,startIter,"to",endIter,".csv",sep="")
    write.csv(dataSet, row.names=FALSE, file=filename)
    return(filename)
  }) # end sapply
}

#
# Write out SAS/IML code for an R matrix
#
rMatrixToSAS = function(name, m) {
  return(paste(c(paste(c(name, " = {"), collapse=""),
                 paste(apply(m, 1, paste, collapse = " " ), collapse=",\n"),
                 "};"), collapse="\n"))
}

#
# Convert the design and hypothesis to corresponding
# SAS/IML code for simulation
#
genSASCode = function(filename, design, glh, mixedCall,
                      dataPrefix) {
  
  sink(filename)
  # header information
  cat("* ! generated SAS code, do not edit ! ;\n")
  cat("* design: ", design@name, ";\n\n")
  cat("%include \"common.sas\";\n")
  
  # write the proc mixed call macro
  cat("* define the mixed model fitting macro;\n")
  cat("* this must contain a 'by setID' statement, but can otherwise;\n")
  cat("* be defined as needed by the model;\n")
  cat("%macro mixedCall(datasetName);\n")
  cat(mixedCall)
  cat("%mend;\n")
  
  # open PROC IML
  cat("* Generate the data sets;\n")
  cat("proc iml;\n")
  cat("%INCLUDE \"&MODULES_DIR\\simulateMixedModel.sxs\"/NOSOURCE2;\n")
  
  # generate X matrix and column names
  X = getXStacked(design)
  cat("XFullColNames={", 
      paste(sapply(names(X), function(n) { 
        return(paste(c("'", n, "'"), collapse=""))}), collapse=" "), "};\n") 
  cat("XModelColNames = {", 
      paste(sapply(names(X)[2:length(names(X))], function(n) { 
        return(paste(c("'", n, "'"), collapse=""))}), collapse=" "), "};\n") 
  
  
  cat(rMatrixToSAS("X", X), "\n")
  
  # generate Sigma
  cat(rMatrixToSAS("Sigma", getSigmaStacked(design)), "\n")
  
  # generate the beta matrix
  cat(rMatrixToSAS("beta", design@beta), "\n")
  
  # generate the fixed effects contrast
  cat(rMatrixToSAS("C", glh@fixedContrast), "\n")
  
  # generate theta null
  cat(rMatrixToSAS("thetaNull", glh@thetaNull), "\n")
  
  # generate the type I error
  cat("alpha = ", glh@alpha, ";\n")
  
  # 
  cat("simlib= \"outData\";\n")
  cat("simprefix = \"", dataPrefix, "\";\n")
  
  
  # generate the simulation call
  cat("call calculateEmpiricalPowerConditional(10000, 1000,\n")
  cat("simlib, simprefix, \"mixedCall\", X, XFullColNames, XModelColNames, Beta, SigmaS,
      powerResults);\n")
  cat("print powerResults;")
  cat("quit;")
      
  sink(NULL)      
}

