\name{calculateApproximatePowerForDesignList}
\alias{calculateApproximatePowerForDesignList}
\title{calculateApproximatePowerForDesignList}
\usage{
calculateApproximatePowerForDesignList(designList, output.data.dir = ".")
}
\arguments{
  \item{designList}{list of pairs of design.glmmFG and glh
  objects}

  \item{output.data.dir}{directory to which data sets are
  written}
}
\value{
design.glmmFG object with a single covariate
}
\description{
Calculate approximate power values for a list of design/glh
pairs.
}
\note{
This function takes about 15-30 minutes to run
}

