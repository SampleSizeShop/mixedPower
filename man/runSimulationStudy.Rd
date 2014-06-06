\name{runSimulationStudy}
\alias{runSimulationStudy}
\title{runSimulationStudy}
\usage{
runSimulationStudy(study.seed = 5896, study.data.dir = ".",
  study.figures.dir = ".", study.runEmpirical = TRUE,
  study.sasEmpirical = FALSE)
}
\arguments{
  \item{study.seed}{the random number seed (defaults to
  7634)}

  \item{study.data.dir}{the directory into which data files
  are written (defaults to current working directory)}

  \item{study.figures.dir}{the directory into which pdf
  figures are written (defaults to the current working
  directory)}

  \item{study.runEmpirical}{if true, empirical power values
  will be recalculated. If false, existing empirical power
  values will be loaded from the R package.}
}
\description{
This function reproduces the simulation study results for
the manuscript:\cr Kreidler, S. M., Muller, K. E., &
Glueck, D. H. A Power Approximation for Longitudinal
Studies Using the Kenward and Roger Wald Test in the Linear
Mixed Model, In review.
}
\note{
The empirical power calculations may take several hours to
run
}

