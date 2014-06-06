\docType{class}
\name{missingDataPattern}
\alias{missingDataPattern}
\title{missingDataPattern}
\description{
Class describing the pattern of observations for a group of
independent sampling units within a mixed model design.
}
\section{Slots}{

\describe{
\item{\code{group}}{The group (i.e. treatment assignment) for the sampling unit.}

\item{\code{observations}}{The list of indices for which}

\item{\code{size}}{The number of independent sampling units with this pattern
   of observations
An optional \code{character} string specifying the name of the study design.}

\item{\code{designMatrix}}{the design matrix for independent sampling units with
   the described pattern of observations.}
}}
\note{
For theoretical details, please see

Kreidler, S. M., Muller, K. E., & Glueck, D. H. A Power
Approximation for Longitudinal Studies Using the Kenward
and Roger Wald Test in the Linear Mixed Model, In review.
}
\seealso{
design.mixed
}

