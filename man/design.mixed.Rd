\docType{class}
\name{design.mixed}
\alias{design.mixed}
\title{design.mixed}
\description{
Class describing a mixed model study design for use in
power analysis
}
\section{Slots}{

\describe{
\item{\code{name}}{An optional \code{character} string specifying the name of the study design.}

\item{\code{description}}{An optional \code{character} string specifying the brief description of the study design.}

\item{\code{xPatternList}}{list of missing data patterns for all planned participants}

\item{\code{beta}}{The \code{matrix} of regression coefficients.}

\item{\code{Sigma}}{The residual covariance for a complete data case \code{matrix}.}
}}
\note{
For theoretical details, please see

Kreidler, S. M., Muller, K. E., & Glueck, D. H. A Power
Approximation for Longitudinal Studies Using the Kenward
and Roger Wald Test in the Linear Mixed Model, In review.
}

