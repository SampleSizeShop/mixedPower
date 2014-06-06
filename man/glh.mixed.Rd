\docType{class}
\name{glh.mixed}
\alias{glh.mixed}
\title{glh.mixed}
\description{
Class describing the general linear hypothesis for fixed
effects in the mixed model
}
\section{Slots}{

\describe{
\item{\code{alpha}}{The Type I error rate}

\item{\code{fixedContrast}}{The \code{matrix} of fixed effects contrasts.}

\item{\code{thetaNull}}{The \code{matrix} of between participant contrasts.}

\item{\code{test}}{A \code{character} string indicating the statistical test. At present,
only the Kenward and Roger Wald test is supported with value "Wald, KR ddf".}
}}
\note{
For theoretical details, please see

Kreidler, S. M., Muller, K. E., & Glueck, D. H. A Power
Approximation for Longitudinal Studies Using the Kenward
and Roger Wald Test in the Linear Mixed Model, In review.
}

