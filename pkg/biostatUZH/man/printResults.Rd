\name{printResults}
\alias{printResults}
\title{Confidence interval for a given effect}
\description{Prints the Wald confidence interval and the p-value for a given effect and standard error.}

\usage{
printResults(theta, se.theta, conf.level=0.95, FUN=identity, digits=3)
}

\arguments{
  \item{theta}{the effect size.}
  \item{se.theta}{the standard error.}
  \item{conf.level}{Confidence level. Default is 95\%.}
  \item{FUN}{Transformation to applied to effect size and confidence
  interval. Default is the identity function.}
  \item{digits}{number of digits for the effect.}
}

\value{
  Prints effect and associated Wald confidence interval and p-value.
}

\author{Leonhard Held}

\examples{
printResults(1.7, 2.1)
printResults(1.7, 2.1, FUN=exp)
}
