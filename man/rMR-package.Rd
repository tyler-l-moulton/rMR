\name{rMR-package}
\alias{rMR-package}
\alias{rMR}
\docType{package}
\title{
Importing data from Loligo Systems softwares, calculating metabolic rates, and critical tensions
}
\description{
This package is written to accompany Loligo Systems respirometry equipment. The package includes a function for loading data output by AutoResp software (get.witrox.data()), functions for calculating metabolic rates over user-specified time intervals, extracting critical points from data using broken stick regressions based on Yeager and Ultsch (1989), easy functions for converting between different units of barometric pressure, etc. Requires "biglm" in order reset RAM memory during computations, allowing for the analysis of large datasets.
}
\details{
\tabular{ll}{
Package: \tab loligoR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-09-07\cr
License: \tab 3\cr
}
}
\author{
Tyler L. Moulton

Maintainer: Tyler L. Moulton <tyler.moulton@mail.mcgill.ca>

}
\references{
Yeager, D. P. and Ultsch, G. R. (1989). Physiological regulation and conformation: a BASIC program for the determination of critical points. Physiological Zoology, 888-907.
}
\keyword{ package }
\seealso{
\code{\link{biglm}}
}

