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
Package: \tab rMR\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2016-09-07\cr
License: \tab 3\cr
}
}
\author{
Tyler L. Moulton

Maintainer: Tyler L. Moulton <tyler.moulton@mail.mcgill.ca>

}
\references{
Gnaiger, Erich, and Hellmuth Forstner, eds. Polarographic oxygen sensors: Aquatic and physiological applications. Springer Science & Business Media, 2012.

Lumley, Thomas. "biglm: bounded memory linear and generalized linear models". 0.9-1. https://CRAN.R-project.org/package=biglm. 2013.

McDonnell, Laura H., and Lauren J. Chapman. "Effects of thermal increase on aerobic capacity and swim performance in a tropical inland fish." Comparative Biochemistry and Physiology Part A: Molecular & Integrative Physiology 199 (2016): 62-70.

Mechtly, E. A., 1973: The International System of Units, Physical Constants and Conversion Factors. NASA SP-7012, Second Revision, National Aeronautics and Space Administration, Washington, D.C.

Roche, Dominique G., et al. "Finding the best estimates of metabolic rates in a coral reef fish." Journal of Experimental Biology 216.11 (2013): 2103-2110.

Yeager, D. P. and Ultsch, G. R. (1989). Physiological regulation and conformation: a BASIC program for the determination of critical points. Physiological Zoology, 888-907.

http://www.loligosystems.com/

https://en.wikipedia.org/wiki/Barometric_formula

http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
}
\keyword{ package }
\seealso{
\code{\link{biglm}}
}

