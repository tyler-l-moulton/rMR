\name{background.resp}
\alias{background.resp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function for determining the background respiration in a respirometer
}
\description{
background.resp takes user-defined start and end times to calculate the backround respiration rate in a respirometer. 
}
\usage{
background.resp(data, DO.var.name,
time.var.name = "std.time", 
start.time, end.time, ...)
}
\arguments{
  \item{data}{
Dataframe to be used, best if formatted by "get.witrox.data" 
}
  \item{DO.var.name}{
Column name of DO variable, formatted as a character string.
}
  \item{time.var.name}{
Column name of time variable as character string. Time column must be formatted as default class for datetime: class = "POSIXct" "POSIXt", strptime format = "\%Y-\%m-\%d \%H:\%M:\%S".
}
  \item{start.time}{
Input start time as character string of strptime format = "\%Y-\%m-\%d \%H:\%M:\%S".
}
  \item{end.time}{
Input endtime as character string of strptime format = "\%Y-\%m-\%d \%H:\%M:\%S".
}
  \item{\dots}{
Passes on arguments to internal functions
}
}
\details{

}
\value{
returns an object of method 'biglm'. The slope of this funtion is the metabolic rate in input units/sec/
}
\references{
Thomas Lumley (2013). biglm: bounded memory linear and generalized linear models. R package version 0.9-1. http://CRAN.R-project.org/package=biglm

}
\author{
Tyler L. Moulton
}
\note{

}


\seealso{
 \code{\link{as.POSIXct}},
 \code{\link{strptime}},
 \code{\link{biglm}}
}
\examples{
data(fish.resp)

bgd.resp <- 
background.resp(d, "DO.mgL", 
           start.time = "2015-07-02 16:05:00",
           end.time = "2015-07-02 16:35:00",
           ylab = "DO (mg/L)", xlab = "time (min)")
}


