\name{names.vec}
\alias{names.vec}
\docType{data}
\title{
header for fishMR data frame
}
\description{
Vector of column names for a data frame imported with 'get.witrox.data' function from AutoResp}
\usage{data("names.vec")}
\format{
  The format is:
 chr [1:6] "Date.time" "timestamp" "Bar.Pressure.hpa" "Phase" ...
}
\details{
Example header for a text file from Auto Resp.
}

\examples{
data(names.vec)

str(names.vec) 
}
\keyword{datasets}
