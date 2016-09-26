\name{Barom.Press}
\alias{Barom.Press}

\title{
Estimate barometric pressure
}
\description{
 Function to estimate barometric pressure based on altitude.
}
\usage{
Barom.Press(elevation.m, units = "atm")
}

\arguments{
  \item{elevation.m}{
Elevation in meters above sea level
}
  \item{units}{
Output units for barometric pressure must be one of: "atm", "kpa", or "mmHg".
}
}
\details{
pretty straightforward
}
\value{
Barometric pressure in specified units
}
\references{
https://en.wikipedia.org/wiki/Atmospheric_pressure
}
\author{
Tyler L. Moulton
}
\note{

}

\seealso{

}
\examples{
bar.pressure1 <- Barom.Press(elevation.m = 1000) # returns "atm"
bar.pressure2 <- Barom.Press(elevation.m = 1000, "kpa")
bar.pressure3 <- Barom.Press(elevation.m = 1000, "mmHg")
}
