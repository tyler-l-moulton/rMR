\name{DO.saturation}
\alias{DO.saturation}

\title{
  Calculate oxygen saturation of water
}
\description{
Calculate the percent saturation of oxygen in water given external temperature, barometric pressure, and recorded DO concentration in mg/L.
}
\usage{
DO.saturation(DO.mgl, temp.C,
elevation.m = NULL, bar.press = NULL,
bar.units = NULL)
}
\arguments{
  \item{DO.mgl}{
Recorded DO concentration in mg/L.
}
  \item{temp.C}{
Temperature in degrees C.
}
  \item{elevation.m}{
Elevation in meters above sea level. EITHER 'elevation.m' or 'bar.press' must be specified.
}
  \item{bar.press}{
Barometric pressure in user defined units (bar.units)--defaults to 'NULL'. EITHER 'elevation.m' or 'bar.press' must be specified.
}
  \item{bar.units}{
Units of barometric pressure, defaults to 'NULL'. must be 'atm', 'kpa' or "mmHg"
}
}
\details{

}
\value{
Returns numeric value of dissolved oxygen saturation.
}
\references{

}
\author{
Tyler L. Moulton
}
\note{
}


\seealso{
\code{\link{Eq.Ox.conc}},
\code{\link{DO.unit.convert}}, 
}
\examples{
DO.sat1 <- DO.saturation(DO.mgl = 5.5,
temp.C = 20, elevation.m = 1000)

DO.sat2 <- DO.saturation(DO.mgl = 5.5,
temp.C = 20, bar.press = 674.1, bar.units = "mmHg")

DO.sat1
DO.sat2

# Will ya look at that...
}

