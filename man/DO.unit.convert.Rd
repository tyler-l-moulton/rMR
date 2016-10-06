\name{DO.unit.convert}
\alias{DO.unit.convert}

\title{
Convert between different common units of DO concentration
}
\description{
This function converts between different different units of DO concentration. Takes into account ambient temperature and pressure.
}
\usage{
DO.unit.convert(x, DO.units.in, DO.units.out, bar.units.in,
bar.press, temp.C, bar.units.out = "mmHg")
}

\arguments{
  \item{x}{
    Value or object of class numeric to be converted.
}
  \item{DO.units.in}{
    Units of dissolved oxygen concentration measured, i.e. to be converted from. Must be "mg/L", "PP" (partial pressure), or "pct" (percent). If "PP", the units of partial pressure must be equal to 'bar.units.in'.
}
  \item{DO.units.out}{
    Units of dissolved oxygen concentration desired, i.e. to be converted to. Must be "mg/L", "PP", or "pct".
}
  \item{bar.units.in}{
    Units of barometric pressure of user specified barometric pressure measurement. Must take value of "atm", "kpa", or "mmHg"
}
  \item{bar.press}{
    Ambient barometric pressure measurement
}
  \item{temp.C}{
    Water temperature measured in degrees C
}
  \item{bar.units.out}{
    Used in internal calculation, only visible if output DO.units.out = "PP". Must take value of "atm", "kpa", or "mmHg"
}
}
\details{

}
\value{
Numeric object representing dissolved oxygen concentration in the units specified by 'DO.units.out'. 
}
\references{
Mechtly, E. A., 1973: The International System of Units, Physical Constants and Conversion Factors. NASA SP-7012, Second Revision, National Aeronautics and Space Administration, Washington, D.C.

https://en.wikipedia.org/wiki/Barometric_formula

http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
}
\author{
Tyler L. Moulton
}
\note{
Use this function on entire data columns to convert them to desired units before analysing with functions like 'MR.loops' and 'get.pcrit'.
}



\seealso{
\code{\link{plot}}, 
\code{\link{plot.raw}}, 
\code{\link{cbind}}, 
\code{\link{Eq.Ox.conc}}, 
\code{\link{DO.saturation}}, 
}
\examples{
## on a single value ##

DO.pct<- DO.unit.convert(x= 125.6863, DO.units.in = "PP",
                DO.units.out = "pct",
                bar.units.in = "mmHg", bar.press = 750, temp.C =15)

## Apply to a column in a 'data.frame' class object ## 

## load data ##
data(fishMR)

head(fishMR)
#note that DO data are in mg/L (DO.mgL) and
#that there is an instantaneous temperature column
#(temp.C) and a pressure column (Bar.Pressure.hpa)

DO.pct.col.a <- DO.unit.convert(fishMR$DO.mgL, DO.units.in = "mg/L",
                DO.units.out = "pct",
                bar.units.in = "kpa", bar.press = 101.3,
                temp.C = fishMR$temp.C,
                bar.units.out = "kpa")
DO.pct.col.b<- DO.unit.convert(fishMR$DO.mgL, DO.units.in = "mg/L",
                DO.units.out = "pct",
                bar.units.in = "kpa", bar.press = 101.3,
                temp.C = fishMR$temp.C)
head(DO.pct.col.a)
head(DO.pct.col.b)

# Now with df #

fishMR2 <- cbind(fishMR, DO.pct.col.a)

par(mfrow = c(1,2))
plot.raw(data = fishMR, DO.var.name = "DO.mgL",
          start.time = "2015-07-03 06:15:00",
          end.time = "2015-07-03 08:05:00",
          main = "DO (mg/L) vs time",
          xlab = "time",
          ylab = "DO (mg/L)")
          

plot.raw(data = fishMR2, DO.var.name = "DO.pct.col.a",
          start.time = "2015-07-03 06:15:00",
          end.time = "2015-07-03 08:05:00",
          main = "DO (percent saturation) vs time",
          xlab = "time",
          ylab = "DO (percent saturation)")
}
