\name{MR.loops}
\alias{MR.loops}

\title{
Calculate metabolic rates from multiple closed respirometry loops
}
\description{
This function calculates the metabolic rates from multiple closed respirometry loops simultaneously. Requires lots of user input, but is easy to manipulate. Returns list of metabolic rates, as well as the average metabolic rate and the standard deviation of the sample of metabolic rates, as well as biglm objects for each section of data used to calculate MRs.
}
\usage{
MR.loops(data, DO.var.name, in.DO.meas = "mg/L",
time.var.name = "std.time", out.DO.meas = "mg/L",
start.idx, stop.idx, system.vol = 1,
background.consumption = NULL,
background.start = NULL, background.end = NULL,
temp.C, elevation.m = NULL, bar.press = NULL, 
bar.units = "atm", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
data must include a time variable in standard as.POSIXct format. eg 
"2016-09-25 15:30:00 EST"
}
  \item{DO.var.name}{
Column name of DO variable as character
}
  \item{time.var.name}{
Column name of time variable (which is in as.POSIXct format) as character. defaults to "std.time" as generated from 'get.witrox.data'
 }
  \item{in.DO.meas}{
Units of DO measurement entered in the DO variable column: must be one of "mg/L" for miligrams/liter, "PP" for partial pressure, "pct" for saturation percent.
}
  \item{out.DO.meas}{
Units of DO measurement returned for metabolic rate: must be one of "mg/L" for miligrams/liter, "PP" for partial pressure (units determined by bar.units), "pct" for saturation percent.
}
  \item{start.idx}{
Character class value or vector matching as.POSIXct object coding for date time. Each element of the vector represents the start time of a new 'loop' for calculation of metabolic rates.
}
  \item{stop.idx}{
Character class value or vector matching as.POSIXct object coding for date time. Each element of the vector represents the stop time of a new 'loop' for calculation of metabolic rates.
}
  \item{system.vol}{
system volume in Liters (defaults to 1 L)
}
  \item{background.consumption}{
Default is NULL. Enter the slope of the background rate of consumption if already known. If not known, enter the start and end times of background respiration rate assay in the background.start and background.end arguments.
}
  \item{background.start}{
Only enter value if background.consumption = NULL. Start time of trial for calculating background consumption if it is in the same dataframe as the rest of your trial. Must be a character class element matching as.POSIXct format.
}
  \item{background.end}{
Only enter value if background.consumption = NULL. Stop time of trial for calculating background consumption if it is in the same dataframe as the rest of your trial. Must be a character class element matching as.POSIXct format.
}
  \item{temp.C}{
Water temperature in degrees C
}
  \item{elevation.m}{
Elevation in m. Only required if bar.press = NULL.
}
  \item{bar.press}{
barometric pressure in units defined by bar.units argument. Only required if elevation.m = NULL
}
  \item{bar.units}{
units of barometric pressure used as input and in output if DO.meas.out = "PP"
}
  \item{\dots}{
arguments passed on to internal functions
}
}
\details{

}
\value{
A list of 2. $MR.summary is a 'data.frame' with 3 columns: $MR.mg.sec (metabolic rate in miligrams/second, this is the same as the slope in each linear model), $sd.slope (standard deviation of slopes calculation), $r.square (adjusted r square value from each model). This second object is a list of biglm objects, each one representing a metabolic "loop" (see McDonnell and Chapman 2016). 
}
\references{
McDonnell, Laura H., and Lauren J. Chapman. "Effects of thermal increase on aerobic capacity and swim performance in a tropical inland fish." Comparative Biochemistry and Physiology Part A: Molecular & Integrative Physiology 199 (2016): 62-70.
}
\author{
Tyler L. Moulton
}
\note{

}

\seealso{
 \code{\link{as.POSIXct}},
 \code{\link{strptime}},
 \code{\link{background.resp}},
 \code{\link{Barom.Press}},
 \code{\link{Eq.Ox.conc}},
 }
\examples{
## load data ##
data(fishMR)

## set vectors of time indices for resp loops##
starts <- c("2015-07-03 01:22:00", "2015-07-03 02:13:00",
            "2015-07-03 03:02:00", "2015-07-03 03:50:00",
            "2015-07-03 04:50:00")

stops <- c("2015-07-03 01:44:00", "2015-07-03 02:35:30",
           "2015-07-03 03:25:00", "2015-07-03 04:16:00",
           "2015-07-03 05:12:00")

metR <- MR.loops(data = fishMR, DO.var.name ="DO.mgL",
         start.idx = starts,
         stop.idx = stops, time.var.name = "std.time",
         temp.C = "temp.C", elevation.m = 1180,
         bar.press = NULL, in.DO.meas = "mg/L", 
         background.start = "2015-07-02 16:05:00",
         background.end = "2015-07-02 16:35:00",
         ylim=c(6,8))

## view output summary ##

metR$MR.summary


## This looks great, but you need to check your start and
## stop vectors, otherwise, you could end up with some
## atrocious loops, e.g.:

starts <- c("2015-07-03 01:15:00", "2015-07-03 02:13:00",
            "2015-07-03 03:02:00", "2015-07-03 03:50:00",
            "2015-07-03 04:50:00")

stops <- c("2015-07-03 01:50:00", "2015-07-03 02:35:30",
           "2015-07-03 03:25:00", "2015-07-03 04:16:00",
           "2015-07-03 05:12:00")

metR <- MR.loops(data = fishMR, DO.var.name ="DO.mgL",
         start.idx = starts,
         stop.idx = stops, time.var.name = "std.time",
         temp.C = "temp.C", elevation.m = 1180,
         bar.press = NULL, in.DO.meas = "mg/L", 
         background.start = "2015-07-02 16:05:00",
         background.end = "2015-07-02 16:35:00",
         ylim=c(6,8))

}
