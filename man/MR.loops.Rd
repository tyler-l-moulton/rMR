\name{MR.loops}
\alias{MR.loops}

\title{
Calculate Metabolic Rates from Multiple Closed Respirometry Loops
}
\description{
This function calculates the metabolic rates from multiple closed respirometry loops simultaneously. Requires lots of user input, but is easy to manipulate. Returns list of metabolic rates, as well as the average metabolic rate and the standard deviation of the sample of metabolic rates, as well as biglm objects for each section of data used to calculate MRs.
}
\usage{
MR.loops(data, DO.var.name, time.var.name = "std.time",
             in.DO.meas = "mg/L", out.DO.meas = "mg/L",
             start.idx, stop.idx, system.vol = 1,
             background.consumption = 0,
             background.indices = NULL,
             temp.C, elevation.m = NULL,
             bar.press = NULL, bar.units = "atm",
             PP.units, time.units = "sec", ...)
}

\arguments{
  \item{data}{
Must include a time variable in standard \code{POSIXct} format. eg 
"2016-09-25 15:30:00 EST".
}
  \item{DO.var.name}{
Column name of DO variable, must be entered as character string.
}
  \item{time.var.name}{
Column name of time variable (which is in POSIXct format) as character. defaults to "std.time" as generated from get.witrox.data().
 }
  \item{in.DO.meas}{
Units of DO measurement entered in the DO variable column: must be one of "mg/L" for miligrams/liter, "PP" for partial pressure, "pct" for saturation percent.
}
  \item{out.DO.meas}{
Units of DO measurement returned for metabolic rate: must be one of "mg/L" for miligrams/liter, "PP" for partial pressure (units determined by bar.units), "pct" for saturation percent.
}
  \item{start.idx}{
Character class value or vector matching {POSIXct} object coding for date time. Each element of the vector represents the start time of a new 'loop' for calculation of metabolic rates.
}
  \item{stop.idx}{
Character class value or vector matching "POSIXct" object coding for date time. Each element of the vector represents the stop time of a new loop for calculation of metabolic rates.
}
  \item{system.vol}{
System volume in Liters (defaults to 1 L).
}
  \item{background.consumption}{
Default = 0. If using a one point calibration for background, simply set "background.consumption" equal to the value of the calculated respiration rate. If using a multi-point calibration, enter a vector of background respiration rates, and enter a corresponding vector for "background.indices". CAUTION: The slope must be entered in raw units  (i.e. those specified in the input data.frame in the "data" argument). For example, if the "DO.var.name" column is recorded in mg/L, and the as.POSIXct format goes to the resolution of seconds, background consumption units would need to be entered in mgO2 / sec. 
}
  \item{background.indices}{
If using a multi-point calibration to set the background respiration rate, enter a vector of times for when the respiration rates were calculated. There should be one time point per corresponding value in "background.consumption". The background respiration rate is continually factored into all calculations of metabolic rate. The elements of the vector must entered as character strings conforming to the POSIXct format specified in the "time.var.name" column.
}
  \item{temp.C}{
Water temperature in degrees C.
}
  \item{elevation.m}{
Elevation in m. Only required if "bar.press" = NULL.
}
  \item{bar.press}{
barometric pressure in units defined by bar.units argument. Only required if "elevation.m" = NULL.
}
  \item{bar.units}{
Units of barometric pressure used as input and in output if DO.meas.out = "PP". Acceptable arguments: "mmHg", "atm", "kpa".
}
  \item{PP.units}{
Units of barometric pressure used for "PP". 
}
  \item{time.units}{
Denominator for metabolic rate, also displayed as units on X-axis. Acceptable arguments: "hr", "min", "sec".
}
  \item{\dots}{
Arguments passed on to internal functions
}
}
\value{
Returns a list of 2. $MR.summary is of class = data.frame with 3 columns: $MR (metabolic rate in user specified units, this is the same as the slope in each linear model), $sd.slope (standard deviation of slopes calculation), $r.square (adjusted r square value from each model). This second object is a list of "biglm" objects, each one representing a metabolic loop (see McDonnell and Chapman 2016). 
}
\references{
McDonnell, Laura H., and Lauren J. Chapman. "Effects of thermal increase on aerobic capacity and swim performance in a tropical inland fish." Comparative Biochemistry and Physiology Part A: Molecular & Integrative Physiology 199 (2016): 62-70.

Roche, Dominique G., et al. "Finding the best estimates of metabolic rates in a coral reef fish." Journal of Experimental Biology 216.11 (2013): 2103-2110.

}
\author{
Tyler L. Moulton
}
\seealso{
 \code{\link{as.POSIXct}},
 \code{\link{strptime}},
 \code{\link{background.resp}},
 \code{\link{Barom.Press}},
 \code{\link{Eq.Ox.conc}},
  \code{\link{biglm}},
 }
\examples{
## load data ##
data(fishMR)


## calc background resp rate
bgd.resp <- 
    background.resp(fishMR, "DO.mgL", 
                    start.time = "2015-07-02 16:05:00",
                    end.time = "2015-07-02 16:35:00",
                    ylab = "DO (mg/L)", xlab = "time (min)")

bg.slope.a <- bgd.resp$mat[2]


starts <- c("2015-07-03 01:15:00", "2015-07-03 02:13:00",
            "2015-07-03 03:02:00", "2015-07-03 03:50:00",
            "2015-07-03 04:50:00")

stops <- c("2015-07-03 01:44:00", "2015-07-03 02:35:30",
           "2015-07-03 03:25:00", "2015-07-03 04:16:00",
           "2015-07-03 05:12:00")

metR <- MR.loops(data = fishMR, DO.var.name ="DO.mgL",
                 start.idx = starts, time.units = "hr",
                 stop.idx = stops, time.var.name = "std.time",
                 temp.C = "temp.C", elevation.m = 1180,
                 bar.press = NULL, in.DO.meas = "mg/L", 
                 background.consumption = bg.slope.a,
                 ylim=c(6, 8))

metR$MR.summary


## now lets assume we ran a control loop for background rate
## before and after we ran the MR loops
## let:

bg.slope.b <-bg.slope.a -0.0001 
metRa <- MR.loops(data = fishMR, DO.var.name ="DO.mgL",
                 start.idx = starts, time.units = "hr",
                 stop.idx = stops, time.var.name = "std.time",
                 temp.C = "temp.C", elevation.m = 1180,
                 bar.press = NULL, in.DO.meas = "mg/L", 
                 background.consumption = c(bg.slope.a, bg.slope.b),
                 background.indices = c("2015-07-02 16:20:00",
                                        "2015-07-03 06:00:00"),
                 ylim=c(6, 8))


metRa$MR.summary

# note that the calculated slopes
# diverge as time increases. This is
# because the background respiration
# rate is increasing.

metR$MR.summary-metRa$MR.summary

## This looks great, but you need to check your start and
## stop vectors, otherwise, you could end up with some
## atrocious loops, e.g.:

starts <- c("2015-07-03 01:15:00", "2015-07-03 02:13:00",
            "2015-07-03 03:02:00", "2015-07-03 03:50:00",
            "2015-07-03 04:50:00")

stops <- c("2015-07-03 01:50:00", "2015-07-03 02:35:30",
           "2015-07-03 03:25:00", "2015-07-03 04:16:00",
           "2015-07-03 05:12:00")

metRb <- MR.loops(data = fishMR, DO.var.name ="DO.mgL",
         start.idx = starts,
         stop.idx = stops, time.var.name = "std.time",
         temp.C = "temp.C", elevation.m = 1180,
         bar.press = NULL, in.DO.meas = "mg/L", 
         background.consumption = bg.slope.a,
         ylim=c(6,8))

}
