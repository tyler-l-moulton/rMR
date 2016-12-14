\name{get.pcrit}
\alias{get.pcrit}

\title{
Calculate critical tension for rate processes
}
\description{
Determines the critical point of a rate process based on the broken stick model featured in Yeager and Ultsch (1989). The two regressions are selected based on the break point which minimizes the total residual sum of squares.
}
\usage{
get.pcrit(data, DO.var.name, MR.var.name = NULL, Pcrit.below,
             time.interval, time.var = NULL,
             start.time, stop.time, time.units = "sec",
             Pcrit.type = "both",...)
}
\arguments{
  \item{data}{
Data to be used.
}
  \item{DO.var.name}{
Variable name of oxygen concentration variable, formatted as character string. To be used for determination of critical point. Requires 'time.interval' to be specified.
}
  \item{MR.var.name}{
Metabolic rate variable name, formatted as character. Default = 'NULL'. Cannot have value unless DO = 'NULL'.
}
  \item{Pcrit.below}{
    DO concentration below which you are confident that Pcrit occurs. Accelerates process by reducing the number of iterations required to find Pcrit. Data points featuring DO conc > 'Pcrit.below' are still used to calculate regressions for model.
}
  \item{time.interval}{
    If MR = NULL, specify interval in seconds over which to calculate instantaneous MR.
}
  \item{time.var}{
    Column name for indexing(time) variable used to calculate instantaneous Metabolic Rate (MR) from oxygen concentration data (specified by DO.var.name). Must be a character string.
}
  \item{start.time}{
    Beginning of time interval over which to evaluate data for Pcrit. Required if MR.var.name = NULL.
}
  \item{stop.time}{
    End of time interval over which to evaluate data for Pcrit. Required if MR.var.name = NULL.
}
  \item{time.units}{
Units of time in MR calculation. Defaults to "sec", must be "sec", "min", or "hr". Required if MR.var.name = NULL.
}
  \item{Pcrit.type}{
Either "lm" to draw a vertical line at the Pcrit as determined by the intersection point of the best fit lines(Yeager and Ultsch 1989) or "midpoint" as determined by the midpoint between the two points on either side of the Pcrit (Yeager and Ultsch 1989). "both" will plot both as vertical lines on the plot. NULL will plot neither. Both values are returned in the output.
}
  \item{\dots}{
Arguments passed on to internal functions.
}
}
\details{
This calculates the critical oxygen tension for a change in metabolic rate. It is a simple broken stick model which evaluates the data at dissolved oxygen values < 'Pcrit.below'. The data of MR and DO are ordered by decreasing DO value. Then, the function iteratively calculates the total residual sum of squares (using tot.RSS) of two linear models, one spanning from 'Pcrit.below' to 'Pcrit.below' - i, the other with a range from the minimum DO value to 'Pcrit.below' - (i + 1). The broken stick model resulting in the lowest total residual sum of squares is selected. The Pcrit is the DO value at the intersection point of the broken stick model. This is indicated by a blue circle on the plot.
}
\value{
Returns a list of 6. $Pcrit.lm is the Pcrit given by the intersection of the two best fit lines. $P$Adj.r2.above gives the adjusted R2 value of the relationship between MR~DO above the critical point, and likewise, $Adj.r2.below gives the R2 below the critical point. The other 24 list elements are from the two linear models (denoted as above and below the pcrit) in the broken stick model.
}
\references{
Yeager, D. P. and Ultsch, G. R. (1989). Physiological regulation and conformation: a BASIC program for the determination of critical points. Physiological Zoology, 888-9
}
\author{
Tyler L. Moulton
}
\seealso{
\code{\link{tot.rss}},
\code{\link{strptime}}, 
\code{\link{as.POSIXct}},
}
\examples{
## set data ##

data(fishMR)

Pcrit1 <-get.pcrit(data = fishMR, DO = "DO.mgL",
                    Pcrit.below = 2,
                    time.var = "std.time",
                    time.interval = 120,
                    start.time = "2015-07-03 06:15:00",
                    stop.time = "2015-07-03 08:05:00")
## MR units in mgO2 / sec
            
## Change time interval ##
Pcrit2 <-get.pcrit(data = fishMR, DO = "DO.mgL",
                    Pcrit.below = 2,
                    time.var = "std.time",
                    time.interval = 60,
                    start.time = "2015-07-03 06:15:00",
                    stop.time = "2015-07-03 08:05:00",
                    time.units = "min")
## MR units in mgO2 / min

Pcrit3 <-get.pcrit(data = fishMR, DO = "DO.mgL",
                    Pcrit.below = 2,
                    time.var = "std.time",
                    time.interval = 60,
                    start.time = "2015-07-03 06:15:00",
                    stop.time = "2015-07-03 08:05:00",
                    time.units = "hr",
                    ylab = "Met Rate (mg O2 / hour)")

## No vertical lines on plot

Pcrit4 <-get.pcrit(data = fishMR, DO = "DO.mgL",
                    Pcrit.below = 2,
                    time.var = "std.time",
                    time.interval = 60,
                    start.time = "2015-07-03 06:15:00",
                    stop.time = "2015-07-03 08:05:00",
                    time.units = "hr",
                    ylab = "Met Rate (mg O2 / hour)",
                    Pcrit.type = "")
                    

}
