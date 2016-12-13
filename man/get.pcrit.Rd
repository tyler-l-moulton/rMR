\name{get.pcrit}
\alias{get.pcrit}

\title{
Calculate critical tension for rate processes
}
\description{
Determines the critical point of a rate process based on the broken stick model featured in Yeager and Ultsch (1989). The two regressions are selected based on the break point which minimizes the total residual sum of squares.
}
\usage{
get.pcrit(data, DO = NULL, MR = NULL,
Pcrit.below, idx.interval = NULL,
index.var = "std.time", start.idx = NULL,
stop.idx = NULL, time.units = "sec", Pcrit.type = "both",...)
}
\arguments{
  \item{data}{
Data to be used.
}
  \item{DO}{
Variable name of dissolved oxygen variable, formatted as character string. To be used for determination of critical point. Default = 'NULL'. Cannot have value unless MR = 'NULL'. Requires 'idx.interval' to be specified.
}
  \item{MR}{
Metabolic rate variable name, formatted as character. Default = 'NULL'. Cannot have value unless DO = 'NULL'.
}
  \item{Pcrit.below}{
    DO concentration below which you are confident that Pcrit occurs. Accelerates process by reducing the number of iterations required to find Pcrit. Data points featuring DO conc > 'Pcrit.below' are still used to calculate regressions for model.
}
  \item{idx.interval}{
    If MR = NULL, specify interval in seconds over which to calculate instantaneous MR.
}
  \item{index.var}{
    Column name for indexing(time) variable used to calculate instantaneous MR.
}
  \item{start.idx}{
    Beginning of time interval over which to evaluate data for Pcrit. Default = 'NULL'. If calculating Pcrit directly from MR (i.e. DO = 'NULL') start.idx will not be evaluated by the function.
}
  \item{stop.idx}{
    End of time interval over which to evaluate data for Pcrit. Default = 'NULL'. If calculating Pcrit directly from MR (i.e. DO = 'NULL') stop.idx will not be evaluated by the function.
}
  \item{time.units}{
Units of time in MR calculation. Defaults to "sec", must be "sec", "min", or "hr".
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
\note{

}

\seealso{
\code{\link{tot.RSS}},
\code{\link{tot.RSS}},
\code{\link{strptime}}, 
\code{\link{as.POSIXct}},
}
\examples{
## set data ##

data(fishMR)

Pcrit1 <-get.pcrit(data = fishMR, DO = "DO.mgL",
                    Pcrit.below = 2,
                    idx.interval = 120,
                    start.idx = "2015-07-03 06:15:00",
                    stop.idx = "2015-07-03 08:05:00")
## MR units in mgO2 / sec
            
## Change time interval ##
Pcrit2 <-get.pcrit(data = fishMR, DO = "DO.mgL",
                    Pcrit.below = 2,
                    idx.interval = 60,
                    start.idx = "2015-07-03 06:15:00",
                    stop.idx = "2015-07-03 08:05:00",
                    time.units = "min")
## MR units in mgO2 / min

Pcrit3 <-get.pcrit(data = fishMR, DO = "DO.mgL",
                    Pcrit.below = 2,
                    idx.interval = 60,
                    start.idx = "2015-07-03 06:15:00",
                    stop.idx = "2015-07-03 08:05:00",
                    time.units = "hr",
                    ylab = "Met Rate (mg O2 / hour)")

## No vertical lines on plot

Pcrit4 <-get.pcrit(data = fishMR, DO = "DO.mgL",
                    Pcrit.below = 2,
                    idx.interval = 60,
                    start.idx = "2015-07-03 06:15:00",
                    stop.idx = "2015-07-03 08:05:00",
                    time.units = "hr",
                    ylab = "Met Rate (mg O2 / hour)",
                    Pcrit.type = NULL)
                    

}
