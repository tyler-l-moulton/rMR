\name{plot.raw}
\alias{plot.raw}

\title{
Plotting data from witrox
}
\description{
A good way to visualize your respiro data to get an idea of where to set up the time intervals in functions like 'MR.loops' or 'get.pcrit'.
}
\usage{
plot.raw(data, DO.var.name, time.var.name = "std.time",
             start.time = data$x[1],
             end.time = data$x[length(data$x)])
}

\arguments{
 \item{data}{
data object for plotting
}
 \item{DO.var.name}{
A character string matching the column header for the DO variable column
}
  \item{time.var.name}{
Column name of time (or x) axis in character class. 
}
  \item{start.time}{
Character string specifying left bound x limit in 'strptime' compatable format
}
  \item{end.time}{
Character string specifying right bound x limit in 'strptime' compatable format
}
}
\details{
start.time and end.time arguments must match the time.var.name column's format for date time.
}
\value{
Plot showing the overall metabolic data
}
\references{

}
\author{
Tyler L. Moulton
}
\note{

}

\seealso{
code{\link{plot}},
code{\link{strptime}},
code{\link{get.pcrit}},
code{\link{MR.loops}},

}
\examples{
## load data ##

data(fishMR)

plot.raw(data = fishMR, DO.var.name = "DO.mgL",
         start.time = "2015-07-03 06:15:00",
         end.time = "2015-07-03 08:05:00")
         
plot.raw(fishMR, DO.var.name = "DO.mgL",
         start.time = "2015-07-03 01:00:00",
         end.time = "2015-07-03 05:12:00")

}
