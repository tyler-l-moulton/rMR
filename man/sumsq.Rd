\name{sumsq}
\alias{sumsq}

\title{
Sum of squares
}
\description{
Internal function for use in package for calculating sum of squares of a vector.
}
\usage{
sumsq(x)
}

\arguments{
  \item{x}{
Numeric vector to be evaluated
}
}
\details{
Internal function for package
}
\value{
The sum of squares of the vector
}
\author{
Tyler L. Moulton
}
\examples{
vec <- sample(c(100:120), 50, replace = TRUE)
sumsq(vec)
}
