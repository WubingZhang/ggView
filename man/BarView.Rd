% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BarView.R
\docType{methods}
\name{BarView}
\alias{BarView}
\title{Bar plot}
\usage{
BarView(
  df,
  x = "x",
  y = "y",
  fill = "#FC6665",
  color = fill,
  bar.width = 0.8,
  position = "dodge",
  dodge.width = 0.8,
  main = NA,
  xlab = NULL,
  ylab = NA,
  ...
)
}
\arguments{
\item{df}{A data frame.}

\item{x}{A character, specifying the x-axis.}

\item{y}{A character, specifying the x-axis.}

\item{fill}{A character, specifying the fill color.}

\item{color}{A character, specifying the outline color.}

\item{bar.width}{A numeric, specifying the width of bar.}

\item{position}{"dodge" (default), "stack", "fill".}

\item{dodge.width}{A numeric, set the width in position_dodge.}

\item{main}{A charater, specifying the figure title.}

\item{xlab}{A character, specifying the title of x-axis.}

\item{ylab}{A character, specifying the title of y-axis.}

\item{...}{Other parameters in geom_bar}
}
\value{
An object created by \code{ggplot}, which can be assigned and further customized.
}
\description{
Bar plot
}
\examples{
mdata = data.frame(group=letters[1:5], count=sample(1:100,5))
BarView(mdata, x = "group", y = "count")
}
\author{
Wubing Zhang
}
