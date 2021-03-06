% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gseaView.R
\name{gseaView}
\alias{gseaView}
\title{gseaView}
\usage{
gseaView(x, geneSetID, title = "", seg_height = 0.08, base_size = 11)
}
\arguments{
\item{x}{gseaResult object}

\item{geneSetID}{gene set ID}

\item{title}{plot title}

\item{seg_height}{height of segments}

\item{base_size}{base font size}
}
\value{
ggplot object.
}
\description{
GSEA plot that mimic the plot generated by broad institute's GSEA software
}
\author{
Wubing Zhang
}
