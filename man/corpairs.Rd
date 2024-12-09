% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corpairs.R
\name{corpairs}
\alias{corpairs}
\title{corpairs: Correlation Pairs Plotting Function}
\usage{
corpairs(
  dt = demo,
  id.col = "Index",
  cor.method = c("pearson", "kendall", "spearman")[3],
  savefile = "outputfile",
  corpairplotsize = c(3 * 3, 2.5 * 3),
  bin = 50,
  logL = FALSE,
  plottype = c("hex", "point")[2],
  pointcolor = c("A", "B", "C", "D", "E")[3],
  pointsize = 1,
  kde2d.n = 50
)
}
\arguments{
\item{dt}{A data frame to be used for plotting. The first column should be an identifier column.}

\item{id.col}{The column name of the identifier column in \code{dt}. Default is "Index".}

\item{cor.method}{The correlation method to use: "pearson", "kendall", or "spearman". Default is "spearman".}

\item{savefile}{The directory where output plots will be saved. Default is "outputfile".}

\item{corpairplotsize}{A numeric vector of length 2 specifying the width and height of the plot in inches. Default is \code{c(9, 7.5)}.}

\item{bin}{Number of bins for hexbin plots. Default is \code{50}.}

\item{logL}{Logical, whether to log-transform the data. Default is \code{FALSE}.}

\item{plottype}{The type of plot for lower panels: "hex" or "point". Default is "point".}

\item{pointcolor}{A vector of color palette names for scatter plots, with a default of "C".}

\item{pointsize}{Size of points in scatter plots. Default is \code{1}.}

\item{kde2d.n}{Grid size for \code{MASS::kde2d} density calculation. Default is \code{50}.}
}
\value{
A \code{GGally::ggpairs} plot object.
}
\description{
This function provides correlation matrix visualizations using \code{GGally::ggpairs}.
It supports plotting scatter plots with density-based color mappings or hexbin plots.
}
\examples{
# Load the demo data included in the package
data("demo", package = "corhex")
corpairs(dt=demo,
id.col="Index",
cor.method=c("pearson", "kendall", "spearman")[3],
savefile="outputfile",
corpairplotsize=c(3*3,2.5*3),#width height
bin=50,
logL=FALSE,
plottype=c("hex","point")[2],
pointcolor=c("A","B","C","D","E")[3],
pointsize=1,
kde2d.n=50)
}
