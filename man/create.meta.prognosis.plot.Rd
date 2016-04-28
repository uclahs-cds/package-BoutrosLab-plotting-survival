\name{create.meta.prognosis.plot}
\alias{create.meta.prognosis.plot}
\title{Make a Kaplan-Meier plot by meta-dichotomizing a feature}
\description{A major plotting function that takes in a meta-analytic dataset and a feature-name and outputs a survival curve based on dataset-wise median dichotomization of the feature.}
\usage{
create.meta.prognosis.plot(expression.data, survival.data, feature.name, filename,
    xlab.label = "Time (Years)", ylab.label = "Fraction of Cohort", xlimits, xat,
    label.times, label.positions, label.cex = 1.2, main = NULL, main.cex = 3.0,
    rounding = 2, covariates = NULL, line.colours = c("blue", "red"), resolution = 2000,
    enable.warnings = FALSE, description = NULL);
}
\arguments{
  \item{expression.data}{A list of matrices, with each list element representing an expression matrix for a separate dataset}
  \item{survival.data}{A list of survival objects, with each list element representing a survobj for a separate dataset}
  \item{feature.name}{The name of the feature to be dichotomized}
  \item{filename}{The filename to be written}
  \item{xlab.label}{X Axis Title.  Defaults to 'Time (Years)'}
  \item{ylab.label}{Y Axis Title.  Defaults to 'Fraction of Cohort'}
  \item{xlimits}{The limits of the x (time) axis}
  \item{xat}{The locations of the x-axis labels}
  \item{label.times}{A vector of times at which to calculate the fraction of patients still at risk}
  \item{label.positions}{A vector of positions at which to print the fraction of patients still at risk}
  \item{label.cex}{Size of patients at risk text}
  \item{main}{The main title for the plot (space is reclaimed if NULL)}
  \item{main.cex}{Size of text for main plot title}
  \item{rounding}{How many decimal places of precision should be shown in printed survival statistics?}
  \item{covariates}{A data-frame (or matrix?) of variables to be controlled in the Cox model. If null, no adjustment is done.  No interactions are fit.}
  \item{line.colours}{A vector of colours to be used in plotting}
  \item{resolution}{Figure resolution, defaults to 2000}
  \item{enable.warnings}{Print warnings if set to TRUE, defaults to FALSE}
  \item{description}{Short description of image/plot, defaults to NULL}
}
\value{Returns a survival curve (written to file) and a four-element vector containing the survival results.}
\author{Paul C. Boutros}
\seealso{
  \code{\link{make.survfit.plot}}, \code{\link{make.survival.plot}}
  }
\examples{ # add examples }
\keyword{hplot}
