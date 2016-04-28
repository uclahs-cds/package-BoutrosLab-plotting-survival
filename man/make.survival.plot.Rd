\name{make.survival.plot}
\alias{make.survival.plot}
\title{Utility function to make complete Kaplan-Meier plots}
\description{A utility function that makes a pretty Kaplan-Meier survival curve given survival data and patient grouping. This is a highly general function that does very little with the survival data itself.}
\usage{
make.survival.plot(survival.object, patient.groups, filename, xlimits, xat,
    label.positions, label.cex = 1.2, at.risk.low, at.risk.high,
    risk.labels = c("Low", "High"), main = NULL, main.cex = 3.0, rounding = 2,
    covariates = NULL, line.colours = c("blue", "red"), key = NULL,
    statistical.method = NA, x.statistic = 0.200, y.statistic = 0.610,
    resolution = 2000, height = 7.007874, width = 7.007874, size.units = 'in',
    enable.warnings = FALSE);
}
\arguments{
  \item{survival.object}{An object of class Surv -- patient ordering needs to be identical as for patient.groups}
  \item{patient.groups}{Grouping of patients (passed directly to coxph, so factors & continuous variables are okay)}
  \item{filename}{The filename to be written}
  \item{xlimits}{A two component vector listing the upper and lower bounds of the x-axis (time dimension)}
  \item{xat}{A vector listing the points at which x-axis labels should be drawn}
  \item{label.positions}{A vector of positions at which to print the fraction of patients still at risk}
  \item{label.cex}{Size of patients at risk text}
  \item{at.risk.low}{A vector containing the number of patients in the low-risk group at each label.position}
  \item{at.risk.high}{A vector containing the number of patients in the high-risk group at each label.position}
  \item{risk.labels}{A character(2) vector that gives the labels for the at-risk rows}
  \item{main}{The main title for the plot (space is reclaimed if NULL)}
  \item{main.cex}{Size of text for main plot title}
  \item{rounding}{How many decimal places of precision should be shown in printed survival statistics?}
  \item{covariates}{A data-frame (or matrix?) of variables to be controlled in the Cox model. If null, no adjustment is done.  No interactions are fit.}
  \item{line.colours}{A vector of colours to be used in plotting}
  \item{key}{A list describing the legend one wants to add to the plot, defaults to no legend}
  \item{statistical.method}{A character indicating what statistical test should be done ('cox', 'logrank', 'ttest'). If NA, tries to choose intelligently}
  \item{x.statistic}{x position for statistical data, defaults to left}
  \item{y.statistic}{y position for statistical data, defaults to low}
  \item{resolution}{Figure resolution in dpi}
  \item{height}{Figure height in size.units}
  \item{width}{Figure width in size.units}
  \item{size.units}{Units of size for the figure}
  \item{enable.warnings}{Print warnings if set to TRUE, defaults to FALSE}
}
\value{Creates a tiff file containing a full Kaplan-Meier curve, including number of patients at risk at multiple time-points.}
\author{Paul C. Boutros}
\seealso{\code{\link{make.survfit.plot}}, \code{\link[lattice]{xyplot}}, \code{\link[lattice]{lattice}} for an overview of the package.}
\examples{

# create some artificial data
survtime <- c(
	sample(1:75, 100, replace = TRUE),
	sample(1:99, 100, replace = TRUE)
	);

survstat <- c(
	sample( c(0,1,1), 100, replace = TRUE),
	sample( c(0,0,1), 100, replace = TRUE)
	)

survobj <- Surv(survtime, survstat);
groups <- as.factor(c( rep('A', 100), rep('B', 100) ));

# set some plotting parameters
xlimits <- c(0,100);
xat <- seq(0,100,25);
label.times <- seq(0,100,20);
label.positions <- 0.167 + 0.1605 * (1:length(label.times) - 1);

make.survival.plot(
	survival.object = survobj,
	patient.groups = groups,
	filename = 'survival_plot_simple.tiff',
	xlimits = xlimits,
	xat = xat,
	label.positions = label.positions,
	at.risk.low = calculate.number.at.risk(
		survobj = survobj[groups == 'B',],
		cut.points = label.times
		),
	at.risk.high = calculate.number.at.risk(
		survobj = survobj[groups == 'A',],
		cut.points = label.times
		),
	risk.labels = c('A = High', 'B = Low'),
	resolution = 100
	);
}
\keyword{hplot}
