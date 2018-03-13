\name{create.km.plot}
\alias{create.km.plot}
\title{Utility function to make complete Kaplan-Meier plots}
\description{A utility function that makes a pretty Kaplan-Meier survival curve given survival data and patient grouping. This is a highly general function that does very little with the survival data itself.}
\usage{
create.km.plot(
	survival.object, 
	patient.groups = NA, 
	filename = NULL, 
	xat = NA, 
	yat = seq(0,1,0.2),
	xlimits = NA, 
	ylimits = c(0,1.03),
	xaxis.cex = 2,
	yaxis.cex = 2,
	xlab.label = 'Time (Months)',
	xlab.cex = 2.75,
	ylab.label = 'Estimated Proportion',
	ylab.cex = 2.75,
        xaxis.fontface = 'bold',
        yaxis.fontface = 'bold',
	risk.labels = NA,
	key.groups.labels = levels(as.factor(patient.groups)),
	key.groups.cex = 2,
	risk.label.pos = NA,
	risk.label.fontface = 'bold',
	key.groups.title = NULL, 
	key.groups.title.cex = 1.4,
	key.stats.cex = 1.5,
	explicit.HR.label = TRUE,
	main = NULL, 
	main.cex = 3.0,
	covariates = NULL,
	stratification.factor = NULL,
	stratification.value = NULL,
	lwd = 2,
	lty = 1,
	censoring.pch.cex = 1.1,
	digits = 2, 
	line.colours = NA, 
	statistical.method = NA,
	predefined.p = NULL,
	predefined.hr = NULL,
	predefined.hr.ci = NULL,
	predefined.p.statistic.type = 'P',
	ph.assumption.check = "warning",
	cox.zph.threshold = 0.1, 
	cox.zph.truncation.thresholds = c(5, 10),
	show.key.groups = NA, 
	show.risktable = TRUE,
	risktable.fontsize = NULL, 
	key.groups.corner = c(0,0), 
	key.groups.x.pos = 0, 
	key.groups.y.pos = 0.01, 
	key.stats.corner = c(1,0), 
	key.stats.x.pos = 1, 
	key.stats.y.pos = 0.01, 
	ylab.axis.padding = 1,
	bottom.padding = 2,  
	top.padding = 0.1,
	right.padding = 0.1,
	left.padding = 0.5,
	return.statistics = FALSE,
	height = 7, 
	width = 7, 
	style = 'BoutrosLab',
	resolution = 1000, 
	size.units = 'in', 
	enable.warnings = TRUE,
	description = NULL,
	use.legacy.settings = FALSE 
	);
}

\arguments{
  \item{survival.object}{An object of class Surv -- patient ordering needs to be identical as for patient.groups}
  \item{patient.groups}{Factor vector representing the grouping of patients.  The first level of the factor is the baseline risk group, used if 'cox' is selected for statistical method.  Note that even if statistical.method = NA, cox method may be selected by the function, so a judicious choice should be made.  If patient.groups is not a factor, it is coerced into a factor variable, and the baseline risk group is taken to be the risk group that comes first alpha-numerically.  In the plot, the baseline risk group appears first in the risk table and in the legend, and its Kaplan-Meier estimate is plotted in blue.  To plot survival curves for a subset of the groups in a factor, only specifies the desired groups in the "levels" argument of factor.}
  \item{filename}{The filename to be written}
  \item{xat}{A vector listing the points at which x-axis labels should be drawn}
  \item{yat}{A vector listing the points at which y-axis labels should be drawn.  Defaults to seq(0,1,0.2)}
  \item{xlimits}{A two component vector listing the upper and lower bounds of the x-axis (time dimension)}
  \item{ylimits}{A two component vector listing the upper and lower bounds of the y-axis.  Defaults to c(0,1.03)}
  \item{xaxis.cex}{Size of x-axis scales, defaults to 2}
  \item{yaxis.cex}{Size of y-axis scales, defaults to 2}
  \item{xlab.label}{X Axis Title}
  \item{xlab.cex}{Size of x-axis labels, defaults to 2.75}
  \item{ylab.label}{Y Axis Title.  Defaults to "Proportion", but should be customized to better describe survival outcome of interest.}
  \item{ylab.cex}{Size of y-axis labels, defaults to 2.75}
  \item{xaxis.fontface}{font face of the xaxis, defaults to bold}
  \item{yaxis.fontface}{font face of the yaxis, defaults to bold}
  \item{risk.labels}{A character vector that gives the labels for risk groups in risk table.  Defaults to full name of each risk group if there is more than one risk group, and to a blank string otherwise.  Useful to replace long names by shorter ones.)}
  \item{key.groups.labels}{A character vector that gives the labels in the legend (listing name and line colour for each group).  Defaults to full name of risk groups.  Useful to replace long names by shorter ones.}
  \item{key.groups.cex}{Size of key groups text}
  \item{risk.label.pos}{A value for the x-axis position of the risk labels (in the risk table below the KM plot).  This value will typically be negative, in order for it to appear to the left of the value '0' on the x axis.}
  \item{risk.label.fontface}{font face of the risk group label.  Defaults to 'bold'.  Set it to NULL if plain font face is wanted for risk group label.}
  \item{key.groups.title}{Title for legend listing group names and line colours.  Defaults to NULL.}
  \item{key.groups.title.cex}{Font size of title for legend listing group names and line colours.  Defaults to 2.2 (default cex for group names is 2).}
  \item{key.stats.cex}{Font size of the statisitcal results key. Defaults to 1.5}
  \item{explicit.HR.label}{If TRUE, prints the label for the comparison risk group in subscript to HR.  Only has an effect if the cox statistical method is used.  Defaults to TRUE.}
  \item{main}{The main title for the plot (space is reclaimed if NULL)}
  \item{main.cex}{Size of text for main plot title}
  \item{covariates}{A data-frame (or matrix?) of variables to be controlled in the Cox model. If null, no adjustment is done.  No interactions are fit.}  
  \item{stratification.factor}{If building a cox model for statistical analysis purposes, the resulting cox model will be stratified on stratification.factor (if it is non-NULL).  User needs to provide a factor vector to specify stratification.factor Defaults to NULL.}
  \item{stratification.value}{A pre-determined value to split the stratification factor (i.e. strata(stratification.factor > stratification.value))}
  \item{lwd}{Specifies line width, defaults to 2.5}
  \item{lty}{Specifies line type, defaults to 1}
  \item{censoring.pch.cex}{Specifies length of tick marks indicating censoring events, defaults to 1.1}
  \item{digits}{How many decimal places of precision should be shown in printed survival statistics}
  \item{line.colours}{A vector of colours to be used in plotting.  Defaults to the survival palette in BoutrosLab.plotting.general::default.colours}
  \item{statistical.method}{A string indicating what statistical test should be done ('cox', 'logrank', 'ttest', 'none'). If NA, tries to choose intelligently}
  \item{predefined.p}{A predefined p value (e.g. from a log likelihood test) to display on the plot IN PLACE OF any other statistical test (above)}
  \item{predefined.hr}{A predefined hazard ratio value to display on the plot IN PLACE OF the hazard ratio calculated from a cox test.This value will only be displayed if the cox test is used. If the user provides a predefined hazard ratio, then they MUST specify a confidence interval (see predefined.hr.ci)}
  \item{predefined.hr.ci}{A vector of two numeric values representing the boundaries of the confidence interval associated with the predefined hazard ratio. This vector should only be provided if a predefined hazard ratio is specified.}
  \item{predefined.p.statistic.type}{String giving the type of predefined p-value statistic to be displayed. Passed directly to display.statistical.result function in BPG. Defaults to 'P'.}
  \item{ph.assumption.check}{A string indicating how the result of the Proportional Hazards assumption check (using the cox.zph function) should be reported.  If "warning", a warning is printed on the plot if the pvalue of cox.zph is smaller than cox.zph.threshold.  If "warning.and.plot", a warning is printed on the plot and a Schoenfeld residual plot is produced and saved if the pvalue is smaller than cox.zph.threshold.  If "logrank", the statistical method is changed to logrank if the pvalue is smaller than cox.zph.threshold.  If "residual.plot", the Schoenfled residual plot is produced and saved as a tiff, no matter what the results of the cox.zph test are.  If "ignore", the cox.zph test is not performed at all.  Defaults to "warning".}
  \item{cox.zph.threshold}{Numeric value, between 0 and 1).  If the pvalue of the cox.zph test is smaller than cox.zph.threshold, warnings will be provided to the user, indicating that the PH assumption may not be valid.  Defaults to 0.1}
  \item{cox.zph.truncation.thresholds}{If ph.assumption.check is not 'ignore', and a variable fails the CPH assumption, a multi-point table will be generated with CPH results for all variables at multiple time points. The truncation points specified here, indicate the time points at which to split the survival data. See BoutrosLab.statistics.general/R/multi.point.HR.table.R}
  \item{show.key.groups}{Defaults to FALSE if there is only one risk group, to TRUE otherwise}
  \item{show.risktable}{Defaults to TRUE if there are optimal display values listed for the corresponding number of risk groups in the optimal.risktable.parameters.txt file in the inst directory (currently for 1 to 10 risk groups).  Otherwise, defaults to FALSE.}
  \item{risktable.fontsize}{font size of the risk group table. Defaults to optimal.risktable.parameters.}
  \item{key.groups.corner}{A two-component vector specifiying which corner of the groups key will be anchored by key.groups.x.pos and key.groups.y.pos.  Defaults to c(0,0), i.e. bottom left corner.}
  \item{key.groups.x.pos}{Specifies the x position of the groups key (on scale of 0-1).  Defaults to 0}
  \item{key.groups.y.pos}{Specifies the y position of the groups key (on scale of 0-1).  Defaults to 0.01}
  \item{key.stats.corner}{A two-component vector specifying which corner of the statistical results key will be anchored by key.stats.x.pos and key.stats.y.pos.  Defaults to c(1,0), i.e. bottom right corner.}
  \item{key.stats.x.pos}{Specifies the x position of the statistical results key (on scale of 0-1).  Defaults to 1}
  \item{key.stats.y.pos}{Specifies the y position of the statistical results key (on scale of 0-1).  Defaults to 0.01}
  \item{ylab.axis.padding}{A number specifying the distance of ylabel to the y-axis, defaults to 1.  This parameter (along with left.padding) can be changed to make more space to display long group names in the risk table.}
  \item{bottom.padding}{A number specifying the distance to the bottom margin, defaults to 0.7}
  \item{top.padding}{A number specifying the distance to the top margin, defaults to 0.1}
  \item{right.padding}{A number specifying the distance to the right margin, defaults to 0.1}
  \item{left.padding}{A number specifying the distance to the left margin, defaults to 0.5.  This parameter (along with ylab.axis.padding) can be changed to make more space to display long group names in the risk table.}
  \item{return.statistics}{If TRUE, function will return a list containing trellis object (or 0/1 indicator) and relevant statistics.  If FALSE, no statistics are returned.  Defaults to FALSE.}
  \item{height}{Figure height in size.units}
  \item{width}{Figure width in size.units}
  \item{style}{defaults to \dQuote{BoutrosLab}, also accepts \dQuote{Nature}, which changes parameters according to Nature formatting requirements}
  \item{resolution}{Figure resolution in dpi}
  \item{size.units}{Units of size for the figure}
  \item{enable.warnings}{Print warnings if set to TRUE, defaults to TRUE}
  \item{description}{short description of image/plot; defaults to null}
  \item{use.legacy.settings}{boolean to set wheter or not to use legacy mode settings (font)}
}

\value{If return.stats is TRUE: if \code{filename} is \code{NULL} then returns the trellis object, otherwise creates a plot and returns a 0/1 success code.  If return.stats is FALSE: returns a list whose first element is the trellis object if \code{filename}) is \code{NULL}, 0/1 otherwise, and whose other elements are relevent statistics.}
\author{Nathalie Moon}
\seealso{\code{\link{make.survival.plot}},\code{\link{make.survfit.plot}},\code{\link{create.scatterplot}}, \code{\link[lattice]{xyplot}}, \code{\link[lattice]{lattice}} for an overview of the package.}


\examples{
set.seed(100);

survtime <- c(
	sample(1:75, 100, replace = TRUE),
	sample(1:99, 100, replace = TRUE)
	);

survstat <- c(
	sample( c(0,0,0,1,1), 100, replace = TRUE),
	sample( c(0,0,0,1), 100, replace = TRUE)
	);

survobj <- Surv(survtime, survstat);
two.groups <- factor(c( rep('A', 100), rep('B', 100)), levels = c('A','B'));
three.groups <- factor(c( rep('A', 100), rep('B', 50), rep('C', 50) ), levels = c('A','B','C'));

# set some plotting parameters
xlimits <- c(0,100);
xat <- seq(0,100,25);

# simple survival plot with 2 risk groups
BoutrosLab.plotting.survival::create.km.plot(
	survival.object = survobj,
	patient.groups = two.groups,
	filename = "simple.survival.plot.tiff",
	xlimits = xlimits,
	xat = xat,
	resolution = 200
	);

# survival plot with 2 risk groups, group B specified as reference group for cox analysis
BoutrosLab.plotting.survival::create.km.plot(
	survival.object = survobj,
	patient.groups = factor(two.groups, levels=c("B","A")),
	filename = "simple.survival.plot2.tiff",
	xlimits = xlimits,
	xat = xat,
	resolution = 200
	);

# survival plot for single group
BoutrosLab.plotting.survival::create.km.plot(
	survival.object = survobj,
	filename = "single.group.survival.plot.tiff",
	xlimits = xlimits,
	xat = xat,
	resolution = 200
	);

# survival plot with three risk groups, customized (longer) group names
BoutrosLab.plotting.survival::create.km.plot(
	survival.object = survobj,
	patient.groups = three.groups,
	filename = "three.group.survival.plot.tiff",
	xlimits = xlimits,
	xat = xat,
	key.groups.labels = c("Group A","Group B","Group C"),
	risk.labels = c("First Group","Second Group","Third Group"),
	ylab.axis.padding = 3.5,
	left.padding = 3.5,
	risk.label.pos = -5,
	resolution = 200
	);

# Creating data for six groups
survtime.six.groups <- c(
	sample(1:75, 105, replace = TRUE),
	sample(1:99, 105, replace = TRUE)
	);

survstat.six.groups <- c(
	sample(c(0,1), 35, replace = TRUE),
	sample(c(0,1,1), 35, replace = TRUE),
	sample(c(0,1,1,1), 35, replace = TRUE),
	sample(c(0,0,1,1,1), 35, replace = TRUE),
	sample(c(0,0,1,1,1,1), 35, replace = TRUE),
	sample(c(0,0,0,0,0,1,1), 35, replace = TRUE)
	);

survobj.six.groups <- Surv(survtime.six.groups, survstat.six.groups);

six.groups <- factor(
	c(
		rep('A',35),
		rep('B',35),
		rep('C',35),
		rep('D',35),
		rep('E',35),
		rep('F',35)
		), 
	levels = c('A','B','C','D','E','F')
	);

# Survival plot with six risk groups, legend and statistical results in 
# different locations
BoutrosLab.plotting.survival::create.km.plot(
	survival.object = survobj.six.groups,
	patient.groups = six.groups,
	filename = "simple.survival.plot.six.groups.tiff",
	xlimits = xlimits,
	xat = xat,
	key.stats.corner = c(0,0), # bottom left corner selected
	key.stats.x.pos = 0,
	key.stats.y.pos = 0.01,
	key.groups.x.pos = 0,
	key.groups.y.pos = 0.09,
	resolution = 200
	);

### Merge levels of a factor and then plot 
survtime.stage <- c(
	sample(1:99, 100, replace = TRUE),
	sample(1:85, 100, replace = TRUE),
	sample(1:60, 100, replace = TRUE),
	sample(1:55, 100, replace = TRUE)
	);

survstat.stage <- c(
	sample( c(0,0,0,0,0,1), 100, replace = TRUE),
	sample( c(0,0,0,1), 100, replace = TRUE),
	sample( c(0,0,1,1), 100, replace = TRUE),
	sample( c(0,1,1,1), 100, replace = TRUE)
	);

survobj.stage <- Surv(survtime.stage, survstat.stage);
groups.stage <- factor(c( rep('I', 100), rep('II', 100), rep('III', 100), rep('IV',100)),
    levels = c('I','II','III','IV'));

# Survival plot with 4 risk groups, one for each stage
BoutrosLab.plotting.survival::create.km.plot(
	survival.object = survobj.stage,
	patient.groups = groups.stage,
	filename = "stage.survival.plot.tiff",
	xlimits = xlimits,
	xat = xat,
	resolution = 200
	);

# merge stages I and II together, and stages II and III together by defining
# new levels for the factor variables
merged.stage.groups <- factor(x = groups.stage);
levels(merged.stage.groups) <- list(
	low.stage = c("I","II"),
	high.stage = c("III","IV")
	);

# Survival plot with merged stage groups
BoutrosLab.plotting.survival::create.km.plot(
	survival.object = survobj.stage,
	patient.groups = merged.stage.groups,
	filename = "merged.stage.survival.plot.tiff",
	xlimits = xlimits,
	xat = xat,
	left.padding = 1.5,
	resolution = 200
	);

# Survival plot with group 3 vs 2 only, and produce Schoenfeld residual plot to check PH assumption
BoutrosLab.plotting.survival::create.km.plot(
	survival.object = survobj.stage,
	patient.groups = factor(groups.stage, levels = c("II","III")),
	filename = "subset.of.stages.survival.plot.tiff",
	ph.assumption.check = "residual.plot",
	xlimits = xlimits,
	xat = xat,
	resolution = 200
	);

# Create data that will fail PH assumption
survtime.ph <- c(
	sample(1:100, 100, replace = TRUE),
	sample(1:60, 50, replace = TRUE),
	sample(40:100, 50, replace = TRUE)
	);

survstat.ph <- c(
	sample( c(0,1), 100, replace = TRUE),
	sample( c(0,0,0,0,1), 50, replace = TRUE),
	sample( c(0,1,1,1,1), 50, replace = TRUE)
	);

survobj.ph <- Surv(survtime.ph, survstat.ph);
groups.ph <- factor(c(rep('A', 100), rep('B', 100)), levels = c('A','B'));
gender.ph <- factor(sample(c("Male","Female"), 200, replace = TRUE), levels = c("Male","Female"));
score.ph <- rgamma(n=200, shape = 20, rate=1); # this is just an arbitrary "continuous" variable

# set some plotting parameters
xlimits <- c(0,100);
xat <- seq(0,100,25);

# Survival plot with PH assumption checking (warning on plot)
BoutrosLab.plotting.survival::create.km.plot(
	survival.object = survobj.ph,
	patient.groups = factor(groups.ph, levels = c("A","B")),
	covariates = data.frame(gender.ph, score.ph),
	filename = "ph.warning.survival.plot.tiff",
	xlimits = xlimits,
	xat = xat,
	ph.assumption.check = "warning.and.plot",
	resolution = 200
	);

# Survival plot with PH assumption checking (switch to logrank test if fails)
BoutrosLab.plotting.survival::create.km.plot(
	survival.object = survobj.ph,
	patient.groups = factor(groups.ph, levels = c("A","B")),
	filename = "ph.logrank.survival.plot.tiff",
	xlimits = xlimits,
	xat = xat,
	ph.assumption.check = "logrank",
	resolution = 200
	);

# Example with and without stratification
set.seed(151)
survtime <- sample(seq(0.1,10,0.1), 100, replace = TRUE);
survstat <- sample(c(0,1), 100, replace = TRUE);
survobj <- Surv(survtime, survstat);
groups <- sample(c('A','B'), 100, replace = TRUE);
gender <- factor(sample(c('Male','Female'), 100, replace = TRUE), levels=c("Male","Female"));
age <- factor(sample(c('Young','MiddleAged','Old'), 100, replace = TRUE),
    levels = c('Young', 'MiddleAged', 'Old'));

BoutrosLab.plotting.survival::create.km.plot(
	survival.object = survobj,
	patient.groups = as.factor(groups),
	filename = "surv.plot.no.stratification.tiff",
	covariates = data.frame(age, gender)
	);

BoutrosLab.plotting.survival::create.km.plot(
	survival.object = survobj,
	patient.groups = as.factor(groups),
	filename = "surv.plot.with.stratification.tiff",
	covariates = data.frame(age),
	stratification.factor = gender
	);

}
\keyword{hplot}
