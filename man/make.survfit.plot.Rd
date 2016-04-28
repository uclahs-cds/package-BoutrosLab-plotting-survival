\name{make.survfit.plot}
\alias{make.survfit.plot}
\title{Utility function to make Kaplan-Meier plots}
\description{A utility function that makes half of each Kaplan-Meier survival curve.  It needs to be called twice, once for patients who had an event and once for the overall curve.}
\usage{
make.survfit.plot(x, xlimits, xat, type = 1, main = NULL, main.cex = 3.0,
    line.colours = c("blue", "red"), enable.warnings = FALSE,...);
}
\arguments{
  \item{x}{An object of class survfit}
  \item{xlimits}{A two component vector listing the upper and lower bounds of the x-axis (time dimension)}
  \item{xat}{A vector listing the points at which x-axis labels should be drawn}
  \item{type}{If type = 1, draws the step curve.  If type = 2, draws the censored points}
  \item{main}{The main title for the plot (space is reclaimed if NULL)}
  \item{main.cex}{Size of text for main plot title}
  \item{line.colours}{A vector of colours to be used in plotting}
  \item{enable.warnings}{Print warnings if set to TRUE, defaults to FALSE}
  \item{...}{Any additional parameters to be passed to xyplot}
}
\details{This function is rarely called directly, but instead is used in larger plotting functions.}
\value{Takes a survfit object and returns a lattice object containing part of a Kaplan-Meier curve. For most practical cases this needs to be called twice, once with type = 1 and once with type = 2.}
\author{Paul C. Boutros}
\seealso{\code{\link[lattice]{xyplot}}, \code{\link[lattice]{lattice}} for an overview of the package.}
\examples{
# create some synthetic survival data
survstat <- sample(c(0,1), 100, TRUE);
survtime <- runif(100) * 100 / 12;
survobj  <- Surv(survtime, survstat);
groups   <- sample(c(0,1), 100, TRUE);

# set the graphics driver
current.type <- getOption("bitmapType");
options(bitmapType = "cairo");

# open output file
tiff(
	filename = "Test Plot.tiff",
	height = 7,
	width = 7,
	units = "in",
	res = 100,
	compression = "lzw"
	);

# plot the step curve
plot(
	x = make.survfit.plot(
		x = survfit(survobj ~ groups),
		xlimits = c(0,9),
		xat = seq(0,9,3),
		type = 1,
		main = 'Test Plot',
		enable.warnings = TRUE
		),
	position = c(0.03,0.2,1.0,1.0),
	newpage = TRUE
	);

# without creating a new plot, add the censored points
plot(
	x = make.survfit.plot(
		x = survfit(survobj ~ groups),
		xlimits = c(0,9),
		xat = seq(0,9,3),
		type = 2,
		main = 'Test Plot',
		enable.warnings = TRUE
		),
	position = c(0.03,0.2,1.0,1.0),
	newpage = FALSE
	);

# finish plot
dev.off();
options(bitmapType = current.type);
}
\keyword{hplot}
