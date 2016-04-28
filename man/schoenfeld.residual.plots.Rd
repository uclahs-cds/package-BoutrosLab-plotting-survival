\name{schoenfeld.residual.plots}
\alias{schoenfeld.residual.plots}
\title{Utility function to produce Schoenfeld residual plots}
\description{A utility function that produces Schoenfeld residual plots (in .tiff files), for each variable in a coxph model.  Relevant pvalues are printed on each plot.}
\usage{
schoenfeld.residual.plots(cox.model, filename);
}

\arguments{
  \item{cox.model}{The coxph object for which the Schoenfeld residual plots will be produced.}
  \item{filename}{String specifiying the root of the filename for the Schoenfeld residual plots.}
}

\value{For each variable in the coxph model, creates a tiff file containing a plot of the Schoenfeld residuals.}

\author{Nathalie Moon}

%\seealso{\code{\link{BoutrosLab.plotting.survival::create.survival.plot}},\code{\link{fit.coxmodel}}}
\seealso{\code{\link{fit.coxmodel}}}


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
three.groups <- factor(c( rep('A', 100), rep('B', 50), rep('C', 50) ), levels = c('A','B','C'));

cox.model <- fit.coxmodel(
	groups = three.groups,
	survobj = survobj,
	return.cox.model = TRUE
	)

schoenfeld.residual.plots(cox.model = cox.model, filename = "testing");

}

\keyword{hplot}
