\name{create.gene.meta.prognosis.plots}
\alias{create.gene.meta.prognosis.plots}
\title{Make Kaplan-Meier plots of many features by meta-dichotomizing each separately in each dataset}
\description{A major plotting function that takes in a meta-analytic dataset and a feature-name and outputs a survival curve based on dataset-wise median dichotomization of the feature.}
\usage{
create.gene.meta.prognosis.plots(file.stem, feature.names, tumour.types = 'all',
    mains = NA, rounding = 2, enable.warnings = TRUE, description = NULL);
}
\arguments{
  \item{file.stem}{The stem of the filename to be used to create each tumour-specific plot}
  \item{feature.names}{A vector giving the name of the feature to be dichotomized}
  \item{tumour.types}{Which tumour types should be run (defaults to all). NOTE: 'all' does not include prostate cancer methylation dataset.}
  \item{mains}{A vector of the main titles for each plot (must be the same for each plot for a feature). Defaults to the empty string}
  \item{rounding}{How many decimal places of precision should be shown in printed survival statistics?}
  \item{enable.warnings}{Print warnings if set to TRUE, defaults to FALSE}
  \item{description}{Short description of image/plot; defaults to NULL}
}
\value{Creates a survival curve for the meta-dichotomization of many feature in all site-specific cancer datasets.}
\author{Paul C. Boutros}
\seealso{
  \code{\link{make.survfit.plot}}, \code{\link{make.survival.plot}}, \code{\link{create.meta.prognosis.plot}}
  }
\examples{
#create.gene.meta.prognosis.plots(
#	file.stem = 'create_gene_meta_prognosis_plot',
#	feature.names = c('100_at', '1000_at'),
#	tumour.types = c('nsclc')
#	);
}
\keyword{hplot}
