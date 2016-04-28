\name{create.gene.meta.prognosis.plot}
\alias{create.gene.meta.prognosis.plot}
\title{Make Kaplan-Meier plots by meta-dichotomizing a feature separately in each dataset}
\description{A major plotting function that takes in a meta-analytic dataset and a feature-name and outputs a survival curve based on dataset-wise median dichotomization of the feature.}
\usage{
create.gene.meta.prognosis.plot(file.stem, feature.name, tumour.types = 'all',
    main = '', rounding = 2, enable.warnings = TRUE, description = NULL);
}
\arguments{
  \item{file.stem}{The stem of the filename to be used to create each tumour-specific plot}
  \item{feature.name}{The name of the feature to be dichotomized}
  \item{tumour.types}{Which tumour types should be run (defaults to all). NOTE: 'all' does not include prostate cancer methylation dataset.}
  \item{main}{Add a main title to each plot (must be the same for each plot). Defaults to the empty string}
  \item{rounding}{How many decimal places of precision should be shown in printed survival statistics?}
  \item{enable.warnings}{Print warnings if set to TRUE, defaults to FALSE}
  \item{description}{Short description of image/plot; defaults to NULL}
}
\value{Creates a survival curve for the meta-dichotomization of a feature in all site-specific cancer datasets.}
\author{Paul C. Boutros}
\seealso{
  \code{\link{make.survfit.plot}}, \code{\link{make.survival.plot}}, \code{\link{create.meta.prognosis.plot}}
  }
\examples{
#create.gene.meta.prognosis.plot(
#	file.stem = 'create_gene_meta_prognosis_plot',
#	feature.name = '100_at',
#	tumour.types = 'ovarian'
#	);

#create.gene.meta.prognosis.plot(
#	file.stem = 'create_gene_meta_prognosis_plot',
#	feature.name = '100_at'
#	);
}
\keyword{hplot}
