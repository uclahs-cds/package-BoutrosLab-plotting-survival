# The BoutrosLab.plotting.survival package is copyright (c) 2012 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

create.gene.meta.prognosis.plots <- function(file.stem, feature.names, tumour.types = 'all', mains = NA, rounding = 2, enable.warnings = TRUE, description = NULL) {

	# handle requests for all tumour types
	if (tumour.types == 'all') {
		tumour.types <- c('breast', 'nsclc', 'ovarian', 'gbm', 'cervix', 'prostate', 'colon', 'breast.miRNA');
		}

	# handle case mismatches
	tumour.types <- tolower(tumour.types);

	# loop through requested tumour types
	for (tumour.type in tumour.types) {

		# load the dataset
		if (tumour.type == 'breast' && requireNamespace("BoutrosLab.datasets.breast.cancer")) {
			expression.data <- BoutrosLab.datasets.breast.cancer::load.breast.cancer.datasets();
			xlimits <- c(0,25);
			xat <- seq(0,25,5);
			label.times <- seq(0,25,5);
			label.positions <- 0.167 + 0.1605 * (1:length(label.times) - 1);
			}
		else if (tumour.type == 'nsclc' && requireNamespace("BoutrosLab.datasets.nsclc")) {
			expression.data <- BoutrosLab.datasets.nsclc::load.nsclc.datasets();
			xlimits <- c(0,20);
			xat <- seq(0,20,4);
			label.times <- seq(0,20,4);
			label.positions <- 0.167 + 0.1605 * (1:length(label.times) - 1);
			}
		else if (tumour.type == 'ovarian' && requireNamespace("BoutrosLab.datasets.ovarian.cancer")) {
			expression.data <- BoutrosLab.datasets.ovarian.cancer::load.ovarian.cancer.datasets();
			xlimits <- c(0,20);
			xat <- seq(0,20,4);
			label.times <- seq(0,20,4);
			label.positions <- 0.167 + 0.1605 * (1:length(label.times) - 1);
			}
		else if (tumour.type == 'gbm' | tumour.type == 'glioblastoma' && requireNamespace("BoutrosLab.datasets.gbm")) {
			expression.data <- BoutrosLab.datasets.gbm::load.gbm.datasets();
			xlimits <- c(0,20);
			xat <- seq(0,20,4);
			label.times <- seq(0,20,4);
			label.positions <- 0.167 + 0.1605 * (1:length(label.times) - 1);
			}
		else if (tumour.type == 'cervix' && requireNamespace("BoutrosLab.datasets.cervix.cancer")) {
			expression.data <- BoutrosLab.datasets.cervix.cancer::load.cervix.cancer.datasets();
			xlimits <- c(0,20);
			xat <- seq(0,20,4);
			label.times <- seq(0,20,4);
			label.positions <- 0.167 + 0.1605 * (1:length(label.times) - 1);
			}
		else if (tumour.type == 'prostate' && requireNamespace("BoutrosLab.datasets.prostate.cancer")) {
			expression.data <- BoutrosLab.datasets.prostate.cancer::load.prostate.cancer.datasets();
			xlimits <- c(0,20);
			xat <- seq(0,20,4);	
			label.times <- seq(0,20,4);
			label.positions <- 0.167 + 0.1605 * (1:length(label.times) - 1);
			}
		else if (tumour.type == 'colon' && requireNamespace("BoutrosLab.datasets.colon.cancer")) {
			expression.data <- BoutrosLab.datasets.colon.cancer::load.colon.cancer.datasets();
			xlimits <- c(0,20);
			xat <- seq(0,20,4);
			label.times <- seq(0,20,4);
			label.positions <- 0.167 + 0.1605 * (1:length(label.times) - 1);
			}
		else if (tumour.type == 'prostate.methylation' && requireNamespace("BoutrosLab.datasets.prostate.cancer.methylation")) {
			expression.data <- BoutrosLab.datasets.prostate.cancer.methylation::load.prostate.cancer.methylation.datasets(with.survival.only = TRUE);
			xlimits <- c(0,10);
			xat <- seq(0,10,2);
			label.times <- seq(0,10,2);
			label.positions <- 0.167 + 0.1605 * (1:length(label.times) - 1);
			}
		else if (tumour.type == 'prostate.CNV' && requireNamespace("BoutrosLab.datasets.prostate.cancer.CNV")) {
			expression.data <- BoutrosLab.datasets.prostate.cancer.CNV::load.prostate.cancer.datasets.CNV(with.survival.only = TRUE);
			xlimits <- c(0,10);
			xat <- seq(0,10,2);
			label.times <- seq(0,10,2);
			label.positions <- 0.167 + 0.1605 * (1:length(label.times) - 1);
			}
		else if (tumour.type == 'breast.miRNA' && requireNamespace("BoutrosLab.datasets.breast.cancer.miRNA")) {
			expression.data <- BoutrosLab.datasets.breast.cancer.miRNA::load.miRNA.breast.cancer.datasets();
			xlimits <- c(0,20);
			xat <- seq(0,20,4);
			label.times <- seq(0,20,4);
			label.positions <- 0.167 + 0.1605 * (1:length(label.times) - 1);
			}
			
		# verify correct loading or skip this tumour type
		if (class(expression.data) != 'list') {
			warning('Loading failed: ', tumour.type);
			next;
			}

		# Create the appropriate plot for each gene
		for (i in 1:length(feature.names)) {
			BoutrosLab.plotting.survival::create.meta.prognosis.plot(
				expression.data = expression.data$all.data,
				survival.data = expression.data$all.survobj,
				feature.name = feature.names[i],
				filename = paste(file.stem, feature.names[i], tumour.type, 'SurvivalCurve.tiff', sep = '_'),
				xlimits = xlimits,
				xat = xat,
				label.times = label.times,
				label.positions = label.positions,
				rounding = rounding,
				main = mains[i],
				enable.warnings = enable.warnings,
				resolution = 20, 
				description = description
				);
			}
		}
	}
