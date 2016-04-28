# The BoutrosLab.plotting.survival package is copyright (c) 2013 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

create.meta.prognosis.plot <- function(
	expression.data,
	survival.data,
	feature.name,
	filename,
	xlab.label = 'Time (Years)',
	ylab.label = 'Fraction of Cohort',
	xlimits,
	xat,
	label.times,
	label.positions,
	label.cex = 1.2,
	main = NULL,
	main.cex = 3.0,
	rounding = 2,
	covariates = NULL,
	line.colours = c("blue", "red"),
	resolution = 2000,
	enable.warnings = FALSE,
	description = NULL
	) {

	# verify that we got appropriate input data
	if ("list" != class(expression.data)) { return( rep(NA,4) ); }
	if ("list" != class(survival.data)) { return( rep(NA,4) ); }
	if (length(expression.data) != length(survival.data)) { return( rep(NA,4) ); }

	# dichotomize the dataset
	dichotomized.data <- BoutrosLab.prognosticsignature.general::dichotomize.meta.dataset(
		feature.name    = feature.name,
		expression.data = expression.data,
		survival.data   = survival.data
		);

	# create locals to store the patient information
	overall.patient.groups   <- dichotomized.data[[1]];
	overall.patient.survival <- dichotomized.data[[2]];
	overall.patient.status   <- dichotomized.data[[3]];

	# check that not all the data is missing (i.e. feature found in at least one dataset)
	if (0 == length(overall.patient.groups)   ) { return( rep(NA,4) ); }
	if (0 == length(overall.patient.survival) ) { return( rep(NA,4) ); }
	if (0 == length(overall.patient.status)   ) { return( rep(NA,4) ); }

	# create an overall survival object
	overall.survival.object <- survival::Surv(overall.patient.survival, overall.patient.status);

	# calculate the number of patients at risk
	at.risk.low  <- BoutrosLab.statistics.survival::calculate.number.at.risk(overall.survival.object[overall.patient.groups == FALSE,], label.times);
	at.risk.high <- BoutrosLab.statistics.survival::calculate.number.at.risk(overall.survival.object[overall.patient.groups == TRUE, ], label.times);

	# create a survival plot
	survival.data <- BoutrosLab.plotting.survival::create.km.plot(
		survival.object = overall.survival.object,
		patient.groups  = factor(x = as.character(overall.patient.groups), levels = c('0','1'), labels = c("Low","High")),
		filename        = filename,
		xlab.label      = xlab.label,
		ylab.label      = ylab.label,
		xlimits         = xlimits,
		xat             = xat,
		main            = main,
		main.cex        = main.cex,
		digits          = rounding,
		covariates      = covariates,
		line.colours    = line.colours,
		resolution      = resolution,
		enable.warnings = enable.warnings,
		description = description
		);

	# try to prevent any memory leakage
	BoutrosLab.utilities::GarbageCollect();

	# return statistics
	return(survival.data);

	}
