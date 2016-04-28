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

create.meta.multifeature.prognosis.plot <- function(
	expression.data,
	survival.data,
	feature.names,
	feature.weights,
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
	rounding = 3,
	covariates = NULL,
	survival.file = NA,
	line.colours = c("blue", "red"),
	resolution = 2000,
	enable.warnings = FALSE,
	description = NULL
	) {

	# verify that we got appropriate input data
	if ('list' != class(expression.data)) { warning('Invalid expression.data'); return(NA); }
	if ('list' != class(survival.data)) { warning('Invalid survival.data'); return(NA); }
	if (length(expression.data) != length(survival.data)) { warning('Invalid survival.data length'); return(NA); }
	if (length(feature.names) != length(feature.weights)) { warning('Invalid feature.weights length'); return(NA); }
	if (length(names(expression.data)) == 1) { if (is.null(names(expression.data))) { warning("Ensure lists are named"); return(NA); } }

	# covariates will *not* currently work unless the meta-object is already unflattened -- this needs fixing!
	if ('list' == class(covariates)) { stop('covariates handling is borked: only works if the meta-analysis object is already flattened'); }
	if (!is.null(covariates)) { warning('covariates behaviour is borked: will only work if the meta-analysis object is already flattened'); }

	# calculate signature performance
	signature.results <- BoutrosLab.prognosticsignature.general::score.multifeature.signature(
		expression.data = expression.data,
		survival.data = survival.data,
		feature.names = feature.names,
		feature.weights = feature.weights
		);

	# verify we got back a proper analysis that we can interpret
	if (1 == length(signature.results) && is.na(signature.results)) { warning("Signature scoring failed"); return( rep(NA,4) ); }

	# if desired print the per-patient annotation to file
	if (!is.na(survival.file)) {
		write.table(
			x = signature.results,
			file = survival.file,
			sep = "\t",
			col.names = NA
			);
		}

	# calculate the number of patients at risk
	survobj <- survival::Surv(signature.results$survtime, signature.results$survstat);
	at.risk.low  <- BoutrosLab.statistics.survival::calculate.number.at.risk(survobj[signature.results$group == FALSE,], label.times);
	at.risk.high <- BoutrosLab.statistics.survival::calculate.number.at.risk(survobj[signature.results$group == TRUE, ], label.times);

	# extract survival statistics and create a survival plot
	survivaldata <- BoutrosLab.plotting.survival::create.km.plot(
		survival.object = survobj,
		patient.groups  = factor(x = as.character(signature.results$group), levels = c('FALSE','TRUE'), labels = c('Low','High')),
		filename        = filename,
		xlab.label      = xlab.label,
		ylab.label      = ylab.label,
		xlimits         = xlimits,
		xat             = xat,
		main            = main,
		main.cex        = main.cex,
		# label.positions = label.positions,
		# label.cex = label.cex,
		# at.risk.low = at.risk.low,
		# at.risk.high = at.risk.high,
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
	return(survivaldata);

	}

