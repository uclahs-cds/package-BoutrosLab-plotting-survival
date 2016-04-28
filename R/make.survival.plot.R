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

make.survival.plot <- function(survival.object, patient.groups, filename, xlimits, xat, label.positions, label.cex = 1.2, at.risk.low, at.risk.high, risk.labels = c("Low", "High"), main = NULL, main.cex = 3.0, rounding = 2, covariates = NULL, line.colours = c("blue", "red"), key = NULL, statistical.method = NA, x.statistic = 0.200, y.statistic = 0.610, resolution = 2000, height = 7.007874, width = 7.007874, size.units = 'in', enable.warnings = FALSE) {

	# warn the user that this function is being deprecated
	warning('This function will be deprecated -- consider using create.survival.plot.R() instead');

	# Set up parameters about the display
	if ('in' == size.units) {
		xpixels <- resolution * width;
		ypixels <- resolution * height;
		}
	else {
		return( rep(NA,4) );
		}

	# verify parameters
	if (!is.na(statistical.method) & (!tolower(statistical.method) %in% c('cox', 'logrank', 'ttest'))) {
		stop('Invalid statistical method selected');
		}

	# set the graphics driver
	current.type <- getOption("bitmapType");
	options(bitmapType = "cairo");

	# open output file
	tiff(
		filename = filename,
		height = height,
		width = width,
		units = size.units,
		res = resolution,
		compression = "lzw"
		);

	# create first plot: to draw the curve only (without censorship marks)
	graphics::plot(
		x = BoutrosLab.plotting.survival::make.survfit.plot(
			x = survival::survfit(survival.object ~ patient.groups),
			xlimits = xlimits,
			xat = xat,
			type = 1,
			line.colours = line.colours,
			main = main,
			main.cex = main.cex,
			key = key
			),
		position = c(0.03,0.2,1.0,1.0),
		newpage = TRUE
		);

	# add the second plot (i.e. censorship marks) if necessary (i.e. if there are any!)
	if (any(survival.object[, "status"] == 0)) {
		plot(
			x = BoutrosLab.plotting.survival::make.survfit.plot(
				x = survival::survfit(survival.object ~ patient.groups),
				xlimits = xlimits,
				xat = xat,
				type = 2,
				line.colours = line.colours,
				main = main,
				main.cex = main.cex, 
				),
			position = c(0.03,0.2,1.0,1.0),
			newpage = FALSE
			);
		}

	# if the user hasn't set it, decide what statistical analysis to use
	if (is.na(statistical.method)) {

		# if no events, must use a t-test
		if (all(survival.object[, "status"] == 1)) { statistical.method <- 'ttest'; }

		# if no censored, must use a log-rank test
		else if (all(survival.object[, "status"] == 0)) { statistical.method <- 'logrank'; }

		# otherwise try a Cox model
		else { statistical.method <- 'cox'; }

		}

	# do Cox-modeling analysis
	if (statistical.method == 'cox') {

		# fit a Cox model
		stats <- BoutrosLab.statistics.survival::fit.coxmodel(
			groups = patient.groups,
			survobj = survival.object,
			rounding = rounding,
			other.data = covariates
			);

		# define statistical parameters
		this.hr   <- stats[1];
		this.95l  <- stats[2];
		this.95u  <- stats[3];
		this.pval <- stats[4];

		# add survival statistics: labels
		ltext(x = x.statistic * xpixels, y = (y.statistic + 0.00) * ypixels, cex = 1.55, adj = c(0,0), label = 'HR:');
		ltext(x = x.statistic * xpixels, y = (y.statistic + 0.05) * ypixels, cex = 1.55, adj = c(0,0), label = 'P:');

		# add survival statistics: values
		ltext(x = (x.statistic + 0.07) * xpixels, y = y.statistic * ypixels, cex = 1.55, adj = c(0,0), label = paste(sprintf(paste("%.", 2, "f", sep = ""), this.hr)," (", sprintf(paste("%.", 2, "f", sep = ""), this.95l), " - ", sprintf(paste("%.", 2, "f", sep = ""), this.95u), ")", sep = ""));
		if (this.pval > 0.01)                      { ltext(x = (x.statistic + 0.07) * xpixels, y = (y.statistic + 0.05) * ypixels, cex = 1.55, adj = c(0,0), label = sprintf(paste("%.", 2, "f", sep = ""), this.pval)); }
		if (this.pval <= 0.01 & this.pval > 0.001) { ltext(x = (x.statistic + 0.07) * xpixels, y = (y.statistic + 0.05) * ypixels, cex = 1.55, adj = c(0,0), label = sprintf(paste("%.", 3, "f", sep = ""), this.pval)); }
		if (this.pval <= 0.001)                    { ltext(x = (x.statistic + 0.07) * xpixels, y = (y.statistic + 0.05) * ypixels, cex = 1.55, adj = c(0,0), label = scientific.notation(x = this.pval, digits = 2)); } 

		ret.stats <- c(this.pval, this.hr, this.95l, this.95u);

		}

	# do t-test analysis
	else if (statistical.method == 'ttest') {

		# perform t-test to assess survival differences
		stats <- t.test(survival.object[,1] ~ patient.groups);
                
		# define statistical parameters
		this.pval <- stats$p.value;
		average.group1 <- stats$estimate[1];
		average.group2 <- stats$estimate[2];

		# add t-test statistics: labels
		ltext(x = x.statistic * xpixels, y = (y.statistic + 0.05) * ypixels, cex = 1.55, adj = c(0,0), label = 'P:');

		# add survival statistics: values
		if (this.pval > 0.01)                      { ltext(x = (x.statistic + 0.07) * xpixels, y = (y.statistic + 0.05) * ypixels, cex = 1.55, adj = c(0,0), label = sprintf(paste("%.", 2, "f", sep = ""), this.pval)); }
		if (this.pval <= 0.01 & this.pval > 0.001) { ltext(x = (x.statistic + 0.07) * xpixels, y = (y.statistic + 0.05) * ypixels, cex = 1.55, adj = c(0,0), label = sprintf(paste("%.", 3, "f", sep = ""), this.pval)); }
		if (this.pval <= 0.001)                    { ltext(x = (x.statistic + 0.07) * xpixels, y = (y.statistic + 0.05) * ypixels, cex = 1.55, adj = c(0,0), label = scientific.notation(x = this.pval, digits = 2)); }

		ret.stats <- c('t-test', this.pval, average.group1, average.group2);

		}

	# do log-rank analysis
	else if (statistical.method == 'logrank') {

		# perform t-test to assess survival differences
		stats <- survdiff(survival.object ~ patient.groups);

		# define statistical parameters
		this.pval <- pchisq(stats$chisq, df = 1, lower.tail = FALSE);
		expected.group1 <- stats$exp[1];
		expected.group2 <- stats$exp[2];

		# add log-rank statistics: label
		ltext(x = x.statistic * xpixels, y = (y.statistic + 0.05) * ypixels, cex = 1.55, adj = c(0,0), label = 'P:');

		# add survival statistics: values
		if (this.pval > 0.01)                      { ltext(x = (x.statistic + 0.07) * xpixels, y = (y.statistic + 0.05) * ypixels, cex = 1.55, adj = c(0,0), label = sprintf(paste("%.", 2, "f", sep = ""), this.pval)); }
		if (this.pval <= 0.01 & this.pval > 0.001) { ltext(x = (x.statistic + 0.07) * xpixels, y = (y.statistic + 0.05) * ypixels, cex = 1.55, adj = c(0,0), label = sprintf(paste("%.", 3, "f", sep = ""), this.pval)); }
		if (this.pval <= 0.001)                    { ltext(x = (x.statistic + 0.07) * xpixels, y = (y.statistic + 0.05) * ypixels, cex = 1.55, adj = c(0,0), label = scientific.notation(x = this.pval, digits = 2)); }

		ret.stats <- c('log-rank', this.pval, expected.group1, expected.group2);

		}

	# this should never happen -- weird settings should have been caught earlier
	else {
		stop('Invalid statistical-method selected');
		}

	# add at-risk annotations: labels
	ltext(x = 0.005 * xpixels, y = 0.85 * ypixels, cex = label.cex + 0.2, label = "Number at Risk", fontfamily = BoutrosLab.plotting.general::get.defaults(property = "fontfamily"), fontface = "bold", adj = c(0, 0.5));
	ltext(x = 0.005 * xpixels, y = 0.90 * ypixels, cex = label.cex, label = risk.labels[1],   fontfamily = BoutrosLab.plotting.general::get.defaults(property = "fontfamily"), fontface = "bold", adj = c(0, 0.5));
	ltext(x = 0.005 * xpixels, y = 0.96 * ypixels, cex = label.cex, label = risk.labels[2],   fontfamily = BoutrosLab.plotting.general::get.defaults(property = "fontfamily"), fontface = "bold", adj = c(0, 0.5));
        
	# add at-risk annotations: values
	for (i in 1:length(label.positions)) {
		ltext(x = label.positions[i] * xpixels, y = 0.90 * ypixels, cex = label.cex - 0.1, fontfamily = BoutrosLab.plotting.general::get.defaults(property = "fontfamily"), fontface = "bold", label = at.risk.low[i]);
		ltext(x = label.positions[i] * xpixels, y = 0.96 * ypixels, cex = label.cex - 0.1, fontfamily = BoutrosLab.plotting.general::get.defaults(property = "fontfamily"), fontface = "bold", label = at.risk.high[i]);
		}

	# finish plot
	dev.off();
	options(bitmapType = current.type);

	# return statistics back to the caller
	return(ret.stats);

	}
