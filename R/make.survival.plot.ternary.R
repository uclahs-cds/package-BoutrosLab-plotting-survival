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

make.survival.plot.ternary <- function(survival.object, patient.groups, filename, xlimits, xat, label.positions, label.cex = 1.2, at.risk.low, at.risk.medium, at.risk.high, risk.labels = c("Low", "Medium", "High"), main = NULL, main.cex = 3.0, rounding = 2, covariates = NULL, line.colours = c("blue", "black", "red"), key = NULL, resolution = 2000, height = 7.007874, width = 7.007874, size.units = 'in', enable.warnings = FALSE) {

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

	# add at-risk annotations: labels
	ltext(x = 0.005 * xpixels, y = 0.84 * ypixels, cex = label.cex + 0.2, label = "Number at Risk", fontfamily = get.defaults(property = "fontfamily"), fontface = "bold", adj = c(0, 0.5));
	ltext(x = 0.005 * xpixels, y = 0.89 * ypixels, cex = label.cex, label = risk.labels[1],   fontfamily = get.defaults(property = "fontfamily"), fontface = "bold", adj = c(0, 0.5));
	ltext(x = 0.005 * xpixels, y = 0.93 * ypixels, cex = label.cex, label = risk.labels[2],   fontfamily = get.defaults(property = "fontfamily"), fontface = "bold", adj = c(0, 0.5));
	ltext(x = 0.005 * xpixels, y = 0.97 * ypixels, cex = label.cex, label = risk.labels[3],   fontfamily = get.defaults(property = "fontfamily"), fontface = "bold", adj = c(0, 0.5));

	# add at-risk annotations: values
	for (i in 1:length(label.positions)) {
		ltext(x = label.positions[i] * xpixels, y = 0.89 * ypixels, cex = label.cex - 0.1, fontfamily = get.defaults(property = "fontfamily"), fontface = "bold", label = at.risk.low[i]);
		ltext(x = label.positions[i] * xpixels, y = 0.93 * ypixels, cex = label.cex - 0.1, fontfamily = get.defaults(property = "fontfamily"), fontface = "bold", label = at.risk.medium[i]);
		ltext(x = label.positions[i] * xpixels, y = 0.97 * ypixels, cex = label.cex - 0.1, fontfamily = get.defaults(property = "fontfamily"), fontface = "bold", label = at.risk.high[i]);
		}

	# finish plot
	dev.off();
	options(bitmapType = current.type);

	}
