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

make.survfit.plot <- function(x, xlimits, xat, type = 1, main = NULL, main.cex = 3.0, line.colours = c("blue", "red"), enable.warnings = FALSE,...) {

	# create groups
	g <- with(x, rep(names(strata), strata));

	# we need to handle the origin, so we'll create a new object and add dummy points
	x <- data.frame(
		strata  = c(g, unique(g)),
		surv    = c(x$surv, rep(1, length(unique(g)))),
		time    = c(x$time, rep(0, length(unique(g)))),
		n.event = c(x$n.event, rep(1, length(unique(g))))
		);

	# set parameters for plotting lines
	if (1 == type) { subset.to.plot = (x$n.event > -1); }

	# set parameters for plotting points
	if (2 == type) { subset.to.plot = (x$n.event == 0); }
        
	str <- "%s = %s";

	# check if graphics device is postscript
	if ("postscript" %in% rownames(as.matrix(dev.cur()))) {
		ps.options(family = "sans");
		}

	# check if graphics device is not set i-e "null device"
	if (enable.warnings == TRUE && dev.cur() == 1) {
		warning("If you wish to print this plot to postscript device, please set family param as: postscript(family=\"sans\")\n");
		}
        
	# make and return the plot
	lattice::xyplot(
		surv ~ time,
		data = x,
		panel = function(...) {
			if (1 == type) { panel.xyplot(type = "s", cex = 1, ...); }
			if (2 == type) { panel.xyplot(type = "p", cex = 1.1, pch = "|",...); }
			},
		groups = as.character(x$strata),
		subset = subset.to.plot,
		main = BoutrosLab.plotting.general::get.defaults(
			property = "fontfamily", 
			add.to.list = list(
				label = main,
				fontface = "bold",
				cex = main.cex
				)
			),
		xlab = BoutrosLab.plotting.general::get.defaults(
			property = "fontfamily", 
			add.to.list = list(
				label = "Time (Years)",
				fontface = "bold",
				cex = 2.75
				)
			),
		ylab = BoutrosLab.plotting.general::get.defaults(
			property = "fontfamily", 
			add.to.list = list(
				label = "Fraction of Cohort",
				fontface = "bold",
				cex = 2.75
				)
			),
		scales = list(
			x = list(
				limits = xlimits,
				at = xat,
				tck = c(1,1),
				axs = "r",
				cex = 2
				),
			y = list(
				limits = c(0,1.05),
				at = seq(0,1,0.25),
				tck = c(0.75,0.75),
				axs = "r",
				cex = 2
				)
			),
		par.settings = list(
			axis.line = list(
				lwd = 2
				),
			layout.heights = list(
				top.padding = 0.1,
				main = if (is.null(main)) { 0.1 } else { 1 },
				main.key.padding = 0.1,
				key.top = 0.1,
				key.axis.padding = 0.1,
				axis.top = 0.5,
				axis.bottom = 1.0,
				axis.xlab.padding = 0.7,
				xlab = 1,
				xlab.key.padding = 0,
				key.bottom = 0,
				key.sub.padding = 0,
				sub = 0,
				bottom.padding = 1
				),
			layout.widths = list(
				left.padding = 0,
				key.left = 0,
				key.ylab.padding = 0,
				ylab = 0.4,
				ylab.axis.padding = 1,
				axis.left = 1,
				axis.right = 1,
				axis.key.padding = 0.1,
				key.right = 0.1,
				right.padding = 0.2
				)
			),
		lwd = 3,
		col = line.colours,
		...
		);

	}
