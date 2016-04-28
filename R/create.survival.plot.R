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

### NOTE ########################################################################################
# This function is deprecated.  Use create.km.plot instead.

### create.survival.plot.R ######################################################################
# Description:
#       Produces a single Kaplan-Meier (KM) plot for any number of risk groups, as well as a
#       corresponding risk table.
# Input variables:
#       
# Output variables: If no filename is specified returns a trellis object.  If a filename is
#       specified, prints the trellis object to a .tiff file
#       
create.survival.plot <- function(
	survival.object, 
	patient.groups = NA, 
	filename = NULL, 
	xat = NA, 
	yat = seq(0,1,0.2),
	xlimits = NA, 
	ylimits = c(0,1.03),
	xaxis.cex = 2,
	yaxis.cex = 2,
	xlab.label = 'Time (Months)',
	xlab.cex = 2.75,
	ylab.label = 'Proportion',
	ylab.cex = 2.75,
	risk.labels = NA,
	key.groups.labels = levels(as.factor(patient.groups)),
	risk.label.pos = NA,
	key.groups.title = NULL, 
	explicit.HR.label = TRUE,
	main = NULL, 
	main.cex = 3.0,
	covariates = NULL,
	stratification.factor = NULL,
	lwd = 2,
	censoring.pch.cex = 1.1,
	digits = 2, 
	line.colours = NA, 
	statistical.method = NA, 
	ph.assumption.check = "warning",
	cox.zph.threshold = 0.1,
	show.key.groups = NA, 
	show.risktable = TRUE, 
	key.groups.corner = c(0,0), 
	key.groups.x.pos = 0, 
	key.groups.y.pos = 0.01, 
	key.stats.corner = c(1,0), 
	key.stats.x.pos = 1, 
	key.stats.y.pos = 0.01, 
	ylab.axis.padding = 1,
	bottom.padding = 0.7,  
	top.padding = 0.1,
	right.padding = 0.1,
	left.padding = 0.5,
	return.statistics = FALSE,
	height = 7, 
	width = 7, 
	resolution = 1000, 
	size.units = 'in', 
	enable.warnings = TRUE,
	description = NULL
	) { 

	create.km.plot <- function(
		survival.object, 
		patient.groups, 
		filename, 
		xat, 
		yat,
		xlimits, 
		ylimits,
		xaxis.cex,
		yaxis.cex,
		xlab.label,
		xlab.cex,
		ylab.label,
		ylab.cex,
		risk.labels,
		key.groups.labels,
		risk.label.pos,
		key.groups.title, 
		explicit.HR.label,
		main, 
		main.cex,
		covariates,
		stratification.factor,
		lwd,
		censoring.pch.cex = 1.1,
		digits, 
		line.colours, 
		statistical.method, 
		ph.assumption.check,
		cox.zph.threshold,
		show.key.groups, 
		show.risktable, 
		key.groups.corner, 
		key.groups.x.pos, 
		key.groups.y.pos, 
		key.stats.corner, 
		key.stats.x.pos, 
		key.stats.y.pos, 
		ylab.axis.padding,
		bottom.padding,  
		top.padding,
		right.padding,
		left.padding,
		return.statistics,
		height, 
		width, 
		resolution, 
		size.units, 
		enable.warnings,
		description
		)

		warning("This function is deprecated.  Use create.km.plot instead.");

	}
