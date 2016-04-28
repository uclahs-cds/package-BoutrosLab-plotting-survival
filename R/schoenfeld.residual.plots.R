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

schoenfeld.residual.plots <- function(cox.model, filename) {
	
	cox.zph.test <- cox.zph(cox.model, transform = "identity");

	# Get number of betas estimated in the cox model
	if ("GLOBAL" %in% row.names(cox.zph.test$table)) {
		number.of.betas <- nrow(cox.zph.test$table) - 1;
		row.number.for.GLOBAL <- nrow(cox.zph.test$table);
		}
	else {
		number.of.betas <- nrow(cox.zph.test$table);
		row.number.for.GLOBAL <- NA;
		}

	# Get vector of all pvalues computed for cox.zph test, one for each beta and one for the overall model (global)
	pvalues.zph <- cox.zph.test$table[,'p'];

	# Generate one Schoenfeld residual plot for each covariate in the cox model, and print global and covariate-specific
	# cox.zph pvalues on each plot.
	if (is.null(filename)){
		stop("The residual plot(s) can only be produced and returned if a filename is specified.")
		}
	else {
		for (i in 1:number.of.betas) {
			# Generate filename for the i-th residual plot
			residual.plot.filename <- BoutrosLab.utilities::generate.filename(
				project.stem = paste(gsub('.tiff', '', filename), 'Schoenfeld.residuals', sep = "_"), 
				file.core = names(pvalues.zph)[i],
				extension = "tiff",
				file.date = FALSE
				);			
			#tiff(paste(gsub('.tiff', '', filename), '_Schoenfeld.residuals.', names(pvalues.zph)[i], '.tiff', sep=''));

			# Open tiff device, produce residual plot, annotate it with relevant p-values, and close the device

			# set the graphics driver
			current.type <- getOption("bitmapType");
			options(bitmapType = "cairo");

			tiff(residual.plot.filename);			
			plot(cox.zph.test[i]);
			title(
				main = ifelse(
					!is.na(row.number.for.GLOBAL),
					yes = BoutrosLab.plotting.general::display.statistical.result(cox.zph.test$table[row.number.for.GLOBAL,3], statistic.type="GLOBAL, P"),
					no = ""
					),
				sub = BoutrosLab.plotting.general::display.statistical.result(cox.zph.test$table[i,3], statistic.type=paste(row.names(cox.zph.test$table)[i], ", P", sep=" "))
				);

			# add a line for the coefficient from the Cox model
			abline(h=summary(cox.model)$coef[i], lty=2, col='red');

			# finish plot
			dev.off();
			options(bitmapType = current.type);
			}
		}
	}
