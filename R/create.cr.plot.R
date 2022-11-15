### create.cr.plot.R ######################################################################
# Description:
#       Produces a competing risk (CR) plot for any number of risk groups, as well as a
#       corresponding risk table.
# Input variables:
#   data: a data.frame contains `time`, `status`, and `group` information
#   formula: time ~ status
#            time: numeric values, representing months or years;
#            status: numeric values, such as 0, 1, 2, 3 ..., representing different groups
#   patient.groups: the colname of group information

# Output variables: If no filename is specified returns a trellis object.  If a filename is
#       specified, prints the trellis object to a plot.pdf file

create.cr.plot <- function(
    data,
    formula,
    cencode = 0,
    patient.groups = NA,
    filename = NULL,
    xlab.label = 'Time',
    ylab.label = 'Cumulative Incidence Function',
    return.statistics = FALSE,
    xaxis.cex = 1,
    yaxix.cex = 1,
    line.colours = NULL,
    height = 7,
    width = 7,
    style = 'BoutrosLab',
    resolution = 1000
    ) {

    # create CIF curve estimates
    fit.time <- data[,toString(formula[[2]])];
    fit.status <- data[,toString(formula[[3]])];
    fit.group <- eval(substitute(patient.groups), data, parent.frame());

    if(cencode == 0 | is.null(cencode)){
        fit.cencode <- 0
    } else {
        cencode <- fit.cencode
    }
    # create cumulative incidence analysis
    fit <- cuminc(
        ftime = fit.time,
        fstatus = fit.status,
        group = fit.group,
        cencode = cencode
    );
    print('Here is the cumulative incidence analysis result:');
    print(fit);

    # extract the list to a data frame
    df <- data.frame();
    suppressMessages(
        for (i in 1:(length(fit) - 1)) {
        simplify.fit <- data.frame(simplify2array(fit[[i]]));
        simplify.fit$type <- names(fit)[i];
        if (dim(df)[1] == 0 ){
            df <-  simplify.fit;
            rm(simplify.fit);
        } else {
            df <-  rbind(df, simplify.fit);
            rm(simplify.fit);
        }
        }
    )
    colnames(df) <- c('Time', 'CIF', 'Variance', 'Group');

    # create competing risks plot
    if (is.null(filename)) {
        filename <- 'plot.pdf';
    }

    # set line colours
    if (is.null(line.colours)) {
        line.colours = default.colours(length(unique(fit.group)), palette.type = 'survival');
    }
    cr.plot <- BoutrosLab.plotting.general::create.scatterplot(
        filename = filename,
        formula =  CIF ~ Time,
        data = df,
        xaxix.tck = c(1,0),
        yaxis.tck = c(1,0),
        groups = Group,
        type = c('l','p'),
        col = as.vector(rbind(line.colours, line.colours)),
        pch = '|',
        lty = as.numeric(gsub('.*[ ](.*)', '\\1', unique(df$Group))),
        lwd = 5,
        right.padding = 2,
        left.padding = 2,
        top.padding = 2,
        bottom.padding = 2,
        width = width,
        height = height,
        style = style,
        key = list(
        text = list(
            lab = unique(df$Group),
            cex = 1,
            col = 'black'
            ),
        lines = list(
            lty = as.numeric(gsub('.*[ ](.*)', '\\1', unique(df$Group))),
            col = as.vector(cbind(line.colours, line.colours)),
            lwd = 5
            ),
        x = 0.04,
        y = 0.95,
        padding.text = 2
        );
    )
    return(cr.plot);
    }
