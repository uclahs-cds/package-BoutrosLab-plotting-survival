# Clone of survival:::plot.cox.zph but for BPG
plot.cox.zph.bpg <- function(
    x,
    resid=TRUE,
    se=TRUE,
    df=4,
    nsmo=40,
    var,
    lty=1:2,
    col=1,
    lwd=1,
    hr = FALSE,
    xlab.label = 'Time (transformation)',
    ylab.label = expression(bold(beta~'(t) for <var>')),
    ...
    ) {
    xx <- x$x
    yy <- x$y
    df <- max(df)     # in case df is a vector
    nvar <- ncol(yy)
    pred.x <- seq(from=min(xx), to=max(xx), length=nsmo)
    temp <- c(pred.x, xx)
    lmat <- splines::ns(temp, df=df, intercept=TRUE)
    pmat <- lmat[1:nsmo,]       # for prediction
    xmat <- lmat[-(1:nsmo),]
    if (!is.logical(hr)) stop('hr parameter must be TRUE/FALSE')

    if (missing(ylab.label)) {
        if (hr) {
            ylab.label <- paste0('HR(t) for ', dimnames(yy)[[2]])
        } else {
            ylab.label <- paste0('\u03b2(t) for ', dimnames(yy)[[2]])
        }
    }
    if (missing(xlab.label)) {
        xlab.label <- sprintf('Time (%s transform)', x$transform);
        }

    if (missing(var)) {
        var <- 1:nvar
        }
    else {
        if (is.character(var)) var <- match(var, dimnames(yy)[[2]])
        if  (any(is.na(var)) || max(var)>nvar || min(var) <1)
            stop('Invalid variable requested')
        }

    # Figure out a 'good' set of x-axis labels.  Find 8 equally spaced
    #    values on the 'transformed' axis.  Then adjust until they correspond
    #    to rounded 'true time' values.  Avoid the edges of the x axis, or
    #    approx() may give a missing value
    if (x$transform == 'log') {
        xx <- exp(xx)
        pred.x <- exp(pred.x)
        }
    else if (x$transform != 'identity') {
        xtime <- x$time;
        indx <- !duplicated(xx);  #avoid a warning message in R
        apr1  <- approx(
            x = xx[indx],
            y = xtime[indx],
            xout = seq(min(xx), max(xx), length=17)[2*(1:8)]
            );
        temp <- signif(apr1$y,2)
        apr2  <- approx(xtime[indx], xx[indx], temp)
        xaxisval <- apr2$y
        xaxislab <- rep('',8)
        for (i in 1:8) xaxislab[i] <- format(temp[i])
        }
    col <- rep(col, length=2)
    lwd <- rep(lwd, length=2)
    lty <- rep(lty, length=2)

    # Now, finally do the work
    rtn.plots <- sapply(var, function(i) {
        #   Since release 3.1-6, yy can have missing values.  If a covariate is
        # constant within a stratum then it's Shoenfeld residual is identially
        # zero for all observations in that stratum.  These 'structural zeros'
        # are marked with an NA.  They contain no information and should not
        # by plotted.  Thus we need to do the spline fit one stratum at a time.
        y <- yy[,i]
        keep <- !is.na(y)
        if (!all(keep)) y <- y[keep]

        qmat <- qr(xmat[keep,])
        if (qmat$rank < df) {
            warning('spline fit is singular, variable ', i, ' skipped')
            next
        }

        yhat <- pmat %*% qr.coef(qmat, y)
        if (resid) yr <-range(yhat, y)
        else       yr <-range(yhat)

        if (se) {
            bk <- backsolve(qmat$qr[1:df, 1:df], diag(df))
            xtx <- bk %*% t(bk)
            seval <- ((pmat%*% xtx) *pmat) %*% rep(1, df)

            temp <- 2* sqrt(x$var[i,i]*seval)
            yup <- yhat + temp
            ylow<- yhat - temp
            yr <- range(yr, yup, ylow)
            }

            plotpar <- list(...);

            yaxis.log <- FALSE;
            if (hr) {
                y <- exp(y);
                yr <- exp(yr);
                yhat <- exp(yhat);
                yup <- exp(yup);
                ylow <- exp(ylow);
                # Not sure I understand fully the log y-scale for BPG or lattice
                # This *should* be
                # yaxis.log <- exp(1);
            }

            # x-axis log works fine
            xaxis.log <- if(x$transform == 'log') exp(1) else FALSE;

            if (x$transform=='identity') {
                plot.data <- data.frame(
                    x = range(xx),
                    y = range(yr)
                    )
                zph.plot <- create.scatterplot(
                    formula = y ~ x,
                    data = plot.data,
                    type = 'n',
                    ylab.label = ylab.label,
                    xlab.label = xlab.label,
                    yaxis.log = yaxis.log,
                    ...
                    )
            } else {
                plot.data <- data.frame(
                    x = range(xx[keep]),
                    y = yr
                    )
                if (is.null(plotpar$xat) && is.null(plotpar$xaxis.lab)) {
                    zph.plot <- create.scatterplot(
                        formula = y ~ x,
                        data = plot.data,
                        type = 'n',
                        xaxis.log = xaxis.log,
                        yaxis.log = yaxis.log,
                        ylab.label = ylab.label,
                        xlab.label = xlab.label,
                        xat = xaxisval,
                        xaxis.lab = xaxislab,
                        ...
                        )
                    } else {
                    zph.plot <- create.scatterplot(
                        formula = y ~ x,
                        data = plot.data,
                        type = 'n',
                        xaxis.log = xaxis.log,
                        yaxis.log = yaxis.log,
                        ylab.label = ylab.label,
                        xlab.label = xlab.label,
                        ...
                        )
                    }
                }

                if (resid) {
                    points.data <- data.frame(
                        x = xx[keep],
                        y = y
                        )
                    zph.plot <- zph.plot + create.scatterplot(
                        formula = y ~ x,
                        data = points.data,
                        xaxis.log = xaxis.log,
                        yaxis.log = yaxis.log,
                        ...
                        )
                    }

                pred.data <- data.frame(
                    y = yhat,
                    yupper = yup,
                    ylower = ylow,
                    x = pred.x
                    )

                zph.plot <- zph.plot + create.scatterplot(
                    formula = y ~ x,
                    data = pred.data,
                    type = 'l',
                    xaxis.log = xaxis.log,
                    yaxis.log = yaxis.log,
                    lty=lty[1],
                    col=col[1],
                    lwd=lwd[1]
                    )

                if (se) {
                    zph.plot <- zph.plot + create.scatterplot(
                        formula = yupper ~ x,
                        data = pred.data,
                        type = 'l',
                        xaxis.log = xaxis.log,
                        yaxis.log = yaxis.log,
                        lty=lty[2],
                        col=col[2],
                        lwd=lwd[2]
                        )
                    zph.plot <- zph.plot + create.scatterplot(
                        formula = ylower ~ x,
                        data = pred.data,
                        type = 'l',
                        xaxis.log = xaxis.log,
                        yaxis.log = yaxis.log,
                        lty=lty[2],
                        col=col[2],
                        lwd=lwd[2]
                        )
                    }
        return(zph.plot);
    }, simplify = FALSE, USE.NAMES = TRUE)

    if (length(rtn.plots) == 1) {
        return(rtn.plots[[1]]);
        } else {
        return(rtn.plots);
        }
    }
