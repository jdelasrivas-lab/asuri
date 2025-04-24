#' Extra functions to graphically represent the results obtained by
#' other functions in the package
#'
# Plot geneSurv() results -------------------------------------------------
plotBoxplot <- function(gs_result) {
    patientExpr <- gs_result$patientExpr
    patientClass <- gs_result$patientClass
    pvalue <- gs_result$wilcox.pvalue
    geneName <- gs_result$geneName
    plot_values <- gs_result$plot_values

    if (plot_values$source == "geneSurv-risk") {
        stop("This function cannot be used when c(type = 'risk').")
    }


    par(cex.lab = 1.5) # is for y-axis
    par(cex.axis = 1.5) # is for x-axis
    boxplot(patientExpr[patientClass == 1],
        patientExpr[patientClass == 2],
        col = c(3, 2),
        ylab = "expression level",
        xlab = "Wilcoxon test p.value for groups:",
        sub = format(pvalue$p.value),
        names = c("LowExpr", "HighExpr")
    )

    title(main = paste0("Expression plot (", geneName, ")"))
    legend("bottomright",
        title = "Nº of samples",
        c(
            as.character(sum(patientClass == 1)),
            as.character(sum(patientClass == 2))
        ),
        fill = c(3, 2),
        horiz = FALSE
    )
}

plotProbClass <- function(gs_result) {
    patientExpr <- gs_result$patientExpr
    patientClass <- gs_result$patientClass
    patientClassProbality <- gs_result$patientClassProbality
    geneName <- gs_result$geneName

    plot_values <- gs_result$plot_values

    if (plot_values$source == "geneSurv-risk") {
        stop("This function cannot be used when c(type = 'risk').")
    }

    colVector <- patientClass[order(patientExpr)]
    colVector <- replace(colVector, colVector == 1, 3)
    patientClassprobs <- patientClassProbality
    plot(patientClassprobs[order(patientExpr)],
        col  = colVector,
        ylab = "Class Probability",
        xlab = "Patients ordered by expression (coloured by groups)"
    )
    title(main = paste0("Patient Class probability (", geneName, ")"))
}

# Plot patientRisk() result -----------------------------------------------
plotLogRank <- function(pr_result) {
    plot_values <- pr_result$plot_values

    par(cex.lab = 1.5) # is for y-axis
    par(cex.axis = 1.5) # is for x-axis

    p.val30a70 <- plot_values$sigmoid_logrank$p.val30a70
    lower10 <- plot_values$sigmoid_logrank$lower10
    lowIndexReal <- plot_values$sigmoid_logrank$lowIndexReal
    highIndexReal <- plot_values$sigmoid_logrank$highIndexReal
    cutPoint <- plot_values$sigmoid_logrank$cutPoint
    nSamples <- plot_values$sigmoid_logrank$nSamples

    lowIndex0.30 <- round(0.3 * nSamples)
    highIndex0.70 <- round(0.7 * nSamples)
    ylimits <- c(0, summary(p.val30a70[lowIndex0.30:highIndex0.70])[3])
    graphics::plot(p.val30a70,
        type = "b",
        pch  = 18,
        col  = "black",
        lty  = 1,
        lwd  = 1,
        ylim = ylimits,
        xlab = "Patients ordered by Risk",
        ylab = "Logrank p.value",
        xpd  = FALSE
    )
    # horizontal black line
    abline(h = lower10, col = "black", lwd = 2, xpd = FALSE) 
    # lower threshold blue
    abline(v = lowIndexReal, col = 4, lwd = 2, xpd = FALSE) 
    # upper threshold red
    abline(v = highIndexReal, col = 2, lwd = 2, xpd = FALSE) 
    # green cutpoint
    abline(v = cutPoint, col = 3, lwd = 2, xpd = FALSE) 
}

plotSigmoid <- function(pr_result) {
    plot_values <- pr_result$plot_values

    lowIndexReal <- plot_values$sigmoid_logrank$lowIndexReal
    highIndexReal <- plot_values$sigmoid_logrank$highIndexReal
    cutPoint <- plot_values$sigmoid_logrank$cutPoint
    orderednormalized_risk <- plot_values$sigmoid_logrank$orderednormalized_risk

    colores <- rep(0, length(orderednormalized_risk))
    colores[seq(1, lowIndexReal)] <- 4 # Low risk
    colores[seq(lowIndexReal + 1, (cutPoint - 1))] <- "skyblue" # Low-like risk
    colores[cutPoint] <- 3
    colores[seq((cutPoint + 1), (highIndexReal - 1))] <- "pink" # High-like risk
    colores[seq(highIndexReal, length(orderednormalized_risk))] <- 2 # High risk

    par(xpd = FALSE)
    graphics::plot(orderednormalized_risk,
        col  = colores,
        xlab = "Patients ordered by Risk",
        ylab = "Risk Score"
    )
    grid(nx = NA, ny = NULL, lty = 2, col = "gray", lwd = 1)

    abline(v = lowIndexReal, col = 4, lwd = 2, xpd = FALSE)
    abline(v = highIndexReal, col = 2, lwd = 2, xpd = FALSE)
    abline(v = cutPoint, col = 3, lwd = 2, xpd = FALSE)
}

plotLambda <- function(pr_result) {
    plot_values <- pr_result$plot_values

    nested.p.vals <- plot_values$lambda$nested.p_vals
    number.features <- plot_values$lambda$number.features

    pvalsByLambda <- apply(nested.p.vals$p.vals, 2, FUN = function(x) min(x))
    ylimits <- c(0, summary(pvalsByLambda)[3])

    names(pvalsByLambda) <- number.features
    graphics::plot(pvalsByLambda,
        xlab = "Number of genes by lambda",
        ylab = "Logrank p.value",
        xaxt = "n",
        type = "b"
    )

    axis(1, at = seq(1, length(number.features)), 
         labels = number.features, cex.axis = 0.8)
}

plotBetas <- function(pr_result) {
    betasPlot <- pr_result$betasplot

    top25betasPlot <- as.data.frame(betasPlot[betasPlot$beta_value != 0, ])
    top25betasPlot$gene <- factor(top25betasPlot$gene, 
                                  levels = rev(top25betasPlot$gene))

    if (nrow(top25betasPlot) >= 25) {
        ngenes <- 25
        top25betasPlot <- head(top25betasPlot, 25)
    } else {
        ngenes <- nrow(top25betasPlot)
    }
    level_to_number <- c("***" = 4, "**" = 3, "*" = 2, "·" = 1, " " = 0)
    top25betasPlot$signif_numeric <- 
      as.integer(level_to_number[top25betasPlot$signif])
    top25betasPlot <- 
      top25betasPlot[order(top25betasPlot$signif_numeric, 
                           abs(top25betasPlot$beta_value), decreasing = TRUE), ]
    top25betasPlot$point_size <- -log10(top25betasPlot$p_value) * 1.1

    h <- ggplot(top25betasPlot, aes(
        x = abs(beta_value),
        y = gene,
        color = beta_group,
        size = point_size,
        shape = p_value < 0.06
    )) +
        scale_shape_manual(name = "p.value < 0.05", values = c(1, 19)) +
        ggtitle(paste0("Top ", ngenes, 
                       " genes with the greatest influence on risk")) +
        geom_point() +
        theme_light() +
        theme(
            legend.position = "right",
            panel.border    = element_blank(),
            axis.text.x     = element_text(size = 12, angle = 45, hjust = 1),
            axis.text.y     = element_text(size = 12)
        ) +
        xlab("Absolute Beta Values (Genes Risk Influence)") +
        ylab("") +
        scale_color_manual(
            name = "Beta Sign",
            breaks = c("type1", "type2"),
            labels = c("Increases risk (+)", "Decreases risk (-)"),
            values = c("type1" = "red", "type2" = 4)
        ) +
        # labs(size = "1 / p.value") +
        scale_size_continuous(
            name = "-log10(p.value)",
            breaks = -log10(c(0.05, 0.01, 0.001, 0.0001)) * 1.1, 
            labels = c("0.05", "0.01", "0.001", "0.0001")
        ) +
        guides(
            color = guide_legend(order = 1), 
            size = guide_legend(order = 2),
            shape = guide_legend(order = 3)
        )
    graphics::plot(h)
}

# Both geneSurv() and patientRisk() result --------------------------------
plotKM <- function(result,
                   col.surv = NULL,
                   col.ci = NULL,
                   par.bot = 6,
                   par.left = 10,
                   par.top = 2,
                   par.right = 4,
                   y.just.legend = .5,
                   x.title.adj = 1.75,
                   mark = 3,
                   simple = FALSE,
                   xaxis.at = c(0:10),
                   xaxis.lab = xaxis.at,
                   lty.surv = 1,
                   lty.ci = 3,
                   lwd = 4,
                   lwd.surv = 4,
                   lwd.ci = 4,
                   # lwd.grid = 6,
                   group.names = "",
                   group.order = seq(length(km$n)),
                   extra.left.margin = 4,
                   label.n.at.risk = FALSE,
                   draw.lines = TRUE,
                   cex.axis = 1.25,
                   main = "",
                   xlim = c(0, 10),
                   ylim = c(0, 1),
                   grid = TRUE,
                   lty.grid = 1,
                   lwd.grid = 1,
                   col.grid = grey(.9),
                   legend = TRUE,
                   loc.legend = "bottomleft",
                   add = FALSE,
                   ... # ... is passed to par()
) {
    par(cex.lab = 1.5) # is for y-axis
    par(cex.axis = 1.5) # is for x-axis
    
    plot_values <- result$plot_values
    nbetas <- nrow(result$betasplot)
    geneName <- result$geneName
    km <- plot_values$km$fitsKM
    p.value <- plot_values$km$p.val
    h.ratio <- plot_values$km$hazardR

    # Define colours based on source function
    if (is.null(col.surv)) {
        if (plot_values$source == "geneSurv-exprs") {
            col.surv <- c(3, 2)
            col.ci <- c(2, 3)
            group.names <- c("LowExpr", "HighExpr")
        } else if (plot_values$source == "patientRisk" || 
                   plot_values$source == "geneSurv-risk") {
            col.surv <- c(4, "red")
            col.ci <- c("red", 4)
            group.names <- c("LowRisk", "HighRisk")
        } else {
            col.surv <- c(3, 2)
            col.ci <- c(2, 3) # Default
        }
    }

    ng0 <- length(km$strata)
    ng <- max(ng0, 1)
    # When only one group...
    if (ng0 == 0) {
        km$strata <- length(km$time)
        names(km$strata) <- "All"
        legend <- draw.lines <- FALSE
    }

    lty.surv <- rep(lty.surv, ng)
    lwd.surv <- rep(lwd.surv, ng)
    col.surv <- rep(col.surv, ng)
    lty.ci <- rep(lty.ci, ng)
    lwd.ci <- rep(lwd.ci, ng)
    col.ci <- rep(col.ci, ng)

    ## group names and error checking
    gr <- c(km$strata)
    if (is.null(group.names) || length(group.names) != length(km$strata)) {
        group.names <- names(km$strata)
    }
    if (length(unique(group.names)) != ng) {
        stop("\n", "length(unique(group.names)) != number of groups.")
    }
    if (length(group.order) != length(seq_len(ng)) || 
        any(sort(group.order) != seq_len(ng))) {
        stop("\n", "Something wrong with group.order.", "\n", 
             "sort(group.order) must equal 1:", ng, ".")
    }
    # to remove unwanted white spaces in group.names.
    group.names <- gsub(" *$", "", group.names) 
    if (ng == 1 & (group.names[1] == "group.names")) {
        group.names <- "N at risk"
        label.n.at.risk <- FALSE
    }

    ## graphic parameters
    if (!add) {
        par(list(oma = c(1, 1, 1, 1), 
                 mar = c(4 + ng, 4 + extra.left.margin, 4, 2) + .1))
        if (simple) par(mar = c(3, 4, 2, 1) + .1)
        par(list(...))
    }

    ## reformat survival estimates
    dat <- data.frame(
        time = km$time, n.risk = km$n.risk, n.event = km$n.event, 
        survival = km$surv, std.err = km$std.err,
        lower = km$lower, upper = km$upper, group = rep(group.names, gr)
    )
    dat.list <- split(dat, f = dat$group)

    ## plot (but not survival curves)
    graphics::plot(0, type = "n", lwd = lwd, xlim = xlim, ylim = ylim, 
                   xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    if (grid) {
        par("xpd" = FALSE)
        abline(v = xaxis.at, lty = lty.grid, 
               lwd = lwd.grid, col = col.grid)
        abline(h = pretty(c(0, 1)), lty = lty.grid, 
               lwd = lwd.grid, col = col.grid)
    }
    axis(side = 2, at = pretty(c(0, 1)), cex.axis = cex.axis)
    axis(side = 1, at = xaxis.at, label = xaxis.lab, line = -0.45, 
         tick = FALSE, cex.axis = cex.axis)
    axis(side = 1, at = xaxis.at, label = rep("", length(xaxis.at)), 
         line = 0, tick = TRUE, cex.axis = cex.axis)
    title(xlab = paste("Overall Survival (years) and p.value: ",
                       format(p.value, 3)), line = 2, adj = .5, ...)
    title(ylab = paste("Survival Probability", "\n Hazard Ratio:", 
                       format(h.ratio$hazard.ratio, 3), " (", 
                       format(h.ratio$lower.ci, 3), ", ", 
                       format(h.ratio$upper.ci, 3), ")"), ...)

    if (!simple) {
        ## write group names
        group.name.pos <- (par()$usr[2] - par()$usr[1]) / -8
        padding <- abs(group.name.pos / 8)
        line.pos <- (seq(1, ng))[order(group.order)] + 3
        mtext(group.names, side = 1, line = line.pos, 
              at = group.name.pos, adj = 1, col = 1, las = 1)

        ## draw matching lines for n at risk.
        if (draw.lines) {
            par("xpd" = TRUE)
            for (i in seq(1, ng)) {
                axis(
                    side = 1,
                    at = c(group.name.pos + padding * 3, 0 - 3 * padding),
                    labels = FALSE,
                    line = line.pos[i] + 0.5,
                    lwd.ticks = 0,
                    col = col.surv[i],
                    lty = lty.surv[i],
                    lwd = 8,
                )
            }
        }

        ## numbers at risk
        kms <- summary(km, times = xaxis.at)
        if (is.null(kms$strata)) kms$strata <- rep(1, length(kms$time))
        d1 <- data.frame(time = kms$time, n.risk = kms$n.risk, 
                         strata = c(kms$strata))
        d2 <- split(d1, f = d1$strata)

        ## Right-justifying the numbers
        ndigits <- lapply(d2, function(x) nchar(x[, 2]))
        # max.len <- max(sapply(ndigits, length))
        max.len <- max(vapply(ndigits, length, FUN.VALUE = integer(1)))
        L <- do.call("rbind", lapply(ndigits, function(z) {
            length(z) <- max.len
            z
        }))
        nd <- apply(L, 2, max, na.rm = TRUE)
        for (i in seq(ng)) {
            this <- d2[[i]]
            w.adj <- strwidth("0", cex = cex.axis, 
                              font = par("font")) / 2 * nd[seq(1, nrow(this))]
            mtext(side = 1, at = this$time + w.adj, 
                  text = this$n.risk, line = line.pos[i], 
                  cex = cex.axis, adj = 1, col = 1, las = 1)
        }
        if (label.n.at.risk) mtext(side = 1, text = "N at risk", 
                                   at = group.name.pos, line = 1.5, 
                                   adj = 1, col = 1, las = 1, cex = cex.axis)
    } ## End of if(!simple)

    # Legend
    rlp <- group.order
    if (legend) {
        bgc <- ifelse(par("bg") == "transparent", "white", par("bg"))
        legend(
            x = loc.legend,
            legend = group.names[rlp],
            # lty = lty.surv[rlp],
            # lwd = lwd.surv[rlp],
            pch = 15,
            col = col.surv[rlp],
            bty = "o",
            cex = cex.axis,
            bg = bgc,
            box.col = "transparent",
            inset = .01
        )
    }

    ## draw confidence intervals
    for (i in seq(1, ng)) {
        this <- dat.list[[i]]
        x <- this$time
        L <- this$lower
        U <- this$upper
        S <- this$survival
        naL <- which(is.na(L))
        L[naL] <- L[naL - 1]
        U[naL] <- U[naL - 1]
        lines(x, L, type = "s", 
              col = alpha(col.ci[i], .5), lty = lty.ci[i], lwd = lwd.ci[i])
        lines(x, U, type = "s", 
              col = alpha(col.ci[i], .5), lty = lty.ci[i], lwd = lwd.ci[i])
    }
    # draw curves
    lines(km, conf.int = FALSE, col = col.surv, lty = lty.surv, 
          lwd = lwd.surv, mark = mark, xmax = xlim[2], ymin = ylim[1])
    box(bty = par("bty"))
    # par(op)
    if (is.null(geneName)) {
        title(main = paste0("Kaplan-Meier plot for ", nbetas, " selected genes"))
    } else {
        title(main = paste0("Kaplan-Meier plot (", geneName, ")"))
    }
}


# KM plot for custom groups -----------------------------------------------
plotKmCustomGroups <- function(genExpr, time, status,
                               group.assignation.vector, group.labels,
                               group.colors, boxplot = FALSE) {
    if (length(genExpr) != length(time)) {
        stop("the expr values and time vector must have the same length")
    }
    if (length(genExpr) != length(status)) {
        stop("the expr values and status vector must have the same length")
    }
    # gene name in vectors time status
    if (!identical(names(genExpr), names(time))) {
        stop("expr values names must match time vector names")
    }
    if (!identical(names(genExpr), names(status))) {
        stop("expr values names must match status vector names")
    }

    mSurv <- cbind(time, status)
    colnames(mSurv) <- c("time", "status")
    mSurv <- as.data.frame(mSurv)
    rownames(mSurv) <- names(time)
    # limit to 10 years and new censored
    mSurv$status[mSurv$time > 10] <- 0
    mSurv$time[mSurv$time > 10] <- 10.1

    n.genes <- 1
    n.samples <- length(genExpr)
    probesets.names <- names(genExpr)

    fits1 <- survfit(Surv(time, status) ~ group.assignation.vector, 
                     data = mSurv)
    log.rank.groups.surv <- 
      survdiff(Surv(time, status) ~ group.assignation.vector, data = mSurv)
    p.val <- (1 - pchisq(log.rank.groups.surv$chisq, 
                         length(log.rank.groups.surv$n) - 1))

    indices <- c(
        table(group.assignation.vector)[[1]],
        table(group.assignation.vector)[[2]]
    )
    group.colors <- rep(group.colors, indices)


    hazardR <- hazard.ratio(
        x = group.assignation.vector,
        surv.time = mSurv$time[match(colnames(genExpr), names(time))],
        surv.event = mSurv$status[match(colnames(genExpr), names(status))]
    )

    names(hazardR) <- c(
        "hazard.ratio", "coef", "se", "lower.ci", "upper.ci", "p.value",
        "n", "coxm", "data"
    )

    if (hazardR$hazard.ratio < 1) {
        stop("1/hazar.ratio was calculated")
        hazardR$hazard.ratio <- format(1 / as.numeric(hazardR$hazard.ratio, 3))
        hazardR$upper <- format(1 / as.numeric(hazardR$lower, 3))
        hazardR$lower <- format(1 / as.numeric(hazardR$upper, 3))
    }

    plot_values <- NULL
    plot_values$km <- list(
        "fitsKM" = fits1,
        "p.val" = p.val,
        "hazardR" = hazardR
    )

    rList <- list("plot_values" = plot_values)

    plotKM(rList,
        xaxis.at = c(0:10),
        col.surv = unique(group.colors),
        group.names = c(group.labels[1], group.labels[2]),
        col.ci = unique(group.colors)
    )
    if (boxplot) {
        # expression and t test
        vector.exprs.j <- as.numeric(genExpr)
        level1 <- vector.exprs.j[group.assignation.vector == 
                                   levels(factor(group.assignation.vector))[1]]
        level2 <- vector.exprs.j[group.assignation.vector == 
                                   levels(factor(group.assignation.vector))[2]]
        Ttest <- wilcox.test(level1, level2)


        mm <- data.frame("values" = vector.exprs.j, "group" = c(
            rep(
                unique(group.assignation.vector)[1],
                indices[1]
            ),
            rep(unique(group.assignation.vector)[2], indices[2])
        ))

        plot <- ggplot(data = mm, aes(x = group, y = values)) +
            stat_boxplot(geom = "errorbar") +
            geom_boxplot(aes(fill = group), coef = 1.5, alpha = .5, 
                         show.legend = FALSE, outlier.alpha = 0, notch = TRUE) +
            geom_jitter(aes(fill = group), color = group.colors, size = 2) +
            scale_fill_manual(labels = unique(mm$group), 
                              values = unique(group.colors)) +
            ggtitle(paste0("Expresion for gene ", rownames(genExpr), 
                           " in ", unique(mm$group)[1], " and ", 
                           unique(mm$group)[2])) +
            ylab("Expression level") +
            xlab(paste0("wilcox-test p.value for groups: ", "\n", 
                        Ttest$p.value)) +
            theme_bw() +
            theme(
                panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                legend.title = element_blank(),
                legend.text = element_text(size = 10),
                legend.key.size = unit(1.5, "cm")
            ) +
            guides(fill = guide_legend(override.aes = list(shape = 22, 
                                                           size = 10)))
        plot
    }
    msg <- paste("P value: ", p.val)
    message(msg)
}
