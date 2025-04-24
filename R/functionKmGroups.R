#' Kaplan-Meier Group Stratification Based on Gene Expression
#'
#' This function stratifies samples into two groups (high/low expression) 
#' based on gene expression levels to identify optimal cut-points that minimize 
#' p-values in survival analysis.
#'
#' @param genExpr A numeric vector representing the expression levels of a 
#' gene across samples.
#' @param mSurv A data.frame containing survival data with the 
#' following columns:
#'   \itemize{
#'     \item \code{time}: Survival time for each sample.
#'     \item \code{status}: Event indicator (1 = event occurred, 0 = censored).
#'   }
#' @param geneName A character vector of gene names corresponding to the
#'  expressions in \code{genExpr}.
#'
#' @details
#' The function uses the log-rank test to assess the difference in survival 
#' between two groups formed by varying cut-points in the expression levels. 
#' It identifies the cut-point that minimizes the p-value for each gene.
#'
#' @return
#' A list with the following elements:
#' \describe{
#'   \item{\code{matrix.groups}}{A matrix where each row corresponds to the 
#'   group assignment for a gene. Columns represent samples.}
#'   \item{\code{best_pvalue}}{The smallest p-value achieved across all 
#'   genes and cut-points.}
#' }
#'
#' @examples
#' # Simulated data example
#' set.seed(123)
#' genExpr <- rnorm(100)
#' mSurv <- data.frame(
#'     time = rexp(100, rate = 0.1),
#'     status = sample(0:1, 100, replace = TRUE)
#' )
#' geneName <- c("Gene1")
#' result <- functionKmGroups(genExpr, mSurv, geneName)
#' print(result$best_pvalue)
#'
#' @export

functionKmGroups <- function(genExpr, mSurv, geneName) {
    if (!is.numeric(genExpr)) {
        stop("'genExpr' must be a numeric vector with the gene expressions.")
    }
    if (!inherits(mSurv, "data.frame")) {
        stop("'mSurv' must be a data.frame with columns “time” and “status”.")
    }
    if (!all(c("time", "status") %in% colnames(mSurv))) {
        stop("'mSurv' must contain the columns “time” ", 
             "(survival time) and “status” (event).")
    }
    if (!is.character(geneName)) {
        stop("'geneName' must be a character vector ",
             "with the names of the genes.")
    }
    if (length(genExpr) != nrow(mSurv)) {
        stop("The length of “genExpr” must match ", 
             "the number of rows in “mSurv”.")
    }

    matrix.groups <- NULL
    n.genes <- 1
    n.samples <- length(genExpr)

    for25 <- round(n.samples * 0.25)
    for75 <- round(n.samples * 0.75)
    p.value.genes.ordering <- rep(0, n.genes)


    for (j in geneName) {
        vector.exprs.j <- as.numeric(genExpr)
        order.vector.exprs.j <- order(vector.exprs.j)

        group.assignation.vector <- rep(0, n.samples)
        p.val <- rep(1, n.samples)
        for (i in seq(for25, for75)) {
            # run for each possible group computing p values to take the lowest
            # group1 <- order.vector.exprs.j[1:i]
            group1 <- order.vector.exprs.j[seq(1, i)]
            group2 <- order.vector.exprs.j[seq((i + 1), n.samples)]
            group.assignation.vector[group1] <- 1
            group.assignation.vector[group2] <- 2
            log.rank.groups.surv <- survdiff(Surv(time, status) ~ 
                                               group.assignation.vector, 
                                             data = mSurv)
            p.val[i] <- 1 - pchisq(log.rank.groups.surv$chisq, 
                                   length(log.rank.groups.surv$n) - 1)
        }

        ordered.pval.indexes <- order(p.val)


        # searching lowest p value and extracting index to order final cut point
        lowest.pvalue.index <- ordered.pval.indexes[1]
        group1 <- order.vector.exprs.j[seq(1, lowest.pvalue.index)]
        group2 <- order.vector.exprs.j[seq((lowest.pvalue.index + 1), 
                                           n.samples)]
        group.assignation.vector[group1] <- 1
        group.assignation.vector[group2] <- 2

        p.value.genes.ordering[j] <- p.val[lowest.pvalue.index]

        if (j == 1) {
            matrix.groups <- group.assignation.vector
        } else {
            matrix.groups <- rbind(matrix.groups, group.assignation.vector)
        }
        # performing the last fit using optimised groups
        fits1 <- survfit(Surv(time, status) ~ group.assignation.vector, 
                         data = mSurv)
    }
    list(as.vector(matrix.groups), p.val[lowest.pvalue.index])
}
