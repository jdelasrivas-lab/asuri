#' Double Nested Cross-Validation for Lambda Optimization in Survival Analysis
#'
#' This function performs double nested cross-validation to optimize the lambda
#' parameter based on the risk matrix.
#' The goal is to find the lambda value that results in the lowest p-value in
#' the Kaplan-Meier survival analysis.
#'
#' @param matrixOfRisks A numeric matrix where each column corresponds to the
#' risk score of each sample for a particular lambda value.
#' @param mSurv A data.frame with survival data with the following columns:
#'  \itemize{
#'  \item \code{time}: Survival time for each sample.
#'  \item \code{status}: Event indicator (1 = event occurred, 0 = censored).
#'   }
#' @param thresholds A numeric vector of two elements that defines the index
#' range (lowIndex, highIndex) for grouping in each lambda 
#' optimization iteration. If \code{NULL}, the function uses the 
#' 30%-70% range (default = NULL).
#'
#' @details
#' This function iterates over different lambda values, and for each lambda,
#' it uses nested cross-validation to find the optimal cutoff for grouping the
#' samples based on their risk scores. It performs Kaplan-Meier survival
#' analysis for each group and calculates p-values using the log-rank test.
#' The lambda value that results in the lowest p-value is considered optimal.
#'
#' @return
#' A list containing the following elements:
#' \describe{
#'   \item{\code{p.vals}}{A matrix where each column corresponds to the
#'   p-values obtained for each lambda value across different cutoffs.}
#' }
#'
#' @examples
#' # Simulated data example
#' set.seed(123)
#' matrixOfRisks <- matrix(rnorm(1000), ncol = 10)
#' mSurv <- data.frame(
#'   time = rexp(100, rate = 0.1),
#'   status = sample(0:1, 100, replace = TRUE)
#' )
#' result <- function_double_nested_lambda(matrixOfRisks, mSurv)
#' print(result$p.vals)
#'
#' @export
function_double_nested_lambda <- function(matrixOfRisks, 
                                          mSurv, thresholds = NULL) {
    nLambdas <- dim(matrixOfRisks)[2]
    nSamples <- dim(matrixOfRisks)[1]
    p.vals <- NULL
    groupsByLambda <- NULL

    # Optimized lambda using double nested CV
    if (is.null(thresholds)) {
        message("Nested Cross Validation: optimizing lambda...\n")
    }

    # for (i in 1:nLambdas) {
    for (i in seq(1, nLambdas)) {
        progress <- paste0("Progress: ", round(i * 100 / nLambdas), "%")
        message("\r", progress)
        flush.console()
        # in each iteration (by lambdas) the risk score is used to optimize
        # separability of KM curves and thus obtain the optimal logrank p.val
        order.riskDefinitive <- matrixOfRisks[, i][order(matrixOfRisks[, i])]
        group.assignation.vector <- rep(0, nSamples)
        p.val <- rep(1, nSamples)
        if (is.null(thresholds)) {
            lowIndex0.3 <- round(0.3 * nSamples)
            highIndex0.7 <- round(0.7 * nSamples)
        } else {
            lowIndex0.3 <- thresholds[1]
            highIndex0.7 <- thresholds[2]
        }
        for (j in lowIndex0.3:highIndex0.7) {
            group.assignation.vector[matrixOfRisks[, i] < 
                                       order.riskDefinitive[j]] <- 1
            group.assignation.vector[matrixOfRisks[, i] >= 
                                       order.riskDefinitive[j]] <- 2
            msurv_data <- mSurv[match(rownames(matrixOfRisks), 
                                      rownames(mSurv)), ]
            log.rank.groups.surv <- survdiff(Surv(time, status) ~ 
                                               group.assignation.vector, 
                                             data = msurv_data)
            
            p.val[j] <- pchisq(log.rank.groups.surv$chisq, 
                               length(log.rank.groups.surv$n) - 1, 
                               lower.tail = FALSE)
        }
        # min p.val is restrained to 40% central interval, in order to observe 
        # relative minima (if they are near the interval limits)
        p.vals <- cbind(p.vals, p.val)
        lowest.p.value.index <- which.min(p.val)

        # cutPoint is used to create the groups
    }

    if (is.null(thresholds)) {
        message("\nRisk predicted!\n")
    }
    # the function returns a vector with optimal p.values for each lambda 
    # value in order to choose a lambda
  rList <- list("p.vals" = p.vals)
  return(rList)
}
