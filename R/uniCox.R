#' uniCox: Fit a Penalized Cox Proportional Hazards Model
#'
#' This function fits a Cox proportional hazards model with L1 regularization 
#' (Lasso) for survival analysis. The regularization parameter, lambda, is 
#' chosen through an iterative procedure using coordinate descent.
#'
#' @param x A numeric matrix of size n by p (n samples, p features). 
#' The input gene expression matrix.
#' @param y A numeric vector of length n. Survival times corresponding to the 
#' n samples.
#' @param status A numeric vector of length n. Binary vector where 1 indicates 
#' the event (death) occurred, and 0 indicates censored data.
#' @param lamlist A numeric vector of lambda values for regularization. 
#' If NULL, a sequence of values is generated.
#' @param nlam The number of lambda values to consider. Default is 20.
#' @param del.thres Convergence threshold for the coefficient update step. 
#' Default is 0.01.
#' @param max.iter Maximum number of iterations for each lambda value. 
#' Default is 5.
#'
#' @return A list of class "uniCoxFit" containing the following components:
#' \item{lamlist}{The list of lambda values used for regularization.}
#' \item{beta}{A matrix of size p by nlam, where each column contains the 
#' coefficients for a different lambda value.}
#' \item{nfeatures}{A vector of length nlam, where each element indicates 
#' the number of non-zero features selected for each lambda value.}
#' \item{mx}{The mean of each feature in the input matrix x.}
#' \item{vx}{The variance of each feature in the input matrix x.}
#' \item{s0}{The scaling parameter, which is the median of the square root of 
#' the variances.}
#' \item{call}{The original function call.}
#'
#' @examples
#' # Example usage
#' # x: Gene expression matrix (samples x genes)
#' # y: Survival times
#' # status: Censoring status (0 = censored, 1 = event)
#' fit <- uniCox(x, y, status)
#' print(fit$beta) # Coefficients for each lambda value
#' print(fit$nfeatures) # Number of selected features for each lambda
#'
#' @importFrom stats scale
#' @export
uniCox <- function(x, y, status, lamlist = NULL, nlam = 20, 
                   del.thres = .01, max.iter = 5) {
    # x is n by p
    this.call <- match.call()
    # Error control: check if input dimensions match
    if (nrow(x) != length(y)) {
        stop("The number of rows in the gene expression matrix 'x' must ",
             "match the length of the survival times 'y'.")
    }
    if (nrow(x) != length(status)) {
        stop("The number of rows in the gene expression matrix 'x' must ",
             "match the length of the status vector.")
    }

    # Error control: check if status vector contains only 0 (censored) 
    # and 1 (event occurred)
    if ((sum(status == 0) + sum(status == 1)) != length(status)) {
        stop("The 'status' vector must contain only 0 (censored) and ",
             "1 (event occurred) values.")
    }

    mx <- colMeans(x)
    x <- scale(x, center = mx, scale = FALSE)
    vx <- 1 / coxvar(t(x), y, status)
    s0 <- quantile(sqrt(vx), .5)
    xs <- scale(x, center = FALSE, scale = sqrt(vx) + s0)

    u0 <- coxscor(t(xs), y, status)$scor
    inffull <- NULL
    if (is.null(lamlist)) {
        lamlist <- seq(0, max(abs(u0)), length = nlam)
    }
    nlam <- length(lamlist)

    p <- ncol(xs)
    beta2 <- matrix(0, nrow = p, ncol = length(lamlist))
    beta0 <- rep(0, p)
    for (k in seq(1, length(lamlist))) {
        message(c("lambda value ", k, "out of ", length(lamlist)), fill = TRUE)
        lam <- lamlist[k]
        om <- abs(u0) > lam
        beta <- rep(0, sum(om))
        if (k > 1) {
            beta <- beta2[om, k - 1, drop = FALSE]
        }
        niter <- 0
        u <- rep(1, p)
        go <- TRUE
        while (niter < max.iter & go & sum(om)) {
            niter <- niter + 1
            offset <- t(scale(xs[, om, drop = FALSE], center = FALSE, 
                              scale = 1 / beta))
            u <- coxscor2(t(xs[, om, drop = FALSE]), y, status, 
                          offset = offset)$scor - lam * sign(beta)
            v <- coxvar2(t(xs[, om, drop = FALSE]), y, status, offset = offset)
            obeta <- beta
            beta <- obeta + (u / v)
            del <- mean(abs(beta - obeta), na.rm = TRUE)
            if (del < del.thres) {
                go <- FALSE
            }
        }
        if (lam == 0) {
            betafull <- beta
            inffull <- v
        }
        beta2[om, k] <- beta
    }

    nfeatures <- colSums(beta2 != 0, na.rm = TRUE)
    nm <- sum(is.na(beta2))
    if (nm > 0) {
        message(c(nm, " betas missing"), fill = TRUE)
    }
    junk <- list(lamlist = lamlist, beta = beta2, 
                 nfeatures = nfeatures, mx = mx, vx = vx, s0 = s0, 
                 call = this.call)
    class(junk) <- "uniCoxFit"
    return(junk)
}
