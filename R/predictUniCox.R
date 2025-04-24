#' Predict Risk Scores Using a UniCox Model
#'
#' This function predicts risk scores based on a pre-trained uniCox model using 
#' the provided input data. It scales the input data based on the model's 
#' parameters (mean and variance) and computes the risk score by multiplying 
#' the scaled input with the model's coefficients.
#'
#' @param object A list representing a pre-trained uniCox model. It should 
#' contain the following elements:
#'   \itemize{
#'     \item \code{mx}: A numeric vector of means for each feature used in 
#'     scaling the input data.
#'     \item \code{vx}: A numeric vector of variances for each feature used 
#'     in scaling the input data.
#'     \item \code{s0}: A numeric value representing the scaling factor for 
#'     the input data.
#'     \item \code{beta}: A numeric vector of coefficients (model parameters)
#'      for the model.
#'   }
#' @param x A numeric matrix or data frame of input data to be used for 
#' prediction. The matrix should have dimensions \code{n x p}, where \code{n} 
#' is the number of samples and \code{p} is the number of features (variables).
#' @param ... Additional arguments passed to or from other methods (not used 
#' in this function).
#'
#' @details
#' The function first scales the input data (\code{x}) using the provided model 
#' parameters (\code{mx}, \code{vx}, \code{s0}). Then, it calculates the 
#' predicted risk score by multiplying the scaled input data with the model 
#' coefficients (\code{beta}). The output is a numeric vector of predicted 
#' risk scores.
#'
#' @return A numeric vector of predicted risk scores for each sample in the 
#' input data.
#'
#' @examples
#' # Example usage with a pre-trained model and input data
#' model <- list(mx = rep(0, 10), vx = rep(1, 10), s0 = 1, beta = rnorm(10))
#' input_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' predicted_risk <- predict.uniCox(model, input_data)
#' print(predicted_risk)
#'
#' @export
predict.uniCox <- function(object, x, ...) {
    # Error control: Check if object contains required elements
    if (is.null(object$mx) || is.null(object$vx) || 
        is.null(object$s0) || is.null(object$beta)) {
        stop("The input model 'object' must contain 'mx', 'vx', 's0', ",
             "and 'beta'.")
    }

    # Error control: Check if input data 'x' is a matrix or data frame
    if (!is.matrix(x) && !is.data.frame(x)) {
        stop("The input 'x' should be a matrix or data frame.")
    }
    # x is n by p
    mx <- object$mx
    vx <- object$vx
    s0 <- object$s0
    beta <- object$beta
    x <- scale(x, center = mx, scale = FALSE)
    xs <- scale(x, center = FALSE, scale = sqrt(vx) + s0)
    yhat <- x %*% beta
    return(yhat)
}
