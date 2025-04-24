#' Print Method for uniCoxFit Objects
#'
#' This function prints a summary of the `uniCoxFit` object, which includes the 
#' call to the function that created the object, the list of lambda values 
#' (regularization parameters) used during the fitting process, and the number 
#' of features (genes) selected.
#'
#' @param x An object of class `uniCoxFit` (from a `uniCox` fitting procedure).
#' @param ... Additional arguments (not used in this method).
#'
#' @details
#' The printed output includes the following:
#'  - The original call that generated the `uniCoxFit` object.
#'  - A table with the lambda values (regularization parameters) of the model.
#'  - The number of features (genes) selected for each lambda value.
#'
#' @return This function prints a summary to the console.
#'
#' @examples
#' # Example usage after running the uniCox model
#' fit <- uniCox(X, y) # X is the predictor matrix, y is the survival outcome
#' print(fit)
#'
#' @export
print.uniCoxFit <- function(x, ...) {
    message("Call:\n")
    dput(x$call)
    mat <- rbind(lambda = format(round(x$lamlist, 3)), 
                 number.of.genes = x$nfeatures)
    dimnames(mat) <- list(dimnames(mat)[[1]], paste(seq(1, ncol(mat))))
    print(t(mat), quote = FALSE)
    invisible()
}
