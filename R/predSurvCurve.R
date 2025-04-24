#' Predict Survival Curve for a Test Patient Using Cox Model Predictions
#'
#' This function computes the survival curve for a test patient based on their 
#' predicted risk score using a Cox proportional hazards model. It estimates 
#' the baseline survival function from training data and predicts the survival 
#' probability for the test patient at specified times.
#'
#' @param cox_pred_training A numeric vector representing the predicted risk 
#' scores (or log-risk scores) from the Cox model for the training set.
#' @param mSurv A data frame with two columns: "time" representing survival 
#' times, and "status" representing the event status (1 for event, 0 for 
#' censored).
#' @param cox_pred_test A numeric vector of length 1 representing the predicted 
#' risk score (or log-risk score) for the test patient.
#' @param eval_surv_times A numeric vector of times at which the survival curve 
#' should be evaluated. If NULL (default), the times will be taken from the 
#' training data up to the maximum event time.
#'
#' @details
#' This function calculates the baseline survival curve using the training data 
#' by first estimating the hazard function and cumulative hazard function.
#' Then, it computes the survival probability for the test patient using the 
#' baseline survival function and the test patient's risk score.
#' If the test patient's risk is associated with a hazard ratio greater than 1, 
#' the survival curve will decrease more rapidly.
#' If `eval_surv_times` is provided, the curve is evaluated at those specific 
#' times; otherwise, the function uses the survival times from the training 
#' data.
#'
#' @return A list containing:
#'   \describe{
#'     \item{\code{eval_times}}{The times at which the survival curve is 
#'     evaluated.}
#'     \item{\code{S_0_t}}{The baseline survival function evaluated at the 
#'     times in \code{eval_surv_times}.}
#'     \item{\code{S_test_patient}}{The survival curve for the test patient at 
#'     the times in \code{eval_surv_times}.}
#'   }
#'
#' @examples
#' # Example usage
#' cox_pred_training <- rnorm(100)
#' mSurv <- data.frame(time = rexp(100, rate = 0.1), 
#'                     status = sample(0:1, 100, replace = TRUE))
#' cox_pred_test <- 0.5
#' result <- predSurvCurve(cox_pred_training, mSurv, cox_pred_test)
#' plot(result$eval_times, result$S_test_patient, type = "l", 
#'      xlab = "Time", ylab = "Survival Probability")
#'
#' @export
predSurvCurve <- function(cox_pred_training, mSurv, 
                          cox_pred_test, eval_surv_times = NULL) {
    # Error control: Check if mSurv has the correct structure
    if (dim(mSurv)[[2]] != 2 || anyNA(match(c("time", "status"), 
                                            colnames(mSurv)))) {
        stop("mSurv should be a data frame with two columns: \"time\" ",
             "(Survival times) and \"status\" (Event status)")
    }

    # Error control: Check if cox_pred_training is a numeric vector
    if (!is.vector(cox_pred_training) || !is.numeric(cox_pred_training)) {
        stop("cox_pred_training should be a numeric vector")
    }

    # Error control: Check if the number of survival times matches the length 
    # of cox_pred_training
    if (dim(mSurv)[[1]] != length(cox_pred_training)) {
        stop("The number of survival times should be equal to the length ",
             "of cox_pred_training")
    }

    # Error control: Check if cox_pred_test is a numeric vector of length 1
    if (!is.vector(cox_pred_test) || length(cox_pred_test) != 1 || 
        !is.numeric(cox_pred_test)) {
        stop("cox_pred_test should be a numeric vector of length 1 ",
             "representing the risk of the test patient")
    }

    surv_times <- mSurv$time
    surv_status <- mSurv$status

    # Calculate hazard rate function
    event_times <- sort(unique(surv_times[surv_status == 1]))
    hazard_0 <- rep(0, length(event_times))

    for (i in seq_along(event_times)) {
        hazard_0[i] <- sum(surv_times[mSurv$status == 1] == event_times[i]) /
            sum(exp(cox_pred_training[surv_times >= event_times[i]]))
    }

    # Computation of cumulative hazard rate function
    Hazard_0 <- cumsum(hazard_0)

    # Baseline survival function
    S_0 <- exp(-Hazard_0)
    if (is.null(eval_surv_times)) {
        eval_surv_times <- surv_times[surv_times <= max(event_times)]
    } else if (max(eval_surv_times) > max(event_times)) {
        message("The Survival curve is only estimated for the ",
                "interval [%.1f, %.1f]", 0, max(event_times))
        message("No observation beyond this point in the training set")
        eval_surv_times <- eval_surv_times[eval_surv_times <= max(event_times)]
    }

    # Smoothed baseline survival function using Friedman super smoother
    S_0_smooth <- supsmu(event_times, S_0)

    # Survival curve estimation for test patient
    S_0_t <- approx(S_0_smooth$x, S_0_smooth$y, xout = eval_surv_times)$y
    S_z_t <- (S_0_t)^exp(cox_pred_test[1])

    # Reordering of survival times to plot the curve
    S_z_t <- S_z_t[order(eval_surv_times)]
    eval_surv_times <- eval_surv_times[order(eval_surv_times)]

    # Plot the survival curve in eval_surv_times for the test patient
    index_median_survival_t <- order(abs(S_z_t - 0.5))[1]
    plot(eval_surv_times, S_z_t, type = "l", ylim = c(0, 1), col = 2, 
         xaxt = "n", yaxt = "n", xpd = FALSE, xlab = "Time", 
         ylab = "Survival probability", 
         main = "Survival curve for the test patient")
    
    ticks_x <- seq(0, floor(max(eval_surv_times)))
    axis(1, at = seq(0, floor(max(eval_surv_times)), by = 1), 
         labels = ticks_x, tck = 1, lty = 2, col = "grey")
    axis(1, at = seq(0, floor(max(eval_surv_times))))
    axis(2, at = seq(0, 1, by = 0.2), labels = c(0.0, "", 0.4, "", 0.8, 1), 
         tck = 1, lty = 2, col = "grey")
    axis(2, at = seq(0, 1, by = 0.2), labels = rep("", 6))
    axis(2, at = S_z_t[index_median_survival_t], 
         labels = round(S_z_t[index_median_survival_t], 1), cex.axis = 0.8)
    segments(-0.4, S_z_t[index_median_survival_t], 
             eval_surv_times[index_median_survival_t], 
             S_z_t[index_median_survival_t], col = "blue", lty = 2)
    
    segments(eval_surv_times[index_median_survival_t], -0.05, 
             eval_surv_times[index_median_survival_t], 
             S_z_t[index_median_survival_t], col = "blue", lty = 2)

    return(list(eval_times = eval_surv_times, S_0_t = S_0_t, 
                S_test_patient = S_z_t))
}
