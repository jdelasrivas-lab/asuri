#' Predict Patient Risk Based on Gene Expression Data
#'
#' This function calculates the risk score for test patients based on the gene 
#' expression data using a pre-fitted model. It scales the risk scores and 
#' classifies patients into risk groups (low, intermediate, or high). 
#' Additionally, it provides a visualization of the risk scores and thresholds 
#' using a sigmoid curve.
#'
#' @param model.fit A list containing the pre-fitted model and necessary 
#' parameters for risk prediction, including the optimal lambda value, risk 
#' thresholds, and plot values.
#' @param mExpr.testData A data frame representing the gene expression data of 
#' test patients, where each row is a gene and each column is a sample.
#'
#' @details
#' This function uses a pre-trained Cox model to predict the risk scores for 
#' test patients. The risk scores are then scaled and categorized into three 
#' risk groups: low, intermediate, and high. The thresholds for the risk groups 
#' are based on the pre-fitted model. The function also allows visualization of 
#' the risk scores. If only one sample is provided, it will print the 
#' normalized risk score and classification (low, intermediate, or high).
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{risk_score}}{The raw risk score for each sample in the 
#'   test dataset.}
#'   \item{\code{scaled_risk_score}}{The normalized risk score 
#'   (scaled from 0 to 100) for each sample.}
#'   \item{\code{plot_values}}{A list containing information for visualizing 
#'   the sigmoid curve and risk thresholds.}
#' }
#'
#' @examples
#' # Example usage with a pre-fitted model and gene expression data
#' model_fit <- load_model("path_to_pretrained_model")
#' test_data <- load_test_data("path_to_test_data")
#' result <- predict.patientRisk(model.fit = model_fit, 
#'                               mExpr.testData = test_data)
#' print(result$scaled_risk_score)
#'
#' @export
predict.patientRisk <- function(model.fit, mExpr.testData) {
    mExpr.testData <- as.data.frame(mExpr.testData)
    # Error control: Check for row names in the gene expression matrix
    if (is.null(rownames(mExpr.testData))) {
        stop("Row names in gene expression matrix should be the gene names")
    }

    # Error control: Check for column names in the gene expression matrix
    if (is.null(colnames(mExpr.testData))) {
        stop("Column names in gene expression matrix should be the ",
             "sample names")
    }

    riskscore_testData <- asuri:::predict.uniCox(model.fit$model.optimalLambda, 
                                                 t(mExpr.testData))

    plot_values <- model.fit$plot_values

    scaled_riskscore_testData <- 
      ((riskscore_testData - model.fit$range.risk[1]) / abs(
        model.fit$range.risk[2] - model.fit$range.risk[1])) * 100
    
    ordered_scaled_riskscore_testData <- 
      scaled_riskscore_testData[order(scaled_riskscore_testData)]
    min_value <- min(ordered_scaled_riskscore_testData)
    ordered_scaled_riskscore_testData <- if (min_value < 0) {
        ordered_scaled_riskscore_testData + abs(min_value)
    } else {
        ordered_scaled_riskscore_testData
    }

    low_risk_threshold <- model.fit$riskThresholds[1, 1]
    high_risk_threshold <- model.fit$riskThresholds[3, 1]
    cutpoint_risk_threshold <- model.fit$riskThresholds[2, 1]

    cutpoints <- c(low_risk_threshold, 
                   high_risk_threshold, 
                   cutpoint_risk_threshold)
    names(cutpoints) <- c("low.cutpoint", 
                          "high.cutpoint", 
                          "cutpoint")

    cutpoints_check <- 
      cutpoints[!cutpoints %in% ordered_scaled_riskscore_testData]
    
    ordered_scaled_riskscore_testData <- 
      c(ordered_scaled_riskscore_testData, cutpoints_check)
    ordered_scaled_riskscore_testData <- 
      ordered_scaled_riskscore_testData[
        order(ordered_scaled_riskscore_testData)]

    plot_values$sigmoid_logrank$orderednormalized_risk <- 
      ordered_scaled_riskscore_testData
    
    plot_values$sigmoid_logrank$lowIndexReal <- 
      which(ordered_scaled_riskscore_testData == low_risk_threshold)
    
    plot_values$sigmoid_logrank$highIndexReal <- 
      which(ordered_scaled_riskscore_testData == high_risk_threshold)
    
    plot_values$sigmoid_logrank$cutPoint <- 
      which(ordered_scaled_riskscore_testData == cutpoint_risk_threshold)

    if (dim(mExpr.testData)[2] == 1) {
        message("Normalized patient Risk (0 100): ", 
                scaled_riskscore_testData, "\n")
        if (scaled_riskscore_testData >= 0 & 
            scaled_riskscore_testData < low_risk_threshold) {
            message("The patient is classified as Low Risk \n")
            message("Low Risk interval: (", 0, ", ", low_risk_threshold, ")")
        } else if (scaled_riskscore_testData < high_risk_threshold) {
            message("The patient is categorized as Intermediate Risk", "\n")
            message("Medium Risk interval (", 
                    c(low_risk_threshold, 
                      high_risk_threshold), ")")
        } else if (scaled_riskscore_testData <= 100) {
            message("The patient is categorized as High Risk \n")
            message("High Risk interval (", 
                    c(high_risk_threshold, 100), ")", "\n")
        } else {
            stop("Risk score is outside the accepted range \n")
        }
    } else {
        # Plot risk score for test patients
        model.fit$plot_values <- plot_values
        plotSigmoid(model.fit)
    }

    riskscore_testData <- t(riskscore_testData)
    scaled_riskscore_testData <- t(scaled_riskscore_testData)
    samplenames <- colnames(riskscore_testData)

    riskscore_testData <- as.vector(riskscore_testData)
    scaled_riskscore_testData <- as.vector(scaled_riskscore_testData)

    names(riskscore_testData) <- samplenames
    names(scaled_riskscore_testData) <- samplenames

    predicted.riskscore <- list(
        risk_score = riskscore_testData,
        scaled_risk_score = scaled_riskscore_testData,
        plot_values = plot_values
    )
}
