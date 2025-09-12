#' Predict Patient Risk Based on Gene Expression Data
#'
#' Function to predict the risk for new patients considering the gene 
#' expression for a subset of genes and the multivariate Cox regression model 
#' trained by function patientRisk, it may be used to predict over a single 
#' patient.
#'
#' @param model.fit A list containing the pre-fitted model and necessary 
#' parameters for risk prediction, including the optimal lambda value, risk 
#' thresholds, and plot values.
#' @param mExpr.testData A data frame representing the gene expression data of 
#' test patients, where each row is a gene and each column is a sample.
#'
#' @details
#' A risk score is estimated for new patients considering the optimal 
#' regularized multivariate Cox regression model trained by the function 
#' patientRisk. The risk score is normalized to be interpretable in the 
#' scale (0-100). The function generates a risk plot for new patients and 
#' stratifies them in three risk groups (low, intermediate, high) considering 
#' the thresholds learned by function patientRisk.
#'
#' @return A list containing the following elements:
#' \itemize{
#'  \item{\code{risk_score}: A vector with the unscaled risk score for new 
#'  patients estimated by a multivariate Cox regression model.}
#'  \item{\code{scaled_risk_score}:  A vector with the risk score for the new 
#'  patients scaled to be interpretable in the range 0-100.}
#'  \item{\code{plot_values}: A list containing information for visualizing the 
#'  sigmoid curve and risk thresholds.}
#'  }
#'  
#' @examples
#' data(seBRCA)
#' 
#' # prefilterSAM ---
#' data(ex_prefilterSAM)
#' 
#' # genePheno ---
#' data(ex_genePheno)
#' 
#' # patientRisk ---
#' data(ex_patientRisk)
#'                                            
#' # Simulate expression data
#' num_samples <- 20
#' geneList <- names(ex_genePheno$genes)
#' set.seed(5)
#' mExprs_testData <- matrix(rnorm(length(geneList) * num_samples, 
#'                           mean = 10, sd = 3),
#'                           nrow = length(geneList), ncol = num_samples)
#' 
#' # Assign row names (genes) and column names (samples)
#' rownames(mExprs_testData) <- geneList
#' colnames(mExprs_testData) <- paste0("Sample", 1:num_samples)
#' 
#' set.seed(5)
#' risk_prediction_validation_set <- predict_PatientRisk(
#'                                                 ex_patientRisk, 
#'                                                 mExprs_testData)
#' 
#' 
#' # Example for single patient prediction: Patient fourth is selected.
#' ex_predictPatientRisk <- data.frame(mExprs_testData[, 4])
#' colnames(ex_predictPatientRisk) <- colnames(mExprs_testData)[4]
#' # Risk prediction for the optimal subset of genes selected by patientRisk() 
#' set.seed(5)
#' risk_prediction_one_patient <- predict_PatientRisk(
#'                                                 ex_patientRisk, 
#'                                                 ex_predictPatientRisk)
#' 
#' @references
#' \itemize{
#'   \item{\insertRef{martinezromero2018}{asuri}} 
#'   \item{\insertRef{BuenoFortes2023}{asuri}}
#' } 
#' 
#' @export
predict_PatientRisk <- function(model.fit, mExpr.testData) {
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

    riskscore_testData <- predict_uniCox(model.fit$model.optimalLambda, 
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
      if (scaled_riskscore_testData >= 0 && 
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
