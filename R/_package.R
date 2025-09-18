#' asuri: \strong{A}nalysis of disease \strong{SU}rvival and patient 
#' \strong{RI}sk prediction based on gene signatures
#'
#' The \pkg{asuri} package provides functions to discovers marker genes that 
#' are related to risk prediction capabilities and to a clinical variable of 
#' interest. It uses two main steps, including subsampling glmnet and unicox. 
#' The package implements robust functions to discover survival markers 
#' related to a clinical phenotype and to predict a risk score, allowing 
#' to study the patient's risk based on the gene signatures. Several plots 
#' are provided to visualise the relevance of the genes, the risk score, 
#' and patient stratification, as well as a robust version of the 
#' Kaplan-Meier curves.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{prefilterSAM}}: Bootstrapped differential expression based on SAM.
#'   \item \code{\link{genePheno}}: Analyze gene-phenotype associations.
#'   \item \code{\link{geneSurv}}: Compute survival analysis for genes.
#'   \item \code{\link{patientRisk}}: Calculate patient risk scores.
#'   \item \code{\link{predict_PatientRisk}}: Predict patient risk for new data.
#'   \item \code{\link{predict_SurvCurve}}: Predict survival curves for patients.
#'   \item Plotting functions: 
#'     \code{\link{plotBoxplot}}, \code{\link{plotProbClass}}, \code{\link{plotLogRank}}, 
#'     \code{\link{plotSigmoid}}, \code{\link{plotLambda}}, \code{\link{plotBetas}}, \code{\link{plotKM}}.
#' }
#' @return This page does not return any value. Provides general usage information.
#' @name asuri
NULL