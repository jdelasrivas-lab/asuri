\name{predict.patientRisk}
\alias{predict.patientRisk}
\title{ 
  predict.patientRisk
}
\description{
  Function to predict the risk for new patients considering the gene expression for a subset of genes and the multivariate Cox regression model trained by function patientRisk, it may be used to predict over a single patient.
}
\usage{
  predict.patientRisk(model.fit, mExpr.testData)
}
\arguments{
  \item{model.fit}{ 
    A list containing the pre-fitted model and necessary parameters for risk prediction, including the optimal lambda value, risk thresholds, and plot values given by *patientRisk()* function.
  }
  \item{mExpr.testData}{
    A data frame for an independent set of patients representing the gene expression data, where each row is a gene and each column is a sample. The genes must be the same subset of genes (as rows) considered in the training of the Cox model using patientRisk function. *colnames(mExpr)* corresponds to the sample names and *rownames(mExpr)* to the gene names.
  }
}
\details{
  A risk score is estimated for new patients considering the optimal regularized multivariate Cox regression model trained by the function patientRisk. The risk score is normalized to be interpretable in the scale (0-100). The function generates a risk plot for new patients and stratifies them in three risk groups (low, intermediate, high) considering the thresholds learned by function patientRisk.
}
\value{
  This function generates one plot with the risk for new patients scaled to be interpretable in the range 0-100. Low risk pattients are shown in blue while high risk patients in red. If only the expression profile of a single patient is provided, the function predicts the risk and classifies the sample into one of the three possible risk groups.
  The function also returns a list with the following objects:
  \item{risk_score}{
	 A vector with the unscaled risk score for new patients estimated by a multivariate Cox regression model.
  }
  \item{scaled_risk_score}{
    A vector with the risk score for the new patients scaled to be interpretable in the range 0-100.
  }
}
\examples{
data(matrixData)
set.seed(5)

# patientRisk function estimates the parameters of a multivariate Cox model
multivariate_risk_predictor <- patientRisk(mExprSelectedGenes, time, status)

# Risk prediction for a subset of new patients. If you don't have the multivariate_risk_predictor object, run the previous line
# mExprs_testData is the gene expression matrix for the test samples (as columns) and the same subset 
# of genes (as rows) considered in the training of the Cox model using patientRisk function.
risk_prediction_validation_set <- predict.patientRisk(multivariate_risk_predictor, mExprs_testData)

# Risk prediction for a single patient. Patient fourth is selected.
mExprs_testSingleData <- data.frame( mExprs_testData[, 4 ] )
colnames( mExprs_testSingleData ) <- colnames(mExprs_testData)[ 4 ]
risk_prediction_one_patient <- predict.patientRisk(multivariate_risk_predictor, mExprs_testSingleData)
}
\references{
  Martinez-Romero, J., Bueno-Fortes, S., Martín-Merino, M. \emph{et al.} (2018) \bold{Survival marker genes of colorectal cancer derived from consistent transcriptomic profiling}. \emph{BMC Genomics} 19 (Suppl 8), 857 . \doi{10.1186/s12864-018-5193-9}

  Santiago Bueno-Fortes, Alberto Berral-Gonzalez, José Manuel Sánchez-Santos, Manuel Martin-Merino, Javier De Las Rivas (2023) \bold{Identification of a gene expression signature associated with breast cancer survival and risk that improves clinical genomic platforms}, \emph{Bioinformatics Advances}, Volume 3, Issue 1, vbad037, \doi{10.1093/bioadv/vbad037}
}
