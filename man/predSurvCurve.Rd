\name{predSurvCurve}
\alias{predSurvCurve}
\title{ 
  predSurvCurve
}
\description{
  Function to predict the expected survival curve for a test patient given the risk score provided by a Cox regression model.
}
\usage{
  predSurvCurve( cox_pred_training, mSurv, cox_pred_test, eval_surv_times = NULL )
}
\arguments{
  \item{cox_pred_training}{ 
    A numeric vector representing the predicted risk scores (or log-risk scores) from the Cox model for the training set. It can be obtained from *predict.patientRisk* function.
  }
  \item{mSurv}{
     A data frame with two columns: "time" representing survival times, and "status" representing the event status (1 for event, 0 for censored) for the training dataset.
  }
    \item{cox_pred_test}{
    A numeric vector of length 1 representing the predicted risk score (or log-risk score) for the test patient. It can be obtained from *predict.patientRisk* function.
  }
   \item{eval_surv_times}{
    A numeric vector of times at which the survival curve should be evaluated. If NULL (default), the times will be taken from the training data up to the maximum event time.
  }
}
\details{
 This function will obtain the expected survival curve for a test patient given the risk score estimated by a multivariate Cox regression model. The Breslow estimator is considered to approximate the baseline Hazard. The cumulative baseline survival function is smoothed using the Friedman's super smoother. The function interpolates linearly the survival probability for time points not in the training set. }
\value{
The function provides a plot for the expected survival curve. It shows also the median survival time, that is, the time at which the survival probability drops to 0.5.
  The function also returns a list with the following objects:
  \item{eval_times}{
Vector of times in years at which the survival curve has been evaluated.
  }
  \item{S_0_t}{
Vector with the baseline survival probability estimated at each evaluation time.
  }
  \item{S_test_patient}{
Vector with the survival probability estimated at each evaluation time for a given test patient.
  }
}
\examples{
  data(predSurvCurve)
  
  # COX prediction for the training set
  set.seed(5)
  cox_pred_training <- predict.patientRisk(multivariate_risk_predictor, mExprSelectedGenes)
  cox_pred_training$risk_score
  
  # COX prediction for the test patient
  set.seed(5)
  cox_pred_test <- predict.patientRisk(multivariate_risk_predictor, mExprs_testSingleData)
  cox_pred_test$risk_score
  
  # Survival curve estimation
  eval_surv_times <- seq(0, max(mPheno$time), by = 0.1)
  set.seed(5)
  surv_curv_cox <- predSurvCurve(cox_pred_training$risk_score, mPheno[, c(2, 3)], cox_pred_test$risk_score, eval_surv_times)
}
\references{
  Martinez-Romero, J., Bueno-Fortes, S., Martín-Merino, M. \emph{et al.} (2018) \bold{Survival marker genes of colorectal cancer derived from consistent transcriptomic profiling}. \emph{BMC Genomics} 19 (Suppl 8), 857 . \doi{10.1186/s12864-018-5193-9}

  Santiago Bueno-Fortes, Alberto Berral-Gonzalez, José Manuel Sánchez-Santos, Manuel Martin-Merino, Javier De Las Rivas (2023) \bold{Identification of a gene expression signature associated with breast cancer survival and risk that improves clinical genomic platforms}, \emph{Bioinformatics Advances}, Volume 3, Issue 1, vbad037, \doi{10.1093/bioadv/vbad037}
}
