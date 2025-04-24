\name{patientRisk}
\alias{patientRisk}
\title{ 
  patientRisk
}
\description{
  This function selects a subset of good risk markers and estimates a multivariate risk score based on the UNICOX algorithm. The patients are stratified into two or more prognostic groups based on the risk score. The Cox regression is trained using a ten-fold double nested crossvalidation strategy to avoid overfitting.
}
\usage{
  patientRisk(mExpr, time, status, 
  method = c("min.pval", "med.pval", "class.probs"), nboot = 50)
}
\arguments{
  \item{mExpr}{ 
    A numeric matrix of gene expression data where each row represents a gene and each column represents a sample.  
  }
  \item{time}{
    A numeric vector representing the survival time for each sample, names(time) should be defined as the sample names.
  }
  \item{status}{
    A numeric vector representing the event status (1 = event, 0 = censored) for each sample, names(status) should be defined as the sample names.
  }
  \item{group.vector}{
    A numeric vector specifying predefined risk groups for the patients. This is optional.
  }
   \item{method}{
   A character string specifying the method for defining risk groups. Possible options are:
   \itemize{
    \item \code{"min.pval"}: Defines risk groups based on the minimum log-rank p-value.
    \item \code{"med.pval"}: Defines risk groups based on the median of the log-rank p-values.
    \item \code{"class.probs"}: Defines risk groups based on the classification probabilities from the model.
   }
   If \code{NULL}, the default method is \code{"class.probs"}. 
  }
  \item{nboot}{
    An integer specifying the number of bootstrap iterations for risk score calculation. Default is 50.
  }
}
\details{
A multivariate Cox regression is trained to select a subset of genes significantly associated with the risk and to estimate a risk score based on these risk markers. 
The algorithm considered is based on UNICOX, a regularized multivariate Cox regression model (see Tibshirani et al., 2009 for more details). In this predictor, the variables are penalized individually using an \eqn{L_1} norm term which allow us to keep more relevant genes correlated with risk than in Lasso. The Lasso model selects only one representative gene randomly from the set of correlated genes. The optimal value for the lambda parameter as well as the risk score are estimated using a double nested crossvalidation strategy. Finally, the risk score allow us to stratify the whole set of patients according to their risks. Three algorithms are implemented to estimate the optimal threshold that classifies the patients in risk groups. "min.pval" determines the optimal threshold by minimization of the log-rank p-value statistics, that is by  maximization of the separability between the K-M curves for the high and low risk groups, see  (Martinez-Romero et al., 2018). When several local minima arise this may be sample dependent and unstable. To avoid this problem, "med.pval" estimates the optimal threshold as the median of the lower 10th percentile logrank p-values. The lower 10th percentile selects the smallest values from the p-value distribution corresponding to intermediate risk patients that are on the boundary between both groups. This interval is more robust than a single minimum and provides good experimental results for a large variety of problems tested. The median threshold in this interval may change from one iteration to another because the distribution of p-values for patients with intermediate risk may change due to sample variations. Finally, "class.probs" implements a bootstrap strategy for the patients corresponding to the lower 10th percentile p-values and estimates a robust threshold to stratify the patients. It estimates also a membership probability of classification.
}
\value{
  The function generate 5 plots:
  \item{•}{Curve of log-rank p-values for each possible splitting of patients in two risk groups considering their risk score. The green line shows the robust minimum selected which provides a threshold to stratify the patients in two risk groups. The blue/red lines show the thresholds to stratify the patients in three risk groups (low, intermediate and high).}
  \item{•}{Risk plot obtained by a multivariate Cox regression model (UNICOX) considering the optimal subset of genes correlated with the risk. The green and blue/red lines show the thresholds to stratify the patients in two or three risk groups respectively.}
  \item{•}{Log-rank p-value versus number of genes selected by the model.}
  \item{•}{Genes ranked by the value of the regression coefficients in the multivariate Cox model (UNICOX). Genes with a significant p-value are marked with asterisks following the standard statistical convention. Similarly, genes with a significance level between 0.1 and 0.05 are marked with a period.}
  \item{•}{Kaplan Meier plot: Blue (low risk) and red (high risk). Risk groups are obtained from the stratification of the risk curve (plot 2).} \cr
  The function also returns an object containing:
  \item{cv_risk_score}{Risk score prediction for the training set using a double nested crossvalidated strategy}
  \item{cv_normalized_risk}{Normalized risk score in the interval (0,100).}
  \item{table_genes_selected}{Data frame with the following columns: The names for the genes selected by the Cox regression, the beta coefficients for the optimal multivariate Cox regression fitted to the training set, the Hazard Ratio for each gene and the p-value for the univariate log-rank statistical test. Genes are shown by descending order of the HR index.}
  \item{table_genes_selected_extended}{Table with the same format as table_genes_selected. A search for local minima within a 5\% range of the selected minimum is performed. The goal is expanding the list of significant genes to improve biological interpretability, since the lasso penalty drastically reduces the number of significant genes.}
  \item{model.optimalLambda}{The fitted model for the optimal regularization parameter.}
  \item{groups}{Vector of classification of patients in two risk groups, high (2) or low (1).}
  \item{riskThresholds}{Thresholds that allows to stratify the test patients in three groups according to the predicted risk score: low, intermediate and high risk.}
  \item{range.risk}{Range of the unscaled risk score in the training set.}
  \item{list.models}{List of models tested for different values of the regularization parameter.}
  \item{evaluation.models}{Data frame that provides several metrics for each model evaluated. The lambda column provides the regularization parameter for the multivariate Cox regression adjusted, the number_features gives the number of genes selected by this model, c.index and se.c.index the concordance index and the standard deviation for the risk prediction and finally, the p_value_c.index and the logrank_p_value give the p-values for the the concordance index and the log-rank statistics respectively. Models are shown by ascending order of the log-rank p-value and the best one is marked with two asterisks.}
  \item{betasplot}{Dataset used to create the plot of genes ranked according to the regression coefficients in the multivariate Cox model (UNICOX).}
  \item{plot_values}{A list containing Kaplan-Meier fit results, logrank p-value, and hazard ratio.}
  \item{membership_prob}{If method "class.probs" is selected a table with two columns is returned. The first one is the probability of classification to the low risk group while the second one is the membership probability to the high risk group.}
}
\examples{
  data(patientRisk)

  # Survival times should be provided in YEARS
  time <- mPheno$time
  names(time) <- rownames(mPheno)
  status <- mPheno$status
  names(status) <- rownames(mPheno)

  # Pred_ER.IHC$genes is the subset of genes to be tested. In our case study,
  # it is the list of genes related to the ER clinical variable that was
  # obtained using the function **genePheno()**.
  geneList <- names(Pred_ER.IHC$genes)
  
  # Next, the expression matrix for the list of genes selected is obtained.
  mExprSelectedGenes <- mExprs[match(geneList, rownames(mExprs)), ]

  # Training of the multivariate COX model. Provide the expression matrix
  # (genes as rows and samples as columns) for the list of genes selected,
  # the time and the status vectors, and the method to stratify the patients
  # (select one of these methods: `min.pval`, `med.pval`, `class.probs`).
  set.seed(5)
  multivariate_risk_predictor <- patientRisk(mExprSelectedGenes, time, status, method = "class.probs")

  # Generate the plots again
  asuri:::plotLogRank(multivariate_risk_predictor)
  asuri:::plotSigmoid(multivariate_risk_predictor)
  asuri:::plotLambda(multivariate_risk_predictor)
  asuri:::plotBetas(multivariate_risk_predictor)
  asuri:::plotKM(multivariate_risk_predictor)
}
\references{
  Martinez-Romero, J., Bueno-Fortes, S., Martín-Merino, M. \emph{et al.} (2018) \bold{Survival marker genes of colorectal cancer derived from consistent transcriptomic profiling}. \emph{BMC Genomics} 19 (Suppl 8), 857 . \doi{10.1186/s12864-018-5193-9}

  Santiago Bueno-Fortes, Alberto Berral-Gonzalez, José Manuel Sánchez-Santos, Manuel Martin-Merino, Javier De Las Rivas (2023) \bold{Identification of a gene expression signature associated with breast cancer survival and risk that improves clinical genomic platforms}, \emph{Bioinformatics Advances}, Volume 3, Issue 1, vbad037, \doi{10.1093/bioadv/vbad037}
}
