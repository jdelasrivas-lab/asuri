#' Example SummarizedExperiment object
#'
#' A `SummarizedExperiment` object created from gene expression and phenotypic
#' data included in the ASURI package. The dataset included in ASURI as a stud
#' case is intended to facilitate the use of the package. The dataset was 
#' obtained from several GSE series from 
#' \href{http://www.ncbi.nlm.nih.gov/geo/}{GEO database} corresponding to 
#' breast cancer (BRCC) samples, where expression and survival time were 
#' reported (Bueno \emph{et al.}, 2023). Each sample corresponds to 
#' genome-wide expression profiles of BRCC primary tumor samples hybridized on 
#' the transcriptomic platform: *Affymetrix* *HGU133 Plus2.0*; which is a 
#' high-density microarray expression platform containing 594,000 
#' oligo-nucleotide probes (organized in probe-sets) that we mapped to ENSEMBL 
#' genes using the *hgu133plus2hsensgcdf* CDF package, obtained from 
#' BRAINARRAY (hgu133plus2hsensgcdf_24.0.0.tar.gz.) The expression signal of 
#' the samples was normalized using fRMA (McCall,MN \emph{et al.}, 2010), or 
#' RMA (Gautier,L \emph{et al.}, 2004) and 
#' Combat (McCall,MN \emph{et al.}, 2010), as described in 
#' Bueno \emph{et al.}, 2023. It integrates 1070 genes measured across 200 
#' breast cancer samples with associated clinical annotations.
#' 
#' @docType data
#' @keywords datasets
#' @name seBRCA
#' @usage data(seBRCA)
#' @return Load a `SummarizedExperiment` object.
#' @format A `SummarizedExperiment` object with:
#' \describe{
#'   \item{assay:}{A matrix of normalized gene expression 
#'   values (1070 genes × 200 samples).}
#'   \item{colData:}{A DataFrame with phenotypic variables for 200 samples.}
#' }
#' 
#' @references
#' \itemize{
#'   \item{\insertRef{schwender2025siggenes}{asuri}}
#'   \item{\insertRef{martinezromero2018}{asuri}} 
#'   \item{\insertRef{BuenoFortes2023}{asuri}}
#'   \item{\insertRef{fRMA}{asuri}}
#'   \item{\insertRef{RMA}{asuri}}
#'   \item{\insertRef{combat}{asuri}}
#' } 
#' @seealso \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' @source Derived from GEO datasets as described in Bueno \emph{et al.}, 2023].
"seBRCA"

#' ex_prefilterSAM
#'
#' Result of running the prefilterSAM function.
#' 
#' @docType data
#' @keywords datasets
#' @name ex_prefilterSAM
#' @usage data(ex_prefilterSAM)
#' @return Load a `character` vector object.
#' @format A `character` vector object with the list of obtained genes.
#' @seealso \code{\link{prefilterSAM}}
"ex_prefilterSAM"

#' ex_genePheno
#'
#' Result of running the genePheno function.
#'   
#' @docType data
#' @keywords datasets
#' @name ex_genePheno
#' @usage data(ex_genePheno)
#' @return Load a `List` object.
#' @format A `List` object with:
#' \itemize{
#'  \item{\code{genes}: A list of genes ranked according to the degree of 
#'  association with the clinical or phenotypic variable tested.}
#'  \item{\code{listCoeff}: A list with the beta regression coefficients and 
#'  the AUC score for each bootstrap iteration.}
#'  \item{\code{stability}: Gene selection probability estimated by bootstrap 
#'  (the number of times discovered over "n" iterations).}
#'  \item{\code{betasMedian}: Median of the beta coefficients over the B 
#'  replicates.}
#'  \item{\code{betasMean}:  Mean of the beta coefficients over the B 
#'  replicates.}
#'  \item{\code{betasTable}: Table of genes ordered by decreasing value of the 
#'  stability coefficient. Contains several metrics: the stability index, 
#'  the mean and the median of the beta coefficients.}
#'  }
#' @seealso \code{\link{genePheno}}  
"ex_genePheno"

#' ex_patientRisk
#'
#' Result of running the patientRisk function.
#'
#' @docType data
#' @keywords datasets
#' @name ex_patientRisk
#' @usage data(ex_patientRisk)
#' @return Load a `List` object.
#' @format  A `List` containing the following elements:
#' \itemize{
#'  \item{\code{cv_risk_score}: Risk score prediction for the training set 
#'  using a double nested crossvalidated strategy.}
#'  \item{\code{cv_normalized_risk}: Normalized risk score in the 
#'  interval (0,100).}
#'  \item{\code{table_genes_selected}: Data frame with the following columns: 
#'  The names for the genes selected by the Cox regression, the beta 
#'  coefficients for the optimal multivariate Cox regression fitted to the 
#'  training set, the Hazard Ratio for each gene and the p-value for the 
#'  univariate log-rank statistical test. Genes are shown by descending order 
#'  of the HR index.}
#'  \item{\code{table_genes_selected_extended}: Table with the same format as 
#'  table_genes_selected. A search for local minima within a 5\% range of the 
#'  selected minimum is performed. The goal is expanding the list of 
#'  significant genes to improve biological interpretability, since the lasso 
#'  penalty drastically reduces the number of significant genes.}
#'  \item{\code{model.optimalLambda}: The fitted model for the optimal 
#'  regularization parameter.}
#'  \item{\code{groups}: Vector of classification of patients in two risk 
#'  groups, high (2) or low (1).}
#'  \item{\code{riskThresholds}: Thresholds that allows to stratify the test 
#'  patients in three groups according to the predicted risk score: low, 
#'  intermediate and high risk.}
#'  \item{\code{range.risk}: Range of the unscaled risk score in the 
#'  training set.}
#'  \item{\code{list.models}: List of models tested for different values of the
#'   regularization parameter.}
#'  \item{\code{evaluation.models}: Data frame that provides several metrics 
#'  for each model evaluated. The lambda column provides the regularization 
#'  parameter for the multivariate Cox regression adjusted, the number_features 
#'  gives the number of genes selected by this model, c.index and se.c.index 
#'  the concordance index and the standard deviation for the risk prediction 
#'  and finally, the p_value_c.index and the logrank_p_value give the p-values 
#'  for the the concordance index and the log-rank statistics respectively. 
#'  Models are shown by ascending order of the log-rank p-value and the best 
#'  one is marked with two asterisks.}
#'  \item{\code{betasplot}: Dataset used to create the plot of genes ranked 
#'  according to the regression coefficients in the multivariate Cox model.}
#'  \item{\code{plot_values}: A list containing Kaplan-Meier fit results, 
#'  logrank p-value, and hazard ratio.}
#'  \item{\code{membership_prob}: If method "class.probs" is selected a table 
#'  with two columns is returned. The first one is the probability of 
#'  classification to the low risk group while the second one is the 
#'  membership probability to the high risk group.}
#'  }
#' @seealso \code{\link{patientRisk}}  
"ex_patientRisk"

#' ex_predictPatientRisk
#'
#' Test matrix created in the PredictPatientRisk example
#' 
#' @docType data
#' @keywords datasets
#' @name ex_predictPatientRisk
#' @usage data(ex_predictPatientRisk)
#' @return Load a `data.frame` object.
#' @format A `data.frame` object with:
#' \describe{
#'   \item{data:}{A matrix of normalized gene expression 
#'   values (77 genes × 1 samples).}
#' }
#' @seealso \code{\link{predict_PatientRisk}}  
"ex_predictPatientRisk"
