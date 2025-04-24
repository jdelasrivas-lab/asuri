## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------
library(formatR)
library(knitr)
library(asuri)
knitr::opts_chunk$set(tidy.opts = list(width.cutoff = 65), tidy = TRUE, 
                      echo = TRUE, message = FALSE, warning = FALSE, results = "hide")

## ----diagramaCrop, fig.cap="Workflow diagram of the package ASURI showing the decision-making lines and the five main functions.", out.width = "100%", fig.align = "center", echo=FALSE, results="asis"----
knitr::include_graphics("figures/diagramaCrop.png")

## -----------------------------------------------------------------------------------------------------------------------------
data(mExprs)
# To visualize the dimensions of the ExpressionSet object
dim(mExprs)

## -----------------------------------------------------------------------------------------------------------------------------
data(mPheno)
# To visualize the dimensions of the phenotypic information object
dim(mPheno)

## ----echo = TRUE, results = 'hide', fig.show='hide'---------------------------------------------------------------------------
data(geneSurvExprs)
genExpr <- mExprs[match("ESR1", rownames(mExprs)), ]
time <- mPheno$time
names(time) <- rownames(mPheno)
status <- mPheno$status
names(status) <- rownames(mPheno)
# The TIME value must be transformed to YEARS
# The gene expression vector must be provided with the NAMES of each sample,
# that should match the time and status NAMES.
set.seed(5)
outputKM <- geneSurv(genExpr, time, status, "ESR1", type = "exprs")

## ----pclassesr1, fig.cap="Patient Class Probability for low exprs (89) and high exprs (111) groups).", out.width = "80%", fig.align = "center", echo=FALSE, results="asis"----
knitr::include_graphics("figures/PclassProbability.png")

## ----boxesr1, fig.cap="Boxplot: high expression vs low expression for ESR1 gene.", out.width = "80%", fig.align = "center", echo=FALSE, results="asis", results="asis"----
knitr::include_graphics("figures/expr.png")

## ----kmesr1, fig.cap="Kaplan Meier curves obtained after the assignment of the patients by bootstrap to high or low survival using the selected feature ESR1.", out.width = "80%", fig.align = "center", echo=FALSE, results="asis"----
knitr::include_graphics("figures/kmESR1.png")

## ----echo = TRUE, results = 'hide', fig.show='hide'---------------------------------------------------------------------------
# If we instead consider to run the function as *type* = risk
data(geneSurvRisk)
genExpr <- mExprs[match("BRCA1", rownames(mExprs)), ]
time <- mPheno$time
names(time) <- rownames(mPheno)
status <- mPheno$status
names(status) <- rownames(mPheno)
set.seed(5)
outputKM.TP53 <- geneSurv(genExpr, time, status, "BRCA1", type = "risk")

## ----kmrisk, fig.cap="Kaplan Meier curves obtained after the assignment of the patients by bootstrap to high or low survival using the estimated risk based on feature BRCA1.", out.width = "80%", fig.align = "center", echo=FALSE, results="asis"----
knitr::include_graphics("figures/kmtyperisk.png")

## ----echo = TRUE, results = 'hide'--------------------------------------------------------------------------------------------
# Bootstrapped differential expression based on SAM.
# Parameters: FDR = 0.05, iter = 100, percentage filter = 80 %
# CAUTION: if the data have a high number of genes this function will take several minutes to compute.
data(prefilterSAM)

set.seed(5)
DE_list_genes <- prefilterSAM(mExprs, mPheno$ER.IHC)

## ----echo = TRUE, results = 'hide'--------------------------------------------------------------------------------------------
data(genePheno)
# Gene expression matrix filtered by SAM differential expression function
mExprsDE <- mExprs[match(DE_list_genes, rownames(mExprs)), ]
dim(mExprsDE)
# [1] 865  200
tmExprsDE <- t(mExprsDE)
# Parameters: Number of bootstrap iterations: 100.
# The object `tmExprsDE` is the gene expression matrix for the subset of genes tested (samples as rows and genes as columns). The object `mPheno$ER.IHC` is a vector indicating the value per-sample of ER (as a binary variable: "0" or "1").
vectorSampleID <- as.character(rownames(mPheno))
vectorGroups <- as.numeric(mPheno$ER.IHC)
Pred_ER.IHC <- genePheno(tmExprsDE, vectorGroups, vectorSampleID)

# Pred_ER.IHC is an output object with the list of genes that show a significant correlation with the clinical variable. Since a bootstrap is performed, the results of how many times across iterations a gene is found significant are reported as *stability* (in relative numbers 0-1, 1=100%) and the *beta values* from the regression across iterations are also provided as *betaMedian* and *betaMean* :

Pred_ER.IHC$stability
Pred_ER.IHC$betasMedian
Pred_ER.IHC$betasMean
Pred_ER.IHC$betasTable

## ----echo = TRUE, results = 'hide'--------------------------------------------------------------------------------------------
names(Pred_ER.IHC)
# [1] "genes" "listCoeff" "stability" "betasMedian" "betasMean" "betasTable"

## ----echo = TRUE, results = 'hide', fig.show='hide'---------------------------------------------------------------------------
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

## ----plot1risk, fig.cap="Curve of log-rank p-values for each possible splitting of patients into 2 risk groups. The X axis represents the patient number ranked from low to high risk. The green line marks the selected cutpoint, i.e., the point that provides the threshold for stratifying patients into 2 risk groups: low risk, in blue; high risk, in red. The blue & red lines define thresholds for classifying patients into 3 risk groups: *low*, *intermediate* and *high*. These lines mark samples that fall in an *intermediate* region, where risk assignment is not very precise and samples may fluctuate in position.", out.width = "80%", fig.align = "center", echo=FALSE, results="asis"----
knitr::include_graphics("figures/pval.png")

## ----plot2risk, fig.cap="Risk Score, between 0 and 100, calculated for each patient, i.e., for each of the 200 patient tumor samples included in this study, obtained with a multivariate UNICOX regression model. The model includes as features the subset of selected genes included in the *geneList*, which are tested to find out their correlation with risk. The green and blue & red lines provide the thresholds estimated by the **patientRisk()** function to split the patients respectively in 2 or 3 groups, as described in the legend of the previous figure. Note that the implemented algorithm is capable of detecting changes in the slope of the risk curve, giving rise to risk groups of differenciated prognosis.", out.width = "80%", fig.align = "center", echo=FALSE, results="asis"----
knitr::include_graphics("figures/risk.png")

## ----plot3risk, fig.cap="Log-rank p-values of the minimum optimal versus the number of genes selected by the model. As the regularization parameter in UNICOX increases, the regression model becomes more especific in the selection of features and, therefore the number of genes that have non-zero (non-null) coefficients is reduced.", out.width = "80%", fig.align = "center", echo=FALSE, results="asis"----
knitr::include_graphics("figures/pvaldistribution.png")

## ----plot4risk, fig.cap="Plot of the genes selected by the multivariate COX model, ranked by the value of the regression coefficients in the model. Genes with *positive* values (red) increase the risk in the samples. While genes with *negative* values (blue) reduce the risk in the samples. Genes with significant p-values are marked with asterisks (** p-value<0.01, * p-value<0.05). Similarly, genes with slightly significance level are marked with a dot (\u00B7 between 0.1 and 0.05). All the analyses carried out in this part are multivariate and take into account additive interactions among the genes.", out.width = "80%", fig.align = "center", echo=FALSE, results="asis"----
knitr::include_graphics("figures/betas.png")

## ----plot5risk, fig.cap="Kaplan Meier plot: Blue (low risk) and red (high risk) of the two groups of samples with different survival obtained from the stratification of the multivariate risk curve (Figure 6).", out.width = "80%", fig.align = "center", echo=FALSE, results="asis"----
knitr::include_graphics("figures/km.png")

## ----echo = TRUE, results = 'hide', fig.show='hide'---------------------------------------------------------------------------
data(predictPatientRisk)

# Generate the validation set, mExprs_testData if necessary.
# Vector of genes (same ones used in Cox model training)
genes <- rownames(mExprSelectedGenes)

# Simulate expression data
num_samples <- 20
set.seed(5)
mExprs_testData <- matrix(rnorm(length(genes) * num_samples, mean = 10, sd = 3),
                          nrow = length(genes), ncol = num_samples)

# Assign row names (genes) and column names (samples)
rownames(mExprs_testData) <- genes
colnames(mExprs_testData) <- paste0("Sample", 1:num_samples)

set.seed(5)
risk_prediction_validation_set <- predict.patientRisk(multivariate_risk_predictor, mExprs_testData)

## ----predprisk, fig.cap="Risk score prediction for an independent set of patients using the optimal multivariate COX model trained by function patientRisk(). Green and blue/red lines correspond to the optimal thresholds that split the patients in two and three groups respectively. These thresholds were learned in training phase.", out.width = "80%", fig.align = "center", echo=FALSE, results="asis"----
knitr::include_graphics("figures/risk_plot_predict_patientRisk.png")

## ----echo = TRUE, results = 'hide'--------------------------------------------------------------------------------------------
data(predictPatientRisk)
# Example for single patient prediction: Patient fourth is selected.
mExprs_testSingleData <- data.frame(mExprs_testData[, 4])
colnames(mExprs_testSingleData) <- colnames(mExprs_testData)[4]
# Risk prediction for the optimal subset of genes selected by patientRisk function
set.seed(5)
risk_prediction_one_patient <- predict.patientRisk(multivariate_risk_predictor, mExprs_testSingleData)

## ----echo = TRUE, results = 'hide'--------------------------------------------------------------------------------------------
# Normalized patient Risk (0 100): 27.9017675117161
# The patient is classified as Low Risk 
# Low Risk interval: (0, 37.0600635839135)

## ----echo = TRUE, results = 'hide', fig.show='hide'---------------------------------------------------------------------------
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

## ----predSurv, fig.cap="Survival curve estimation for a single patient based on the risk score predicted by a multivariate COX model. The blue line shows the median survival time at which the survival probability drops to 0.5.", out.width = "80%", fig.align = "center", echo=FALSE, results="asis"----
knitr::include_graphics("figures/surv_plot_breslow.png")

## ----echo = TRUE, results = 'hide'--------------------------------------------------------------------------------------------
devtools::session_info()

