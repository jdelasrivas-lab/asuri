#' Identify Predictive Genes for a Phenotype
#'
#' This function implements robust algorithms to obtain a list of genes 
#' associated to a given clinical variable. It is based on the elastic net 
#' algorithm and the robustness and reproducibility of the subset of genes is 
#' improved using a bootstrap strategy combined with ensemble methods.
#'
#' @param seData SummarizedExperiment object with the normalized expression 
#' data and the phenotypic data in colData.
#' @param DEgenes Vector containing the genes to be used. Expected to be in the 
#' same format as the rows of the assay(seData). Usually this vector is the 
#' result of running prefilterSAM().
#' @param vectorGroups Clinical variable or phenotypic variable tested. It must 
#' be provided as a numeric binary vector.
#' @param vectorSampleID Vector containing the sample names in the same order 
#' as in assay(seData).
#' @param iter Number of bootstrap iterations (default: 100, should be 
#' changed if the function takes too long to execute).
#' @param numberOfFolds Number of folds to implement nested cross-validation. 
#' By default 5.

#'
#' @details
#' This function implements a robust version of the elastic net algorithm 
#' proposed by Tibshirani (Tibshirani et al., 2009). This algorithm considers a 
#' penalty term to avoid overfitting that is a convex combination of the 
#' \out{L<sub>2</sub>} norm (ridge regression) and \out{L<sub>1</sub>} 
#' (Lasso regression). When the alpha parameter is 1, the regularization term 
#' perfoms similarly to Lasso and minimizes the number of non-null 
#' coefficients. If a subset of features are slightly correlated Lasso selects 
#' only one of them randomly. To avoid this extreme behavior the alpha 
#' parameter is set up to 0.75 that includes more relevant variables than 
#' Lasso and improves the prediction accuracy. Besides, this choice will help 
#' to improve the stability and to reduce the variance in the feature selection 
#' process. In order to improve the robustness and reproducibility of the gene 
#' signature discovered, a bootstrap strategy is implemented. The patients are 
#' resampled with replacement giving rise to B replicates. For each replicate, 
#' a gene signature is obtained using double nested cross-validation to avoid 
#' overfitting. The final gene list is built as an ensemble of lists, 
#' considereing several metrics that evaluate the stability, the robustness 
#' and the predictive power of each gene. See (Martinez-Romero et al., 2018) 
#' for more details.
#'
#' @return
#' A list containing the following elements:
#'  - \code{genes}: A list of genes ranked according to the degree of 
#'  association with the clinical or phenotypic variable tested.
#'  - \code{listCoeff}: A list with the beta regression coefficients and the 
#'  AUC score for each bootstrap iteration.
#'  - \code{stability}: Gene selection probability estimated by bootstrap 
#'  (the number of times discovered over "n" iterations).
#'  - \code{betasMedian}: Median of the beta coefficients over the B replicates.
#'  - \code{betasMean}:  Mean of the beta coefficients over the B replicates.
#'  - \code{betasTable}: Table of genes ordered by decreasing value of the 
#'  stability coefficient. Contains several metrics: the stability index, 
#'  the mean and the median of the beta coefficients.
#'
#' @examples
#' data(seBRCA)
#' 
#' # prefilterSAM ---
#' groupsVector <- SummarizedExperiment::colData(seBRCA)$ER.IHC
#' set.seed(5)
#' DE_list_genes <- prefilterSAM(seBRCA, groupsVector)
#' 
#' # genePheno ---
#' vectorSampleID <- rownames(SummarizedExperiment::colData(seBRCA))
#' vectorGroups <- SummarizedExperiment::colData(seBRCA)$ER.IHC
#' 
#' Pred_ER.IHC <- genePheno(seBRCA, DE_list_genes, vectorGroups, vectorSampleID)
#' 
#' # Pred_ER.IHC is an output object with the list of genes that show a 
#' # significant correlation with the clinical variable. Since a bootstrap is 
#' # performed, the results of how many times across iterations a gene is found 
#' # significant are reported as *stability* (in relative numbers 0-1, 1=100%) 
#' # and the *beta values* from the regression across iterations are also 
#' # provided as *betaMedian* and *betaMean* :
#' 
#' names(Pred_ER.IHC)
#' # [1] "genes" "listCoeff" "stability" "betasMedian" "betasMean" "betasTable"
#'
#' @references 
#' \insertRef{martinezromero2018}{asuri}
#' \insertRef{BuenoFortes2023}{asuri}
#' 
#' @export

genePheno <- function(seData, DEgenes, vectorGroups, vectorSampleID, 
                        iter = 100, numberOfFolds = 5) {
    # Validating input types
    if (!is(seData, "SummarizedExperiment")) {
      stop("SEdata must be a 'SummarizedExperiment'.")
    }
    if (!is.character(DEgenes)) {
      stop("'DEgenes' must be a character vector.")
    }
    if (!is.numeric(vectorGroups)) {
        stop("'vectorGroups' must be a numeric vector.")
    }
    if (!is.character(vectorSampleID)) {
        stop("'vectorSampleID' must be a character vector.")
    }
    if (!is.numeric(iter) || iter <= 0 || iter %% 1 != 0) {
        stop("'iter' must be a positive integer.")
    }
    if (!is.numeric(numberOfFolds) || numberOfFolds <= 1 || 
        numberOfFolds %% 1 != 0) {
      stop("'numberOfFolds' must be an integer greater than 1.")
    }

    mExpr <- assay(seData)
    mExpr <- mExpr[match(DEgenes, rownames(mExpr)), ]
    mExpr <- t(mExpr)
    
    # Validating dimensions and consistency
    if (nrow(mExpr) != length(vectorGroups)) {
        stop("The length of 'vectorGroups' must equal the number of rows ", 
             "in 'mExpr' (samples).")
    }
    if (nrow(mExpr) != length(vectorSampleID)) {
        stop("The length of 'vectorSampleID' must equal the number of rows ",
             "in 'mExpr' (samples).")
    }
    if (is.null(colnames(mExpr)) || is.null(rownames(mExpr))) {
        stop("'mExpr' must have row names (sample IDs) and ", 
             "column names (gene names).")
    }

    
    # n.genes <- dim(mExpr)[2]
    n.samples <- dim(mExpr)[1]

    list <- NULL
    outp <- NULL
    pb <- txtProgressBar(min = 0, max = iter,  style = 3, 
                         width = 50, char = "=") 
    init <- numeric(iter)
    end <- numeric(iter)
    for (i in seq(1, iter)) {
      init[i] <- Sys.time()
      sampl <- sample(seq(1, n.samples), size = n.samples, replace = TRUE)
      NOsampl <- setdiff(seq(1, n.samples), unique(sampl))
      # for each time a sample is taken
      mExpr_i <- mExpr[sampl, ]
      vectorGroups_i <- vectorGroups[sampl]
      # vectorSampleID_i <- vectorSampleID[sampl]
      # outersect items
      mExpr_o <- mExpr[NOsampl, ]
      # vectorGroups_o <- vectorGroups[NOsampl]
      # vectorSampleID_o <- vectorSampleID[NOsampl]

      # calling predictor: training
      object_cv_glmnet_train <- cv.glmnet(
          x = mExpr_i,
          y = vectorGroups_i,
          nfolds = numberOfFolds,
          type.measure = "auc",
          alpha = 0.75,
          family = "binomial"
      )
      # calling predictor: predicting, samples not used in bootstrap
      object_cv_glmnet_coeff <- predict(
          object = object_cv_glmnet_train,
          newmat = mExpr_o,
          type = "coeff",
          s = object_cv_glmnet_train$lambda.1se
      )
      # AUC measuring predictive power
      object_cv_glmnet_response <- predict(
          object = object_cv_glmnet_train,
          newx = mExpr_o,
          type = "response",
          s = object_cv_glmnet_train$lambda.1se
      )
      x <- NULL
      # cumulative matrix of beta values picturing each probeset 
      # predictive power. AUC value is also stored
      x$coeff <- object_cv_glmnet_coeff[object_cv_glmnet_coeff[, 1] != 0, ]
      x$coeff <- x$coeff[2:length(x$coeff)]
      outp$genes <- c(outp$genes, names(x$coeff))

      auc_prediction <- prediction(
          as.double(object_cv_glmnet_response[, 1]),
          vectorGroups[match(rownames(object_cv_glmnet_response), 
                             as.character(vectorSampleID))]
      )

      x$auc <- as.numeric((performance(auc_prediction, "auc"))@y.values)
      list[[i]] <- x
      end[i] <- Sys.time()
      setTxtProgressBar(pb, i)
      # time <- round(lubridate::seconds_to_period(sum(end - init)), 0)
      # 
      # # Estimated remaining time based on the
      # # mean time that took to run the previous iterations
      # est <- iter * (mean(end[end != 0] - init[init != 0])) - time
      # remainining <- round(lubridate::seconds_to_period(est), 0)
      # text_msg <- paste(" // Execution time:", time, 
      #                   " // Estimated time remaining:", remainining)
      # message(text_msg, "")
    }
    close(pb)
    
    outp$listCoeff <- list
    outp$genes <- table(outp$genes)
    outp$genes <- outp$genes[outp$genes > 0.1 * iter]
    outp$stability <- outp$genes / iter

    betasM <- matrix(NA, nrow = length(outp$genes), ncol = iter)
    colnames(betasM) <- seq(1, iter)
    rownames(betasM) <- names(outp$genes)
    for (i in seq(1, iter)) {
        matching <- match(names(list[[i]]$coeff), rownames(betasM))
        index <- as.vector(na.omit(matching))
        betasM[index, i] <- as.vector(list[[i]]$coeff[!is.na(matching)])
    }
    betasMedian <- apply(betasM, 1, function(x) median(na.omit(x)))
    # confidence interval beta as col in table
    outp$betasMedian <- betasMedian

    betasMean <- apply(betasM, 1, function(x) mean(na.omit(x)))
    outp$betasMean <- betasMean

    betasTable <- cbind(outp$stability, betasMedian, betasMean)
    colnames(betasTable) <- c("stability", "betasMedian", "betasMean")
    outp$betasTable <- betasTable[order(betasTable[, 1], decreasing = TRUE), ]

    outp
}
