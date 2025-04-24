#' Identify Predictive Genes for a Phenotype
#'
#' This function implements a robust elastic net algorithm to identify a 
#' subset of genes that characterize a given phenotype. It uses bootstrap 
#' resampling and cross-validation to ensure reliable feature selection 
#' and avoids overfitting.
#'
#' @param mExpr A matrix with normalized gene expression data. Rows correspond 
#' to samples, and columns correspond to genes. `rownames(mExpr)` should be set 
#' to sample names, and `colnames(mExpr)` should be set to gene names.
#' @param vectorGroups A numeric vector indicating the phenotype groups. 
#' Positive samples should be marked as `1` and negative samples as `0`.
#' @param vectorSampleID A character vector of sample IDs corresponding 
#' to the rows of `mExpr`.
#' @param iter Integer. The number of bootstrap iterations for model training 
#' and validation. Defaults to 100.
#' @param numberOfFolds Integer. The number of folds for cross-validation 
#' during training. Defaults to 5.
#'
#' @details
#' The `genePheno` function applies the elastic net regularization method to 
#' predict genes that are significantly associated with a binary phenotype. 
#' During each iteration, a bootstrap sample is used for training, and the 
#' remaining samples are used for validation. The function computes stability
#' metrics and estimates confidence intervals for selected gene coefficients.
#'
#' @return
#' A list containing the following elements:
#' \describe{
#'   \item{\code{listCoeff}}{A list of coefficient matrices for each iteration.}
#'   \item{\code{genes}}{A table of gene frequencies across iterations, 
#'   filtered by a stability threshold.}
#'   \item{\code{stability}}{A named vector indicating the stability of each 
#'   selected gene (frequency of selection / total iterations).}
#'   \item{\code{betasMedian}}{A vector of median beta coefficients for the 
#'   selected genes across iterations.}
#'   \item{\code{betasMean}}{A vector of mean beta coefficients for the 
#'   selected genes across iterations.}
#'   \item{\code{betasTable}}{A matrix summarizing stability, median, and mean 
#'   beta coefficients for the selected genes, ordered by stability.}
#' }
#'
#' @examples
#' # Simulated data example
#' mExpr <- matrix(rnorm(1000), nrow = 50, ncol = 20)
#' rownames(mExpr) <- paste0("Sample_", 1:50)
#' colnames(mExpr) <- paste0("Gene_", 1:20)
#' vectorGroups <- c(rep(1, 25), rep(0, 25))
#' vectorSampleID <- rownames(mExpr)
#'
#' result <- genePheno(mExpr, vectorGroups, vectorSampleID)
#' print(result$betasTable)
#'
#' @export

genePheno <- function(mExpr, 
                      vectorGroups, 
                      vectorSampleID, 
                      iter = 100, 
                      numberOfFolds = 5) {
    # Validating input types
    if (!is.matrix(mExpr)) {
        stop("'mExpr' must be a matrix with genes as columns ", 
               "and samples as rows.")
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

    n.genes <- dim(mExpr)[2]
    n.samples <- dim(mExpr)[1]

    list <- NULL
    outp <- NULL
    message(Sys.time())
    for (i in seq(1, iter)) {
        if (i %% 10 == 0) {
            progress <- paste0("Progress: ", i, " %")
            message("\r", progress)
            flush.console()
        }
        # printing iter
        sampl <- sample(seq(1, n.samples), size = n.samples, replace = TRUE)
        NOsampl <- setdiff(seq(1, n.samples), unique(sampl))
        # for each time a sample is taken
        mExpr_i <- mExpr[sampl, ]
        vectorGroups_i <- vectorGroups[sampl]
        vectorSampleID_i <- vectorSampleID[sampl]
        # outersect items
        mExpr_o <- mExpr[NOsampl, ]
        vectorGroups_o <- vectorGroups[NOsampl]
        vectorSampleID_o <- vectorSampleID[NOsampl]

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
    }

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
