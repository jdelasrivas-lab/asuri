#' patientRisk
#' @description Ten fold crossvalidated estimation of the multivariate risk 
#' score given a group of candidate marker genes for microarray or RNAseq
#' @export
#' @param mExpr A numeric matrix of gene expression data where each row 
#' represents a gene and each column represents a sample.
#' @param time A numeric vector representing the survival time for each sample.
#' @param status A numeric vector representing the event status (1 = event, 
#' 0 = censored) for each sample.
#' @param group.vector A numeric vector specifying predefined risk groups for 
#' the patients. This is optional.
#' @param method A character string specifying the method for defining risk 
#' groups. Possible options are:
#'   \itemize{
#'     \item \code{"min.pval"}: Define risk groups based on the minimum p-value.
#'     \item \code{"med.pval"}: Define risk groups based on the median p-value.
#'     \item \code{"class.probs"}: Defines risk groups based on the 
#'     classification probabilities from the model.
#'   }
#'   If \code{NULL}, the default method is \code{"class.probs"}.
#' @param nboot An integer specifying the number of bootstrap iterations for 
#' risk score calculation. Default is 50.
#' @details \code{patientRisk} This function estimates a multivariate risk 
#' score based on univariate Cox Regression by ten fold crossvalidation.
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{cv_risk_score}}{Risk score prediction for the training set 
#'   using a double nested crossvalidated strategy.}
#'   \item{\code{cv_normalized_risk}}{Normalized risk score in the 
#'   interval (0,100).}
#'   \item{\code{table_genes_selected}}{Data frame with the following columns: 
#'   The names for the genes selected by the Cox regression, the beta 
#'   coefficients for the optimal multivariate Cox regression fitted to the 
#'   training set, the Hazard Ratio for each gene and the p-value for the 
#'   log-rank statistical test. Genes are shown by ascending order of the 
#'   HR index.}
#'   \item{\code{table_genes_selected_extended}}{Table with the same format as 
#'   table_genes_selected. A search for local minima within a 5\% range of the 
#'   selected minimum is performed. The goal is expanding the list of 
#'   significant genes to improve biological interpretability, since the lasso 
#'   penalty drastically reduces the number of significant genes.}
#'   \item{\code{model.optimalLambda}}{The fitted model for the optimal 
#'   regularization parameter.}
#'   \item{\code{groups}}{Vector of classification of patients in two risk 
#'   groups, high (2) or low (1).}
#'   \item{\code{riskThresholds}}{Thresholds that allows to stratify the test 
#'   patients in three groups according to the predicted risk score: low, 
#'   intermediate and high risk.}
#'   \item{\code{range.risk}}{Range of the unscaled risk score in the 
#'   training set.}
#'   \item{\code{list.models}}{List of models tested for different values of 
#'   the regularization parameter.}
#'   \item{\code{evaluation.models}}{Data frame that provides several metrics 
#'   for each model evaluated. The lambda column provides the regularization 
#'   parameter for the multivariate Cox regression adjusted, the 
#'   number_features gives the number of genes selected by this model, 
#'   c.index and se.c.index the concordance index and the standard deviation 
#'   for the risk prediction and finally, the p_value_c.index and the 
#'   logrank_p_value give the p-values for the the concordance index and the 
#'   log-rank statistics respectively. Models are shown by ascending order of 
#'   the log-rank p-value and the best one is marked with two asterisks.}
#'   \item{\code{betasplot}}{Dataset used to create the plot of genes ranked 
#'   according to the regression coefficients in 
#'   the multivariate Cox model (UNICOX).}
#'   \item{\code{plot_values}}{A list containing Kaplan-Meier fit results, 
#'   logrank p-value, and hazard ratio.}
#'   \item{\code{membership_prob}}{If method "class.probs" is selected a table 
#'   with two columns is returned. The first one is the probability of 
#'   classification to the low risk group while the second one is the 
#'   membership probability to the high risk group.}
#' }
#'
#' @examples
#' data(matrixData)
#' set.seed(5)
#'
#' # The TIME value must be transformed to YEARS
#' ## Names of time and status vectors are the sample names
#' time <- mPheno$time
#' names(time) <- rownames(mPheno)
#' status <- mPheno$status
#' names(status) <- rownames(mPheno)
#'
#' # Pred_Er.IHC$genes is the subset of genes to be tested. In our case study 
#' # it is the list of genes related to the ER clinical variable obtained by 
#' # function genePheno. It can be any subset of genes to analyze the 
#' # relationship with survival.
#' geneList <- names(Pred_ER.IHC$genes)
#' # Next, the expression matrix for the list of genes selected is obtained.
#' mExprSelectedGenes <- mExprs[match(geneList, rownames(mExprs)), ]
#'
#' # Training of the multivariate Cox model.
#' # You have to provide: the expression matrix (genes as rows and samples 
#' # as columns) for the list of genes selected, the time #' and status vectors.
#' multivariate_risk_predictor <- patientRisk(mExprSelectedGenes, time, status, 
#'                                            method = "class.probs")
#'
#' # Generate the plots again
#' asuri:::plotLogRank(multivariate_risk_predictor)
#' asuri:::plotSigmoid(multivariate_risk_predictor)
#' asuri:::plotLambda(multivariate_risk_predictor)
#' asuri:::plotBetas(multivariate_risk_predictor)
#' asuri:::plotKM(multivariate_risk_predictor)
#'
patientRisk <- function(mExpr,
                        time,
                        status,
                        group.vector,
                        method = NULL,
                        nboot = 50) {
    # Error control: Verify the validity of input dimensions
    if (dim(mExpr)[2] != length(time)) {
        stop("The gene expression matrix and the time vector must have ",
             "the same length.")
    }
    if (dim(mExpr)[2] != length(status)) {
        stop("The gene expression matrix and the status vector must have ",
             "the same length.")
    }

    if (!identical(colnames(mExpr), names(time))) {
        stop("The column names of mExpr must match the names ",
             "in the time vector.")
    }
    if (!identical(colnames(mExpr), names(status))) {
        stop("The column names of mExpr must match the names ",
             "in the status vector.")
    }

    # Set the method for risk analysis (default if not provided)
    if (is.null(method)) {
        warning("No method selected for risk group division. ",
                "'class.probs' has been selected as the default method.\n")
        method <- "class.probs"
    } else {
        if (!method %in% c("min.pval", "med.pval", "class.probs")) {
            stop("Incorrect method for risk group division. Choose one of the ",
                 "following methods: 'min.pval', 'med.pval', 'class.probs'.")
        }
    }

    # Verify the status vector (must have exactly two unique values)
    if (length(unique(status)) != 2) {
        stop("The status vector must contain exactly two values",
             "(event and censored).")
    }
    # defining pData matrix (train)
    mSurv <- data.frame(
        "time" = as.numeric(time),
        "status" = as.numeric(status),
        row.names = names(time),
        stringsAsFactors = FALSE
    )

    folds <- sample(cut(seq(1, length(mExpr[1, ])),
        breaks = 10,
        labels = FALSE
    )) # vector is randomized with sample()
    # Perform 10 fold cross validation
    betasMatrix <- NULL

    riskDefinitive <- matrix(data = 0, nrow = length(mExpr[1, ]), ncol = 30)
    rownames(riskDefinitive) <- names(mExpr[1, ])

    # prefit model to choose lambda range 30 lambdas 
    # min number of genes = 1/4 of total genes
    capture.output(fit <- uniCox(t(mExpr), mSurv[, 1], mSurv[, 2], 
                                 nlam = 30, del.thres = .05, max.iter = 20))
    lambMax <- max(fit$lamlist[fit$nfeatures >= ceiling(
      length(mExpr[, 1]) * 0.15)])
    # lambMax is the limit to perform the lambda optimal net of values
    listFits <- NULL
    message("Nested ten fold cross validation: Predicting the risk for ", 
            "each lambda...\n")

    ### Grid of values for the regularization parameter.
    list.of.lambdas <- seq(from = 0, to = lambMax, length.out = 30)
    for (i in seq(1, 10)) {
        progress <- paste0("Progress: ", i * 10, " %")
        message("\r", progress)
        flush.console()
        # splitting the data using folds
        testIndexes <- which(folds == i, arr.ind = TRUE)
        testData <- t(mExpr[, testIndexes])
        trainData <- t(mExpr[, -testIndexes])
        mSurvUniCox <- mSurv
        capture.output(fit <- uniCox(trainData,
            mSurvUniCox[-testIndexes, 1],
            mSurvUniCox[-testIndexes, 2],
            lamlist = list.of.lambdas,
            del.thres = .01,
            max.iter = 20
        ))
        # beta values are stored, their means will be the outcome used to 
        # order the genes by their global influence in the model
        riskDefinitive[testIndexes, ] <- as.matrix(predict.uniCox(fit,testData))
        listFits[[i]] <- list("fit" = fit, "betas" = fit$beta)
    }
    message("Risk predicted!")
    # Control of risk-matrix

    vars_near_zero <- function(x) {
        if (length(unique(x)) < 0.25 * length(x)) {
            if (max(table(x) / length(x)) >= 0.9) {
                return(0)
            }
        }
        return(1)
    }

    lambdas.risk.constant <- apply(riskDefinitive, 
                                   MARGIN = 2, FUN = "vars_near_zero")
    list.of.lambdas <- list.of.lambdas[lambdas.risk.constant == 1]
    riskDefinitive <- riskDefinitive[, lambdas.risk.constant == 1]

    # list.of.lambdas <- list.of.lambdas[colSums(riskDefinitive) != 0]
    # riskDefinitive  <- riskDefinitive[, colSums(riskDefinitive) != 0]

    # Estimation of index errors for each lambda
    # Logrank p.value
    nested.p.vals <- function_double_nested_lambda(riskDefinitive, mSurv)
    optimal.p.vals <- apply(nested.p.vals$p.vals, 2, FUN = function(x) min(x))
    # after CV optimal lambda is stored and used to re-compute groups and pvalue
    index.optimalLambda <- which.min(optimal.p.vals)

    range_optimal.p.vals_5 <- abs(
      range(optimal.p.vals)[2] - range(optimal.p.vals)[1])
    threshold_optimal.p.vals_5 <- 
      optimal.p.vals[index.optimalLambda] + 0.2 * range_optimal.p.vals_5
    optimal.p.vals_temp <- optimal.p.vals[seq(1, index.optimalLambda)]

    if (!any(optimal.p.vals_temp >= threshold_optimal.p.vals_5)) {
        index.optimalLambda_5 <- which.min(
          abs(optimal.p.vals - threshold_optimal.p.vals_5))
    } else {
        min_optimal.p.vals_temp <- min(
          optimal.p.vals_temp[optimal.p.vals_temp >= 
                                threshold_optimal.p.vals_5])
        index.optimalLambda_5 <- which.min(
          abs(optimal.p.vals - min_optimal.p.vals_temp))
        if (min(optimal.p.vals[seq(1, (index.optimalLambda_5 - 1))]) < 
            threshold_optimal.p.vals_5) {
            index.optimalLambda_5 <- which.min(
              optimal.p.vals[seq(1, (index.optimalLambda_5 - 1))])
        }
    }

    # Concordance index
    fConcordanceIndex <- function(risk_vector) {
        concordance.index(
            x = risk_vector,
            surv.time = mSurv$time,
            surv.event = mSurv$status,
            method = "noether"
        )
    }

    concordance_risk <- lapply(
      as.data.frame(riskDefinitive), FUN = "fConcordanceIndex")
    concordance_risk$c.index <- lapply(
      concordance_risk[seq(1, length(list.of.lambdas))], 
      FUN = function(x) x$c.index)
    concordance_risk$Cse <- lapply(
      concordance_risk[seq(1, length(list.of.lambdas))], 
      FUN = function(x) x$se)
    concordance_risk$Cp.value <- lapply(
      concordance_risk[seq(1, length(list.of.lambdas))], 
      FUN = function(x) x$p.value)

    # Estimation of regression coefficients for each lambda
    capture.output(fit.Lambdas <- uniCox(t(mExpr),
        mSurv[, 1],
        mSurv[, 2],
        lamlist = list.of.lambdas,
        del.thres = .05,
        max.iter = 20
    ))

    # Number of features for each lambda
    number_features <- fit.Lambdas$nfeatures

    table_models <- data.frame(
        lambda = list.of.lambdas,
        number_features = number_features,
        c.index = unlist(concordance_risk$c.index),
        se.c.index = unlist(concordance_risk$Cse),
        p_value_c.index = unlist(concordance_risk$Cp.value),
        logrank_p_value = optimal.p.vals
    )

    #### Plots for the Cox model that maximizes the log-rank test
    riskDefinitive <- riskDefinitive[, index.optimalLambda]
    p.val30a70 <- nested.p.vals$p.vals[, index.optimalLambda]

    ### Computing betas for the Cox model that maximizes log-rank test
    betas <- NULL
    capture.output(fit.optimalLambda <- uniCox(t(mExpr),
        mSurv[, 1],
        mSurv[, 2],
        lamlist = list.of.lambdas[index.optimalLambda],
        del.thres = .05,
        max.iter = 20
    ))
    risk.predict.optimal.lambda <- predict.uniCox(fit.optimalLambda, t(mExpr))
    betas <- NULL
    betas <- as.vector(fit.optimalLambda$beta)
    names(betas) <- rownames(mExpr)

    # Table list of genes selected, regression coefficients, p-values and HR 
    # for the optimal multivariate Cox model
    betas_norm0 <- betas / (sqrt(fit.optimalLambda$vx) + fit.optimalLambda$s0)
    betas_norm <- betas_norm0[abs(betas_norm0) > 0]
    # Get logrank-pvalue for each gene with non-null beta_i
    genes_tested <- rownames(mExpr[abs(betas_norm0) > 0, ])
    logrank_pvalue_genes_tested_uni <- NULL
    logrank_pvalue_genes_uni <- rep(NA, length(betas_norm0))
    names(logrank_pvalue_genes_uni) <- rownames(mExpr)

    for (gk in genes_tested) {
        univariate_cox <- coxph(Surv(time, status) ~ get(gk), 
                                data = as.data.frame(t(mExpr)))
        logrank_pvalue_genes_tested_uni <- c(logrank_pvalue_genes_tested_uni, 
                                             summary(univariate_cox)$sctest[3])
    }
    names(logrank_pvalue_genes_tested_uni) <- genes_tested
    which_betas <- abs(betas_norm0) > 0
    logrank_pvalue_genes_uni[which_betas] <- logrank_pvalue_genes_tested_uni

    # The LASSO penalty reduces the number of important genes significantly. 
    # A search for other nearby minima within 5% of the global minimum has been 
    # conducted to broaden the list of genes and improve biological 
    # interpretability.
    capture.output(fit.optimalLambda_5 <- uniCox(t(mExpr),
        mSurv[, 1],
        mSurv[, 2],
        lamlist = list.of.lambdas[index.optimalLambda_5],
        del.thres = .05,
        max.iter = 20
    ))

    betas_5 <- NULL
    betas_5 <- as.vector(fit.optimalLambda_5$beta)
    names(betas_5) <- rownames(mExpr)

    # Table list of genes selected, regression coefficients, p-values and HR 
    # for the optimal multivariate Cox model
    betas_norm0_5 <- betas_5 / (
      sqrt(fit.optimalLambda_5$vx) + fit.optimalLambda_5$s0)
    betas_norm_5 <- betas_norm0_5[abs(betas_norm0_5) > 0]
    # Get logrank-pvalue for each gene with non-null beta_i
    genes_tested_5 <- rownames(mExpr[abs(betas_norm0_5) > 0, ])
    logrank_pvalue_genes_tested_uni_5 <- NULL
    logrank_pvalue_genes_uni_5 <- rep(NA, length(betas_norm0_5))
    names(logrank_pvalue_genes_uni_5) <- rownames(mExpr)

    for (gk in genes_tested_5) {
        univariate_cox_5 <- coxph(Surv(time, status) ~ get(gk), 
                                  data = as.data.frame(t(mExpr)))
        logrank_pvalue_genes_tested_uni_5 <- 
          c(logrank_pvalue_genes_tested_uni_5, 
            summary(univariate_cox_5)$sctest[3])
    }
    names(logrank_pvalue_genes_tested_uni_5) <- genes_tested_5
    which_norm0 <- abs(betas_norm0_5) > 0
    logrank_pvalue_genes_uni_5[which_norm0] <- logrank_pvalue_genes_tested_uni_5

    # Significance codes p-vals
    significance_logrank_pvalue_uni <- symnum(logrank_pvalue_genes_uni,
        cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
        symbols   = c("***", "**", "*", "·", " "),
        na        = FALSE
    )

    table_genes_optimallambda <- data.frame(
        gene.name = names(betas_norm),
        betas = betas_norm,
        HR = exp(betas_norm),
        p_val_log_rank_uni = logrank_pvalue_genes_tested_uni,
        stringsAsFactors = FALSE
    )
    rownames(table_genes_optimallambda) <- NULL

    table_genes_optimalLambda_5 <- data.frame(
        gene.name = names(betas_norm_5),
        betas = betas_norm_5,
        HR = exp(betas_norm_5),
        p_val_log_rank_uni = logrank_pvalue_genes_tested_uni_5,
        stringsAsFactors = FALSE
    )

    # Optimal threshold selection
    nSamples <- dim(mExpr)[2]
    lowIndex0.30 <- round(0.3 * nSamples)
    highIndex0.70 <- round(0.7 * nSamples)

    orderedRiskDefinitive <- riskDefinitive[order(riskDefinitive)]
    orderednormalized_risk <- rescale(as.numeric(orderedRiskDefinitive), 
                                      to = c(0, 100))
    names(orderednormalized_risk) <- names(orderedRiskDefinitive)

    lower10 <- quantile(p.val30a70, probs = 0.1)
    lower10index <- p.val30a70 <= lower10

    lowIndexReal <- min(which(lower10index == TRUE))
    highIndexReal <- max(which(lower10index == TRUE))

    lower10index[c(lowIndexReal, highIndexReal)] <- FALSE

    if (method %in% "min.pval") {
        cutPoint <- which.min(p.val30a70)
    }
    if (method %in% "med.pval") {
        med.index <- which.min(
          abs(p.val30a70[lower10index] - median(p.val30a70[lower10index])))
        cutPoint <- which(p.val30a70 == p.val30a70[lower10index][med.index])
    }
    if (method %in% "class.probs") {
        cutPointLow <- orderednormalized_risk[lowIndexReal]
        cutPointHigh <- orderednormalized_risk[highIndexReal]

        thresholds_boot <- rep(0, 2)
        ordered_mSurv <- mSurv[order(riskDefinitive), ]
        med_threshold <- rep(0, nboot)

        membership_prob <- matrix(
            data = 0,
            ncol = 2,
            nrow = (highIndexReal - lowIndexReal) + 1
        )
        colnames(membership_prob) <- c("Low risk", "High risk")
        rownames(membership_prob) <- as.character(lowIndexReal:highIndexReal)


        lower10_normalized_risk <- 
          orderednormalized_risk[lowIndexReal:highIndexReal]

        for (i in seq(1, nboot)) {
            bsample <- unique(sort(sample(seq(1, nSamples), nSamples, 
                                          replace = TRUE)))
            orderednormalized_risk_bsample <- orderednormalized_risk[bsample]
            thresholds_boot[1] <- which.max(
              orderednormalized_risk_bsample >= cutPointLow)
            thresholds_boot[2] <- which.max(
              orderednormalized_risk_bsample >= cutPointHigh)
            nested.p.vals.b <- function_double_nested_lambda(
                data.frame(orderednormalized_risk_bsample),
                ordered_mSurv[bsample, ],
                thresholds_boot
            )
            # Threshold estimation
            nested.p.vals_lower10b <- 
              nested.p.vals.b$p.vals[thresholds_boot[1]:thresholds_boot[2]]
            lower10b <- quantile(nested.p.vals_lower10b, probs = 0.1)
            median_lower10b <- 
              median(nested.p.vals_lower10b[nested.p.vals_lower10b <= lower10b])
            med_index <- 
              which.min(abs(nested.p.vals.b$p.vals - median_lower10b))
            med_threshold[i] <- orderednormalized_risk_bsample[med_index]
            
            membership_prob[lower10_normalized_risk < med_threshold[i], 1] <- 
              membership_prob[lower10_normalized_risk < med_threshold[i], 1]+ 1
            
            membership_prob[lower10_normalized_risk >= med_threshold[i], 2] <- 
              membership_prob[lower10_normalized_risk >= med_threshold[i], 2]+ 1
        }

        membership_prob <- membership_prob / nboot
        lower10_boot_group <- 
          apply(membership_prob, MARGIN = 1, FUN = "which.max")
        
        max_membership_prob <- apply(membership_prob, MARGIN = 1, FUN = "max")
        lower10_boot_group[max_membership_prob < 0.8] <- 8
        boot_med_threshold <- mean(med_threshold)
        cutPoint <- which.min(abs(boot_med_threshold - orderednormalized_risk))
    }

    # cutPoint is used to create the groups
    definitiveGroups <- rep(0, length(riskDefinitive))
    definitiveGroups[riskDefinitive < orderedRiskDefinitive[cutPoint]] <- 1
    definitiveGroups[riskDefinitive >= orderedRiskDefinitive[cutPoint]] <- 2

    # Normalized coefficients by the Fisher Information Matrix
    orderedBetas <- betas_norm0[order(betas_norm0)]
    betasPlot <- data.frame(x = names(orderedBetas), y = orderedBetas)
    mycolor <- ifelse(orderedBetas > 0, "type1", "type2")
    betasPlot <- data.frame(betasPlot, mycolor = mycolor)
    colnames(betasPlot) <- c("gene", "beta_value", "beta_group")

    pvals <- logrank_pvalue_genes_tested_uni
    siglvl <- significance_logrank_pvalue_uni[names(pvals)]
    index.match <- match(betasPlot$gene, names(pvals))
    betasPlot$p_value <- pvals[index.match]
    betasPlot$signif <- as.character(siglvl[index.match])
    betasPlot$signif[is.na(betasPlot$signif)] <- " "
    betasPlot <- betasPlot[betasPlot$beta_value != 0, ]
    betasPlot <- betasPlot[with(betasPlot, order(-abs(beta_value), p_value)), ]

    hazardR <- hazard.ratio(
        x = definitiveGroups,
        surv.time = mSurv$time[match(colnames(mExpr), names(time))],
        surv.event = mSurv$status[match(colnames(mExpr), names(status))]
    )

    names(hazardR) <- c(
        "hazard.ratio", "coef", "se", "lower.ci", "upper.ci", "p.value",
        "n", "coxm", "data"
    )

    if (hazardR$hazard.ratio < 1) {
        # message("1/hazar.ratio was calculated")
        hazardR$hazard.ratio <- format(1 / as.numeric(hazardR$hazard.ratio), 3)
        hazardR$upper.ci <- format(1 / as.numeric(hazardR$lower.ci), 3)
        hazardR$lower.ci <- format(1 / as.numeric(hazardR$upper.ci), 3)
    }

    if (missing(group.vector)) {
        group.vector <- definitiveGroups
    } else {
        group.vector <- group.vector
    }

    mSurvKM <- mSurv[match(colnames(mExpr), rownames(mSurv)), ]
    mSurvKM$time[mSurvKM$time > 10] <- 10
    mSurvKM$status[mSurvKM$time > 10] <- 0

    fitsKM <- survfit(Surv(time, status) ~ group.vector, data = mSurvKM)
    log.rank.groups.surv <- survdiff(Surv(time, status) ~ group.vector, 
                                     data = mSurvKM)

    p.val <- pchisq(log.rank.groups.surv$chisq, 
                    length(log.rank.groups.surv$n) - 1, lower.tail = FALSE)

    normalized_risk <- rescale(as.numeric(riskDefinitive), to = c(0, 100))

    risk.predict.optimal.lambda.scaled <- 
      rescale(sort(risk.predict.optimal.lambda), to = c(0, 100))
    riskThresholds <- data.frame(
        Threshold = c(
            risk.predict.optimal.lambda.scaled[lowIndexReal],
            risk.predict.optimal.lambda.scaled[cutPoint],
            risk.predict.optimal.lambda.scaled[highIndexReal]
        ),
        Sample = c(lowIndexReal, cutPoint, highIndexReal)
    )
    rownames(riskThresholds) <- c("LowriskThreshold", 
                                  "riskCutpoint", 
                                  "HighriskThreshold")
    # Printing indexex for Cox models obtained by diff regularization parameters
    rownames(table_models) <- NULL
    rownames(table_models)[index.optimalLambda] <- 
      paste0(index.optimalLambda, " **")
    table_genes_optimallambda <- 
      table_genes_optimallambda[order(table_genes_optimallambda$HR, 
                                      decreasing = TRUE), ]
    table_genes_optimalLambda_5 <- 
      table_genes_optimalLambda_5[order(table_genes_optimalLambda_5$HR, 
                                        decreasing = TRUE), ]

    plot_values <- NULL
    plot_values$lambda <- list(
        "nested.p_vals" = nested.p.vals,
        "number.features" = number_features
    )

    plot_values$sigmoid_logrank <- list(
        "p.val30a70" = p.val30a70,
        "lower10" = lower10,
        "lowIndexReal" = lowIndexReal,
        "highIndexReal" = highIndexReal,
        "cutPoint" = cutPoint,
        "nSamples" = nSamples,
        "orderednormalized_risk" = orderednormalized_risk
    )

    plot_values$km <- list(
        "fitsKM" = fitsKM,
        "p.val" = p.val,
        "hazardR" = hazardR
    )
    plot_values$source <- "patientRisk"
    
    eval_models <- table_models[order(table_models$logrank_p_value),]
    
    rList <- list(
        "cv_risk_score" = riskDefinitive,
        "cv_normalized_risk" = normalized_risk,
        "table_genes_selected" = table_genes_optimallambda,
        "table_genes_selected_extended" = table_genes_optimalLambda_5,
        "model.optimalLambda" = fit.optimalLambda,
        "groups" = definitiveGroups,
        "riskThresholds" = riskThresholds,
        "range.risk" = range(risk.predict.optimal.lambda),
        "list.models" = fit.Lambdas,
        "evaluation.models" = eval_models,
        "betasplot" = betasPlot,
        "plot_values" = plot_values
    )


    plotLogRank(rList)

    # patient ordered by risk value
    plotSigmoid(rList)

    # plot p values by lambda
    plotLambda(rList)

    plotBetas(rList)

    plotKM(rList)

    if (method == "class.probs") {
        rList$membership_prob <- membership_prob
    }
    betasPlot
    return(rList)
}
