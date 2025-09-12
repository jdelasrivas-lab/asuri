#' Kaplan-Meier Survival Analysis Based on Gene Expression or Risk Score
#'
#' This function analyzes the ability of a gene to mark survival based on a 
#' robust version of the KM curves. The robust K-M estimator is obtained by a 
#' bootstrap strategy.
#' 
#' @param seData SummarizedExperiment object with the normalized expression 
#' data and the phenotypic data in colData. Phenotypic colData must contain 
#' the samples name in the first column and two columns with the time and the 
#' status.
#' @param time SummarizedExperiment colData column name containing the survival 
#' time in years for each sample in numeric format.
#' @param status SummarizedExperiment colData column name containing the status 
#' (censored 0 and not censored 1) for each sample.
#' @param geneName A character string with the name of the gene being analyzed.
#' @param boxplot A logical value indicating whether to generate a boxplot of 
#' gene expression by survival group (default = TRUE).
#' @param iter The number of iterations (bootstrap resampling) for calculating 
#' optimal group cutoffs (default = 100).
#' @param type Defines if the KM curve groups are computed using risk ("risk") 
#' or gene expression (default "exprs").
#' @param cut_time A numeric value specifying the cutoff time (in years) for 
#' survival analysis. All events beyond this time are treated as censored 
#' (default = 10 years).
#'
#' @details
#' This function improves the stability and robustness of the K-M estimator 
#' using a bootstrap strategy. Patients are resampled with replacement giving 
#' rise to B replicates. The K-M estimator is obtained based on the replicates 
#' as well as the confidence intervals. The patients are stratified in two risk 
#' groups by an expression threshold that optimizes the log-rank statistics, 
#' that is the separability between the Kaplan-Meier curves for each group. 
#' This function implements a novel method to find the optimal threshold 
#' avoiding the problems of instability and unbalanced classes that suffer 
#' other implementations. Besides, a membership probability for each risk group 
#' is estimated from the classification of each sample in the replicates. 
#' This membership probability allow us to reclassify patients around the gene 
#' expression threshold in a more robust way.
#' The function provides a robust estimation of the log-rank p-value and the 
#' Hazard ratio that allow us to evaluate the ability of a given gene 
#' to mark survival.
#'
#' @return Depending on the type run, the output changes:
#' \itemize{ 
#'  \item{For \code{type = exprs}, a Kaplan-Meier plot based on expression groups, a 
#'  differential expression boxplot and a plot with the membership probability 
#'  for each risk group.
#'  Additionally, an object with the following components:}
#'  \itemize{
#'    \item{\code{geneName}: A character string with the selected name of 
#'    the gene to analyze.}
#'    \item{\code{patientExpr}: The expression level of each patient for 
#'    the gene.}
#'    \item{\code{patientClass}: Vector of group classification according 
#'    to the gene expression level: 2 = high expression, 
#'    and 1 = low expression level.}
#'    \item{\code{patientClassProbality}: Vector of membership probabilities 
#'    for the classification.}
#'    \item{\code{wilcox.pvalue}: The p-value from the Wilcoxon test comparing 
#'    the two expression groups.}
#'    \item{\code{plot_values}: A list containing Kaplan-Meier fit results, 
#'    log-rank p-value, and hazard ratio.}
#'    }
#'  \item{For \code{type = risk}, a Kaplan-Meier plot based on risk groups. 
#'  Additionally, an object with the following components:}
#'  \itemize{
#'    \item{\code{geneName}: A character string with the selected name of the 
#'    gene to analyze.}
#'    \item{\code{patientExpr}: The expression level of each patient 
#'    for the gene.}
#'    \item{\code{risk_score_predicted}: A numeric vector of predicted relative 
#'    risk scores for each patient.}
#'    \item{\code{plot_values}: A list containing Kaplan-Meier fit results,
#'    log-rank p-value, and hazard ratio.}
#'    }
#'  }
#'
#' @examples
#' data(seBRCA)
#' time <- "time"
#' status <- "status"
#' geneName <- "ESR1"
#' # The TIME value must be transformed to YEARS
#' # The gene expression vector must be provided with the NAMES of each sample,
#' # that should match the time and status NAMES.
#' set.seed(5)
#' outputKM <- geneSurv(seBRCA, time, status, geneName, type = "exprs")
#' 
#' # Generate the plots again
#' ## Plots for c(type = exprs)
#' plotBoxplot(outputKM)
#' plotProbClass(outputKM)
#' plotKM(outputKM)
#' 
#' # If we instead consider to run the function as *type* = risk
#' 
#' geneName <- "BRCA1"
#' set.seed(5)
#' outputKM.TP53 <- geneSurv(seBRCA, time, status, geneName, type = "risk")
#' 
#' ## Plots for c(type = risk)
#' plotKM(outputKM.TP53)
#' 
#' @references
#' \itemize{
#'   \item{\insertRef{martinezromero2018}{asuri}} 
#'   \item{\insertRef{BuenoFortes2023}{asuri}}
#' } 
#' 
#' @importFrom SummarizedExperiment assay colData
#' @importFrom survival coxph survdiff survfit
#' @importFrom stats median pchisq predict wilcox.test
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom survcomp hazard.ratio
#' @importFrom lubridate seconds_to_period
#' 
#' @export

geneSurv <- function(seData, time, status, geneName, boxplot = TRUE, 
                     iter = 100, type = c("exprs", "risk"), cut_time = 10) {
  
  if (!is(seData, "SummarizedExperiment")) {
    stop("SEdata must be a 'SummarizedExperiment'.")
  }
  
  
  genExpr <- assay(seData)
  genExpr <- genExpr[rownames(genExpr) %in% geneName,]
  time    <- colData(seData)$time
  status  <- colData(seData)$status
  
  names(time)   <- colData(seData)[,1]
  names(status) <- colData(seData)[,1]
  
  time   <- time[names(genExpr)] 
  status <- status[names(genExpr)]  
  
  # Error control: Ensure vectors have the same length
  if (length(genExpr) != length(time)) {
    stop("Different number of expression samples and time values. Please",
         "check the SummarizedExperiment object.")
  }
  if (length(genExpr) != length(status)) {
    stop("Different number of expression samples and status values. Please",
         "check the SummarizedExperiment object.")
  }
  
  # Ensure that the names of the vectors match
  if (!identical(names(genExpr), names(time))) {
    stop("Different sample names in the expression data and time vector.",
         "Please check the SummarizedExperiment object.")
  }
  if (!identical(names(genExpr), names(status))) {
    stop("Different sample names in the expression data and status vector.",
         "Please check the SummarizedExperiment object.")
  }
  
  # Warn if time is not in years
  if (max(time) > 30) {
    sprintf("CAUTION: Time may need to be converted to years. ",
            "Otherwise, this may result in errors.")
  }
  
  # Default to "exprs" if no type is specified
  if (is.null(type)) {
    type <- "exprs"
    warning("Type not specified. Defaulting to 'exprs'.")
  } else {
    if (!type %in% c("exprs", "risk")) {
      stop("Invalid type selected. Please choose 'exprs' or 'risk'.")
    }
  }
  
  # Check that the status vector contains only two unique values
  if (length(unique(status)) != 2) {
    stop("The status vector must contain exactly two unique values ",
         "(1 for event, 0 for censored).")
  }
  # genExpr <- as.matrix(genExpr)
  geneName <- geneName
  mSurv <- cbind(time, status)
  colnames(mSurv) <- c("time", "status")
  mSurv <- as.data.frame(mSurv)
  rownames(mSurv) <- names(time)
  
  mSurv$status[mSurv$time > cut_time] <- 0
  mSurv$time[mSurv$time > cut_time] <- cut_time + 0.1
  n.samples <- length(genExpr)
  
  if (type == "exprs") {
    # for25 <- round(n.samples * 0.25)
    # for75 <- round(n.samples * 0.75)
    
    vector.exprs <- as.numeric(genExpr)
    # order.vector.exprs <- order(vector.exprs)
    
    # matrix to fill with results
    matrixgr <- matrix(0, nrow = n.samples, ncol = iter)
    rownames(matrixgr) <- names(genExpr)
    pb <- txtProgressBar(min = 0, max = iter,  style = 3, 
                         width = 50, char = "=")
    init <- numeric(iter)
    end <- numeric(iter)
    
    for (i in seq(1, iter)) { #################
      init[i] <- Sys.time()
      # sampl <- sample(seq(1, n.samples), size = n.samples, replace = TRUE)
      
      muestra <- sample(seq(1, n.samples),
                        size = n.samples,
                        replace = TRUE
      )
      
      genExpr2 <- genExpr[muestra]
      mSurv2 <- mSurv[muestra, ]
      
      g <- functionKmGroups(genExpr2, mSurv2, geneName)
      # print(g[[2]])
      g <- g[[1]]
      names(g) <- names(genExpr2)
      matrixgr[match(names(g), rownames(matrixgr)), i] <- as.vector(g)
      end[i] <- Sys.time()
      setTxtProgressBar(pb, i)
      # timer <- round(lubridate::seconds_to_period(sum(end - init)), 0)
      # 
      # # Estimated remaining time based on the
      # # mean time that took to run the previous iterations
      # est <- iter * (mean(end[end != 0] - init[init != 0])) - timer
      # remainining <- round(lubridate::seconds_to_period(est), 0)
      # 
      # text_msg <- paste(" // Execution time:", time, 
      #                   " // Estimated time remaining:", remainining)
      # message(text_msg, "")
    }
    close(pb)
    group.assignation.vector <- NULL
    group.assignation.vector$assigned_group <- apply(
      matrixgr, 1, function(x) round(mean(x[x != 0])))
    group.assignation.vector$probability <- 
      group.assignation.vector$assigned_group
    for (k in seq(1, n.samples)) {
      group.assignation.vector$probability[k] <- sum(
        matrixgr[k, ] == group.assignation.vector$assigned_group[k]) / 
        sum(matrixgr[k, ] != 0)
    }
    p.val <- NULL
    message(table(is.na(group.assignation.vector$probability)))
    log.rank.grupos.surv <- survdiff(
      Surv(time, status) ~ group.assignation.vector$assigned_group, 
      data = mSurv)
    p.val <- 1 - pchisq(log.rank.grupos.surv$chisq, 
                        length(log.rank.grupos.surv$n) - 1)
    
    
    fits1 <- survfit(
      Surv(time, status) ~ group.assignation.vector$assigned_group, 
      data = mSurv)
    names(fits1)
    
    # library hazardR
    hazardR <- hazard.ratio(
      x = (as.numeric(group.assignation.vector$assigned_group)),
      surv.time = time,
      surv.event = status
    )
    
    names(hazardR) <- c(
      "hazard.ratio", "coef", "se", "lower.ci", "upper.ci", "p.value",
      "n", "coxm", "data"
    )
    
    if (hazardR$hazard.ratio < 1) {
      # message("1/hazar.ratio was calculated")
      hazardR$hazard.ratio <- format(1/as.numeric(hazardR$hazard.ratio),3)
      hazardR$upper.ci <- format(1 / as.numeric(hazardR$lower.ci), 3)
      hazardR$lower.ci <- format(1 / as.numeric(hazardR$upper.ci), 3)
    }
    
    Ttest <- wilcox.test(
      vector.exprs[group.assignation.vector$assigned_group == 1],
      vector.exprs[group.assignation.vector$assigned_group == 2]
    )
    
    plot_values <- NULL
    plot_values$km <- list(
      "fitsKM" = fits1,
      "p.val" = p.val,
      "hazardR" = hazardR
    )
    plot_values$source <- "geneSurv-exprs"
    
    output <- NULL
    
    output$geneName <- geneName
    output$patientExpr <- vector.exprs
    output$patientClass <- group.assignation.vector$assigned_group
    output$patientClassProbality <- group.assignation.vector$probability
    output$wilcox.pvalue <- Ttest
    output$plot_values <- plot_values
    
    
    if (boxplot) {
      plotBoxplot(output)
    }
    plotProbClass(output)
    plotKM(output)
    return(output)
  }
  
  if (type == "risk") {
    fit <- coxph(Surv(mSurv$time, mSurv$status) ~ genExpr)
    fitPredict <- predict(fit, type = "risk")
    riskVal <- as.integer(fitPredict > median(fitPredict))
    fits1 <- survfit(Surv(mSurv$time, mSurv$status) ~ riskVal, data = mSurv)
    p.val <- summary(fit)$sctest[3]
    
    hazardR <- NULL
    hazardR$hazard.ratio <- exp(-fit$coefficients)
    hazardR$upper.ci <- 1 / summary(fit)$conf.int[3]
    hazardR$lower.ci <- 1 / summary(fit)$conf.int[4]
    
    plot_values <- NULL
    plot_values$km <- list(
      "fitsKM" = fits1,
      "p.val" = p.val,
      "hazardR" = hazardR
    )
    
    plot_values$source <- "geneSurv-risk"
    
    output <- NULL
    output$geneName <- geneName
    output$patientExpr <- genExpr
    output$risk_score_predicted <- fitPredict
    output$plot_values <- plot_values
    
    plotKM(output)
    
    return(output)
  }
}
