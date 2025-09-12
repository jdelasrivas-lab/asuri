#' Pre-filter Genes using SAM (Significance Analysis of Microarrays)
#'
#' This function performs pre-filtering of genes based on Significance Analysis 
#' of Microarrays (SAM). It runs a bootstrapping procedure to select 
#' significant genes based on the False Discovery Rate (FDR) and allows 
#' filtering by the percentage of times a gene is selected across multiple 
#' iterations.
#'
#' @param seData SummarizedExperiment object with the normalized expression 
#' data and the phenotypic data in colData.
#' @param groupsVector A binary vector indicating group assignment the samples.
#' @param FDRfilter A numeric value indicating the FDR threshold for selecting 
#' significant genes. Default is 0.05.
#' @param iter The number of iterations for bootstrapping. Default is 100.
#' @param percentageFilter A numeric value indicating the percentage of 
#' iterations a gene must appear in to be considered significant. Default is 80.
#'
#' @details
#' This function implements SAM (Schwender H., 2022) robust diferential 
#' expression analysis based on bootstrap . It helps to remove noisy genes 
#' reducing the computational complexity of further analysis. The function 
#' uses a bootstrapping approach, where in each iteration a random sample is 
#' drawn from the input data (with replacement), and the SAM algorithm is 
#' applied to select significant genes. The genes that appear as significant 
#' in a specified percentage of iterations (controlled by `percentageFilter`) 
#' are retained.
#'
#' @return An ordered vector with the names of differentially expressed genes 
#' between the categories of the grouping vector. A list of DE genes ordered 
#' by SAM d.value and filtered by percentageFilter.
#'
#' @examples
#' data(seBRCA)
#' 
#' # Bootstrapped differential expression based on SAM.
#' # Parameters: FDR = 0.05, iter = 100, percentage filter = 80 %
#' # CAUTION: if the data have a high number of genes this function will take 
#' # several minutes to compute.
#' 
#' groupsVector <- SummarizedExperiment::colData(seBRCA)$ER.IHC
#' 
#' set.seed(5)
#' ex_prefilterSAM <- prefilterSAM(seBRCA, groupsVector, iter = 25)
#' 
#' # NOTE: For consistent results with the vignettes and example data, use 
#' # default parameters (e.g., iter = 100).
#' 
#' @references
#' \itemize{
#'   \item{\insertRef{schwender2025siggenes}{asuri}}
#'   \item{\insertRef{martinezromero2018}{asuri}} 
#'   \item{\insertRef{BuenoFortes2023}{asuri}}
#' } 
#' 
#' @importFrom SummarizedExperiment assay
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom siggenes sam findDelta list.siggenes
#' @importFrom spsUtil quiet
#' @importFrom lubridate seconds_to_period
#' 
#' @export
prefilterSAM <- function(seData, groupsVector, FDRfilter = 0.05, 
                         iter = 100, percentageFilter = 80) {
    if (!is(seData, "SummarizedExperiment")) {
      stop("SEdata must be a 'SummarizedExperiment'.")
    }
    mExpr <- assay(seData)  
    
    # Error control: Ensure matrix dimensions match the group vector length
    if (dim(mExpr)[2] != length(groupsVector)) {
        stop("Different number of samples in SummarizedExperiment object and",
             "groupsVector. Please check their lengths.")
    }

    # Error control: Check for presence of column and row names in mExpr
    if (is.null(colnames(mExpr))) {
        stop("The SummarizedExperiment does not have column names.")
    }
    if (is.null(rownames(mExpr))) {
        stop("The SummarizedExperiment does not have row names.")
    }

    # Warning if any rownames are NA, and remove rows with NA gene names
    if (any(is.na(rownames(mExpr)))) {
        warning("Some row names are NA. Genes with NA names have been ",
                "removed. Please assign valid names if you don't want ",
                "them removed.", immediate. = TRUE)
        mExpr <- mExpr[!is.na(rownames(mExpr)), ]
    }
    ## bootstrap 100 samples
    # n.genes <- dim(mExpr)[1]
    n.samples <- dim(mExpr)[2]

    list.genes <- NULL
    message(Sys.time())
    # lista <- NULL
    #
    pb <- txtProgressBar(min = 0, max = iter,  style = 3, 
                         width = 50, char = "=")
    init <- numeric(iter)
    end <- numeric(iter)
    
    for (i in seq(1, iter)) {
      init[i] <- Sys.time()
      sampl <- sample(seq(1, n.samples), size = n.samples, replace = TRUE)
      # checking iterations
      # list of 500 vectors with relevant names
      # using a restrictive delta
      mExpr2 <- mExpr[, sampl]
      groupsVector2 <- groupsVector[sampl]

      samR <- try(
        sam(mExpr2, groupsVector2, method = d.stat, var.equal = FALSE),
        silent = TRUE
      )
      if (inherits(samR, "try-error")) next
      
      # extracting best genes by FDR
      delta <- try(
        spsUtil::quiet(findDelta(samR, fdr = FDRfilter)),
        silent = TRUE
      )
      if (inherits(delta, "try-error")) next
      
      delta <- unlist(delta)[1]

      new_genes <- try(
        siggenes::list.siggenes(samR, delta),
        silent = TRUE
      )
      
      if (inherits(new_genes, "try-error")) next
      list.genes <- c(list.genes, new_genes)
      # incidence as number of times it shows as significative value, 
      # is what it's retourned as a table
      end[i] <- Sys.time()
      setTxtProgressBar(pb, i)
      # time <- round(lubridate::seconds_to_period(sum(end - init)), 0)
      # 
      # # Estimated remaining time based on the
      # # mean time that took to run the previous iterations
      # est <- iter * (mean(end[end != 0] - init[init != 0])) - time
      # remainining <- round(lubridate::seconds_to_period(est), 0)
      # 
      # text_msg <- paste(" // Execution time:", time, 
      #                   " // Estimated time remaining:", remainining)
      # message(text_msg, "")
    }
    close(pb)
    message(Sys.time())
    list.genes <- factor(list.genes, levels = unique(list.genes))
    result <- names(table(list.genes)[table(list.genes) >= 
                                        (iter * (percentageFilter / 100))])
    return(result)
}
