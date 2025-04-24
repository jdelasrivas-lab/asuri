#' Pre-filter Genes using SAM (Significance Analysis of Microarrays)
#'
#' This function performs pre-filtering of genes based on Significance Analysis 
#' of Microarrays (SAM). It runs a bootstrapping procedure to select 
#' significant genes based on the False Discovery Rate (FDR) and allows 
#' filtering by the percentage of times a gene is selected across multiple 
#' iterations.
#'
#' @param mExpr A matrix where rows represent genes and columns represent 
#' samples. Each entry is the expression value of a gene in a sample.
#' @param groups_vector A vector indicating group assignment for each sample.
#' @param FDRfilter A numeric value indicating the FDR threshold for selecting 
#' significant genes. Default is 0.05.
#' @param iter The number of iterations for bootstrapping. Default is 100.
#' @param percentageFilter A numeric value indicating the percentage of 
#' iterations a gene must appear in to be considered significant. Default is 80.
#'
#' @details
#' The function uses a bootstrapping approach, where in each iteration a 
#' random sample is drawn from the input data (with replacement), and the SAM 
#' algorithm is applied to select significant genes. The genes that appear as 
#' significant in a specified percentage of iterations (controlled by 
#' `percentageFilter`) are retained.
#'
#' @return A character vector of significant gene names.
#'
#' @examples
#' # Bootstrapped differential expression based on SAM.
#' # Parameters: FDR = 0.05, iter = 100, percentage filter = 80 %
#' # CAUTION: if the data have a high number of genes this function will take several minutes to compute.
#' data(prefilterSAM)
#' 
#' set.seed(5)
#' DE_list_genes <- prefilterSAM(mExprs, mPheno$ER.IHC)
#'
#' @export
prefilterSAM <- function(mExpr, groups_vector, FDRfilter = 0.05, 
                         iter = 100, percentageFilter = 80) {
  
    # Error control: Ensure matrix dimensions match the group vector length
    if (dim(mExpr)[2] != length(groups_vector)) {
        stop("The number of columns in mExpr (samples) must match the length ",
             "of the groups_vector.")
    }

    # Error control: Check for presence of column and row names in mExpr
    if (is.null(colnames(mExpr))) {
        stop("The matrix should have sample (columns) names in ",
             "colnames(mExpr).")
    }
    if (is.null(rownames(mExpr))) {
        stop("The matrix should have gene names in rownames(mExpr).")
    }

    # Warning if any rownames are NA, and remove rows with NA gene names
    if (any(is.na(rownames(mExpr)))) {
        warning("Some row names are NA. Genes with NA names have been ",
                "removed. Please assign valid names if you don't want ",
                "them removed.", immediate. = TRUE)
        mExpr <- mExpr[!is.na(rownames(mExpr)), ]
    }
    ## bootstrap 100 samples
    n.genes <- dim(mExpr)[1]
    n.samples <- dim(mExpr)[2]

    list.genes <- NULL
    message(Sys.time())
    lista <- NULL
    #
    for (i in seq(1, iter)) { ##########
        progress <- paste0("Progress: ", i, "/", iter)
        message("\r", progress)
        flush.console()
        sampl <- sample(seq(1, n.samples), size = n.samples, replace = TRUE)
        # checking iterations
        # list of 500 vectors with relevant names
        # using a restrictive delta
        mExpr2 <- mExpr[, sampl]
        groups_vector2 <- groups_vector[sampl]

        samR <- try(
            sam(mExpr2, groups_vector2, method = d.stat, var.equal = FALSE),
            silent = TRUE
        )
        if (inherits(samR, "try-error")) next
        
        # extracting best genes by FDR
        delta <- try(
            findDelta(samR, fdr = FDRfilter),
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
    }
    list.genes <- factor(list.genes, levels = unique(list.genes))
    result <- 
      names(table(list.genes)[table(list.genes) >= 
                                (iter * (percentageFilter / 100))])
}
