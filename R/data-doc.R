#' Example SummarizedExperiment object
#'
#' A `SummarizedExperiment` object created from gene expression and phenotypic
#' data included in the ASURI package. The dataset included in ASURI as a stud
#' case is intended to facilitate the use of the package. The dataset was 
#' obtained from several GSE series from 
#' [GEO database](http://www.ncbi.nlm.nih.gov/geo/) corresponding to breast 
#' cancer (BRCC) samples, where expression and survival time were 
#' reported [@Bueno2023]. Each sample corresponds to genome-wide expression 
#' profiles of BRCC primary tumor samples hybridized on the transcriptomic 
#' platform: *Affymetrix* *HGU133 Plus2.0*; which is a high-density microarray 
#' expression platform containing 594,000 oligo-nucleotide probes (organized 
#' in probe-sets) that we mapped to ENSEMBL genes using the 
#' *hgu133plus2hsensgcdf* CDF package, obtained from BRAINARRAY( 
#' hgu133plus2hsensgcdf_24.0.0.tar.gz.) The expression signal of the samples 
#' was normalized using fRMA [@fRMA] or RMA [@RMA] and Combat [@combat], 
#' as described in [@Bueno2023]. It integrates 1070 genes measured across 200 
#' breast cancer samples with associated clinical annotations.
#'
#' @format A `SummarizedExperiment` object with:
#' \describe{
#'   \item{assay:}{A matrix of normalized gene expression 
#'   values (1070 genes × 200 samples).}
#'   \item{colData:}{A DataFrame with phenotypic variables for 200 samples.}
#' }
#' @usage data(seBRCA)
#' @seealso \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' @source Derived from GEO datasets as described in [@Bueno2023].
"seBRCA"