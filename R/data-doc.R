#' Example gene expression matrix
#'
#' The dataset included in ASURI as a study case is intended to facilitate the 
#' use of the package. The dataset was obtained from several GSE series from 
#' [GEO database](http://www.ncbi.nlm.nih.gov/geo/) corresponding to breast 
#' cancer (BRCC) samples, where expression and survival time were 
#' reported [@Bueno2023]. Each sample corresponds to genome-wide expression 
#' profiles of BRCC primary tumor samples hybridized on the transcriptomic 
#' platform: *Affymetrix* *HGU133 Plus2.0*; which is a high-density microarray 
#' expression platform containing 594,000 oligo-nucleotide probes (organized 
#' in probe-sets) that we mapped to ENSEMBL genes using the 
#' *hgu133plus2hsensgcdf* CDF package, obtained from 
#' [BRAINARRAY](http://mbni.org/customcdf/24.0.0/ensg.download/hgu133plus2hsensgcdf_24.0.0.tar.gz). 
#' The expression signal of the samples was normalized using fRMA [@fRMA] 
#' or RMA [@RMA] and Combat [@combat], as described in [@Bueno2023].
#' @usage data(mExprs)
#' @format A matrix with 1070 genes as rows and 200 samples as columns.
#' @source Derived from several GSE series from GEO as explained in [@Bueno2023].
"mExprs" 

#' Example phenotypic expression matrix
#'
#' Phenotypic matrix that matches de expression dataset included in ASURI 
#' package as a study case. This matrix contains the time, status, age, and
#' different immunohistochemistry markers used in clinical practice to diagnose
#' breast cancer. It is also included which GSE series each samples is from.
#' @usage data(mPheno)
#' @format A matrix with 200 rows as samples and 8 features as columns.
#' @source Derived from several GSE series from GEO as explained in [@Bueno2023].
"mPheno"