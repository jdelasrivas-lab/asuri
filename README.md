<!-- README.md is generated from README.Rmd. Please edit that file -->

# asuri

<!-- badges: start -->

[![](https://img.shields.io/badge/devel%20version-0.99.12-green.svg)](https://github.com/jdelasrivas-lab/asuri)
[![](https://img.shields.io/badge/download-NA/total-blue.svg)](https://bioconductor.org/packages/stats/bioc/asuri)
[![](https://img.shields.io/badge/download-NA/month-blue.svg)](https://bioconductor.org/packages/stats/bioc/asuri)
<!-- badges: end -->

The ASURI (Analysis of SUrvival and patients RIsk prediction based on
gene signatures) package discovers marker genes that are related to risk
prediction capabilities and to a clinical variable of interest. It uses
two main steps, including subsampling glmnet and unicox. The package
implements robust functions to discover survival markers related to a
clinical phenotype and to predict a risk score, allowing to study the
patient’s risk based on the gene signatures. Several plots are provided
to visualise the relevance of the genes, the risk score, and patient
stratification, as well as a robust version of the Kaplan-Meier curves.

## Rationale

Modern medicine based on *omic* technologies provides a new approach to
the study of diseases because it renders a way to interrogate about the
role of the genes as biomolecular *markers* associated with the risk,
prognosis, or outcome of the patients. In this regard, the discovery and
validation of *survival biomarkers* associated with a given phenotype or
with a specific clinical variable is a critical step to achieve better
disease prognosis and prediction. Furthermore, accurate risk prediction
and patients stratification based on such predictions will help to
advance in personalized treatments and precision medicine. Currently,
there are not many tools available to discover genetic survival markers
or to assess the prognostic capacity of specific gene signatures.
Moreover, it is common that gene signatures discovered are not
reproducible and robust, and cannot be correlated well with clinical
phenotypes, or with the stages and outcome of the disease. Besides,
there are not easy integrated tools providing molecular based assessment
of patients *risk* and *survival*.

## Methodology

The *asuri* package provides an integrated set of functions to analyse
disease *SURVIVAL* and provide patient *RISK* predictions based on gene
signatures. The tool allows: (1) Discovery of robust and reproducible
gene lists associated with disease survival based on gene expression or
on another gene-related activity signal (**geneSurv**); (2) Discovery of
gene markers by identification of the significant association of gene
expression (or other gene-related signal) with clinical variables or
phenotypic characteristics (**genePheno**); (3) Construction of robust
patient risk predictors based on gene signatures using univariate and
multivariate approaches (**patientRisk**).

## :writing_hand: Authors

A. Berral-Gonzalez, S. Bueno-Fortes, N. Alonso-Moreda, J.M.
Sanchez-Santos, M. Martin-Merino Acera and J. De Las Rivas

-Bioinformatics and Functional Genomics Group- Cancer Research Centre
(CiC-IBMCC, USAL/CSIC/IBSAL) Salamanca (Spain)

Learn more at <http://bioinfow.dep.usal.es/>.

If you use asuri in published research, please cite:

- Bueno-Fortes S, Berral-Gonzalez A, Sanchez-Santos J, Martin-Merino M,
  De Las Rivas J (2023). [Identification of a gene expression signature
  associated with breast cancer survival and risk that improves clinical
  genomic platforms.](https://doi.org/10.1093/BIOADV/VBAD037).
  *Bioinformatics Advances*, *3*(1), vbad037. ISSN 26350041,
  <doi:10.1093/BIOADV/VBAD037>

## :arrow_double_down: Installation

You can install the stable version from Bioconductor:

``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
## BiocManager::install("BiocUpgrade") ## you may need this
BiocManager::install("asuri")
```

or the development version of asuri from
[GitHub](https://github.com/jdelasrivas-lab/asuri) with:

``` r
# install.packages("remotes")
remotes::install_github("jdelasrivas-lab/asuri")
#> Downloading GitHub repo jdelasrivas-lab/asuri@HEAD
#> 
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/tmp/RtmpAsgDCa/remotes371b6c53917da0/jdelasrivas-lab-asuri-64b8b25/DESCRIPTION’ ... OK
#> * preparing ‘asuri’:
#> * checking DESCRIPTION meta-information ... OK
#> * installing the package to process help pages
#> Loading required namespace: asuri
#> * saving partial Rd database
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * looking to see if a ‘data/datalist’ file should be added
#> * building ‘asuri_0.99.12.tar.gz’
#> Installing package into '/tmp/Rtmp8lYYVN/temp_libpath370180768e9bd'
#> (as 'lib' is unspecified)
```
