\name{genePheno}
\alias{genePheno}
\title{
  genePheno
}
\description{
This function implements robust algorithms to obtain a list of genes associated to a given clinical variable. It is based on the elastic net algorithm and the robustness and reproducibility of the subset of genes is improved using a bootstrap strategy combined with ensemble methods.
}
\usage{
  genePheno(mExpr, vectorGroups, vectorSampleID, iter = 100, numberOfFolds = 10)
}
\arguments{
  \item{mExpr}{
    A matrix with normalized gene expression data. Rows correspond to samples, and columns correspond to genes. `rownames(mExpr)` should be set to sample names, and `colnames(mExpr)` should be set to gene names.
  }
  \item{vectorGroups}{
    Clinical variable. It must be provided as a numeric binary vector.
  }
  \item{vectorSampleID}{
    Vector containing the sample names ordered as in the expression matrix.
  }
  \item{iter}{
    Number of bootstrap iterations (default: 100, should be changed if the function takes too long to execute).
  }
  \item{numberOfFolds}{
    Number of folds to implement nested cross-validation. By default 10.
  }}
\details{
  This function implements a robust version of the elastic net algorithm proposed by Tibshirani (Tibshirani et al., 2009).
This algorithm considers a penalty term to avoid overfitting that is a convex combination of the \out{L<sub>2</sub>} norm (ridge regression) and \out{L<sub>1</sub>} (Lasso regression). When the alpha parameter is 1, the regularization term perfoms similarly to Lasso and minimizes the number of non-null coefficients. If a subset of features are slightly correlated Lasso selects only one of them randomly. To avoid this extreme behavior the alpha parameter is set up to 0.75 that includes more relevant variables than Lasso and improves the prediction accuracy. Besides, this choice will help to improve the stability and to reduce the variance in the feature selection process.
In order to improve the robustness and reproducibility of the gene signature discovered, a bootstrap strategy is implemented. The patients are resampled with replacement giving rise to B replicates.
  For each replicate a gene signature is obtained using double nested cross-validation to avoid overfitting. The final gene list is built as an ensemble of lists, considereing several metrics that evaluate the stability, the robustness and the predictive power of each gene. See (Martinez-Romero et al., 2018) for more details.
}
\value{
  A list is provided with the following objects:
  \item{genes}{
    A list of genes ranked according to the degree of association with the clinical or phenotypic variable tested.
  }
  \item{listCoeff}{
    A list with the beta regression coefficients and the AUC score for each bootstrap iteration.
  }
  \item{stability}{
    Gene selection probability estimated by bootstrap (the number of times discovered over "n" iterations).
  }
  \item{betasMedian}{
    Median of the beta coefficients over the B replicates.
  }
  \item{betasMean}{
    Mean of the beta coefficients over the B replicates.
  }
  \item{betasTable}{
    Table of genes ordered by decreasing value of the stability coefficient. Contains several metrics: the stability index, the mean and the median of the beta coefficients.
  }
}
\examples{ 
  
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

  names(Pred_ER.IHC)
  # [1] "genes" "listCoeff" "stability" "betasMedian" "betasMean" "betasTable"
}
\references{
  Martinez-Romero, J., Bueno-Fortes, S., Martín-Merino, M. \emph{et al.} (2018) \bold{Survival marker genes of colorectal cancer derived from consistent transcriptomic profiling}. \emph{BMC Genomics} 19 (Suppl 8), 857 . \doi{10.1186/s12864-018-5193-9}

  Santiago Bueno-Fortes, Alberto Berral-Gonzalez, José Manuel Sánchez-Santos, Manuel Martin-Merino, Javier De Las Rivas (2023) \bold{Identification of a gene expression signature associated with breast cancer survival and risk that improves clinical genomic platforms}, \emph{Bioinformatics Advances}, Volume 3, Issue 1, vbad037, \doi{10.1093/bioadv/vbad037}
}
