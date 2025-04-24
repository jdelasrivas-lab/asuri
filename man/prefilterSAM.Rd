\name{prefilterSAM}
\alias{prefilterSAM}
\title{
  prefilterSAM
}
\description{
  R function to pre-filter noisy genes, keeping the most relevant ones for microarrays or RNAseq data. These genes can be used later for further analysis.
}
\usage{
  prefilterSAM(mExpr, groups_vector, FDRfilter = 0.05, percentageFilter = 80, iter = 100)
}
\arguments{
  \item{mExpr}{
    A matrix where rows represent genes and columns represent samples. Each entry is the expression value of a gene in a sample.
  }
  \item{groups_vector}{
    A binary vector indicating the group assignment for each sample.
  }
  \item{FDRfilter}{
    A numeric value indicating the FDR threshold for selecting significant genes. Default is 0.05.
  }
  \item{iter}{
    The number of iterations for bootstrapping. Default is 100.
  }
  \item{percentageFilter}{
    A numeric value indicating the percentage of iterations a gene must appear in to be considered significant. Default is 80.
  }
  }
\details{
  This function implements SAM (Schwender H., 2022) robust diferential expression analysis based on bootstrap . It helps to remove noisy genes reducing the computational complexity of further analysis.
}
\value{
   An ordered vector with the names of differentially expressed genes between the categories of the grouping vector.
  \item{result}{
    A list of DE genes ordered by SAM d.value and filtered by percentageFilter.
  }  
}
\references{
  Schwender H (2022). \bold{_siggenes: Multiple Testing using SAM and Efron's Empirical Bayes Approaches_.} R package version 1.72.0.
  
  Martinez-Romero, J., Bueno-Fortes, S., Martín-Merino, M. \emph{et al.} (2018) \bold{Survival marker genes of colorectal cancer derived from consistent transcriptomic profiling}. \emph{BMC Genomics} 19 (Suppl 8), 857 . \doi{10.1186/s12864-018-5193-9}

  Santiago Bueno-Fortes, Alberto Berral-Gonzalez, José Manuel Sánchez-Santos, Manuel Martin-Merino, Javier De Las Rivas (2023) \bold{Identification of a gene expression signature associated with breast cancer survival and risk that improves clinical genomic platforms}, \emph{Bioinformatics Advances}, Volume 3, Issue 1, vbad037, \doi{10.1093/bioadv/vbad037}
}
\examples{ 
  # Bootstrapped differential expression based on SAM.
  # Parameters: FDR = 0.05, iter = 100, percentage filter = 80 %
  # CAUTION: if the data have a high number of genes this function will take several minutes to compute.
  data(prefilterSAM)

  set.seed(5)
  DE_list_genes <- prefilterSAM(mExprs, mPheno$ER.IHC)
}
