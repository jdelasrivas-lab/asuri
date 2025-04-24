\name{geneSurv}
\alias{geneSurv}
\title{
  geneSurv
}
\description{ 
This function analyzes the ability of a gene to mark survival based on a robust version of the KM curves. The robust K-M estimator is obtained by a bootstrap strategy.
}
\usage{
  geneSurv(genExpr, time, status, geneName, boxplot = T, iter = 100, type = c("exprs", "risk"), cut_time = 10)
}
\arguments{
  \item{genExpr}{
    Vector with the normalized gene expression for each sample. names(genExpr) should contains the sample names.
  }
  \item{time}{
    Numeric vector containing the survival time for each sample in years, including the sample name in names(time).
  }
  \item{status}{
    Numeric vector with the status (censored 0 and not censored 1) for each sample. names(status) shoud include the sample names.
  }
  \item{geneName}{
    A character string with the name of the gene being analyzed.
  }
   \item{boxplot}{
    A logical value indicating whether to generate a boxplot of gene expression by survival group (default = TRUE).
  }
   \item{iter}{
    The number of iterations (bootstrap resampling) for calculating optimal group cutoffs (default = 100).
  }
   \item{type}{
    Defines if the KM curve groups are computed using risk ("risk") or gene expression (default "exprs").
  }
   \item{cut_time}{
    A numeric value specifying the cutoff time (in years) for survival analysis. All events beyond this time are treated as censored (default = 10 years).
    }
}
\details{
This function improves the stability and robustness of the K-M estimator using a bootstrap strategy. Patients are resampled with replacement giving rise to B replicates. The K-M estimator is obtained based on the replicates as well as the confidence intervals. The patients are stratified in two risk groups by an expression threshold that optimizes the log-rank statistics, that is the separability between the Kaplan-Meier curves for each group. This function implements a novel method to find the optimal threshold avoiding the problems of instability and unbalanced classes that suffer other implementations.
Besides, a membership probability for each risk group is estimated from the classification of each sample in the replicates. This membership probability allow us to reclassify patients around the gene expression threshold in a more robust way.
The function provides a robust estimation of the log-rank p-value and the Hazard ratio that allow us to evaluate the ability of a given gene to mark survival.
}
\value{
Depending on the type run, the output changes. 
  \describe{
    \item{For "type = exprs":}{
  A Kaplan-Meier plot based on expression groups, a differential expression boxplot and a plot with the membership probability for each risk group. Additionally, an object with the following components:
    \describe{
      \item{geneName}{
        A character string with the name of the gene being analyzed.
      }
      \item{patientExpr}{
        The expression level of each patient for the given gene.
      }
      \item{patientClass}{
        Vector of group classification according to the gene expression level: 2 = high expression and 1 = low expression level.
      }
      \item{patientClassProbality}{
        Vector of membership probabilities for the classification.
      }
      \item{wilcox.pvalue}{
        The p-value from the Wilcoxon test comparing the two expression groups.
      }
      \item{plot_values}{
        A list containing Kaplan-Meier fit results, log-rank p-value, and hazard ratio.
      }
    }
    }
    \item{For "type = risk":}{
    A Kaplan-Meier plot based on risk groups. Additionally, an object with the following components:
    \describe{
        \item{geneName}{
          A character string with the name of the gene being analyzed.
        }
        \item{patientExpr}{
          The expression level of each patient for the given gene.
        }
        \item{risk_score_predicted}{
          A numeric vector of predicted relative risk scores for each patient.
        }
        \item{plot_values}{
          A list containing Kaplan-Meier fit results, log-rank p-value, and hazard ratio.
        }
      }
    }
  }
}
\examples{ 
  data(geneSurvExprs)
  
  genExpr <- mExprs[match("ESR1", rownames(mExprs)), ]
  time <- mPheno$time
  names(time) <- rownames(mPheno)
  status <- mPheno$status
  names(status) <- rownames(mPheno)
  
  # The TIME value must be transformed to YEARS
  # The gene expression vector must be provided with the NAMES of each sample,
  # that should match the time and status NAMES.
  set.seed(5)
  outputKM <- geneSurv(genExpr, time, status, "ESR1", type = "exprs")
                            
  # Generate the plots again
  ## Plots for c(type = exprs)
  asuri:::plotBoxplot(outputKM)
  asuri:::plotProbClass(outputKM)
  asuri:::plotKM(outputKM)
  
  # If we instead consider to run the function as *type* = risk
  data(geneSurvRisk)
  
  genExpr <- mExprs[match("BRCA1", rownames(mExprs)), ]
  time <- mPheno$time
  names(time) <- rownames(mPheno)
  status <- mPheno$status
  names(status) <- rownames(mPheno)
  set.seed(5)
  outputKM.TP53 <- geneSurv(genExpr, time, status, "BRCA1", type = "risk")

  ## Plots for c(type = risk)
  asuri:::plotKM(outputKM)
}
\references{
  Martinez-Romero, J., Bueno-Fortes, S., Martín-Merino, M. \emph{et al.} (2018) \bold{Survival marker genes of colorectal cancer derived from consistent transcriptomic profiling}. \emph{BMC Genomics} 19 (Suppl 8), 857 . \doi{10.1186/s12864-018-5193-9}

  Santiago Bueno-Fortes, Alberto Berral-Gonzalez, José Manuel Sánchez-Santos, Manuel Martin-Merino, Javier De Las Rivas (2023) \bold{Identification of a gene expression signature associated with breast cancer survival and risk that improves clinical genomic platforms}, \emph{Bioinformatics Advances}, Volume 3, Issue 1, vbad037, \doi{10.1093/bioadv/vbad037}
}

