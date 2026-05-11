#' Stein's loss
#'
#' This function calculates Stein's loss between true and estimated positive-definite matrices,
#' which equals two times the Kullback-Leibler divergence, and is defined as follows:
#' \deqn{\mathrm{tr}(\widehat{\mathbf{M}}\mathbf{M}^{-1}) - \log |\widehat{\mathbf{M}}\mathbf{M}^{-1}| - p}
#' 
#' @param true_matrix The true \code{p} by \code{p} covariance or precision matrix.
#' @param est_matrix The estimate of the \code{p} by \code{p} covariance or precision matrix.
#'
#' @return Scalar, which is calculated Stein's loss between two matrices.
#' @references James, W. and Stein, C. (1961). Estimation with quadratic loss. In \emph{Proceedings
#'  of the Fourth Berkeley Symposium on Mathematical Statistics and Probability}, volume 1, pages
#'  361--379. University of California Press.
#' @export
#'
#' @examples
#' m1 <- matrix(c(1,0,0,1), 2, 2)
#' m2 <- matrix(c(0.8, 0, 0, 0.8), 2, 2)
#' stein_loss(m1, m2)
stein_loss <- function(true_matrix, est_matrix){
  p <- nrow(true_matrix)
  P <- est_matrix %*% chol2inv(chol(true_matrix))
  return(sum(diag(P)) - determinant(P)$modulus[1] - p)
}

#' Confusion matrix for network estimation
#' 
#' Construct confusion matrix using the true adjacency and estimated adjacency matrices of
#' the network. Can be used both directed and undirected networks.
#'
#' @param true_adj The true \code{p} by \code{p} adjacency matrix of the network.
#' @param est_adj The estimated \code{p} by \code{p} adjacency matrix of the network.
#' @param margins Boolean, add marginal sums if \code{TRUE}. The default value is \code{FALSE}.
#' @param normalize Boolean, use normalized values (proportions) instead of absolute values if
#'   \code{TRUE}. The default value is \code{FALSE}.
#' @param undirected Boolean, if \code{TRUE} (default) assumes that the adjacency matrices are
#'   symmetric.
#'
#' @return A data frame, which size depends if marginal sums are displayed or not. The default size
#'  is 2 by 2, and if marginal sums added, then 3 by 3.
#' @export
#'
#' @examples
#' true_adj <- matrix(c(0, 0, 1, 0, 0,
#'                      0, 0, 1, 1, 0,
#'                      1, 1, 0, 0, 0,
#'                      0, 1, 0, 0, 1,
#'                      0, 0, 0, 1, 0), 5, 5)
#' est_adj <- matrix(c(0, 1, 1, 0, 0,
#'                     1, 0, 0, 1, 0,
#'                     1, 0, 0, 0, 1,
#'                     0, 1, 0, 0, 1,
#'                     0, 0, 1, 1, 0), 5, 5)
#' conf_matrix(true_adj, est_adj)
conf_matrix <- function(true_adj, est_adj, margins = FALSE, normalize = FALSE,
                        undirected = TRUE) {
  same_edges <- true_adj * est_adj
  diff <- true_adj - est_adj
  summ <- true_adj + est_adj
  p <- dim(true_adj)[1]
  max_edges <- (p^2 - p)
  tp <- sum(same_edges)
  tn <- sum(summ == 0) - p
  fp <- sum(diff == -1)
  fn <- sum((same_edges - true_adj) == -1)
  P <- sum(true_adj)
  N <- max_edges - P
  EP <- sum(est_adj)
  EN <- max_edges - EP
  cm <- matrix(c(tp, fn, fp, tn), nrow = 2, byrow = TRUE,
               dimnames = list(c("True P", "True N"), c("Estim. P", "Estim. N")))
  if (margins) {
    cm <- matrix(c(tp, fn, P, fp, tn, N, EP, EN, max_edges), nrow = 3, byrow = TRUE,
                 dimnames = list(c("True P", "True N", "Sum"), c("Estim. P", "Estim. N", "Sum")))
  }
  if (undirected) {
    cm <- cm * 0.5
  }
  if (normalize) {
    cm <- matrix(c(tp/P, fn/P, fp/N, tn/N), nrow = 2, byrow = TRUE,
                 dimnames = list(c("True P", "True N"), c("Estim. P", "Estim. N")))
  }
  return(cm)
}

#' Performance scores for network estimation
#' 
#' This function calculates performance scores for network estimation. The scores are same as binary
#' classification scores as the adjacency matrices are binary matrices.
#' 
#' @param cm The 2 by 2 confusion matrix calculated using the function [conf_matrix()] or
#'   manually created 2 by 2 matrix using following order:
#'   \tabular{cc}{
#'      \eqn{TP} \tab \eqn{FN}\cr
#'      \eqn{FP} \tab \eqn{TN}
#'   }
#'
#' @return
#' A data frame, which contains following scores:
#' \item{ACC}{
#'    Accuracy, which is \eqn{(TP + TN) / (TP + TN + FN + FP)}.
#' }
#' \item{ACC_bal}{
#'    Balanced accuracy, which is \eqn{(TPR + TNR) / 2}.
#' }
#' \item{MCC}{
#'    Matthews correlation coefficient, which is \eqn{(TP \cdot TN - FP \cdot FN) / \sqrt{(TP + FP) (TP + FP) (TN + FP) (TN + FN)}}.
#' }
#' \item{F1}{
#'    \eqn{F_1}-score, which is \eqn{2(PPV \cdot TPR) / (PPV + TPR)}.
#' }
#' \item{TPR}{
#'    True positive rate (or recall or sensitivity), which is \eqn{TP / (TP + TN)}.
#' }
#' \item{TNR}{
#'    True negative rate (or specificity or selectivity), which is \eqn{TN / (TN + FP)}.
#' }
#' \item{PPV}{
#'    Positive predictive value (or precision), which is \eqn{TP / (TP + FP)}.
#' }
#' \item{NPV}{
#'    Negative predictive value, which is \eqn{TN / (TN + FN)}.
#' }
#' \item{FPR}{
#'    False positive rate (type I error), which is \eqn{1 - TNR}.
#' }
#' \item{FNR}{
#'    False negative rate (type II error), which is \eqn{1 - TPR}.
#' }
#' \item{FDR}{
#'    False discovery rate, which is \eqn{1 - PPV}.
#' }
#' \item{FOR}{
#'    False omission rate, which is \eqn{1 - NPV}.
#' }
#' \item{PT}{
#'    Prevalence threshold, which is \eqn{(\sqrt{TPR \cdot FPR} - FPR) / (TPR - FPR)}.
#' }
#' \item{TS}{
#'    Threat score (or Jaccard index or critical success index (CSI)), which is \eqn{TP / (TP + FN + FP)}.
#' }
#' \item{FM}{
#'    Fowlkes-Mallows index, which is \eqn{\sqrt{PPV \cdot TPR}}.
#' }
#' \item{MK}{
#'    Markedness, which is \eqn{PPV + NPV - 1}.
#' }
#' \item{LRp}{
#'    Positive likelihood ratio, which is \eqn{TPR / FPR}.
#' }
#' \item{LRn}{
#'    Negative likelihood ratio, which is \eqn{FNR / TNR}.
#' }
#' @export
#' @references Sammut, C., & Webb, G. I. (2017). \emph{Encyclopedia of machine learning and data mining.}
#'  Springer Publishing Company, Incorporated.
#'  
#'  Powers, D. M. (2020). Evaluation: from precision, recall and F-measure to ROC, informedness,
#'  markedness and correlation. \emph{arXiv preprint arXiv:2010.16061.}
#'  
#'  Chicco, D., & Jurman, G. (2020). The advantages of the Matthews correlation coefficient (MCC)
#'  over F1 score and accuracy in binary classification evaluation. \emph{BMC genomics, 21,} 1-13.
#'  
#'  Balayla, J. (2020). Prevalence threshold (\eqn{\phi e}) and the geometry of screening curves.
#'  \emph{Plos one, 15}(10), e0240215.
#'
#' @examples
#' true_adj <- matrix(c(0, 0, 1, 0, 0,
#'                      0, 0, 1, 1, 0,
#'                      1, 1, 0, 0, 0,
#'                      0, 1, 0, 0, 1,
#'                      0, 0, 0, 1, 0), 5, 5)
#' est_adj <- matrix(c(0, 1, 1, 0, 0,
#'                     1, 0, 0, 1, 0,
#'                     1, 0, 0, 0, 1,
#'                     0, 1, 0, 0, 1,
#'                     0, 0, 1, 1, 0), 5, 5)
#' cm <- conf_matrix(true_adj, est_adj)
#' calculate_scores(cm)
calculate_scores <- function(cm) {
  tp <- cm[1,1]
  tn <- cm[2,2]
  fp <- cm[2,1]
  fn <- cm[1,2]
  tpr <- tp / (tp + fn)
  tnr <- tn / (tn + fp)
  ppv <- tp / (tp + fp)
  npv <- tn / (tn + fn)
  fnr <- 1 - tpr
  fpr <- 1 - tnr
  fdr <- 1 - ppv
  FOR <- 1 - npv
  lr_plus <- tpr / fpr
  lr_neg <- fnr / tnr
  pt <- (sqrt(tpr * fpr) - fpr) / (tpr - fpr)
  ts <- tp / (tp + fn + fp)
  fm <- sqrt(ppv * tpr)
  mk <- ppv + npv - 1
  acc <- (tp + tn) / (tp + tn + fn + fp)
  bal_acc <- (tpr + tnr) / 2
  F1_score <- 2 * (ppv * tpr) / (ppv + tpr)
  mcc <- (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  results <- data.frame(ACC = acc, ACC_bal = bal_acc, MCC = mcc, F1 = F1_score,
                        TPR = tpr, TNR = tnr, PPV = ppv, NPV = npv, FPR = fpr,
                        FNR = fnr, FDR = fdr, FOR = FOR, PT = pt, TS = ts, FM = fm, MK = mk,
                        LRp = lr_plus, LRn = lr_neg)
  return(results)
}