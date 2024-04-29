#' @title  Effect of stir-frying on \emph{L. monocytogenes} in frozen vegetables
#'
#' @description
#' The function [stirFry_MildHeat()] describes the effect of stir-frying on the numbers of \emph{L. monocytogenes} 
#' present in vegetables such as broccoli, mushroom, onions, peas and pepper. The function is based on the Bigelow
#' model, which describes the decimal reduction time (`D`) as a function of temperature, with parameters `z` and reference `D` (`Dref`) 
#' at 70 \eqn{^\circ C}.
#' This function assumes that the extent of reduction is global for the aforementioned vegetables, and depends on the duration and
#' temperature of stir-frying; Whereas the algorithm considers the `z` value as fixed.
#' 
#' @author Ursula Gonzales-Barron \email{ubarron@ipb.pt}
#' 
#' @examples
#'                   
#' data <- LotGen(nLots=1000, sizeLot=1000, unitSize = 250,
#'               P=0.1, C0MeanLog=0.8,C0SdLog=0.3)
#' res <- stirFry_MildHeat(data,
#'                     tempBlanch = 68,
#'                     timeBlanch = 1.5,
#'                     logDrefMean=-1.78,
#'                     logDrefSd=0.251,
#'                     zT=6.06)
#' res$P
#' hist(c(res$N),xlim=c(0,5000))
#' hist(c(log10(res$N)), xlim=c(0,4))
#'
#'
stirFry_MildHeat <- function(data=list(),
                        nLots=NULL,
                        tempBlanch,
                        timeBlanch,
                        logDrefMean=-1.78,
                        logDrefSd=0.251,
                        zT=6.06){

  ifelse (exists("nLots", data)==TRUE, nLots <- data$nLots, nLots <- nrow(data$N))
#  if (exists("nLots", data)==TRUE) {
#    nLots <- data$nLots
#    } else {
#     print("Add 'nLots=#' to function arguments")
#   }
#  nLots <- nrow(data$N)
#  sizeLot <- ncol(data$N)
  N0 <- data$N
  
  #Draw heat inactivation values at batch level
  log10Dref <- stats::rnorm(nLots,
                            logDrefMean,
                            logDrefSd)
  D <- timeBlanch*(log10Dref-((tempBlanch-70)/zT))
  pSurvive <- 10^D
  
  # Trick to get only positive lots
  # apply the log-reduction to the sum and redistribute, proportionally to the original distribution
  # Watch out: overall equivalent, but not serving to serving (may lead to an increase for a given portion)
  # but not a problem here
  sumN <- rowSums(data$N)
  sumN1 <- extraDistr::rtbinom(n=nLots, size=sumN, prob=pSurvive, a=0)
  # To avoid special cases with extremely low pSurvive
  # (check extraDistr::rtbinom(n= 10, size=10, prob = 1E-18, a=0))
  sumN1[sumN1 == sumN & pSurvive < 1E-5] <- 1
  
  N1 <- mc2d::rmultinomial(n=nLots, size = sumN1, prob = data$N)
  
  # Evaluate the probabilities of at least one survive
  atLeastOneSurvive <- 1-(1-pSurvive)^sumN
  
  ProbUnitPos <- data$ProbUnitPos * atLeastOneSurvive 

  P1 <- data$P * mean(atLeastOneSurvive)
  
  data$nLots <- nLots
  data$N <- N1
  data$P <- P1
  data$ProbUnitPos <- ProbUnitPos
  data$pSurviveBlanching <- pSurvive
  return(data)
}
