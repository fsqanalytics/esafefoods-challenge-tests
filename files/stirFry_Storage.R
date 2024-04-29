#' @title Storage of left overs of stir fried vegetables at home
#'
#' @description
#' The [stirFry_Storage()] function simulates the growth of \emph{L. monocytogenes} in left overs of vegetables that were stir fried. .
#' The function employs the lognormal distribution to represent the variability (at the portion level) in exponential growth rate at 5 \eqn{^\circ C\ log10\ CFU/g/h}
#' of \emph{L. monocytogenes} in heat-treated vegetables including corn, green peas, carrots, broccoli, beans and asparagus stored in air,
#' as determined in \insertCite{EFSA2020;textual}{JEMRA}. A log-linear growth model is employed, without any lag phase duration.
#'
#' @param data a list of:
#'   \describe{
#'      \item{`N`}{(`CFU`)  A matrix of size `nLots` lots by `sizeLot` units containing the numbers of \emph{L. monocytogenes} per
#'          portion of frozen vegetables to be defrosted;}
#'      \item{`ProbUnitPos`}{Lot-specific probability of contaminated portions or servings before defrosting (vector);}
#'      \item{`P`}{Mean prevalence of contaminated portions or servings (scalar).}
#'      }
#' @param nLots Number of lots sampled or size of the Monte Carlo simulation (scalar).
#' @param sizeLot Number of units or portions produced in a lot (scalar).
#' @param time (`h`) time that a portion of frozen vegetables is kept for defrosting (scalar or vector).
#' @param Temp (\eqn{^\circ C}) temperature at which frozen vegetables are defrosted (scalar or vector).
#' @param MPD (`log10 CFU/g`) is the  maximum population density of \emph{L. monocytogenes} in blanched vegetables (scalar or vector).
#' @param servingSize (`g`) weight of frozen vegetables portioned from the pack or serving size (scalar or vector).
#' @param pDefrost probability of defrosting frozen vegetables (scalar).
#' @param meanEGR5 (\eqn{log10\ CFU/g/h}) mean of the exponential growth rate of \emph{L. monocytogenes} in blanched vegetables at
#'  5 \eqn{^\circ C} (scalar or vector).
#' @param sdEGR5 (\eqn{log10\ CFU/g/h}) standard deviation of the exponential growth rate of \emph{L. monocytogenes} in
#'   blanched vegetables at 5 \eqn{^\circ C} (scalar or vector).
#' @param Tmin (\eqn{^\circ C}) minimum temperature for growth of \emph{L. monocytogenes} in blanched vegetables at
#'  5 \eqn{^\circ C} (scalar or vector).
#' @return A list of three elements:
#'     \describe{
#'        \item{`N`}{(`CFU`) A matrix of size `nLots` lots by `sizeLot` units representing the numbers of \emph{L. monocytogenes}
#'        in the portions of defrosted or non-defrosted vegetables.}
#'        \item{`ProbUnitPos`}{Lot-specific probability of contaminated portions or servings, defrosted or not (vector).}
#'        \item{`P`}{Mean prevalence of contaminated portions of defrosted or non-defrosted vegetables (scalar).}
#'        }
#'
#' @author Ursula Gonzales-Barron \email{ubarron@ipb.pt}
#'
#' @examples
#' Tmin <- -1.18
#' meanEGR5 <- 0.0117
#' sdEGR5 <- 0.00816
#' data <- LotGen(nLots=1000, sizeLot=1000, unitSize = 250,
#'               P=0.1, C0MeanLog=-1.5,C0SdLog=0.1)
#' res <- stirFry_Storage(data,
#'                        Temp = 5,
#'                        time = 16,
#'                        MPD = 8,
#'                        Tmin = -1.18,
#'                        meanEGR5 = meanEGR5,
#'                        sdEGR5 = sdEGR5, 
#'                        servingSize = 100,
#'                        pDefrost = 1)
#' hist(c(res$N))
#' hist(c(log10(res$N)), xlim=c(0,4))# histogram of microbial cells in contaminated defrosted/non-defrosted servings

stirFry_Storage <- function(data = list(),
                      nLots=NULL,
                      sizeLot=NULL,
                      Temp,
                      time,
                      MPD,
                      Tmin,
                      meanEGR5,
                      sdEGR5,
                      servingSize=NULL,
                      pDefrost) {
  if (pDefrost == 0) {
    return(data)
  } # Shortcut

  N_m <- data$N

  ifelse (exists("nLots", data)==TRUE, nLots <- data$nLots, nLots <- nrow(N_m))
  ifelse (exists("sizeLot", data)==TRUE, sizeLot <- data$sizeLot, sizeLot <- ncol(N_m))
  if (missing(servingSize)) servingSize <- data$servingSize # test if servingSize was defined
  if (is.null(servingSize)) warning("Add 'servingSize=#' to function arguments") # test again if servingSize is defined
  
#  nLots <- nrow(N_m)
#  sizeLot <- ncol(N_m)

  # Applied similarly to all lots to better compare lot impact
  portions_defrosted <- matrix(as.logical(stats::rbinom(n = sizeLot, size = 1, prob = pDefrost)),
    nrow = nLots, ncol = sizeLot, byrow = TRUE
  )

  lnorm_a <- log(meanEGR5^2 / sqrt(sdEGR5^2 + meanEGR5^2))
  lnorm_b <- sqrt(log(1 + (sdEGR5^2 / meanEGR5^2)))
  # Same EGR5 to better compare lot impact
  EGR5 <- matrix(stats::rlnorm(sizeLot, lnorm_a, lnorm_b),
    nrow = nLots, ncol = sizeLot, byrow = TRUE
  )

  # Evaluate on all and apply only on defrost
  delta <- (EGR5 * ((Temp - Tmin) / (5 - Tmin))^2) * time # log10 CFU/g
  Nf <- 10^(log10(N_m / servingSize) + delta) * servingSize
  Nf[N_m == 0] <- 0
  Nf <- pmin(Nf, servingSize * (10^MPD)) # checking for N values higher than MPD*servingSize

  data$N[portions_defrosted] <- round(Nf)[portions_defrosted]
   return(data)
}
