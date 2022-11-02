#' Mixture of 2 Gauss Distributions $4.2.1
#' @format A list of numeric vectors of Bimodal Gaussian Mixture Model.
#' \describe{
#'   \item{n200}{200 elements numeric vector}
#'   \item{n400}{400 elements numeric vector}
#'   \item{n600}{600 elements numeric vector}
#'   \item{n800}{800 elements numeric vector}
#'   \item{n1000}{1000 elements numeric vector}
#'   \item{n1200}{1200 elements numeric vector}
#' }
#' @source \code{\link{mix2gauss_gen}}
"mix2gauss"

#' Mixture of 3 Gaussian Distributions $4.2.2
#' @format A list of numeric vectors of Unimodal Gaussian Asymmetric Mixture
#'  Model.
#' \describe{
#'   \item{n200}{200 elements numeric vector}
#'   \item{n400}{400 elements numeric vector}
#'   \item{n600}{600 elements numeric vector}
#'   \item{n800}{800 elements numeric vector}
#'   \item{n1000}{1000 elements numeric vector}
#'   \item{n1200}{1200 elements numeric vector}
#' }
#' @source \code{\link{mix3gauss_gen}}
"mix3gauss"

#' Mixture of Exponential Distribution and Gamma Distribution
#' @format A list of numeric vectors of Mixture of Exponential and Gamma
#'  distribution Model.
#' \describe{
#'   \item{n200}{200 elements numeric vector}
#'   \item{n400}{400 elements numeric vector}
#'   \item{n600}{600 elements numeric vector}
#'   \item{n800}{800 elements numeric vector}
#'   \item{n1000}{1000 elements numeric vector}
#'   \item{n1200}{1200 elements numeric vector}
#' }
#' @source \code{\link{mixexpgamma_gen}}
"mixexpgamma"

#' Mixture of 2 Gauss Distributions $4.2.1, histogram version
#' @format A list of numeric vectors of Bimodal Gaussian Mixture Model.
#' \describe{
#'   \item{n200}{200 elements numeric vector}
#'   \item{n200p}{histogram sample data with bins 25}
#'   \item{n200f}{histogram frequency data with bins 25}
#'   \item{n400}{400 elements numeric vector}
#'   \item{n400p}{histogram sample data with bins 50}
#'   \item{n400f}{histogram frequency data with bins 50}
#'   \item{n800}{800 elements numeric vector}
#'   \item{n800p}{histogram sample data with bins 100}
#'   \item{n800f}{histogram frequency data with bins 100}
#' }
#' @source \code{\link{mix2gauss_gen}} \code{\link{databinning}}
"mix2gaussHist"


#' Mixture of Exponential Distribution and Gamma Distribution,
#' histogram version.
#' @format A list of numeric vectors of Mixture of Exponential and Gamma
#'  distribution Model.
#' \describe{
#'   \item{n200}{200 elements numeric vector}
#'   \item{n200p}{histogram sample data with bins 25}
#'   \item{n200f}{histogram frequency data with bins 25}
#'   \item{n400}{400 elements numeric vector}
#'   \item{n400p}{histogram sample data with bins 50}
#'   \item{n400f}{histogram frequency data with bins 50}
#'   \item{n800}{800 elements numeric vector}
#'   \item{n800p}{histogram sample data with bins 100}
#'   \item{n800f}{histogram frequency data with bins 100}
#' }
#' @source \code{\link{mix2gauss_gen}} \code{\link{databinning}}
"mixExpGammaHist"