#' echoice2
#'
#' Choice Models with economic foundations
#'
#' @docType package
#' @author Nino Hardt
#' @import Rcpp RcppArmadillo
#' @importFrom Rcpp evalCpp
#' @useDynLib echoice2
#' @name echoice2
NULL


#' voldata_demo
#'
#' @name voldata_demo
#' @docType data
#' @details Simulated data mimicing an empirical dataset
NULL

#' ICdata
#'
#' @name ICdata
#' @docType data
#' @details Data from volumetric conjoint analysis in the ice cream category. 601 respondents total. 
#' Volumetric demand in units of 4 ounces each. Attributes include brand name, flavor, and container size. 
#' 8 brand intercepts, vector of 1s for baseline estimation.
NULL

#' ICdata
#'
#' @name ICdataUL
#' @docType data
#' @details Data from volumetric conjoint analysis in the ice cream category. 601 respondents total. 
#' Volumetric demand in units of 4 ounces each. Attributes include brand name, flavor, and container size. 
#' 8 brand intercepts, vector of 1s for baseline estimation.
#' Additionally, covariates are added, including 'kids in household' and certain 'Needs', such as look for a 'treat'.
#' 
#' }
NULL