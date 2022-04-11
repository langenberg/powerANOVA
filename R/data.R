#' @title Example data from Janczyk and Leuthold (2017).
#'
#'@references{
#'   \insertRef{Janczyk2017a}{powerANOVA}
#' }
#'
#' @description
#' The \code{simon} data set originates from a Simon experiment conducted by Janczyk and Leuthold (2017, second experiment). N = 36 subjects participated in the experiment. Three repeated measures factors were manipulated: (1) congruency in task n, (2) congruency in task n-1, and (3) effector switch. See the original study for more information.
#'
#' @format
#' \describe{
#'     \item{id}{Subject ID.}
#'     \item{incongruent_incongruent_repitition}{The first part of the variable name refers to congruency in trial n. The second part of the variable name refers to congruency in trial n-1. The last part of the variable indicates whether the effector system was repeated or switched.}
#'     \item{...}{More variables in the above format.}
#' }
#'
#'@importFrom Rdpack reprompt
#'
"simon"
