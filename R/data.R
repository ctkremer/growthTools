#' Microbial abundance time series
#'
#' A data set tracking the abundances of experimental phytoplankton populations over time. Fluorescence 
#' is used as a proxy for numerical abundance/biomass.
#'
#' @format A data frame with 46046 rows and 19 variables:
#' \describe{
#'   \item{isolate.id}{name of a given phytoplankton isolate}
#'   \item{experiment.ID}{unique identifier combining isolate, temperature, and dilution number}
#'   \item{dilution}{dilution number, ranging from 1-3. Populations were periodically diluted to provide new nutrients for growth}
#'   \item{replicate}{experimental replicate (4 replicates per temperature/isolate)}
#'   \item{plate.file}{name of data file produced by plate reader, containing fluorescence readings}
#'   \item{date}{date of measurement, yyyy-mm-dd format}
#'   \item{time}{time of day of measurement, hh:mm}
#'   \item{media}{growth media use for culturing populations}
#'   \item{well}{location of specific population in the 24-well plate}
#'   \item{fluorescence}{measured fluorescence of population, a standard proxy for abundance}
#'   \item{temperature}{temperature, in Celsius}
#'   \item{dtime}{time since first measure of a population's abundance, in hours}
#'   \item{log.fluor}{natural log of fluorescence}
#'   \item{genus of phytoplankton population}
#'   \item{species}{species epithet of phytoplankton population}
#'   \item{group}{functional group phytoplankton species belongs to}
#'   \item{source}{location where phytoplankton species were isolate (NB = Narragansett Bay)}
#'   \item{plate.temperature}{temperature of well plate during fluorescence readings(?)}
#'   \item{id}{unique identifier for single experimental phytoplankton populations (experimental units)}
#' }
#' @source Elena Litchman/Tatiana Severin, W.K. Kellogg Biological Station, Michigan State University
"abundance"

#' Effect of temperature on phytoplankton growth rates
#'
#' A data set containing estimated growth rates of two phytoplankton isolates across a range
#' of temperatures. Growth rates were estimated from abundance time series (abundance data set) 
#' using methods in the growthTools package.
#'
#' @format A data frame with 160 rows and 8 variables:
#' \describe{
#'   \item{experiment.ID}{unique identifier combining isolate, temperature, and dilution number}
#'   \item{isolate.id}{name of a given phytoplankton isolate}
#'   \item{temperature}{temperature, in Celsius}
#'   \item{dilution}{dilution number, ranging from 1-3}
#'   \item{replicate}{experimental replicate (4 replicates per temperature/isolate)}
#'   \item{mu}{estimated growth rate}
#'   \item{best.model}{time series model generating estimated growth rate}
#'   \item{id}{unique identifier combining isolate id and dilution only}
#' }
#' @source Elena Litchman/Tatiana Severin, W.K. Kellogg Biological Station, Michigan State University
"example_TPC_data"



