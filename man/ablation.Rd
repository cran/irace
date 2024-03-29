% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ablation.R
\name{ablation}
\alias{ablation}
\title{Performs ablation between two configurations (from source to target).}
\usage{
ablation(
  iraceResults,
  src = 1L,
  target = NULL,
  ab.params = NULL,
  type = c("full", "racing"),
  n_instances = 1L,
  seed = 1234567,
  ablationLogFile = "log-ablation.Rdata",
  ...
)
}
\arguments{
\item{iraceResults}{(\code{list()}|\code{character(1)}) \code{iraceResults} object created by \pkg{irace} and typically saved in the log file \code{irace.Rdata}. If a character string is given, then it is interpreted as the path to the log file from which the \code{iraceResults} object will be loaded.}

\item{src, target}{Source and target configuration IDs. By default, the first configuration ever evaluated (ID 1) is used as \code{src} and the best configuration found by irace is used as target.}

\item{ab.params}{Specific parameter names to be used for the ablation. They must be in \code{parameters$names}. By default, use all parameters.}

\item{type}{Type of ablation to perform: \code{"full"} will execute each configuration on all \code{n_instances} to determine the best-performing one; \code{"racing"} will apply racing to find the best configurations.}

\item{n_instances}{(\code{integer(1)}) Number of instances used in \code{"full"} ablation will be \code{n_instances * scenario$firstTest}.}

\item{seed}{(\code{integer(1)}) Integer value to use as seed for the random number generation.}

\item{ablationLogFile}{(\code{character(1)}) Log file to save the ablation log. If \code{NULL}, the results are not saved to a file.}

\item{...}{Further arguments to override scenario settings, e.g., \code{debugLevel}, \code{parallel}, etc.}
}
\value{
A list containing the following elements:
\describe{
\item{configurations}{Configurations tested in the ablation.}
\item{instances}{A matrix with the instances used in the experiments. First column has the
instances IDs from \code{iraceResults$scenario$instances}, second column the seed assigned to the instance.}
\item{experiments}{A matrix with the results of the experiments (columns are configurations, rows are instances).}
\item{scenario}{Scenario object with the settings used for the experiments.}
\item{trajectory}{IDs of the best configurations at each step of the ablation.}
\item{best}{Best configuration found in the experiments.}
}
}
\description{
Ablation is a method for analyzing the differences between two configurations.
}
\examples{
\donttest{
logfile <- system.file(package="irace", "exdata", "sann.rda")
# Execute ablation between the first and the best configuration found by irace.
ablog <- ablation(logfile, ablationLogFile = NULL)
plotAblation(ablog)
# Execute ablation between two selected configurations, and selecting only a
# subset of parameters, directly reading the setup from the irace log file.
ablog <- ablation(logfile, src = 1, target = 10,
                  ab.params = c("temp"), ablationLogFile = NULL)
plotAblation(ablog)
}

}
\references{
C. Fawcett and H. H. Hoos. Analysing differences between algorithm
configurations through ablation. Journal of Heuristics, 22(4):431–458, 2016.
}
\seealso{
\code{\link[=plotAblation]{plotAblation()}}
}
\author{
Leslie Pérez Cáceres and Manuel López-Ibáñez
}
