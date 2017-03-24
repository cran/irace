## ----include=FALSE-------------------------------------------------------
library(knitr)

## ----exampleload,eval=TRUE,include=FALSE------------------
library("irace")
load("examples.Rdata")
load("irace-output.Rdata")
options(width = 60)

## ----R_irace_install, prompt=FALSE, eval=FALSE------------
#  install.packages("irace")

## ----R_irace_launch,eval=FALSE, prompt=FALSE--------------
#  library("irace")
#  q() # To exit R

## ----install_win1,eval=FALSE, prompt=FALSE----------------
#  # Replace <package> with the path to the downloaded file.
#  install.packages("<package>", repos = NULL)

## ----install_win2,eval=FALSE, prompt=FALSE----------------
#  # Replace <package> with the path to the downloaded file.
#  # Replace <R_LIBS_USER> with the path used for installation.
#  install.packages("<package>", repos = NULL, lib = "<R_LIBS_USER>")
#  # Tell R where to find R_LIBS_USER.
#  # This must be executed for every new session.
#  .libPaths(c("<R_LIBS_USER>", .libPaths()))

## ----R_irace_test1, prompt=FALSE, eval=FALSE--------------
#  # Load the package
#  library("irace")
#  # Obtain the installation path
#  system.file(package = "irace")

## ----windows_irace_help ,eval=FALSE, prompt=FALSE---------
#  library("irace")
#  irace.cmdline("--help")

## ----irace_R_exe0, eval=FALSE, prompt=FALSE---------------
#  library("irace")
#  parameters <- readParameters("parameters.txt")
#  scenario <- readScenario(filename = "scenario.txt",
#                           scenario = defaultScenario())
#  irace(scenario = scenario, parameters = parameters)

## ----irace_R_check, eval=FALSE, prompt=FALSE--------------
#  library("irace")
#  parameters <- readParameters("parameters.txt")
#  scenario <- readScenario(filename = "scenario.txt",
#                           scenario = defaultScenario())
#  checkIraceScenario(scenario = scenario, parameters = parameters)

## ----irace_R_exe, eval=FALSE, prompt=FALSE----------------
#  library("irace")
#  # Go to the directory containing the scenario files
#  setwd("~/tuning")
#  # Create the R objects scenario and parameters
#  parameters <- readParameters("parameters.txt")
#  scenario <- readScenario(filename = "scenario.txt",
#                           scenario = defaultScenario())
#  irace(scenario = scenario, parameters = parameters)

## ----runexample2,prompt=FALSE,eval=FALSE------------------
#  library("irace")
#  setwd("~/tuning/")
#  parameters <- readParameters("parameters-acotsp.txt")
#  scenario <- readScenario(filename = "scenario.txt",
#                           scenario = defaultScenario())
#  irace(scenario = scenario, parameters = parameters)

## ----readParameters,prompt=FALSE, eval=FALSE--------------
#  parameters <- readParameters(file = "parameters.txt")

## ----parameterlist,eval=TRUE, prompt=TRUE, size='normalsize', comment=""----
print(parameters)

## ----targetRunner,prompt=FALSE, eval=FALSE----------------
#  targetRunner(experiment, scenario)

## ----experimentlist,eval=TRUE,size='normalsize', prompt=TRUE,  comment=""----
print(experiment)

## ----targetEvaluator, prompt=FALSE, eval=FALSE------------
#  targetEvaluator(experiment, num.configurations, all.conf.id,
#                  scenario, target.runner.call)

## ----instance1, prompt=FALSE, eval=FALSE------------------
#  scenario$instances <- c("rosenbrock_20", "rosenbrock_40",
#                          "rastrigin_20", "rastrigin_40")
#  scenario$instances.extra.params <-
#    c("--function=12 --nvar 20", "--function=12 --nvar 30",
#      "--function=15 --nvar 20", "--function=15 --nvar 30")

## ----repairEx,prompt=FALSE, eval=FALSE--------------------
#  repairConfiguration <- function (configuration, parameters, digits)
#  {
#    isreal <- parameters$type[colnames(configuration)] %in% "r"
#    configuration[isreal] <- configuration[isreal] / sum(configuration[isreal])
#    return(configuration)
#  }

## ----targetRunnerParallel,prompt=FALSE, eval=FALSE--------
#    targetRunnerParallel(experiments, exec.target.runner, scenario)

## ----targetRunnerParallel2,prompt=FALSE, eval=FALSE-------
#  targetRunnerParallel <- function(experiments, exec.target.runner, scenario)
#  {
#    return (lapply(experiments, exec.target.runner, scenario = scenario))
#  }

## ----targetRunnerParallel3,prompt=FALSE, eval=TRUE--------
   print(output)

## ----testing_r, prompt=FALSE, eval=FALSE------------------
#  testing.main(logFile = "./irace.Rdata")

## ----change_recover, prompt=TRUE, eval=FALSE--------------
#  load ("~/tuning/irace.Rdata")
#  new.path <- "~/experiments/tuning/instances/"
#  iraceResults$scenario$instances <-
#     paste (new.path,
#           basename(iraceResults$scenario$instances),
#           sep="")
#  save (iraceResults, file="~/tuning/irace.Rdata")

## ----load_rdata, prompt=FALSE, eval=FALSE-----------------
#  load("irace-output.Rdata")

## ----show_version, prompt=TRUE, eval=TRUE, comment=""-----
iraceResults$irace.version
iraceResults[["irace.version"]]

## ----show_configurations, prompt=TRUE, eval=TRUE, comment=""----
head(iraceResults$allConfigurations)

## ----show_idelites, prompt=TRUE, eval=TRUE, comment=""----
print(iraceResults$allElites)

## ----get_elites, prompt=TRUE, eval=TRUE, comment=""-------
getFinalElites(irace.logFile = "irace-output.Rdata", n = 0)

## ----show_iditelites, prompt=TRUE, eval=TRUE, comment=""----
print(iraceResults$iterationElites)

## ----get_elite, prompt=TRUE, eval=TRUE, comment=""--------
last <- length(iraceResults$iterationElites)
id <- iraceResults$iterationElites[last]
getConfigurationById(irace.logFile = "irace-output.Rdata", ids = id)

## ----get_experiments, prompt=TRUE, eval=TRUE, comment=""----
# As an example, we use the best configuration found
best.config <- getFinalElites(iraceResults = iraceResults, n = 1)
id <- best.config$.ID.
# Obtain the configurations using the identifier
# of the best configuration
all.exp <- iraceResults$experiments[,as.character(id)]
all.exp[!is.na(all.exp)]

## ----get_instance_seed, prompt=TRUE, eval=TRUE, comment=""----
# As an example, we get seed and instance of the experiments
# of the best candidate.
# Get index of the instances
pair.id <- names(all.exp[!is.na(all.exp)])
index <- iraceResults$state$.irace$instancesList[pair.id,"instance"]
# Obtain the instance names
iraceResults$scenario$instances[index]
# Get the seeds
iraceResults$state$.irace$instancesList[index,"seed"]

## ----get_model, prompt=TRUE, eval=TRUE, comment=""--------
# As an example, we get the model probabilities for the
# localsearch parameter.
iraceResults$state$model["localsearch"]
# The order of the probabilities corresponds to:
iraceResults$parameters$domain$localsearch

## ----get_test_exp, prompt=TRUE, eval=TRUE, comment=""-----
# Get the experiments of the testing
iraceResults$testing$experiments

## ----get_test_seeds, prompt=TRUE, eval=TRUE, comment=""----
# Get the experiments of the testing
iraceResults$testing$seeds

## ----plot_test, fig.pos="H", fig.align="center", out.width='0.75\\textwidth', fig.cap="Boxplot of the testing results of the best configurations.", prompt=TRUE, eval=TRUE, comment=""----
results <- iraceResults$testing$experiments
# Wilcoxon paired test
conf <- gl(ncol(results), # number of configurations
           nrow(results), # number of instances
           labels = colnames(results))
pairwise.wilcox.test (as.vector(results), conf, paired = TRUE, p.adj = "bonf")
# Plot the results
boxplot (iraceResults$testing$experiments,
         ylab = "Solution cost",
         xlab = "Configuration ID")

## ----freq, fig.pos="H", fig.cap="Parameters sampling frequency.", prompt=TRUE, eval=TRUE, comment=""----
parameterFrequency(iraceResults$allConfigurations, iraceResults$parameters)

## ----parcord, fig.pos="H", fig.align="center", out.width='0.75\\textwidth', fig.cap="Parallel coordinate plots of the parameters of the configurations in the last two iterations of a run of \\irace.", prompt=FALSE, eval=TRUE, comment=""----
# Get last iteration number
last <- length(iraceResults$iterationElites)
# Get configurations in the last two iterations
conf <- getConfigurationByIteration(iraceResults = iraceResults,
                                    iterations = c(last - 1, last))
parallelCoordinatesPlot (conf, iraceResults$parameters,
                         param_names = c("algorithm", "alpha",
                                         "beta", "rho", "q0"),
                         hierarchy = FALSE)

## ----faq3, eval=FALSE-------------------------------------
#  library(Rmpi)
#  mpi.spawn.Rslaves(nslaves = 10)
#  x <- mpi.applyLB(1:10, function(x) {
#    library(irace)
#    return(path.package("irace")) })
#         print(x)

## ----R_irace_home2, prompt=FALSE, eval=FALSE--------------
#  system.file(package = "irace")

