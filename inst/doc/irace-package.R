## ----include=FALSE------------------------------------------------------------
library(knitr)

## ----exampleload,eval=TRUE,include=FALSE----------------------------
library("irace")
load("examples.Rdata")
iraceResults <- irace::read_logfile("irace-acotsp.Rdata")
log_ablation_file <- file.path(system.file(package="irace"), "exdata", "log-ablation.Rdata")
load(log_ablation_file)
options(width = 70)

## ----R_irace_install, prompt=FALSE, eval=FALSE----------------------
#  install.packages("irace")

## ----R_irace_launch,eval=FALSE, prompt=FALSE------------------------
#  library("irace")
#  q() # To exit R

## ----install_win1,eval=FALSE, prompt=FALSE--------------------------
#  install.packages("<package>", repos = NULL)

## ----install_win2,eval=FALSE, prompt=FALSE--------------------------
#  # Replace <package> with the path to the downloaded file.
#  # Replace <R_LIBS_USER> with the path used for installation.
#  install.packages("<package>", repos = NULL, lib = "<R_LIBS_USER>")
#  # Tell R where to find R_LIBS_USER.
#  # This must be executed for every new session.
#  .libPaths(c("<R_LIBS_USER>", .libPaths()))

## ----R_irace_test1, prompt=FALSE, eval=FALSE------------------------
#  # Load the package
#  library("irace")
#  # Obtain the installation path
#  system.file(package = "irace")

## ----windows_irace_help ,eval=FALSE, prompt=FALSE-------------------
#  library("irace")
#  irace.cmdline("--help")

## ----irace_R_check, eval=FALSE, prompt=FALSE------------------------
#  library("irace")
#  scenario <- readScenario(filename = "scenario.txt",
#                           scenario = defaultScenario())
#  checkIraceScenario(scenario = scenario)

## ----irace_R_exe, eval=FALSE, prompt=FALSE--------------------------
#  library("irace")
#  # Go to the directory containing the scenario files
#  setwd("~/tuning")
#  scenario <- readScenario(filename = "scenario.txt",
#                           scenario = defaultScenario())
#  irace.main(scenario = scenario)

## ----runexample2,prompt=FALSE,eval=FALSE----------------------------
#  library("irace")
#  setwd("~/tuning/")
#  irace.cmdline()

## ----readParameters,prompt=FALSE, eval=FALSE------------------------
#  parameters <- readParameters(file = "parameters.txt")

## ----parametersetup,eval=TRUE,include=FALSE-------------------------
# Setup example
parameters <- irace::readParameters(text='
algorithm       "--"                 c		(as,mmas,eas,ras,acs)
ants            "--ants "            i	  	(5, 100)
q0              "--q0 "              r  	(0.0, 1.0) 		| algorithm %in% c("acs")
')

## ----parameterlist,eval=TRUE, prompt=TRUE, size='normalsize', comment=""----
str(parameters, vec.len = 10)

## ----targetRunner,prompt=FALSE, eval=FALSE--------------------------
#  targetRunner(experiment, scenario)

## ----experimentlist,eval=TRUE,size='normalsize', prompt=TRUE,  comment=""----
print(experiment)

## ----targetEvaluator, prompt=FALSE, eval=FALSE----------------------
#  targetEvaluator(experiment, num.configurations, all.conf.id,
#                  scenario, target.runner.call)

## ----instance1, prompt=FALSE, eval=FALSE----------------------------
#  scenario$instances <- c("rosenbrock_20 --function=12 --nvar 20",
#                          "rosenbrock_40 --function=12 --nvar 30",
#                          "rastrigin_20 --function=15 --nvar 20",
#                          "rastrigin_40 --function=15 --nvar 30")

## ----repairEx,prompt=FALSE, eval=FALSE------------------------------
#  repairConfiguration = function (configuration, parameters, digits)
#  {
#    isreal <- parameters$type[colnames(configuration)] %in% "r"
#    configuration[isreal] <- configuration[isreal] / sum(configuration[isreal])
#    configurations[isreal] <- round(configuration[isreal], digits)
#    return(configuration)
#  }

## ----repairEx2,prompt=FALSE, eval=FALSE-----------------------------
#  repairConfiguration = function (configuration, parameters, digits)
#  {
#   columns <- c("p1","p2","p3")
#   # cat("Before"); print(configuration)
#   configuration[columns] <- sort(configuration[columns])
#   # cat("After"); print(configuration)
#   return(configuration)
#  }

## ----targetRunnerParallel,prompt=FALSE, eval=FALSE------------------
#  targetRunnerParallel(experiments, exec.target.runner, scenario, target.runner)

## ----targetRunnerParallel2,prompt=FALSE, eval=FALSE-----------------
#  targetRunnerParallel <- function(experiments, exec.target.runner, scenario)
#  {
#    return (lapply(experiments, exec.target.runner, scenario = scenario,
#                   target.runner = target.runner))
#  }

## ----targetRunnerParallel3,prompt=FALSE, eval=TRUE------------------
   print(output)

## ----testing_r, prompt=FALSE, eval=FALSE----------------------------
#  testing_fromlog(logFile = "./irace.Rdata", testNbElites = 1)

## ----change_recover, prompt=TRUE, eval=FALSE------------------------
#  load ("~/tuning/irace.Rdata")
#  new.path <- "~/experiments/tuning/instances/"
#  iraceResults$scenario$instances <-
#     paste (new.path,
#           basename(iraceResults$scenario$instances),
#           sep="")
#  save (iraceResults, file="~/tuning/irace.Rdata")

## ----load_rdata, prompt=FALSE, eval=FALSE---------------------------
#  load("irace-acotsp.Rdata")

## ----show_version, prompt=TRUE, eval=TRUE, comment=""---------------
iraceResults$irace.version
iraceResults[["irace.version"]]

## ----show_configurations, prompt=TRUE, eval=TRUE, comment=""--------
head(iraceResults$allConfigurations)

## ----show_idelites, prompt=TRUE, eval=TRUE, comment=""--------------
print(iraceResults$allElites)

## ----get_elites, prompt=TRUE, eval=TRUE, comment=""-----------------
getFinalElites("irace-acotsp.Rdata", n = 0)

## ----show_iditelites, prompt=TRUE, eval=TRUE, comment=""------------
print(iraceResults$iterationElites)

## ----get_elite, prompt=TRUE, eval=TRUE, comment=""------------------
last <- length(iraceResults$iterationElites)
id <- iraceResults$iterationElites[last]
getConfigurationById(iraceResults, ids = id)

## ----get_experiments, prompt=TRUE, eval=TRUE, comment=""------------
# As an example, we use the best configuration found
best.config <- getFinalElites(iraceResults = iraceResults, n = 1)
id <- best.config$.ID.
# Obtain the configurations using the identifier
# of the best configuration
all.exp <- iraceResults$experiments[,as.character(id)]
all.exp[!is.na(all.exp)]

## ----get_instance_seed, prompt=TRUE, eval=TRUE, comment=""----------
# As an example, we get seed and instance of the experiments
# of the best candidate.
# Get index of the instances
pair.id <- which(!is.na(all.exp))
index <- iraceResults$state$.irace$instancesList[pair.id, "instance"]
# Obtain the instance names
iraceResults$scenario$instances[index]
# Get the seeds
iraceResults$state$.irace$instancesList[pair.id, "seed"]

## ----get_model, prompt=TRUE, eval=TRUE, comment=""------------------
# As an example, we get the model probabilities for the
# localsearch parameter.
iraceResults$state$model["localsearch"]
# The order of the probabilities corresponds to:
iraceResults$parameters$domain$localsearch

## ----get_test_exp, prompt=TRUE, eval=TRUE, comment=""---------------
# Get the results of the testing
iraceResults$testing$experiments

## ----get_test_seeds, prompt=TRUE, eval=TRUE, comment=""-------------
# Get the seeds used for testing
iraceResults$testing$seeds

## ----wilcox_test,prompt=TRUE, eval=TRUE, comment=""-----------------
results <- iraceResults$testing$experiments
# Wilcoxon paired test
conf <- gl(ncol(results), # number of configurations
           nrow(results), # number of instances
           labels = colnames(results))
pairwise.wilcox.test (as.vector(results), conf, paired = TRUE, p.adj = "bonf")

## ----conc, prompt=TRUE, eval=TRUE, comment=""-----------------------
irace:::concordance(iraceResults$testing$experiments)

## ----testEvo, fig.pos="tb", fig.align="center", out.width='0.7\\textwidth', fig.cap="Testing set performance of the best-so-far configuration over number of experiments. Label of each point is the configuration ID.", prompt=FALSE, eval=TRUE, comment=""----
# Get number of iterations
iters <- unique(iraceResults$experimentLog[, "iteration"])
# Get number of experiments (runs of target-runner) up to each iteration
fes <- cumsum(table(iraceResults$experimentLog[,"iteration"]))
# Get the mean value of all experiments executed up to each iteration
# for the best configuration of that iteration.
elites <- as.character(iraceResults$iterationElites)
values <- colMeans(iraceResults$testing$experiments[, elites])
stderr <- function(x) sqrt(var(x)/length(x))
err <- apply(iraceResults$testing$experiments[, elites], 2, stderr)
plot(fes, values, type = "s",
     xlab = "Number of runs of the target algorithm",
     ylab = "Mean value over testing set", ylim=c(23000000,23500000))
points(fes, values, pch=19)
arrows(fes, values - err, fes, values + err, length=0.05, angle=90, code=3)
text(fes, values, elites, pos = 1)

## ----ablation, prompt=FALSE, eval=FALSE-----------------------------
#  ablog <- ablation("irace.Rdata", src = 1, target = 60)
#  plotAblation(ablog)

## ----testAb, fig.pos="htb!", fig.align="center", out.width="0.75\\textwidth", fig.cap="Example of plot generated by \\code{plotAblation()}.", prompt=FALSE, eval=TRUE, echo=FALSE----
logfile <- file.path(system.file(package="irace"), "exdata", "log-ablation.Rdata")
plotAblation(logfile)

## ----ablation_cmdline, prompt=FALSE, eval=TRUE,echo=FALSE-----------
ablation_cmdline("--help")

## ----postsel, prompt=FALSE, eval=FALSE------------------------------
#  # Execute all elite configurations in the iterations
#  psRace(iraceLogFile="irace.Rdata", elites=TRUE)
#  # Execute a set of configurations IDs providing budget
#  psRace(iraceLogFile="irace.Rdata",
#               conf.ids=c(34, 87, 102, 172, 293),
#               max.experiments=500)

## ----faq3, eval=FALSE-----------------------------------------------
#  library(Rmpi)
#  mpi.spawn.Rslaves(nslaves = 10)
#  paths <- mpi.applyLB(1:10, function(x) {
#    library(irace); return(path.package("irace")) })
#  print(paths)

## ----R_irace_home2, prompt=FALSE, eval=FALSE------------------------
#  system.file(package = "irace")

