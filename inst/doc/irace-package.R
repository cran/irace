## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
knit_hooks$set(inline = function(x) {
  if (is.numeric(x)) return(knitr:::format_sci(x, 'latex'))
  highr::hi_latex(x)
})

## ----include=FALSE------------------------------------------------------------
library(knitr)

## ----exampleload,eval=TRUE,include=FALSE----------------------------
library("irace")
load("examples.Rdata") # loads "experiment" and "output"
iraceResults <- irace::read_logfile(system.file(package="irace", "exdata", "irace-acotsp.Rdata", mustWork=TRUE))
log_ablation_file <- system.file(package="irace", "exdata", "log-ablation.Rdata", mustWork = TRUE)
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
#  irace_cmdline("--help")

## ----irace_R_check, eval=FALSE, prompt=FALSE------------------------
#  library("irace")
#  scenario <- readScenario(filename = "scenario.txt",
#                           scenario = defaultScenario())
#  checkIraceScenario(scenario = scenario)

## ----irace_R_exe, eval=FALSE, prompt=FALSE--------------------------
#  library("irace")
#  # Go to the directory containing the scenario files
#  setwd("./tuning")
#  scenario <- readScenario(filename = "scenario.txt",
#                           scenario = defaultScenario())
#  irace_main(scenario = scenario)

## ----runexample2,prompt=FALSE,eval=FALSE----------------------------
#  library("irace")
#  setwd("./tuning/")
#  irace_cmdline()

## ----readParameters,prompt=FALSE, eval=FALSE------------------------
#  parameters <- readParameters(file = "parameters.txt")

## ----parametersetup,eval=TRUE,include=FALSE-------------------------
# Setup example
parameters <- readParameters(text='
algorithm       "--"                 c		(as,mmas,eas,ras,acs)
ants            "--ants "            i	        (5, 100)
q0              "--q0 "              r          (0.0, 1.0)              | algorithm %in% c("acs")
')

## ----parameterlist,eval=TRUE, prompt=TRUE, size='normalsize', comment=""----
str(parameters, vec.len = 10)

## ----targetRunner,prompt=FALSE, eval=FALSE--------------------------
#  targetRunner(experiment, scenario)

## ----experimentlist,eval=TRUE,size='normalsize', prompt=TRUE,  comment=""----
print(experiment)

## ----targetEvaluator, eval=FALSE------------------------------------
#  targetEvaluator(experiment, num_configurations, all_conf_id, scenario,
#                  target_runner_call)

## ----instance1, prompt=FALSE, eval=FALSE----------------------------
#  scenario$instances <- c("rosenbrock_20 --function=12 --nvar 20",
#                          "rosenbrock_40 --function=12 --nvar 30",
#                          "rastrigin_20 --function=15 --nvar 20",
#                          "rastrigin_40 --function=15 --nvar 30")

## ----repairEx,prompt=FALSE, eval=FALSE------------------------------
#  repairConfiguration = function(configuration, parameters)
#  {
#    isreal <- names(which(parameters$types[colnames(configuration)] == "r"))
#    # This ignores 'digits'
#    c_real <- unlist(configuration[isreal])
#    c_real <- c_real / sum(c_real)
#    configuration[isreal] <- c_real
#    return(configuration)
#  }

## ----repairEx2,prompt=FALSE, eval=FALSE-----------------------------
#  repairConfiguration = function(configuration, parameters)
#  {
#   columns <- c("p1","p2","p3")
#   # cat("Before"); print(configuration)
#   configuration[columns] <- sort(unlist(configuration[columns], use.names=FALSE))
#   # cat("After"); print(configuration)
#   return(configuration)
#  }

## ----targetRunnerParallel,prompt=FALSE, eval=FALSE------------------
#  targetRunnerParallel(experiments, exec_target_runner, scenario, target_runner)

## ----targetRunnerParallel2,prompt=FALSE, eval=FALSE-----------------
#  targetRunnerParallel <- function(experiments, exec_target_runner, scenario,
#                                   target_runner)
#  {
#    lapply(experiments, exec_target_runner, scenario = scenario,
#           target_runner = target_runner)
#  }

## ----targetRunnerParallel3,prompt=FALSE, eval=TRUE------------------
   print(output)

## ----testing_r, prompt=FALSE, eval=FALSE----------------------------
#  testing_fromlog(logFile = "./irace.Rdata", testNbElites = 1)

## ----change_recover, prompt=TRUE, eval=FALSE------------------------
#  iraceResults <- read_logfile("./tuning/irace.Rdata")
#  new_path <- "./experiments/tuning/instances/"
#  iraceResults$scenario$instances <-
#     paste0(new_path, basename(iraceResults$scenario$instances))
#  save(iraceResults, file="./tuning/irace.Rdata")

## ----load_rdata, prompt=FALSE, eval=FALSE---------------------------
#  logfile <- system.file(package="irace", "exdata", "irace-acotsp.Rdata", mustWork=TRUE)
#  iraceResults <- read_logfile(logfile)

## ----show_version, prompt=TRUE, eval=TRUE, comment=""---------------
iraceResults$irace_version

## ----show_configurations, prompt=TRUE, eval=TRUE, comment=""--------
head(iraceResults$allConfigurations)

## ----show_idelites, prompt=TRUE, eval=TRUE, comment=""--------------
print(iraceResults$allElites)

## ----get_elites, prompt=TRUE, eval=TRUE, comment=""-----------------
logfile <- system.file(package="irace", "exdata", "irace-acotsp.Rdata", mustWork=TRUE)
getFinalElites(logfile, n = 0)

## ----show_iditelites, prompt=TRUE, eval=TRUE, comment=""------------
print(iraceResults$iterationElites)

## ----get_elite, prompt=TRUE, eval=TRUE, comment=""------------------
last <- length(iraceResults$iterationElites)
id <- iraceResults$iterationElites[last]
getConfigurationById(iraceResults, ids = id)

## ----get_experiments, prompt=TRUE, eval=TRUE, comment=""------------
# As an example, we use the best configuration found
best_config <- getFinalElites(iraceResults, n = 1)
best_id <- as.character(best_config$.ID.)
# Obtain the results of the best configuration
all_exp <- iraceResults$experiments[, best_id]
# all_exp is a vector and names(all_exp) is the (instance,seed) index.
all_exp
# Obtain the results of the first and best configurations
all_exp <- iraceResults$experiments[, c("1", best_id)]
# all_exp is a matrix: colnames(all_exp) is configurationID and
# rownames(all_exp) is the (instance,seed) index.
all_exp

## ----get_instance_seed, prompt=TRUE, eval=TRUE, comment=""----------
# As an example, we get instanceID, seeds and instances of the experiments
# of the best configuration.
# We could get the indexes of the instances on which at least one
# configuration was executed:
pair_index <- which(apply(!is.na(all_exp), 1L, any))
# or the instances on which all configurations were executed:
pair_index <- which(apply(!is.na(all_exp), 1L, all))
# but in this example we get the indexes of the instances executed for
# the best configuration.
pair_index <- which(!is.na(all_exp[, best_id]))
instanceID <- get_instanceID_seed_pairs(iraceResults)[["instanceID"]][pair_index]
# or get the seeds
get_instanceID_seed_pairs(iraceResults)[["seed"]][pair_index]
# or obtain the actual instances.
iraceResults$scenario$instances[instanceID]
# If the instances are of atomic type (integers, floating-point numbers or
# character strings), the above is similar to:
get_instanceID_seed_pairs(iraceResults, index = pair_index, instances=TRUE)

## ----get_model, prompt=TRUE, eval=TRUE, comment=""------------------
# As an example, we get the model probabilities for the
# localsearch parameter.
iraceResults$state$model["localsearch"]
# The order of the probabilities corresponds to:
iraceResults$scenario$parameters$domains$localsearch

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
# Get summary data from the logfile.
irs <- irace_summarise(iraceResults)
# Get number of iterations
iters <- irs$n_iterations
# Get number of experiments (runs of target-runner) up to each iteration
fes <- cumsum(table(iraceResults$state$experiment_log[["iteration"]]))
# Get the mean value of all experiments executed up to each iteration
# for the best configuration of that iteration.
elites <- as.character(iraceResults$iterationElites)
values <- colMeans(iraceResults$testing$experiments[, elites])
stderr <- function(x) sqrt(var(x)/length(x))
err <- apply(iraceResults$testing$experiments[, elites], 2L, stderr)
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
logfile <- system.file(package="irace", "exdata", "log-ablation.Rdata", mustWork=TRUE)
plotAblation(logfile)

## ----ablation_cmdline, prompt=FALSE, eval=TRUE,echo=FALSE, comment=""----
ablation_cmdline("--help")

## ----postsel, prompt=FALSE, eval=FALSE------------------------------
#  # Execute all elite configurations in the iterations
#  psRace("irace.Rdata", max_experiments = 0.5, elites=TRUE)
#  # Execute a set of configurations IDs providing budget
#  psRace("irace.Rdata", conf_ids = c(34, 87, 102, 172, 293), max_experiments = 500)

## ----targetCmdline, prompt=FALSE, eval=FALSE------------------------
#  targetRunner="./real_target_runner.py"
#  targetRunnerLauncher="python"
#  targetCmdLine="-m {targetRunner} {configurationID} {instanceID}\
#   --seed {seed} -i {instance} --cutoff {bound} {targetRunnerArgs}"

## ----faq3, eval=FALSE-----------------------------------------------
#  library(Rmpi)
#  mpi.spawn.Rslaves(nslaves = 10)
#  paths <- mpi.applyLB(1:10, function(x) {
#    library(irace); return(path.package("irace")) })
#  print(paths)

## ----R_irace_home2, prompt=FALSE, eval=FALSE------------------------
#  system.file(package = "irace")

