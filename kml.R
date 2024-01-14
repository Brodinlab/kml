#!/usr/bin/env Rscript

# Set the library path
.libPaths("env")

# Set options
options(repos = c(CRAN = "https://cran.r-project.org/"),
	rgl.useNULL = TRUE
	)


# Function to check and install packages
check_and_install_packages <- function(package_names) {
  installed_packages <- rownames(installed.packages())
  
  for (pkg in package_names) {
    if (!pkg %in% installed_packages) {
      message("Installing package: ", pkg)
      install.packages(pkg, lib = "env")
    } else {
      message("Package already installed: ", pkg)
    }
  }
}

# Read the package names from the environment file
env_file_path <- "env/kml.env"
package_names <- readLines(env_file_path)

# Check and install packages
check_and_install_packages(package_names)

# Load the packages
lapply(package_names, library, character.only = TRUE)

#######################
# SETTING UP LOGGING #
######################
timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
logfile <- paste0("logfile_", timestamp, ".log")
lf = log_open(logfile)

options("logr.compact" = TRUE)

#######################
# SETTING UP ARGPARSER #
########################

sep("parsing cli args")

# Add the --input argument to the parser
parser <- ArgumentParser(description = "Command Line Argument Parsing in R") # parser
parser$add_argument("--input", type = "character", required = TRUE, 
                    help = "Input .rds with count matrices per taxa")
parser$add_argument("--metadata", type = "character", required = TRUE, 
                    help = "Input metadata path")
parser$add_argument("--nclusters", type = "integer", required = TRUE, 
                    help = "Number of trajectories to generate per taxa.")
parser$add_argument("--algorithm", type = "character", required = TRUE, 
                    help = "Clustering algorithm to use. If you're unsure use parWithEuclidean_rndm")
parser$add_argument("--subset", type = "character", required = TRUE, 
                    help = "Naming variable only. Is this data a particular subset?")
parser$add_argument("--transformation", type = "character", required = TRUE, 
                    help = "Naming variable only. Is this data transformed?")
parser$add_argument("--level", type = "character", required = TRUE, 
                    help = "Naming variable only. What taxonomic level does the data represent?")

# Parse the command line arguments
args <- parser$parse_args()

# Access the input argument
input_arg <- args$input
metadata_arg <- args$metadata
ncluster_arg <- args$nclusters
algorithm_arg <- args$algorithm
subset_arg <- args$subset
transformation_arg <- args$transformation
level_arg <- args$level

options("logr.notes" = FALSE)
put(paste("user provided input: ", input_arg, sep = ""))
put(paste("user provided metadata: ", metadata_arg, sep = ""))
put(paste("user provided ncluster: ", ncluster_arg, sep = ""))
put(paste("user provided algorithm: ", algorithm_arg, sep = ""))
put(paste("user provided subset: ", subset_arg, sep = ""))
put(paste("user provided transformation: ", transformation_arg, sep = ""))
put(paste("user provided level: ", level_arg, sep = ""))
options("logr.notes" = TRUE)

##########################
# Setting kml Parameters #
##########################

sep("setting upp kml parameters")

kml_parameters = list()
kml_parameters[["parWithMinkowski"]] = parALGO(distanceName="minkowski", startingCond="randomAll")
kml_parameters[["parWithMax"]] = parALGO(distanceName = "maximum", startingCond = "maxDist")
kml_parameters[["parWithEuclidean"]] = parALGO(distanceName = "euclidean", startingCond = "maxDist")
kml_parameters[["parWithEuclidean_rndm"]] = parALGO(distanceName = "euclidean", startingCond = "all")
kml_parameters[["parWithMedian"]] = parALGO(distanceName = "maximum", "startingCond" = "maxDist", centerMethod=function(x){median(x,na.rm=TRUE)})

put("kml parameters loaded")

#####################
# Loading functions #
#####################

sep("loading functions")

source("scripts/deconvoluteRobustness.R")
source("scripts/generateKmlData.R")
source("scripts/plotRobustness.R")
source("scripts/runRobustness.R")
source("scripts/statisticsRobustness.R")

put("functions loaded successfully")

######################
# RUNNING ROBUSTNESS #
######################

sep("running robustness")

input_rds = readRDS(input_arg)
metadata = read.csv(metadata_arg, row.names = "X")
rownames(metadata) = metadata[["subject"]]

kmlSeeds = run_robustness(seed_range = c(1:5), 
                          taxa_rds=input_rds,
                          nclusters=ncluster_arg,
                          parAlgorithm=kml_parameters[[algorithm_arg]],
                          robustnessdirname="output",
                          subset=subset_arg,
                          transformation=transformation_arg,
                          level=level_arg)

put("robustness finished")

#########################
# RUNNING DECONVOLUTION #
#########################

sep("running deconvolution")

kmlSeedsDeconvoluted = deconvolute_robustness(kmlSeeds=kmlSeeds, 
                                              subset=subset_arg, 
                                              transformation=transformation_arg)

put("deconvolution finished")

##########################
# SUMMARISING ROBUSTNESS #
##########################

sep("summarising robustness")

seedStatistics = robustness_statistics(kml_deconvoluted=kmlSeedsDeconvoluted, 
                                       metadata=metadata)

put("summarising finished")

log_close()
