#!/usr/bin/env Rscript

# Set the library path
.libPaths("env")

# Set options
options(repos = c(CRAN = "https://cran.r-project.org/"),
    rgl.useNULL = TRUE
    )

######################
# SET UP ENVIRONMENT # 
######################

# will install packages in env/

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

env_file_path <- "env/kml.env"
package_names <- readLines(env_file_path)

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

source("scripts/generateKmlData.R")
source("scripts/runRobustness.R")
source("scripts/deconvoluteRobustness.R")
source("scripts/plotRobustnessLogreg.R")
source("scripts/statisticsRobustnessLogreg.R")
#source("scripts/plotRobustness.R")
#source("scripts/statisticsRobustness.R")
#source("scripts/topSeedsRobustness.R")

put("functions loaded successfully")

######################
# RUNNING ROBUSTNESS #
######################

sep("running kml across seeds")

input_rds = readRDS(input_arg)
metadata = read.csv(metadata_arg, row.names = "X")
rownames(metadata) = metadata[["subject"]]

robustness_path = paste("output", "/", subset_arg, "/", transformation_arg, "/", level_arg, "/", sep = "")
plot_path = paste("output", "/", subset_arg, "/", transformation_arg, "/", level_arg, "/plots/", sep = "")
check_dirs(robustness_path)
check_dirs(plot_path)

kmlSeeds = run_robustness(seed_range = c(1:100), 
                          taxa_rds=input_rds,
                          nclusters=ncluster_arg,
                          parAlgorithm=kml_parameters[[algorithm_arg]],
                          robustnessdirname="output",
                          subset=subset_arg,
                          transformation=transformation_arg,
                          level=level_arg)

kmlSeeds_path = paste(robustness_path, "kmlSeeds.rds", sep = "")
saveRDS(kmlSeeds, kmlSeeds_path)

put("kml finished")

#########################
# RUNNING DECONVOLUTION #
#########################

sep("running deconvolution")

kmlSeedsDeconvoluted = deconvolute_robustness(kml_seeds=kmlSeeds, 
                                              subset=subset_arg, 
                                              transformation=transformation_arg,
                                              level=level_arg)

kmlSeedsDeconvoluted_path = paste(robustness_path, "kmlSeedsDeconvoluted.rds", sep = "")
saveRDS(kmlSeedsDeconvoluted, kmlSeedsDeconvoluted_path)

put("deconvolution finished")

###################################
# RUNNING STATISTICS ON ALL SEEDS #
###################################

sep("running per seed statistics")

kmlSeedStatistics = seed_statistics_logreg(kml_deconvoluted=kmlSeedsDeconvoluted, 
                                          metadata=metadata)

kmlSeedStatistics_path = paste(robustness_path, "kmlSeedStatistics.rds", sep = "")
saveRDS(kmlSeedStatistics, kmlSeedStatistics_path)


taxaSeedsSummed = do.call(rbind, kmlSeedStatistics[["taxaSeedsSummarised"]])

taxaSeedsSummed_path = paste(robustness_path, "taxaSeedsSummed.csv", sep = "")
write.csv(taxaSeedsSummed, taxaSeedsSummed_path)

put("per seed statistics finished")

####################################
# PLOTTING STATISTICS ON ALL SEEDS #
####################################


put("plotting volcano plots for all the seed scorings per taxa")

for (taxa in names(kmlSeeds)) {

    put(paste("plotting volcano: ", taxa, sep=""))

    a = plot_TaxaSeedStats(kmlSeedStatistics[["taxaSeedStats"]][[taxa]], taxa)
    
    volcano_plot_path = paste(plot_path, taxa, "_volcano.pdf", sep="")
    ggsave(volcano_plot_path, a, height = 10, width = 15)
}

put("plotting trajectories deconvoluted with log2 OR per taxa")
for (taxa in names(kmlSeeds)) {

    put(paste("plotting deconvoluted trajectories: ", taxa, sep=""))

    a = plotDeconvolutedMeanTraj_logreg(kmlSeedsDeconvoluted[[taxa]][["meltedMeanTraj"]], kmlSeedStatistics[["taxaSeedStats"]][[taxa]], taxa)
    
    deconvoluted_plot_path = paste(plot_path, taxa, "_deconvoluted.pdf", sep="")
    ggsave(deconvoluted_plot_path, a, height = 7, width = 24)
}

put("plotting volcano for final subject assignments")

taxaSeedsSummed_volcano = plot_TaxaSeedsSummed(taxaSeedsSummed)
taxaSeedsSummed_volcano_path = paste(robustness_path, "taxaSeedsSummed_volcano.pdf", sep = "")
ggsave(taxaSeedsSummed_volcano_path, taxaSeedsSummed_volcano, height = 10, width = 15)


put("plotting trajectories with final subject assignments per taxa with log2 OR")

for (taxa in names(kmlSeeds)) {

    put(paste("plotting final trajectories: ", taxa, sep=""))

    a = plot_subjectassignments_with_meanTraj(input_rds, kmlSeedStatistics, kmlSeedsDeconvoluted, taxa)
    
    trajectory_plot_path = paste(plot_path, taxa, "_FinalTrajectories.pdf", sep="")
    ggsave(trajectory_plot_path, a, height = 7, width = 24)
}


sep("RUN FINISHED CONGRATULATIONS")

log_close()
