#!/usr/bin/env Rscript


########################
# SETTING UP LIBRARIES #
########################

library("argparse")
library("dplyr")
library("phyloseq")
library("vegan")
library("logr")


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

parser <- ArgumentParser(description = "script to run a permanova for trajectory assignments") 
parser$add_argument("--input", type = "character", required = TRUE, 
                    help = "Path to CLR tx rds.")
parser$add_argument("--metadata", type = "character", required = TRUE, 
                    help = "Path metadata")
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
tax_arg <- args$tax
subset_arg <- args$subset
transformation_arg <- args$transformation
level_arg <- args$level

options("logr.notes" = FALSE)
put(paste("user provided input: ", input_arg, sep = ""))
put(paste("user provided metadata: ", metadata_arg, sep = ""))
put(paste("user provided taxtable: ", tax_arg, sep = ""))
put(paste("user provided subset: ", subset_arg, sep = ""))
put(paste("user provided transformation: ", transformation_arg, sep = ""))
put(paste("user provided level: ", level_arg, sep = ""))
options("logr.notes" = TRUE)

#############################
# SETTING UP BASE FUNCTIONS #
#############################

check_dirs <- function(subdir_path) {
    
        if (!dir.exists(subdir_path)) {
        dir.create(subdir_path, recursive = TRUE)
        put(paste("Directory", subdir_path, "created."))
    } else {
        put(paste("Directory", subdir_path, "already exists."))
    }
}


#######################
# SETTING UP PATHINGS #
#######################

sep("preparations to iterate adonis2")

output_path = paste("output", "/", subset_arg, "/", transformation_arg, "/", level_arg, sep = "")
check_dirs(output_path)

#######################
# SETTING UP METADATA #
#######################

path_to_subject_assignments = paste(output_path, "/kmlSeedStatistics.rds", sep = "")
subject_assignments = readRDS(path_to_subject_assignments)[["taxaSubjectAssignments"]]

meged_subjects_assignments <- subject_assignments %>%
    lapply(function(df) {dplyr::select(df, -count, -prop)}) %>% # Remove columns
    purrr::reduce(dplyr::full_join, by = "subject") %>% 
    distinct(subject, .keep_all = TRUE) # Merge on 'subject'


metadata = read.csv(metadata_arg, row.names = "X")
adonis2_md = merge(metadata, meged_subjects_assignments, by = "subject")
rownames(adonis2_md) = metadata$Sample


############################
# SETTING UP ADONIS2 INPUT #
############################

adonis2_input = readRDS(input_arg)

# running adonis2
sep("running adoni2")

adonis2_output = list()

for (level in c("species", "genus", "family", "order", "class")) {
    
    put(paste("working on level: ", level, sep = ""))
    
    adonis2_output[[level]] = list()
    
    for (taxa in names(subject_assignments)) {
        
        put(paste("working on taxa: ", taxa, sep = ""))
        
        adonis2_output[[level]][[taxa]] = list()
        
        taxa_columns = paste(taxa, "_clusters", sep = "")
        
        adonis2_taxa_md = adonis2_md %>% dplyr::rename(taxa_var = paste(taxa_columns))
        rownames(adonis2_taxa_md) = adonis2_taxa_md$Sample
        adonis2_taxa_md = adonis2_taxa_md[rownames(adonis2_input[[level]]), ]
        
        if (!identical(rownames(adonis2_taxa_md), rownames(adonis2_input[[level]]))) {
            
            put(taxa)
            put("DATAFRAMES ARE SORTED INACCURATELY")
            put("RESULTS ARE NOT REPRESENTATIVE")
            
        } else {
        
            adonis2_model = adonis2(adonis2_input[[level]] ~ taxa_var, data = adonis2_taxa_md, permutations = 999, method = "euclidean", by = "terms", strata = adonis2_taxa_md$subject)
            adonis2_model[["taxa"]] = taxa
            adonis2_model[["level"]] = level
            adonis2_model[["terms"]] = rownames(adonis2_model)
            adonis2_output[[level]][[taxa]] = adonis2_model
        
        }
    
    }
    
}

saveRDS(adonis2_output, paste(output_path, "/adonis2_output.rds", sep = ""))