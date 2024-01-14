#!/usr/bin/env Rscript

########################################################
# Setting ROBUSTNESS RUNNING (KML ITERATION) FUNCTIONS #
########################################################

check_dirs <- function(subdir_path) {
    
        if (!dir.exists(subdir_path)) {
        dir.create(subdir_path, recursive = TRUE)
        message(paste("Directory", subdir_path, "created."))
    } else {
        message(paste("Directory", subdir_path, "already exists."))
    }
}

run_robustness <- function(seed_range = c(1:2), taxa_rds, nclusters, parAlgorithm, robustnessdirname, subset, transformation, level, subsetting = TRUE) {
    
    # nclusters = number of clusters to generate. Will do the same number of clusters for all the taxa. Maybe implement a way to change this. or run them separately.
    # parAlgorithm = the parameters for running kml, prespecified

    
    taxa_list = names(taxa_rds) # gets all taxa names
    kml_seeds <- list()

    for (taxa in taxa_list) {
        print(taxa)
        subdir_path = paste(robustnessdirname, "/", subset, "/", transformation, "/", level, "/", taxa, "/",  sep="")
        check_dirs(subdir_path)
        kml_seeds[[taxa]][["subdir_path"]] = subdir_path
    
        for (seed_val in seed_range) {
            print(seed_val)
            set.seed(seed_val)
        
            if (subsetting == TRUE) {
                
                taxa_wide_subset = slice_sample(taxa_rds[[taxa]], prop = 0.9)
                
            } else {taxa_wide_subset = taxa_rds[[taxa]]}
            
            
            taxa_kml = run_kml(taxa_wide_subset, nclusters, 1, parAlgorithm, taxa, subset, transformation)
            taxa_kml$subjectAssignment[["seed"]] = seed_val
            taxa_kml$taxaMeanTraj[["seed"]] = seed_val
            
            kml_seeds[[taxa]][[paste("seed: ", seed_val, sep = "")]] = taxa_kml
            
            taxa_kml_path = paste(subdir_path, "/", taxa, "_", seed_val, ".rds", sep = "")
            saveRDS(taxa_kml, taxa_kml_path)
    
            
        }
    }
    
    rds_path = paste(robustnessdirname, "/", subset, "/", transformation, "/", level, "/kml_seeds.rds", sep = "")
    
    saveRDS(kml_seeds, rds_path)
    
    return(kml_seeds)
}
