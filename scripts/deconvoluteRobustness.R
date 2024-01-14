#!/usr/bin/env Rscript

#############################################################################
# Setting ROBUSTNESS DECONVOLUTION (REASSIGNING TRAJECTORY NAMES) FUNCTIONS #
#############################################################################

deconvolute_robustness <- function(kml_seeds, subset, transformation) {
    
    taxa_list = names(kml_seeds)
    kml_deconvoluted <- list()
    
    for (taxa in taxa_list) {
        
        print(taxa)
        cluster_column = paste(taxa, "_clusters", sep = "")
        
#        meanTraj_allSeeds <- do.call(rbind, lapply(robustness_genus[[taxa]], function(x) {x[["taxaMeanTraj"]]})) # rbind all meanTrajectories
        meanTraj_allSeeds <- do.call(rbind, lapply(robustness_genus[[taxa]], function(x) {
            if ("taxaMeanTraj" %in% names(x)) {
                return(x[["taxaMeanTraj"]])
            } else {invisible()}}))

        
        meanTraj_allSeeds$trajSeedIdentifier = paste(meanTraj_allSeeds[[cluster_column]], meanTraj_allSeeds$seed, sep = "_") # preserve original trajectory and seed label
        
#        subjectAssignments_allSeeds <- do.call(rbind, lapply(robustness_genus[[taxa]], function(x) {x[["subjectAssignment"]]})) # rbind all subjectAssignments
        
        subjectAssignments_allSeeds <- do.call(rbind, lapply(robustness_genus[[taxa]], function(x) {
            if ("subjectAssignment" %in% names(x)) {
                return(x[["subjectAssignment"]])
            } else {invisible()}}))
        
        
        subjectAssignments_allSeeds$trajSeedIdentifier = paste(subjectAssignments_allSeeds[[cluster_column]], subjectAssignments_allSeeds$seed, sep = "_") # preserve original trajectory and seed label
        
        ### prepping another wide table for running final kml to deconvolute labels
        
        meanTraj_allSeeds_wide = meanTraj_allSeeds %>% 
            select(trajSeedIdentifier, value, variable) %>% # variable is the timepoint ...
            pivot_wider(names_from = variable, values_from = value)
        
        nClusters = length(unique(meanTraj_allSeeds[[cluster_column]])) # number of clusters        
        meanTraj_allSeeds_wide_kml = run_kml(meanTraj_allSeeds_wide, nClusters, 1, parWithEuclidean_rnd, taxa, subset, transformation) # deconvolution kml (dictionary)
        
        ### relabeling individuals 
        subjectAssignments_allSeeds = subjectAssignments_allSeeds %>% select(subject, trajSeedIdentifier)
        subjectAssignments_deconvoluted = merge(subjectAssignments_allSeeds, meanTraj_allSeeds_wide_kml$subjectAssignment, by = "trajSeedIdentifier")
        

        ### for plotting the meanTrajs
        
        
        kml_deconvolution_table = merge(meanTraj_allSeeds_wide_kml$subjectAssignment, meanTraj_allSeeds_wide, by = "trajSeedIdentifier")
        kml_deconvolution_table_melted = reshape2::melt(kml_deconvolution_table, id.vars = c("trajSeedIdentifier", cluster_column))
        kml_deconvolution_table_melted$variable_numeric <- as.numeric(gsub("t", "", kml_deconvolution_table_melted$variable))
        kml_deconvolution_table_melted$original_label <- gsub("_.*", "", kml_deconvolution_table_melted$trajSeedIdentifier) # remove seed label and get original label
        
        kml_deconvoluted[[taxa]][["meltedMeanTraj"]] = kml_deconvolution_table_melted
        kml_deconvoluted[[taxa]][["wideMeanTraj"]] = meanTraj_allSeeds_wide
        kml_deconvoluted[[taxa]][["kmlDictionary"]] = meanTraj_allSeeds_wide_kml
        kml_deconvoluted[[taxa]][["subjectAssignmentsDeconvoluted"]] = subjectAssignments_deconvoluted

    }
    
    return(kml_deconvoluted)
}
