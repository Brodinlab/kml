#!/usr/bin/env Rscript

######################################################################
# Setting ROBUSTNESS STATISTICS (AGGREGATING ACROSS SEEDS) FUNCTIONS #
######################################################################

tile_seed_stats <- function(subjectAssignmentsDeconvoluted, taxa) {
    
    seed_stats = list()
    subjectAssignmentsDeconvoluted$seed <- gsub("*._", "", subjectAssignmentsDeconvoluted$trajSeedIdentifier)
    seeds = unique(subjectAssignmentsDeconvoluted$seed)
    
    cluster_column = paste(taxa, "_clusters", sep = "")
    clusters = unique(subjectAssignmentsDeconvoluted[[cluster_column]])
    
    for (seed_val in seeds) {
        
        seedDeconvoluted = subjectAssignmentsDeconvoluted %>% 
            subset(seed == seed_val & reactivity != "_nodata") %>% 
            group_by(.data[[cluster_column]], reactivity) %>% summarise(count = n()) %>% 
            pivot_wider(names_from = reactivity, values_from = count) %>% data.frame(check.names = FALSE)
        
        
        clusters = unique(seedDeconvoluted[[cluster_column]])
        
        for (cluster in clusters) {
            
            seedDeconvolutedCluster = seedDeconvoluted
            
            seedDeconvolutedCluster[[cluster_column]] <- ifelse(seedDeconvolutedCluster[[cluster_column]] == cluster, cluster, "X")
            seedDeconvolutedCluster = seedDeconvolutedCluster %>% 
                group_by(.data[[cluster_column]]) %>% 
                summarise(across(everything(), sum, na.rm = TRUE)) %>% 
                arrange(desc(.data[[cluster_column]])) %>%
                data.frame(check.names = FALSE)
            
            seedDeconvolutedCluster = seedDeconvolutedCluster[, c(cluster_column, "_NR", "_R")]
            rownames(seedDeconvolutedCluster) = seedDeconvolutedCluster[[cluster_column]]
            seedDeconvolutedCluster[[cluster_column]] <- NULL
            
            fisher_test = fisher.test(seedDeconvolutedCluster)
            
            seedDeconvolutedCluster[["trajectory"]] = rownames(seedDeconvolutedCluster)
            seedDeconvolutedCluster[["taxa"]] = taxa
            seedDeconvolutedCluster[["seed"]] = paste("seed: ", seed_val, sep = "")
            seedDeconvolutedCluster[["comparison"]] = paste(cluster, "v", "X", sep = "")
            seedDeconvolutedCluster[["fisher_pval"]] = fisher_test$p.value
            seedDeconvolutedCluster[["fisher_OR"]] = fisher_test$estimate
            seedDeconvolutedCluster[["fraction_reactive"]] = seedDeconvolutedCluster$`_R` / seedDeconvolutedCluster$`_NR`
            
            seed_stats[[paste(taxa, "_seed_", seed_val, "_", cluster, "v", "X", sep = "")]] = seedDeconvolutedCluster

        }
        
    }
    
    return(seed_stats)
}

robustness_statistics <- function(kml_deconvoluted, metadata) {
    # in metadata one entry per subject
    
    taxa_list = names(kml_deconvoluted)
    
    robustness_statistics = list()
    robustness_statistics[["taxaSeedStats"]] = list()
    robustness_statistics[["taxaSeedsMax"]] = list()
    robustness_statistics[["taxaSeedsMin"]] = list()
    robustness_statistics[["taxaSeedsSummarised"]] = list()    

    for (taxa in taxa_list) {
        
	put(taxa)
        taxa_deconvoluted = merge(kml_deconvoluted[[taxa]][["subjectAssignmentsDeconvoluted"]], metadata, by = "subject")
        robustness_statistics[["taxaSeedStats"]][[taxa]] = do.call(rbind, tile_seed_stats(taxa_deconvoluted, taxa))
	put("tiling completed, now summarising")
        
        # for each trajectory, returns a dataframe with mean OR across all seeds
        avg = robustness_statistics[["taxaSeedStats"]][[taxa]] %>% 
                        subset(trajectory != 'X') %>% 
                        group_by(trajectory) %>% 
                        summarise(OR_mean = mean(fisher_OR), mean_NR = mean(.data[["_NR"]]), mean_R = mean(.data[["_R"]]))
        
        # for each seed, returns a dataframe with the trajectory that had the highest OR
        robustness_statistics[["taxaSeedsMax"]][[taxa]] = robustness_statistics[["taxaSeedStats"]][[taxa]] %>% 
                        subset(trajectory != 'X') %>% 
                        arrange(seed, desc(fisher_OR)) %>% 
                        distinct(.data[["seed"]], .keep_all = TRUE)
        
        # for each seed, returns a dataframe with the trajectory that had the lowest OR
        robustness_statistics[["taxaSeedsMin"]][[taxa]] = robustness_statistics[["taxaSeedStats"]][[taxa]] %>% 
                        subset(trajectory != 'X') %>% 
                        arrange(seed, fisher_OR) %>% 
                        distinct(.data[["seed"]], .keep_all = TRUE)
            
	uniqueTrajectories = unique(avg[["trajectory"]]) # remembers the unique trajectories for this taxa
        
        ### completes a chisq test comparing max OR cluster distributions vs expected distribution (i.e. numberUniqueSeeds / numberUniqueClusters) which is random
        max_grouped = robustness_statistics[["taxaSeedsMax"]][[taxa]] %>% group_by(trajectory) %>% summarise(count = n()) %>% data.frame(check.names = FALSE)
	max_grouped = max_grouped %>% complete(trajectory = uniqueTrajectories, fill = list("count" = 0)) # this code autocompletes the 2x2 matrix and adds missing trajectories with a count of 0
        max_grouped$p.val = chisq.test(max_grouped$count)$p.val
        max_grouped$comparison = 'max'
        
        ### completes a chisq test comparing min OR cluster distributions vs expected distribution (i.e. numberUniqueSeeds / numberUniqueClusters) which is random
        min_grouped = robustness_statistics[["taxaSeedsMin"]][[taxa]] %>% group_by(trajectory) %>% summarise(count = n()) %>% data.frame(check.names = FALSE)
	min_grouped = min_grouped %>% complete(trajectory = uniqueTrajectories, fill = list("count" = 0)) # this code autocompletes the 2x2 matrix and adds missing trajectories with a count of 0
        min_grouped$p.val = chisq.test(min_grouped$count)$p.val
        min_grouped$comparison = 'min'

        ### aggregates data with mean OR across all seeds, not just min or max
        min_max_avg = merge(rbind(max_grouped, min_grouped), avg, by = "trajectory")
        min_max_avg$filter <- ifelse(min_max_avg$comparison == "max" & min_max_avg$OR_mean > 1, "keep",
                            ifelse(min_max_avg$comparison == "min" & min_max_avg$OR_mean < 1, "keep", "filter"))
        
        
        robustness_statistics[["taxaSeedsSummarised"]][[taxa]] = min_max_avg
            
    }
    
    return(robustness_statistics)
}
