#!/usr/bin/env Rscript

######################################################################
# Setting TOP SEEDS AGGREGATION (AGGREGATING ACROSS SEEDS) FUNCTIONS #
######################################################################

robustness_top_seeds <- function(kmlSeedStatistics) {
    
    taxaTrajTopSeeds = list()
    taxa_list = names(kmlSeedStatistics[["taxaSeedsSummarised"]])
    
    for (taxa in taxa_list) {
        
        taxaTrajTopSeeds[[taxa]] = list()
        
        taxa_summarised_max = kmlSeedStatistics[["taxaSeedsSummarised"]][[taxa]] %>% subset(filter == 'keep' & comparison == 'max')
        taxa_summarised_min = kmlSeedStatistics[["taxaSeedsSummarised"]][[taxa]] %>% subset(filter == 'keep' & comparison == 'min')
        
        max_trajectories = unique(taxa_summarised_max[["trajectory"]])
        min_trajectories = unique(taxa_summarised_min[["trajectory"]])
        
        for (traj in max_trajectories) {
            
            taxaTrajTopSeeds[[taxa]][[traj]] = list()
            taxaSeedStatistics = kmlSeedStatistics[["taxaSeedStats"]][[taxa]] %>% subset(trajectory == traj) %>% arrange(desc(fisher_OR)) %>% head(20) # grab 20 best seeds
            taxaTrajTopSeeds[[taxa]][[traj]] = taxaSeedStatistics[, c("trajectory", "seed", "fisher_OR", "fisher_pval")]
            
        }
        
        for (traj in min_trajectories) {
            
            taxaTrajTopSeeds[[taxa]][[traj]] = list()
            taxaSeedStatistics = kmlSeedStatistics[["taxaSeedStats"]][[taxa]] %>% subset(trajectory == traj) %>% arrange(fisher_OR) %>% head(20) # grab 20 best seeds
            taxaTrajTopSeeds[[taxa]][[traj]] = taxaSeedStatistics[, c("trajectory", "seed", "fisher_OR", "fisher_pval")]
            
        }
        
    }
    
    return(taxaTrajTopSeeds)
}

summarise_top_seeds <- function(kmlSeedStatistics, kmlSeedsDeconvoluted, metadata, threshold) {
    # give threshold like a fraction e.g. .75%. needs to be > 0.5 at least.
    
    topSeedsSummarised = list()
    topSeedsSummarised[["topSubjectAssginments"]] = list()
    topSeedsSummarised[["topSeedsSummarised"]] = list()
    
    topSeedsSummarised[["taxaTrajTopSeeds"]] = robustness_top_seeds(kmlSeedStatistics)
    taxa_list = names(topSeedsSummarised[["taxaTrajTopSeeds"]])
    
    for (taxa in taxa_list) {
        
        topSeedsSummarised[["topSubjectAssginments"]][[taxa]] = list()
        topSeedsSummarised[["topSeedsSummarised"]][[taxa]] = list()
        
        cluster_column = paste(taxa, "clusters", sep="_")
        traj_list = names(topSeedsSummarised[["taxaTrajTopSeeds"]][[taxa]])
        
        for (traj in traj_list) {
            
            topSeedsSummarised[["topSubjectAssginments"]][[taxa]][[traj]] = list()
            topSeedsSummarised[["topSeedsSummarised"]][[taxa]][[traj]] = list()
            
            topSeeds = unique(topSeedsSummarised[["taxaTrajTopSeeds"]][[taxa]][[traj]][["seed"]])
            
            subjectAssignments = kmlSeedsDeconvoluted[[taxa]][["subjectAssignmentsDeconvoluted"]]
            subjectAssignments[["seed"]] = paste("seed:", gsub("*._", "", subjectAssignments[["trajSeedIdentifier"]]))
            
            subjectAssignmentsTopSeeds = subjectAssignments %>% subset(seed %in% topSeeds)
            subjectAssignmentsTopSeedsRobustness = subjectAssignmentsTopSeeds %>% group_by(subject, .data[[cluster_column]]) %>% summarise(count = n()) %>% mutate(prop = count / sum(count)) %>% subset(prop > threshold)
            
            topSeedsSummarised[["topSubjectAssginments"]][[taxa]][[traj]] = merge(subjectAssignmentsTopSeedsRobustness, metadata, by = "subject")
            
#            return(topSeedsSummarised)
            
            ### TILING TOP SEEDS
            
            topSeedsStatistics = topSeedsSummarised[["topSubjectAssginments"]][[taxa]][[traj]] %>% 
                subset(reactivity != "_nodata") %>% 
                group_by(.data[[cluster_column]], reactivity) %>% summarise(count = n()) %>% 
                pivot_wider(names_from = reactivity, values_from = count) %>% data.frame(check.names = FALSE)
            
            topSeedsStatistics[[cluster_column]] <- ifelse(topSeedsStatistics[[cluster_column]] == traj, traj, "X")
            
            topSeedsStatistics = topSeedsStatistics %>% # summarises all other trajectories to X
                group_by(.data[[cluster_column]]) %>% 
                summarise(across(everything(), sum, na.rm = TRUE)) %>% 
                arrange(desc(.data[[cluster_column]])) %>%
                data.frame(check.names = FALSE)
            
            topSeedsStatistics = topSeedsStatistics[, c(cluster_column, "_NR", "_R")]
            rownames(topSeedsStatistics) = topSeedsStatistics[[cluster_column]]
            topSeedsStatistics[[cluster_column]] <- NULL
            
            fisher_test = fisher.test(topSeedsStatistics)
            
            topSeedsStatistics[["trajectory"]] = rownames(topSeedsStatistics)
            topSeedsStatistics[["taxa"]] = taxa
            topSeedsStatistics[["comparison"]] = paste(traj, "v", "X", sep = "")
            topSeedsStatistics[["fisher_pval"]] = fisher_test$p.value
            topSeedsStatistics[["fisher_OR"]] = fisher_test$estimate
            topSeedsStatistics[["fraction_reactive"]] = topSeedsStatistics$`_R` / topSeedsStatistics$`_NR`
            
            topSeedsSummarised[["topSeedsSummarised"]][[taxa]][[traj]] = topSeedsStatistics

            
        }
    }
    
    return(topSeedsSummarised)
}
