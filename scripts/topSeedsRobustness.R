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
            taxaSeedStatistics = kmlSeedStatistics[["taxaSeedStats"]][[taxa]] %>% 
                subset(trajectory == traj & fisher_OR > 1) %>% # only keep seeds where OR > 1 because hypothesis is max
                arrange(fisher_pval) %>% # arranges by pvalues (increasing)
                head(40) # grab 20 best seeds
            
            taxaTrajTopSeeds[[taxa]][[traj]] = taxaSeedStatistics[, c("trajectory", "seed", "fisher_OR", "fisher_pval")]
            
        }
        
        for (traj in min_trajectories) {
            
            taxaTrajTopSeeds[[taxa]][[traj]] = list()
            taxaSeedStatistics = kmlSeedStatistics[["taxaSeedStats"]][[taxa]] %>% 
                subset(trajectory == traj & fisher_OR < 1) %>% # only keep seeds where OR < 1 because hypothesis is min
                arrange(fisher_pval) %>% 
                head(40) # grab 20 best seeds
            
            taxaTrajTopSeeds[[taxa]][[traj]] = taxaSeedStatistics[, c("trajectory", "seed", "fisher_OR", "fisher_pval")]
            
        }
        
    }
    
    return(taxaTrajTopSeeds)
}

summarise_top_seeds <- function(kmlSeedStatistics, kmlSeedsDeconvoluted, metadata, threshold=0.6) {
    # give threshold like a fraction e.g. .75%. needs to be > 0.5 at least.
    # higher threshold does not necessarily produce better results because individuals who are split between two very similar trajectories are filtered away completely
    # In some comparisons it doesn't matter to which of these they belong to e.g. bacteroidia in relab data. 
    
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
            subjectAssignmentsTopSeeds[["deconvoluted_trajectory"]] = subjectAssignmentsTopSeeds[[cluster_column]]
            subjectAssignmentsTopSeeds[[cluster_column]] <- ifelse(subjectAssignmentsTopSeeds[[cluster_column]] == traj, traj, "X")
            subjectAssignmentsTopSeedsRobustness = subjectAssignmentsTopSeeds %>% 
                group_by(subject, .data[[cluster_column]]) %>% 
                summarise(count = n()) %>% 
                mutate(prop = count / sum(count)) %>% 
                subset(prop > threshold)
            
            subjectAssignmentsTopSeedsDeconvolutedLabel <- subjectAssignmentsTopSeeds %>% # attempts to assign individuals to A, B or C instead of A v X etc
                group_by(subject, deconvoluted_trajectory) %>% 
                summarise(count = n()) %>% 
                mutate(prop = count / sum(count)) %>%
                group_by(subject) %>%  # Regroup by subject only to compare within each subject
                filter(prop == max(prop))
            
            subjectAssignmentsTopSeedsRobustness = merge(subjectAssignmentsTopSeedsRobustness, subjectAssignmentsTopSeedsDeconvolutedLabel, by = "subject", all.x = TRUE)
            topSeedsSummarised[["topSubjectAssginments"]][[taxa]][[traj]] = merge(subjectAssignmentsTopSeedsRobustness, metadata, by = "subject")
            

            ### TILING TOP SEEDS
            
            topSeedsStatistics = topSeedsSummarised[["topSubjectAssginments"]][[taxa]][[traj]] %>% 
                subset(reactivity != "_nodata") %>% 
                group_by(.data[[cluster_column]], reactivity) %>% 
                summarise(count = n()) %>% 
                pivot_wider(names_from = reactivity, values_from = count, values_fill = 0) %>% 
                data.frame(check.names = FALSE) %>%
                arrange(desc(.data[[cluster_column]]))
            
            
            topSeedsStatistics = topSeedsStatistics[, c(cluster_column, "_NR", "_R")]
            rownames(topSeedsStatistics) = topSeedsStatistics[[cluster_column]]
            topSeedsStatistics[[cluster_column]] <- NULL
            topSeedsStatistics = topSeedsStatistics + 1
            
            fisher_test = fisher.test(topSeedsStatistics)
            
            topSeedsStatistics[["trajectory"]] = rownames(topSeedsStatistics)
            topSeedsStatistics[["taxa"]] = taxa
            topSeedsStatistics[["comparison"]] = paste(traj, "v", "X", sep = "")
            topSeedsStatistics[["fisher_pval"]] = fisher_test$p.value
            topSeedsStatistics[["fisher_OR"]] = fisher_test$estimate
            topSeedsStatistics[["fisher_OR_reciprocal"]] = 1/fisher_test$estimate
            topSeedsStatistics[["proportion_reactive"]] = topSeedsStatistics$`_R` / (topSeedsStatistics$`_NR` + topSeedsStatistics$`_R`)
            
            topSeedsSummarised[["topSeedsSummarised"]][[taxa]][[traj]] = topSeedsStatistics

        }
    }
    
    return(topSeedsSummarised)
}