#!/usr/bin/env Rscript

######################################################################
# Setting ROBUSTNESS STATISTICS (AGGREGATING ACROSS SEEDS) FUNCTIONS #
######################################################################
########################## LOG REG ###################################
######################################################################

tile_seed_stats_logreg <- function(subjectAssignmentsDeconvoluted, taxa) {
    
    seed_stats = list()
    subjectAssignmentsDeconvoluted$seed <- paste("seed:", gsub("*._", "", subjectAssignmentsDeconvoluted$trajSeedIdentifier))
    seeds = unique(subjectAssignmentsDeconvoluted$seed)
    
    cluster_column = paste(taxa, "_clusters", sep = "")
    clusters = unique(subjectAssignmentsDeconvoluted[[cluster_column]])
    
    for (seed_val in seeds) {
        
        seedDeconvoluted = subjectAssignmentsDeconvoluted %>% 
            subset(seed == seed_val & reactivity != "_nodata") %>%
            data.frame(check.names = FALSE)
        
        
        clusters = unique(seedDeconvoluted[[cluster_column]])
        
        if (length(clusters) == 1) {
            put(paste("SKIPPING SEED: ", seed_val, " only one trajectory"), sep = "")
            next
            
        } else {
        
            for (cluster in clusters) {
                
                nSubjects = subsetseedDeconvoluted %>% group_by(.data[[cluster_column]]) %>% filter(.data[[cluster_column]] == cluster) %>% dim()                                       
                cluster_comparison = paste(cluster, "vX", sep = "")
                seedDeconvoluted[[cluster_comparison]] <- ifelse(seedDeconvoluted[[cluster_column]] == cluster, cluster, "X")
                seedDeconvoluted[[cluster_comparison]] = factor(seedDeconvoluted[[cluster_comparison]], c("X", cluster))
                
                formula_ = paste("reactivity_coded ~ ", cluster_comparison, sep = "")
                seed_logmodel <- glm(formula_, family=binomial(link='logit'), data=seedDeconvoluted)
                seed_logmodel_data = plot_model(seed_logmodel)$data
                seed_logmodel_data[["term"]] = paste(taxa, seed_logmodel_data[["term"]], sep = "_")
                seed_logmodel_data[["seed"]] = seed_val
                seed_logmodel_data[["taxa"]] = taxa
                seed_logmodel_data[["cluster"]] = cluster
                seed_logmodel_data[["nSubjects"]] = nSubjects[1]
                
                
                seed_stats[[paste(taxa, cluster, seed_val, sep = "_")]] = seed_logmodel_data
            
            }
            
        }
    }
    
    return(seed_stats)
}

retrieve_subject_assignments <- function(subjectAssignmentsDeconvoluted, taxa) {
    
    cluster_column = paste(taxa, "_clusters", sep = "")
    
    subjectAssignmentsTaxa = subjectAssignmentsDeconvoluted %>% 
        group_by(subject, .data[[cluster_column]]) %>%
        summarise(count = n()) %>% 
        mutate(prop = count / sum(count)) %>%
        group_by(subject) %>%  # Regroup by subject only to compare within each subject
        filter(prop == max(prop))
    
    return(subjectAssignmentsTaxa)
}

subject_assignments_logreg <- function(subjectAssignmentsTaxa, metadata, taxa) {
    
    cluster_column = paste(taxa, "_clusters", sep = "")
    clusters = unique(subjectAssignmentsTaxa[[cluster_column]])
    
    subjectAssignmentsTaxa = merge(subjectAssignmentsTaxa, metadata, by = "subject") %>% subset(reactivity != "_nodata")
    subjectAssignmentsTaxa[["reactivity_coded"]] <- ifelse(subjectAssignmentsTaxa[["reactivity"]] == "_R", 1,
                                                    ifelse(subjectAssignmentsTaxa[["reactivity"]] == "_NR", 0, 2))
    
    taxa_stats = list()
    
    
    for (cluster in clusters) {
        
        nSubjects = subjectAssignmentsTaxa %>% group_by(.data[[cluster_column]]) %>% filter(.data[[cluster_column]] == cluster) %>% dim()                                       
        cluster_comparison = paste(cluster, "vX", sep = "")
        subjectAssignmentsTaxa[[cluster_comparison]] <- ifelse(subjectAssignmentsTaxa[[cluster_column]] == cluster, cluster, "X")
        subjectAssignmentsTaxa[[cluster_comparison]] = factor(subjectAssignmentsTaxa[[cluster_comparison]], c("X", cluster))
        
        formula_ = paste("reactivity_coded ~ ", cluster_comparison, sep = "")
        taxa_logmodel <- glm(formula_, family=binomial(link='logit'), data=subjectAssignmentsTaxa)
        taxa_logmodel_data = plot_model(taxa_logmodel)$data
        taxa_logmodel_data[["term"]] = paste(taxa, taxa_logmodel_data[["term"]], sep = "_")
        taxa_logmodel_data[["taxa"]] = taxa
        taxa_logmodel_data[["cluster"]] = cluster
        taxa_logmodel_data[["nSubjects"]] = nSubjects[1]
        
        
        taxa_stats[[paste(taxa, cluster, sep = "_")]] = taxa_logmodel_data
        
        }
    
    taxa_stats = do.call(rbind, taxa_stats)
    
    return(taxa_stats)
    
}

seed_statistics_logreg <- function(kmlSeedsDeconvoluted, metadata) {
    
    taxa_list = names(kmlSeedsDeconvoluted)
    
    robustness_statistics = list()
    robustness_statistics[["taxaSeedStats"]] = list()
    robustness_statistics[["taxaSeedsSummarised"]] = list()
    
    for (taxa in taxa_list) {
        
        taxa_deconvoluted = merge(kmlSeedsDeconvoluted[[taxa]][["subjectAssignmentsDeconvoluted"]], metadata, by = "subject")
        taxa_deconvoluted[["reactivity_coded"]] <- ifelse(taxa_deconvoluted[["reactivity"]] == "_R", 1, 
                                                   ifelse(taxa_deconvoluted[["reactivity"]] == "_NR", 0, 2))
        
        robustness_statistics[["taxaSeedStats"]][[taxa]] = do.call(rbind, tile_seed_stats_logreg(taxa_deconvoluted, taxa))
        
        robustness_statistics[["taxaSubjectAssignments"]][[taxa]] = retrieve_subject_assignments(taxa_deconvoluted, taxa)
        
        robustness_statistics[["taxaSeedsSummarised"]][[taxa]] = subject_assignments_logreg(robustness_statistics[["taxaSubjectAssignments"]][[taxa]], metadata, taxa)

        
    }
    
    return(robustness_statistics)
}
