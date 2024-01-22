#!/usr/bin/env Rscript

####################################################################
# Setting ROBUSTNESS PLOTTING (AGGREGATING ACROSS SEEDS) FUNCTIONS #
####################################################################

plot_taxaSeedsSummarised <- function(seed_stats, plot = TRUE) {
    
    taxaSeedsSummarised = do.call(rbind, seed_stats[["taxaSeedsSummarised"]]) %>% subset(filter == 'keep')
    taxaSeedsSummarised$taxa <- paste(sapply(strsplit(rownames(taxaSeedsSummarised), ".", fixed = TRUE), function(x) {unlist(x)[2]}), taxaSeedsSummarised$trajectory, sep = "_")

    if (plot == TRUE) {
        
        a = ggplot() +
        geom_point(data=taxaSeedsSummarised, aes(x=log2(OR_mean), 
                                                 y= count / 100, 
                                                 fill = ifelse(OR_mean > 1, "higher risk of atopy", "lower risk of atopy"),
                                                 size=(.data[["mean_NR"]] + .data[["mean_R"]])),
                   shape = 21, alpha=0.8) +

        geom_hline(yintercept = 1.3, linetype = 5, lwd = 1.1, alpha=0.7) + 

        geom_text_repel(data = taxaSeedsSummarised, 
                        aes(x = log2(OR_mean),
                            y= count / 100, 
                            label = taxa),
                        vjust = -1, 
                        size = 5) +

        scale_fill_manual(values = c("#fe9929", "#41ab5d")) +
        
        theme_minimal() + 
        theme(
            panel.border = element_rect(color = "black", fill = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"), 
            axis.text.x = element_text(size = 14, face = "bold"), 
            axis.text.y = element_text(size = 14, face = "bold")) +
        xlab("log2(mean_OR) across all seeds for this trajectory") + 
        ylab("% seeds this cluster performed best / worst") + 
        ggtitle("trajectory consistency across seeds")
        
        return(a)
        
    } else {
        
        return(taxaSeedsSummarised)
    }
    
}

plotDeconvolutedMeanTraj <- function(kml_deconvoluted_taxa_melted, kmlSeedStatistics_taxa, taxa) { 
    
    cluster_column = paste(taxa, "clusters", sep="_")
    
    meanTraj = kml_deconvoluted_taxa_melted
    meanTraj[["seed"]] = paste("seed:", gsub("*._", "", meanTraj[["trajSeedIdentifier"]]))
    
    kmlSeedStatistics_taxa[["OR_capped"]] = ifelse(kmlSeedStatistics_taxa[["fisher_OR"]] >= 4, 4, 
                                                   ifelse(kmlSeedStatistics_taxa[["fisher_OR"]] <= 0.25, 0.25, kmlSeedStatistics_taxa[["fisher_OR"]]))
    
    seed_metadata = kmlSeedStatistics_taxa %>% 
        subset(trajectory != 'X') %>% 
        select(trajectory, seed, OR_capped, fisher_OR, fisher_pval) %>%
        dplyr::rename(!!cluster_column := "trajectory")
    
    kmlMeanTraj_metadata = merge(meanTraj, seed_metadata, by = c(cluster_column, "seed"))
    
    if (dim(meanTraj)[1] != dim(kmlMeanTraj_metadata)[1]) {
        put(paste(taxa, ": DIMENSIONS DONT ADD UPP AFTER MERGING"))
    }
    
    kml_plot_theme = theme_minimal() + 
        theme(
            panel.border = element_rect(color = "black", fill = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"), 
            axis.text.x = element_text(size = 14, face = "bold"), 
            axis.text.y = element_text(size = 14, face = "bold"), legend.position = "top")
    
    microshades_settings = scale_color_gradient2(
          name="log2(OR) capped at -2, 2",
          low = '#41ab5d', 
          mid = "white", # This color is at the midpoint, which is 1:1 (log2 = 0) in this case.
          high = '#fe9929', 
          midpoint = 0, limits = c(-2, 2), na.value = 'red',
          space = "Lab",
          guide = "colourbar"
        )
        
    a <- ggplot() + # the original labeling
        geom_line(data = kmlMeanTraj_metadata, aes(x = .data[["variable_numeric"]], y = .data[["value"]], col = .data[["original_label"]], group = .data[["trajSeedIdentifier"]])) +
        kml_plot_theme + xlab("timepoint") + ylab("abundance (clr or relab)")
    
    b <- ggplot() + # deconvoluted labeling
        geom_line(data = kmlMeanTraj_metadata, aes(x = .data[["variable_numeric"]], y = .data[["value"]], col = .data[[cluster_column]], group = .data[["trajSeedIdentifier"]])) +
        kml_plot_theme + xlab("timepoint") + ylab("abundance (clr or relab)")
    
    c <- ggplot() + # deconvoluted labeling
        geom_line(data = kmlMeanTraj_metadata, aes(x = .data[["variable_numeric"]], y = .data[["value"]], col = log2(.data[["OR_capped"]]), group = .data[["trajSeedIdentifier"]])) +
        kml_plot_theme + xlab("timepoint") + ylab("abundance (clr or relab)") + 
        microshades_settings
    
    
    
    return(arrangeGrob(a,b,c, ncol = 3))
}

################################################################
# Setting ROBUSTNESS PLOTTING FOR FINAL TRAJECTORIES FUNCTIONS #
################################################################

plot_topSeedsSummarised <- function(topSeedsSummarised) {
    
    topSeedsSummarised$taxa = paste(topSeedsSummarised$taxa, topSeedsSummarised$trajectory, sep = "_")
    topSeedsSummarised$color_column <- ifelse(topSeedsSummarised$fisher_pval <= 0.05 & topSeedsSummarised$fisher_OR <= 1, '#238b45',
                                              ifelse(topSeedsSummarised$fisher_pval <= 0.05 & topSeedsSummarised$fisher_OR >= 1, '#ff7f00', '#525252'))

    color_vector <- setNames(topSeedsSummarised$color_column, topSeedsSummarised$taxa)
    
    a = ggplot() +
        geom_point(data=topSeedsSummarised, aes(x=log2(fisher_OR), y=-log10(fisher_pval), 
                                    fill = taxa, 
                                    size=.data[["_NR"]] + .data[["_R"]]), 
                   shape= 21, alpha=0.6) +

        geom_hline(yintercept = 1.3, linetype = 5, lwd = 1.1, alpha=0.7) + 

        geom_text_repel(data = topSeedsSummarised, 
                        aes(x = log2(fisher_OR), 
                            y=-log10(fisher_pval), 
                            label = taxa),
                        vjust = -1, 
                        size = 5) +

        scale_fill_manual(values = color_vector, aesthetics = c("colour", "fill")) +
        
        theme_minimal() + 
        theme(
            panel.border = element_rect(color = "black", fill = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"), 
            axis.text.x = element_text(size = 14, face = "bold"), 
            axis.text.y = element_text(size = 14, face = "bold"),
            legend.position = "none") +
        xlab("log2(OR)") + ylab("-log10(pvalue)") + ggtitle("atopy in association to longitudinal microbial trajectories")
    
    return(a)
}

#######################################################################
# Setting ROBUSTNESS PLOTTING FOR FINAL SUBJECT ASSIGNMENTS FUNCTIONS #
#######################################################################

retrieve_proportion_reactive <- function(topsubjectAssignmentTaxaTraj, taxa) {
    
    cluster_column = paste(taxa, "clusters", sep = "_")
    
    topsubjectAssignmentTaxaTraj_wide = topsubjectAssignmentTaxaTraj %>% 
        group_by(deconvoluted_trajectory, reactivity) %>% 
        tally %>% 
        pivot_wider(names_from = reactivity, values_from = n) %>%
        rename(!!cluster_column := deconvoluted_trajectory)

    topsubjectAssignmentTaxaTraj_wide[["deconvoluted_proportion_reactive"]] = topsubjectAssignmentTaxaTraj_wide$`_R` / (topsubjectAssignmentTaxaTraj_wide$`_NR` + topsubjectAssignmentTaxaTraj_wide$`_R`)

    return(topsubjectAssignmentTaxaTraj_wide)
}

process_subjects <- function(inputrds, topSeeds, taxa, traj) {
    
    abundance = melt(inputrds[[taxa]], id.vars = "subject") %>% subset(is.na(value) == FALSE)
    abundance[["variable_numeric"]] = gsub("t", "", abundance[["variable"]])
    abundance[["variable_numeric"]] = as.numeric(abundance[["variable_numeric"]])
    
    abundance_grouped = merge(topSeeds[["topSubjectAssginments"]][[taxa]][[traj]], abundance)
    
    cluster_column = paste(taxa, "clusters", sep = "_")
    abundance_grouped[["trajectory"]] = abundance_grouped[[cluster_column]]

    return(abundance_grouped)   
}


process_meanTraj <- function(topSeeds, kmlSeedsDeconvoluted, taxa, traj) {
    
    cluster_column = paste(taxa, "clusters", sep = "_")
    
    topSeedsTaxaTraj = topSeeds[["taxaTrajTopSeeds"]][[taxa]][[traj]]
    
    trajectoriesAllSeeds = kmlSeedsDeconvoluted[[taxa]][["meltedMeanTraj"]]
    trajectoriesAllSeeds[["seed"]] = paste("seed:", gsub("*._", "", trajectoriesAllSeeds[["trajSeedIdentifier"]]))
    trajectoriesTopSeeds = merge(trajectoriesAllSeeds, topSeedsTaxaTraj, by = "seed")
    meanTrajTopSeeds = trajectoriesTopSeeds %>% group_by(.data[[cluster_column]], variable_numeric) %>% summarise(mean_value = mean(value))
    
    
    meanTrajTopSeeds[["trajectory"]] <- ifelse(meanTrajTopSeeds[[cluster_column]] == traj, traj, "X")
    meanTrajTopSeeds_OR = merge(meanTrajTopSeeds, topSeeds[["topSeedsSummarised"]][[taxa]][[traj]], by = c("trajectory"))
    
    deconvoluted_reactive_prop = retrieve_proportion_reactive(topSeeds[["topSubjectAssginments"]][[taxa]][[traj]], taxa)
    
    meanTrajTopSeeds_OR_reactive_prop = merge(meanTrajTopSeeds_OR, deconvoluted_reactive_prop, by = cluster_column, all.x = TRUE)

    return(meanTrajTopSeeds_OR_reactive_prop)
    
}


plot_subjectassignments_with_meanTraj <- function(inputrds, topSeeds, kmlSeedsDeconvoluted, taxa, traj) {
    
    subjectTrajectories = process_subjects(inputrds, topSeeds, taxa, traj)
    meanTrajTopSeeds = process_meanTraj(topSeeds, kmlSeedsDeconvoluted, taxa, traj)
    
    meanTrajTopSeeds[["OR"]] <- ifelse(meanTrajTopSeeds[["trajectory"]] == 'X', meanTrajTopSeeds[["fisher_OR_reciprocal"]], meanTrajTopSeeds[["fisher_OR"]])
    
    cluster_column = paste(taxa, "clusters", sep = "_")
    
    meanTrajTopSeeds_labels = meanTrajTopSeeds %>% distinct(trajectory, .keep_all = TRUE)
    
    a <- ggplot() +
#        geom_point(data=subjectTrajectories, aes(x=.data[["variable_numeric"]], y=.data[["value"]]), alpha = 0.2) + 
        geom_line(data=subjectTrajectories, aes(x=.data[["variable_numeric"]], y=.data[["value"]], group=subject), lwd=1.1, alpha=0.2) +
    
        geom_line(data=meanTrajTopSeeds, aes(x=.data[["variable_numeric"]], y=.data[["mean_value"]], col = .data[["deconvoluted_proportion_reactive"]], group = .data[[cluster_column]]), lwd = 4.5, alpha = 0.8) +
        geom_line(data=meanTrajTopSeeds, aes(x=.data[["variable_numeric"]], y=.data[["mean_value"]], col = .data[["deconvoluted_proportion_reactive"]], group = .data[[cluster_column]]), lwd = 7, alpha = 0.5) +
    
        geom_label_repel(data = meanTrajTopSeeds_labels, 
                        aes(x = .data[["variable_numeric"]], 
                            y= .data[["mean_value"]], 
                            label = paste(paste("OR ", traj, "vX: " , round(.data[["OR"]], digits = 3), sep = ""), 
                                          paste("pval ", traj, "vX: " round(.data[["fisher_pval"]], digits = 6), sep = ""),
                                          sep = "\n")),
                        vjust = 0, 
                        size = 6,
                        box.padding = 0.5) +
    
        scale_color_gradient2(
          name="proportion atopic",
          low = microshades_palettes$micro_green[4], 
          mid = "white", # This color is at the midpoint, which is 1.5 in this case.
          high = microshades_palettes$micro_orange[4], 
          midpoint = 0.5, limits = c(0, 1), na.value = "red",
          space = "Lab",
          guide = "colourbar"
        ) +
    
        facet_wrap(~trajectory) + 
        theme(
            panel.background = element_blank(),
            panel.border = element_rect(color = "black", fill = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"), 
            axis.text.x = element_text(size = 14, face = "bold"), 
            axis.text.y = element_text(size = 14, face = "bold"), 
            legend.position = "top") + ggtitle(paste(taxa, " ", traj, sep = ""))
    
    
    return(a)
}