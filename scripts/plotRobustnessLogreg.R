#!/usr/bin/env Rscript

######################################################################
# Setting ROBUSTNESS STATISTICS (AGGREGATING ACROSS SEEDS) FUNCTIONS #
######################################################################
########################## LOG REG ###################################
######################################################################

plot_TaxaSeedStats <- function(taxaSeedStats, taxa) {
    
    a <- ggplot() +
    geom_point(data=taxaSeedStats, aes(x=log2(estimate), 
                                             y= -log10(p.value), 
                                             fill = term),
    #                                         size=(.data[["mean_NR"]] + .data[["mean_R"]])),
               shape = 21, size = 7, alpha=0.8) +
    geom_hline(yintercept = 1.3, linetype = 5, lwd = 1.1, alpha=0.7) + 
    
    #geom_text_repel(data = taxaSeedStats,
    #                aes(x = log2(estimate),
    #                    y= -log10(p.value), 
    #                    label = term),
    #                vjust = -1, 
    #                size = 5) +
    #scale_fill_manual(values = c("#fe9929", "#41ab5d")) +
    
    theme_minimal() + 
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(size = 14, face = "bold"), 
        axis.text.y = element_text(size = 14, face = "bold"),
        legend.position = "top") + labs(y = "-log10(pval)", x = "log2 OR (from logreg)", title = paste(taxa, sep = ""))
    
    return(a)
    
}

plot_TaxaSeedsSummed <- function(taxaSeedsSummed) {
    
    taxaSeedsSummed_toplot = taxaSeedsSummed
    taxaSeedsSummed_toplot$estimate <- ifelse(taxaSeedsSummed_toplot$estimate < 0.25, 0.25, taxaSeedsSummed_toplot$estimate)
    taxaSeedsSummed_toplot$estimate <- ifelse(taxaSeedsSummed_toplot$estimate > 4, 4, taxaSeedsSummed_toplot$estimate)
    
    a <- ggplot() +
    geom_point(data=taxaSeedsSummed_toplot, aes(x=log2(estimate), 
                                             y= -log10(p.value)), 
#                                             fill = term),
    #                                         size=(.data[["mean_NR"]] + .data[["mean_R"]])),
               shape = 21, size = 7, alpha=0.8) +
    geom_hline(yintercept = 1.3, linetype = 5, lwd = 1.1, alpha=0.7) + 
    
    geom_text_repel(data = taxaSeedsSummed_toplot,
                    aes(x = log2(estimate),
                        y= -log10(p.value), 
                        label = term),
                    vjust = -1, 
                    size = 5) +
    
    theme_minimal() + 
    theme(
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(size = 14, face = "bold"), 
        axis.text.y = element_text(size = 14, face = "bold"),
        legend.position = "top") + labs(y = "-log10(pval)", x = "log2 OR (from logreg)", title = paste(taxa, sep = ""))
    
    return(a)
    
}

plotDeconvolutedMeanTraj_logreg <- function(kml_deconvoluted_taxa_melted, kmlSeedStatistics_taxa, taxa) { 
    
    cluster_column = paste(taxa, "clusters", sep="_")
    
    meanTraj = kml_deconvoluted_taxa_melted
    meanTraj[["seed"]] = paste("seed:", gsub("*._", "", meanTraj[["trajSeedIdentifier"]]))

    seed_metadata = kmlSeedStatistics_taxa
    seed_metadata[[cluster_column]] = substr(seed_metadata$term, nchar(seed_metadata$term), nchar(seed_metadata$term))
    
#    seed_metadata$term = as.character(seed_metadata$term)
#    seed_metadata[[cluster_column]] = substr(seed_metadata$term, nchar(seed_metadata$term), nchar(seed_metadata$term))
    
    kmlMeanTraj_metadata = merge(meanTraj, seed_metadata, by = c(cluster_column, "seed"))
    
    print(dim(meanTraj))
    print(dim(kmlMeanTraj_metadata))
    
    if (dim(meanTraj)[1] != dim(kmlMeanTraj_metadata)[1]) {
        print(paste(taxa, ": DIMENSIONS DONT ADD UPP AFTER MERGING"))
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
          name="log2(OR)",
          low = "#41ab5d", 
          mid = "white", # This color is at the midpoint, which is 1.5 in this case.
          high = "#fe9929",
          midpoint = 0, limits = c(-2, 2), na.value = "red",
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
        geom_line(data = kmlMeanTraj_metadata, aes(x = .data[["variable_numeric"]], y = .data[["value"]], col = log2(.data[["estimate"]]), group = .data[["trajSeedIdentifier"]])) +
        kml_plot_theme + xlab("timepoint") + ylab("abundance (clr or relab)") + 
        microshades_settings
    
    return(arrangeGrob(a,b,c, ncol = 3))
}

process_subjects <- function(inputrds, seedStats, taxa) {
    
    abundance = melt(inputrds[[taxa]], id.vars = "subject") %>% subset(is.na(value) == FALSE)
    abundance[["variable_numeric"]] = gsub("t", "", abundance[["variable"]])
    abundance[["variable_numeric"]] = as.numeric(abundance[["variable_numeric"]])
    
    abundance_grouped = merge(seedStats[["taxaSubjectAssignments"]][[taxa]], abundance, by = "subject")
    
    cluster_column = paste(taxa, "clusters", sep = "_")
    abundance_grouped[["trajectory"]] = abundance_grouped[[cluster_column]]

    return(abundance_grouped)   
}

process_meanTraj <- function(seedStats, kmlSeedsDeconvoluted, taxa) {
    
    cluster_column = paste(taxa, "clusters", sep = "_")
    
    taxaStats = seedStats[["taxaSeedsSummarised"]][[taxa]]
    taxaStats[[cluster_column]] = substr(taxaStats$term, nchar(taxaStats$term), nchar(taxaStats$term))
    
    trajectoriesAllSeeds = kmlSeedsDeconvoluted[[taxa]][["meltedMeanTraj"]]
    meanTrajAllSeeds = trajectoriesAllSeeds %>% group_by(.data[[cluster_column]], variable_numeric) %>% summarise(mean_value = mean(value))
    
    meanTrajAllSeeds_Stats = merge(meanTrajAllSeeds, taxaStats, by = cluster_column)
    meanTrajAllSeeds_Stats[["trajectory"]] = meanTrajAllSeeds_Stats[[cluster_column]]
    
    return(meanTrajAllSeeds_Stats)
    
}


plot_subjectassignments_with_meanTraj <- function(inputrds, seedStats, kmlSeedsDeconvoluted, taxa) {
    
    subjectTrajectories = process_subjects(inputrds, seedStats, taxa)
    meanTrajTopSeeds = process_meanTraj(seedStats, kmlSeedsDeconvoluted, taxa)
    
    cluster_column = paste(taxa, "clusters", sep = "_")
    
    meanTrajTopSeeds_labels = meanTrajTopSeeds %>% distinct(term, .keep_all = TRUE)
    
    a <- ggplot() +
#        geom_point(data=subjectTrajectories, aes(x=.data[["variable_numeric"]], y=.data[["value"]]), alpha = 0.2) + 
        geom_line(data=subjectTrajectories, aes(x=.data[["variable_numeric"]], y=.data[["value"]], group=subject), lwd=1.1, alpha=0.2) +
    
        geom_line(data=meanTrajTopSeeds, aes(x=.data[["variable_numeric"]], y=.data[["mean_value"]], col = log2(.data[["estimate"]]), group = .data[[cluster_column]]), lwd = 4.5, alpha = 0.8) +
        geom_line(data=meanTrajTopSeeds, aes(x=.data[["variable_numeric"]], y=.data[["mean_value"]], col = log2(.data[["estimate"]]), group = .data[[cluster_column]]), lwd = 7, alpha = 0.5) +
    
        geom_label_repel(data = meanTrajTopSeeds_labels, 
                        aes(x = .data[["variable_numeric"]], 
                            y= .data[["mean_value"]], 
                            label = paste(paste("OR: ", round(.data[["estimate"]], digits = 3), sep = ""), 
                                          paste("pval: ", round(.data[["p.value"]], digits = 6), sep = ""),
                                          sep = "\n")),
                        vjust = 0, 
                        size = 6,
                        box.padding = 0.5) +
    
        scale_color_gradient2(
          name="log2(OR)",
          low = '#41ab5d', 
          mid = "white", # This color is at the midpoint, which is 1.5 in this case.
          high = '#fe9929',
          midpoint = 0, limits = c(-2, 2), na.value = "red",
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
            legend.position = "top") + ggtitle(paste(taxa, sep = ""))
    
    
    return(a)
}