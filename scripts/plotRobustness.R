#!/usr/bin/env Rscript

####################################################################
# Setting ROBUSTNESS PLOTTING (AGGREGATING ACROSS SEEDS) FUNCTIONS #
####################################################################

plot_taxaSeedsSummarised <- function(seed_stats, plot = TRUE) {
    
    taxaSeedsSummarised = do.call(rbind, seed_stats[["taxaSeedsSummarised"]]) %>% subset(filter == 'keep')
    taxaSeedsSummarised$taxa <- paste(sapply(strsplit(rownames(taxaSeedsSummarised), ".", fixed = TRUE), function(x) {unlist(x)[2]}), taxaSeedsSummarised$trajectory, sep = "_")
    
    if (plot == TRUE) {
        
        a = ggplot() +
        geom_point(data=taxaSeedsSummarised, aes(x=OR_mean, y=-log10(p.val), fill = ifelse(OR_mean > 1, "higher risk of atopy", "lower risk of atopy")), size=8, shape = 21, alpha=0.8) +
        geom_hline(yintercept = 1.3, linetype = 5, lwd = 1.1, alpha=0.7) + 

        geom_text_repel(data = taxaSeedsSummarised, 
                        aes(x = OR_mean, 
                            y=-log10(p.val), 
                            label = taxa),
                        vjust = -1, 
                        size = 5) +

        scale_fill_manual(values = c("grey", "blue")) +
        
        theme_minimal() + 
        theme(
            panel.border = element_rect(color = "black", fill = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"), 
            axis.text.x = element_text(size = 14, face = "bold"), 
            axis.text.y = element_text(size = 14, face = "bold")) +
        xlab("OR for being atopic in specific trajectory") + ylab("-log10(pvalue)") + ggtitle("atopy in association to longitudinal microbial trajectories")
        
        return(a)
    } else {
        
        return(taxaSeedsSummarised)
    }
    
}

plotDeconvolutedMeanTraj <- function(kml_deconvoluted_taxa_melted, kmlSeedStatistics_taxa, taxa) { 
    
    cluster_column = paste(taxa, "clusters", sep="_")
    
    meanTraj = kml_deconvoluted_taxa_melted
    meanTraj[["seed"]] = paste("seed:", gsub("*._", "", meanTraj[["trajSeedIdentifier"]]))
    
    seed_metadata = kmlSeedStatistics_taxa %>% 
        subset(trajectory != 'X') %>% 
        select(trajectory, seed, fraction_reactive, fisher_OR, fisher_pval) %>%
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
          name="Fraction Reactive",
          low = microshades_palettes$micro_green[4], 
          mid = "white", # This color is at the midpoint, which is 1.5 in this case.
          high = microshades_palettes$micro_orange[4], 
          midpoint = 1.5, limits = c(0.5, 2.5), na.value = microshades_palettes$micro_orange[4],
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
        geom_line(data = kmlMeanTraj_metadata, aes(x = .data[["variable_numeric"]], y = .data[["value"]], col = .data[["fraction_reactive"]], group = .data[["trajSeedIdentifier"]])) +
        kml_plot_theme + xlab("timepoint") + ylab("abundance (clr or relab)") + 
        microshades_settings
    
    
    
    return(arrangeGrob(a,b,c, ncol = 3))
}
