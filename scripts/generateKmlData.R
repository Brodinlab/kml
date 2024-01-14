#!/usr/bin/env Rscript

##########################
# Setting kml Parameters #
##########################

kml_parameters = list()
kml_parameters[["parWithMinkowski"]] = parALGO(distanceName="minkowski", startingCond="randomAll")
kml_parameters[["parWithMax"]] = parALGO(distanceName = "maximum", startingCond = "maxDist")
kml_parameters[["parWithEuclidean"]] = parALGO(distanceName = "euclidean", startingCond = "maxDist")
kml_parameters[["parWithEuclidean_rndm"]] = parALGO(distanceName = "euclidean", startingCond = "all")
kml_parameters[["parWithMedian"]] = parALGO(distanceName = "maximum", "startingCond" = "maxDist", centerMethod=function(x){median(x,na.rm=TRUE)})

#########################################
# Setting KML DATA EXTRACTION FUNCTIONS #
#########################################

get_CLD <- function(df, subject, tpcol, parAlgo_) {
    
    dfCLD <<- cld(data.frame(df), timeInData = tpcol:ncol(df), maxNA = 5)
    capture.output(kml(dfCLD, nbRedrawing=100, parAlgo=parAlgo_), file="trash")
    
    return(dfCLD)
    
}

get_meanTraj <- function(dfCLD, nclustersStr, clusterrank) {
    
    meanTraj <<- calculTrajMean(dfCLD["traj"], dfCLD[nclustersStr][[clusterrank]]['clusters'])
    meanTraj <- meanTraj %>% 
        as.data.frame()
    
    meanTraj$cluster <- row.names(meanTraj)
    meanTraj <- reshape2::melt(meanTraj, id.var = c('cluster'))

    return(meanTraj)
}

get_IdClusterData <- function(dfCLD, nclusters, clusterrank, idcol, antigencol) {

    list_ID <<- list(dfCLD['idAll'])
    list_cluster <<- list(getClusters(dfCLD, nbCluster=nclusters, clusterRank = clusterrank, asInteger = FALSE))
    IdCluster_merged <- do.call(rbind, Map(data.frame, subject=list_ID, antigen=list_cluster))
    
    colnames(IdCluster_merged) <- list(idcol, antigencol)

    return(IdCluster_merged)
    
}

get_subjectTraj <- function(df, tpcol, dfCLD, IdCluster_merged, subject) {


    subjectTraj <<- merge(x=df, y=IdCluster_merged, by=subject)
    colnames(subjectTraj) <- append(list(subject), append(colnames(dfCLD['traj']), "cluster")) # lol at this solution
    subjectTraj <- reshape2::melt(subjectTraj, id.var = c("cluster", subject))
    
    return(subjectTraj)
}

################################
# Setting KML RUNNING FUNCTION #
################################

run_kml <- function(taxa_wide, nclusters, clusterrank, parAlgo_, taxa, subset, transformation) {
    
    # first column is the unique identifier column
    
    nclusters_str = paste("c", nclusters, sep = "")
    identifier_col = colnames(taxa_wide)[1]
    
    taxaCLD = get_CLD(taxa_wide, identifier_col, tpcol=2, parAlgo_=parAlgo_)
    taxaMeanTraj = get_meanTraj(taxaCLD, nclusters_str, clusterrank) # mean trajectory path
    taxaIDCluster = get_IdClusterData(taxaCLD, nclusters, clusterrank, identifier_col, "cluster") # subject cluster assignment. 
    taxaSubTraj = get_subjectTraj(taxa_wide, tpcol=2, taxaCLD, taxaIDCluster, identifier_col) # don't think we need this. 
    
    
    colnames(taxaIDCluster)[2] = paste(taxa, "clusters", sep="_")
    colnames(taxaMeanTraj)[1] = paste(taxa, "clusters", sep="_")
    
    return(
        list(
            subjectAssignment = taxaIDCluster, 
            taxaMeanTraj = taxaMeanTraj)
    )   
}

