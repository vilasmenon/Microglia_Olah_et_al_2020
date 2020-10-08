####Code workflow for clustering of cells in Olah et al. (2020)#####
####Note that this code requires an older version of Seurat (v2), and takes a long time to run in its entirety###

require(Seurat)
require(randomForest)
require(Matrix)
memory.limit(80000)
source("code_for_submission_functions.R")

###Step 1: load data, remove genes often associated with technical artifacts, and filter by total UMI count###
alldat = read.csv("counts_matrix.csv",as.is=T,row.names=1,header=T)
alldat = alldat[grep("^LOC|^MT-|^RP[0-9]|^BC[0-9]|-PS",rownames(alldat),invert=T),]
alldat = alldat[,which(colSums(alldat)>=1000)]
batchval = do.call(rbind,strsplit(colnames(alldat),"_"))

###Step 2: run first-level clustering, ranging over parameters
###Note that the cluster IDs are saved as the names of the list, not as the values
top_clustering_over_parameters = cluster_over_parameters(alldat,batchval,keepcols=1:ncol(alldat),pcrange=5:15,
                                                       resrange=seq(from=0.2,to=0.8,by=0.2),krange=c(5,10,20,30),batchreg=1)

###Step 3: assess cluster robustness for every parameter set and select the one with the highest number of clusters
###with minimum prediction score >=0.75
cluster_numbers=rep(0,length(top_clustering_over_parameters))
predmatrices=list
for (parval in 1:length(top_clustering_over_parameters)) {
  clusterids=strsplit(names(top_clustering_over_parameters)[parval],",")[[1]]
  prediction_matrix=cluster_robustness(distweights=alldat,clids=clusterids,num_iterations=20)
  minrows=apply(prediction_matrix,1,min)
  mincols=apply(prediction_matrix,2,min)
  cluster_numbers[parval]=length(which(minrows>=0.75 & mincols>=0.75))
  predmatrices[[parval]]=prediction_matrix
}
optimal_pars=which(cluster_numbers==max(cluster_numbers))[1]
top_level_clusters=strsplit(names(top_clustering_over_parameters)[optimal_pars],",")
prediction_matrix=predmatrices[[optimal_pars]]
minrows=apply(prediction_matrix,1,min)
mincols=apply(prediction_matrix,2,min)
distinct_clusters=rownames(prediction_matrix)[which(minrows>=0.75 & mincols>=0.75)]
merge_cells=which(top_level_clusters %in% setdiff(rownames(prediction_matrix),distinct_clusters))
top_level_clusters[merge_cells]="merged1"

###show distribution of top level clusters
table(top_level_clusters)

###Step 4: for each top level cluster, re-run steps 2 and 3
subcluster_list=list()
unique_top_level_clusters=unique(top_level_clusters)
for (topcl in unique_top_level_clusters) {
  keepcells=which(top_level_clusters==topcl)
  sub_clustering_over_parameters = cluster_over_parameters(alldat,batchval,keepcols=keepcells,pcrange=5:15,
                                                           resrange=seq(from=0.2,to=0.8,by=0.2),krange=c(5,10,20,30),batchreg=1)
  cluster_numbers=rep(0,length(sub_clustering_over_parameters))
  predmatrices=list
  for (parval in 1:length(sub_clustering_over_parameters)) {
    clusterids=strsplit(names(sub_clustering_over_parameters)[parval],",")[[1]]
    prediction_matrix=cluster_robustness(distweights=alldat[,keepcells],clids=clusterids,num_iterations=20)
    minrows=apply(prediction_matrix,1,min)
    mincols=apply(prediction_matrix,2,min)
    cluster_numbers[parval]=length(which(minrows>=0.75 & mincols>=0.75))
    predmatrices[[parval]]=prediction_matrix
  }
  if (max(cluster_numbers)>1) {
    optimal_pars=which(cluster_numbers==max(cluster_numbers))[1]
    sub_level_clusters=strsplit(names(sub_clustering_over_parameters)[optimal_pars],",")
    prediction_matrix=predmatrices[[optimal_pars]]
    minrows=apply(prediction_matrix,1,min)
    mincols=apply(prediction_matrix,2,min)
    distinct_clusters=rownames(prediction_matrix)[which(minrows>=0.75 & mincols>=0.75)]
    merge_cells=which(sub_level_clusters %in% setdiff(rownames(prediction_matrix),distinct_clusters))
    sub_level_clusters[merge_cells]="merged1"
    subcluster_list[[as.character(topcl)]]=sub_level_clusters
  } else {
    subcluster_list[[as.character(topcl)]]="No_subdivisions"
  }
}


###compile top level and sub-level cluster names###
two_level_clusters=as.character(top_level_clusters)
for (ii in names(subcluster_list)) {
  if (subcluster_list[[ii]][1]!="No_subdivisions") {
    renamecells=which(two_level_clusters==ii)
    two_level_clusters[renamecells]=paste0(two_level_clusters[renamecells],"_",subcluster_list[[ii]])
  }
}
table(two_level_clusters)


