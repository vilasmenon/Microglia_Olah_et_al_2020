###functions called by main script code_for_submission.R###


require(Seurat)   ###note: Seurat version 2 was used for the clustering analysis in the accompanying paper Olah et al. (2020)
require(randomForest)

###Cluster over parameters (number of PCs and resolutions)###
cluster_over_parameters=function(dat,batchval,keepcols,pcrange=5:15,resrange=seq(from=0.2,to=0.8,by=0.2),krange=c(5,10,20,30),batchreg=1) {
  fullclusterings=list()
  if (length(keepcols)>5) {
    datobj_0 <- new("seurat", raw.data = data.frame(dat[,keepcols]))
    datobj_0 <- CreateSeuratObject(raw.data = data.frame(as.matrix(alldat[,keepcols])), min.cells = 3, min.genes = 1)
    datobj_0 <- AddMetaData(object = datobj_0, metadata = batchval[keepcols,])
    datobj_0 = NormalizeData(datobj_0)
    if (batchreg==1 & (length(unique(datobj_0@meta.data$orig.ident))>1)) {
      datobj_0 = ScaleData(datobj_0,vars.to.regress=c("nUMI","batchid"))
    } else {
      datobj_0 = ScaleData(datobj_0,vars.to.regress=c("nUMI"))
    }
    vargenes=rownames(datobj_0@raw.data)[which(apply(datobj_0@raw.data,1,var)>apply(datobj_0@raw.data,1,mean))]
    datobj_0 <- RunPCA(object = datobj_0, pc.genes = vargenes, do.print = FALSE)
    comboval="0"
    maxcls=1
    maxpcs=0
    maxident=rep("0",length(keepcols))
    for (num_pcs in pcrange) {
      for (resval in resrange) {
        for (k.param in krange) {
          datobj_0 <- FindClusters(datobj_0, dims.use = 1:num_pcs, reduction.type="pca", force.recalc=TRUE, algorithm=algval, print.output = 0, resolution = resval, save.SNN = F)
          clutstids=as.numeric(datobj_0@ident)
          clustids=datobj_0@active.ident
          print(paste("running parameter combination ",num_pcs,resval,k.param,length(unique(clustids)),sep=" "))
          if (length(unique(clustids))>1) {
            outnam=paste(clustids,collapse=",")
            if (outnam %in% names(fullclusterings)) {
              fullclusterings[[outnam]]=c(fullclusterings[[outnam]],paste(num_pcs,resval,k.param))
            } else {
              fullclusterings[[outnam]]=paste(num_pcs,resval,k.param)
            }
          }
        }
      }
    }
  }
  return(fullclusterings)
}


####function to assess cluster robustness - returns a matrix of minimum prediction values for each cluster over all pairs of clusters###
cluster_robustness=function(distweights,clids,num_iterations=20) {
  allcl=unique(clids)
  allcl=allcl[order(allcl)]
  clmat=matrix(1,nrow=length(allcl),ncol=length(allcl))
  rownames(clmat)=allcl
  colnames(clmat)=allcl
  for (cl1 in 1:(length(allcl)-1)) {
    clinds1=which(clids==allcl[cl1])
    samp1=round(length(clinds1)*0.5)
    for (cl2 in (cl1+1):length(allcl)) {
      clinds2=which(clids==allcl[cl2])
      samp2=round(length(clinds2)*0.5)
      prediction_vals=c(1,1)
      for (rf_iteration in 1:num_iterations) {
        numcells=min(c(10,round(samp1/2),round(samp2/2)))
        set.seed(10000*cl1+100*cl2+rf_iteration)
        sampcells1=sample(clinds1,samp1)
        sampcells2=sample(clinds2,samp2)
        remcells1=setdiff(clinds1,sampcells1)
        remcells2=setdiff(clinds2,sampcells2)
        rfmod1=randomForest(x=distweights[c(sampcells1,sampcells2),],y=as.factor(clids[c(sampcells1,sampcells2)]))
        rfout1=predict(rfmod1,distweights[c(remcells1,remcells2),])
        ###confusion matrix###
        confmat=table(rfout1,clids[c(remcells1,remcells2)])
        prediction_vals[1]=min(prediction_vals[1],confmat[1,1]/(confmat[1,1]+confmat[2,1]))
        prediction_vals[2]=min(prediction_vals[2],confmat[2,2]/(confmat[1,2]+confmat[2,2]))
      }
      clmat[cl1,cl2]=prediction_vals[1]
      clmat[cl2,cl1]=prediction_vals[2]
    }
  }
  return(clmat)
}
