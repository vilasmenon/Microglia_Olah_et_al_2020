###functions called by main script code_for_submission.R###


require(Seurat)   ###note: Seurat version 2 was used for the clustering analysis in the accompanying paper Olah et al. (2020)
require(randomForest)
require(edgeR)

###Cluster over parameters (number of PCs and resolutions)###
cluster_over_parameters=function(dat,batchval,keepcols,pcrange=5:15,resrange=seq(from=0.2,to=0.8,by=0.2),krange=c(5,10,20,30),batchreg=1) {
  fullclusterings=list()
  if (length(keepcols)>5) {
    datobj_0 <- new("seurat", raw.data = data.frame(alldat[,keepcols]))
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


####random forest classification for each cell####
cell_by_cell_prediction=function(dataset,clusterids,crossval=4,iterations=100,filename="rf_pred.rda") {
  allf=unique(clusterids)
  predlist=list()
  for (ii in 1:(length(allf)-1)) {
    testcols1=rownames(dataset)[which(clusterids==allf[ii])]
    for (jj in (ii+1):length(allf)) {
      testcols2=rownames(dataset)[which(clusterids==allf[jj])]
      outmat1=matrix(0,nrow=100,ncol=length(testcols1))
      colnames(outmat1)=testcols1
      outmat2=matrix(0,nrow=100,ncol=length(testcols2))
      colnames(outmat2)=testcols2
      numtrain=min(round(length(testcols1)*3/4),round(length(testcols2)*3/4))
      alldone1=c()
      alldone2=c()
      set.seed(ii*iterations+jj)
      for (kk in 1:iterations) {
        sampids1=sample((1:length(testcols1))%%crossval)
        sampids2=sample((1:length(testcols2))%%crossval)
        for (mm in unique(sampids1)) {
          trainset=c(sample(testcols1[sampids1!=mm],min(numtrain,length(which(sampids1!=mm)))),sample(testcols2[sampids2!=mm],min(numtrain,length(which(sampids2!=mm)))))
          testset=c(testcols1[sampids1==mm],testcols2[sampids2==mm])
          ttt=as.factor(rep(c(allf[ii],allf[jj]),times=c(min(numtrain,length(which(sampids1!=mm))),min(numtrain,length(which(sampids2!=mm))))))
          predval=randomForest(as.matrix(dataset[trainset,]),ttt)
          outpred=predict(predval,as.matrix(dataset[testset,]))
          names(outpred)=testset
          outmat1[kk,testcols1[sampids1==mm]]=as.character(outpred[testcols1[sampids1==mm]])
          outmat2[kk,testcols2[sampids2==mm]]=as.character(outpred[testcols2[sampids2==mm]])
        }
        print(c(ii,jj,kk))
      }
      nam=paste0(allf[ii],"-",allf[jj])
      predlist[[nam]]=list()
      predlist[[nam]][[1]]=outmat1
      predlist[[nam]][[2]]=outmat2
      save(predlist,file=filename)
    }
  }
}


####Function to return differential gene expression across each pair of clusters####
pairwise_differential_expression=function(counts,clusterids) {
  allf=unique(clusterids)
  allf=allf[order(allf)]
  allfits=list()
  cpm=sweep(counts,2,colSums(counts),"/")*10^6
  for (ii in 1:length(allf)) {
    keepcols1=which(clusterids==ii)
    for (jj in setdiff(1:length(allf),ii)) {
      keepcols2=which(clusterids==jj)
      subvec=factor(c(rep(1,length(keepcols1)),rep(2,length(keepcols2))))
      e_design=model.matrix(~subvec)
      y2 = DGEList(counts=counts[,c(keepcols1,keepcols2)])
      y2 = estimateDisp(y2, e_design)
      fit = glmQLFit(y2, e_design)
      qlf.2vs1 <- glmQLFTest(fit, coef=2)
      outval=topTags(qlf.2vs1,n=nrow(col),p.value=1)
      mean1=apply(cpm[rownames(outval$table),keepcols1],1,mean)
      mean2=apply(cpm[rownames(outval$table),keepcols2],1,mean)
      frac1=rowSums(cpm[rownames(outval$table),keepcols1]>0)/length(keepcols1)
      frac2=rowSums(cpm[rownames(outval$table),keepcols2]>0)/length(keepcols2)
      print(c(ii,jj))
      outval$table=cbind(outval$table,mean1,mean2,frac1,frac2)
      allfits[[paste0(allf[ii],"-",allf[jj])]]=outval
    }
  }
  return(allfits)
}

####Function to return differential gene expression for each cluster when compared to all other clusters####
one_versus_all_differential_expression=function(counts,clusterids) {
  allf=unique(clusterids)
  allf=allf[order(allf)]
  allfits=list()
  cpm=sweep(counts,2,colSums(counts),"/")*10^6
  for (ii in 1:length(allf)) {
    keepcols1=which(clusterids==ii)
    keepcols2=which(clusterids!=ii)
    subvec=factor(c(rep(1,length(keepcols1)),rep(2,length(keepcols2))))
    e_design=model.matrix(~subvec)
    y2 = DGEList(counts=counts[,c(keepcols1,keepcols2)])
    y2 = estimateDisp(y2, e_design)
    fit = glmQLFit(y2, e_design)
    qlf.2vs1 <- glmQLFTest(fit, coef=2)
    outval=topTags(qlf.2vs1,n=nrow(col),p.value=1)
    mean1=apply(cpm[rownames(outval$table),keepcols1],1,mean)
    mean2=apply(cpm[rownames(outval$table),keepcols2],1,mean)
    frac1=rowSums(cpm[rownames(outval$table),keepcols1]>0)/length(keepcols1)
    frac2=rowSums(cpm[rownames(outval$table),keepcols2]>0)/length(keepcols2)
    outval$table=cbind(outval$table,mean1,mean2,frac1,frac2)
    allfits[[paste0(allf[ii],"-all")]]=outval
  }
  return(allfits)
}
