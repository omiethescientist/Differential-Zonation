library(DESeq2)
library(stringr)
library(foreach)

home_dir <-  "~/Documents/projects/DifferentialZonation/"
setwd(home_dir)

rm(list = ls())

# load("../Differential-Zonation/data/simulatedData.RData")
# simData

df <-  read.csv("./data/20230430_ExpSums_dl.csv", row.names = 1)
counts <- data.matrix(df)
G <- dim(counts)[1]
top500 <- order(rowSums(counts))[9500:G]

counts <- counts[top500, ]
# counts <- counts[1:10, 1:dim(counts)[2]]
mode(counts) <- 'integer'

dose_layer_sample <- as.data.frame(str_split_fixed(colnames(df), "_", 3))
colnames(dose_layer_sample) <- c("Dose", "Layer") 
dose_layer_sample
data <- list(countData = counts, group = dose_layer_sample$Dose, zone = dose_layer_sample$Layer)

countData = data[["countData"]]
group     = data[["group"]]
zone      = data[["zone"]]
sample_name=sample_name=colnames(countData)

sel       = order(group,zone)
zone      = as.numeric(zone[sel])
group     = group[sel]
countData = countData[,sel]
sample_name = sample_name[sel]
# 
# countData = countData[rowSums(countData)!=0,]

p1 <- zone
p2 <- 0.5*(3*zone^2-1)

conds  = cbind(group,p1,p2)
colnames(conds) = c("group","p1","p2")

colData <- data.frame(row.names=colnames(countData), conds)
N=length(unique(group))

nbt = function(x){
  l=length(which(x) == TRUE)
  l
}

create_matrix_list = function(t, conds, n.co){
  require(combinat)
  
  my_matrix = list()
  
  c <- t
  s <- 0.5*(3*t^2-1)
  
  MAT <- cbind(rep(1,length(t)),c[1:length(t)],s[1:length(t)])
  GMAT <- matrix(NA,ncol=3*n.co, nrow =length(t))
  rownames(GMAT) <- conds
  colnames(GMAT) <- c(paste(c('u','a','b'),rep(1:n.co,each =3), sep = "."))
  
  it <- 1
  for(i in unique(rownames(GMAT))){
    GMAT[rownames(GMAT)==i,grep(paste0('.',it,'$'),colnames(GMAT))] = MAT[rownames(GMAT)==i,]
    it=it+1
  }
  
  vn = rep(F,n.co)
  for(i in 1:n.co){
    g = rep(F,n.co)
    g[1:i] = TRUE
    p = unique(combinat::permn(g))
    v = matrix(unlist(p),ncol = n.co,byrow = TRUE)
    vn = rbind(vn,v)
    
  }
  
  
  vn = vn[,rep(1:n.co,each=3)]
  vn[,seq(1,3*n.co,3)] = TRUE
  vn = data.frame(vn,row.names= NULL)
  vn[,dim(vn)[2] + 1]=(apply(vn,1,nbt)-n.co)/2
  colnames(vn) = c(paste(c('u','a','b'),rep(1:n.co,each =3), sep = "."),'nb_cycl')
  
  model = 1
  for(g in 0:n.co){
    
    
    nb_cycl =g
    com = expand.grid(rep(list(1:nb_cycl),nb_cycl))
    simply = apply(com,1,simply_it)
    poss =match(unique(simply),simply)
    com_l = com[poss,]
    pos = which(vn$nb_cycl == g)
    
    for(k in pos){
      if(g > 1){
        for(v in 1:nrow(com_l)){
          gmat = GMAT[,unlist(vn[k,-dim(vn)[2]])]
          ve = as.numeric(com_l[v,])
          id =1
          sa = ve
          while(length(ve) !=0){
            
            poc = which(sa == ve[1])
            po = which(ve ==ve[1])
            if(length(poc) !=1){
              poch =c(2*poc-1,2*poc)
              poch =poch[order(poch)]
              he = grep("[ab]",colnames(gmat))
              he = he[poch]
              pp=0
              for(z in 1:((length(he)-2)/2)){
                repl1 = which(gmat[,he[2*z+1]]!='NA')
                repl2 = which(gmat[,he[2*z+2]]!='NA')
                gmat[repl1,he[1]] = gmat[repl1,he[2*z+1]]
                gmat[repl2,he[2]] = gmat[repl2,he[2*z+2]]
                colnames(gmat)[he[1]]= paste(colnames(gmat)[he[1]],colnames(gmat)[he[2*z+1]],sep=',')
                colnames(gmat)[he[2]]= paste(colnames(gmat)[he[2]],colnames(gmat)[he[2*z+2]],sep=',')
                gmat[repl1,he[2*z+1]] =NA
                gmat[repl2,he[2*z+2]]=NA
                pp = pp+2
              }
              id = id+1
              ve = ve[-po]
            }else{
              ve = ve[-1]
            }
            
          }
          gmat[is.na(gmat)] =0
          del=which(apply(gmat,2,function(x) length(which(x == 0))) == length(t))
          if(length(del)!=0){
            gmat = gmat[,-del]
          }
          my_matrix[[model]] = gmat
          model = model + 1
        }
      }else{
        gmat = GMAT[,unlist(vn[k,-dim(vn)[2]])]
        gmat[is.na(gmat)] =0
        del =which(apply(gmat,2,function(x) length(which(x == 0))) == length(t))
        if(length(del)!=0){
          gmat = gmat[,-del]
        }
        my_matrix[[model]] = gmat
        model = model +1
      }
    }
    
  }
  
  
  
  return(my_matrix)
}

simply_it = function(x){
  a = 0
  for(i in x) {
    a= paste(a,paste(which(x == as.numeric(i)), collapse = "",sep = ""), collapse = "", sep = "")
  }
  a
}

models <-  create_matrix_list(zone,group, N)

models = lapply(models, function(l) l[,c(grep("u",colnames(l)),grep("a|b",colnames(l)))]   )

for (i in 1:length(models)){
  rownames(models[[i]]) = rownames(colData)}

# if (length(unique(batch))>1) {
#   # add the batch effect
#   model_b = as.matrix(model.matrix(~  batch),contrasts.arg=NULL)[,2:length(unique(batch)),drop=F]
#   colnames(model_b)=paste0("BATCH_",unique(batch)[-1])
#   models = lapply(models, function(l) cbind(model_b,l))
#   models = lapply(models, function(l) l[,c(grep("u",colnames(l)),grep("BATCH",colnames(l)),grep("a|b",colnames(l)))]   )
# }

countData
dds = DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = colData,  design=~1)
dds.full = DESeq2::DESeq(dds, full=models[[length(models)]], betaPrior = F, fitType = "parametric", test = "Wald", parallel =T, quiet = F)

deviances = sapply(models[-length(models)], function(m){
  dds.x = DESeq2::nbinomWaldTest(dds.full, modelMatrix= m, betaPrior = F, quiet = T)
  return(mcols(dds.x)$deviance)
})

deviances

message("computing BICW (Zonation)")


# calculate the BIC
BIC = as.data.frame(sapply(1:ncol(deviances), function(i) { deviances[,i] + log(ncol(countData)) * ncol(models[[i]] )}   ))

compute_BICW = function(x){
  x = as.numeric(x)
  BIC_min = min(x)
  test = exp(-0.5*(x-BIC_min))/sum(exp(-0.5*(x-BIC_min)))
  return(test)
}
#calculate the BICW
BICW               = t(apply(BIC,1,compute_BICW))
chosen_model      = apply(BIC,1,which.min)
chosen_model_BICW = apply(BICW,1,max)

message("fitting mean models")

create_matrix_list_mean = function(N,group){
  com = expand.grid(rep(list(1:N),N))
  simply = as.data.frame(t(apply(com,1,simply_it.2)))
  simply = do.call("paste",simply)
  poss =match(unique(simply),simply)
  com_l = com[poss,]
  names(com_l)=unique(group)
  com_l=com_l[order(apply(com_l,1,function(x) length(unique(x))),apply(com_l,1,function(x) length(which(x==max(x))))),]
  rownames(com_l)=1:nrow(com_l)
  
  com_l=com_l[,match(group,names(com_l))]
  p=list()
  for(j in 1:nrow(com_l)){
    if(j==1){
      p[[j]] = as.matrix(rep(1,ncol(com_l)))
      
    }else{
      p[[j]]= model.matrix(~0+ factor(as.numeric(com_l[j,])))
    }
  }
  p
}

simply_it.2 = function(x){
  
  a = match(x,x)
}

model_mean_cond=create_matrix_list_mean(N,group)

annotate_matrix = function(m,group){
  if(ncol(m)==1){
    colnames(m)=paste("u",1:length(unique(group)),sep=".",collapse=".")
  }else{
    pos_ind= match(unique(group),group)
    m=as.matrix(m)
    l=list()
    for(k in 1:ncol(m)){
      l[[k]]=as.numeric(which(m[pos_ind,k]==1))
    }
    colnames(m)=sapply(l,function(x) paste("u",x,sep=".",collapse="."))
  }
  m
}
model_mean_cond=lapply(model_mean_cond,annotate_matrix,group)

for (i in 1:length(model_mean_cond)){
  rownames(model_mean_cond[[i]]) = rownames(colData)}

for (i in 1:length(model_mean_cond)){
  rownames(model_mean_cond[[i]]) = rownames(colData)}

DDS_dev =  foreach (i = 1:length(models)) %dopar% {
  sel = which(chosen_model==i)
  gene = rownames(dds.full)[sel]
  
  if(length(gene)>0){
    M=models[[i]]
    #build the gene specific model from the rhythmic point of view
    gene_specific_mean_models = lapply(model_mean_cond,
                                       function(x) cbind(x,M[,-grep("u",colnames(M))]))
    
    dev <- lapply(gene_specific_mean_models,function(m){
      dds.m <- dds.full # Copying the full model
      dds.m <- DESeq2::nbinomWaldTest(dds.m[gene], modelMatrix= as.matrix(m), betaPrior = F) # Re-run wald test
      return(list(dds.m, mcols(dds.m)$deviance)) # Returning deviances (-2 * log likelihood) // https://support.bioconductor.org/p/107472/
      
    })
  }
  
  if(length(gene)==0){dev = list (NA, NA)}
  
  return(dev)
}

deviance_mean = NULL
for (cm_r in 1:length(models)){
  
  if(!is.na(DDS_dev[[cm_r]][1])){
    deviance_mean.x            = rbind(sapply(1:length(model_mean_cond),function(x) {DDS_dev[[cm_r]][[x]][[2]]}))
    rownames(deviance_mean.x)  = rownames(DDS_dev[[cm_r]][[1]][[1]])
    deviance_mean              = rbind(deviance_mean, deviance_mean.x)}
}

deviance_mean = deviance_mean[rownames(countData),]

message("computing BICW (mean)")

# calculate the BIC
BIC_mean = as.data.frame(sapply(1:ncol(deviance_mean), function(i) { deviance_mean[,i] + log(ncol(countData)) * ncol(model_mean_cond[[i]] )}   ))

#calculate the BICW
BICW_mean = t(apply(BIC_mean,1,compute_BICW))

chosen_model_mean = apply(BIC_mean,1,which.min)

chosen_model_mean_BICW = apply(BICW_mean,1,max)

message("extracting rhythmic parameters")

parameters=NULL

compute_param = function(dds, gene, N){
  
  dds = dds[gene,]
  param = c(paste(rep(c('u','a','b'),each=N),rep(1:N,3), sep = "."))
  
  paramout = rep(NA,N*3)
  
  for(i in 1:N){
    
    u=coef(dds)[grep(paste(param[i],"Intercept",sep="|"), colnames(coef(dds)))]
    a=coef(dds)[grep(param[i+N], colnames(coef(dds)))]
    b=coef(dds)[grep(param[i+N*2], colnames(coef(dds)))]
    
    if(length(u) ==0) u=NA
    if(length(a) ==0) a=NA
    if(length(b) ==0) b=NA
    
    # phase=period/(2*pi)*atan2(b,a)
    # amp =2*sqrt(a^2+b^2)
    # relamp=0.5*amp/u
    # if(!is.na(phase)){
    #   #if(phase<0) phase=phase+period
    #   #if(phase>period) phase=phase-period
    #   phase=phase%%period
    # }
    paramout[(1:3 + 3*(i-1))] = c(u,a,b)
  }
  
  #names(paramout) = c(paste(c('mean','a','b','amp','relamp','phase'),rep(1:N,each =6), sep = "_"))
  paramout
}

parameters =  foreach (i = 1:nrow(deviance_mean)) %dopar% {
  gene = rownames(deviance_mean)[i]
  cm_r = chosen_model[i]
  cm_m = chosen_model_mean[i]
  dds= DDS_dev[[cm_r]][[cm_m]][[1]]
  out = compute_param(dds, gene ,N)
  return(data.frame(row.names= gene, t(matrix(out)))           )
}

parameters            = data.frame(do.call(rbind.data.frame, parameters))
colnames(parameters)  = c(paste(c('mean','a','b'),rep(unique(group),each =3), sep = "_"))
parameters            = parameters[rownames(countData),]

# Generate all the count and expression data
# raw counts
counts_RF        =  counts(dds.full, normalized = FALSE)

#vst stabilized counts
vsd <- DESeq2::varianceStabilizingTransformation(dds.full)
vsd <- assay(vsd)

#normalized counts
ncounts_RF       = counts(dds.full, normalized = TRUE)

# generate a table summarizing the analysis
complete_parameters = cbind(parameters,chosen_model,chosen_model_BICW, chosen_model_mean, chosen_model_mean_BICW)
global_table = merge(ncounts_RF,complete_parameters, by="row.names")
rownames(global_table) = global_table$Row.names
global_table_df  = global_table[,-grep("Row.names",colnames(global_table))]

global_table_df = global_table_df[rownames(countData),]

out = list()

out[["zone"]]        = zone
out[["group"]]       = group
out[["results"]]     = global_table_df
out[["BICW_zonation"]] = BICW
out[["BICW_mean"]]   = BICW_mean
out[["vsd"]]         = vsd
out[["ncounts"]]     = ncounts_RF
out[["counts"]]      = counts_RF
out[["parameters"]]  = complete_parameters
out[["cook"]]        = assays(dds.full)[["cooks"]]
out[["dds.full"]]    = dds.full

message("finished!")
return(out)
