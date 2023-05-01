library(DESeq2)

home_dir <-  "/mnt/home/kanaomar/sudin_pipeline/projects/TCDD_Timeseries/"
setwd(home_dir)

rm(list = ls())

load("../Differential-Zonation/data/simulatedData.RData")
simData

df <-  read.csv("./Data/20230430_ExpSums_dl.csv", row.names = 1)
counts <- data.matrix(df)
G <- dim(counts)[1]
top500 <- order(rowSums(counts))[500:G]

counts <- counts[top500, ]
# counts <- counts[1:10, 1:dim(counts)[2]]
mode(counts) <- 'integer'

library(stringr)

dose_layer_sample <- as.data.frame(str_split_fixed(colnames(df), "_", 3))
colnames(dose_layer_sample) <- c("Dose", "Layer", "Sample") 
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

countData = countData[rowSums(countData)!=0,]

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

