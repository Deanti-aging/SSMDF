SSMDC_20220928 = function()
{
  install.packages("dplyr")
  install.packages("ggplot2")
  install.packages("tidytext")
  install.packages("pacman") 
  install.packages("ggthemes")
  install.packages("readr")
  install.packages("showtext")
  install.packages("patchwork")
  install.packages("cluster")
  install.packages("kmed")
  install.packages("conclust") 
  install.packages("clv") 
  install.packages("SSLR")
  install.packages("dbscan")
  install.packages("plot3D")
  install.packages("wesanderson")
  library(wesanderson)
  library(plot3D)
  library(dplyr)
  library(cluster)
  library(ggplot2)
  library(tidytext)
  library(patchwork)
  library(kmed)
  library(conclust)
  library(clv)
  library(SSLR)
  library(pacman)
  library(ggthemes)
  library(readr)
  library(showtext)
  library(dbscan)
  datalist = data_initialize()
  tround = 10
  dataset_name = Link_num = mustLink_num = cantLink_num = clustering_method = ARI = sd = minptsd = minptsh = eps = c()
  
  for(i in c(1:10)){
    dataset = datalist[[i]]
    if(length(dataset$idcat)==0 && length(dataset$idbin)==0){
      data1 = data2 = dataset$data
    }else{
      data1 = dataset$data[, dataset$idnum]
      data2 = dataset$data
    }
    
    datalabel = dataset$datalabel
    ncluster = length(unique(dataset$datalabel))
    iterate = 100
    idnum = dataset$idnum
    idbin = dataset$idbin
    idcat = dataset$idcat
    consnum_list = round(seq(from = 50, to =1000, length.out = 5)) 
    comb = t(combn(c(1:length(datalabel)),2))
    
    for (j in c(1:length(consnum_list))){
      constrainedLink_num = consnum_list[j]
      resfastkmed = reshdbscan = resdbscan = resfastkmedo = resccls = reslcvqe = resckmeans = c()
      
      for (p in c(1:tround)) {
        if(constrainedLink_num>0){
          set.seed(p)
          labelpair_idx = sample(c(1:nrow(comb)), size = constrainedLink_num, replace = FALSE)
          allink = comb[labelpair_idx,]
          mustLink = NULL
          cantLink = NULL
          for(k in c(1:nrow(allink))){
            if(datalabel[allink[k,1]]==datalabel[allink[k,2]]){mustLink = rbind(mustLink, allink[k,])}
            if(datalabel[allink[k,1]]!=datalabel[allink[k,2]]){cantLink = rbind(cantLink, allink[k,])}
          }
        }
        
        resckmeans = tryCatch(
          {
            ckmeans_result = ckmeans(data1, k=ncluster, mustLink, cantLink, maxIter = 100)
            resckmeans = c(resckmeans, indicator(datalabel, ckmeans_result))
          },
          error = function(e){-0.1}
        )
        cat("Processed", constrainedLink_num,"COP-kmeans", " .\n")
        
        ccls_result = ccls(data1, k=ncluster, mustLink, cantLink, maxIter=100, tabuIter = 100, tabuLength = 20)
        resccls = c(resccls, indicator(datalabel, ccls_result))
        cat("Processed", constrainedLink_num,"CCLS", " .\n")
        
        ssmdc_result = ssmdc(data2, datalabel, ncluster, iterate, idnum, idbin, idcat, mustLink, cantLink)
        resfastkmed = c(resfastkmed,ssmdc_result$SSMDK)
        resdbscan = c(resdbscan,ssmdc_result$SSMDD)
        resfastkmedo = c(resfastkmedo,ssmdc_result$SSMDKO)
        reshdbscan = c(reshdbscan,ssmdc_result$SSMDH)
        minptsd = c(minptsd,ssmdc_result$minptsd)
        minptsh = c(minptsh,ssmdc_result$minptsh)
        eps = c(eps,ssmdc_result$eps)
        cat("Processed", constrainedLink_num,"SSMDC in round", p, " .\n")
        
        
        if (length(idnum)==1){
          reslcvqe = -0.1
        }else{
          lcvqe_result = lcvqe(data1, k=ncluster, mustLink, cantLink, maxIter=iterate)
          reslcvqe = c(reslcvqe, indicator(datalabel, lcvqe_result))
        }
        cat("Processed", constrainedLink_num,"LCVQE", " .\n")
        
      }
      dataset_name = c(dataset_name, names(datalist[i]))
      Link_num = c(Link_num, constrainedLink_num)
      mustLink_num = c(mustLink_num, nrow(mustLink))
      cantLink_num = c(cantLink_num, nrow(cantLink))
      clustering_method = c(clustering_method, "COP-kmeans")
      ARI = c(ARI, mean(resckmeans))
      sd = c(sd, sd(resckmeans))
      
      dataset_name = c(dataset_name, names(datalist[i]))
      Link_num = c(Link_num, constrainedLink_num)
      mustLink_num = c(mustLink_num, nrow(mustLink))
      cantLink_num = c(cantLink_num, nrow(cantLink))
      clustering_method = c(clustering_method, "SSMDD")
      ARI = c(ARI, mean(resdbscan))
      sd = c(sd, sd(resdbscan))
      
      dataset_name = c(dataset_name, names(datalist[i]))
      Link_num = c(Link_num, constrainedLink_num)
      mustLink_num = c(mustLink_num, nrow(mustLink))
      cantLink_num = c(cantLink_num, nrow(cantLink))
      clustering_method = c(clustering_method, "SSMDK")
      ARI = c(ARI, mean(resfastkmed))
      sd = c(sd, sd(resfastkmed))
      
      dataset_name = c(dataset_name, names(datalist[i]))
      Link_num = c(Link_num, constrainedLink_num)
      mustLink_num = c(mustLink_num, nrow(mustLink))
      cantLink_num = c(cantLink_num, nrow(cantLink))
      clustering_method = c(clustering_method, "SSMDH")
      ARI = c(ARI, mean(reshdbscan))
      sd = c(sd, sd(reshdbscan))
      
      dataset_name = c(dataset_name, names(datalist[i]))
      Link_num = c(Link_num, constrainedLink_num)
      mustLink_num = c(mustLink_num, nrow(mustLink))
      cantLink_num = c(cantLink_num, nrow(cantLink))
      clustering_method = c(clustering_method, "SSMDKO")
      ARI = c(ARI, mean(resfastkmedo))
      sd = c(sd, sd(resfastkmedo))
      
      dataset_name = c(dataset_name, names(datalist[i]))
      Link_num = c(Link_num, constrainedLink_num)
      mustLink_num = c(mustLink_num, nrow(mustLink))
      cantLink_num = c(cantLink_num, nrow(cantLink))
      clustering_method = c(clustering_method, "CCLS")
      ARI = c(ARI, mean(resccls))
      sd = c(sd, sd(resccls))
      
      dataset_name = c(dataset_name, names(datalist[i]))
      Link_num = c(Link_num, constrainedLink_num)
      mustLink_num = c(mustLink_num, nrow(mustLink))
      cantLink_num = c(cantLink_num, nrow(cantLink))
      clustering_method = c(clustering_method, "LCVQE")
      ARI = c(ARI, mean(reslcvqe) )
      sd = c(sd, sd(reslcvqe))
      
    }
    #picture
    out = data.frame(
      dataset_name,
      Link_num,
      mustLink_num,
      cantLink_num,
      clustering_method,
      ARI,
      sd
    )
    real_out = list(
      dataset_name,
      Link_num,
      mustLink_num,
      cantLink_num,
      clustering_method,
      ARI,
      sd,
      minptsh,
      minptsd,
      eps
    )
    
    out_now = out[which(out$dataset_name==names(datalist[i])),]
    p1 = ggplot(out_now, aes(x=Link_num, y=ARI,linetype = clustering_method,  color= clustering_method, group = clustering_method, shape = clustering_method))#
    p1 = p1 + geom_errorbar(aes(ymin=ARI-sd, ymax=ARI+sd), width = 6, size = 1, position = position_dodge(100)) 
    p1 = p1 + geom_line(size = 1.1)
    p1 = p1 + geom_point(size = 3)
    p1 = p1 + labs(x="Constraints Number",
                   y = "Rand Index",
                   title = names(datalist[i]))
    p1 = p1 + theme_bw() 
    p1 = p1 + theme(legend.position = "bottom",
                    legend.direction = "horizontal",
                    legend.title = element_blank(),
                    axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 45),
                    axis.text.y = element_text(size = 12, vjust = 0.5, hjust = 0.5, angle = 45),
                    panel.grid = element_blank(),
                    plot.title = element_text(hjust = 0.5)
    )
    p1 = p1 + scale_x_continuous(breaks = consnum_list)
    print(p1)
    ggsave(p1, filename = paste(names(datalist[i])," Dateset.PNG"))
  }
  scriptpath = rstudioapi::getSourceEditorContext()$path
  scriptpathn = nchar(scriptpath)
  suppath = substr(scriptpath, 1, scriptpathn-8)
  save(real_out, file=paste(suppath,names(datalist[i])))
  return(real_out)
}

indicator = function(datalabel, resultlabel){
  datalabel = as.integer(datalabel+1)
  resultlabel = as.integer(resultlabel+1)
  std = std.ext(datalabel, resultlabel)
  res = clv.Rand(std)
  return(res)
}

ssmdc = function(data, datalabel, ncluster, iterate = 100, idnum = NULL, idbin = NULL, idcat = NULL, mustLink = NULL, cantLink = NULL){
  #get density peak
  if(is.null(idbin) && is.null(idcat)){
    distdata = as.matrix(dist(data, method = "euclidean"))
  } else {
    distdata = distmix(data, method = "gower", idnum = idnum, idbin = idbin, idcat = idcat)
  }
  densitymaxlist = rowSums(distdata)
  densitypoint = c()
  for (i in c(1:ncluster)) {
    densitypoint = c(densitypoint, which(densitymaxlist == max(densitymaxlist))[1])
    densitymaxlist[which(densitymaxlist == max(densitymaxlist))[1]] = min(densitymaxlist)
  }
  
  #pro part
  n = nrow(distdata)
  m = nrow(mustLink)
  c = nrow(cantLink)
  promatrix = matrix(0,n,n)
  for (i in c(1:m)) {promatrix[mustLink[i,1], mustLink[i,2]] = promatrix[mustLink[i,2], mustLink[i,1]] = 1}
  for (i in c(1:c)) {promatrix[cantLink[i,1], cantLink[i,2]] = promatrix[cantLink[i,2], cantLink[i,1]] = 2}
  mustpoint = unique(c(mustLink))
  mustLink = matrix(,0,2)
  cantLink = matrix(,0,2)
  while (length(mustpoint) > 0) {
    pro = used = subi = c(mustpoint[1]) 
    mustpoint = mustpoint[-1]
    cant = c()
    while (length(pro) > 0) {
      for (i in c(1:n)) {
        if (promatrix[pro[1],i]==1 && !i%in%used) {
          pro = c(pro,i); used = c(used,i); subi = c(subi,i); mustpoint=mustpoint[-which(mustpoint==i)]
        }else if (promatrix[pro[1],i]==2) {
          cant = c(cant,i)
        }
      }
      pro = pro[-1]
      if (length(pro)==0 && length(subi)>0){
        tmustLink = t(combn(subi, 2))
        mustLink = rbind(mustLink, tmustLink)
      }
      if (length(pro)==0 && length(subi)>0 && length(cant)>0){
        tcantLink = matrix(unlist(expand.grid(subi,cant)), ncol=2)
        cantLink = rbind(cantLink, tcantLink)
      }
    }
  }
  
  
  #constrained part
  m = nrow(mustLink)
  c = nrow(cantLink)
  distdata_old = distdata_new = distdata
  maxdist = distdata_new[which(distdata==0 ,arr.ind = TRUE)] = max(distdata_new)
  mindist = min(distdata_new)
  
  if(m>0){
    for(i in c(1:m)){
      distdata[mustLink[i,2], mustLink[i,1]] = distdata[mustLink[i,1], mustLink[i,2]] = mindist/n
    }
  }
  
  if(c>0){
    for(i in c(1:c)){
      distdata[cantLink[i,2], cantLink[i,1]] = distdata[cantLink[i,1], cantLink[i,2]] = maxdist*n
    }
  }
  
  minptslist = seq(from = 2, to = n/2)
  maxresult1 = -1
  for (i in minptslist){
    hdbscanresult = indicator(datalabel, hdbscan(distdata, i)$cluster)
    if(hdbscanresult > maxresult1){
      maxresult1 = hdbscanresult
      minptsh = i
    }
  }
  
  maxvalue = max(distdata) 
  distdata0 = distdata
  distdata0[which(distdata0==0)] = maxvalue
  minvalue = min(distdata0)
  minptslist = seq(from = 2, to = n/2)
  epslist = seq(from = minvalue, to = maxvalue, length.out = n)
  maxresult2 = -1
  for (i in minptslist){
    for (j in epslist){
      dbscanresult = indicator(datalabel, dbscan(distdata, j, i, weights = NULL, borderPoints = TRUE)$cluster)
      if(dbscanresult > maxresult2){
        maxresult2 = dbscanresult
        minptsd = i
        eps = j
      }
    }
  }
  resultlabel = list()
  resultlabel$SSMDK = tryCatch(
    {
      resultlabel$SSMDK = indicator(datalabel, fastkmed(distdata, ncluster, iterate, init = densitypoint )$cluster)#
    },
    error = function(e){-0.1}
  )
  resultlabel$eps = eps
  resultlabel$SSMDH = maxresult1
  resultlabel$SSMDD = maxresult2
  resultlabel$minptsd = minptsd
  resultlabel$minptsh = minptsh
  resultlabel$SSMDKO = tryCatch(
    {
      resultlabel$SSMDKO = indicator(datalabel, fastkmed(distdata, ncluster, iterate)$cluster)#
    },
    error = function(e){-0.1}
  )
  return(resultlabel)
}

data_initialize = function(){
  scriptpath = rstudioapi::getSourceEditorContext()$path
  scriptpathn = nchar(scriptpath)
  suppath = substr(scriptpath, 1, scriptpathn-8)
  datasetpath = paste(suppath,"dataset/",sep="")
  load(file=paste(datasetpath, "overlapping_clusters", sep="")) 
  load(file=paste(datasetpath, "noise20_100_clusters", sep=""))
  load(file=paste(datasetpath, "halfring_clusters", sep=""))
  load(file=paste(datasetpath, "ring_clusters", sep=""))
 
  data = overlapping_clusters
  data[which(data=="cluster 1", arr.ind = TRUE)] = 1
  data[which(data=="cluster 2", arr.ind = TRUE)] = 2
  data[,3] = as.numeric(data[,4])
  datalabel = data[,4]
  data = data[,1:3]
  idnum = c(1,2)
  idbin = c(3)
  idcat = c()
  OVERLAPPING = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat
  )
  
  data = noise20_100_clusters
  data[which(data=="cluster 1", arr.ind = TRUE)] = 1
  data[which(data=="cluster 2", arr.ind = TRUE)] = 2
  data[which(data=="noise", arr.ind = TRUE)] = 3
  data[,3] = as.numeric(data[,3])
  datalabel = data[,3]
  data = data[,1:2]
  idnum = c(1,2)
  idbin = c()
  idcat = c()
  NOISE = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat
  )
  
  data = ring_clusters
  data[which(data=="cluster 1", arr.ind = TRUE)] = 1
  data[which(data=="cluster 2", arr.ind = TRUE)] = 2
  data[,3] = as.numeric(data[,3])
  datalabel = data[,3]
  data = data[,1:2]
  idnum = c(1,2)
  idbin = c()
  idcat = c()
  RING = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat
  )
  
  data = halfring_clusters
  data[which(data=="cluster 1", arr.ind = TRUE)] = 1
  data[which(data=="cluster 2", arr.ind = TRUE)] = 2
  data[,3] = as.numeric(data[,3])
  datalabel = data[,3]
  data = data[,1:2]
  idnum = c(1,2)
  idbin = c()
  idcat = c()
  HALFRING = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat
  )
  
  data = iris
  data$Species = as.numeric(data$Species)
  datalabel = data$Species
  data = data[,1:4]
  idnum = c(1:4)
  idbin = NULL
  idcat = NULL
  IRIS = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat
  )
  
  data = wine
  datalabel = as.numeric(data$Wine) 
  data = data[,1:13]
  idnum = c(1:13)
  idbin = NULL
  idcat = NULL
  WINE = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat
  )
  
  fertility_Diagnosis = read.csv(paste(datasetpath, "fertility_Diagnosis.txt", sep=""), header=F)
  data = fertility_Diagnosis
  data[which(data=="N", arr.ind = TRUE)] = 1
  data[which(data=="O", arr.ind = TRUE)] = 2
  datalabel = as.numeric(data[,10])
  data = data[,1:9]
  idnum = c(1,2,9)
  idbin = c(3,4,5)
  idcat = c(6,7,8)
  FERTILITY = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat
  )
  
  zoo = read.csv(paste(datasetpath, "zoo.data", sep=""), header=F)
  data = zoo[,2:18]
  datalabel = as.numeric(data[,17])
  data = data[,1:16]
  idnum = c(13)
  idbin = c(1:12,14:16)
  idcat = NULL
  ZOO = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat
  )
  
  tae = read.csv(paste(datasetpath, "tae.data", sep=""), header=F)
  data = tae
  datalabel = as.numeric(data[,6])
  data = data[,1:5]
  idnum = c(5)
  idbin = c(1,4)
  idcat = c(2,3)
  TAE = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat
  )
 
  flag = read.csv(paste(datasetpath, "flag.data", sep=""), header=F)
  data = flag[,c(2:ncol(flag))]
  type18 = unique(flag[,18])
  type29 = unique(flag[,29])
  type30 = unique(flag[,30])
  for(i in c(1:length(type18))){
    data[which(data==type18[i], arr.ind = TRUE)] = i
  }
  for(i in c(1:length(type29))){
    data[which(data==type29[i], arr.ind = TRUE)] = i
  }
  for(i in c(1:length(type30))){
    data[which(data==type30[i], arr.ind = TRUE)] = i
  }
  datalabel = data[,6] + 1
  for(i in c(1:29)){data[,i] = as.numeric(data[,i])}
  idnum = c(3,4,7,8,9,18,19,20,21,22,28,29)
  idbin = c(10,11,12,13,14,15,16,23,24,25,26,27)
  idcat = c(1,2,5,17)
  FLAG = list(
    data = data,
    datalabel = datalabel,
    idnum = idnum,
    idbin = idbin,
    idcat = idcat
  )
  
  datalist = list(
    NOISE=NOISE,
    OVERLAPPING=OVERLAPPING,
    RING=RING,
    HALFRING=HALFRING,
    
    WINE=WINE,
    IRIS=IRIS,
 
    FLAG=FLAG,
    ZOO=ZOO,
    TAE=TAE,
    FERTILITY=FERTILITY
  )
  return(datalist)
}
