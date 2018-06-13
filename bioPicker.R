library(caret)
library()

filedir = './data/brca_data.csv'

eset<-read.csv(filedir,header = TRUE)


#将csv文件读取并处理后返回结果,行样本，列基因名
csvToEset<-function(fileDir){
  eset<-read.csv(fileDir,header = TRUE)
  tryCatch(
    {
      rownames(eset)<-eset[,1]
      eset<-eset[,2:ncol(eset)]
    },
    error=function(e){cat("fail to read goal csv file",conditionMessage(e),"\n\n")},
    finally={}
  )
  return(eset)
}

eset = csvToEset('./data/brca_data.csv')




#remove the outliers and return the remain samples and labels
removeOutliers<-function(eset,label){
  sampleTree<-hclust(dist(eset),method="average")
  sizeGrWindow(12,9)
  par(cex=0.6)
  par(mar=c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
  #get user input
  again=TRUE
  while(again){
    tryCatch(
      { print("please input a height(numeric) to cut the tree")
        height<-scan("",what=numeric(0),nlines=1)
        abline(h=height,col="red")
        print("please input the min size to recognize a module")
        minSize<-scan("",what=integer(0),nlines=1)
        again=FALSE},
      error=function(e){print("there are some error，please input again...");again=TRUE},
      finally={}
    )
  }
  #cut tree
  clust = cutreeStatic(sampleTree, cutHeight = height, minSize = minSize)
  rm(sampleTree)
  again<-TRUE
  removeSample<-c(1:nrow(eset)) #init
  while(again){
    tryCatch(
      { print("please input the index of cluster to delete, press enter to complete:")
        removeClust<-scan("",what=integer(0),nlines = length(table(clust)))
        ifelse(
          #there is a bug to fix
          all(removeClust%in%c(0:length(table(clust))))==TRUE&length(removeClust)!=0,
          #表达式拼???
          {eps<-sapply(removeClust,function(x){paste("clust==",x,sep="")})
          eps<-parse(text=paste(eps,"",sep = "",collapse = "|"))
          removeSample<-eval(eps)
          rm(eps)
          again=FALSE},
          {print("some index is valid or empty input,please input again..")
            again=TRUE}
        )
      },
      error=function(e){print("输入有误，请重新输入...");again=TRUE},
      finally={}
    )
    
  }
  eset2<-eset[-which(removeSample),]
  label2<-label[-which(removeSample)]
  rm(eset,removeSample,again,clust)#移除删除后的异常点和样本
  e_l<-list(eset=eset2,label=label2)
  return(e_l)
}



#a function to detect modules
moduleDetect<-function(eset,dissTOM,MEDissThres=0.25){
  #基于WGCNA算法中的TOM矩阵，利用相异性矩阵作为层次聚类的距离进行聚类。
  
  geneTree <- hclust(as.dist(dissTOM), method = "average")#聚类，通过1-TOM作为聚类的距离
  minSize<-floor(dim(eset)[2]/100)#
  if(minSize<4)sessionInfo()
    minSize=4
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minSize)
  dynamicColors = labels2colors(dynamicMods)#将对应的模块标记颜色
  # MEDissThres<-0.25
  # Call an automatic merging function
  merge = mergeCloseModules(eset, dynamicColors, cutHeight = MEDissThres, verbose = 0)
  # The merged module colors
  mergedColors = merge$colors;
  # Rename to moduleColors
  moduleColors = mergedColors
  return(moduleColors)
}



#十折交叉验证
testbioPicker<-function(eset,label,model="SVM",cor1=0.85,k=1,MEDissThres=0.25,type="acc"){
  # set.seed(100)
  result<-sapply(1:10,function(x){
    #十折交叉
    folds<-createFolds(y=label,k=5)#创建5折交叉验证
    result1<-sapply(folds,function(x){
      #每次先选好训练集和测试集
      trainset<-eset[-x,]
      testset<-eset[x,]
      trainlabel<-label[-x]
      testlabel<-label[x]
      re<-wgcnaPredict(trainset,trainlabel,cor1=cor1,k=k,MEDissThres=MEDissThres)#返回的预测集标签
      result2<-trainModel(testset[,re],testlabel,type=type)$result
      result2_2<-trainModel3(testset[,re],testlabel,k=k,type=type)
      result2_3<-trainModel3(testset[,re],testlabel,k=k)
      print(paste(result2,"---------",result2_2,"-------------",result2_3,"--------"))
      result2
    })
    mean(result1)
  })
  return(mean(result))
}



#eset行样本，列特征
wgcnaPredict<-function(eset,label,stop_acc=1,cor1=0.85,k=1,MEDissThres=0.25){
  eset<-prepareData(eset,label,cor1)
  eset2<-scale(eset)#eset2是标准化的eset,仅用于聚类
  dissTOM<-1-cor(eset2)#相似矩阵化为相异矩阵，用于层次聚类
  moduleColors<-moduleDetect(eset2,dissTOM,MEDissThres)#获取簇
  colors<-table(moduleColors)
  removeColors<-names(which(colors==1)) #移除只有一个元素的簇，否则会出错
  colors<-setdiff(names(colors),removeColors) #剩下模块的名字
  #和colors顺序一致
  colors<-lapply(colors,function(x){x1<-as.data.frame(eset[,which(moduleColors==x)])})
  #获取各个模块和label的cor的降序基因名,命名为color_dec
  #顺序与colors一致
  
  colors_dec<-lapply(colors,function(x){
    cor1<-cor(x,label)
    cor1<-t(cor1)
    cor1<-as.data.frame(cor1)
    cor1<-sort(abs(cor1),decreasing = T)
    cor1<-colnames(cor1)
  })
  rm(colors,eset2,dissTOM,removeColors)
  gc()
  #将各个模块cor第一名放进去，颜色顺序与colors同
  first<-sapply(colors_dec,function(x){x[1]})
  first<-unlist(first)
  
  #迭代替换基因
  first<-replaceGene(first,colors_dec,eset,label,stop_acc,k=k)

  #print("Starting remove least contribution feature ... ")
  #result2<-removeWF(eset[,first],label,k=k)
  #用于测试目标基因在模块中的位置
  #pos<-showCorPos(eset,moduleColors,result2$genelist,label)
  #li <- list(result=result2,pos=pos)
  #
  print("Start adding gene..")
  first = addGene(first,eset,label,end=stop_acc,k=k)
  return(first)
}

#replace geneVector genes with gene in coresponding module,to reach the highest predict score
#first 由各个colors_dec第一个元素组成的基因列表
#colors_dec 某个模块基因降序排列（与label的cor）
#fast

replaceGene<-function(first,colors_dec,eset,label,end=1,model="NB",k=1){
  cl.cores <- detectCores()#
  cl <- makeCluster(cl.cores)
  clusterEvalQ(cl,library(caret))
  clusterEvalQ(cl,library(nnet))
  clusterEvalQ(cl,library(e1071))
  clusterEvalQ(cl,library(class))
  clusterExport(cl,c("trainModel3","removeWF","replaceGene","getAcc"))
  #
  genelist<-as.character(first)
  acc<-trainModel3(eset[,genelist],label,k)
  # index<-sample(c(1:length(genelist)))
  # index<-c(14,3,1,12,6,9,13,8,5,7,10,11,2,4,15)
  # print(index)
  for(i in 1:length(genelist)){
    # for(i in index){
    if(length(colors_dec[[i]])==1){next}
    #if the max acc bigger than end,stop the loop
    if(acc>=end)
      break
    accs<-parSapply(cl,colors_dec[[i]][-1],function(x){
      genelist[i]<-x
      acc2<-trainModel3(eset[,genelist],label,k)
    })
    #acc提升，替换
    if(max(accs)>acc){
      max_index<-1+which(accs==max(accs))[1]
      genelist[i]<-colors_dec[[i]][max_index]
      acc<-max(accs)
    }
  }
  stopCluster(cl)
  first<-genelist
  return(first)
}

addGene<-function(first,eset,label,end=1,model="NB",k=1){
  cl.cores <- detectCores()
  cl <- makeCluster(cl.cores)
  clusterEvalQ(cl,library(caret))
  clusterEvalQ(cl,library(nnet))
  clusterEvalQ(cl,library(e1071))
  clusterEvalQ(cl,library(class))
  clusterExport(cl,c("trainModel3","removeWF","replaceGene","getAcc"))
  genelist<-as.character(first)
  acc<-trainModel3(eset[,genelist],label,k)
  allgene<-colnames(eset)
  continue = T
  while(continue){
    backup<-setdiff(allgene,genelist)
    accs<-parSapply(cl,backup,function(x){
      genelist2<-c(genelist,x)
      acc2<-trainModel3(eset[,genelist2],label,k)
    })
    #acc提升，替换
    if(max(accs)>acc){
      max_index<-which(accs==max(accs))[1]
      genelist<-c(genelist,backup[max_index])
      acc<-max(accs)
      acc<-max(accs)
      print("-----------------ti huan ----------------")
      print(acc)
    }else{
      continue = F
    }
    if(acc>=end){
      continue = F
    }
    # if(acc==max(accs))
    #   continue = F
    
  }
  stopCluster(cl)
  return(genelist)
}



