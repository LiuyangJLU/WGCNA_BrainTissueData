#原始脑部数据的处理
#step1:将原始的txt文件写入csv文件，打开文件，删除不必要的行和列
#step2:将数据框的第一列作为列名
#step3:删除第一列
#step4:mode()查看数据类型，如果是list，转换我Numeric()
#eset<-apply(eset,1,function(x){as.numeric(x)})
#使相应的行名等于列名
CRBL_rawset  = read.table("E:/R_Code/tissue/BrainTissue/expr_CRBL.txt")
CRBL_csv = write.csv(CRBL_rawset,file="E:/R_Code/tissue/BrainTissue/expr_CRBL.csv")

FCTX_rawset = read.table("E:/R_Code/tissue/BrainTissue/expr_aveALL.txt")
FCTX_csv = write.csv(FCTX_rawset,file="E:/R_Code/tissue/BrainTissue/expr_FCTX.csv")

WHMT_rawset = read.table("E:/R_Code/tissue/BrainTissue/expr_WHMT.txt")
WHMT_csv = write.csv(WHMT_rawset,file="E:/R_Code/tissue/BrainTissue/expr_WHMT.csv")

MEDU_rawset = read.table("E:/R_Code/tissue/BrainTissue/expr_MEDU.txt")
MEDU_csv = write.csv(MEDU_rawset,file="E:/R_Code/tissue/BrainTissue/expr_MEDU.csv")

HIPP_rawset = read.table("E:/R_Code/tissue/BrainTissue/expr_HIPP.txt")
HIPP_csv = write.csv(HIPP_rawset,file="E:/R_Code/tissue/BrainTissue/expr_HIPP.csv")

OCTX_rawset = read.table("E:/R_Code/tissue/BrainTissue/expr_OCTX.txt")
OCTX_csv = write.csv(OCTX_rawset,file="E:/R_Code/tissue/BrainTissue/expr_OCTX.csv")

PUTM_rawset = read.table("E:/R_Code/tissue/BrainTissue/expr_PUTM.txt")
PUTM_csv = write.csv(PUTM_rawset,file="E:/R_Code/tissue/BrainTissue/expr_PUTM.csv")

SNIG_rawset = read.table("E:/R_Code/tissue/BrainTissue/expr_SNIG.txt")
SNIG_csv = write.csv(SNIG_rawset,file="E:/R_Code/tissue/BrainTissue/expr_SNIG.csv")

TCTX_rawset = read.table("E:/R_Code/tissue/BrainTissue/expr_TCTX.txt")
TCTX_csv = write.csv(TCTX_rawset,file="E:/R_Code/tissue/BrainTissue/expr_TCTX.csv")

THAL_rawset = read.table("E:/R_Code/tissue/BrainTissue/expr_THAL.txt")
THAL_csv = write.csv(THAL_rawset,file="E:/R_Code/tissue/BrainTissue/expr_THAL.csv")

#统一读取表中txt文件，转换为csv文件并统一写入文件夹中
#step1: 首先，函数形参传入原始数据
#原始数据
file1 = "E:/R_Code/tissue/BrainTissue/expr_CRBL.txt"
file2 = "E:/R_Code/tissue/BrainTissue/expr_WHMT.txt"
file3 = "E:/R_Code/tissue/BrainTissue/expr_FCTX.txt"
file4 = "E:/R_Code/tissue/BrainTissue/expr_HIPP.txt"
file5 = "E:/R_Code/tissue/BrainTissue/expr_MEDU.txt"
file6 = "E:/R_Code/tissue/BrainTissue/expr_OCTX.txt"
file7 = "E:/R_Code/tissue/BrainTissue/expr_PUTM.txt"
file8 = "E:/R_Code/tissue/BrainTissue/expr_SNIG.txt"
file9 = "E:/R_Code/tissue/BrainTissue/expr_TCTX.txt"
file10 = "E:/R_Code/tissue/BrainTissue/expr_THAL.txt"

str1 ="CRBL"
str2 ="WHMT"
str3 = "FCTX"
str4 = "HIPP"
str5 = "MEDU"
str6 = "OCTX"
str7 = "PUTM"
str8 = "SNIG"
str9 = "TCTX"
str10 = "THAL"




Data_progress<-function(file,str)
{
  writer = read.table(paste(file,1:10,seq="_",str,1:10))
  
  filename<-as.character(x = )
  
  for(i in 10)
  {
    
     writer = read.table(rawset,file="E:\R_code\tissue\expr_'i'.txt")
     
     
  }

  filename<-str()
  eset = write.csv(rawset,file="E:\R_Code\tissue\expr_'filename'.csv")
  
}


#DataProcess
#eset = read.table("E:/R_Code/tissue/expr_aveALL.txt")

P<-function(rawset)
{
  rownames(rawset)<-rawset[,1]
  rawset<-rawset[,-1]
  eset<-rawset
  eset<-apply(eset,1,function(x){as.numeric(x)})
  rownames(eset)<-colnames(rawset)
  return(eset)
}

Rawset = read.csv("E:/R_Code/tissue/rawset.csv")
dim(rawset)#查看行和列
head(rawset[1:4,1:4])#查看前4行和4列
rownames(rawset)<-rawset[,1]#将第一列赋值给行名
rawset<-rawset[,-1]#删除多余的第一列
eset<-rawset#重新赋值eset
eset<-apply(eset,1,function(x){as.numeric(x)})#将数据类型转换为数字类型
mode(eset)
rownames(eset)<-colnames(rawset)

#Coefficient of Variation
eset1 <- eset[,order(apply(eset,2,mad), decreasing = T)[1:4000]]
eset<-eset1



library(WGCNA)
#detect sampleOutlier
sampleTree<-hclust(dist(eset),method="average")
sizeGrWindow(12,9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree,main="Sample clustering to detect outliers",sub="",xlab="",cex.lab=1.5,cex.axis=1.5,cex.main=2)



#Network Construction

softPower =12
adjacency = adjacency(eset,power = softPower)
TOM_P<-TOMsimilarity(adjacency)
dissTOM<-1-TOM_P

geneTree<-hclust(as.dist(dissTOM),method="average")
plot(geneTree,xlab="",sub="",main = "Gene clustering TOM-based dissimilarity",labels=FALSE,hang=0.04)
#这里的最小的聚类数量依赖于基因数
dynamicMods = cutreeDynamic(dendro = geneTree,distM = dissTOM,minClusterSize = 30)
dynamicColors = labels2colors(dynamicMods)



#module eigengenes defined as the first principal component of the expression of the gene within the module
MElist <- moduleEigengenes(eset,colors = dynamicColors)
MEs<-MElist$eigengenes
MEdiss<-1-cor(MEs)
METree<-hclust(as.dist(MEdiss),method = "average")

#Graphical the result
sizeGrWindow(7,6)
plot(METree,main="Clustering of module eigengenes(Normal)")
MEDissThres = 0.1
abline(h=MEDissThres,col="red")
#合并模块的目的是看主要的模块之间的差异
merge<-mergeCloseModules(eset,dynamicColors,cutHeight = MEDissThres)
#The merged module colors
mergeColors = merge$
  cellToExps<-functio





moduleDetect<-function(eset,dissTOM)
{
  geneTree<-hclust(as.dist(dissTOM),method="average")
  plot(geneTree,xlab="",sub="",main = "Gene clustering TOM-based dissimilarity",labels=FASLE,hang=0.04)
  dynamicMods = cutreeDynamic(dendro = geneTree,distM = dissTOM,minClusterSize = 30)
  dynamicColors = labels2colors(dynamicMods)
  
}


