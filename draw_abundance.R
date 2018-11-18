# !/usr/bin/env Rscript
# 输入样本相对丰度表，以及meta文件
# 功能1: 生成样本之间的相关性矩阵即热图（基于皮尔逊相关系数）
# 功能2: 为每一个样本计算top20的物种丰富度并绘制bar图
# 功能3: 将所有样本集中绘制物种丰度组成的堆叠图

options(warn=-1)
# 依赖关系检查，安装和加载
package_list <- c("optparse","reshape2","ggplot2","ggpubr","RColorBrewer","gplots","dplyr")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p,character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    suppressWarnings(suppressMessages(library(p, character.only=TRUE, quietly=TRUE, warn.comflicts=FALSE)))
  }}

# 清理工作环境 clean environment object
rm(list=ls())
# 加载依赖关系
library('optparse')
library('reshape2')
library('ggplot2')
library('dplyr')
library('ggpubr')
library("gplots")
library("RColorBrewer")

# 解析命令行
if(TRUE){
  option_list <- list(
    make_option(c("-i","--input"),type="character",default="genus_relative_abundance.csv",help="input abundance file to read[default %default]"),
    make_option(c("-m","--meta"),type="character", default="meta.csv", help="meta file[default %defailt]"),
    make_option(c("-t","--filtered_thread"), type="numeric", default=0.01, help="绘制物种丰富度堆叠图时候过滤菌的阈值[default %defailt]")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
}

# 读取abundance文件以及meta文件
data <- read.csv(opts$input, header = T)
meta <- read.csv(opts$meta, header=T, sep=',')
# 删除data表中的第一列
data <- data[2:ncol(data)]
# 将run_id作为meta的行名
meta2 = meta[2:ncol(meta)]
rownames(meta2) <- meta[,1]
# 对meta2表按照data表过滤，有时候有可能meta表中具有data中所没有的run_id
# 此时SRR*顺序也变为一致
meta2 <- meta2[colnames(data)[2:ncol(data)],]
# 将NA值补成0
data[is.na(data)]=0

# 修改data表头名字，相当于将样本后加上疾病状态
for(i in 2:ncol(data)){
  colnames(data)[i]=paste(meta2$sample_name[i-1],meta2$host_age[i-1],sep="_")
  }

# ------------------------------功能模块1------------------------------
# 计算所有样品之间的相关性矩阵，这里使用的是皮尔逊相关系数
print("执行功能1:...")
sim = cor(data[2:ncol(data)], method = "pearson")
# 将缺失值变为0
sim[is.na(sim)]=0

# 使用热图可视化,并保存
# 绘制样本相关性热图
pdf(file=paste("heat_cor_samples.pdf",sep=""))
heatmap.2(sim, Rowv=TRUE, Colv=TRUE, dendrogram='both', trace='none', srtCol=55 , 
          margins=c(6,6), col=rev(colorRampPalette(brewer.pal(11, "RdYlGn"))(256)),density.info="none") 
dev.off()

# -------------------------功能模块2----------------------
# 功能2: 为每一个样本计算top20的物种丰富度并绘制bar图
print("执行功能2...")
for(i in 2:ncol(data)){
  # 筛选出最多的20个丰度
  df <- data[order(desc(data[i])),][1:20,][,c(1,i)]
  # 对第一列改名
  colnames(df)[1]="Genus"
  ggplot(df,aes(reorder(Genus,df[,colnames(df)[2]]),df[,colnames(df)[2]],
                  label=paste(round(df[,colnames(df)[2]]*100,2),"%",sep="")))+
    geom_bar(stat="identity")+
    geom_text(stat="identity",hjust=0, nudge_y=0.01) +
    coord_flip() + xlab("Genus")+ylab("Relative abundance(%)") + 
    labs(title=paste(colnames(df)[2],"genus relative abundance",sep="_"), size=0.5) + 
    ylim(0, max(df[2])*1.2)+
    theme(plot.title = element_text(hjust = 0.5,size=10))
    ggsave(file=paste(paste(colnames(df)[2],".pdf",sep=''),sep=""),width=5, height=3)
}


#-------------------功能模块3---------------------
# h绘制物种丰富度堆叠图
print("执行功能3...")
n = nrow(data) # 原先具有的菌数量
print(paste("原来具有的菌种数量为:",n,sep=""))
# 将小于0.01的菌过滤掉，所有样本的丰富度总和小于0.01的菌
# 过滤后剩余菌数量为
data2 <- data[rowSums(data[2:ncol(data)])>opts$filtered_thread,]
filtered_n <- nrow(data2)
print(paste("过滤后剩余菌种类数量为:",filtered_n,sep=""))
# 将data2按照meta中的disease顺序排序的样本顺序一致的排序
data3 <- data2[2:ncol(data2)]
data3 <- data3[,order(meta2$disease)]
data3['Genus'] <- data2[1]


data4 <- melt(data3)
# 改名
data4$variable <- as.factor(data4$variable)
# 绘图,所有样本的物种组成堆叠图，并隐藏图例
ggplot(data4,aes(variable, value, fill=factor(Genus))) + 
  geom_bar(position='stack',stat="identity")+
  # scale_fill_discrete(guide=FALSE)+
  coord_flip() +
  xlab("genus")+
  ylab("Samples")
ggsave("samples_genus_abundance.pdf")

