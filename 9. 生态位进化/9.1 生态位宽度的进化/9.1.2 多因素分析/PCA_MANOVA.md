#### 9.1.2 多维度因子分析

```R
说明:
## 进行多维度分析,可以从两个多维度出发,一个是多近缘物种对角度,另外一个是多维环境因子的角度;
## 上述两个角度在进行分析时,可能采用的方法和分析手段是不同的;
##　一般pca/cca和t-sne以及其他线性模型都比较适合于从环境维度分析物种建模结果;而pcoa则比较适合从物种维度出发来聚类评估物种多样性特征;
```

##### 9.1.2.1 降维排序分析

###### 9.1.2.1.1PCA概念

```r
## 说明:
pca分析是一种探索式分析方法;
### 在pca分析中对主轴解释更为重要的是特征根(Eigenvalue)的选择是非常重要的；当然也会生成对应 的参数占比；
### 选择保留多少个排序周并没有统一的标准，通常很随意；也有方法如Kaiser-guttman准则，先计算所有与轴的特征根平均值，然后选择保留特征根超过平均值的轴。
### pca的中部区域是方差变化最小的区域，往往代表数据差异性较小；而pca的边缘区域则往往代表着最大的方差变异量；
### pca是一种展示样方之间的欧式距离的线性方法；因此有时在反映环境因子(往往是单峰模型)特征是不合适的,建议选择ca或者cca;

## pca在构建时是原始数据重新投影到新的坐标系,会改变数据原有位置排序,因此不能直观的反映数据在空间维度的上的特征偏好;因此对于pca的单物种或者多物种分布偏好,需要谨慎对待;
```

```r
## 关于pca的一些补充：
1、pca原始假设要求数据符合多元正态分布，但在生态学领域只要数据偏离正态不太离谱，pca对于数据是否正态并不敏感；
2、pca分析必须从相同的量纲表出发；原因是需要将变量总方差分配给特征根；
3、pca计算的直接原理是分解线性协方差矩阵的特征分解，找垂直特征根作为第一主轴；
4、矩阵数据不能倒置；
5、二元数据也可以应用于pca分析，但需要进行Hellinger转化或者弦转化对数据进行处理；
```

```R
## PCA-双序图：
## pca双序图，也称之为样方和变量双序图(bioplots of sites and variables)
## pca双序图可以直接调用vegan包里的bioplot.rda()函数；在bioplot分析时会有两种设置参数：
其一是scaling =1（关注对象）；其二是scaling =2（关注变量）；

# 第一个变量分析图反映的是物种分布在环境梯度上的变化关系；可以将物种分布数据划分为不同组的结果；然后在pca的中部的变量不具有环境特征性；如果物种分布数据主要分布在主轴中西周围则说明这些分布点对应pca1和pca2的贡献较小；

# 第二个变量分析图反映的环境不同分组之间的差异，以及这些变量之间的相关关系；如果两个变量之间呈较小的角度差异则反映两个变量具有正相关性；如果两个变量近乎垂直相交则说明两个变量几乎没有关系；如果两个变量的之间的夹角大于90度，则说明两个变量之间为负相关关系；
##  在因子图中的圆圈，称作平衡贡献圆。它的半径是根号下d/p;其中d是双序图中轴数量(通常d=2),p是pca的维度(即变量的个数)。平衡贡献圆的半径代表变量的向量长度对排序的平均贡献率。如果某个变量的箭头长度长于圆的半径，则代表它对这个排序空间的贡献大于所有变量的平均贡献；箭头的长度反映对两个主轴的贡献量；
```

###### 简化pca实现

```{r eval =F}
## 用于pca分析的数据构建:
## 简化pca实现

## 养成习惯,设置工作路径:
rm(list = ls())
setwd("D:\\zhang")
localdata <- read.csv("bjpoints.csv")
tiffs <- list.files("./preasc",pattern = "asc",full.names = TRUE)
library(raster)
ascs <- stack(tiffs)
data <- extract(ascs,localdata[,-1])
###  为画图方便,仅选择其中6种环境变量;
data <- as.data.frame(data[,-7])

## 初级pca实现:
pca2 <- princomp(data,cor = T)
pca1 <- prcomp(data,scale = T)
##查看pca提取结果:
#可以查看个主轴的特征根和贡献以及各因子在因子的负载(相关关系)；
summary(pca2,loadings =T)

##简单绘图:
## 注意使用的biplot();不是bioplot()
## black为点的颜色,red为各项主轴的颜色;
biplot(pca2, col = c("black", "red"))

## 获取降维后的数据矩阵:
pca2_data <- predict(pca2)
## 
```

###### 优化pca实现

```{r eval =F}
## 加载环境:
library(ade4) ## dudi.pca()
library(factoextra) ## 可视化
library(FactoMineR) ## PCA()函数

## 使用来自FactoMineR::PCA()
## PCA(X, scale.unit = TRUE, ncp = 5, graph = TRUE)
## ncp为结果保留的维数；graph为结果打印图表；
library("FactoMineR")
res.pca <- PCA(decathlon2.active, graph = FALSE)

## pca结果可视化：
library("factoextra")
## 关于pca双标图见上面关于pca双序图的解释：
# get_eigenvalue(res.pca)：提取主成分的特征值/方差
# fviz_eig(res.pca)：可视化特征值
fvizeig(res.pca, addlabels = TRUE, ylim = c(0, 50))
# get_pca_ind(res.pca)，get_pca_var(res.pca)：分别提取个体和变量的结果。
# fviz_pca_ind(res.pca)，fviz_pca_var(res.pca)：分别可视化结果个体和变量。
# fviz_pca_biplot(res.pca)：制作个人和变量的双重图

## 查看不同维度的贡献（特征根）：----
get_eigenvalue(res.pca) ## 提取主成分的特征值/方差
## 可视化特征根的信息压缩：
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

```

###### PCA双序图----环境因子、个体

```{r eval =F}
#### PCA双序图----环境因子图：#####

## 查看单物种不同变量的相关关系(##等价于在主轴上相关)----
var <- get_pca_var(res.pca)
var$contrib
### 绘制 变量相关图(##加变量平均圈)
fvizpcavar(res.pca, col.var = "black")
## 投影cos2到轴上可视化；
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
## 绘制不同变量在主轴上对相关性：
library("corrplot")
## 这里涉及到的cos2概念，解释起来很麻烦：
## 理论上应该是方差的转换，用于表示单维度因子的重要性；
corrplot(var$cos2, is.corr=FALSE)


## 查看变量对PC主轴的贡献
## head(var$contrib, 4)
library("corrplot")
corrplot(var$contrib, is.corr=FALSE)
## 将贡献投影到轴上，与cos2差别不大；
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)


##根据变量数量分类着色：
set.seed(123)
my.cont.var <- rnorm(10)
# Color variables by the continuous variable
fviz_pca_var(res.pca, col.var = my.cont.var,
             gradient.cols = c("blue", "yellow", "red"),
             legend.title = "Cont.Var")


## 分组聚类可视化：(##无因子轴)-----
## 可利用kmeans等将聚类结果传到pca中：
res.km <- kmeans(var$coord, centers = 3, nstart = 25)
grp <- as.factor(res.km$cluster)
# Color variables by groups
fviz_pca_var(res.pca, col.var = grp, 
             palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
             legend.title = "Cluster")

## 查看变量对主成分的贡献来显示变量；
## 并给出变量的相关性可靠程度；
res.desc <- dimdesc(res.pca, axes = c(1,2), proba = 0.05)
# Description of dimension 1
res.desc$Dim.1


## 
#### PCA双序图---个体图：#####
## 流程与pca因子图类似，区别在于因子图与pcoa类似，
## 可以反映物种间在因子梯度上的分布关系：
ind <- get_pca_ind(res.pca)
ind
## 绘制个体图重要性和贡献图：
## 绘制简图：
fvizpcaind(res.pca)
## 根据cos2着色，可以反映不同分布点之间的重要性差异；
fviz_pca_ind(res.pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

## 基于个体分组，添加置信圈（无变量图）；
fviz_pca_ind(iris.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = iris$S
             
             pecies, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)

```

###### pca-双标图

```{r eval =F}


###  pca-双标图（因子图加个体图）：
##  将变量和个体图同时绘制；
fviz_pca_biplot(iris.pca, 
                col.ind = iris$Species, palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Species")

## 优化后：
# col1 = colorRampPalette(c("lightpink", "indianred"),bias=1,alpha=0.4)(n=20)
# col2 = colorRampPalette(c("lightyellow", "orange"),bias=1)(n=20)

fviz_pca_biplot(pca3, arrowsize=1,pointsize=5,labelsize=12,
                col.ind = sdm3$group, palette = c("lightpink","khaki","#FF3366","gold"), 
                addEllipses = FALSE, label = "var",
                col.var = "black", repel = TRUE
                )

fviz_pca_var(pca3 ,col.var = "black", arrowsize=1,pointsize=5,labelsize=12)+theme_bw() + theme(axis.line=element_blank(), 
                                                         axis.text.x=element_blank(), 
                                                         axis.text.y=element_blank(), 
                                                         axis.ticks=element_blank(), 
                                                         axis.title.x=element_blank(), 
                                                         axis.title.y=element_blank(), 
                                                         legend.position="none", 
                                                         panel.background=element_blank(), 
                                                         panel.border=element_blank(), 
                                                         panel.grid.major=element_blank(), 
                                                         panel.grid.minor=element_blank(), 
                                                         plot.background=element_blank()) +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
```

```{r eval =F}
## 将pca所有结果导出到指定的txt文件：
## Export into a TXT file
write.infile(res.pca, "pca.txt", sep = "\t")
# Export into a CSV file
write.infile(res.pca, "pca.csv", sep = ";")

## 将所有绘制结果图批量导出：
## 使用pppubr::ggexport()导出成pdf文件
library(ggpubr)
##一张图一页；
ggexport(plotlist = list(scree.plot, ind.plot, var.plot), 
         filename = "PCA.pdf")
## 一张图多页；
ggexport(plotlist = list(scree.plot, ind.plot, var.plot), 
         nrow = 2, ncol = 2,
         filename = "PCA.pdf")

##将图表导出到png文件。如果指定绘图列表，则将自动创建多个png文件以保存每个绘图
ggexport(plotlist = list(scree.plot, ind.plot, var.plot),
         filename = "PCA.png")
```

###### 9.1.2.1.2 CA/CCA分析

```

```

###### 9.1.2.1.3 无监督聚类

**KMEANS**

```{r eval =F}
##数据标准化
## data1的数据结构包含一个物种名，其他为环境因子，在kmeans构造中将物种名作为行名称；
rown <- as.character(data1[,1])
data2 <- data1[,-1]
rownames(data2) <- rown
## head(data1)
data3 <- scale(data2, center = TRUE, scale = TRUE)
##绘制分类预测限制，确定最佳聚类水平
## kmeans一般适合较低维度，此结果可能不具参考价值；
set.seed(123)
library(factoextra)
fviz_nbclust(data3,kmeans,method="wss")+geom_vline(
  xintercept = 4,linetype=2)
## 根据滚石图，最佳聚类可能为6,7,8
km_res6 <- kmeans(data3,6,nstart = 10,iter.max = 10)
cluster6 <- km_res6$cluster
## 下图的结果也充分说明用此种办法建模效果并不好；
fviz_cluster(km_res6, data = data3,
             ellipse.type = "euclid",
             star.plot = TRUE, 
             repel = TRUE,
             ggtheme = theme_minimal(),
             main="6次kmeans"
)

```

**层次聚类**

```{r eval =F}
### 分类-层次分类
#### 层次分析聚类：####
result <- dist(data3, method = "euclidean")
#产生层次结构
result_hc <- hclust(d = result, method = "ward.D2")
#进行初步展示
fviz_dend(result_hc, cex = 0.6,horiz=TRUE)


fviz_dend(result_hc, k = 8, 
          cex = 0.5, 
          color_labels_by_k = TRUE, 
          rect = TRUE,horiz=TRUE          
)
```

**T-SNE**

```{r eval =F}
### 无监督聚类-tsne
## install.packages("Rtsne")
library(Rtsne)
## tem 第一列有数据标签
## Choose the train.csv file downloaded from the link above  

## Curating the database for analysis with both t-SNE and PCA

train <- tem <- data1
labels <- train$names
train$label<-as.factor(train$names)
## for plotting
colors = rainbow(length(unique(train$label)))
names(colors) = unique(train$label)

## Executing the algorithm on curated data
tsne<- Rtsne(train[,-1], dims = 2, perplexity=8, verbose=TRUE, max_iter = 500)
## exeTimeTsne<- system.time(Rtsne(train[,-1], dims = 2, perplexity=10, verbose=TRUE, max_iter = 500))

## Plotting
plot(tsne$Y, t='n', main="tsne")
text(tsne$Y, labels=train$label, col=colors[train$label])


### 重新绘制：(优化)
# tsnedata <- data.frame(tsne$Y,type= train$label)
# names(tsnedata) <- c("tsne_1","tsne_2","class")
# ##head(tsnedata)
# library(ggpubr)
# ggscatter(tsnedata, x="tsne_1",y="tsne_2",
#           color ="class",size =0.5,label = "class",font.label = c(7, "plain"),
#           main="tSNE plot")
```

###### 9.1.2.1.4 多方差置换分析

```{r eval =F}
## 多元置换法方差分析的目的：
比较不同数据矩阵之间的相关关系和差异情况；
1）自变量的变化是否对因变量有显着影响？
2）因变量之间的关系是什么？
3）自变量之间有什么关系？
## 传统manova面对诸多检验和预设模型的问题，因此可以使用以下两种办法替换解决的：
## 方法1：多元方差置换分析：
### R实现：
## adonis(formula, data, permutations = 999, method = "bray", strata = NULL,contr.unordered = "contr.sum", contr.ordered = "contr.poly") 
## 
library(vegan)
sdmdata =cbind(data$species,data1[,2:10])
perm <- how(nperm = 999)

##需要注意是关于矩阵计算距离的方法选择，参见vegan包中的vegdist()函数，该公式默认的是“bray”布雷距离，但是该距离矩阵中不能包含负值。因此上述的代码中选择使用欧式距离矩阵进行计算
## 这里可能是指将矩阵的行名称改为speies
## 注意这里是formula为前半部分dsmma为环境因子矩阵，species为物种分类矩阵；
sdmenv为物种分类矩阵；
adonis2 (dsmma ~ species, data = sdm.env, permutations = perm,method = "euclidean")
## 打印结果
sink("C:/Users/chengshanmei/Desktop/多因素方差置换检验.txt")
##公式计算结果
##out.avop
sink(NULL)

##方法2：RDA-MNOVA：
## 只需要满足一个前提：方差-协方差矩阵齐性通过检验；可以通过vegan包的betadisper()函数进行检验；
## 在本文中举了一个二元回归置换分析（双因素方差分析）可以解释不同因子对样方的显著性；也可以用于评估因子间的交互作用是否显著；
## 这一块还需要继续补充；

```

