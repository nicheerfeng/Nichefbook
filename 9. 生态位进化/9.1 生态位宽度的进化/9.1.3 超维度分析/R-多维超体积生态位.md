### R-多维超体积生态位

###  niceOverPlot

```R
source("niceOverPlot_code.R")

niceOverPlot()

 
# niceOverPlot function is based of this two posts:
# http://stackoverflow.com/questions/20474465/using-different-scales-as-fill-based-on-factor
# http://rforpublichealth.blogspot.com.es/2014/02/ggplot2-cheatsheet-for-visualizing.html

# niceOverPlot function can be used in several ways. See example above to learn the basic use. Different 
# approaches will be posted as soon as possible.

niceOverPlot<-function(sc1,sc2=NULL,n1=NULL,n2=NULL, plot.axis = TRUE, bw = NULL, b=NULL, a1cont=NULL, a2cont=NULL){
  
  # prepare the data, depending of the type of input ("pca"/"dudi" object or raw scores)
  if (is.null(sc2))
  {sc_1<-sc1
  sc_2<-sc1
  sc1<- sc_1$li[1:n1,]
  sc2<- sc_1$li[(n1+1):(n1+n2),]
  }
  
  if (class(sc1)==c("pca","dudi") && class(sc2)==c("pca","dudi")) 
  {sc_1<-sc1
  sc_2<-sc1
  sc1<- sc1$li
  sc2<- sc2$li}
  
  # recognize both species
  scores<-rbind(sc1,sc2)
  g<-c(rep(0,nrow(sc1)),rep(1,nrow(sc2)))
  df<-data.frame(cbind(scores$Axis1,scores$Axis2,g))
  names(df)<-c("x","y","g")
  df$g<-as.factor(df$g)
  
  # establish an empty plot to be placed at top-right corner (X)
  empty <- ggplot()+geom_point(aes(1,1), colour="white") +
    theme(                              
      plot.background = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    )
  # sp1
  p1 <- ggplot(data = df, aes(x, y,color = as.factor(g))) +
    stat_density2d(aes(fill = ..level..), alpha = 0.2, bins=b, geom = "polygon", h=c(bw,bw)) +
    scale_fill_continuous(low = "#fdae61", high = "#d7191c", space = "Lab", name = "sp1") +
    scale_colour_discrete(guide = FALSE) + scale_x_continuous(name = "axis1", limits= c(min(df$x)-100, max(df$x)+100))+
    scale_y_continuous(name = "axis2", limits= c(min(df$y)-100, max(df$y)+100))+
    theme(legend.position="none")
  # sp2
  p2 <- ggplot(data = df, aes(x, y, color = as.factor(g))) +
    stat_density2d(aes(fill = ..level..), alpha = 0.2, bins=b, geom = "polygon", h=c(bw,bw)) +
    scale_fill_continuous(low = "#abd9e9", high = "#2b83ba", space = "Lab", name = "sp2") +
    scale_colour_discrete(guide = FALSE) +  scale_x_continuous(name = "axis1", limits=c(min(df$x)-100, max(df$x)+100))+
    scale_y_continuous(name = "axis2", limits= c(min(df$y)-100, max(df$y)+100))+
    theme(legend.position="none")
  
  pp1 <- ggplot_build(p1)
  ppp1 <- ggplot_build(p1 + aes(alpha=0.15) + theme_classic() + theme(legend.position="none") + theme(text = element_text(size=15)) + xlab("axis1") + ylab("axis2") + xlim(c(min(pp1$data[[1]]$x)-0.5,max(pp1$data[[1]]$x)+0.5)) + ylim(c(min(pp1$data[[1]]$y)-0.5,max(pp1$data[[1]]$y)+0.5)))
  pp2 <- ggplot_build(p2 + aes(alpha=0.15) + theme_classic() + theme(legend.position="none")+ xlab("axis1") + ylab("axis2")  + xlim(c(min(pp1$data[[1]]$x)-0.5,max(pp1$data[[1]]$x)+0.5)) + ylim(c(min(pp1$data[[1]]$y)-0.5,max(pp1$data[[1]]$y)+0.5)))$data[[1]]
  
  ppp1$data[[1]]$fill[grep(pattern = "^2", pp2$group)] <- pp2$fill[grep(pattern = "^2", pp2$group)]
  
  grob1 <- ggplot_gtable(ppp1)
  grob2 <- ggplotGrob(p2)
  grid.newpage()
  grid.draw(grob1)
  
  #marginal density of x - plot on top
  
  if (class(sc_1)==c("pca","dudi") && class(sc_2)==c("pca","dudi")) 
  {plot_top <- ggplot(df, aes(x, y=..scaled..,fill=g)) + 
    geom_density(position="identity",alpha=.5) +
    scale_x_continuous(name = paste("Contribution ",(round((sc_1$eig[1]*100)/sum(sc_1$eig),2)),"%",sep=""), limits=c(min(pp1$data[[1]]$x)-0.5,max(pp1$data[[1]]$x)+0.5))+
    scale_fill_brewer(palette = "Set1") + 
    theme_classic() + theme(legend.position = "none")
  }
  
  else {
    
    if(is.null(a1cont)) plot_top <- ggplot(df, aes(x, y=..scaled..,fill=g)) + 
        geom_density(position="identity",alpha=.5) +
        scale_x_continuous(name = "axis1", limits=c(min(pp1$data[[1]]$x)-0.5,max(pp1$data[[1]]$x)+0.5))+
        scale_fill_brewer(palette = "Set1") + 
        theme_classic() + theme(legend.position = "none")  
    
    
    
    else plot_top <- ggplot(df, aes(x, y=..scaled..,fill=g)) + 
        geom_density(position="identity",alpha=.5) +
        scale_x_continuous(name = paste("Contribution ",a1cont,"%",sep=""), limits=c(min(pp1$data[[1]]$x)-0.5,max(pp1$data[[1]]$x)+0.5))+
        scale_fill_brewer(palette = "Set1") +
        theme_classic() + theme(legend.position = "none")
    
  }
  #marginal density of y - plot on the right
  
  if (class(sc_1)==c("pca","dudi") && class(sc_2)==c("pca","dudi")) 
  {plot_right <- ggplot(df, aes(y, y=..scaled.., fill=g)) + 
    geom_density(position="identity",alpha=.5) + 
    scale_x_continuous(name = paste("Contribution ",(round((sc_1$eig[2]*100)/sum(sc_1$eig),2)),"%",sep=""), limits= c(min(pp1$data[[1]]$y)-0.5,max(pp1$data[[1]]$y)+0.5)) +
    coord_flip() + 
    scale_fill_brewer(palette = "Set1") + 
    theme_classic() + theme(legend.position = "none") 
  }
  
  else {
    
    if(is.null(a2cont)) plot_right <- ggplot(df, aes(y, y=..scaled.., fill=g)) + 
        geom_density(position="identity",alpha=.5) + 
        scale_x_continuous(name = "axis2", limits= c(min(pp1$data[[1]]$y)-0.5,max(pp1$data[[1]]$y)+0.5)) +
        coord_flip() + 
        scale_fill_brewer(palette = "Set1") + 
        theme_classic() + theme(legend.position = "none") 
    
    
    else plot_right <- ggplot(df, aes(y, y=..scaled.., fill=g)) + 
        geom_density(position="identity",alpha=.5) + 
        scale_x_continuous(name = paste("Contribution ",a2cont,"%",sep=""), limits= c(min(pp1$data[[1]]$y)-0.5,max(pp1$data[[1]]$y)+0.5)) +
        coord_flip() + 
        scale_fill_brewer(palette = "Set1") +
        theme_classic() + theme(legend.position = "none") 
    
  }
  
  if (plot.axis == TRUE) grid.arrange(plot_top, empty , grob1, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
  else grid.draw(grob1)
  
}
```

###  dynRB

```R
## 学习使用dynRB:
# https://cran.r-project.org/web/packages/dynRB/vignettes/TutorialdynRB.html

library(dynRB)
library(ggplot2)
library(reshape2)
library(vegan)
library(RColorBrewer)

##  使用的数据来自加拉帕斯群岛的达尔文雀；
data(finch)
head(finch)[,1:5]

# 构建超体积：
# r的输出结果中大小（vol）和重叠（port）
r <- dynRB_VPa(finch)

##  输出的结果具体意思还需要参见原网页中的论述；
res <- r$result
res[1:6, 1:5]

##  量化每个维度的重叠
r <- dynRB_Pn(finch)
head(r$result)

## 量化每个生态位(多维超体积)的大小
r <- dynRB_Vn(finch)

## 数据转换：
r <- dynRB_VPa(finch)
# aggregation method product  
om1 <- reshape(r$result[,1:3], direction="wide", idvar="V1", timevar="V2") 
#aggregation method mean  
om2 <- reshape(r$result[,c(1,2,4)], direction="wide", idvar="V1", timevar="V2") 
#aggregation method gmean  
om3 <- reshape(r$result[,c(1,2,5)], direction="wide", idvar="V1", timevar="V2") 


## 评估重叠不对称矩阵之间的相关性；
# 如果mental测试支持两者是相似的，则说明无论是 对称性评估是可靠的；
r <- dynRB_VPa(finch)
#aggregation method mean  
om <- reshape(r$result[,c(1,2,4)], direction="wide", idvar="V1", timevar="V2") 
mantel(as.dist(om[2:ncol(om)]), as.dist(t(om[2:ncol(om)])), permutations = 1000)


## 使用热图可视化重叠
r <- dynRB_VPa(finch)  #main function of dynRB to calculate size and overlap of hypervolumes
result <- r$result
Overlap <- as.numeric(ifelse(result$V1 == result$V2, "NA", result$port_prod))  


is.numeric(Overlap)
Result2 <- cbind(result, Overlap)
ggplot(Result2, aes(x = V1, y = V2)) +
  geom_tile(data = subset(Result2, !is.na(Overlap)), aes(fill = Overlap), color="black") +
  geom_tile(data = subset(Result2,  is.na(Overlap)), fill = "lightgrey", color="black")


is.numeric(Overlap)
Result2<-cbind(result, Overlap)
breaks <- seq(min(Overlap, na.rm=TRUE),max(Overlap, na.rm=TRUE),  
              by=round(max(Overlap, na.rm=TRUE)/10, digits=3))
col1 <- colorRampPalette(c("white", "navyblue")) #define color gradient
ggplot(Result2, aes(x = V1, y = V2)) +
  geom_tile(data = subset(Result2, !is.na(Overlap)), aes(fill = Overlap), color="black") +
  geom_tile(data = subset(Result2,  is.na(Overlap)), fill = "lightgrey", color="black") +
  scale_fill_gradientn(colours=col1(8), breaks=breaks, guide="colorbar",  
                       limits=c(min(Overlap, na.rm=TRUE),max(Overlap, na.rm=TRUE))) 
```

9.1.3.3 ellipsenm

```R
library(ellipsenm)
library(rgl) # for prettier plots
library(raster) # for processing raster layers 

## 设置工作目录：
## 加载物种分布数据：

## 加载环境变量和物种分布数据：

xh_na <- read.csv("D:/XH/xh/na.csv")[,1:3]
na <- as.data.frame()
xh_as <- read.csv("D:/XH/xh/as.csv")[,1:3]
xh_eu <- read.csv("D:/XH/xh/eu.csv")[,1:3]
xh_au <- read.csv("D:/XH/xh/au.csv")[,1:3]
xh_sa <- read.csv("D:/XH/xh/sa.csv")[,1:3] 
xh_ynbd <- read.csv("D:/XH/xh/ynbd.csv")[,1:3] 
head(xh_na)

####  加载环境数据：#####
tifs <- list.files(path ="D:\\XH\\third_env\\envsV4asc",pattern="asc",full.names =T)
tiffs <- stack(tifs)
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  
proj4string(tiffs) <- crs.geo  
envs <-  tiffs


niche_sa <- overlap_object(xh_sa, species =  "name", longitude = "longitude", 
                           latitude = "latitude", method = "mve1", level = 95, 
                           variables =envs)

niche_as <- overlap_object(xh_as, species =  "name", longitude = "longitude", 
                           latitude = "latitude", method = "mve1", level = 95, 
                           variables =envs)

niche_au <- overlap_object(xh_au, species =  "name", longitude = "longitude", 
                           latitude = "latitude", method = "mve1", level = 95, 
                           variables =envs)

niche_na <- overlap_object(xh_na, species =  "name", longitude = "longitude", 
                           latitude = "latitude", method = "mve1", level = 95, 
                           variables =envs)

niche_eu <- overlap_object(xh_eu, species =  "name", longitude = "longitude", 
                           latitude = "latitude", method = "mve1", level = 95, 
                           variables =envs)

overlap_sa_as <- ellipsoid_overlap(niche_sa, niche_as, overlap_type = "back_union", 
                                 significance_test = TRUE, replicates = 10, 
                                 confidence_limit = 0.05)

overlap_sa_au <- ellipsoid_overlap(niche_sa, niche_as, overlap_type = "back_union", 
                                   significance_test = TRUE, replicates = 10, 
                                   confidence_limit = 0.05)

overlap_sa_na <- ellipsoid_overlap(niche_sa, niche_as, overlap_type = "back_union", 
                                   significance_test = TRUE, replicates = 10, 
                                   confidence_limit = 0.05)

overlap_sa_eu <- ellipsoid_overlap(niche_sa, niche_as, overlap_type = "back_union", 
                                   significance_test = TRUE, replicates = 10, 
                                   confidence_limit = 0.05)
```

### ECOSPAT

```
## 参见9.2 ecospat包的可视化脚本；
```

