### PCA维度得分可视化1

```
## 工作思路：
是先构建ecospat：再计算复杂pca/简单的pca；
此后计算pca轴的得分，然后利用ggplot包绘制得分；
library(ecospat)
library(ade4)
library(raster)
library(rworldmap)

## 加载物种分布数据：#####

as_occ <- read.csv("D:/XH/xh/as.csv")[,2:3]
sa_occ  <- read.csv("D:/XH/xh/sa.csv")[,2:3] 

## 加载环境背景数据：
# 加载bg采样点：
#### 方法1： bg_enm D:\XH\fifth_bg2\bg_buffers

as_e_bg <- read.csv("D:/XH/fifth_bg2/bg_buffers/bg_as0.5.csv")[,2:3]
sa_e_bg <- read.csv("D:/XH/fifth_bg2/bg_buffers/bg_sa0.5.csv")[,2:3]

## 加载栅格背景数据：
env <- stack(list.files("D:/XH/third_env/envsV5tif",pattern = "tif",full.names = T))

### 构建数据组：

## 南美与亚洲：
env_as_occ <- data.frame(raster::extract(env,as_occ))
env_sa_occ <- data.frame(raster::extract(env,sa_occ))


## 方法1：bg_enm：
env_as_e <- data.frame(raster::extract(env,as_e_bg))
env_sa_e <- data.frame(raster::extract(env,sa_e_bg))


## 构建nat和inv：
## 构建nat_sa：
spocc <- c(rep(1,nrow(sa_occ)),rep(0,nrow(env_sa_e )))
envs_enm <- rbind(env_sa_occ,env_sa_e)
names(sa_e_bg) <- names(sa_occ)
sa_bg_occ <- rbind(sa_occ,sa_e_bg)
nat_sa <- na.exclude(cbind(sa_bg_occ,envs_enm ,spocc))


## 构建inv_as:
spocc3 <- c(rep(1,nrow(as_occ)),rep(0,nrow(env_as_e )))
envs_enm3 <- rbind(env_as_occ,env_as_e)
names(as_e_bg) <- names(as_occ)
as_bg_occ3 <- rbind(as_occ,as_e_bg)
inv_as <- na.exclude(cbind(as_bg_occ3,envs_enm3 ,spocc3))
names(inv_as) <- names(nat_sa)


### sa_as_enmbg_PCA ####

pca.env <-dudi.pca(rbind(nat_sa,inv_as)[3:13], center = T, scale = T, scannf = F, nf = 2)

ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)

## PCA结果可视化：
library(ade4) 
library(factoextra)
library(FactoMineR)
## 查看对应轴的贡献：
fviz_cos2(pca.env, choice = "var", axes = 1)
fviz_cos2(pca.env, choice = "var", axes = 2)

## 重新构建矩阵可视化分布点：
saenm <- cbind(nat_sa[1:nrow(na.exclude(env_sa_occ)),3:13],"sa_occ")
asenm <- cbind(inv_as[1:nrow(na.exclude(env_as_occ)),3:13],"as_occ")
sabgenm <- cbind(nat_sa[-c(1:nrow(na.exclude(env_sa_occ))),3:13],"sa_bg")
asbgenm <- cbind(inv_as[-c(1:nrow(na.exclude(env_as_occ))),3:13],"as_bg")
names(saenm) <- names(asenm) <- names(sabgenm) <- names(asbgenm)
sdm3 <- rbind(sabgenm,asbgenm,saenm,asenm)
colnames(sdm3)[12] <- "group"
sdm3$group <- factor(sdm3$group,levels=c("sa_bg","as_bg","sa_occ","as_occ")	)
pca3 <- PCA(sdm3[,1:11],graph =FALSE)


## 导出数据尺度：1000*1000 tif
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

### COMP_NICEH ####

##  predict the scores on the PCA axes
scores.globclim <- pca.env$li
# PCA scores for the species native distribution
scores.sp.nat <- suprow(pca.env,nat_sa[which(nat_sa[,14]==1),3:13])$li
# PCA scores for the species invasive distribution
scores.sp.inv <- suprow(pca.env,inv_as[which(inv_as[,14]==1),3:13])$li
# PCA scores for the whole native study area
scores.clim.nat <- suprow(pca.env,nat_sa[,3:13])$li
# PCA scores for the whole invaded study area
scores.clim.inv <- suprow(pca.env,inv_as[,3:13])$li

# calculation of occurence density
z1<- ecospat.grid.clim.dyn(glob=scores.globclim,
                           glob1=scores.clim.nat,
                           sp=scores.sp.nat, R=100,
                           th.sp=0.05)
z2 <- ecospat.grid.clim.dyn(glob=scores.globclim,
                            glob1=scores.clim.inv,
                            sp=scores.sp.inv, R=100,
                            th.sp=0.05)
                            
library(tidyverse )
library(PMCMR)
col=c("gray20","gray40")
ss1 <- c(rep(1,length(scores.sp.nat$Axis1 )),rep(2,length(scores.sp.inv$Axis1)))
yy1 <- c(scores.sp.nat$Axis1,scores.sp.inv$Axis1) %>% scale(.,center = TRUE)
pc1 <- cbind(ss1,yy1) %>% data.frame(.)
pc1$ss1 <- factor(pc1$ss1)
posthoc.kruskal.nemenyi.test(yy1~ss1,data=pc1, dist="Tukey")
## 计算结果显著：=
ggplot(pc1,aes(x=yy1 ,fill=ss1))+ geom_density(alpha =0.4)+ylab("")+scale_fill_manual(values=col)+guides(fill= FALSE)+xlab("PC1(SCALED)")+
  theme_bw() + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=NULL,expand=c(0,0))+theme(axis.line.x=element_line(linetype=1,color="gray20",size=1))+ 
  theme(axis.text.x = element_text(size = 25,color="black"))+xlim(-2.5,2.5)## 计算显著性：

ss2 <- c(rep(1,length(scores.sp.nat$Axis2 )),rep(2,length(scores.sp.inv$Axis2)))
yy2 <- c(scores.sp.nat$Axis2,scores.sp.inv$Axis2) %>% scale(.,center = TRUE)
pc2 <- cbind(ss2,yy2) %>% data.frame(.)
pc2$ss2 <- factor(pc2$ss2)
posthoc.kruskal.nemenyi.test(yy2~ss2,data=pc2, dist="Tukey")
## 计算结果不显著：
ggplot(pc2,aes(x=yy2 ,fill=ss2))+ geom_density(alpha =0.3)+ylab("")+scale_fill_manual(values=col)+guides(fill= FALSE)+xlab("PC2(SCALED)")+
  theme_bw() + theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=NULL,expand=c(0,0))+theme(axis.line.x=element_line(linetype=1,color="gray20",size=1))+ 
  theme(axis.text.x = element_text(size = 25,color="black"))+xlim(-2.5,2.5)


```



### PCA维度得分可视化2

```
 
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

