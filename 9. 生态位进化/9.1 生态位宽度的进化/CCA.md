```R


library(vegan)
library(ggplot2)

##读取数据
#读入物种数据OTU 水平丰度表
otu <- read.delim('otu_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu <- data.frame(t(otu))

#读取环境数据
env <- read.delim('env_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)


# 2.调用vegan包中cca()命令执行CCA，有以下两种格式；

#调用格式 1，cca(Y, X, W)：Y，响应变量矩阵；X，解释变量矩阵；W，协变量矩阵（偏 CCA 时使用）
#这里无协变量矩阵，所以直接输入响应变量矩阵和解释变量矩阵
#otu_cca <- cca(otu, env)

#或者格式 2，cca(Y~var1+var2+var3+factorA+var2*var3+Condition(var4))
#var1、var2 等，数值型解释变量；factorA，因子型解释变量；var2*var3，考虑变量间的交互作用；Condition(var4)，变量 4 作为协变量
#Y~. 是 Y~var1+var2+...+varn 的简写，不涉及交互作用及协变量
otu_cca <- cca(otu~., env)
# 3.使用coef() 命令提取 CCA 典范系数；

##coef() 提取 CCA 典范系数
cca_coef <- coef(otu_cca)
# 4.使用RsquareAdj()校正R2并约束轴特征值；

##R2 校正
#RsquareAdj() 提取 R2，详情 ?RsquareAdj() 
r2 <- RsquareAdj(otu_cca)
otu_cca_noadj <- r2$r.squared   #原始 R2
otu_cca_adj <- r2$adj.r.squared #校正后的 R2

#关于约束轴承载的特征值或解释率，应当在 R2 校正后重新计算
otu_cca_exp_adj <- otu_cca_adj * otu_cca$CCA$eig/sum(otu_cca$CCA$eig)
otu_cca_eig_adj <- otu_cca_exp_adj * otu_cca$tot.chi
# 5.对约束轴的进行置换检验和p值校正；

##置换检验
#所有约束轴的置换检验，即全局检验，基于 999 次置换，详情 ?anova.cca
otu_cca_test <- anova.cca(otu_cca, permutations = 999)

#各约束轴逐一检验，基于 999 次置换
otu_cca_test_axis <- anova.cca(otu_cca, by = 'axis', permutations = 999)

#p 值校正（Bonferroni 为例）
otu_cca_test_axis$`Pr(>F)` <- p.adjust(otu_cca_test_axis$`Pr(>F)`, method = 'bonferroni')
# 6.计算方差膨胀因子，进行CCA的变量选择；

##变量选择
#计算方差膨胀因子，详情 ?vif.cca
vif.cca(otu_cca)

#前向选择，这里以 ordiR2step() 的方法为例，基于 999 次置换检验，详情 ?ordiR2step
otu_cca_forward_pr <- ordiR2step(cca(otu~1, env), scope = formula(otu_cca), R2scope = TRUE, direction = 'forward', permutations = 999)
# 7.可视化示例；

#以上述前向选择后的简约模型 otu_cca_forward_pr 为例作图展示前两轴

#计算校正 R2 后的约束轴解释率
exp_adj <- RsquareAdj(otu_cca_forward_pr)$adj.r.squared * otu_cca_forward_pr$CCA$eig/sum(otu_cca_forward_pr$CCA$eig)
cca1_exp <- paste('CCA1:', round(exp_adj[1]*100, 2), '%')
cca2_exp <- paste('CCA2:', round(exp_adj[2]*100, 2), '%')
# 8.使用ggplot2包绘图；

#下面是 ggplot2 方法
#提取样方和环境因子排序坐标，前两轴，I 型标尺
otu_cca_forward_pr.scaling1 <- summary(otu_cca_forward_pr, scaling = 1)
otu_cca_forward_pr.site <- data.frame(otu_cca_forward_pr.scaling1$sites)[1:2]
otu_cca_forward_pr.env <- data.frame(otu_cca_forward_pr.scaling1$biplot)[1:2]

#添加分组
otu_cca_forward_pr.env$name <- rownames(otu_cca_forward_pr.env)
#读取分组文件按
map<- read.delim('group.txt', row.names = 1, sep = '\t')
otu_cca_forward_pr.site$name <- rownames(otu_cca_forward_pr.site)
otu_cca_forward_pr.site$group <- map$group
#merged2<-merge(map,otu,by="row.names",all.x=TRUE)

#ggplot2 作图
library(ggrepel)    #用于 geom_label_repel() 添加标签

color=c( "#3C5488B2","#00A087B2", 
         "#F39B7FB2","#91D1C2B2", 
         "#8491B4B2", "#DC0000B2", 
         "#7E6148B2","yellow", 
         "darkolivegreen1", "lightskyblue", 
         "darkgreen", "deeppink", "khaki2", 
         "firebrick", "brown1", "darkorange1", 
         "cyan1", "royalblue4", "darksalmon", 
         "darkgoldenrod1", "darkseagreen", "darkorchid")

p <- ggplot(otu_cca_forward_pr.site, aes(CCA1, CCA2)) +
  geom_point(size=1,aes(color = group,shape = group)) +
  stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE, linetype = 2) +
  scale_color_manual(values = color[1:length(unique(map$group))]) +
  theme(panel.grid.major = element_line(color = 'gray', size = 0.1), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), plot.title = element_text(hjust = 0.5)) + 
  labs(x = cca1_exp, y = cca2_exp) +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  geom_segment(data = otu_cca_forward_pr.env, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), arrow = arrow(length = unit(0.2, 'cm')), size = 0.3, color = 'blue') +
  geom_text(data = otu_cca_forward_pr.env, aes(CCA1 * 1.2, CCA2 * 1.2, label = name), color = 'blue', size = 3) +
  geom_label_repel(aes(label = name, color = group), size = 3, box.padding = unit(0, 'lines'), show.legend = FALSE)

p

```

