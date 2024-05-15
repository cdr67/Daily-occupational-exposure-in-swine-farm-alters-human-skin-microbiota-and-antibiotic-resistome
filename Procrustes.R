library(vegan)
library(dplyr)
library(ggsci)
group <- read.csv('skingroup.CSV', sep = ",", header = T)
colnames(group)[1]<-"clade"
tax <- read.csv('skin.csv', sep = ",", header = T,row.names = 1)
tax.dist<-vegdist(tax[1:96,1:849],method='bray')
mds.t <- monoMDS(tax.dist)
arg <- read.csv('argskinpcoa.csv', sep = ",", header = T,row.names = 1)
arg.dist<-vegdist(tax[1:96,1:484],method='bray')
mds.a <- monoMDS(arg.dist)
pro.t.a <- procrustes(mds.t,mds.a, symmetric = TRUE)
summary(pro.t.a)
plot(pro.t.a, kind = 2)
residuals(pro.t.a)
set.seed(1)
pro.t.a_t <- protest(mds.t,mds.a, permutations = 999)
pro.t.a_t
pro.t.a_t$ss
pro.t.a_t$signif
library(ggplot2)
Pro_Y <- cbind(data.frame(pro.t.a$Yrot), data.frame(pro.t.a$X))
Pro_X <- data.frame(pro.t.a$rotation)
Pro_Y[,5]<-rownames(Pro_Y)
colnames(Pro_Y)[5]<-"clade"
Pro_Y<-dplyr::left_join(Pro_Y, group, by = "clade")
colnames(Pro_Y)[6]<-"type"
rownames(Pro_Y)<-Pro_Y$clade
Pro_Y<-Pro_Y[,-5]
Pro_Y$type[which(Pro_Y$type %in% "T2_skin")] <- "T2"
Pro_Y$type[which(Pro_Y$type %in% "T1_skin")] <- "T1"
Pro_Y$type[which(Pro_Y$type %in% "T-1_skin")] <- "T-1"
Pro_Y$type=factor(Pro_Y$type,levels = c("T-1","T1","T2","Dust","Faeces"))

p=ggplot(data=Pro_Y) +
  geom_segment(aes(x = X1, y = X2,
                   xend = MDS1, yend = MDS2,color=type),
               arrow = arrow(length = unit(0, 'cm')),
              size = 1) +
  geom_point(aes(X1, X2,color=type), size = 5,shape = 16) +
    
  geom_point(aes(MDS1, MDS2,color=type),size = 5, shape = 17) +
  scale_color_npg()+
        panel.background = element_rect(color = 'black',
                                        fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size=14),
        axis.title.y=element_text(colour='black', size=14),
        axis.text=element_text(colour='black',size=12)) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  labs(title="Correlation between microbial community and ARGs") +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.3) +
  geom_abline(intercept = 0, slope = Pro_X[1,2]/Pro_X[1,1], size = 0.3) +
  geom_abline(intercept = 0, slope = Pro_X[2,2]/Pro_X[2,1], size = 0.3) +
  annotate('text', label = 'Procrustes analysis:\n
                    M2 = 0.2010, p-value = 0.001',
           x = -0.2, y = 0.15, size = 4,hjust = 0) +
  theme(plot.title = element_text(size=14,colour = "black",
                                  hjust = 0.5,face = "bold"))



