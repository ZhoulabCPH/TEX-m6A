# Figure 1B
setwd("D:\\TEX_FY\\中间数据")
load("all.Rdata")
exp<-as.data.frame(t(all))
exp<-log2(exp+1)
rm(all)
mRNA<-read.table("mRNA_GRCh38.104.gtf.txt")#读取参考基因组
exp_gene<-mRNA[mRNA[,1]%in%rownames(exp),]
exp<-cbind(row.names(exp),exp)
colnames(exp)[1]<-"gene_id"
exp<-merge(exp_gene,exp,by="gene_id")
rm(mRNA,exp_gene)
m6a<-read.table("m6A.txt",header = T)
m6a_exp<-exp[exp$gene_id%in%m6a$ENSEMBL,]
PD1<-exp[exp$gene_name=="PDCD1",]
library(dplyr)
PD1<-PD1[,-c(1,2)]
pd1<-t(PD1)
pd1<-as.data.frame(pd1)
pd1<-cbind(rownames(pd1),pd1)
colnames(pd1)<-c("sample","expression")
pd1$expression<-as.numeric(pd1$expression)
pd1<-arrange(pd1,expression)
low_sample<-pd1[1:133,]
high_sample<-pd1[398:530,]
middle_sample<-pd1[134:397,]
pd1_exp<-m6a_exp[,c(low_sample$sample,middle_sample$sample,high_sample$sample)]
rownames(pd1_exp)<-m6a_exp$gene_name
rm(pd1,low_sample,high_sample,middle_sample,m6a_exp)
eraser<-pd1_exp[c("ALKBH5","FTO"),]
writer<-pd1_exp[c("METTL5","RBM15","METTL3","CBLL1","RBM15B","WTAP",
                  "METTL14","ZCCHC4","ZC3H13","VIRMA","METTL16"),]
reader<-pd1_exp[c("EIF3A","HNRNPC","ELAVL1","G3BP1","HNRNPA2B1","YTHDF2",
                  "PRRC2A","YTHDC1","YTHDF1","G3BP2","YTHDC2","YTHDF3","LRPPRC"),]
pd1_exp<-rbind(eraser,writer,reader)
label<-data.frame(sample=c(colnames(pd1_exp[1:133]),colnames(pd1_exp[398:530])),subtype=c(rep("PD1-",133),rep("PD1+",133)))
library(tidyr)
pd1_exp_1<-cbind(pd1_exp[,1:133],pd1_exp[,398:530])
sample<-colnames(pd1_exp_1)
pd1_exp_1 <- cbind(rownames(pd1_exp_1),pd1_exp_1)
colnames(pd1_exp_1)[1]<-"gene"
pd1_exp_1$gene<-factor(pd1_exp_1$gene,levels = pd1_exp_1$gene)#按该顺序排序
a<-pivot_longer(pd1_exp_1,cols = c(sample),names_to = c("sample"),values_to = "Expression")
a<-merge(a,label,by="sample")
m6a<-pd1_exp_1[,1]
b<-data.frame(gene=m6a)
for(i in 1:26){
  my_data<-a[a$gene==m6a[i],]
  b[i,2]<-t.test(Expression~subtype,data = my_data)$p.value
}
colnames(b)[2]<-"pvalue"
b$p.adj<-p.adjust(b$pvalue,method = "fdr",n=length(b$pvalue))
b$signiff<-ifelse(b$p.adj<0.05,"<0.05",">0.05")
b$signiff2<-ifelse(b$p.adj<0.001,"***",
                   ifelse(b$p.adj<0.01,"**",
                          ifelse(b$p.adj<0.05,"*","ns")))
rm(a,my_data,i,label,sample,pd1_exp_1,m6a)
annotation_col<-data.frame(subtype=c(rep("PD1-",133),rep("middle",264),rep("PD1+",133)))
rownames(annotation_col)<-colnames(pd1_exp)
annotation_row<-data.frame(Type=factor(c(rep("eraser",2),rep("writer",11),rep("reader",13))),gene=rownames(pd1_exp))
annotation_row<-cbind(annotation_row,b)
annotation_row<-as.data.frame(annotation_row[,-2])
pd1_exp<-cbind(annotation_row,pd1_exp)
pd1_exp<-pd1_exp[,-2]
library(dplyr)
pd1_exp<-arrange(pd1_exp,Type,signiff)
pd1_exp<-rbind(pd1_exp[1:2,],pd1_exp[16:26,],pd1_exp[3:15,])
annotation_row<-pd1_exp[,1:4]
annotation_row<-annotation_row[,-c(2,3)]
colnames(annotation_row)[2]<-"p.adjust"
pd1_exp<-pd1_exp[,-(1:5)]
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
library(ggsci)
ann_colors = list(subtype=c(`PD1-`="#96CEB4",`middle`="#E4D1D3",`PD1+`="#EA7070"),
                  Type=c(eraser="#CCCC99",writer="#2E94B9",reader="#F78D3F"),
                  p.adjust=c(`<0.05`="#D5A4CF",`>0.05`="#F4F0E6"))
pheatmap(pd1_exp, 
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         cluster_rows = F, cluster_cols = F,
         show_rownames = T,show_colnames = F,
         annotation_colors = ann_colors,
         scale="row",
         fontsize = 8,
         gaps_col = c(133,397),
         gaps_row = c(2,13),
         cutree_cols=3,
         cutree_rows = 4,
         color = c(colorRampPalette(colors = c("#0000A1","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","#f1404b"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk)

# Figure 1C
rm(ann_colors,annotation_col,bk)
pd1_exp<-as.data.frame(cbind(pd1_exp[,1:133],pd1_exp[,398:530]))
label<-data.frame(condition = factor(c(rep("PD1-",133),rep("PD1+",133)), levels = c('PD1-', 'PD1+')))
label<-cbind(colnames(pd1_exp),label)
colnames(label)<-c("sample","subtype")
library(tidyr)
sample<-colnames(pd1_exp)
pd1_exp<-cbind(rownames(pd1_exp),pd1_exp)
colnames(pd1_exp)[1]<-"gene"
median<-data.frame(median=apply(pd1_exp,1,median))
pd1_exp<-cbind(median,annotation_row,pd1_exp)
pd1_exp<-arrange(pd1_exp,Type,median)
pd1_exp<-rbind(pd1_exp[1:2,],pd1_exp[16:26,],pd1_exp[3:15,])
pd1_exp$gene<-factor(pd1_exp$gene,levels = pd1_exp$gene)
a<-pivot_longer(pd1_exp,cols = c(sample),names_to = c("sample"),values_to = "Expression")
a<-merge(a,label,by="sample")
rm(label,annotation_row)
library(ggplot2)
library(ggsci)
library(ggsignif)
library(ggpubr)
compared_list = list(c("PD1-", "PD1+"))  
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          text = element_text(size = 8),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y=element_text(face = "bold",size = 14),
          axis.text = element_text(face = "bold",size = 12), 
          axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 90),
          axis.ticks = element_line(color='black'),
          legend.title=element_blank(),
          legend.position=c(0.75, 0.93),
          legend.direction = "horizontal",
          legend.text = element_text(face = "bold",size = 2,margin = margin(r=8)),
          legend.background = element_rect( linetype="solid",colour ="black"),
    )
}
ggplot(a, aes(x = gene, y = Expression))+ 
  geom_boxplot(aes(fill = subtype),position=position_dodge(0.9),width=0.8)+ 
  scale_fill_manual(values = c("#39BAE8","#f1404b"))+theme_zg()+
  stat_compare_means(aes(group=subtype,label = "p.signif"))


# Figure 1D
load("all.Rdata")
exp<-as.data.frame(t(all))
exp<-log2(exp+1)
rm(all)
mRNA<-read.table("mRNA_GRCh38.104.gtf.txt")#读取参考基因组
exp_gene<-mRNA[mRNA[,1]%in%rownames(exp),]
exp<-cbind(row.names(exp),exp)
colnames(exp)[1]<-"gene_id"
exp<-merge(exp_gene,exp,by="gene_id")
rm(mRNA,exp_gene)
m6a<-read.table("m6A.txt",header = T)
m6a_exp<-exp[exp$gene_id%in%m6a$ENSEMBL,]
rm(m6a)
PD1<-exp[exp$gene_name=="PDCD1",]
TCF1<-exp[exp$gene_name=="TCF7",]
library(dplyr)
PD1<-PD1[,-c(1,2)]
TCF1<-TCF1[,-c(1,2)]
pd1<-t(PD1)
a<-quantile(pd1)
c_index<-which(pd1>a[4])
tcf1<-t(TCF1)                                           
d<-quantile(tcf1)
e_index<-which(tcf1<d[2])
f_index<-which(tcf1>d[4])
pd1_tcf1_high<-intersect(c_index,f_index)
pd1_tcf1_low<-intersect(e_index,c_index)
pt_exp<-m6a_exp[,c(pd1_tcf1_low+2,pd1_tcf1_high+2)]
rownames(pt_exp)<-m6a_exp$gene_name
TCF1<-TCF1[,match(colnames(pt_exp),colnames(TCF1))]
tcf1<-t(TCF1)
tcf1<-as.data.frame(tcf1)
tcf1<-cbind(rownames(tcf1),tcf1)
colnames(tcf1)<-c("sample","expression")
tcf1$expression<-as.numeric(tcf1$expression)
tcf1<-arrange(tcf1,expression)
pt_exp<-pt_exp[,match(tcf1$sample,colnames(pt_exp))]
rm(a,c_index,d,e_index,f_index,pd1,PD1,TCF1,tcf1)
eraser<-pt_exp[c("ALKBH5","FTO"),]
writer<-pt_exp[c("METTL5","RBM15","METTL3","CBLL1","RBM15B","WTAP",
                 "METTL14","ZCCHC4","ZC3H13","VIRMA","METTL16"),]
reader<-pt_exp[c("EIF3A","HNRNPC","ELAVL1","G3BP1","HNRNPA2B1","YTHDF2",
                 "PRRC2A","YTHDC1","YTHDF1","G3BP2","YTHDC2","YTHDF3","LRPPRC"),]
pt_exp<-rbind(eraser,writer,reader)
label<-data.frame(sample=colnames(pt_exp),subtype=c(rep("PD1+TCF1-",27),rep("PD1+TCF1+",44)))
library(tidyr)
sample<-colnames(pt_exp)
pt_exp<- cbind(rownames(pt_exp),pt_exp)
colnames(pt_exp)[1]<-"gene"
pt_exp$gene<-factor(pt_exp$gene,levels = pt_exp$gene)
a<-pivot_longer(pt_exp,cols = c(sample),names_to = c("sample"),values_to = "Expression")
a<-merge(a,label,by="sample")
m6a<-pt_exp[,1]
b<-data.frame(gene=m6a)
for(i in 1:26){
  my_data<-a[a$gene==m6a[i],]
  b[i,2]<-t.test(Expression~subtype,data = my_data,)$p.value
}
colnames(b)[2]<-"pvalue"
b$p.adj<-p.adjust(b$pvalue,method = "fdr",n=length(b$pvalue))
b$signiff<-ifelse(b$p.adj<0.05,"<0.05",">0.05")
b$signiff2<-ifelse(b$p.adj<0.001,"***",
                   ifelse(b$p.adj<0.01,"**",
                          ifelse(b$p.adj<0.05,"*","ns")))
rm(a,my_data,i,label,sample,m6a)
annotation_col<-data.frame(subtype=c(rep("PD1+TCF1-",27),rep("PD1+TCF1+",44)))
rownames(annotation_col)<-colnames(pt_exp[,2:72])
annotation_row<-data.frame(Type=factor(c(rep("eraser",2),rep("writer",11),rep("reader",13))),gene=rownames(pt_exp))
annotation_row<-cbind(annotation_row,b)
annotation_row<-as.data.frame(annotation_row[,-2])
pt_exp<-cbind(annotation_row,pt_exp)
pt_exp<-pt_exp[,-2]
library(dplyr)
pt_exp<-arrange(pt_exp,Type,signiff)
pt_exp<-rbind(pt_exp[1:2,],pt_exp[16:26,],pt_exp[3:15,])
annotation_row<-pt_exp[,1:4]
annotation_row<-annotation_row[,-c(2,3)]
colnames(annotation_row)[2]<-"p.adjust"
pt_exp<-pt_exp[,-(1:6)]
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
library(ggsci)
ann_colors = list(subtype=c(`PD1+TCF1-`="#96C0CE",`PD1+TCF1+`="#C25B56"),
                  Type=c(eraser="#CCCC99",writer="#2E94B9",reader="#F78D3F"),
                  p.adjust=c(`<0.05`="#D5A4CF",`>0.05`="#F4F0E6"))

pheatmap(pt_exp, 
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         cluster_rows = F, cluster_cols = F,
         show_rownames = T,show_colnames = F,
         annotation_colors = ann_colors,
         border_color=NA,
         scale="row",
         fontsize = 8,
         gaps_col = 27,
         gaps_row = c(2,13),
         cutree_cols=3,
         cutree_rows = 3,
         color = c(colorRampPalette(colors = c("#0000A1","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","#f1404b"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk)

rm(bk,pd1_tcf1_high,pd1_tcf1_low,ann_colors,annotation_col)

#  Figure 1E
label<-data.frame(condition = factor(c(rep("PD1+TCF1-",27),rep("PD1+TCF1+",44)), levels = c('PD1+TCF1-', 'PD1+TCF1+')))
label<-cbind(colnames(pt_exp),label)
colnames(label)<-c("sample","subtype")
library(tidyr)
sample<-colnames(pt_exp)
pt_exp<-cbind(rownames(pt_exp),pt_exp)
colnames(pt_exp)[1]<-"gene"
median<-data.frame(median=apply(pt_exp[,2:72],1,median))
pt_exp<-cbind(median,annotation_row,pt_exp)
pt_exp<-arrange(pt_exp,Type,median)#排序
pt_exp<-rbind(pt_exp[1:2,],pt_exp[16:26,],pt_exp[3:15,])
pt_exp$gene<-factor(pt_exp$gene,levels = pt_exp$gene)
a<-pivot_longer(pt_exp,cols = c(sample),names_to = c("sample"),values_to = "Expression")
a<-merge(label,a,by="sample")
library(ggplot2)
library(ggsci)
library(ggsignif)
library(ggpubr)
compared_list = list(c("PD1+TCF1-", "PD1+TCF1+"))  
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          text = element_text(size = 8),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 65),
          axis.title.y=element_text(face = "bold",size = 14),
          axis.text = element_text(face = "bold",size = 12),        
          axis.ticks = element_line(color='black'),
          legend.title=element_blank(),
          legend.position=c(0.75, 0.93),
          legend.direction = "horizontal",
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=8)),
          legend.background = element_rect( linetype="solid",colour ="black")
    )
}
ggplot(a, aes(x = gene, y = Expression))+ 
  geom_boxplot(aes(fill = subtype),position=position_dodge(0.9),width=0.8)+ 
  scale_fill_manual(values = c("#39BAE8","#f1404b"))+theme_zg()+
  stat_compare_means(aes(group=subtype,label = ..p.format..),
                     method = "t.test")


# Figure 2C
setwd("D:\\TEX_FY\\中间数据")
cd8<-readRDS("cd8t.Rds")
library(dplyr)
library(Seurat)
cluster.averages <- AverageExpression(cd8,slot = "data", return.seurat = F,)
cluster.averages<-as.data.frame(cluster.averages)
markers.1<-c("PDCD1","TOX","CXCL13","TIGIT",
             "CTLA4","TNFRSF9","HAVCR2",#Exhaustion markers
             "CST7","GZMK","GZMA","NKG7","IFNG",
             "PRF1","GZMB","GNLY",#Cytokines and effector molecules
             "CD28","EOMES","CCR5","CCL4","CCL5",
             "CD200R1","TMEM155","CD27",#Higher in GZMK+ Tex
             "TOX2","TSHZ2","BATF","GEM","CD200",
             "SNX9","ENTPD1","LAYN","ETV1","PRDM1",
             "CSF1","MYO7A","MYO1E","GNG4",#Higher in terminal Tex
             "BTLA","FOXP3","TNFRSF4","TCF7","EEF1A1",
             "SELL","CCR7","IL6R","IGFBP4","IGFL2"#Higher in TCF7+ Tex
)
markers.2<-c("CD3D","CD3E","CD3G","CD8A","CD8B",#CD8T
             "PDCD1","TOX","CXCL13","TIGIT",
             "CTLA4","TNFRSF9","HAVCR2",#Exhaustion markers
             #Pre-exhausted
             "GNLY","NKG7","GZMB","GZMA"#Toxicity molecules
)
VlnPlot(cd8,features = c("ZNF683","CXCR6"),pt.size = 0)
bubble.df=as.matrix(cd8[["RNA"]]@data[markers.2,])
bubble.df=t(bubble.df)
bubble.df=as.data.frame(scale(bubble.df))
bubble.df$CB=rownames(bubble.df)
cd8@meta.data$CB=rownames(cd8@meta.data)
bubble.df=merge(bubble.df,cd8@meta.data[,c("CB","seurat_clusters")],by = "CB")
bubble.df$CB=NULL
seurat_clusters_v=c()
gene_v=c()
mean_v=c()
ratio_v=c()
for (i in unique(bubble.df$seurat_clusters)) {
  bubble.df_small=bubble.df%>%filter(seurat_clusters==i)
  for (j in markers.2) {
    exp_mean=mean(bubble.df_small[,j])
    exp_ratio=sum(bubble.df_small[,j] > min(bubble.df_small[,j])) / length(bubble.df_small[,j])
    seurat_clusters_v=append(seurat_clusters_v,i)
    gene_v=append(gene_v,j)
    mean_v=append(mean_v,exp_mean)
    ratio_v=append(ratio_v,exp_ratio)
  }
}
plotdf=data.frame(
  seurat_clusters=seurat_clusters_v,
  gene=gene_v,
  exp=mean_v,
  ratio=ratio_v
)
library(RColorBrewer)
library(ggplot2)
plotdf$seurat_clusters<-as.numeric(plotdf$seurat_clusters)
plotdf<-arrange(plotdf,seurat_clusters)
a<-c(rep(1:14,each = 16))
a<-paste("c",a,sep = "")
plotdf$seurat_clusters<-a
plotdf$seurat_clusters=factor(plotdf$seurat_clusters,
                              levels = c(c("c3","c12"),
                                         c("c9","c4","c10","c11"),
                                         c("c1","c2","c5","c6","c7","c8","c13","c14")))
plotdf$gene=factor(plotdf$gene,levels = rev(as.character(markers.2)))
plotdf$exp=ifelse(plotdf$exp>2,2,plotdf$exp)
plotdf$exp=ifelse(plotdf$exp<0,0,plotdf$exp)

mid=0

plotdf%>%
  ggplot(aes(x=seurat_clusters,y=gene,size=ratio,color=exp))+geom_point()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_gradient2(midpoint = mid,low = "#0000A1",
                        mid = "#fff7f7", high = "#e42c64")+
  #scale_colour_gradient(low = "#0000A1",high = "#e42c64")+
  scale_size_continuous(limits = c(0,1))+theme_bw()+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 90)
  )
exhausted_T=c(3,8,9,10)
toxicity_T=c(2,11)
other_T=c(0,1,4,5,6,7,12,13)

current.cluster.ids <- c(exhausted_T,
                         toxicity_T,
                         other_T)
new.cluster.ids <- c(rep("exhausted_T",length(exhausted_T)),
                     rep("toxicity_T",length(toxicity_T)),
                     rep("other_T",length(other_T)))
cd8@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(cd8@meta.data$seurat_clusters)), 
                                          from = current.cluster.ids, to = new.cluster.ids)
cd8<-readRDS("cd8_ann.Rds")
pre_Tex=c(8)
term_T=c(10)
mid_T=c(3,9)
toxicity_T=c(2,11)
other_T=c(0,1,4,5,6,7,12,13)
current.cluster.ids <- c(pre_Tex,
                         mid_T,
                         term_T,
                         toxicity_T,
                         other_T)
new.cluster.ids <- c(rep("pre_Tex",length(pre_Tex)),
                     rep("mid_T",length(mid_T)),
                     rep("term_T",length(term_T)),
                     rep("toxicity_T",length(toxicity_T)),
                     rep("other_T",length(other_T)))
cd8@meta.data$celltype2 <- plyr::mapvalues(x = as.integer(as.character(cd8@meta.data$seurat_clusters)), 
                                           from = current.cluster.ids, to = new.cluster.ids)
library(ggpubr)
library(ggsci)
x<-as.data.frame(cd8@reductions$umap@cell.embeddings) 
x$CB <-rownames(x)
y <-data.frame(cd8@meta.data[,c('CB','seurat_clusters','celltype','celltype2')])
lab <-merge(x,y,barcode='CB')
lab$seurat_clusters<-as.character(lab$seurat_clusters)
for(i in 0:8){
  lab[lab$seurat_clusters==paste("",i,sep = ""),4]<-paste("c0",i+1,sep = "")
  i+1
}
for(i in 9:13){
  lab[lab$seurat_clusters==paste("",i,sep = ""),4]<-paste("c",i+1,sep = "")
  i+1
}
lab<-arrange(lab,seurat_clusters)
lab$seurat_clusters <-factor(lab$seurat_clusters,levels=c(unique(lab$seurat_clusters)))   ##按照cluster id进行排序
library(ggsci)
color<-c("#96C3D8","#5F9BBE","#F5B375","#C0937E","#67A59B",
         "#A5D38F","#499D47","#F19294","#E45D61","#3277A9",
         "#BDA7CB","#684797","#8D75AF","#CD9C9B")#每个cluster
cell_type_med <- lab %>%
  group_by(seurat_clusters) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
library(ggrepel)
library(ggplot2)
lab$seurat_clusters<-as.character(lab$seurat_clusters)
ggplot(lab,aes(x= UMAP_1 , y = UMAP_2 ,color = seurat_clusters)) +  
  geom_point(size = 0.1 , alpha =0.3 )  +  scale_color_manual(values = color)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))+
  theme(
    legend.title = element_blank(), 
    legend.key=element_rect(fill='white'), 
    legend.text = element_text(size=20), 
    legend.key.size=unit(1,'cm') ) +  
  guides(color = guide_legend(override.aes = list(size=5)))+ 
  geom_label_repel(aes(label=seurat_clusters), fontface="bold",data = cell_type_med,
                   point.padding=unit(0.5, "lines"))+
  theme(legend.position = "none")



#  Figure 3A
a<-readRDS("texp.Rds")
cluster<-read.table("k3_class.txt",sep = "\t",header = T)
cluster<-arrange(cluster,Cluster)
cluster$Cluster<-as.factor(cluster$Cluster)
a<-a[,match(cluster$sample,colnames(a))]
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
library(ggsci)
annotation_col<-cluster[,2]%>%as.data.frame()
rownames(annotation_col)<-colnames(a)
colnames(annotation_col)<-"Cluster"
annotation_col$Cluster<-as.factor(annotation_col$Cluster)
ann_colors = list(Cluster=c(C1="#F6C097",C2="#7C94BE",C3="#ea7070"))
pheatmap(a, 
         annotation_col =  annotation_col,
         annotation_colors = ann_colors,
         cluster_rows = F, cluster_cols = F,
         show_rownames = T,show_colnames = F,
         scale="row",
         fontsize = 8,
         cutree_cols = 3,
         color = c(colorRampPalette(colors = c("#00818A","#F6F6F6"))(length(bk)/2),
                   colorRampPalette(colors = c("#F6F6F6","#F9CE00"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk)


#  Figure 3B
library(tidyr)
sample.1<-colnames(a)
a <- cbind(rownames(a),a)
colnames(a)[1]<-"pathway"
a<-as.data.frame(a)
b<-pivot_longer(a,cols = c(sample.1),names_to = c("sample"),values_to = "Score")
b<-merge(b,cluster,by="sample")
library(ggplot2)
library(ggsignif)
library(ggpubr)
compared_list = list(c("C1", "C2","C3"))  #设定比较组
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          text = element_text(size = 8),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y=element_text(face = "bold",size = 14),
          axis.text = element_text(face = "bold",size = 12),       
          axis.ticks = element_line(color='black'),
          legend.title=element_blank(),
          legend.position=c(0.7, 0.1),
          legend.direction = "horizontal",
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=8)),
          legend.background = element_rect( linetype="solid",colour ="black")
    )
}
b$Score<-as.numeric(b$Score)
b$Cluster<-as.factor(b$Cluster)
b$pathway<-factor(b$pathway,levels = c("IFN_GAMMA","TNF","IL2","CTL","M6A"))
ggplot(b, aes(x = reorder(pathway,desc(pathway)), y = Score))+ 
  geom_boxplot(aes(fill = Cluster),position=position_dodge(0.9),width=0.8)+ 
  scale_fill_manual(values = c("#F6C097","#7C94BE","#ea7070"))+theme_zg()+
  stat_compare_means(aes(group=Cluster),label = "p.signif",
                     method = "kruskal.test")+
  labs(y = "Score", x = "Pathway")


# Figure 3C
clin<-readRDS("clin.Rds")
cluster<-read.table("k3_class.txt",sep = "\t",header = T)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(stringi)
sample<-clin[,c(1,5)]
sample<-merge(cluster,sample,by="sample")
sample$Cluster<-as.factor(sample$Cluster)
sample<-arrange(sample,Cluster)
sample$project_id<-stri_sub(sample$project_id,6)
library(survival)
library(survminer)
library(ggsci)
clin_1<-clin[,c(1,8,10)]
clin_1<-merge(clin_1,cluster,by="sample")
clin_1$OS.time<-clin_1$OS.time/365
library(TSHRC)
survobj <- with(clin_1, Surv(OS.time,OS)) # 这里主要是指定时间和生存状态
fit1 <- survfit(survobj~Cluster,data=clin_1) # 在公式里指定“sex”就可以比较性别差异了
summary(fit1)
ggsurvplot(fit1,data=clin_1,
           surv.median.line = "hv",
           pval=T,pval.method=T,
           risk.table = T,
           risk.table.col="strata",
           xlab = "Survival time(years)",
           legend.labs = c("C1", "C2","C3"), 
           palette = c(C1="#F6C097",C2="#7C94BE",C3="#ea7070"))











# Figure 3D, F
clin<-readRDS("clin.Rds")#临床数据
cluster<-read.table("k3_class.txt",sep = "\t",header = T)#分型
library(ggplot2)
library(dplyr)
library(stringi)
library(openxlsx)
library(tidyr)
sample<-clin[,c(1,5)]
sample<-merge(cluster,sample,by="sample")
sample$Cluster<-as.factor(sample$Cluster)
sample<-arrange(sample,Cluster)
sample$project_id<-stri_sub(sample$project_id,6)
table(sample$Cluster)
sample.1<-sample[sample$Cluster=="C1",]%>%
  group_by(project_id)%>%
  dplyr::summarise(n=(dplyr::n()/2242)*100)
sample.2<-sample[sample$Cluster=="C2",]%>%
  group_by(project_id)%>%
  dplyr::summarise(n=(dplyr::n()/3153)*100)
sample.3<-sample[sample$Cluster=="C3",]%>%
  group_by(project_id)%>%
  dplyr::summarise(n=(dplyr::n()/4092)*100)
data<-data.frame(class=c(rep("C1",29),rep("C2",29),rep("C3",27)))
sample.4<-rbind(sample.1,sample.2,sample.3)
data<-cbind(data,sample.4)
rm(sample.1,sample.2,sample.3,sample.4)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
colnames(data)<-c("Cluster","Cancer_Type","Proportion")
ggplot(data=data,aes(Cluster,Proportion,fill=Cancer_Type))+
  geom_bar(stat="identity",position="fill", color="black", width=0.7,size=0.25)+
  scale_fill_manual(values = c(ACC="#c75250",BLCA="#cddaef",BRCA="#766955",CESC="#e2c9a1",CHOL="#717777",
                               COAD="#d8e7e7",ESCA="#61bb90",GBM="#aed049",HNSC="#eda815",KICH="#6bc4cb",
                               KIRC="#e34f32",KIRP="#ddba1a",LGG="#ce2115",LIHC="#3065af",LUAD="#555656",
                               LUSC="#d8e7d3",MESO="#4fb133",OV="#d88997",PAAD="#c6e5e6",PRAD="#78bf6d",
                               READ="#4d51a0",SARC="#824f9b",SKCM="#c32b89",STAD="#e0bd96",TGCT="#4c509f",
                               THCA="#7bb9de",THYM="#e76d2f",UCEC="#df9044",UCS="#766b4f"))+
  labs(x = "",y = "Proportion%")+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25), 
        axis.line=element_line(colour="black",size=0.25),
        axis.title=element_text(size=13,color="black"), 
        axis.text = element_text(size=12,color="black") 
  )
sample.1<-data.frame()
Cancer_Type<-c()
for(i in unique(sample$project_id)){
  sample.2<-sample[sample$project_id==i,]%>%group_by(Cluster)
  sample.2<-dplyr::summarise(sample.2,n=(dplyr::n()/nrow(sample.2)))
  sample.1<-rbind(sample.1,sample.2)
  Cancer<-c(rep(i,nrow(sample.2)))
  Cancer_Type<-c(Cancer_Type,Cancer)
}
data<-cbind(sample.1,Cancer_Type)
rm(sample.1,sample.2,Cancer,Cancer_Type)
colnames(data)[2]<-"Proportion"
data$Cancer_Type<-factor(data$Cancer_Type)
ggplot(data=data,aes(Cancer_Type,Proportion,fill=Cluster))+
  geom_bar(aes(fill = Cluster), position = position_stack(reverse = TRUE),stat="identity", color="black", width=0.7,size=0.25)+
  scale_fill_manual(values = c(C1="#ce524b",C2="#ecbf00",C3="#0271c2"))+
  labs(x = "",y = "Proportion%")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(limits=rev(levels(data$Cancer_Type)))+
  theme_classic()+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25),
        axis.line=element_line(colour="black",size=0.25), 
        axis.title=element_text(size=13,color="black"), 
        axis.text = element_text(size=12,color="black"), 
        legend.position = "top"#legend.position= c(3,2) 
  )+coord_flip()


# Figure 4A
C1_high<-read.table("C1_high.txt",sep = "\t",header = T)
C2_high<-read.table("C2_high.txt",sep = "\t",header = T)
C3_high<-read.table("C3_high.txt",sep = "\t",header = T)
C1_high$label<-c(rep("TexHm6AH_Special",170))
C2_high$label<-c(rep("TexLm6AL_Special",88))
C3_high$label<-c(rep("TexLm6AH_Special",499))
library(data.table)
library(dplyr)
load("pancancer_exp.Rdata")
mRNA<-fread("mRNA_GRCh38.104.gtf.TXT")
mRNA<-mRNA[,-1]
exp<-pancancer_exp[pancancer_exp$Ensembl_ID%in%mRNA$gene_id,]
rm(pancancer_exp)
colnames(mRNA)[1]<-"Ensembl_ID"
exp<-merge(mRNA,exp,by="Ensembl_ID")
exp<-exp[,-1]
C1_exp<-exp[exp$gene_name%in%C1_high$x,]
C2_exp<-exp[exp$gene_name%in%C2_high$x,]
C3_exp<-exp[exp$gene_name%in%C3_high$x,]
exp.1<-rbind(C1_exp,C2_exp,C3_exp)%>%as.data.frame()
rownames(exp.1)<-exp.1$gene_name
exp.1<-exp.1[,-1]
rm(exp,mRNA)
cluster<-read.table("k3_class.txt",sep = "\t",header = T)
cluster<-arrange(cluster,Cluster)
exp.1<-exp.1[,match(cluster$sample,colnames(exp.1))]
library(ggsci)
annotation_col<-data.frame(Cluster=cluster$Cluster)
rownames(annotation_col)<-cluster$sample
annotation_row<-data.frame(label=c(C1_high$label,C2_high$label,C3_high$label))
rownames(annotation_row)<-c(C1_high$x,C2_high$x,C3_high$x)
exp.1<-exp.1[match(rownames(annotation_row),rownames(exp.1)),]
ann_colors<-list(Cluster=c(C1="#F6C097",C2="#7C94BE",C3="#ea7070"),
                 label=c(TexHm6AH_Special="#f9eca9",TexLm6AL_Special="#87bfc9",
                         TexLm6AH_Special="#87bfc9"))
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
library(pheatmap)
pheatmap(exp.1, 
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         show_rownames = F,show_colnames = F,
         cluster_cols = T,cluster_rows = F,
         annotation_colors = ann_colors,
         scale = "row",
         fontsize = 8,
         cutree_rows = 3,
         gaps_row = c(170,258),
         cutree_cols = 3,
         gaps_col = c(2242,5395),
         legend = F,
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk)
dev.off()


# Figure 3H
library(stringi)
clin1<-readRDS("clin.Rds")
clin1$project_id<-stri_sub(clin1$project_id,6)
clin1$sample<-stri_sub(clin1$sample,1,15)
clin2<-read.csv("PANCAN.allsurvival.csv",stringsAsFactors = T)##临床数据2
clin2$X<-as.character(clin2$X)
clin2$X<-gsub("\\.","-",clin2$X)
clin2<-clin2[match(clin1$sample,clin2$X),]
colnames(clin2)[1]<-"sample"
cluster<-read.table("k3_class.txt",sep = "\t",header = T)#分型
cluster$sample<-stri_sub(cluster$sample,1,15)
clin<-merge(clin1[,c(1,5,6,7)],clin2[,c(1,6)],by="sample")
clin<-merge(clin,cluster,by="sample")
clin$Age<-ifelse(clin$Age.at.Diagnosis.in.Years<="19","young",
                 ifelse(clin$Age.at.Diagnosis.in.Years<="39","adult",
                        ifelse(clin$Age.at.Diagnosis.in.Years<="59","pre-aging","old")))
colnames(clin)[3]<-"Ages"
colnames(clin)[2]<-c("Cancer_Type")
colnames(clin)[5]<-c("Stage")
clin<-merge(clin,clin1[,c(1,8,10)],by="sample")
clin$OS.time<-clin$OS.time/365
clin$Age<-factor(clin$Age,levels = c("young","adult","pre-aging","old"))

clin$Stage<-as.character(clin$Stage)
clin$Stage<-ifelse(clin$Stage=="Stage I"|clin$Stage=="Stage IA"|clin$Stage=="Stage IB","Stage I",
                   ifelse(clin$Stage=="Stage II"|clin$Stage=="Stage IIA"|clin$Stage=="Stage IIB"|clin$Stage=="Stage IIC","Stage II",
                          ifelse(clin$Stage=="Stage III"|clin$Stage=="Stage IIIA"|clin$Stage=="Stage IIIB"|clin$Stage=="Stage IIIC","Stage III",
                                 ifelse(clin$Stage=="Stage IV"|clin$Stage=="Stage IVA"|clin$Stage=="Stage IVB"|clin$Stage=="Stage IVC","Stage IV","Unknown"))))
clin$Stage<-factor(clin$Stage,levels = c("Stage I","Stage II","Stage III","Stage IV","Unknown"))

clin$Gender<-as.character(clin$Gender)
clin$Gender<-factor(clin$Gender,levels = c("Female","Male"))

clin$Cluster<-ifelse(clin$Cluster=="C1","TexHm6AH",
                     ifelse(clin$Cluster=="C2","TexLm6AL","TexLm6AH"))
clin$Cluster<-factor(clin$Cluster,levels = c("TexLm6AL","TexHm6AH","TexLm6AH"))

clin$Stage2<-ifelse(clin$Stage=="Stage I"|clin$Stage=="Stage II","Stage I/II",
                    ifelse(clin$Stage=="Stage III"|clin$Stage=="Stage IV","Stage III/IV","Uknown"))
library(tableone)  
library(survival)
library(forestplot)
library(stringr)
precox<-clin
result1=data.frame()
index<-which(precox$Stage2=="Uknown")
precox<-precox[-index,]
precox$Stage2<-as.factor(precox$Stage2)
for(i in colnames(precox[,c(10,4,3)])){
  cox <- coxph(Surv(OS.time, OS) ~ get(i), data = precox)
  coxSummary = summary(cox)
  multi1<-as.data.frame(round(coxSummary$conf.int[, c(1, 3, 4)],3))%>%t()
  rownames(multi1)<-i
  multi2<-ShowRegTable(cox, 
                       exp=TRUE, 
                       digits=3, 
                       pDigits =3,
                       printToggle = TRUE, 
                       quote=FALSE, 
                       ciFun=confint)
  result<-cbind(multi1,multi2)
  result1<-rbind(result1,result)
}
cox <- coxph(Surv(OS.time, OS) ~ Cluster, data = precox)
coxSummary = summary(cox)
multi1<-as.data.frame(round(coxSummary$conf.int[, c(1, 3, 4)],3))
multi2<-ShowRegTable(cox, 
                     exp=TRUE, 
                     digits=3, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
result<-cbind(multi1,multi2)
result1<-rbind(result1,result)
result1<-tibble::rownames_to_column(result1, var = "Characteristics")
result1
result1$Characteristics<-str_remove(result1$Characteristics,"s|2|Cluster|Cluster")
ins <- function(x) {c(x, rep(NA, ncol(result)))}
for(i in 5:6) {result1[, i] = as.character(result1[, i])}
result1<-rbind(c("Characteristics", NA, NA, NA, "HR(95%CI)","P-value"),
               result1[1, ],
               ins("Stage III/IV vs I/II"),
               result1[2,],
               ins("Male vs Female"),
               ins("Cluster:"),
               result1[4:5,],
               result1[3,])
for(i in 2:4) {result1[, i] = as.numeric(result1[, i])}
fig<-forestplot(result1[,c(1,5,6)], 
                mean=result1[,2],
                lower=result1[,3],
                upper=result1[,4],
                zero=1,
                boxsize=0.1,
                graph.pos= 2 ,

                graphwidth = unit(.4,"npc"),
                xlab="",
                xticks=c(1,1.5,2,2.5) ,
                txt_gp=fpTxtGp(
                  label=gpar(cex=1),
                  ticks=gpar(cex=1),
                  xlab=gpar(cex=1.5),
                  title=gpar(cex=1)),
                lwd.zero=1,
                lwd.ci=0.5,
                lwd.xaxis=1,
                lty.ci=1.5,
                ci.vertices =T,
                ci.vertices.height=0.1,
                clip=c(0.1,2.5),
                ineheight=unit(8, 'mm'),
                line.margin=unit(8, 'mm'),
                colgap=unit(6, 'mm'),
                fn.ci_norm="fpDrawNormalCI",
                mar=unit(rep(1.25,times=4),"cm"),
                title="Univariate Cox",
                col=fpColors(box ='black', 
                             lines ='black', 
                             zero = "black",
                ))
fig


# Figure 4E
library(org.Hs.eg.db)
library(clusterProfiler)
deg<-read.table("C1vsC2_all.txt",header = T)
deg<-read.table("C1vsC3_all.txt",header = T)
deg<-read.table("C2vsC3_all.txt",header = T)
gene<-deg$GENEs
gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
colnames(deg)[4]<-"SYMBOL"
deg<-deg[,c(3,4)]
gene<-merge(deg,gene,by="SYMBOL")
library(dplyr)
gene<-arrange(gene,desc(logFC))
geneList<-gene$logFC
names(geneList)<-gene$ENTREZID
hallmarker_gmt <- read.gmt("h.all.v7.5.1.entrez.gmt") #读gmt文件
gsea <- GSEA(geneList,TERM2GENE = hallmarker_gmt,
             pvalueCutoff = 1,minGSSize = 0,nPermSimple = 2000000) #GSEA分析
gsea <- as.data.frame(gsea)
name<-read.table("hallmark_name_50.txt")
name<-name$x
data<-gsea[,c(1,5,7)]
data$col<-ifelse(data$NES>0&data$p.adjust<0.05,'TexHm6AH',
                 ifelse(data$NES<0&data$p.adjust<0.05,'TexLm6AL',"normal"))
data$col<-ifelse(data$NES>0&data$p.adjust<0.05,'TexHm6AH',
                 ifelse(data$NES<0&data$p.adjust<0.05,'TexLm6AH',"normal"))
data$col<-ifelse(data$NES>0&data$p.adjust<0.05,'TexLm6AH',
                 ifelse(data$NES<0&data$p.adjust<0.05,'TexLm6AL',"normal"))
data<-data[match(name,data$ID),]
library(ggplot2)
ggplot(data = data,aes(x = ID, y = NES))+
  coord_flip()+
  geom_bar(stat = 'identity',
           fill = ifelse(data$col=="TexLm6AH","#DE6C6D",
                         ifelse(data$col=="TexLm6AL","#788EB6","#dadada")),
           width = 0.8)+  
  scale_x_discrete(limits=rev(data$ID))+#坐标轴顺序确定
  theme_classic()+
  theme(panel.background=element_rect(fill="white",colour="black",size=0.25), 
        axis.line=element_line(colour="black",size=0.25), 
        axis.title=element_text(size=12,color="black"), 
        axis.text = element_text(size=7,color="black") 
  )


# Figure 4F
Fges<-read.csv("29Fges_marker_PMID_34951955.csv")
load("pancancer_exp.Rdata")
mRNA<-fread("mRNA_GRCh38.104.gtf.TXT")
mRNA<-mRNA[,-1]
exp<-pancancer_exp[pancancer_exp$Ensembl_ID%in%mRNA$gene_id,]
rm(pancancer_exp)
exp$Ensembl_ID<-mRNA[match(exp$Ensembl_ID,mRNA$gene_id),2]
gene<-exp$Ensembl_ID
exp<-exp[,-1]
exp<-as.matrix(exp)
rownames(exp)<-gene$gene_name
rm(mRNA)
geneSets<-c()
Fges_p<-Fges$name
for(i in 1:nrow(Fges)){
  df<-as.data.frame(t(Fges[i,]))
  df<-df[!apply(df == "", 1, all),]
  list<-list(as.character(df[3:length(df)]))
  names(list)<-Fges_p[i]
  geneSets <- c(geneSets,list)
}
a<-gsva(exp,geneSets,method="ssgsea",kcdf="Gaussian")
cluster<-read.table("k3_class.txt",sep = "\t",header = T)
cluster<-arrange(cluster,Cluster)
a<-a[,match(cluster$sample,colnames(a))]
cluster$Cluster<-ifelse(cluster$Cluster=="C1","TexHm6AH",
                        ifelse(cluster$Cluster=="C2","TexLm6AL","TexLm6AH"))
protumor<-c(Fges_p[9:16])
antitummor<-c(Fges_p[17:27])
fibroblast<-c(Fges_p[1:8])
other<-c(Fges_p[28:29])
annotation_row<-data.frame(Pathway=c(protumor,antitummor,fibroblast,other),
                           Type=c(rep("protumor",length(protumor)),rep("antitummor",length(antitummor)),
                                  rep("fibroblast",length(fibroblast)),rep("other",length(other))))
a<-a[match(annotation_row$Pathway,rownames(a)),]
annotation_row<-annotation_row$Type%>%as.data.frame()
colnames(annotation_row)<-"Immune_type"
rownames(annotation_row)<-rownames(a)
annotation_col<-cluster[,2]%>%as.data.frame()
rownames(annotation_col)<-colnames(a)
colnames(annotation_col)<-"Cluster"
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))#热图表格标签颜色间断
library(ggsci)
ann_colors = list(Cluster=c(TexHm6AH="#E8B891",TexLm6AL="#7488AE",TexLm6AH="#D2696A"),
                  Immune_type=c(protumor="#f17d80",antitummor="#737495",
                                fibroblast="#68a8ad",other="#c4d4af"))
pheatmap(a, 
         annotation_col =  annotation_col,
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         cluster_rows = F, cluster_cols = F,
         show_rownames = T,show_colnames = F,
         scale="row",
         fontsize = 8,
         gaps_col = c(2242,5395),
         cutree_cols = 3,
         gaps_row = c(8,19,27),
         cutree_rows = 2,
         color = c(colorRampPalette(colors = c("#0000A1","white"))(length(bk)/2),
                   colorRampPalette(colors = c("white","#f1404b"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),
         breaks=bk)


# Figure 5A
setwd("D:\\TEX_FY\\中间数据\\免疫治疗")
k_class<-read.table("Gide_k3.txt",header = T)
clin<-read.csv("Gide_clinical_all.csv")
clin<-clin[clin$Patient%in%k_class$Patient,]
clin<-clin[,c(1,11,12)]
clin$OS_time<-clin$OS_time/365
clin<-na.omit(clin)
clin_1<-merge(k_class,clin,by="Patient")
library(survival)
library(survminer)
library(ggsci)
library(TSHRC)
survobj <- with(clin_1, Surv(OS_time,OS))
fit1 <- survfit(survobj~Cluster,data=clin_1) 
ggsurvplot(fit1,data=clin_1,
           pval=T,pval.method=T,
           risk.table = T,
           risk.table.col="strata",
           xlab = "Survival time(years)",
           legend.labs = c("C1", "C2","C3"), 
           palette = c(C1="#D25565",C2="#3DC7BE",C3="#2E94B9")
)



# Figure 5C
class1<-read.table("Au_k3.txt",header = T)
clin<-read.table("Au_standardization_clinical_infor.txt",sep = "\t",header = T)
clin<-clin[,c(1,12,13)]
clin<-na.omit(clin)
clin1<-merge(class1,clin,by="Patient")
class2<-read.table("Hu_k3.txt",header = T)
clin<-read.table("Hu_standardization_clinical_infor.txt",sep = "\t",header = T)
clin<-clin[,c(1,8,9)]
clin<-na.omit(clin)
clin2<-merge(class2,clin,by="Patient")
clin2<-clin2[,c(1,2,4,3)]
clin<-rbind(clin1,clin2)
clin$OS_time<-clin$OS_time/365
survobj <- with(clin, Surv(OS_time,OS)) 
fit1 <- survfit(survobj~Cluster,data=clin) 
ggsurvplot(fit1,data=clin,
           pval=T,pval.method=T,
           risk.table = T,
           risk.table.col="strata",
           xlab = "Survival time(years)",
           legend.labs = c("C1", "C2","C3"), 
           palette = c(C1="#D25565",C2="#3DC7BE",C3="#2E94B9")
)






















