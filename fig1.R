###Fig1
###
###

#####
###Fig1 B
#####
sc = read.csv("~/data/all_sc.csv",row.names = 1) # loading coordinate of cells
sc = merge(sc,metatable,by = "Well_ID")

show_col(pal_npg("nrc")(10))

major = c("DC","macrophage","monocyte","progenitor","granulocyte","MK","ILC","T","NK","B","erythrocyte")
major_color = c("#E64B35FF","#F39B7FFF","#B24745FF","#8491B4FF","#7E6148FF","#B09C85FF","#95CC5EFF","#00A087FF", "#91D1C2FF","#F79D1EFF","#69C8ECFF")
names(major_color) = major
save(major_color,file = "~/data/major_color.Rdata") # save color scheme for major lineages
sc = sc[sc$time != "Adult",]
p1 = ggplot(sc)+
  geom_point(aes(x = x,y = y,fill = major),shape = 21,color = "grey40",size = 3,stroke = 0.2)+
  theme_bw()+
  xlab("")+ylab("")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        legend.position = "none")+
  scale_fill_manual(values = as.character(major_color),breaks = names(major_color))
ggsave("~[folders]/2d_noadult.pdf",p1,width = 18,height = 18,dpi = 300) # save the image


#####
###Fig1 C
#####

load("~/data/metadatable.rdata",row.names = 1) # loading single cell metatable 

metatable2 = metatable
metatable2$major[metatable2$subtype == "microglial"] = "microglia"

tmp = metatable2$major
names(tmp) = metatable2$Well_ID

save(tmp,file = "~/data/lfp_major.Rdata")


metatable2 = metatable2[metatable2$major != "erythrocyte",]

marker = read.csv("~/data/all_marker.csv")
marker = marker$gene


a_reduced?

use.gene = c("HMGA2","CD34","MYB","MPO","LMO4","CPA3","TPSAB1","GATA2","PF4","GP9","GP1BB",
             "S100A9","CLEC10A","CD1C","HLA-DPB1","HLA-DRA","CD163","MRC1","C1QA","LYVE1","DAB2","P2RY12","C3","CD79B",
             "IGHM","VPREB3","NKG7","GNLY","GZMA","CD3D","KLRB1","TRBC2","GATA3","IL7R")
marker = use.gene
tmp1 = a_reduced[marker,intersect(colnames(a_reduced),metatable$Well_ID)]
tmp1 = as.matrix(tmp1)

gene_exp = as.data.frame(t(tmp1))

gene_exp$wellid = rownames(gene_exp)

tmp = metatable2[,c("Well_ID","major")]
colnames(tmp) = c("wellid","subtype")

gene_exp = merge(gene_exp,tmp,by = "wellid",all.x = T)
gene_exp = gene_exp[,-1]

gene_exp = melt(gene_exp)

gene_exp2 = gene_exp[,c("subtype","variable")]
gene_exp2 = gene_exp2[!duplicated(gene_exp2),]

gene_exp2$mean = 100
gene_exp2$ratio = 100
gene_exp2 = na.omit(gene_exp2)
gene_exp[is.na(gene_exp)] = 0
gene_exp2 =  gene_exp2[gene_exp2$variable %in% use.gene,]
for (i in unique(gene_exp2$subtype)){
  print(i)
  for (j in unique(gene_exp2$variable[gene_exp2$subtype == i])) {
    gene_exp2$mean[gene_exp2$subtype == i & gene_exp2$variable == j] = mean(gene_exp$value[gene_exp$subtype == i & gene_exp$variable == j])
    gene_exp2$ratio[gene_exp2$subtype == i & gene_exp2$variable == j] = nrow(gene_exp[gene_exp$subtype == i & gene_exp$variable == j & gene_exp$value >= 2,])/nrow(gene_exp[gene_exp$subtype == i & gene_exp$variable == j,])
  }
  
}

lfp_major = read.csv("~/data/lfp_major.csv")
colnames(lfp_major)[1] = "gene"
lfp_major = reshape2::melt(lfp_major)
lfp_major = as.data.frame(lfp_major)

gene_exp2 = as.data.frame(gene_exp2)
gene_exp2$lfp = NA
for (i in unique(gene_exp2$variable)){
  for (j in unique(gene_exp2$subtype)) {
    gene_exp2$lfp[gene_exp2$subtype == j & gene_exp2$variable == i] = lfp_major$value[lfp_major$gene == i & lfp_major$variable == j]
  }
}

gene_exp3 = gene_exp2
gene_exp3$subtype = factor(gene_exp3$subtype,levels = c("ILC","T","NK","B","microglia","macrophage",
                                                        "DC","monocyte","MK","granulocyte","progenitor"))
gene_exp3$variable = factor(gene_exp3$variable,levels = use.gene,ordered = T)
gene_exp3$mean = log2(gene_exp3$mean+1)
gene_exp3$mean[gene_exp3$mean >2] = 2
gene_exp3$lfp[gene_exp3$lfp >3] = 3
gene_exp3$lfp[gene_exp3$lfp < 0] = 0

p1 = ggplot(gene_exp3,aes(x = variable , y= subtype,fill = lfp, size =ratio))+
  geom_point(shape = 21,color = "grey30" )+
  scale_fill_distiller(palette = "Reds",direction = 1)+
  #scale_fill_gradient(low = "blue",high = "red")+
  xlab("")+ylab("")+
  theme_bw()+
  theme(axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1))
ggsave("~[folders]/marker_newnormalize.pdf",p1,width = 8,height = 4) # save the image


#####
###Fig1 D
#####

load("~/data/metatable_0929.Rdata") # loading single cell metatable 

the_metatable = metatable
the_metatable = the_metatable[the_metatable$major != "erythrocyte",]
the_metatable$new = paste0(the_metatable$tissue,the_metatable$major)
for (new in unique(the_metatable$new)) {
  if (nrow(the_metatable[the_metatable$new == new,])<10){
    the_metatable= the_metatable[-which(the_metatable$new == new),]
  }
}

the_metatable = the_metatable[,c("time","major","Well_ID","tissue")]
the_metatable = the_metatable[-which(the_metatable$time == "Adult"),]
the_metatable = the_metatable[-which(the_metatable$time == "week19"),]
the_metatable$time <- factor(the_metatable$time,levels=time_order,ordered = TRUE)

tmp2 = the_metatable[,c("time","major")]

tmp2 = tmp2[!duplicated(tmp2),]
rownames(tmp2) = NULL
tmp2$rate = 0
for(t in rownames(tmp2)){
  tmp2[t,"rate"] = nrow(the_metatable[the_metatable$time == tmp2[t,"time"] & the_metatable$major == tmp2[t,"major"],])/nrow(the_metatable[the_metatable$time == tmp2[t,"time"],])
}

p1 = ggplot(tmp2,aes(x = time, 
                     stratum = major, 
                     alluvium = major,
                     y = rate,
                     fill = major,
                     label = major)) +
  scale_y_continuous(label = scales::percent_format(),
                     expand=c(0,0)) +
  scale_fill_manual(values = as.character(major_color),breaks = names(major_color)) +
  geom_flow() +
  geom_stratum(width = 0.6,colour = "grey40") +
  theme(axis.title=element_blank(),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.title=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(),
        axis.ticks=element_line(),
        legend.position = "none",
        axis.text = element_text(size = 12))
ggsave("E:\\fig1/ratio_legend.pdf",p1,width =4,height =3,dpi = 300)

###
figre2_major_ratio_per_tissue = function(the_metatable,the_dir,thew , theh){
  tmp = the_metatable
  tmp2 = tmp
  
  for(i in unique(tmp2$tissue)){
    tmp = tmp2[,c("tissue","time","subtype")]
    tmp = tmp[tmp$tissue == i,]
    tmp = tmp[,c("time","subtype")]
    other = setdiff(time_order,c(unique(tmp$time),"Adult","week19"))
    other = as.data.frame(other)
    other$blank = "blank"
    colnames(other) = colnames(tmp)
    tmp = rbind(tmp,other)
    tmp$time <- factor(tmp$time,levels=time_order,ordered = TRUE)
    p2 = ggplot(tmp,aes(x= time,fill = subtype))+
      geom_bar(stat="count",position = "fill",width = 0.8,color = "black")+
      xlab("")+ylab("")+ggtitle(paste(i," total cell number ",nrow(tmp),sep = ""))+
      scale_fill_manual(values = c(as.character(mac.color),"white"),breaks = c(names(mac.color),"blank"))+
      scale_y_continuous(expand=c(0,0))+
      theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
            axis.text = element_text(size = 12),
            axis.title=element_blank(),
            legend.title=element_blank(),
            panel.background=element_blank(),
            axis.line=element_line(),
            axis.ticks=element_line(),
            legend.position = "none")
    ggsave(paste(the_dir,i,"_tissue.pdf",sep = ""),p2,width = thew,height = theh) # save the image
  }
}

#####
###Fig1 E
#####

library(ggsignif)
library(ggplot2)
source("~/Funcs_signif_fig2.R")
setwd("./path")
#load("mac_color.Rdata")
load("~/major_color.Rdata") # loading color scheme for major lineages
load("~/data/metatable_0929.Rdata") # loading single cell metatable 
#######################################################################################
# major cell type proportion significance test
metatable <- metatable %>% filter(time!="Adult")
metatable <- metatable %>% filter(major!="erythrocyte") # remove erythrocyte
metatable <- metatable %>% filter(embryo!="embryo 10") # embryo 10 contain 1 cell only
table(metatable$major) %>% sort
metatable %>% filter(subtype=="discard") %>% pull(major) %>% table 
metatable <- metatable %>% filter(subtype!="discard") # drop discard annotation
metatable$major <- as.factor(metatable$major)

### divide all time into 3 periods
CS_stage = c("cs11","cs12","cs13","cs14","cs18","cs19","cs21","cs23")
`9-13 PCW` = c("week9","week10","week11","week12","week13")
`After 16 PCW` = c("week16","week19","week20","week23","week27")
metatable$period <- metatable$time
metatable$period[metatable$period%in%CS_stage] <- "CS_stage"
metatable$period[metatable$period%in%`9-13 PCW`] <- "9-13 PCW"
metatable$period[metatable$period%in%`After 16 PCW`] <- "After 16 PCW"
### filter out adult cells
meta_fetal <- metatable %>% filter(time!="Adult")
table(meta_fetal$period)
### compute ratio of subtypes in different time periods
ratio <- list()
for(i in unique(meta_fetal$embryo)){
  df <- meta_fetal %>% filter(embryo==i)
  ratio[[i]] <- tapply(df$major,df$period,function(x){prop.table(table(x))}) %>% do.call(cbind,.) %>% as.data.frame()
  ratio[[i]] <- ratio[[i]] %>% mutate(major=rownames(.),period=colnames(.),embryo=i)
  colnames(ratio[[i]])[1] <- "ratio"
}
ratio <- do.call(rbind,ratio)
ratio$period <- factor(ratio$period,levels = c("CS_stage","9-13 PCW","After 16 PCW"))
### draw median horizontal line for each group; control segment length
median <- tapply(ratio$ratio,list(ratio$major,ratio$period),median) %>% melt()
colnames(median) <- c("major","period","median")
median$x_start <- as.numeric(median$period)-0.3
median$x_end <- as.numeric(median$period)+0.3
my_comp <- list(c("CS_stage","9-13 PCW"),c("9-13 PCW","After 16 PCW"),c("CS_stage","After 16 PCW"))
color_grp <- c("#A6A8A9FF","#484A4BFF","#141414FF")
names(color_grp) <- sort(unique(ratio$period))

plt_Frac_signif(ratio, "./signif/") # save the image




#####
###Fig1 F
#####

load("`/data/metatable_0929.Rdata") # loading single cell metatable 


data1<-metatable

##############################
# remove Adult

data1<-data1[!data1$time%in%c("Adult"),]

majorType<-unique(data1$major)

majoryType.cellsum<-vector()

for(i in 1:length(majorType)){
  
  majoryType.cellsum[i]<-nrow(subset(data1,major==majorType[i]))
  
  
}


names(majoryType.cellsum)<-majorType


#####################################

minorType<-unique(data1$subtype)

minorType.cellsum<-vector()

for(i in 1:length(minorType)){
  
  minorType.cellsum[i]<-nrow(subset(data1,subtype==minorType[i]))
  
  
}

names(minorType.cellsum)<-minorType

tissue.num<-length(unique(data1$tissue))
tissue<-unique(data1$tissue)

tissue_order<-read.csv("~/data/tissue_order.csv",header=T)

subtype_order<-read.csv("~/data/subtypes_order.csv",header=T)

heatmap.data<-matrix(data=NA, ncol = tissue.num , nrow =length(minorType), byrow = FALSE, dimnames = NULL)


colnames(heatmap.data)<-as.matrix(tissue_order)
rownames(heatmap.data)<-as.matrix(subtype_order$subtype)
minorType_order<-subtype_order$subtype
tissue_order<-as.character(as.matrix(tissue_order))


for(i in 1:tissue.num){
  
  for(j in 1:length(minorType) ){
    
    
    numerator = nrow(subset(data1,subtype == minorType_order[j] & tissue == tissue_order[i]))
    
    denominator = nrow(subset(data1,subtype == minorType_order[j]))
    
    ratio = numerator/denominator
    
    heatmap.data[j,i]<-ratio
    
  }
}

####

heatmap.data<-t(heatmap.data)

annotation_col<-subtype_order
rownames(annotation_col)<-annotation_col$subtype
annotation_col<-subset(annotation_col,select = -subtype)

ann_colors = list(
  major = major_color
)

bk <- c(seq(0,0.2,by=0.01),seq(0.21,0.6,by=0.01))
p<-pheatmap(heatmap.data,breaks = bk,gaps_col = c(6,9,12,18,23,29,32,35,41),
            color = c(colorRampPalette(colors = c("white","firebrick"))(length(bk)/2),colorRampPalette(colors = c("firebrick","Maroon"))(length(bk)/2)),show_colnames=T,display_numbers = F,cluster_cols =F,cluster_rows = F,fontface="bold",annotation_colors = ann_colors,annotation_col = annotation_col)


ggsave(p,filename = "~[folders]/fig1F.pdf",width = 20, height = 8, units = "in", device='pdf') # save the image

