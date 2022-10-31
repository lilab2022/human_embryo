####Fig3
####write by Feng Ruoqing & Wang Hao & Li Muxi
####Date 2022-10-4

##Fig3B 

marker.violin <- c("RAN","S100B","ACY3",
                   "AFP", "TTR","STMN1","AGR2","AXL",
                   "MRC1",
                   "LIPA","HMOX1","TIMD4","CD5L",
                   "RNASE1","F13A1",
                   "IL1B","CD83","CXCL8","C5AR1",
                   "MMP9",
                   "CD74","HLA-DQA1",
                   "MS4A6A",
                   "DNASE1L3","CD209","ADAMDEC1",
                   "CD207","CD1A","PTGS2",
                   "SIGLEC15","ACP5",
                   "C3","P2RY12")
marker = marker.violin
tmp1 = a_reduced[marker,intersect(colnames(a_reduced),mac.metatable$Well_ID)]
tmp1 = as.matrix(tmp1)

gene_exp = as.data.frame(t(tmp1))

gene_exp$wellid = rownames(gene_exp)

tmp = mac.metatable[,c("Well_ID","subtype")]
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
gene_exp2 =  gene_exp2[gene_exp2$variable %in% marker.violin,]
for (i in unique(gene_exp2$subtype)){
  print(i)
  for (j in unique(gene_exp2$variable[gene_exp2$subtype == i])) {
    gene_exp2$mean[gene_exp2$subtype == i & gene_exp2$variable == j] = mean(gene_exp$value[gene_exp$subtype == i & gene_exp$variable == j])
    gene_exp2$ratio[gene_exp2$subtype == i & gene_exp2$variable == j] = nrow(gene_exp[gene_exp$subtype == i & gene_exp$variable == j & gene_exp$value >= 3,])/nrow(gene_exp[gene_exp$subtype == i & gene_exp$variable == j,])
  }
  
}

lfp_major = read.csv("E:\\mac_new/subtype_lfp.csv",row.names = 1)
colnames(lfp_major) = mac.sudomc$subtype
lfp_major$gene = rownames(lfp_major)
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
gene_exp3$subtype = factor(gene_exp3$subtype,levels = rev(c("HdpM","YsdM_AFP_high","YsdM_AFP_low","pre_microglia",
                                                            "red_pulp","Kupffer_cell","pre_PraM","PraM","gonad_macrophage",
                                                            "Adrenalgland_macrophage","intestine CD209+ Mφ",
                                                            "intestine CD207+ Mφ","langerhans","osteoclast","microglia")))
gene_exp3$variable = factor(gene_exp3$variable,levels = marker.violin,ordered = T)
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
ggsave("E:\\mac_new/mac_marker_new_normalize.pdf",p1,width = 9,height =4)

###Fig3 C

#distance

# Define the distance formula, such as Euclidean distance
euclidean_dist <- function(x, y) sqrt(sum((x - y)^2))



data<-read.csv("G://embryo/20220801/lfp_mf_high_v2.csv")
rownames(data)<-data$X
data<-data[,-1]

#data<-scale(data)


## meta information 
load("G://embryo/20220929/metatable_0929.Rdata")
load("G://embryo/20220929/mac.color2.rdata")

mc.vs.subtype<-metatable %>% group_by(mc) %>% filter (! duplicated(mc))


mc.vs.subtype.mf<-mc.vs.subtype[mc.vs.subtype$mc %in% colnames(data),]

######################################################

## calculate centroid 

subtype.set<-unique(mc.vs.subtype.mf$subtype)

centroid.matrix<-matrix(data=NA,nrow=length(unique(mc.vs.subtype.mf$subtype)),ncol=nrow(data),byrow=FALSE,dimnames=NULL)
colnames(centroid.matrix)<-rownames(data)
rownames(centroid.matrix)<-subtype.set

for(i in 1:length(unique(mc.vs.subtype.mf$subtype))){
  
  
  tmp<-mc.vs.subtype.mf[mc.vs.subtype.mf$subtype%in%subtype.set[i],]
  
  data.tmp<-data[,colnames(data)%in%tmp$mc]
  
  tmp.coor2<-rowMeans(data.tmp)
  
  
  centroid.matrix[i,]<-tmp.coor2
  
  
}

########

dist.matrix<-matrix(data=NA,nrow=nrow(centroid.matrix),ncol=nrow(centroid.matrix),byrow=FALSE,dimnames=NULL)

rownames(dist.matrix)<-rownames(centroid.matrix)
colnames(dist.matrix)<-rownames(centroid.matrix)



for(i in 1:nrow(centroid.matrix)){
  
  for(j in 1:nrow(centroid.matrix)){
    
    x<-centroid.matrix[i,]
    y<-centroid.matrix[j,]
    
    dist.matrix[i,j]<-euclidean_dist(x,y)
    
    
  }
  
  
}

p2<-pheatmap(dist.matrix,color =colorRampPalette(c("#1B519C", "white", "#A00021"))(50),
             display_numbers = F,angle_col = 45,fontsize_row = 12, 
             fontsize_col = 12,cluster_rows = T, cluster_cols = T,cutree_rows = 6,annotation_colors = ann_colors,
             show_colnames = F,annotation_row = annotation_col, treeheight_col=100,legend=F
)




ggsave(p2,filename="G://embryo/figures/fig3/mf_subtype_distance_20220930.pdf",width=15,height=10,units="in",device='pdf')
ggsave(p2,filename="G://embryo/figures/fig3/mf_subtype_distance_20220930.png",width=15,height=10,units="in",device='png')


###Fig3E
load("G://embryo/20220929/metatable_0929.Rdata")

load("G://embryo/20220929/mac.color2.rdata")


data1<-metatable

data1<-data1[!data1$time%in%c("Adult"),]



majorType<-unique(data1$major)
majoryType.cellsum<-vector()
for(i in 1:length(majorType)){
  
  majoryType.cellsum[i]<-nrow(subset(data1,major==majorType[i]))
  
  
}

names(majoryType.cellsum)<-majorType
data1<-data1[data1$major%in%"macrophage",]

minorType<-unique(data1$subtype)

minorType.cellsum<-vector()

for(i in 1:length(minorType)){
  
  minorType.cellsum[i]<-nrow(subset(data1,subtype==minorType[i]))
  
  
}

names(minorType.cellsum)<-minorType


tissue_order<-read.csv("G://embryo/20220822/tissue_order2.csv",header=T)
subtype_order<-read.csv("G://embryo/20220929/subtypes_new2_mf.csv",header=T)


tissue.num<-length(tissue_order$tissue)
tissue<-unique(data1$tissue)


data1<-data1[data1$tissue%in%tissue_order$tissue,]

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

heatmap.data<-t(heatmap.data)

annotation_col<-subtype_order
rownames(annotation_col)<-annotation_col$subtype
annotation_col<-subset(annotation_col,select = -major)

ann_colors = list(
  subtype = mac.color
)

bk <- c(seq(0,0.2,by=0.01),seq(0.21,0.6,by=0.01))
p<-pheatmap(heatmap.data,breaks = bk,
            color = c(colorRampPalette(colors = c("white","firebrick"))(length(bk)/2),colorRampPalette(colors = c("firebrick","Maroon"))(length(bk)/2)),
            show_rownames = T,show_colnames=F,display_numbers = F,cluster_cols =F,cluster_rows = F,fontface="bold",annotation_colors = ann_colors,annotation_col = annotation_col)

ggsave(p,filename = "G://embryo/figures/fig3/fig3E.png",width = 12, height = 8, units = "in", device='png')
ggsave(p,filename = "G://embryo/figures/fig3/fig3E.pdf",width = 12, height = 8, units = "in", device='pdf')


###Fig3 F

the_metatable = metatable[metatable$major == "macrophage",]
the_metatable$new = paste0(the_metatable$tissue,the_metatable$subtype)
for (new in unique(the_metatable$new)) {
  if (nrow(the_metatable[the_metatable$new == new,])<10){
    the_metatable= the_metatable[-which(the_metatable$new == new),]
  }
}

the_metatable = the_metatable[,c("time","subtype","Well_ID","tissue")]
the_metatable = the_metatable[-which(the_metatable$time == "Adult"),]
the_metatable = the_metatable[-which(the_metatable$time == "week19"),]
the_metatable$time <- factor(the_metatable$time,levels=time_order,ordered = TRUE)

tmp2 = the_metatable[,c("time","subtype")]

tmp2 = tmp2[!duplicated(tmp2),]
rownames(tmp2) = NULL
tmp2$rate = 0
for(t in rownames(tmp2)){
  tmp2[t,"rate"] = nrow(the_metatable[the_metatable$time == tmp2[t,"time"] & the_metatable$subtype == tmp2[t,"subtype"],])/nrow(the_metatable[the_metatable$time == tmp2[t,"time"],])
}


the_p1 = ggplot(tmp2,aes(x = time, 
                         stratum = subtype, 
                         alluvium = subtype,
                         y = rate,
                         fill = subtype,
                         label = subtype)) +
  scale_y_continuous(label = scales::percent_format(),
                     expand=c(0,0)) +
  scale_fill_manual(values = as.character(mac.color),breaks = names(mac.color)) +
  geom_flow() +
  geom_stratum(width = 0.6,colour = "grey30") +
  theme(axis.title=element_blank(),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.title=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(),
        axis.ticks=element_line(),
        axis.text = element_text(size = 12),
        legend.position = "none")
ggsave("E:\\mac_new/macrophage_ratio.pdf",the_p1,width =8,height =6,dpi = 300)

###Fig3 H
tmp.meta = metatable

tmp2 =  tmp.meta[tmp.meta$major == "macrophage",]
tmp2 = tmp2[tmp2$time != "Adult",]
tmp2 = tmp2[tmp2$time != "week19",]
the_dir = "E:\\mac_new/mac_bar2_re4/"
for(i in unique(tmp2$tissue)){
  print(i)
  tmp = tmp2[,c("tissue","time","subtype")]
  tmp = tmp[tmp$tissue == i,]
  tmp = tmp[,c("time","subtype")]
  tmp$new = paste(tmp$time,tmp$subtype,sep = "_")
  for (s in unique(tmp$time)){
    if(nrow(tmp[tmp$time == s,])<25){
      tmp = tmp[tmp$time != s,]
    }
  }
  tmp = tmp[,c("time","subtype")]
  other = setdiff(time_order,c(unique(tmp$time),"Adult","week19"))
  other = as.data.frame(other)
  other$blank = "blank"
  colnames(other) = colnames(tmp)
  tmp = rbind(tmp,other)
  tmp$time <- factor(tmp$time,levels=time_order,ordered = TRUE)
  tmp$subtype <- factor(tmp$subtype,levels=c("HdpM","YsdM_AFP_high","YsdM_AFP_low","pre_microglia",
                                             "red_pulp","Kupffer_cell","pre_PraM","PraM","gonad_macrophage",
                                             "Adrenalgland_macrophage","intestine CD209+ Mφ",
                                             "intestine CD207+ Mφ","langerhans","osteoclast","microglia","blank"),ordered = TRUE)
  p2 = ggplot(tmp,aes(x= time,fill = subtype))+
    geom_bar(stat="count",position = "fill",width = 0.8,color = "grey30")+
    xlab("")+ylab("")+#ggtitle(i)+
    #ggtitle(paste(i," total cell number ",nrow(tmp),sep = ""))+
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
  ggsave(paste(the_dir,i,"_tissue.png",sep = ""),p2,width = 6.8,height = 6)
}
