####Fig5 sup
#### writer： Feng Ruoqing
#### Date：2022-10-4

##fig5sA

#####
prep<-function(A,B,a_reduced,FC_cutoff,mf_new_mc){
  
  
  all_cellids_A = mf_new_mc[mf_new_mc$subtype %in% A,]
  data_A<-a_reduced[,colnames(a_reduced)%in%all_cellids_A$Well_ID]
  
  all_cellids_B = mf_new_mc[mf_new_mc$subtype %in% B,]
  data_B<-a_reduced[,colnames(a_reduced)%in%all_cellids_B$Well_ID]
  
  
  A_row_mean = data_calc_row_mean(data_A,A)
  B_row_mean = data_calc_row_mean(data_B,B)
  
  data_prep <-merge(A_row_mean,B_row_mean,by='row.names',all=TRUE)
  rownames(data_prep)<-data_prep$Row.names
  data_prep<-data_prep[,-1]
  
  
  data_prep<-log2(data_prep+1)
  
  log2FC = (data_prep[,1])-(data_prep[,2])
  
  data_prep_FC<-cbind(data_prep,log2FC)
  
  
  
  gene_tag<-tag(log2FC,as.numeric(FC_cutoff))
  
  data_prep2<-cbind(data_prep_FC,gene_tag)
  
  
  return(data_prep2)
  
}

scalar_plot<-function(data_prep2,A,B,A_1,B_1,genelist,A_color,B_color,path){
  
  
  
  data_prep2$label = ""
  data_prep2$gene = rownames(data_prep2)
  data_prep2$label[data_prep2$gene%in%genelist]=data_prep2$gene[data_prep2$gene%in%genelist]
  
  
  
  corr<-cor(data_prep2[,1],data_prep2[,2])
  corr<-round(corr,digits=4)
  
  colnames(data_prep2)[1:2] = c(A,B)
  p<-ggplot(data_prep2,aes(x = get(A),y = get(B),label = label))+
    geom_point(aes(color = gene_tag), size = 1)+xlim(0,7)+ylim(0,7)+
    geom_text_repel(data = data_prep2[data_prep2 != "",],size = 6,max.overlaps = 100)+
    theme_bw() + 
    
    labs(x = A_1, y = B_1, color = '',subtitle = paste("Correlation:",corr,sep="")) + #鍧愭爣杞存爣棰樿缃?
    
    geom_abline(intercept = 0.585, slope = 1, col = 'black', linetype = 'dashed', size = 0.5) + #杩?3鍙ョ敤浜庢坊鍔? |log2FC|>1 鐨勯槇鍊肩嚎
    
    geom_abline(intercept = -0.585, slope = 1, col = 'black', linetype = 'dashed', size = 0.5) +
    
    geom_abline(intercept = 0, slope = 1, col = 'black', linetype = 'dashed', size = 0.8)+
    
    scale_color_manual(values = c(A_color, 'gray', B_color)) +theme(legend.position = 'none')+theme(text = element_text(size = 20))  
  
  ggsave(p,filename = paste(path,A,"_VS_",B,"_scalar_plot.png",sep=""),width = 8, height = 8, units = "in", device='png')
  ggsave(p,filename = paste(path,A,"_VS_",B,"_scalar_plot.pdf",sep=""),width = 8, height = 8, units = "in", device='pdf')
  
}

data_calc_row_mean<-function(data,subtype){
  
  data<-as.matrix(data)
  
  data = apply(data,1,mean)
  
  data = as.data.frame(data)
  
  colnames(data)<-subtype
  
  return(data)
  
}

tag<-function(data,cutoff){
  
  color_tag<-vector()
  
  for(i in 1:length(data)){
    
    if(data[i]>=cutoff){
      
      color_tag[i]<-"up"
      
    }else if(data[i]< (-cutoff) ){
      
      color_tag[i]<-"down"
      
    }else {
      color_tag[i]<-"none"
    }
  }
  
  return(color_tag)
  
}

load("E:\\fig5/mormalized_umitab/a_reduced.rdata")
#load("E:\\fig5/mormalized_umitab/mf_new_mc.rdata")
A<-"PraM"
B<-"pre_PraM"

A_color<-'#2CA02CFF'
B_color<-'#FF7F0EFF'

FC_cutoff<-"0.585"
genelist<-unique(c("CD83","IL1B","SOD2","PTGS2","VEGFA","ICAM1","BTG1","C5AR1","DAB2",
            c("C5AR1","CXCL3","CXCL8","DAB2","PTGS2","SOD2","TNF","IL1B","HMGA2","CD36","CD83","IL1A",
              "VEGFA","GLUL","ICAM1")))

data_prep2<-prep(A,B,a_reduced,FC_cutoff,mf_new_mc)

A_1<-"PraM"
B_1<-"pre_PraM"

path<-"E:\\fig5/"

scalar_plot(data_prep2,A,B,A_1,B_1,genelist,A_color,B_color,path)



###Fig5s C

mac1.metatable = mac12.metatable[mac12.metatable$subtype == "PraM",]
mac1.metatable = mac1.metatable[mac1.metatable$time != "Adult",]
mac1.metatable = mac1.metatable[mac1.metatable$tissue %in% c("Brain","Skin","Kidney","Female gonad","Lung","Adrenalgland","Male gonad","Heart"),]

tmp = mac1.metatable
#tmp$tissue = tmp$subtype

tmp = tmp[,colnames(metatable)]

tmp2 = metatable[metatable$tissue == "Brain" & metatable$subtype == "microglia",]
tmp2
tmp2$tissue = "brain_microglial"
tmp = rbind(tmp,tmp2)

tmp2 = metatable[metatable$tissue == "Skin" & metatable$subtype == "microglia",]
tmp2
tmp2$tissue = "skin_microglial"
tmp = rbind(tmp,tmp2)

tmp2 = metatable[metatable$tissue == "Skin" & metatable$subtype == "langerhans",]
tmp2
tmp2$tissue = "langerhans"
tmp = rbind(tmp,tmp2)

tmp2 = metatable[metatable$tissue == "Male gonad" & metatable$subtype == "osteoclast",]
tmp2
tmp2$tissue = "osteoclast"
tmp = rbind(tmp,tmp2)

tmp2 = metatable[metatable$tissue == "Male gonad" & metatable$subtype == "gonad_macrophage",]
tmp2
tmp2$tissue = "gonad_macrophage"
tmp = rbind(tmp,tmp2)

tmp2 = metatable[metatable$tissue == "Adrenalgland" & metatable$subtype == "Adrenalgland_macrophage",]
tmp2
tmp2$tissue = "Adrenalgland_macrophage"
tmp = rbind(tmp,tmp2)

mac1.metatable = tmp

mac1.seurat = a[discard_gene2(rownames(a)),tmp$Well_ID]
mac1.seurat = CreateSeuratObject(counts = mac1.seurat)
mac1.seurat <- NormalizeData(mac1.seurat, normalization.method =  "LogNormalize")
mac1.seurat <- FindVariableFeatures(mac1.seurat,selection.method = "vst", nfeatures = 5000)
mac1.seurat <- ScaleData(mac1.seurat)
mac1.seurat <- RunPCA(object = mac1.seurat, pc.genes = VariableFeatures(mac1.seurat))
mac1.seurat <- FindNeighbors(mac1.seurat, dims = 1:30)
mac1.seurat <- FindClusters(mac1.seurat, resolution = 0.5)

tmp3 = tmp$tissue
names(tmp3) = tmp$Well_ID
Idents(mac1.seurat) = tmp3

i = "Brain"
final_marker = FindMarkers(mac1.seurat,ident.1 = i,ident.2 = c("Brain","Skin","Kidney","Female gonad","Lung","Adrenalgland","Male gonad","Heart"))
final_marker$gene = rownames(final_marker)
final_marker$sample = "discard"

for (i in unique(mac1.metatable$tissue)){
  print(i)
  tmp = FindMarkers(mac1.seurat,ident.1 = i,ident.2 = c("Brain","Skin","Kidney","Female gonad","Lung","Adrenalgland","Male gonad","Heart"))
  tmp$gene = rownames(tmp)
  tmp$sample = i
  final_marker = rbind(final_marker,tmp)
}

final_marker2 = final_marker
final_marker= final_marker2
final_marker = final_marker[final_marker$sample != "discard",]
final_marker$sig = "non"
final_marker$sig[final_marker$avg_log2FC >= 1 & final_marker$p_val_adj < 0.05] = "up"
final_marker$sig[final_marker$avg_log2FC <= -1 & final_marker$p_val_adj < 0.05] = "down"

rownames(final_marker) = NULL
final_marker$sample = factor(final_marker$sample,levels = c("Female gonad","Male gonad","Lung","Kidney","Heart",
                                                            "Adrenalgland","Skin","Brain",
                                                            "Adrenalgland_macrophage",
                                                            "gonad_macrophage","osteoclast",
                                                            "langerhans","skin_microglial",
                                                            "brain_microglial"))
p1 = ggplot(final_marker)+
  geom_jitter(aes(x = sample,y = avg_log2FC,color = sig))+
  theme_bw()+xlab('')+
  scale_color_manual(breaks = c("non","down","up"),values = c("grey66","blue","red"))+
  theme(legend.position = 'none',axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
ggsave("E:\\fig5/deg.pdf",p1,width = 8,height = 6,dpi = 300)

###Fig5s D
mac1.metatable = mac12.metatable[mac12.metatable$time != "Adult",]
mac1.metatable = mac1.metatable[,colnames(metatable)]
mac1.metatable = mac1.metatable[mac1.metatable$tissue %in% c("Kidney","Lung","Adrenalgland",
                                                             "Heart","Brain","Skin","Female gonad","Male gonad"),]
mac1.metatable = mac1.metatable[mac1.metatable$subtype == "PraM",]
skin = metatable[metatable$tissue == "Skin" & metatable$subtype %in% c("microglia","langerhans"),]
skin$tissue = skin$subtype
new.mac1 = rbind(mac1.metatable,skin)
brain = metatable[metatable$tissue == "Brain" & metatable$subtype == "microglia",]
brain$tissue = "Brain_microglia"
new.mac1 = rbind(new.mac1,brain)
gonad_male = metatable[metatable$tissue == "Male gonad" & metatable$subtype %in% c("gonad_macrophage","osteoclast"),]
gonad_male$tissue = gonad_male$subtype
new.mac1 = rbind(new.mac1,gonad_male)
Adrenalgland = metatable[metatable$tissue == "Adrenalgland" & metatable$subtype == "Adrenalgland_macrophage",]
Adrenalgland$tissue = Adrenalgland$subtype
new.mac1 = rbind(new.mac1,Adrenalgland)

new.mac1 = new.mac1[new.mac1$time %in% after9,]

new.mac1.seurat = a[discard_gene2(rownames(a)),new.mac1$Well_ID]
new.mac1.seurat = CreateSeuratObject(counts = new.mac1.seurat)
new.mac1.seurat <- NormalizeData(new.mac1.seurat, normalization.method =  "LogNormalize")
new.mac1.seurat <- FindVariableFeatures(new.mac1.seurat,selection.method = "vst", nfeatures = 2000)
new.mac1.seurat <- ScaleData(new.mac1.seurat)
new.mac1.seurat <- RunPCA(object = new.mac1.seurat, pc.genes = VariableFeatures(mac1.seurat))
new.mac1.seurat <- FindNeighbors(new.mac1.seurat, dims = 1:30)
new.mac1.seurat <- FindClusters(new.mac1.seurat, resolution = 0.5)

tmp = new.mac1$tissue
names(tmp) = new.mac1$Well_ID
Idents(new.mac1.seurat) = tmp
DimPlot(new.mac1.seurat)

mac1.pca = as.data.frame(new.mac1.seurat@reductions$pca@cell.embeddings)
mac1.pca = mac1.pca[,c(1,2)]
mac1.center = colMeans(mac1.pca[new.mac1$Well_ID[new.mac1$subtype == "PraM"],])

dis.embryo.subtype = matrix(NA,nrow = length(unique(new.mac1$embryo)),ncol = length(unique(new.mac1$tissue)))
colnames(dis.embryo.subtype) = unique(new.mac1$tissue)
rownames(dis.embryo.subtype) = unique(new.mac1$embryo)
for (i in unique(unique(new.mac1$embryo))){
  for (j in unique(new.mac1$tissue[new.mac1$embryo == i])){
    if(nrow(new.mac1[new.mac1$embryo == i & new.mac1$tissue == j,])>10){
      tmp = mac1.pca[new.mac1$Well_ID[new.mac1$embryo == i & new.mac1$tissue == j],]
      tmp["center",] = mac1.center
      tmp.dist = as.matrix(dist(tmp,method = "euclidean"))
      tmp = as.data.frame(tmp.dist["center",])
      colnames(tmp) = "value"
      tmp = tmp[setdiff(rownames(tmp),"center"),]
      q1 = quantile(tmp,0.001)
      q99 = quantile(tmp,0.999)
      tmp = tmp[tmp>q1]
      tmp = tmp[tmp<q99]
      dis.embryo.subtype[i,j] = mean(tmp)
    }
  }
}
dis.embryo.subtype = as.data.frame(dis.embryo.subtype)
dis.embryo.subtype$embryo = rownames(dis.embryo.subtype)
dis.embryo.subtype2 = reshape2::melt(dis.embryo.subtype)

dis.embryo.subtype2$subtype = "PraM"
dis.embryo.subtype2$subtype[dis.embryo.subtype2$variable %in% c("Adrenalgland_macrophage")] = "Adrenalgland_macrophage"
dis.embryo.subtype2$subtype[dis.embryo.subtype2$variable %in% c("gonad_macrophage")] = "gonad_macrophage"
dis.embryo.subtype2$subtype[dis.embryo.subtype2$variable %in% c("osteoclast")] = "osteoclast"
dis.embryo.subtype2$subtype[dis.embryo.subtype2$variable %in% c("langerhans")] = "langerhans"
dis.embryo.subtype2$subtype[dis.embryo.subtype2$variable %in% c("microglia","Brain_microglia")] = "microglia"
dis.embryo.subtype2$variable = factor(dis.embryo.subtype2$variable,levels = c("Lung","Kidney","Heart","Adrenalgland",
                                                                              "Female gonad","Male gonad","Skin","Brain",
                                                                              "Adrenalgland_macrophage","gonad_macrophage",
                                                                              "osteoclast","langerhans","microglia",
                                                                              "Brain_microglia"))
p1 = ggplot(dis.embryo.subtype2,aes(x = variable, y = value))+
  geom_boxplot(aes(fill = subtype),color = "grey30",outlier.shape = NA)+
  geom_jitter(shape = 21,color = "grey30",fill = NA,size = 1.5)+
  geom_signif(comparisons = list(c("Adrenalgland","Adrenalgland_macrophage"),
                                 c("Male gonad","osteoclast"),
                                 c("Male gonad","gonad_macrophage"),
                                 c("Skin","microglia"),
                                 c("Skin","langerhans"),
                                 c("Brain","Brain_microglia")), 
              map_signif_level = T,step_increase =0.1,test = wilcox.test)+
  theme_bw()+theme_classic()+
  scale_fill_manual(values = as.character(mac.color),breaks = names(mac.color))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
ggsave("E:\\fig5/diversity_score.pdf",p1,width = 5,height = 6,dpi = 300)
