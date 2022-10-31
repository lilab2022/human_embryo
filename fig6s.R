###Figure 6S
###writer: Dr. Wang Hao
###Date: 2022-10-4

###Figure6S A
average_running<-function(test){
  
  test2<-matrix(nrow=nrow(test),ncol=ncol(test))
  
  j=1
  for(i in 1:nrow(test)){
    
    
    
    test2[j,] = roll_mean(test[i,], n = 5, align = "right", fill = NA)
    
    test2[j,] = format(test2[j,], digits = 3)
    
    # print(test2[,j])
    
    j=j+1
  }
  
  rownames(test2)<-rownames(test)
  colnames(test2)<-colnames(test)
  
  test2<-as.data.frame(test2)
  
  test2[,1:4]<-test[,1:4]
  
  
  test2=apply(test2,2,as.numeric)
  rownames(test2)<-rownames(test)
  
  return(test2)
  
  
}

discard_gene2 = function(the_str){
  gene_all <-the_str
  a_relate <- str_subset(gene_all, pattern ='^AC[0-2]|^AE[0-9]|^AF[0-9]|^AJ[0-9]|^AL[0-9]|^AP[0-9]|^FP[0-9]')
  ex9 <- str_subset(gene_all, pattern = '-OT$|-OT[0-9]$|-OT[0-9]\\.[0-9]$|^CDK-AS$|-AS[0-9]$|-AS[0-9]\\.[0-9]$|^FTH|^FTL|^HSP')
  
  ex2 <- str_subset(gene_all, pattern ='^BX[0-9]|^CH17-|^CITF|^CR[0-9]|^CT009|^CT476828|^CT86797|^CT9786|^TUB|^FOS|^CCL|^DNAJA|^DNAJB|^HB')
  ex3 <- str_subset(gene_all, pattern ='^CTA-|^CTB-|^CTC-|^CTD-|^CU[0-9]|^D86|^D87|^DAQB-|^DASS-|^DKFZP|^FAM|^GADD|^HSP|^CDC|^MYH|^XIST')
  ex4 <- str_subset(gene_all, pattern ='^GHc-|^GS1-|^HIST|^hsa-mir|^IGDCC|^IGSF|^MYL|^NFKB|^COL|^MCM|^CXCL|^HB|^Y_|^Metazoa_|^CX3CR|^CCL')
  ex5 <- str_subset(gene_all, pattern ='^KB-|^KIF|^L29074|^L34079|^LA16c-|^LINC|^LL[0-9][0-9]|^LLNL|^MIR|^MRP|^SLC|^TNF|^MKI|^CXCR')
  ex6 <- str_subset(gene_all, pattern ='^MT|^NCRNA|^NDUFA|^NDUF|^NPM|^RN7|^RNA5|^RNU[0-9]|^RNVU|^RNY|^RP[0-9]|^RPL|^RPS|^PDE|^EGR|^IFI|^TGF|^USP|^PPP|^ATP|^DUSP|^GBP|^DNAJ|^ATF|^MALAT|^TOP|^DDIT|^NFKB')
  ex7 <- str_subset(gene_all, pattern ='^SCARNA|^sno|^SNOR|^SNRP|^TRAV|^TRBV|^U[0-9]|^WI2-|^XX-|^XXbac-|^XXcos-|^XXyac-|^Y-RNA|^Z[0-9]')
  ex8 <- str_subset(gene_all, pattern = '-IT$|-IT[0-9]$|-IT[0-9]\\.[0-9]$|-SLIT[0-9]$|-SLIT[0-9]\\.[0-9]$')
  
  exclude_all <- c(a_relate,ex2,ex3,ex4,ex5,ex6,ex7,ex8,ex9)
  gene_all = setdiff(gene_all,exclude_all)
  return(gene_all)
} 

a<-read.csv("G://embryo/20220329/cleanumi.csv",header=T,row.names = "X")
load("G://embryo/20220929/metatable_0929.Rdata")

adult_mc<-c("X2056","X2350","X282","X871","X891")
a<-a[,!colnames(a)%in%adult_mc]



excluded_genes<-discard_gene2(rownames(a))

a_reduced<-a[rownames(a)%in%excluded_genes,]

############global normalization

seob_data <- CreateSeuratObject(counts = a_reduced)
seob_data <- NormalizeData(seob_data, normalization.method =  "LogNormalize") 

a_reduced <- seob_data@assays$RNA@data


### extract umi-tab based on metacell IDs obtained in the last step.

mcmc1<-metatable

adult.mc<-c("X2056","X2350","X282","X871","X891")

mcmc1 <- mcmc1[!mcmc1$mc%in%adult.mc,]
all_mac = c(unique(mcmc1$mc[mcmc1$subtype %in% c("pre_PraM","PraM","YsdM_AFP_high","YsdM_AFP_low")]))

################################

data<-a_reduced

data_sub<-data[,colnames(data)%in%all_mac]


row_var_genes<-rowVars(as.matrix(data_sub[,1:ncol(data_sub)]))
data_sub<-as.data.frame(data_sub)
data_sub$row_var <- row_var_genes


data_sub_ordered<-data_sub[order(data_sub$row_var,decreasing = T),]
data_sub_ordered<-data_sub_ordered[1:2000,] # using the 2000 genes having the largest variance
data_sub_ordered<-data_sub_ordered[,-ncol(data_sub_ordered)] # removing the column of "row_var"



Time<-read.csv("G://embryo/20220930/subtypes_mf_20220930.csv")
rownames(Time)<-Time$MC

realtime<-Time[Time$MC%in%all_mac,]


data.t<-t(data_sub_ordered)
data.t<-cbind(rownames(data.t),data.t)
colnames(data.t)[1]<-"MC"
data.t<-as.data.frame(data.t)
data.t<-merge(realtime,data.t,by="MC")
data.t<-data.t[order(data.t$the_median,decreasing = F),]

meta.info<-data.t[,c(1,2)]

data.t.t<-t(data.t)
colnames(data.t.t)<-data.t$MC
data.t.t<-data.t.t[-c(1,2),]
data.t.t<-as.data.frame(data.t.t)

anno.2 = mcmc1[,c("mc","subtype","major")]
anno.2 = anno.2[!duplicated(anno.2),]
anno.2 = anno.2[anno.2$mc %in% all_mac,]
rownames(anno.2) = anno.2$mc
anno.2 = anno.2[,-1]


write.csv(data.t.t,"G://embryo/20220930/PraM/data_time_order_PraM.csv") 
write.csv(anno.2,"G://embryo/20220930/PraM/anno2_PraM.csv")

## transfer to monocle3 analysis to obtain pseudo time

data = read.csv("G://embryo/20220930/PraM/data_time_order_PraM.csv")
data = data[-c(1:5),]
rownames(data) = data[,1]

data = data[,-1]
data = as.matrix(data)

anno = read.csv("G://embryo/20220930/PraM/anno2_PraM.csv",row.names = 1)
anno = anno[colnames(data),]

gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

data.time.pseudo<-read.csv("G://embryo/20220930/slingshot_bifurcation_pseudo_20220930.csv",header=T) # obtain pseudotime produced by slingshot

colnames(data.time.pseudo)[1]<-"MC"
data.time.pseudo<-merge(data.time.pseudo,data.t,by="MC")
colnames(data.time.pseudo)[2]<-"pseudo.time"
data.time.pseudo<-data.time.pseudo[order(data.time.pseudo$pseudo.time,decreasing = F),]
meta.info2<-data.time.pseudo[,c(1:3)]


data.pseudo.t<-t(data.time.pseudo)
colnames(data.pseudo.t)<-data.pseudo.t[1,]
data.pseudo.t<-data.pseudo.t[-c(1,2,3),]
data.pseudo.t<-as.data.frame(data.pseudo.t)

########

rownames(meta.info2)<-meta.info2$MC

anno.2<-cbind(anno.2,rownames(anno.2))
colnames(anno.2)[3]<-"MC"

meta.info2<-merge(meta.info2,anno.2,by="MC")
meta.info2<-meta.info2[order(meta.info2$pseudo.time,decreasing = F),]

rownames(meta.info2)<-meta.info2$MC
meta.info2<-meta.info2[,c(4,3)]
colnames(meta.info2)<-c("Subtype","Time")


## obtain module score

Time<-read.csv("G://embryo/20220930/")
rownames(Time)<-Time$MC

realtime<-Time[Time$MC%in%all_mac,]

realtime<-realtime[realtime$MC%in%rownames(meta.info2),]
realtime<-realtime[match(rownames(meta.info2),realtime$MC),]

meta.info2$Time<-realtime$the_median

proliferation<-read.table("G://embryo/20220804/G2M_mc.csv")
proliferation<-cbind(proliferation,rownames(proliferation))
proliferation<-as.data.frame(proliferation)  


differentiation<-read.csv("G://embryo/20220808/cytotrace.value.bymc.csv",row.names = "MC")



#proliferation<-proliferation[rownames(proliferation)%in%rownames(meta.info2),]
differentiation<-differentiation[rownames(differentiation)%in%rownames(meta.info2),]

#prolifearation<-as.data.frame(proliferation)
differentiation<-as.data.frame(differentiation)

#proliferation<-proliferation[match(rownames(meta.info2),rownames(proliferation)),]
differentiation<-differentiation[match(rownames(meta.info2),rownames(differentiation)),]


#gene_modules<-read.table("G://embryo/20220808/fp_module(1).csv",header=T)
gene_modules<-read.table("G://embryo/20220930/fp_module(2).csv",header=T)
gene_modules<-t(gene_modules)
gene_modules<-as.data.frame(gene_modules)
gene_modules<-gene_modules[rownames(gene_modules)%in%rownames(meta.info2),]
gene_modules<-gene_modules[match(rownames(meta.info2),rownames(gene_modules)),]


monocle3_pseudo<-read.csv("G://embryo/20220930/monocle3_pse2_PraM.csv",header=T)
monocle3_pseudo<-monocle3_pseudo[monocle3_pseudo$X%in%rownames(meta.info2),]
monocle3_pseudo<-monocle3_pseudo[match(rownames(meta.info2),monocle3_pseudo$X),]

meta.info2$G2M_pct<-proliferation$G2M_pct
meta.info2$cytotrace<-differentiation$cytotrace
meta.info2$pseudoTime<-data.time.pseudo$pseudo.time
meta.info2$gene.module<-gene_modules$PraM
#meta.info2$slingshot<-slingshot_pseudo$slingshot_pseudo
meta.info2$Monocle3<-monocle3_pseudo$pseudotime.cds.

########################
## color schemes
colors_schemes<-colorRampPalette(c("white", "green"))(length(unique(meta.info2$Time))) 

tmp1<-(unique(meta.info2$Time))
tmp1<-as.data.frame(tmp1)
tmp1<-tmp1[mixedorder(tmp1$tmp1),]
tmp1<-as.vector(tmp1)

names(colors_schemes)<-tmp1

#colors_schemes2<-colorRampPalette(c("white", "blue"))(length(unique(data.time.pseudo$pseudo.time))) 
colors_schemes2<-colorRampPalette(c("white", "blue"))(length(unique(meta.info2$pseudoTime))) 

tmp2<-(unique(meta.info2$pseudoTime))
tmp2<-as.data.frame(tmp2)
tmp2<-tmp2[mixedorder(tmp2$tmp2),]
tmp2<-as.vector(tmp2)

names(colors_schemes2)<-tmp2

colors_schemes22<-colorRampPalette(c("white", "#17BECFFF"))(length(unique(meta.info2$Monocle3))) 

tmp22<-(unique(meta.info2$slingshot))
tmp22<-as.data.frame(tmp22)
tmp22<-tmp22[mixedorder(tmp22$tmp22),]
tmp22<-as.vector(tmp22)

names(colors_schemes22)<-tmp22



colors_schemes3<-colorRampPalette(c("white", "darkgreen"))(length(unique(meta.info2$cytotrace))) 

tmp3<-(unique(meta.info2$cytotrace))
tmp3<-as.data.frame(tmp3)
tmp3<-tmp3[mixedorder(tmp3$tmp3),]
tmp3<-as.vector(tmp3)

names(colors_schemes3)<-tmp3

## proliferation

colors_schemes4<-colorRampPalette(c("white", "red"))(length(unique(meta.info2$G2M_pct))) 

tmp4<-(unique(meta.info2$G2M_pct))
tmp4<-as.data.frame(tmp4)
tmp4<-tmp4[mixedorder(tmp4$tmp4),]
tmp4<-as.vector(tmp4)

names(colors_schemes4)<-tmp4

## gene.module

colors_schemes5<-colorRampPalette(c("white", "purple"))(length(unique(meta.info2$gene.module))) 

tmp5<-(unique(meta.info2$gene.module))
tmp5<-as.data.frame(tmp5)
tmp5<-tmp5[mixedorder(tmp5$tmp5),]
tmp5<-as.vector(tmp5)

names(colors_schemes5)<-tmp5

bk <- c(seq(-3,-0.1,by=0.01),seq(0,3,by=0.01))

groups = tribble(~group_id,~color,"PraM","#FF7F0EFF","YsdM_AFP_high","#9467BDFF","YsdM_AFP_low","#E377C2FF","pre_PraM","#2CA02CFF")

value<-groups$color
names(value)<-groups$group_id

ann_colors = list(
  Time = colors_schemes,
  Slingshot=colors_schemes2,
  Monocle3=colors_schemes22,
  cytotrace = colors_schemes3,
  G2M_pct = colors_schemes4,
  gene.module = colors_schemes5,
  Subtype=value
)



#meta.info2<-cbind(meta.info2,data.time.pseudo$pseudo.time)
meta.info2<-as.data.frame(meta.info2)
#colnames(meta.info2)[3]<-"PseudoTime"
#colnames(meta.info2)[4]<-"scores"

annotation_col = meta.info2

##
data.pseudo.t<-data.pseudo.t[-c(1:5),]
##

test=apply(data.pseudo.t,2,as.numeric)
rownames(test)<-rownames(data.pseudo.t)


####

test2<-average_running(test)

colnames(annotation_col)[5]<-"Slingshot"

p3<-pheatmap(test2, cluster_rows = T, cluster_cols = F,scale = "row",color = c(colorRampPalette(colors = c("#1B519C","white"))(length(bk)/2),colorRampPalette(colors = c("white","#A00021"))(length(bk)/2)),
             show_rownames = F,show_colnames = F,cutree_rows = 30,breaks = bk,annotation_col = annotation_col, annotation_colors = ann_colors,treeheight_row=0)


ggsave("G://embryo/20220930/PraM/fig6_PraM_real_split2.png",p3,width = 10,height = 8,dpi = 300,units="in",device="png")
ggsave("G://embryo/20220930/PraM/fig6_PraM_real_split2.pdf",p3,width = 10,height = 8,dpi = 300,units="in",device="pdf")


data.time.ann<-test2

row_cluster <- cutree(p3$tree_row, k=30)
newOrder <- data.time.ann[p3$tree_row$order,]
newOrder<-cbind(newOrder,row_cluster[match(rownames(newOrder), names(row_cluster))])
colnames(newOrder)[ncol(newOrder)]="Cluster"
write.csv(newOrder, "G://embryo/20220930/PraM/cluster30_anno2_pseudo_PraM.csv")


### re-order

newOrder<-as.data.frame(newOrder)
save(newOrder,file="G://embryo/20220930/PraM/newOrder_PraM.rdata")
genes<-data.frame(Genes=rownames(newOrder),Cluster=newOrder$Cluster)

new_order<-c("11","14","12","10","30","28","17","7","5","6","4","2","15")

template<-matrix(,nrow=1,ncol=2)
colnames(template)<-c("Genes","Cluster")
template<-as.data.frame(template)

for(i in 1:length(unique(genes$Cluster))){
  
  
  new_i<-new_order[i]
  
  new_row<-subset(genes,Cluster==new_i)
  
  template<-rbind(template,new_row)
  
}

template<-as.data.frame(template)

template<-template[-1,]

write.csv(template,"G://embryo/20220930/PraM/gene_reorder_PraM.csv")


test3<-test2[rownames(test2)%in%template$Genes,]
test3<-test3[match(template$Genes,rownames(test3)),]

p4<-pheatmap(test3, cluster_rows = F, cluster_cols = F,scale = "row",color = c(colorRampPalette(colors = c("#1B519C","white"))(length(bk)/2),colorRampPalette(colors = c("white","#A00021"))(length(bk)/2)),
             show_rownames = F,show_colnames = F,cutree_rows = 1,breaks = bk,annotation_col = annotation_col, annotation_colors = ann_colors,treeheight_row=0)

ggsave("G://embryo/20220907/fig6_PraM_real_complete.png",p4,width = 10,height = 8,dpi = 300,units="in",device="png")

col_fun = colorRamp2(c(-3, 0, 3), c("#1B519C", "white", "#A00021"))
col_fun(seq(-3, 3))

genes<-read.csv("G://embryo/20220907/gene_list_PraM.csv",header=F)
genes <- as.data.frame(genes)
colnames(genes)<-c("genes")


col_time = colorRamp2(c(min(meta.info2$Time),max(meta.info2$Time)), c("white", "green"))
#col_G2M = colorRamp2(c(min(meta.info2$G2M_pct),max(meta.info2$G2M_pct)), c("white", "red"))
col_pseudotime = colorRamp2(c(min(meta.info2$Monocle3),max(meta.info2$Monocle3)), c("white", "blue"))
col_cytotrace <- colorRamp2(c(min(meta.info2$cytotrace),max(meta.info2$cytotrace)), c("white", "darkgreen"))
col_score<- colorRamp2(c(min(meta.info2$gene.module),max(meta.info2$gene.module)), c("white", "purple"))
col_subtype <- value
col_slingshot<-colorRamp2(c(min(meta.info2$pseudoTime),max(meta.info2$pseudoTime)), c("white", "#17BECFFF"))

column_ann <- HeatmapAnnotation(
  Slingshot = meta.info2$pseudoTime,
  Monocle3 = meta.info2$Monocle3,
  cytotrace = meta.info2$cytotrace,
  Score = meta.info2$gene.module,
  Time = meta.info2$Time,
  Subtype = meta.info2$Subtype,
  
  show_legend = F,
  #col = list(Time = col_time, G2M = col_G2M, Pseudotime = col_pseudotime, cytotrace = col_cytotrace, Subtype = col_subtype)
  col = list(  Slingshot=col_slingshot,Monocle3 = col_pseudotime,cytotrace = col_cytotrace,Score = col_score,Time = col_time, Subtype = col_subtype)
)
gene.local<-c()


for(i in 1:nrow(genes)){
  
  tmp.local<-which(rownames(test3) %in% genes$genes[i])
  gene.local<-c(gene.local,tmp.local)
}

gene.local<-cbind(gene.local,genes)


gene.local2<-gene.local[order(gene.local$gene.local),]

B <- Heatmap( t(scale(t(test3))),
              cluster_rows = F,
              cluster_columns = F,
              col = col_fun,
              show_row_names = F,
              show_column_names = F,
              show_heatmap_legend = F,
              top_annotation = column_ann
)+ rowAnnotation(link = anno_mark(at = gene.local2$gene.local, 
                                  labels = gene.local2$genes, labels_gp = gpar(fontsize = 10)))

#pdf(file="G://embryo/20220903/fig6_PraM_complexheatmap_v4¡ª2.pdf",width = 15,height = 12)
pdf(file="G://embryo/figures/figS6A.pdf",width = 15,height = 12)
dev.off()

###Figure6S B
TFs_all<-read.csv("G://embryo/20220527/TFs_all.csv",header=T)
###

a<-read.csv("G://embryo/20220329/cleanumi.csv",header=T,row.names = "X")
load("G://embryo/20220929/metatable_0929.Rdata")

mf_id<-read.csv("G://embryo/20220930/subtypes_mf_20220930.csv") # loading new version of macrophage meta information
adult_mc<-c("X2056","X2350","X282","X871","X891")
mf_id<-mf_id[!mf_id$MC%in%adult_mc,]


excluded_genes<-discard_gene2(rownames(a))

a_reduced<-a[rownames(a)%in%excluded_genes,]


############global normalization

seob_data <- CreateSeuratObject(counts = a_reduced)
seob_data <- NormalizeData(seob_data, normalization.method =  "LogNormalize") # without log transformation

a_reduced <- seob_data@assays$RNA@data


### extract umi-tab based on metacell IDs obtained in the last step.

all_mac = c(unique(mf_id$MC[mf_id$subtype %in% c("intestine CD209+ M¦Õ")]))
#all_mac = c(unique(mf_id$MC[mf_id$subtype %in% c("osteoclast")]))
#all_mac = c(unique(mf_id$MC[mf_id$subtype %in% c("PraM")]))
#all_mac = c(unique(mf_id$MC[mf_id$subtype %in% c("microglia")]))
#all_mac = c(unique(mf_id$MC[mf_id$subtype %in% c("langerhans")]))
#all_mac = c(unique(mf_id$MC[mf_id$subtype %in% c("Kupffer_cell")]))


data_sub<-a_reduced[,colnames(a_reduced)%in%all_mac]
row_mean_genes<-rowMeans2(as.matrix(data_sub[,1:ncol(data_sub)]))
data_sub<-as.data.frame(data_sub)
data_sub$row_var <- row_mean_genes

data_sub_ordered<-data_sub[order(data_sub$row_var,decreasing = T),]
data_sub_ordered<-data_sub_ordered[1:2000,]
data_sub_ordered<-data_sub_ordered[rownames(data_sub_ordered)%in%TFs_all$Genes,] #only TFs


TFs<-data_sub_ordered[,-ncol(data_sub_ordered)]
TFs<-t(TFs)


module<-read.table("G://embryo/20220930/fp_module(2).csv",header=T)
module<-t(module)

module<-module[rownames(module)%in%rownames(TFs),]
module<-module[match(rownames(TFs),rownames(module)),]
module<-as.data.frame(module)

module.score<-module$gut
module.score<-as.numeric(module.score)

value<-c()

for(i in 1:ncol(TFs)){
  
  tmp<-cor(as.numeric(TFs[,i]),module.score)
  
  value = c(value,tmp)
  
}

names(value)<-colnames(TFs)

value<-as.data.frame(value)
value<-na.omit(value)


write.csv(value,"G://embryo/20221002/TFs_corr_gut_v2.csv")

############################

TFs.data<-cbind(TFs,module.score)

write.csv(TFs.data,"G://embryo/20221002/TFs.data_gut_ML_v2.csv")

##################

data<-read.csv("G://embryo/20221002/TFs_corr_gut_v2.csv",header=T)
colnames(data)<-c("TFs","correlation")

data<-data[order(data$correlation,decreasing = T),]

data<-data[1:30,]

data$TFs=factor(data$TFs,levels = rev(data$TFs),ordered = T)

p<-ggplot(data, aes(x=TFs, y=correlation)) + geom_bar(stat="identity", fill="#AEC7E8FF")+coord_flip()+
  theme_bw()+theme(panel.grid=element_blank())+theme(plot.title = element_text(size = 60, face = "bold"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("G://embryo/20221002/gut_TFs2.png",p,width = 8,height = 10,units = "in",dpi = 150,device = "png")
ggsave("G://embryo/20221002/gut_TFs2.pdf",p,width = 8,height = 10,units = "in",device = "pdf")

###Figure6S C

# Kupffer
data<-read.csv("G://embryo/20221002/TFs_corr_kupffer_v2.csv",header=T)
colnames(data)<-c("TFs","correlation")

data<-data[order(data$correlation,decreasing = T),]

data<-data[1:30,]

data$TFs=factor(data$TFs,levels = rev(data$TFs),ordered = T)

p<-ggplot(data, aes(x=TFs, y=correlation)) + geom_bar(stat="identity", fill="#FF9896FF")+coord_flip()+
  theme_bw()+theme(panel.grid=element_blank())+theme(plot.title = element_text(size = 60, face = "bold"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("G://embryo/20221002/kupffer_TFs2.png",p,width = 8,height = 10,units = "in",dpi = 150,device = "png")
ggsave("G://embryo/20221002/kupffer_TFs2.pdf",p,width = 8,height = 10,units = "in",device = "pdf")


# Langerhans
data<-read.csv("G://embryo/20221002/TFs_corr_langerhans_v2.csv",header=T)
colnames(data)<-c("TFs","correlation")

data<-data[order(data$correlation,decreasing = T),]

data<-data[1:30,]

data$TFs=factor(data$TFs,levels = rev(data$TFs),ordered = T)

p<-ggplot(data, aes(x=TFs, y=correlation)) + geom_bar(stat="identity", fill="#DBDB8DFF")+coord_flip()+
  theme_bw()+theme(panel.grid=element_blank())+theme(plot.title = element_text(size = 60, face = "bold"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("G://embryo/20221002/Langerhans_TFs2.png",p,width = 8,height = 10,units = "in",dpi = 150,device = "png")
ggsave("G://embryo/20221002/Langerhans_TFs2.pdf",p,width = 8,height = 10,units = "in",device = "pdf")


# Ostoclast
data<-read.csv("G://embryo/20220904/TFs_corr_osteoclast_v2.csv",header=T)
colnames(data)<-c("TFs","correlation")

data<-data[order(data$correlation,decreasing = T),]

data<-data[1:30,]

data$TFs=factor(data$TFs,levels = rev(data$TFs),ordered = T)

p<-ggplot(data, aes(x=TFs, y=correlation)) + geom_bar(stat="identity", fill="#C49C94FF")+coord_flip()+
  theme_bw()+theme(panel.grid=element_blank())+theme(plot.title = element_text(size = 60, face = "bold"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("G://embryo/20221002/Ostoclast_TFs2.png",p,width = 8,height = 10,units = "in",dpi = 150,device = "png")
ggsave("G://embryo/20221002/Ostoclast_TFs2.pdf",p,width = 8,height = 10,units = "in",device = "pdf")

# Microglia
data<-read.csv("G://embryo/20221002/TFs_corr_microglia_v2.csv",header=T)
colnames(data)<-c("TFs","correlation")

data<-data[order(data$correlation,decreasing = T),]

data<-data[1:30,]

data$TFs=factor(data$TFs,levels = rev(data$TFs),ordered = T)

p<-ggplot(data, aes(x=TFs, y=correlation)) + geom_bar(stat="identity", fill="#1F77B4FF")+coord_flip()+
  theme_bw()+theme(panel.grid=element_blank())+theme(plot.title = element_text(size = 60, face = "bold"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("G://embryo/20221002/Microglia_TFs2.png",p,width = 8,height = 10,units = "in",dpi = 150,device = "png")
ggsave("G://embryo/20221002/Microglia_TFs2.pdf",p,width = 8,height = 10,units = "in",device = "pdf")


# PraM
data<-read.csv("G://embryo/20221002/TFs_corr_PraM_v2.csv",header=T)
colnames(data)<-c("TFs","correlation")

data<-data[order(data$correlation,decreasing = T),]

data<-data[1:30,]

data$TFs=factor(data$TFs,levels = rev(data$TFs),ordered = T)

p<-ggplot(data, aes(x=TFs, y=correlation)) + geom_bar(stat="identity", fill="#FF7F0EFF")+coord_flip()+
  theme_bw()+theme(panel.grid=element_blank())+theme(plot.title = element_text(size = 60, face = "bold"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("G://embryo/20221002/PraM_TFs2.png",p,width = 8,height = 10,units = "in",dpi = 150,device = "png")
ggsave("G://embryo/20221002/PraM_TFs2.pdf",p,width = 8,height = 10,units = "in",device = "pdf")

##Figure6S C
dat<-read.csv("G://embryo/20220904/TFs.data_microglia_ML_v2.csv",header=T,,row.names = "X")
TFs_list<-read.csv("G://embryo/20220904/TFs_corr_microglia_v2.csv")

colnames(TFs_list)<-c("TFs","correlation")
TFs_list<-TFs_list[order(TFs_list$correlation,decreasing = T),]
TFs_list<-TFs_list[c(1:30),]

datx<-dat[,colnames(dat)%in%TFs_list$TFs,]
datx<-as.data.frame(datx)



datx$module.score<-dat$module.score
datx<-as.data.frame(datx)
#colnames(datx)<-"module.score"

#glimpse(dat)
glimpse(datx)

set.seed(1001) 

index = sample(1:nrow(datx), nrow(datx)) 

k<-3 # k-cross validation
index_assigned<-cross_validaton(k,index)



PreProcess_Predict(index_assigned,datx,30,k)

PreProcess_Predict_regularization(index_assigned,datx,30,k)


PreProcess_Predict<-function(index_assigned,dat,var_num,num){
  
  num<-length(index_assigned)
  
  coefficient_all<-matrix(nrow=(var_num+1),ncol=num)
  
  coefficient_all<-as.data.frame(coefficient_all)
  
  predict_value<-c()
  
  for(i in 1:num){
    
    train = dat[-index_assigned[[i]],]  
    test = dat[index_assigned[[i]],]
    
    cols = colnames(dat)
    cols = cols[1:var_num]
    
    pre_proc_val <- preProcess(train[,cols], method = c("center", "scale"))
    pre_proc_val_test <- preProcess(test[,cols], method = c("center", "scale"))
    
    train[,cols] = predict(pre_proc_val, train[,cols])
    test[,cols] = predict(pre_proc_val_test, test[,cols])
    
    summary(train)
    
    lr = lm(module.score ~ ., data = train)
    summary(lr)
    
    col_names<-paste("round",i)
    print(col_names)
    
    predictions = predict(lr, newdata = train)
    result1<-eval_metrics(lr, train, predictions, target = 'module.score')
    
    
    predictions_test = predict(lr, newdata = test)
    result2<-eval_metrics(lr, test, predictions_test, target = 'module.score')
    
    predict_value<-c(predict_value,predictions_test)
    
    coefficient<-lr$coefficients
    
    coefficient_all[,i]<-coefficient
   
    ##    write.csv(as.data.frame(coefficient),paste("embryo/20220402/monocle_3/TFs_coefficient_",col_names,".csv"))
    
  }
  
  rownames(coefficient_all)<-names(lr$coefficients)
  
  mean_vector<-rowMeans(coefficient_all)
  names(mean_vector)<-names(lr$coefficients)
  
  
  write.csv(coefficient_all,"G://embryo/20220904/ML/TFs_coefficient_all.csv")
  write.csv(mean_vector,"G://embryo/20220904/ML/TFs_coefficent_mean.csv")
  write.csv(as.data.frame(predict_value),"G://embryo/20220904/ML/predictedValues.csv")
}


eval_metrics = function(model, df, predictions, target){
  resids = df[,target] - predictions
  resids2 = resids**2
  N = length(predictions)
  r2 = as.character(round(summary(model)$r.squared, 2))
  adj_r2 = as.character(round(summary(model)$adj.r.squared, 2))
  print(adj_r2) #Adjusted R-squared
  print(as.character(round(sqrt(sum(resids2)/N), 2))) #RMSE
  
  result<-c(adj_r2,as.character(round(sqrt(sum(resids2)/N), 2)))
  
  return(result)
}

##### adding regularizaiton


PreProcess_Predict_regularization<-function(index_assigned,dat,var_num,num){
  
  num<-length(index_assigned)
  
  coefficient_all<-matrix(nrow=(var_num),ncol=num)
  
  coefficient_all<-as.data.frame(coefficient_all)
  
  coefficient_ridge<-coefficient_all
  coefficient_lasso<-coefficient_all
  
  
  
  predict_value_ridge<-c()
  predict_value_lasso<-c()
  
  for(i in 1:num){
    
    train = dat[-index_assigned[[i]],]  
    test = dat[index_assigned[[i]],]
    
    
    cols_reg<-colnames(dat)
    
    
    dummies <- dummyVars(module.score ~ ., data = dat[,cols_reg]) 
    
    train_dummies = predict(dummies, newdata = train[,cols_reg])
    
    test_dummies = predict(dummies, newdata = test[,cols_reg])
    
    print(dim(train_dummies)); print(dim(test_dummies))
    
    x_train = as.matrix(train_dummies)
    y_train = train$module.score
    
    x_test = as.matrix(test_dummies)
    y_test = test$module.score
    
    lambdas <- 10^seq(2, -3, by = -.1)
    ridge_reg = glmnet(x_train, y_train, nlambda = 25, alpha = 0, family = 'gaussian', lambda = lambdas)
    
    summary(ridge_reg)
    
    ### find the best lambda parameter
    cv_ridge <- cv.glmnet(x_train, y_train, alpha = 0, lambda = lambdas)
    optimal_lambda <- cv_ridge$lambda.min
    
    print("Ridge regularization results")
    
    # Prediction and evaluation on train data
    predictions_train_ridge <- predict(ridge_reg, s = optimal_lambda, newx = x_train)
    print(eval_results(y_train, predictions_train_ridge, train))
    
    # Prediction and evaluation on test data
    predictions_test_ridge <- predict(ridge_reg, s = optimal_lambda, newx = x_test)
    print(eval_results(y_test, predictions_test_ridge, test))
    
    ridge_reg1 = glmnet(x_train, y_train, nlambda = 25, alpha = 0, family = 'gaussian', lambda = optimal_lambda)
    
    ################# lasso regularization
    print("Lasso regularization results")
    #lambdas <- 10^seq(2, -3, by = -.1)
    
    # Setting alpha = 1 implements lasso regression
    lasso_reg <- cv.glmnet(x_train, y_train, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
    
    # Best 
    lambda_best <- lasso_reg$lambda.min 
    
    lasso_model <- glmnet(x_train, y_train, alpha = 1, lambda = lambda_best, standardize = TRUE)
    
    predictions_train_lasso <- predict(lasso_model, s = lambda_best, newx = x_train)
    print(eval_results(y_train, predictions_train_lasso, train))
    
    predictions_test_lasso <- predict(lasso_model, s = lambda_best, newx = x_test)
    print(eval_results(y_test, predictions_test_lasso, test))
    
    predict_value_ridge<-c(predict_value_ridge,predictions_test_ridge)
    predict_value_lasso<-c(predict_value_lasso,predictions_test_lasso)
    
    tmp_x<-as.data.frame(ridge_reg1$beta)
    tmp_y<-as.data.frame(lasso_model$beta)
    
    coefficient_ridge[,i]<-tmp_x$s0
    coefficient_lasso[,i]<-tmp_y$s0
    
  }
  
  rownames(coefficient_ridge)<-rownames(ridge_reg1$beta)
  rownames(coefficient_lasso)<-rownames(lasso_model$beta)
  
  mean_vector_ridge<-rowMeans(as.matrix(coefficient_ridge))
  names(mean_vector_ridge)<-names(ridge_reg1$beta)
  
  mean_vector_lasso<-rowMeans(coefficient_lasso)
  names(mean_vector_lasso)<-names(lasso_model$beta)
  
  
  #  write.csv(coefficient_ridge,"G://embryo/20220407/ML/TFs_coefficient_ridge.csv")
  #  write.csv(coefficient_ridge,"G://embryo/20220407/ML/TFs_coefficient_lasso.csv")
  #  write.csv(as.data.frame(predict_value_ridge),"G://embryo/20220407/ML/predictedValues_ridge.csv")
  #  write.csv(as.data.frame(predict_value_lasso),"G://embryo/20220407/ML/predictedValues_lasso.csv")
  
  #  write.csv(mean_vector_ridge,"G://embryo/20220407/ML/mean_vector.ridge.csv")
  #  write.csv(mean_vector_lasso,"G://embryo/20220407/ML/mean_vector.lasso.csv")
  
  write.csv(coefficient_ridge,"G://embryo/20220904/ML/TFs_coefficient_ridge.csv")
  write.csv(coefficient_ridge,"G://embryo/20220904/ML/TFs_coefficient_lasso.csv")
  write.csv(as.data.frame(predict_value_ridge),"G://embryo/20220904/ML/predictedValues_ridge.csv")
  write.csv(as.data.frame(predict_value_lasso),"G://embryo/20220904/ML/predictedValues_lasso.csv")
  
  write.csv(mean_vector_ridge,"G://embryo/20220904/ML/mean_vector.ridge.csv")
  write.csv(mean_vector_lasso,"G://embryo/20220904/ML/mean_vector.lasso.csv")
}


# Compute R^2 from true and predicted values
eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  
  
  # Model performance metrics
  data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  )
  
}


cross_validaton<-function(k,index){
  
  k_1<-(length(index)%%k)
  quant<-((length(index)-k_1)/k)
  
  quant_x<-quant
  
  samplelist<-list()
  
  start<-1
  
  for(i in 1:k){
    
    print(paste(i,",",start,",",quant,sep=""))
    
    if(i==k){
      
      samplelist[[i]]<-index[start:(quant+k_1)]
      
    }else{
      
      samplelist[[i]]<-index[start:quant]
    }
    
    start<-(quant+1)
    quant<-(quant+quant_x)
    
  }
  
  return(samplelist)
  
}

#######################################################


# 
subtype="gut"

## prediction_VS_obs

#pred.<-read.csv("G://embryo/20220810/ML/predictedValues.csv",header=T)
pred.<-read.csv("G://embryo/20220904/gut/predictedValues.csv",header=T)
#pred.<-read.csv("G://embryo/20220810/ML/predictedValues_lasso.csv",header=T)
#obs.<-read.csv("G://embryo/20220810/TFs.data_Kupffer_ML_v2.csv",header=T)
obs.<-read.csv("G://embryo/20220904/TFs.data_gut_ML_v2.csv",header=T)
obs.<-obs.[,c(1,ncol(obs.))]

vs<-merge(pred.,obs.,by="X")

corr<-cor(vs$predict_value,vs$module.score)
corr<-round(corr,digits=4)

p<-ggplot(vs) + geom_point(aes(x=predict_value, y=module.score),shape = 21,stroke = 0.2,size = 3,fill="#AEC7E8FF")+
  theme_bw()+theme(panel.grid=element_blank())+theme(plot.title = element_text(size = 60, face = "bold"))+
  xlab("predicted score")+ylab("observated score")+geom_smooth(method="lm",aes(x=predict_value, y=module.score),se = F)+
  labs(subtitle = paste("Correlation:",corr,sep=""))

ggsave(p,filename = paste("G://embryo/20220904/",subtype,"_TFs_plot.png",sep=""),width = 10, height = 8, units = "in", device='png')
ggsave(p,filename = paste("G://embryo/20220904/",subtype,"_TFs_plot.pdf",sep=""),width = 10, height = 8, units = "in", device='pdf')


subtype="pram"

## prediction_VS_obs


pred.<-read.csv("G://embryo/20220904/PraM/predictedValues.csv",header=T)

obs.<-read.csv("G://embryo/20220904/TFs.data_PraM_ML_v2.csv",header=T)
obs.<-obs.[,c(1,ncol(obs.))]

vs<-merge(pred.,obs.,by="X")

corr<-cor(vs$predict_value,vs$module.score)
corr<-round(corr,digits=4)

p<-ggplot(vs) + geom_point(aes(x=predict_value, y=module.score),shape = 21,stroke = 0.2,size = 3,fill="#FF7F0EFF")+
  theme_bw()+theme(panel.grid=element_blank())+theme(plot.title = element_text(size = 60, face = "bold"))+
  xlab("predicted score")+ylab("observated score")+geom_smooth(method="lm",aes(x=predict_value, y=module.score),se = F)+
  labs(subtitle = paste("Correlation:",corr,sep=""))

ggsave(p,filename = paste("G://embryo/20220904/",subtype,"_TFs_plot.png",sep=""),width = 10, height = 8, units = "in", device='png')
ggsave(p,filename = paste("G://embryo/20220904/",subtype,"_TFs_plot.pdf",sep=""),width = 10, height = 8, units = "in", device='pdf')



subtype="microglia"

## prediction_VS_obs


pred.<-read.csv("G://embryo/20220904/microglia/predictedValues.csv",header=T)

obs.<-read.csv("G://embryo/20220904/TFs.data_microglia_ML_v2.csv",header=T)
obs.<-obs.[,c(1,ncol(obs.))]

vs<-merge(pred.,obs.,by="X")

corr<-cor(vs$predict_value,vs$module.score)
corr<-round(corr,digits=4)

p<-ggplot(vs) + geom_point(aes(x=predict_value, y=module.score),shape = 21,stroke = 0.2,size = 3,fill="#1F77B4FF")+
  theme_bw()+theme(panel.grid=element_blank())+theme(plot.title = element_text(size = 60, face = "bold"))+
  xlab("predicted score")+ylab("observated score")+geom_smooth(method="lm",aes(x=predict_value, y=module.score),se = F)+
  labs(subtitle = paste("Correlation:",corr,sep=""))

ggsave(p,filename = paste("G://embryo/20220904/",subtype,"_TFs_plot.png",sep=""),width = 10, height = 8, units = "in", device='png')
ggsave(p,filename = paste("G://embryo/20220904/",subtype,"_TFs_plot.pdf",sep=""),width = 10, height = 8, units = "in", device='pdf')

## coefficient

subtype<-"Microglia"

coeffi.<-read.csv("G://embryo/20220904/microglia/TFs_coefficent_mean.csv")
coeffi.<-coeffi.[-1,]

colnames(coeffi.)<-c("TFs","coefficient")

coeffi.<-coeffi.[order(coeffi.$coefficient,decreasing = T),]

coeffi.$TFs=factor(coeffi.$TFs,levels = rev(coeffi.$TFs),ordered = T)

p<-ggplot(coeffi., aes(x=TFs, y=coefficient)) + geom_bar(stat="identity", fill="#1F77B4FF")+coord_flip()+
  theme_bw()+theme(panel.grid=element_blank())+theme(plot.title = element_text(size = 60, face = "bold"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p,filename = paste("G://embryo/20220904/",subtype,"_TFs_coefficient_plot.png",sep=""),width = 10, height = 8, units = "in", device='png')
ggsave(p,filename = paste("G://embryo/20220904/",subtype,"_TFs_coefficient_plot.pdf",sep=""),width = 10, height = 8, units = "in", device='pdf')


subtype<-"PraM"

coeffi.<-read.csv("G://embryo/20220904/PraM/TFs_coefficent_mean.csv")
coeffi.<-coeffi.[-1,]

colnames(coeffi.)<-c("TFs","coefficient")

coeffi.<-coeffi.[order(coeffi.$coefficient,decreasing = T),]

coeffi.$TFs=factor(coeffi.$TFs,levels = rev(coeffi.$TFs),ordered = T)

p<-ggplot(coeffi., aes(x=TFs, y=coefficient)) + geom_bar(stat="identity", fill="#FF7F0EFF")+coord_flip()+
  theme_bw()+theme(panel.grid=element_blank())+theme(plot.title = element_text(size = 60, face = "bold"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p,filename = paste("G://embryo/20220904/",subtype,"_TFs_coefficient_plot.png",sep=""),width = 10, height = 8, units = "in", device='png')
ggsave(p,filename = paste("G://embryo/20220904/",subtype,"_TFs_coefficient_plot.pdf",sep=""),width = 10, height = 8, units = "in", device='pdf')


subtype<-"gut"

coeffi.<-read.csv("G://embryo/20220904/gut/TFs_coefficent_mean.csv")
coeffi.<-coeffi.[-1,]

colnames(coeffi.)<-c("TFs","coefficient")

coeffi.<-coeffi.[order(coeffi.$coefficient,decreasing = T),]

coeffi.$TFs=factor(coeffi.$TFs,levels = rev(coeffi.$TFs),ordered = T)

p<-ggplot(coeffi., aes(x=TFs, y=coefficient)) + geom_bar(stat="identity", fill="#AEC7E8FF")+coord_flip()+
  theme_bw()+theme(panel.grid=element_blank())+theme(plot.title = element_text(size = 60, face = "bold"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p,filename = paste("G://embryo/20220904/",subtype,"_TFs_coefficient_plot.png",sep=""),width = 10, height = 8, units = "in", device='png')
ggsave(p,filename = paste("G://embryo/20220904/",subtype,"_TFs_coefficient_plot.pdf",sep=""),width = 10, height = 8, units = "in", device='pdf')
