####Fig5 sup
#### writer： Feng Ruoqing
#### Date：2022-10-4

###Fig5s sup
mc_sc3 = mac12.metatable
mc_sc3$tissue = mc_sc3$subtype
mc_sc3$major_type = "macrophage"
n <- 100
ntime <- 150

pop <- 'macrophage'
tissueA <- "pre_PraM"; cutA <- 30
tissueB <- "PraM"; cutB <- 30

get_n_cells_per_pop_per_disease(mc_sc3, n, pop, tissueA, cutA, diag =F)
get_n_cells_per_pop_per_disease(mc_sc3, n, pop, tissueB, cutB, diag =F)

ind_diff <- diff_calc(mc_sc3, a_reduced,pop, n, ntime, tissueA, tissueB, cutA, cutB)
ind_diff_stat <- diff_stat(ind_diff, ntime)

f_cut <- 1.5
q_cut <- 0.05
pct_cut <- 0.05
ind_diff_stat_good <- ind_diff_stat[ind_diff_stat$pctA.mean > pct_cut | ind_diff_stat$pctB.mean > pct_cut,]
f <- log2((ind_diff_stat_good$avgB.mean + 0.01)/(ind_diff_stat_good$avgA.mean + 0.01))
q <- log(-(log(ind_diff_stat_good$wilcox.p.mean, 10)) + 1)
dat <- data.frame(g = rownames(ind_diff_stat_good), f = f, q = q)
#dat$f = -dat$f
up <- dat[dat$f > log2(f_cut) & dat$q  > log(-log10(q_cut) + 1),]
message('up: ', nrow(up))
dn <- dat[dat$f < -log2(f_cut) & dat$q  > log(-log10(q_cut) + 1),]
message('dn: ', nrow(dn))
rownames(dat) <- rownames(ind_diff_stat_good)
diff_genes <- rownames(dat)[abs(dat$f) > log2(f_cut) & dat$q > log(-log10(q_cut) + 1)]
message('diff: ', length(diff_genes))
dat$sig <- ""
dat[diff_genes,]$sig <- diff_genes
dat <- dat[order(dat$sig),]
dat$group = "non"
dat$group[dat$f > log2(f_cut) & dat$q > log(-log10(q_cut) + 1)] = "B"
dat$group[dat$f < -log2(f_cut) & dat$q > log(-log10(q_cut) + 1)] = "A"

deg = unique(dat$sig)
deg = c("C5AR1","CXCL3","CXCL8","DAB2","PTGS2","SOD2","TNF","IL1B","HMGA2","CD36","CD83","IL1A",
        "VEGFA","GLUL","ICAM1")

dat$sig = ""
for(i in deg){
  dat$sig[dat$g == i] = i
}

p1 = ggplot(dat, aes(f, q, label = sig,fill = group)) +
  geom_point( size=3,shape= 21,color = "grey66")+
  scale_fill_manual(values=c('#2CA02CFF','#FF7F0EFF', "grey88"))+
  geom_text_repel(data = dat[dat$sig != "",],max.overlaps = 100,size = 4) +
  geom_vline(xintercept = log2(f_cut) , linetype="dashed", color="grey66") + 
  geom_vline(xintercept = -log2(f_cut), linetype="dashed", color="grey66") + 
  ylab('log((-log10 q) + 1)') + xlab('log2 f')+
  annotate(geom="text", x=-2, y=3, label="pre_PraM",  fontface="bold",colour='#2CA02CFF', size=5)+
  annotate(geom="text", x=2, y=3, label="PraM",  fontface="bold",colour='#FF7F0EFF', size=5)+
  theme_bw()+theme_classic()+xlab("log2(fold change)")+ylab("q value")+
  theme(legend.position = "none",
        axis.text = element_text(size=14))
ggsave("E:\\fig5/deg_volcano_newnormalize.pdf",p1,width = 5,height = 5,dpi = 300)


###Fig5s B
figre5_mac_ratio = function(the_metatable,the_path,thew , theh){
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
  ggsave(the_path,the_p1,width =thew,height =theh,dpi = 300)
}
figre5_mac_ratio(the_metatable = mac12.metatable,the_path = "E:\\fig5/mac12_dynamic.pdf",thew = 8,theh = 5)

###Fig5s D

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

###Fig5s E
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