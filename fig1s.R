###fig1s

#####
###Fig1s A
#####
human_development = human_development[, metatable$Cell[metatable$time != "Adult"]]

process_cells_annotation(human_development)

plot_cells_annotation(human_development,type="histogram")
plot_cells_annotation(human_development,type="boxplot")
plot_UMIs_vs_Detected_genes(human_development)

#####
###Fig1s B
#####
mc_umi = read.csv("E:\\singlecellaR/cleanumi.csv",row.names = 1)
mc_umi = mc_umi[,unique(metatable$mc)]

the_mat = mc_umi[discard_gene2(rownames(mc_umi)),colnames(mc_umi)]
the_mat = CreateSeuratObject(counts = the_mat)
the_mat <- NormalizeData(the_mat, normalization.method =  "LogNormalize")
the_mat <- FindVariableFeatures(the_mat,selection.method = "vst", nfeatures = 2000)
the_mat <- ScaleData(the_mat)
the_mat <- RunPCA(object = the_mat, pc.genes = VariableFeatures(the_mat))
the_mat <- FindNeighbors(the_mat, dims = 1:30)
the_mat <- FindClusters(the_mat, resolution = 0.5)
the_mat <- RunUMAP(object = the_mat, dims = 1:30, do.fast = TRUE)

mc_anno = metatable[,c("mc","major")]
mc_anno = mc_anno[!duplicated(mc_anno),]

mc.id = mc_anno$major
names(mc.id) = mc_anno$mc
Idents(the_mat) = mc.id
mc.marker = FindAllMarkers(the_mat)
top20 <- mc.marker %>% group_by(cluster) %>% top_n(n =30, wt = avg_log2FC)
write.csv(top20,"E:\\sup_tab/all.csv")

top20 <- mc.marker %>% group_by(cluster) %>% top_n(n =20, wt = avg_log2FC)
p1 = DoHeatmap(the_mat,features = top20$gene)

tmp = as.data.frame(unique(metatable$mc))
colnames(tmp) = "mc"

tmp$sudomc = rownames(tmp)
all.lfp = read.csv("E:\\all/lfp_all_clean_re.csv",row.names = 1)
colnames(all.lfp) = tmp$mc

all.lfp2 = all.lfp[unique(top20$gene),unique(metatable$mc)]

col_order = c()
for(i in levels(p1$data$Identity)){
  col_order = c(col_order,unique(data1$mc[data1$major == i]))
}

tmp = metatable[,c("mc","major")]
tmp = tmp[!duplicated(tmp),]

ann_col.major = as.data.frame(tmp$major)
colnames(ann_col.major) = "major"
rownames(ann_col.major) = tmp$mc
p2 = pheatmap::pheatmap(all.lfp2[rev(levels(p1$data$Feature)),col_order],breaks = bk,
                        show_colnames = F,cluster_cols = F,border_color = "grey30",
                        cluster_rows = F,show_rownames = F,
                        annotation_col = ann_col.major,annotation_colors = list(major = major_color),
                        color = colorRampPalette(c("blue","white", "red"))(length(bk)))
ggsave("E:\\fig1/hc_all_mc.pdf",p2,width = 4,height = 4)

#####
### Fig1s C D E F & G
#####
mc.anno = metatable[,c("mc","subtype","major")]
mc.anno = mc.anno[!duplicated(mc.anno),]

####mono_dc
monodc_mc.umi = mc_umi[discard_gene2(rownames(mc_umi)),mc.anno$mc[mc.anno$major %in% c("monocyte","DC")]]
monodc_mc.umi = CreateSeuratObject(counts = monodc_mc.umi)
monodc_mc.umi <- NormalizeData(monodc_mc.umi, normalization.method =  "LogNormalize")
monodc_mc.umi <- FindVariableFeatures(monodc_mc.umi,selection.method = "vst", nfeatures = 2000)
monodc_mc.umi <- ScaleData(monodc_mc.umi)
monodc_mc.umi <- RunPCA(object = monodc_mc.umi, pc.genes = VariableFeatures(monodc_mc.umi))
monodc_mc.umi <- FindNeighbors(monodc_mc.umi, dims = 1:30)
monodc_mc.umi <- FindClusters(monodc_mc.umi, resolution = 0.5)
monodc_mc.umi <- RunUMAP(object = monodc_mc.umi, dims = 1:30, do.fast = TRUE)

monodc_mc.anno = mc.anno[mc.anno$major %in% c("monocyte","DC"),]
monodc_mc.ident = monodc_mc.anno$subtype
names(monodc_mc.ident) = monodc_mc.anno$mc

Idents(monodc_mc.umi) = monodc_mc.ident
Idents(monodc_mc.umi) = factor(Idents(monodc_mc.umi),levels = c("CD163+CD14+ DC","cycling DC","cDC1","cDC2","mDC","pDC",
                                                                "monocyte","cycling monocyte","non-classical monocyte"))

monodc_mc.markers = FindAllMarkers(monodc_mc.umi)
top20 <- monodc_mc.markers %>% group_by(cluster) %>% top_n(n =30, wt = avg_log2FC)

write.csv(top20,"E:\\sup_tab/monodc_mc_markers.csv")

p1 = DoHeatmap(monodc_mc.umi,top20$gene)

mcs = as.data.frame(levels(p1$data$Cell))
rownames(mcs) = mcs$`levels(p1$data$Cell)`
mcs = mcs$`levels(p1$data$Cell)`[mcs$`levels(p1$data$Cell)` %in% monodc_mc.anno$mc]

rownames(monodc_mc.anno) = monodc_mc.anno[,1]
monodc_mc.anno = monodc_mc.anno[,-1]

p2 = pheatmap::pheatmap(monodc.lfp[rev(levels(p1$data$Feature)),mcs],
                        cluster_rows = F,cluster_cols = F,
                        breaks = bk,
                        annotation_col = monodc_mc.anno,
                        annotation_colors = list(major = major_color,subtype = mono_dc.color),
                        show_rownames = F,show_colnames = F,legend = F,annotation_legend = F,
                        color = colorRampPalette(c("navy","white", "firebrick3"))(length(bk)))
ggsave("E:\\fig1/monodc_marker.pdf",p2,width = 6,height = 4,dpi = 300)

####T NK ILC
tnkilc.umi = mc_umi[discard_gene2(rownames(mc_umi)),mc.anno$mc[mc.anno$major %in% c("T","NK","ILC")]]
tnkilc.umi = CreateSeuratObject(counts = tnkilc.umi)
tnkilc.umi <- NormalizeData(tnkilc.umi, normalization.method =  "LogNormalize")
tnkilc.umi <- FindVariableFeatures(tnkilc.umi,selection.method = "vst", nfeatures = 2000)
tnkilc.umi <- ScaleData(tnkilc.umi)
tnkilc.umi <- RunPCA(object = tnkilc.umi, pc.genes = VariableFeatures(tnkilc.umi))
tnkilc.umi <- FindNeighbors(tnkilc.umi, dims = 1:30)
tnkilc.umi <- FindClusters(tnkilc.umi, resolution = 0.5)
tnkilc.umi <- RunUMAP(object = tnkilc.umi, dims = 1:30, do.fast = TRUE)

tnkilc_mc.anno = mc.anno[mc.anno$major %in% c("T","NK","ILC"),]
tnkilc_mc.ident = tnkilc_mc.anno$subtype
names(tnkilc_mc.ident) = tnkilc_mc.anno$mc

Idents(tnkilc.umi) = tnkilc_mc.ident
Idents(tnkilc.umi) = factor(Idents(tnkilc.umi),levels = c("naive T_S100A4+","naive T_CCR9+","naive T_S100A4-CCR9-","cytotoxic CD8","Treg","thymocyte",
                                                          "ILC_CXCR5+","ILC_ITGAE+","ILC_proliferating",
                                                          "NK_S100A4+S100A5+","NK_S100A4-S100A5-","NK_proliferating"))
tnkilc.markers = FindAllMarkers(tnkilc.umi)
top20 <- tnkilc.markers %>% group_by(cluster) %>% top_n(n =30, wt = avg_log2FC)
write.csv(top20,"E:\\sup_tab/tnkilc_markers.csv")


p1 = DoHeatmap(tnkilc.umi,unique(top20$gene))

mcs = as.data.frame(levels(p1$data$Cell))
rownames(mcs) = mcs$`levels(p1$data$Cell)`
mcs = mcs$`levels(p1$data$Cell)`[mcs$`levels(p1$data$Cell)` %in% tnkilc_mc.anno$mc]

rownames(tnkilc_mc.anno) = tnkilc_mc.anno[,1]
tnkilc_mc.anno = tnkilc_mc.anno[,-1]
tnkilc.lfp = read.csv("E:\\t/lfp_tilcnk.csv",row.names = 1)
load("E:\\t/tilcnk_sudomc.Rdata")
colnames(tnkilc.lfp) = tilcnk.sudomc$mc
names(t.color)[3] = "naive T_S100A4-CCR9-"
names(t.color)[6] = "naive T_CCR9+"
names(t.color)[4] = "naive T_S100A4+"
names(t.color)[11] = "ILC_ITGAE+"


p2 = pheatmap::pheatmap(tnkilc.lfp[rev(levels(p1$data$Feature)),mcs],
                        cluster_rows = F,cluster_cols = F,
                        breaks = bk,
                        annotation_col = tnkilc_mc.anno,
                        annotation_colors = list(major = major_color,subtype = t.color),
                        show_rownames = F,show_colnames = F,legend = F,annotation_legend =F,
                        color = colorRampPalette(c("navy","white", "firebrick3"))(length(bk)))
ggsave("E:\\fig1/tilc_marker.pdf",p2,width = 6,height = 4,dpi = 300)

###grani mk

granu_mk.umi = mc_umi[discard_gene3(rownames(mc_umi)),mc.anno$mc[mc.anno$major %in% c("granulocyte","MK")]]
granu_mk.umi = CreateSeuratObject(counts = granu_mk.umi)
granu_mk.umi <- NormalizeData(granu_mk.umi, normalization.method =  "LogNormalize")
granu_mk.umi <- FindVariableFeatures(granu_mk.umi,selection.method = "vst", nfeatures = 4000)
granu_mk.umi <- ScaleData(granu_mk.umi)
granu_mk.umi <- RunPCA(object = granu_mk.umi, pc.genes = VariableFeatures(granu_mk.umi))
granu_mk.umi <- FindNeighbors(granu_mk.umi, dims = 1:30)
granu_mk.umi <- FindClusters(granu_mk.umi, resolution = 0.5)
granu_mk.umi <- RunUMAP(object = granu_mk.umi, dims = 1:30, do.fast = TRUE)

granu_mk_mc.anno = mc.anno[mc.anno$major %in% c("granulocyte","MK"),]
granu_mk_mc.ident = granu_mk_mc.anno$subtype
names(granu_mk_mc.ident) = granu_mk_mc.anno$mc

Idents(granu_mk.umi) = granu_mk_mc.ident
Idents(granu_mk.umi) = factor(Idents(granu_mk.umi),levels = c("granu_CXCL2+CXCL3+","granu_CXCL2-CXCL3-","mast cell",
                                                              "eosinophil","basophil","neutrophil","MK_GATA2+FCER1A+ITGA4+",
                                                              "MK_GP5+GP6+ITGB3+","MK_ID2-GP5-GP6-ITGA4+"))
granu_mk.markers = FindAllMarkers(granu_mk.umi)
top20 <- granu_mk.markers %>% group_by(cluster) %>% top_n(n =30, wt = avg_log2FC)
write.csv(top20,"E:\\sup_tab/granu_mk_markers.csv")



p1 = DoHeatmap(granu_mk.umi,unique(top20$gene))
p1

mcs = as.data.frame(levels(p1$data$Cell))
rownames(mcs) = mcs$`levels(p1$data$Cell)`
mcs = mcs$`levels(p1$data$Cell)`[mcs$`levels(p1$data$Cell)` %in% granu_mk_mc.anno$mc]

rownames(granu_mk_mc.anno) = granu_mk_mc.anno[,1]
granu_mk_mc.anno = granu_mk_mc.anno[,-1]

granu_mk.lfp = read.csv("E:\\fig1/mk_granu/lfp_mk_granu.csv",row.names = 1)
colnames(granu_mk.lfp) = mk_granu.sudomc$mc

names(granu.color)[3] = "granu_CXCL2-CXCL3-"

p2 = pheatmap::pheatmap(granu_mk.lfp[rev(levels(p1$data$Feature)),mcs],
                        cluster_rows = F,cluster_cols = F,
                        breaks = bk,
                        annotation_col = granu_mk_mc.anno,
                        annotation_colors = list(major = major_color,subtype = c(mk.color,granu.color,eos)),
                        show_rownames = F,show_colnames = F,legend = F,annotation_legend =F,
                        color = colorRampPalette(c("navy","white", "firebrick3"))(length(bk)))
ggsave("E:\\fig1/granu_mk_marker.pdf",p2,width = 6,height = 4,dpi = 300)

###prog B

prog_b.umi = mc_umi[discard_gene2(rownames(mc_umi)),mc.anno$mc[mc.anno$major %in% c("B","progenitor")]]
prog_b.umi = CreateSeuratObject(counts = prog_b.umi)
prog_b.umi <- NormalizeData(prog_b.umi, normalization.method =  "LogNormalize")
prog_b.umi <- FindVariableFeatures(prog_b.umi,selection.method = "vst", nfeatures = 4000)
prog_b.umi <- ScaleData(prog_b.umi)
prog_b.umi <- RunPCA(object = prog_b.umi, pc.genes = VariableFeatures(prog_b.umi))
prog_b.umi <- FindNeighbors(prog_b.umi, dims = 1:30)
prog_b.umi <- FindClusters(prog_b.umi, resolution = 0.5)
prog_b.umi <- RunUMAP(object = prog_b.umi, dims = 1:30, do.fast = TRUE)

prog_b_mc.anno = mc.anno[mc.anno$major %in% c("B","progenitor"),]
prog_b_mc.ident = prog_b_mc.anno$subtype
names(prog_b_mc.ident) = prog_b_mc.anno$mc

Idents(prog_b.umi) = prog_b_mc.ident
Idents(prog_b.umi) = factor(Idents(prog_b.umi),levels = c("HSC/MPP","MPP","Mye Prog",
                                                          "Mk/E Prog","DC-MG Prog","Lym Prog",
                                                          "Pro_B", "Pre_B","Immature_B","Mature_B","plasma"))
prog_b.markers = FindAllMarkers(prog_b.umi)
top20 <- prog_b.markers %>% group_by(cluster) %>% top_n(n =30, wt = avg_log2FC)
write.csv(top20,"E:\\sup_tab/prog_b_markers.csv")



p1 = DoHeatmap(prog_b.umi,unique(top20$gene))
p1

mcs = as.data.frame(levels(p1$data$Cell))
rownames(mcs) = mcs$`levels(p1$data$Cell)`
mcs = mcs$`levels(p1$data$Cell)`[mcs$`levels(p1$data$Cell)` %in% prog_b_mc.anno$mc]

rownames(prog_b_mc.anno) = prog_b_mc.anno[,1]
prog_b_mc.anno = prog_b_mc.anno[,-1]

prog_b.lfp = read.csv("E:\\fig1/b_prog/lfp_bprog.csv",row.names = 1)
colnames(prog_b.lfp) = b_prog.sudomc$mc

p2 = pheatmap::pheatmap(prog_b.lfp[rev(levels(p1$data$Feature)),mcs],
                        cluster_rows = F,cluster_cols = F,
                        breaks = bk,
                        annotation_col = prog_b_mc.anno,
                        annotation_colors = list(major = major_color,subtype = c(b.color,prog.color)),
                        show_rownames = F,show_colnames = F,legend = F,annotation_legend =F,
                        color = colorRampPalette(c("navy","white", "firebrick3"))(length(bk)))
ggsave("E:\\fig1/prog_b.pdf",p2,width = 6,height = 4,dpi = 300)

###mac

mac.umi = mc_umi[discard_gene2(rownames(mc_umi)),mc.anno$mc[mc.anno$major == "macrophage"]]
mac.umi = CreateSeuratObject(counts = mac.umi)
mac.umi <- NormalizeData(mac.umi, normalization.method =  "LogNormalize")
mac.umi <- FindVariableFeatures(mac.umi,selection.method = "vst", nfeatures = 4000)
mac.umi <- ScaleData(mac.umi)
mac.umi <- RunPCA(object = mac.umi, pc.genes = VariableFeatures(mac.umi))
mac.umi <- FindNeighbors(mac.umi, dims = 1:30)
mac.umi <- FindClusters(mac.umi, resolution = 0.5)
mac.umi <- RunUMAP(object = mac.umi, dims = 1:30, do.fast = TRUE)

mac.umi.anno = mc.anno[mc.anno$major == "macrophage",]
mac.umi.ident = mac.umi.anno$subtype
names(mac.umi.ident) = mac.umi.anno$mc

Idents(mac.umi) = mac.umi.ident
Idents(mac.umi) = factor(Idents(mac.umi),levels = c("HdpM","YsdM_AFP_high","YsdM_AFP_low","pre_microglia",
                                                    "red_pulp","Kupffer_cell","pre_PraM","PraM","gonad_macrophage",
                                                    "Adrenalgland_macrophage","intestine CD209+ Mφ",
                                                    "intestine CD207+ Mφ","langerhans","osteoclast","microglia"))
mac.umi.markers = FindAllMarkers(mac.umi)
top20 <- mac.umi.markers %>% group_by(cluster) %>% top_n(n =30, wt = avg_log2FC)
write.csv(top20,"E:\\sup_tab/mac_markers.csv")


p1 = DoHeatmap(mac.umi,top20$gene)
p1

mcs = as.data.frame(levels(p1$data$Cell))
rownames(mcs) = mcs$`levels(p1$data$Cell)`
mcs = mcs$`levels(p1$data$Cell)`[mcs$`levels(p1$data$Cell)` %in% mac.umi.anno$mc]

rownames(mac.umi.anno) = mac.umi.anno[,1]
mac.umi.anno = mac.umi.anno[,-1]
mac.umi.anno2 = as.data.frame(mac.umi.anno$subtype)
rownames(mac.umi.anno2) = rownames(mac.umi.anno)

mac.lfp = read.csv("E:\\mac_new/lfp.csv",row.names = 1)
colnames(mac.lfp) = mac.new.sudomc$mc
write.csv(mac.lfp,"E:\\mac_new/lfp_real_mc.csv")
bk = seq(-2,2,0.5)
p2 = pheatmap::pheatmap(mac.lfp[rev(levels(p1$data$Feature)),mcs],
                        cluster_rows = F,cluster_cols = F,
                        breaks = bk,
                        annotation_col = mac.umi.anno,
                        annotation_colors = list(major = major_color,subtype = mac.color),
                        show_rownames = F,show_colnames = F,
                        legend = F,annotation_legend =F,
                        color = colorRampPalette(c("navy","white", "firebrick3"))(length(bk)))
ggsave("E:\\fig1/mac.pdf",p2,width = 6,height = 4,dpi = 300)

pheatmap::pheatmap(mac.lfp[rev(levels(p1$data$Feature)),mcs],
                   cluster_rows = F,cluster_cols = F,
                   breaks = bk,
                   annotation_col = mac.umi.anno,
                   annotation_colors = list(major = major_color,subtype = mac.color),
                   show_rownames = T,show_colnames = F,
                   legend = F,annotation_legend =F,
                   color = colorRampPalette(c("navy","white", "firebrick3"))(length(bk)))


#####
### Fig1s H
#####
library(randomcoloR)
embryo.color = randomColor(length(unique(sc$embryo)))
show_col(embryo.color)
tmp = sc

p1 = ggplot(tmp)+
  #geom_point(aes(x = x, y = y),color = "grey88",size = 2)+
  geom_point(aes(x = x, y = y,color = embryo),size = 2)+
  theme_bw()+
  scale_color_manual(values = embryo.color,breaks = unique(sc$embryo))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),legend.position = 'none')
ggsave(paste0("E:\\fig1/qc/","time_","all",".png"),p1,width = 18,height = 18,dpi = 300)

p1 = ggplot(tmp)+
  #geom_point(aes(x = x, y = y),color = "grey88",size = 2)+
  geom_point(aes(x = x, y = y,color = tissue),size = 2)+
  theme_bw()+
  scale_color_manual(values = tissue.color,breaks = unique(sc$tissue))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),legend.position = 'none')
ggsave(paste0("E:\\fig1/qc/","tissue_","all",".png"),p1,width = 18,height = 18,dpi = 300)

tmp$gender = NA
tmp$gender[tmp$embryo %in% c("embryo 2","embryo 28","embryo 23","embryo 1","embryo 5","embryo 7","embryo 9",
                             "embryo 11","embryo 13","embryo 16","embryo 17","embryo 18","embryo 19","embryo 22",
                             "embryo 25","embryo 26","embryo 27","embryo 35","embryo 40","embryo 45","embryo 54",
                             "embryo 58","embryo57","embryo 55")] = "female"
tmp$gender[tmp$embryo %in% c("embryo 4","embryo 8","embryo 6","embryo 12","embryo 14","embryo 15","embryo 20",
                             "embryo 21","embryo 24","embryo 29","embryo 30","embryo 31","embryo 32","embryo 33",
                             "embryo 34","embryo 36","embryo 37","embryo 38","embryo 39","embryo 41","embryo 43",
                             "embryo 44","embryo 46","embryo 47","embryo 50","embryo 52","embryo 56")] = "male"

tmp = na.omit(tmp)

p1 = ggplot(tmp)+
  #geom_point(aes(x = x, y = y),color = "grey88",size = 2)+
  geom_point(aes(x = x, y = y,color = gender),size = 2)+
  theme_bw()+
  #scale_color_manual(values = sb.color,breaks = unique(tmp$Seq_batch_ID))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),legend.position = 'none')
ggsave(paste0("E:\\fig1/qc/","gender_","all",".png"),p1,width = 18,height = 18,dpi = 300)


#####
### Fig1s I
#####
tmp = tmp[tmp$time != "Adult",]
for (i in unique(tmp$subtype)){
  dir.create(paste0("E:\\fig1/qc/qc2/",i))
  for (j in unique(metatable$time[metatable$subtype == i])) {
    dir.create(paste0("E:\\fig1/qc/qc2/",i,"/",j))
    tmp$x1 = NA
    tmp$y1 = NA
    
    tmp$x1[tmp$subtype == i & tmp$time == j] = tmp$x[tmp$subtype == i & tmp$time == j]
    tmp$y1[tmp$subtype == i & tmp$time == j] = tmp$y[tmp$subtype == i & tmp$time == j]
    p1 = ggplot(tmp)+
      geom_point(aes(x = x, y = y),color = "grey88",size = 2)+
      geom_point(data = tmp[tmp$subtype == i & tmp$time == j,],aes(x = x1, y = y1,color = embryo),size = 2)+
      theme_bw()+
      scale_color_manual(values = embryo.color,breaks = unique(sc$embryo))+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.border = element_blank(),
            panel.grid=element_blank(),legend.position = 'none')
    ggsave(paste0("E:\\fig1/qc/qc2/",i,"/",j,"/",i,"_",j,".png"),p1,width = 18,height = 18,dpi = 300)
    p1 = ggplot(tmp)+
      geom_point(aes(x = x, y = y),color = "grey88",size = 2)+
      geom_point(data = tmp[tmp$subtype == i & tmp$time == j,],aes(x = x1, y = y1,color = embryo),size = 2)+
      theme_bw()+
      scale_color_manual(values = embryo.color,breaks = unique(sc$embryo))+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.border = element_blank(),
            panel.grid=element_blank())
    ggsave(paste0("E:\\fig1/qc/qc2/",i,"/",j,"/",i,"_",j,".pdf"),p1,width = 10,height = 10,dpi = 300)
  }
}
tmp = tmp[tmp$time != "Adult",]
for (i in unique(tmp$tissue)){
  dir.create(paste0("E:\\fig1/qc/qc3/",i))
  for (j in unique(metatable$time[metatable$tissue == i])) {
    dir.create(paste0("E:\\fig1/qc/qc3/",i,"/",j))
    tmp$x1 = NA
    tmp$y1 = NA
    
    tmp$x1[tmp$tissue == i & tmp$time == j] = tmp$x[tmp$tissue == i & tmp$time == j]
    tmp$y1[tmp$tissue == i & tmp$time == j] = tmp$y[tmp$tissue == i & tmp$time == j]
    p1 = ggplot(tmp)+
      geom_point(aes(x = x, y = y),color = "grey88",size = 2)+
      geom_point(data = tmp[tmp$tissue == i & tmp$time == j,],aes(x = x1, y = y1,color = embryo),size = 2)+
      theme_bw()+
      scale_color_manual(values = embryo.color,breaks = unique(sc$embryo))+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.border = element_blank(),
            panel.grid=element_blank(),legend.position = 'none')
    ggsave(paste0("E:\\fig1/qc/qc3/",i,"/",j,"/",i,"_",j,".png"),p1,width = 18,height = 18,dpi = 300)
    p1 = ggplot(tmp)+
      geom_point(aes(x = x, y = y),color = "grey88",size = 2)+
      geom_point(data = tmp[tmp$tissue == i & tmp$time == j,],aes(x = x1, y = y1,color = embryo),size = 2)+
      theme_bw()+
      scale_color_manual(values = embryo.color,breaks = unique(sc$embryo))+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.border = element_blank(),
            panel.grid=element_blank())
    ggsave(paste0("E:\\fig1/qc/qc3/",i,"/",j,"/",i,"_",j,".pdf"),p1,width = 10,height = 10,dpi = 300)
  }
}




#####
### Fig1s J
#####
tmp2 = matrix(NA,ncol = length(unique(metatable$mc)),nrow = length(unique(metatable$embryo)))
colnames(tmp2) = unique(metatable$mc)
rownames(tmp2) = unique(metatable$embryo)

for (i in unique(metatable$embryo)){
  for (j in unique(metatable$mc)){
    tmp2[i,j] = nrow(metatable[metatable$mc == j & metatable$embryo == i,])/nrow(metatable[metatable$mc == j,])
  }
}

tmp2 = tmp2[setdiff(rownames(tmp2),c("Adult")),]
write.csv(tmp2,"E:\\fig1/qc/mc_embryo_ratio.csv")
tmp2 = read.csv("E:\\fig1/qc/mc_embryo_ratio.csv",row.names = 1)

embryo_order = c()
for (i in setdiff(time_order,c("Adult"))) {
  embryo_order = c(embryo_order,unique(metatable$embryo[metatable$time == i]))
}
mc_order = read.csv("E:\\allcells_realtime_noadult.csv")
mc_order = mc_order[mc_order$mc %in% colnames(tmp2),]
mc_order = mc_order[order(mc_order$the_median),]

bk =seq(0,0.1,0.01)
p1 = pheatmap::pheatmap(tmp2[embryo_order,mc_order$mc],show_colnames = F,breaks = bk,
                        cluster_rows = F,cluster_cols = F,border_color = "black",
                        color = colorRampPalette(c("white", "red"))(length(bk)))

ggsave("E:\\fig1/qc/mc_embryo_ratio.pdf",p1,width = 8,height =6.5,dpi = 300)