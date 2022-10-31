###Fig1
###
###

#####
###Fig1 B
#####
sc = read.csv("E:\\all/all_sc.csv",row.names = 1)
sc = merge(sc,metatable,by = "Well_ID")

show_col(pal_npg("nrc")(10))

major = c("DC","macrophage","monocyte","progenitor","granulocyte","MK","ILC","T","NK","B","erythrocyte")
major_color = c("#E64B35FF","#F39B7FFF","#B24745FF","#8491B4FF","#7E6148FF","#B09C85FF","#95CC5EFF","#00A087FF", "#91D1C2FF","#F79D1EFF","#69C8ECFF")
names(major_color) = major
save(major_color,file = "E:\\singlecellaR/major_color.Rdata")
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
ggsave("E:\\fig1/big_size/2d_noadult.pdf",p1,width = 18,height = 18,dpi = 300)



#####
###Fig1 C
#####
metatable2 = metatable
metatable2$major[metatable2$subtype == "microglial"] = "microglia"

tmp = metatable2$major
names(tmp) = metatable2$Well_ID

save(tmp,file = "E:\\tmp/lfp_major.Rdata")


metatable2 = metatable2[metatable2$major != "erythrocyte",]

marker = read.csv("E:\\tmp/all_marker.csv")
marker = marker$gene

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

lfp_major = read.csv("E:\\fig1/lfp_major.csv")
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
ggsave("E:\\fig1/marker_newnormalize.pdf",p1,width = 8,height = 4)

#####
###Fig D
#####
###T & ILC & NK
t.sc = read.csv("E:\\t/tilcnk_sc.csv")
t.sc = merge(t.sc,metatable,by = "Well_ID",all.x = T)

show_col(t.color)

t.subtype = c("cytotoxic CD8","thymocyte","naive T_IL32-S100A4-CCR9-","naive T_IL32+S100A4+","Treg","naive T_CCR9+",
              "NK_S100A4+S100A5+","NK_S100A4-S100A5-","NK_proliferating",
              "ILC_proliferating","ILC_ITGAE+TRDC+","ILC_CXCR5+")
t.color = c("#008EA0FF","#5A9599FF","#84D7E1FF","#ADE2D0FF","#3182BDFF","#46B8DAFF",
            "#A73030FF","#FF6348FF","#FF95A8FF",
            "#374E55FF","#80796BFF","#3D3B25FF")
names(t.color) = t.subtype
save(t.color,file = "E:\\fig1/t_color.Rdata")
t.sc = t.sc[t.sc$time != "Adult",]
p1 = ggplot(t.sc)+
  geom_point(aes(x = x , y = y, fill = subtype),shape = 21,color = "grey30",stroke = 0.2,size = 3)+
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
  scale_fill_manual(values = as.character(subcolor),breaks = names(subcolor))
ggsave("E:\\fig1/big_size/tilcnk.pdf",p1,width = 18,height = 18,dpi = 300)

###mono & dc 

mono_dc = read.csv("E:\\mono_dc/mono_dc_2d_transfer.csv",row.names = 1)
mono_dc = merge(mono_dc,metatable,by = "Well_ID")
mono_dc = mono_dc[mono_dc$time != "Adult",]

mono_dc.subtype = c("non-classical monocyte","CLEC4E-AQP9- monocyte","monocyte",
                    "pre DC","DC precursor","pDC","cDC1","mDC","cDC2")
mono_dc.color = c("#84BD00FF","#00AF66FF","#526E2DFF",
                  "#CC0C00FF","#5C88DAFF","#FFCD00FF","#7C878EFF","#00B5E2FF","#E89242FF")
names(mono_dc.color) = mono_dc.subtype
save(mono_dc.color,file = "E:\\fig1/monodc_color.Rdata")

show_col(pal_startrek()(7))


p1 = ggplot(mono_dc)+
  geom_point(aes(x = x , y = y, fill = subtype),shape = 21,color = "grey30",stroke = 0.2,,size = 3)+
  theme_bw()+
  xlab("")+ylab("")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),legend.position = 'none')+
  scale_fill_manual(values = as.character(subcolor),breaks = names(subcolor))
ggsave("E:\\fig1/big_size/monodc.pdf",p1,width = 18,height = 18,dpi = 300)

###B & prog
b_prog.sc = read.csv("E:\\fig1/b_prog/sc_bprog.csv",row.names = 1)
b_prog.sc = merge(b_prog.sc,metatable,by = "Well_ID",all.x = T)

b_prog.sudomc = b_prog[,c("mc","sudomc")]
b_prog.sudomc = b_prog.sudomc[!duplicated(b_prog.sudomc),]

subcolor = c(t.color,b.color,prog.color,mono_dc.color,mk.color,granu.color)

p1 = ggplot(b_prog.sc)+
  geom_point(aes(x = x,y = y,fill = subtype),shape = 21,color = "grey40",size  = 3,alpha = 0.8)+
  theme_bw()+
  xlab("")+ylab("")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank())+
  scale_fill_manual(values = as.character(subcolor),breaks = names(subcolor))
ggsave("E:\\fig1/b_prog/b_prog.png",p1,width = 8,height = 6.5)

###mk & granulocyte
mk_granu.sc = read.csv("E:\\fig1/mk_granu/sc_mk_granu.csv",row.names = 1)
mk_granu.sc = merge(mk_granu.sc,metatable,by = "Well_ID",all.x = T)

eos = "#FFB5A7FF"
names(eos) = "eosinophil"
subcolor = c(subcolor,eos)
mk_granu.sc = na.omit(mk_granu.sc)

mk_granu.sc = mk_granu.sc[mk_granu.sc$time != "Adult",]
p1 = ggplot(mk_granu.sc)+
  geom_point(aes(x = x,y = y,fill = subtype),shape = 21,color = "grey40",size  = 3,stroke = 0.2)+
  theme_bw()+
  xlab("")+ylab("")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),legend.position = 'none')+
  scale_fill_manual(values = as.character(subcolor),breaks = names(subcolor))
ggsave("E:\\fig1/big_size/mk_granu.pdf",p1,width = 18,height = 18,dpi = 300)

### mac
sc =  read.csv("E:\\mac_new/sc.csv",row.names = 1)
sc = merge(sc,metatable,by = "Well_ID")
sc = sc[sc$time != "Adult",]
p1 = ggplot(sc)+
  geom_point(aes(x = x, y = y,fill = subtype),size =2,shape =21 ,color = "grey30",stroke = 0.2)+
  geom_point(data = sc[sc$subtype == "Kupffer_cell",],aes(x = x, y = y,fill = subtype),size =2,shape =21 ,color = "grey30",stroke = 0.2)+
  geom_point(data = sc[sc$subtype == "YsdM_AFP_low",],aes(x = x, y = y,fill = subtype),size =2,shape =21 ,color = "grey30",stroke = 0.2)+
  geom_point(data = sc[sc$subtype %in% c("pre_PraM"),],aes(x = x, y = y,fill = subtype),size =2,shape =21 ,color = "grey30",stroke = 0.2)+
  geom_point(data = sc[sc$subtype == "pre_microglia",],aes(x = x, y = y,fill = subtype),size =2,shape =21 ,color = "grey30",stroke = 0.2)+
  geom_point(data = sc[sc$subtype %in% c("Adrenalgland_macrophage","PraM"),],aes(x = x, y = y,fill = subtype),size =2,shape =21 ,color = "grey30",stroke = 0.2)+
  geom_point(data = sc[sc$mc %in% c("X1253","X1258","X1266","X1272","X2858"),],aes(x = x, y = y,fill = subtype),size =2,shape =21 ,color = "grey30",stroke = 0.2)+
  theme_bw()+
  xlab("")+ylab("")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),legend.position = 'none')+
  scale_fill_manual(values = as.character(mac.color),breaks = names(mac.color))
p1
ggsave("E:\\mac_new/mac_2d.png",p1,width = 18,height =18,dpi = 300)
ggsave("E:\\fig1/big_size/mac_2d2.png",p1,width = 18,height =18,dpi = 300)