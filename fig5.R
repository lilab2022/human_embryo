###pram & pre-pram
###writen by Feng Ruoqing
###date: 2022-10-4

unique(mac.metatable$time)
mac.metatable = metatable[metatable$major =="macrophage",]
mac.metatable = mac.metatable[mac.metatable$time != "Adult",]

mac12.metatable = metatable[metatable$subtype %in% c("pre_PraM","PraM"),]

mac12.mc.anno = mac12.metatable[,c("mc","subtype")]
mac12.mc.anno = mac12.mc.anno[!duplicated(mac12.mc.anno),]

keep_cells = mac12.metatable$Well_ID
save(keep_cells,file = "E:\\fig5/mac12.Rdata")

tmp = as.data.frame(unique(mac12.metatable$mc))
tmp$sudomc = rownames(tmp)
colnames(tmp)[1] = "mc"
write.csv(tmp,"E:\\fig5/mac12_sudomc.csv")

mac12.metatable = merge(mac12.metatable,tmp,by = "mc",all.x = T)
mac12.metatable.mcmc = as.integer(mac12.metatable$sudomc)
names(mac12.metatable.mcmc) = mac12.metatable$Well_ID
save(mac12.metatable.mcmc,file = "E:\\fig5/mac12.metatable_mcmc.Rdata")
##### Fig5 A
mac12.sc = read.csv("E:\\fig5/sc_test.csv",row.names = 1)
mac12.sc = merge(mac12.sc,metatable,by = "Well_ID")
mac12.sc = mac12.sc[mac12.sc$time != "Adult",]
p1 = ggplot(mac12.sc,aes(x = x, y = y,fill = subtype))+
  geom_point(shape = 21,stroke = 0.2,size = 3,color = "grey30")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        legend.position = "none")+
  scale_fill_manual(values = as.character(mac.color),breaks = names(mac.color))
ggsave("E:\\fig5/mac12.png",p1,width = 18,height = 18,dpi = 300)

#####Fig5 C
mac.lfp = read.csv("E:\\mac_new/lfp_real_mc.csv",row.names = 1)

mac12.lfp_maclevel = mac.lfp[,unique(mac12.metatable$mc)]
mac12.lfp_maclevel = mac12.lfp_maclevel[c(module$gene[module$category == "PraM"]),]
mac12.lfp_maclevel["total",] = colMeans(mac12.lfp_maclevel)

for (i in  module$gene[module$category == "PraM"]){
tmp = as.data.frame(t(mac12.lfp_maclevel[c(i,"total"),]))
colnames(tmp) = c("gene","total")
tmp$mc = rownames(tmp)
tmp = merge(tmp,mac12.mc.anno,by = "mc")

p1 = ggplot(tmp,aes(x = total, y = gene,fill = subtype))+
  geom_point(size = 4,shape = 21,color = "grey30",stroke = 0.2)+
  theme_bw()+theme_classic()+
  xlab("")+ylab(i)+
  theme(panel.border = element_blank(),
        panel.grid=element_blank(),
        legend.position = "none")+
  scale_fill_manual(values = as.character(mac.color),breaks = names(mac.color))
ggsave(paste0("E:\\fig5/mac1_module/",i,".png"),p1,width = 4.5,height = 4,dpi = 300)
}

####
mc_sc3 = mac12.metatable
mc_sc3$tissue = mc_sc3$subtype
save(mc_sc3,file = "E:\\fig5/mc_sc3.Rdata")

##### Fig5 B
##2d gene exp projection

use.gene = c(module$gene[module$category == "PraM"],"CD209")
tmp = as.data.frame(t(as.matrix(a_reduced[intersect(use.gene,rownames(a_reduced)),intersect(colnames(a_reduced),mac12.metatable$Well_ID)])))
tmp$Well_ID = rownames(tmp)

mac1_2_2d = merge(mac12.sc,tmp,by = "Well_ID",all.x = T)
#mac1_2_2d[is.na(mac1_2_2d)] = 0

for (i in intersect(use.gene,rownames(a_reduced))){
  print(i)
  mac1_2_2d$tmp = mac1_2_2d[,i]
  mac1_2_2d$tmp = log2(mac1_2_2d$tmp+1)
  mac1_2_2d$tmp[mac1_2_2d$tmp >6] = 6
  mac1_2_2d[46048,"tmp"] = 6
  p1 = ggplot(mac1_2_2d)+
    geom_point(aes(x = x,y = y),color = "grey88",size = 3)+
    geom_point(data = mac1_2_2d[mac1_2_2d$tmp >0,],aes(x = x,y = y,color = tmp),size = 3)+
    scale_color_gradientn(values = seq(0,1,0.2),colours = c("pink", "orange", "red", "brown", "black"))+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank(),legend.position = 'none')
  ggsave(paste0("E:\\fig5/gene_exp_2d_new_normalize/",i,".png"),p1,width = 18,height = 18,dpi = 300)
}

use.gene = module$gene[module$category == "PraM"]
tmp = as.data.frame(t(as.matrix(a_reduced[intersect(use.gene,rownames(a_reduced)),intersect(colnames(a_reduced),mac12.metatable$Well_ID)])))
tmp$Well_ID = rownames(tmp)

mac1_2_2d = merge(mac12.sc,tmp,by = "Well_ID",all.x = T)

mac1_2_2d$sum = rowMeans(log2(mac1_2_2d[,intersect(use.gene,rownames(a_reduced))]+1))
mac1_2_2d$sum[mac1_2_2d$sum >3] = 3
mac1_2_2d$sum[mac1_2_2d$sum <=1] = 1
p1 = ggplot(mac1_2_2d)+
  geom_point(aes(x = x,y = y),,color = "grey66",size = 3)+
  geom_point(data = mac1_2_2d[mac1_2_2d$sum >0,],aes(x = x,y = y,color = sum),size = 3)+
  theme_bw()+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("pink", "orange", "red", "brown", "black"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        legend.position = "none")
ggsave("E:\\fig5/pram_module_new_normalize.png",p1,width = 18,height = 18,dpi = 300)

p1 = ggplot(mac1_2_2d)+
  geom_point(aes(x = x,y = y),,color = "grey66",size = 3)+
  geom_point(data = mac1_2_2d[mac1_2_2d$sum >0,],aes(x = x,y = y,color = sum),size = 3)+
  theme_bw()+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("pink", "orange", "red", "brown", "black"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank())
ggsave("E:\\fig5/pram_module_new_normalize_l.pdf",p1,width = 18,height = 18,dpi = 300)


##### Fig5 H
##FACS marker
target.gene = c("MRC1","P2RY12","MMP9","FCGR3A")

mac1_2 = metatable[metatable$subtype %in% c("PraM","pre_PraM","microglia","osteoclast","gonad_macrophage","Adrenalgland_macrophage"),]

sub.ord = c("PraM","pre_PraM","microglia","osteoclast","gonad_macrophage","Adrenalgland_macrophage")

tmp = a_reduced[target.gene,]
tmp = as.data.frame(t(as.matrix(tmp)))
tmp$Well_ID = rownames(tmp)
mac1_2 = merge(mac1_2,tmp,by = "Well_ID",all.x = T)
mac1_2[is.na(mac1_2)] = 0
mac1_2 = mac1_2[,c("subtype",target.gene)]
mac1_2 = reshape2::melt(mac1_2)
mac1_2 = mac1_2[mac1_2$subtype != "0",]
mac1_2$subtype = factor(mac1_2$subtype,levels = sub.ord,ordered = T)
p1 = ggplot(mac1_2)+
  geom_violin(aes(x = subtype, y = log2(value+1),fill = subtype,color = subtype),scale = "width")+
  facet_grid(variable ~ . , scales = "free")+
  xlab("")+ylab("")+theme_bw()+
  scale_x_discrete("")+
  scale_y_continuous(expand = c(0, 0), labels = function(x)
    c(rep(x = '', times = length(x) - 2), x[length(x) - 1], '')) +
  scale_fill_manual(values = as.character(mac.color), breaks = names(mac.color))+
  scale_color_manual(values = as.character(mac.color), breaks = names(mac.color))+
  theme(strip.background = element_rect(color = "black", fill = "white"),
        #axis.text.x = element_blank(),
        panel.spacing = unit(x = 0, units = 'lines'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(size = 12,angle = 360, face = 'bold'),
        legend.position = "none") 
ggsave("E:\\fig5/facs_marker_new_normalize.pdf",p1,width = 5,height = 4,dpi = 300)



###BAR CHART Fig5 G

i = "Male gonad"
tmp = mac.metatable[mac.metatable$tissue == i,]
tmp = tmp[tmp$time %in% after9,]

tmp2 = as.data.frame(table(tmp$subtype))
tmp2 = tmp2[order(tmp2$Freq,decreasing = T),]
tmp2$Var1 = as.character(tmp2$Var1)
tmp2 = tmp2[tmp2$Freq >100,]
tmp2$ratio = round(tmp2$Freq/sum(tmp2$Freq),2) 
tmp2$sample = i
tissue_ratio = tmp2

for (i in c("Lung","Kidney","Heart","Adrenalgland","Female gonad","Skin","Brain")) {
  tmp = mac.metatable[mac.metatable$tissue == i,]
  tmp = tmp[tmp$time %in% after9,]
  
  tmp2 = as.data.frame(table(tmp$subtype))
  tmp2 = tmp2[order(tmp2$Freq,decreasing = T),]
  tmp2 = tmp2[tmp2$Freq >100,]
  rownames(tmp2) = NULL
  tmp2$Var1 = as.character(tmp2$Var1)
  
  if(nrow(tmp2) > 3){
    tmp2_1 = as.data.frame(tmp2[c(1:3),])
    tmp2_1$Var1 = as.character(tmp2_1$Var1)
    tmp2_1[4,"Var1"] = "other"
    tmp2_1[4,"Freq"] = sum(tmp2$Freq[4:nrow(tmp2)])
    tmp2 = tmp2_1
  }
  
  tmp2$ratio = round(tmp2$Freq/sum(tmp2$Freq),2) 
  tmp2$sample = i
  tissue_ratio = rbind(tissue_ratio,tmp2)
}

tissue_ratio = tissue_ratio[tissue_ratio$sample %in% c("Lung","Kidney","Heart","Adrenalgland","Female gonad",
                                                       "Male gonad","Skin","Brain"),]
tissue_ratio$sample = factor(tissue_ratio$sample,levels = c("Female gonad","Kidney","Heart","Lung","Skin","Male gonad","Adrenalgland","Brain"),ordered = T)
tissue_ratio$Var1 = factor(tissue_ratio$Var1,levels = c("PraM","pre_PraM","gonad_macrophage","osteoclast","microglia",
                                                        "YsdM_AFP_low","Adrenalgland_macrophage","other"),ordered = T)
p1 = ggplot(tissue_ratio)+
  geom_bar(aes(x= sample,y=ratio,fill=Var1),stat="identity",position = "fill",width = 0.8,color = "grey30")+
  xlab("")+ylab("")+
  scale_y_continuous(expand=c(0,0),limits = c(0,1))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.text = element_text(size = 12),
        axis.title=element_blank(),
        axis.text.y = element_blank(),
        legend.title=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(),
        legend.position = "none",
        axis.ticks=element_line())+
  scale_fill_manual(values = c(as.character(mac.color),"grey66"),breaks = c(names(mac.color),"other"))
ggsave("E:\\fig5/pram_prepram_ratio.pdf",p1,width = 5,height = 4,dpi = 300)


#####
data = read.csv("E:\\fig5/counting_tab/ovary.csv",row.names = 1)
colnames(data) = c("CD45+MRC1+","CD45+MRC1-")

data2 = data
data2$range = rownames(data2)
data2 = reshape2::melt(data2)

p1 = ggplot(data2,aes(x = range,y = value,group = variable))+
  geom_point(aes(fill = variable),shape = 21,size = 4,color = "grey66")+
  geom_smooth(aes(color = variable),method = "lm",se = F,size = 1.5)+
  stat_poly_eq(
    aes(color = variable,label = paste(..eq.label..,  sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 4, #公式字体大小
    label.x = 0.2,  #位置 ，0-1之间的比例
    #label.y = 0.95
  )+
  theme_bw()+theme_classic()+
  xlab("")+ylab("")+
  scale_color_manual(values = c("navy","firebrick3"))+
  scale_fill_manual(values = c("navy","firebrick3"))+
  theme(legend.position = "none")

p2 = ggplot(data2,aes(x = range,y = value,color = variable,group = variable))+
  geom_line(size = 1.5)+
  theme_bw()+theme_classic()+
  xlab("")+ylab("")+
  scale_color_manual(values = c("navy","firebrick3"))+
  theme(legend.position = "none",
        axis.text.x=element_blank())

vp <- viewport(width = 0.5, height = 0.35, x = 0.95,y = 0.82,just=c("right","top"))

pdf("E:\\fig5/counting_fig/ovary.pdf",width = 5,height = 6)
print(p1)
print(p2,vp=vp)
dev.off()

#####
data = read.csv("E:\\fig5/counting_tab/lung.csv",row.names = 1)
colnames(data) = c("CD45+MRC1+","CD45+MRC1-")

data2 = data
data2$range = rownames(data2)
data2 = reshape2::melt(data2)

p1 = ggplot(data2,aes(x = range,y = value,group = variable))+
  geom_point(aes(fill = variable),shape = 21,size = 4,color = "grey66")+
  geom_smooth(aes(color = variable),method = "lm",se = F,size = 1.5)+
  stat_poly_eq(
    aes(color = variable,label = paste(..eq.label..,  sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 4, #公式字体大小
    label.x = 0.2,  #位置 ，0-1之间的比例
    #label.y = 0.95
  )+
  theme_bw()+theme_classic()+
  xlab("")+ylab("")+
  scale_color_manual(values = c("navy","firebrick3"))+
  scale_fill_manual(values = c("navy","firebrick3"))+
  theme(legend.position = "none")

p2 = ggplot(data2,aes(x = range,y = value,color = variable,group = variable))+
  geom_line(size = 1.5)+
  theme_bw()+theme_classic()+
  xlab("")+ylab("")+
  scale_color_manual(values = c("navy","firebrick3"))+
  theme(legend.position = "none",
        axis.text.x=element_blank())

vp <- viewport(width = 0.5, height = 0.35, x = 1,y = 1,just=c("right","top"))

pdf("E:\\fig5/counting_fig/lung.pdf",width = 5,height = 6)
print(p1)
print(p2,vp=vp)
dev.off()

#####
data = read.csv("E:\\fig5/counting_tab/adrenalgland.csv",row.names = 1)
colnames(data) = c("CD45+MRC1+","CD45+MRC1-")

data2 = data
data2$range = rownames(data2)
data2 = reshape2::melt(data2)

p1 = ggplot(data2,aes(x = range,y = value,group = variable))+
  geom_point(aes(fill = variable),shape = 21,size = 4,color = "grey66")+
  geom_smooth(aes(color = variable),method = "lm",se = F,size = 1.5)+
  stat_poly_eq(
    aes(color = variable,label = paste(..eq.label..,  sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 4, #公式字体大小
    label.x = 0.2,  #位置 ，0-1之间的比例
    #label.y = 0.95
  )+
  theme_bw()+theme_classic()+
  xlab("")+ylab("")+
  scale_color_manual(values = c("navy","firebrick3"))+
  scale_fill_manual(values = c("navy","firebrick3"))+
  theme(legend.position = "none")

p2 = ggplot(data2,aes(x = range,y = value,color = variable,group = variable))+
  geom_line(size = 1.5)+
  theme_bw()+theme_classic()+
  xlab("")+ylab("")+
  scale_color_manual(values = c("navy","firebrick3"))+
  theme(legend.position = "none",
        axis.text.x=element_blank())

vp <- viewport(width = 0.5, height = 0.35, x = 1,y = 1,just=c("right","top"))

pdf("E:\\fig5/counting_fig/adrenalgland.pdf",width = 5,height = 6)
print(p1)
print(p2,vp=vp)
dev.off()

#####
data = read.csv("E:\\fig5/counting_tab/kidney.csv",row.names = 1)
colnames(data) = c("CD45+MRC1+","CD45+MRC1-")

data2 = data
data2$range = rownames(data2)
data2 = reshape2::melt(data2)

p1 = ggplot(data2,aes(x = range,y = value,group = variable))+
  geom_point(aes(fill = variable),shape = 21,size = 4,color = "grey66")+
  geom_smooth(aes(color = variable),method = "lm",se = F,size = 1.5)+
  stat_poly_eq(
    aes(color = variable,label = paste(..eq.label..,  sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 4, #公式字体大小
    label.x = 0.2,  #位置 ，0-1之间的比例
    #label.y = 0.95
  )+
  theme_bw()+theme_classic()+
  xlab("")+ylab("")+
  scale_color_manual(values = c("navy","firebrick3"))+
  scale_fill_manual(values = c("navy","firebrick3"))+
  theme(legend.position = "none")

p2 = ggplot(data2,aes(x = range,y = value,color = variable,group = variable))+
  geom_line(size = 1.5)+
  theme_bw()+theme_classic()+
  xlab("")+ylab("")+
  scale_color_manual(values = c("navy","firebrick3"))+
  theme(legend.position = "none",
        axis.text.x=element_blank())

vp <- viewport(width = 0.5, height = 0.35, x = 1,y = 1,just=c("right","top"))

pdf("E:\\fig5/counting_fig/kidney.pdf",width = 5,height = 6)
print(p1)
print(p2,vp=vp)
dev.off()

#####
data = read.csv("E:\\fig5/counting_tab/skin.csv",row.names = 1)
colnames(data) = c("CD45+MRC1+","CD45+MRC1-")

data2 = data
data2$range = rownames(data2)
data2 = reshape2::melt(data2)

p1 = ggplot(data2,aes(x = range,y = value,group = variable))+
  geom_point(aes(fill = variable),shape = 21,size = 4,color = "grey66")+
  geom_smooth(aes(color = variable),method = "lm",se = F,size = 1.5)+
  stat_poly_eq(
    aes(color = variable,label = paste(..eq.label..,  sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 4, #公式字体大小
    label.x = 0.2,  #位置 ，0-1之间的比例
    #label.y = 0.95
  )+
  theme_bw()+theme_classic()+
  xlab("")+ylab("")+
  scale_color_manual(values = c("navy","firebrick3"))+
  scale_fill_manual(values = c("navy","firebrick3"))+
  theme(legend.position = "none")

p2 = ggplot(data2,aes(x = range,y = value,color = variable,group = variable))+
  geom_line(size = 1.5)+
  theme_bw()+theme_classic()+
  xlab("")+ylab("")+
  scale_color_manual(values = c("navy","firebrick3"))+
  theme(legend.position = "none",
        axis.text.x=element_blank())

vp <- viewport(width = 0.5, height = 0.35, x = 1,y = 1,just=c("right","top"))

pdf("E:\\fig5/counting_fig/skin.pdf",width = 5,height = 6)
print(p1)
print(p2,vp=vp)
dev.off()

#####
data = read.csv("E:\\fig5/counting_tab/Testicle.csv",row.names = 1)
colnames(data) = c("CD45+MRC1+","CD45+MRC1-")

data2 = data
data2$range = rownames(data2)
data2 = reshape2::melt(data2)

p1 = ggplot(data2,aes(x = range,y = value,group = variable))+
  geom_point(aes(fill = variable),shape = 21,size = 4,color = "grey66")+
  geom_smooth(aes(color = variable),method = "lm",se = F,size = 1.5)+
  stat_poly_eq(
    aes(color = variable,label = paste(..eq.label..,  sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 4, #公式字体大小
    label.x = 0.2,  #位置 ，0-1之间的比例
    #label.y = 0.95
  )+
  theme_bw()+theme_classic()+
  xlab("")+ylab("")+
  scale_color_manual(values = c("navy","firebrick3"))+
  scale_fill_manual(values = c("navy","firebrick3"))+
  theme(legend.position = "none")

p2 = ggplot(data2,aes(x = range,y = value,color = variable,group = variable))+
  geom_line(size = 1.5)+
  theme_bw()+theme_classic()+
  xlab("")+ylab("")+
  scale_color_manual(values = c("navy","firebrick3"))+
  theme(legend.position = "none",
        axis.text.x=element_blank())

vp <- viewport(width = 0.5, height = 0.35, x = 1,y = 1,just=c("right","top"))

pdf("E:\\fig5/counting_fig/Testicle.pdf",width = 5,height = 6)
print(p1)
print(p2,vp=vp)
dev.off()

#####
data = read.csv("E:\\fig5/counting_tab/brain.csv",row.names = 1)
colnames(data) = c("CD45+MRC1+","CD45+MRC1-")

data2 = data
data2$range = rownames(data2)
data2 = reshape2::melt(data2)
data2$dis = rep(c(10,20,30,40,50,60,70,80,90,100),2)
p1 = ggplot(data2,aes(x = range,y = value,group = variable))+
  geom_point(aes(fill = variable),shape = 21,size = 4,color = "grey66")+
  geom_smooth(aes(color = variable),method = "lm",se = F,size = 1.5,ymin = 0)+
  stat_poly_eq(
    aes(color = variable,label = paste(..eq.label..,  sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 4, #公式字体大小
    label.x = 0.2,  #位置 ，0-1之间的比例
    #label.y = 0.95
  )+
  #xlim(0,300)+
  theme_bw()+theme_classic()+
  xlab("")+ylab("")+ylim(0,300)+
  scale_color_manual(values = c("navy","firebrick3"))+
  scale_fill_manual(values = c("navy","firebrick3"))+
  theme(legend.position = "none")

p2 = ggplot(data2,aes(x = range,y = value,color = variable,group = variable))+
  geom_line(size = 1.5)+
  theme_bw()+theme_classic()+
  xlab("")+ylab("")+
  scale_color_manual(values = c("navy","firebrick3"))+
  theme(legend.position = "none",
        axis.text.x=element_blank())

vp <- viewport(width = 0.5, height = 0.35, x = 1,y = 1,just=c("right","top"))

pdf("E:\\fig5/counting_fig/brain.pdf",width = 5,height = 6)
print(p1)
print(p2,vp=vp)
dev.off()

#####
data = read.csv("E:\\fig5/counting_tab/heart.csv",row.names = 1)
colnames(data) = c("CD45+MRC1+","CD45+MRC1-")

data2 = data
data2$range = rownames(data2)
data2 = reshape2::melt(data2)

p1 = ggplot(data2,aes(x = range,y = value,group = variable))+
  geom_point(aes(fill = variable),shape = 21,size = 4,color = "grey66")+
  geom_smooth(aes(color = variable),method = "lm",se = F,size = 1.5)+
  stat_poly_eq(
    aes(color = variable,label = paste(..eq.label..,  sep = '~~~~')),
    formula = y ~ x,  parse = TRUE,
    size = 4, #公式字体大小
    label.x = 0.8,  #位置 ，0-1之间的比例
    #label.y = 0.95
  )+
  theme_bw()+theme_classic()+
  xlab("")+ylab("")+
  #xlim(0,250)+
  scale_color_manual(values = c("navy","firebrick3"))+
  scale_fill_manual(values = c("navy","firebrick3"))+
  theme(legend.position = "none")

p2 = ggplot(data2,aes(x = range,y = value,color = variable,group = variable))+
  geom_line(size = 1.5)+
  theme_bw()+theme_classic()+
  xlab("")+ylab("")+
  scale_color_manual(values = c("navy","firebrick3"))+
  theme(legend.position = "none",
        axis.text.x=element_blank())

vp <- viewport(width = 0.5, height = 0.35, x = 0.1,y = 0.1,just=c("left","bottom"))

pdf("E:\\fig5/counting_fig/heart.pdf",width = 5,height = 6)
print(p1)
print(p2,vp=vp)
dev.off()


ggsave("E:\\fig5/counting_fig/legend.pdf",p1,width = 4,height = 3,dpi = 300)

###fig5 J
library(metacell)
library(stringr)
library(foreach)
library(pheatmap)
library(tgconfig)
library(ggsignif)
library(Seurat)

set_param("scm_spike_regexp","^ERCC-","metacell")
set_param("scm_mc_mark_k_per_clust",100,"metacell") #default: 5
set_param("scm_mc_mark_min_gene_cov",0.3,"metacell") # default: 0.25
set_param("scm_mc_mark_min_gene_fold",2,"metacell") # default: 1.5
set_param("mc_plot_device",'pdf', "metacell")

set_param("mcell_mc2d_K",30,"metacell") # default: 20
set_param("mcell_mc2d_T_edge",0.02,"metacell") # default: 0.05
set_param("mcell_mc2d_max_confu_deg",8,"metacell") # default: 5
set_param("mcell_mc2d_edge_asym",FALSE,"metacell") # default: TRUE
set_param("mcell_mc2d_proj_blur",0.02,"metacell") # default: 0.02

dbpath<-"/home/asus/Desktop/invitro_20230308/"
mat_id<-"invitro"  #name of the matrix (umi-tab)

if(!dir.exists(dbpath)) dir.create(dbpath)
scdb_init(dbpath, force_reinit=T)

mcell_import_multi_mars(mat_id,dataset_table_fn = "/home/asus/Desktop/invitro_20230308/MARS_Batches.txt",
                        base_dir = "/home/asus/Desktop/invitro_20230308/umi_tab/",
                        patch_cell_name = T) #loading the amp_batches, patch_cell_name, means, add the name of batch onto the well_id


mat<-scdb_mat(mat_id) # transfer the umi-tab to variabe "mat"

all_mat_seurat<-CreateSeuratObject(mat@mat,min.cells = 40, min.features = 100)
all_mat_seurat[["percentage.mt"]]=PercentageFeatureSet(all_mat_seurat,pattern="MT-")

all_mat_seurat@meta.data$orig.ident = "human_proj" # to avoid the errors of cell id
Idents(all_mat_seurat) = "orig.ident"



VlnPlot(all_mat_seurat,features = c("nCount_RNA","nFeature_RNA", "percentage.mt"))


####
fgpath<-"/home/asus/Desktop/invitro_20230308/"

if(!dir.exists(fgpath)) dir.create(fgpath)
scfigs_init(fgpath)

mcell_plot_umis_per_cell(mat_id)

####


all_mat_seurat <- subset(all_mat_seurat, nCount_RNA > 200 & nCount_RNA < 10000 & percentage.mt< 50) # filter



setwd("/home/asus/Desktop/invitro_20230308/file3/")
library(Matrix)
write.table(data.frame(rownames(all_mat_seurat@assays$RNA@counts),rownames(all_mat_seurat@assays$RNA@counts)),file = 'genes.tsv',
            quote = F,sep = '\t',
            col.names = F,row.names = F)
write.table(colnames(all_mat_seurat@assays$RNA@counts),file = 'barcodes.tsv',quote = F,
            col.names = F,row.names = F)

counts <- as(all_mat_seurat@assays$RNA@counts, "sparseMatrix")
writeMM(counts, file = "matrix.mtx") 

####-----------------------------------------------------------

data<-Read10X("~/Desktop/invitro_20230308/file3/")


## construct Seurat object
seob_data <- CreateSeuratObject(counts = data)
seob_data <- NormalizeData(seob_data, normalization.method =  "LogNormalize")


seob_data <- FindVariableFeatures(seob_data, 
                                  selection.method = "vst", 
                                  nfeatures = 2000) 


seob_data <- ScaleData(seob_data, features = rownames(seob_data))

seob_data <- RunPCA(seob_data)
ElbowPlot(seob_data, ndims = 50)
## UMAP
seob_data <- RunUMAP(seob_data, dims = 1:30)


seob_data <- FindNeighbors(seob_data,
                           dims = 1:30)
seob_data <- FindClusters(seob_data,
                          resolution = 0.3, # 值越大，cluster 越多
                          random.seed = 1) 

umap.coor<-as.data.frame(seob_data@reductions$umap@cell.embeddings)
umap.coor$Well_ID<-rownames(umap.coor)

#meta<-read.csv("~/Desktop/invitro_20230308/wells_cells.csv",header=T)
meta<-read.csv("~/Desktop/invitro_20230308/wells_cells2.csv",header=T)
meta<-merge(umap.coor,meta.orig,by="Well_ID")


pram.score<-read.csv("~/Desktop/invitro_20230308/proangiogenic_score.csv")

#pram.score.matrix<-data[rownames(data)%in%pram.score$Proangiogenic.module,]
mat2<-mat@mat

pram.score.matrix.mat2<-mat2[rownames(mat2)%in%pram.score$Proangiogenic.module,]

pram.score.matrix.mat2<-as.matrix(pram.score.matrix.mat2)


col_mean_genes<-colMeans(as.matrix(pram.score.matrix.mat2[,1:ncol(pram.score.matrix.mat2)]))

col_mean_genes<-as.data.frame(col_mean_genes)
col_mean_genes$Well_ID<-rownames(col_mean_genes)

meta<-merge(meta,col_mean_genes,by="Well_ID")


###########---------------------------------pheatmap
## 2&4 YS
## 5&6 monocyte

meta.sub<-meta[meta$seurat_clusters%in%c("1","2","5","6"),]

data<-seob_data@assays$RNA@data  ## normalized data
data<-data[rownames(data)%in%c("S100A8","S100A9","LYZ", ## monocyte module
                               "MRC1","LYVE1","DAB2","CD163","RNASE1",## core macrophage module
                               "CD83","VEGFA","IL1B","PLAUR","CXCL8","BTG2","ICAM1","TNF","PLAUR"# proangiogenic module
),]
data<-data[,colnames(data)%in%meta.sub$Well_ID]

###
YS_vEC_0h.id<-subset(meta.sub, info == "vEC" & tissue == "Yolksac" & time == "0h")
YS_vEC_60h.id<-subset(meta.sub, info == "vEC" & tissue == "Yolksac" & time == "60h")
monocyte_vEC_0h.id<-subset(meta.sub, info == "vEC" & tissue == "monocyte" & time == "0h")
monocyte_vEC_60h.id<-subset(meta.sub, info == "vEC" & tissue == "monocyte" & time == "60h")

YS_uvEC_0h.id<-subset(meta.sub, info == "uvEC" & tissue == "Yolksac" & time == "0h")
YS_uvEC_60h.id<-subset(meta.sub, info == "uvEC" & tissue == "Yolksac" & time == "60h")
monocyte_uvEC_0h.id<-subset(meta.sub, info == "uvEC" & tissue == "monocyte" & time == "0h")
monocyte_uvEC_60h.id<-subset(meta.sub, info == "uvEC" & tissue == "monocyte" & time == "60h")



YS_vEC_0h<-data[,colnames(data)%in%YS_vEC_0h.id$Well_ID]
YS_vEC_60h<-data[,colnames(data)%in%YS_vEC_60h.id$Well_ID]
YS_uvEC_0h<-data[,colnames(data)%in%YS_uvEC_0h.id$Well_ID]
YS_uvEC_60h<-data[,colnames(data)%in%YS_uvEC_60h.id$Well_ID]


monocyte_vEC_0h<-data[,colnames(data)%in%monocyte_vEC_0h.id$Well_ID]
monocyte_vEC_60h<-data[,colnames(data)%in%monocyte_vEC_60h.id$Well_ID]
monocyte_uvEC_0h<-data[,colnames(data)%in%monocyte_uvEC_0h.id$Well_ID]
monocyte_uvEC_60h<-data[,colnames(data)%in%monocyte_uvEC_60h.id$Well_ID]


calc_mean<-function(tmp){
  
  tmp.mean<-rowMeans(as.matrix(tmp[,1:ncol(tmp)]))
  
  return(tmp.mean)
}

YS_vEC_0h_mean<-calc_mean(YS_vEC_0h)
YS_vEC_60h_mean<-calc_mean(YS_vEC_60h)
YS_uvEC_0h_mean<-calc_mean(YS_uvEC_0h)
YS_uvEC_60h_mean<-calc_mean(YS_uvEC_60h)

monocyte_vEC_0h_mean<-calc_mean(monocyte_vEC_0h)
monocyte_vEC_60h_mean<-calc_mean(monocyte_vEC_60h)
monocyte_uvEC_0h_mean<-calc_mean(monocyte_uvEC_0h)
monocyte_uvEC_60h_mean<-calc_mean(monocyte_uvEC_60h)

x<-rbind(YS_vEC_0h_mean,YS_vEC_60h_mean,YS_uvEC_0h_mean,YS_uvEC_60h_mean,monocyte_vEC_0h_mean,monocyte_vEC_60h_mean,monocyte_uvEC_0h_mean,monocyte_uvEC_60h_mean)
x<-t(x)

###

YS_vEC_0h_pram<-calc_pram(col_mean_genes, YS_vEC_0h)
YS_vEC_60h_pram<-calc_pram(col_mean_genes, YS_vEC_60h)
YS_uvEC_0h_pram<-calc_pram(col_mean_genes, YS_uvEC_0h)
YS_uvEC_60h_pram<-calc_pram(col_mean_genes, YS_uvEC_60h)

monocyte_vEC_0h_pram<-calc_pram(col_mean_genes, monocyte_vEC_0h)
monocyte_vEC_60h_pram<-calc_pram(col_mean_genes, monocyte_vEC_60h)
monocyte_uvEC_0h_pram<-calc_pram(col_mean_genes, monocyte_uvEC_0h)
monocyte_uvEC_60h_pram<-calc_pram(col_mean_genes, monocyte_uvEC_60h)

xx<-c(YS_vEC_0h_pram,YS_vEC_60h_pram,YS_uvEC_0h_pram,YS_uvEC_60h_pram,monocyte_vEC_0h_pram,monocyte_vEC_60h_pram,monocyte_uvEC_0h_pram,monocyte_uvEC_60h_pram)

#######

calc_pram<-function(col_mean_genes,tmp){
  
  
  tmp2<-col_mean_genes[col_mean_genes$Well_ID%in%colnames(tmp),]
  tmp2<-mean(tmp2$col_mean_genes)
  return(tmp2)
}

x<-rbind(xx,x)

rownames(x)[1]<-"Proangiogenic score"


bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))

annotation_col = data.frame(condition = c("vEC","vEC","uvEC","uvEC","vEC","vEC","uvEC","uvEC"),source=c("YS","YS","YS","YS","monocyte","monocyte","monocyte","monocyte"),time=c("0h","60h","0h","60h","0h","60h","0h","60h"))
rownames(annotation_col)<-colnames(x)

pheatmap(x, cluster_rows = F, cluster_cols = F,scale = "row",annotation_col = annotation_col,
         color = c(colorRampPalette(colors = c("#1B519C","white"))(length(bk)/2),colorRampPalette(colors = c("white","#A00021"))(length(bk)/2)),breaks = bk,cutree_rows = 3, show_colnames = F)
