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
#####
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

#####
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

#####
##pram vs pre-pram with new tab

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

write.csv(dat,"E:\\fig5/dat.csv")
#####
##volcano

dat = read.csv("E:\\fig5/compare_pre_PraM_vs_PraM.csv",row.names = 1)
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
  annotate(geom="text", x=-1, y=3, label="pre_PraM",  fontface="bold",colour='#2CA02CFF', size=5)+
  annotate(geom="text", x=2, y=3, label="PraM",  fontface="bold",colour='#FF7F0EFF', size=5)+
  theme_bw()+theme_classic()+xlab("log2(fold change)")+ylab("q value")+
  theme(legend.position = "none",
        axis.text = element_text(size=14))
ggsave("E:\\fig5/deg_volcano.png",p1,width = 5,height = 5,dpi = 300)

#####
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


#####
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

#####
##pie


i = "Male gonad"
tmp = mac.metatable[mac.metatable$tissue == i,]
tmp = tmp[tmp$time %in% after9,]

tmp2 = as.data.frame(table(tmp$subtype))
tmp2 = tmp2[order(tmp2$Freq,decreasing = T),]
tmp2 = tmp2[tmp2$Freq >100,]
tmp2$ratio = round(tmp2$Freq/sum(tmp2$Freq),2) 

pdf(paste0("E:\\fig5/pie_marker/",i,"_pie.pdf"),width = 5,height = 5)
pie(tmp2$ratio,labels = tmp2$ratio,
    col = c(mac.color[as.character(tmp2$Var1)],
            rep("grey88",2)),)
#legend("topright",legend = tmp2$Var1, cex=0.8, fill=mac.color[as.character(tmp2$Var1)])
dev.off()

for (i in setdiff(unique(mac.metatable$tissue),c("AGM","Embryo","Primitive_gut","Yolksac","Male gonad",
                                                 "Gastrointestinal tract","Thymus","Limb","Head"))) {
  tmp = mac.metatable[mac.metatable$tissue == i,]
  tmp = tmp[tmp$time %in% after9,]
  
  tmp2 = as.data.frame(table(tmp$subtype))
  tmp2 = tmp2[order(tmp2$Freq,decreasing = T),]
  tmp2$ratio = round(tmp2$Freq/sum(tmp2$Freq),2) 
  
  pdf(paste0("E:\\fig5/pie_marker/",i,"_pie.pdf"))
  pie(tmp2$ratio,labels = tmp2$ratio[1:3],
      col = c(mac.color[as.character(tmp2$Var1[1:3])],
              rep("grey88",(nrow(tmp2)-3))),)
  #legend("topright",legend = tmp2$Var1[1:3], cex=0.8, fill=mac.color[as.character(tmp2$Var1[1:3])])
  dev.off()
  
}

###BAR CHART

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
figre5_mac_sub_ratio_per_tissue = function(the_metatable,the_dir,thew , theh){
  tmp = the_metatable
  
  tmp= tmp[-which(tmp$time == "Adult"),]
  
  tmp2 = tmp
  
  for(i in unique(tmp2$tissue)){
    tmp = tmp2[,c("tissue","time","subtype")]
    tmp = tmp[tmp$tissue == i,]
    tmp = tmp[,c("time","subtype")]
    tmp$new = paste(tmp$time,tmp$subtype,sep = "_")
    tmp3 = as.data.frame(table(tmp$new))
    #tmp3 = tmp3[tmp3$Freq >= 20,]
    if(nrow(tmp3)>0){
      tmp = tmp[tmp$new %in% tmp3$Var1,]
      tmp = tmp[,c("time","subtype")]
      tmp$count = 0
      for (a in unique(tmp$time)) {
        
        tmp$count[tmp$time == a ] = nrow(tmp[tmp$time == a,])
        
      }
      other = setdiff(time_order,c(unique(tmp$time),"Adult"))
      other = as.data.frame(other)
      other$blank = "blank"
      other$count = 0
      colnames(other) = colnames(tmp)
      tmp = rbind(tmp,other)
      tmp$time <- factor(tmp$time,levels=time_order,ordered = TRUE)
      p2 = ggplot(tmp)+
        geom_bar(aes(x= time,fill = subtype),stat="count",position = "fill",width = 0.9,color = "black")+
        geom_text(data = tmp[!duplicated(tmp),],aes(x = time,y = 0.95,label = count))+
        xlab("")+ylab("")+
        ggtitle(paste(i," total cell number ",nrow(tmp),sep = ""))+
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
      ggsave(paste(the_dir,i,"_tissue.pdf",sep = ""),p2,width = thew,height = theh,dpi = 300)
    }
  }
}
figre5_mac_sub_ratio_per_tissue(the_metatable = mac12.metatable,
                                the_dir = "E:\\fig5/pertissue3/",thew = 8,theh = 5)

#####
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

#####

cs_stage = c("cs11","cs12","cs13","cs14","cs18","cs19","cs21","cs23")
`week_9-13` = c("week9","week10","week11","week12","week13")
after_week16 = c("week16","week19","week20","week23","week27")

mac12.metatable.stat = mac12.metatable

mac12.metatable.stat$period <- mac12.metatable.stat$time
mac12.metatable.stat$period[mac12.metatable.stat$period%in%cs_stage] <- "cs_stage"
mac12.metatable.stat$period[mac12.metatable.stat$period%in%`week_9-13`] <- "week_9-13"
mac12.metatable.stat$period[mac12.metatable.stat$period%in%after_week16] <- "after_week16"

use.embryo = as.data.frame(table(mac12.metatable.stat$embryo))
use.embryo = use.embryo$Var1[use.embryo$Freq >10]
mac12.metatable.stat = mac12.metatable.stat[mac12.metatable.stat$embryo  %in% use.embryo,]
mac12.metatable.stat = mac12.metatable.stat[mac12.metatable.stat$embryo != "embryo 6",]
mac12.metatable.stat = mac12.metatable.stat[mac12.metatable.stat$embryo != "embryo 41",]
mac12.metatable.stat = mac12.metatable.stat[mac12.metatable.stat$embryo != "embryo 4",]

### filter out adult cells
meta_fetal <- mac12.metatable.stat %>% filter(time!="Adult")
table(meta_fetal$period)
### compute ratio of subtypes in different time periods
ratio <- list()
for(i in unique(meta_fetal$embryo)){
  print(i)
  df <- meta_fetal %>% filter(embryo==i)
  ratio[[i]] <- tapply(df$subtype,df$period,function(x){prop.table(table(x))}) %>% do.call(cbind,.) %>% as.data.frame()
  ratio[[i]] <- ratio[[i]] %>% mutate(subtype=rownames(.),period=colnames(.),embryo=i)
  colnames(ratio[[i]])[1] <- "ratio"
}
ratio <- do.call(rbind,ratio)
ratio$period <- factor(ratio$period,levels = c("cs_stage","week_9-13","after_week16"))
### draw mean horizontal line for each group; control segment length
mean <- tapply(ratio$ratio,list(ratio$subtype,ratio$period),mean) %>% melt()
colnames(mean) <- c("subtype","period","mean")
mean$x_start <- as.numeric(mean$period)-0.1
mean$x_end <- as.numeric(mean$period)+0.1
###
plt_Frac <- function(ratio, path){
  for(i in unique(ratio$subtype)){
    color <- mac.color[i]
    if(is.na(color)){color <- "grey"}
    ratio_sub <- ratio %>% filter(subtype==i)
    mean <- tapply(ratio_sub$ratio,list(ratio_sub$subtype,ratio_sub$period),mean) %>% melt()
    colnames(mean) <- c("subtype","period","mean")
    mean$x_start <- as.numeric(mean$period)-0.1
    mean$x_end <- as.numeric(mean$period)+0.1
    
    pdf(paste(path,i,"_Fraction_change",".pdf",sep=""),width =5,height = 8)
    print(ggplot(ratio_sub,aes(x = period, y = ratio))+
            geom_jitter(shape=21,width=0.2,colour="black",size=7, fill=color)+
            geom_segment(data=mean,aes(x=x_start, y=mean, xend=x_end, yend=mean))+
            xlab("")+ylab("")+
            theme_bw()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_text(size=30,angle = 45,vjust = 0.7,hjust=0.8,colour = "black"),
                  axis.text.y = element_text(size=30,colour = "black"),
                  title = element_text(size=30))+
            ggtitle(paste0(i,"  ***"))+
            stat_compare_means(aes(label = ..p.signif.., group=period),vjust=-2,size=6))
    dev.off()
  }
}

plt_Frac(ratio, "E:\\fig5/stat_ratio/")


meta_mf <- metatable %>% filter(time!="Adult") %>% filter(major=="macrophage")
meta_mf <- meta_mf %>% filter(!embryo%in%c("embryo 54","embryo 12")) # contain 2/44 cell respectively

### divide all time into 3 periods 119995 121152
cs_stage = c("cs11","cs12","cs13","cs14","cs18","cs19","cs21","cs23")
`week_9-13` = c("week9","week10","week11","week12","week13")
after_week16 = c("week16","week19","week20","week23","week27")
meta_mf$period <- meta_mf$time
meta_mf$period[meta_mf$period%in%cs_stage] <- "CS stage"
meta_mf$period[meta_mf$period%in%`week_9-13`] <- "9-13 PCW"
meta_mf$period[meta_mf$period%in%after_week16] <- "After 16 PCW"
table(meta_mf$period)

### compute ratio of subtypes in different time periods
ratio <- list()
for(i in unique(meta_mf$embryo)){
  df <- meta_mf %>% filter(embryo==i)
  ratio[[i]] <- tapply(df$subtype,df$period,function(x){prop.table(table(x))}) %>% do.call(cbind,.) %>% as.data.frame()
  ratio[[i]] <- ratio[[i]] %>% mutate(subtype=rownames(.),period=colnames(.),embryo=i)
  colnames(ratio[[i]])[1] <- "ratio"
}

### zero not present in tables, complete
ratio_complete <- function(ratio){
  for(i in names(ratio)){
    subtype_full <- unique(meta_mf$subtype)
    subtype_complete <- setdiff(subtype_full,unique(ratio[[i]]$subtype))
    if(!is.null(subtype_complete)){
      df <- data.frame(subtype=subtype_complete)
      df$ratio <- rep(0,nrow(df))
      df$embryo <- rep(i,nrow(df))
      df$period <- rep(unique(ratio[[i]]$period),nrow(df))
      df <- df[,c("ratio","subtype","period","embryo")]
      ratio[[i]] <- rbind(ratio[[i]],df)
    }
  }
  
  return(ratio)
}
ratio <- ratio_complete(ratio)

ratio <- do.call(rbind,ratio)
ratio$period <- factor(ratio$period,levels = c("CS stage","9-13 PCW","After 16 PCW"))

color_grp <- c("#A6A8A9FF","#484A4BFF","#141414FF")
names(color_grp) <- sort(unique(ratio$period))

my_comp <- list(c("CS stage","9-13 PCW"),c("9-13 PCW","After 16 PCW"),c("CS stage","After 16 PCW"))


#####
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

#####
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

write.csv(final_marker,"E:\\fig5/degofspecific.csv")


#final_marker = final_marker[-which(final_marker$gene %in% c("SOX11","TTN","NR4A1","TNNT2","JUN","ANXA1","HMGA2","NR4A2","JUNB","NR4A4",
#                                                            "HLA-DRB1","HLA-DRB5","HLA-DMA",mhc2)),]
#final_marker$sample[final_marker$sample == "langerhas"] = "langerhans"
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

#####
use.gene = c("IL1B","CXCL8","TNF","DAB2","MRC1","VEGFA","HMGA2")

tmp = metatable[metatable$subtype %in% c("PraM","pre_PraM"),]
#tmp = tmp[,c("Well_ID","subtype")]
tmp_tab = as.data.frame(t(as.matrix(a_reduced[use.gene,])))
tmp_tab$Well_ID = rownames(tmp_tab)

tmp_tab = merge(tmp,tmp_tab,by = "Well_ID",all.x = T)

tmp_tab = tmp_tab[tmp_tab$time != "Adult",]
tmp_tab = tmp_tab[,c("subtype",use.gene)]
tmp_tab = reshape2::melt(tmp_tab)

ggplot(tmp_tab)+
  geom_violin(aes(x = subtype, y = log2(value+1),fill = variable),scale = "width",color = "grey30")+
  facet_grid(variable ~ .,scales = "free")+
  xlab("")+ylab("")+theme_bw()+
  scale_x_discrete("")+
  scale_y_continuous(expand = c(0, 0), labels = function(x)
    c(rep(x = '', times = length(x) - 2), x[length(x) - 1], '')) +
  theme(strip.background = element_rect(color = "black", fill = "white"),
        panel.spacing = unit(x = 0, units = 'lines'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=12,angle = 45,vjust = 0.7,hjust=0.8,colour = "black"),
        strip.text.y = element_text(size = 12,angle = 360, face = 'bold'),
        legend.position = "none")
