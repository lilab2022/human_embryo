###vs adult

library(ggsignif)

#####
##amp10 metatable & umi.tab

adult.amp = list.files("E:\\adult/meta2/")
adult.umi = read.table(paste0("E:\\adult/meta2/",adult.amp[1]))

for (i in adult.amp[2:10]){
  tmp = read.table(paste0("E:\\adult/meta2/",i))
  adult.umi = cbind(adult.umi,tmp)
}

tmp = colnames(adult.umi)
tmp2 = c()

for (i in tmp){
  s = unlist(strsplit(i,"_"))[2]
  tmp2 = c(tmp2,s)
}

adult.metatable = as.data.frame(colnames(adult.umi))
colnames(adult.metatable) = "Well_ID"
adult.metatable$Amp_batch_ID = tmp2
adult.metatable$time = "Adult"
adult.metatable$tissue = "Liver"
adult.metatable$tissue[adult.metatable$Amp_batch_ID %in% c("AB2789","AB2805","AB3069")] = "Spleen"
adult.metatable$tissue[adult.metatable$Amp_batch_ID %in% c("AB3112")] = "Brain"
adult.metatable$tissue[adult.metatable$Amp_batch_ID %in% c("AB3145")] = "Intestine"

adult.mcmc = read.table("E:\\adult/mcmc.txt")
colnames(adult.mcmc) = "metacell"
adult.mcmc$metacell = paste0("X",adult.mcmc$metacell)
adult.mcmc$Well_ID = rownames(adult.mcmc)
adult.metatable = merge(adult.metatable,adult.mcmc,by = "Well_ID")
adult.metatable$subtype = "discard"
adult.metatable$subtype[adult.metatable$metacell == "X11"] = "Kupffer_cell"
adult.metatable$subtype[adult.metatable$metacell %in% c("X9","X10")] = "PraM"
adult.metatable = adult.metatable[adult.metatable$subtype != "discard",]

the_mat = a[setdiff(rownames(a),str_subset(rownames(a), pattern ='^MT-')),metatable$Well_ID]
the_mat = CreateSeuratObject(counts = the_mat)
the_mat <- NormalizeData(the_mat, normalization.method =  "RC")
a_reduced = the_mat@assays$RNA@data
save(a_reduced,file = "E:\\singlecellaR/all_cell_noemalized_tab.Rdata")

the_mat.adult = adult.umi[setdiff(rownames(a),str_subset(rownames(a), pattern ='^MT-')),adult.metatable$Well_ID]
the_mat.adult = CreateSeuratObject(counts = the_mat.adult)

the_mat.combine = merge(x = the_mat, y = the_mat.adult)
the_mat.combine <- NormalizeData(the_mat.combine, normalization.method =  "RC")
the_mat.combine = the_mat.combine@assays$RNA@data

adult.metatable = adult.metatable[,c("Well_ID","Amp_batch_ID","time","tissue","subtype")]
tmp = metatable[,c("Well_ID","Amp_batch_ID","time","tissue","subtype")]

adult.metatable = rbind(adult.metatable,tmp)
#####
## brain microglia
gene_list <- read.csv("D:\\all_v2/vs/human.signature.genes.v1.gmt",sep="\t",header = F)

use.gene = module$gene
save(use.gene,file = "E:\\adult/use_gene.Rdata")

table1 = adult.metatable[adult.metatable$subtype == "microglia" & adult.metatable$tissue == "Brain",]
table1 = table1[table1$time %in% c("Adult","week27"),]
brain = the_mat.combine[intersect(use.gene,rownames(the_mat.combine)),table1$Well_ID]
brain = as.matrix(brain)
brain = as.data.frame(brain)

for (i in unique(module$category)){
  brain[i,] = colMeans(brain[module$gene[module$category == i],])
}

brain["g2m",] = colMeans(brain[intersect(rownames(brain),as.character(gene_list[1,])),])
brain = brain[c("g2m",unique(module$category)),]

time = adult.metatable[,c("time","Well_ID")]

brain = as.data.frame(t(brain))
brain$Well_ID = rownames(brain)
brain = melt(brain)
brain = merge(brain,time,by = "Well_ID",all.x = T)
brain$time[brain$time != "Adult"] = "Embryo"
tmp = brain

for(i in c(unique(module$category))){
brain = tmp
brain = brain[brain$variable == i,]

brain$time[brain$time == "Embryo"] = "26PCW"
brain$time[brain$time == "Adult"] = "Adult stage"
brain$time = factor(brain$time, levels=c('26PCW','Adult stage'))
compaired <- list(c('26PCW','Adult stage'))

brain$time = factor(brain$time, levels=c('26PCW','Adult stage'))

fc = mean(brain$value[brain$time == "Adult stage"])/mean(brain$value[brain$time == "26PCW"]+1)
p1 = ggplot(brain, aes(x= time, y = log2(value+1),fill = time))+
  geom_boxplot(outlier.color = NA)+
  stat_boxplot(geom ="errorbar",width=0.15)+
  geom_signif(comparisons = compaired,step_increase = 0.1,test =wilcox.test)+
  xlab("")+ylab("")+ylim(0,8)+
  ggtitle(paste0("fc ",round(fc,2)))+
  theme_classic()+
  scale_fill_manual(values = c("#546de5", "#ff4757"))+
  theme(legend.position = 'none')
ggsave(paste0("E:\\adult/brain_microglia/",i,".pdf"),p1,width = 2,height = 3.5,dpi = 300)
}

table1$embryo = table1$time
table1$major_type = "macrophage"
table1$tissue = table1$time

pop <- 'macrophage'
tissueA <- "week27"; cutA <- 15
tissueB <- "Adult"; cutB <- 15

n <- 50
ntime <- 100

get_n_cells_per_pop_per_disease(table1, n, pop, tissueA, cutA, diag =F)
get_n_cells_per_pop_per_disease(table1, n, pop, tissueB, cutB, diag =F)

ind_diff <- diff_calc(table1, the_mat.combine[,table1$Well_ID], pop, n, ntime, tissueA, tissueB, cutA, cutB)
ind_diff_stat <- diff_stat(ind_diff, ntime)

f_cut <- 1.5
q_cut <- 0.05
pct_cut <- 0.05
ind_diff_stat_good <- ind_diff_stat[ind_diff_stat$pctA.mean > pct_cut | ind_diff_stat$pctB.mean > pct_cut,]
f <- log2((ind_diff_stat_good$avgB.mean + 0.01)/(ind_diff_stat_good$avgA.mean + 0.01))
q <- log(-(log(ind_diff_stat_good$wilcox.p.mean, 10)) + 1)
dat <- data.frame(g = rownames(ind_diff_stat_good), f = f, q = q)
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

dat$label = ""
dat$module = ""
for (i in module$gene[module$category == "MHC-II"]) {
  dat$label[dat$g == i] = dat$g[dat$g == i]
  dat$module[dat$g == i] = "MHC-II"
}
for (i in module$gene[module$category == "microglia"]) {
  dat$label[dat$g == i] = dat$g[dat$g == i]
  dat$module[dat$g == i] = "microglia"
}


dat = dat[discard_gene2(rownames(dat)),]
dat$label = ""
for (i  in dat$g[dat$group != "non" & dat$module != ""]) {
 dat$label[dat$g == i] = i 
}

p1 = ggplot(dat, aes(f, q,fill = module,label = label)) +
  geom_point(shape = 21,size=2,color = "grey66")+
  geom_point(data = dat[dat$module != "",],shape = 21,size=4,color = "grey66")+
  geom_text_repel(data = dat[dat$label != "",])+
  scale_fill_manual(values=c( "#1F77B4FF","#E64B35FF","grey88"),breaks = c("microglia","MHC-II",""))+
  geom_vline(xintercept = log2(2) , linetype="dashed", color="grey77") + 
  geom_vline(xintercept = -log2(2), linetype="dashed", color="grey77") + 
  geom_hline(yintercept = log(-log10(q_cut) + 1), linetype="dashed", color="grey30") +
  ylab('log((-log10 q) + 1)') + xlab('log2 f')+
  theme_bw()+theme_classic()+
  theme(legend.position = c(0.9,0.15),
        axis.text = element_text(size=14))
ggsave("E:\\adult/microglia_embryo_adult.png",p1,width = 5,height = 5,dpi = 300)
write.csv(dat,"E:\\adult/microglia_embryo_adult.csv")
write.csv(ind_diff_stat,"E:\\adult/kupffer_embryo_adult_ind_diff_stat.csv")





#####
table1 = adult.metatable[adult.metatable$subtype == "PraM" & adult.metatable$tissue == "Skin",]
table1 = table1[table1$time %in% c("Adult","week23","week20"),]

skinpram = the_mat.combine[intersect(use.gene,rownames(the_mat.combine)),table1$Well_ID]
skinpram = as.matrix(skinpram)
skinpram = as.data.frame(skinpram)
       
for (i in unique(module$category)){
  skinpram[i,] = colMeans(skinpram[module$gene[module$category == i],])
}
       
skinpram["g2m",] = colMeans(skinpram[intersect(rownames(skinpram),as.character(gene_list[1,])),])
skinpram = skinpram[c("g2m",unique(module$category)),]
       

skinpram = as.data.frame(t(skinpram))
skinpram$Well_ID = rownames(skinpram)
skinpram = melt(skinpram)
skinpram = merge(skinpram,time,by = "Well_ID",all.x = T)
skinpram$time[skinpram$time != "Adult"] = "Embryo"
tmp = skinpram
       
for(i in c(unique(module$category))){
  skinpram = tmp
  skinpram = skinpram[skinpram$variable == i,]
         
  skinpram$time[skinpram$time == "Embryo"] = "PCW20-23"
  skinpram$time[skinpram$time == "Adult"] = "Adult stage"
  skinpram$time = factor(skinpram$time, levels=c('PCW20-23','Adult stage'))
  compaired <- list(c('PCW20-23','Adult stage'))
         
  skinpram$time = factor(skinpram$time, levels=c('PCW20-23','Adult stage'))
  fc = mean(skinpram$value[skinpram$time == "Adult stage"])/mean(skinpram$value[skinpram$time == "PCW20-23"])
  p1 = ggplot(skinpram, aes(x= time, y = log2(value+1),fill = time))+
           geom_boxplot(outlier.color = NA)+
           stat_boxplot(geom ="errorbar",width=0.15)+
           geom_signif(comparisons = compaired,step_increase = 0.1,test =wilcox.test)+
           xlab("")+ylab("")+ylim(0,8)+
           ggtitle(paste0("fc ",round(fc,2)))+
           theme_classic()+
           scale_fill_manual(values = c("#546de5", "#ff4757"))+
           theme(legend.position = 'none')
         ggsave(paste0("E:\\adult/skin_mac1/",i,".pdf"),p1,width = 2,height = 3.5,dpi = 300)
}

table1$time[table1$time != "Adult"] = "PCW20-23"
table1$embryo = table1$time
table1$major_type = "macrophage"
table1$tissue = table1$time

pop <- 'macrophage'
tissueA <- "PCW20-23"; cutA <- 15
tissueB <- "Adult"; cutB <- 15

n <- 50
ntime <- 100

get_n_cells_per_pop_per_disease(table1, n, pop, tissueA, cutA, diag =F)
get_n_cells_per_pop_per_disease(table1, n, pop, tissueB, cutB, diag =F)

ind_diff <- diff_calc(table1, the_mat.combine[,table1$Well_ID], pop, n, ntime, tissueA, tissueB, cutA, cutB)
ind_diff_stat <- diff_stat(ind_diff, ntime)

f_cut <- 1.5
q_cut <- 0.05
pct_cut <- 0.05
ind_diff_stat_good <- ind_diff_stat[ind_diff_stat$pctA.mean > pct_cut | ind_diff_stat$pctB.mean > pct_cut,]
f <- log2((ind_diff_stat_good$avgB.mean + 0.01)/(ind_diff_stat_good$avgA.mean + 0.01))
q <- log(-(log(ind_diff_stat_good$wilcox.p.mean, 10)) + 1)
dat <- data.frame(g = rownames(ind_diff_stat_good), f = f, q = q)
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

dat$label = ""
dat$module = ""
for (i in module$gene[module$category == "MHC-II"]) {
  dat$label[dat$g == i] = dat$g[dat$g == i]
  dat$module[dat$g == i] = "MHC-II"
}
for (i in module$gene[module$category == "PraM"]) {
  dat$label[dat$g == i] = dat$g[dat$g == i]
  dat$module[dat$g == i] = "PraM"
}


dat = dat[discard_gene2(rownames(dat)),]
dat$label = ""
for (i  in dat$g[dat$group != "non" & dat$module != ""]) {
  dat$label[dat$g == i] = i 
}
dat = dat[dat$g != "SPP1",]
p1 = ggplot(dat, aes(f, q,fill = module,label = label)) +
  geom_point(shape = 21,size=2,color = "grey66")+
  geom_point(data = dat[dat$module != "",],shape = 21,size=4,color = "grey66")+
  geom_text_repel(data = dat[dat$label != "",])+
  scale_fill_manual(values=c( "#FF7F0EFF","#E64B35FF","grey88"),breaks = c("PraM","MHC-II",""))+
  geom_vline(xintercept = log2(2) , linetype="dashed", color="grey77") + 
  geom_vline(xintercept = -log2(2), linetype="dashed", color="grey77") + 
  geom_hline(yintercept = log(-log10(q_cut) + 1), linetype="dashed", color="grey30") +
  ylab('log((-log10 q) + 1)') + xlab('log2 f')+
  theme_bw()+theme_classic()+
  theme(legend.position = c(0.9,0.15),
        axis.text = element_text(size=14))
ggsave("E:\\adult/skin_pram_embryo_adult2.pdf",p1,width = 5,height = 5,dpi = 300)
write.csv(dat,"E:\\adult/skin_pram_embryo_adult.csv")
#####
table1 = adult.metatable[adult.metatable$subtype == "langerhans" & adult.metatable$tissue == "Skin",]
table1 = table1[table1$time %in% c("Adult","week11"),]

skinlangerhans = the_mat.combine[intersect(use.gene,rownames(the_mat.combine)),table1$Well_ID]
skinlangerhans = as.matrix(skinlangerhans)
skinlangerhans = as.data.frame(skinlangerhans)

for (i in unique(module$category)){
  skinlangerhans[i,] = colMeans(skinlangerhans[module$gene[module$category == i],])
}

skinlangerhans["g2m",] = colMeans(skinlangerhans[intersect(rownames(skinlangerhans),as.character(gene_list[1,])),])
skinlangerhans = skinlangerhans[c("g2m",unique(module$category)),]


skinlangerhans = as.data.frame(t(skinlangerhans))
skinlangerhans$Well_ID = rownames(skinlangerhans)
skinlangerhans = melt(skinlangerhans)
skinlangerhans = merge(skinlangerhans,time,by = "Well_ID",all.x = T)
skinlangerhans$time[skinlangerhans$time != "Adult"] = "Embryo"
tmp = skinlangerhans

for(i in c(unique(module$category))){
  skinlangerhans = tmp
  skinlangerhans = skinlangerhans[skinlangerhans$variable == i,]
  
  skinlangerhans$time[skinlangerhans$time == "Embryo"] = "PCW12"
  skinlangerhans$time[skinlangerhans$time == "Adult"] = "Adult stage"
  skinlangerhans$time = factor(skinlangerhans$time, levels=c('PCW12','Adult stage'))
  compaired <- list(c('PCW12','Adult stage'))
  
  skinlangerhans$time = factor(skinlangerhans$time, levels=c('PCW12','Adult stage'))
  fc = mean(skinlangerhans$value[skinlangerhans$time == "Adult stage"])/mean(skinlangerhans$value[skinlangerhans$time == "PCW12"])
  
  p1 = ggplot(skinlangerhans, aes(x= time, y = log2(value+1),fill = time))+
    geom_boxplot(outlier.color = NA)+
    stat_boxplot(geom ="errorbar",width=0.15)+
    geom_signif(comparisons = compaired,step_increase = 0.1,test =wilcox.test)+
    xlab("")+ylab("")+ylim(0,8)+
    ggtitle(paste0("fc ",round(2^abs(fc),2)))+
    theme_classic()+
    scale_fill_manual(values = c("#546de5", "#ff4757"))+
    theme(legend.position = 'none')
  ggsave(paste0("E:\\adult/skinlangerhans/",i,".png"),p1,width = 2,height = 3.5,dpi = 300)
}


table1$embryo = table1$time
table1$major_type = "macrophage"
table1$tissue = table1$time

pop <- 'macrophage'
tissueA <- "week11"; cutA <- 15
tissueB <- "Adult"; cutB <- 15

n <- 20
ntime <- 50

get_n_cells_per_pop_per_disease(table1, n, pop, tissueA, cutA, diag =F)
get_n_cells_per_pop_per_disease(table1, n, pop, tissueB, cutB, diag =F)

ind_diff <- diff_calc(table1, the_mat.combine[,table1$Well_ID], pop, n, ntime, tissueA, tissueB, cutA, cutB)
ind_diff_stat <- diff_stat(ind_diff, ntime)

f_cut <- 1.5
q_cut <- 0.05
pct_cut <- 0.05
ind_diff_stat_good <- ind_diff_stat[ind_diff_stat$pctA.mean > pct_cut | ind_diff_stat$pctB.mean > pct_cut,]
f <- log2((ind_diff_stat_good$avgB.mean + 0.01)/(ind_diff_stat_good$avgA.mean + 0.01))
q <- log(-(log(ind_diff_stat_good$wilcox.p.mean, 10)) + 1)
dat <- data.frame(g = rownames(ind_diff_stat_good), f = f, q = q)
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

dat$label = ""
dat$module = ""
for (i in module$gene[module$category == "MHC-II"]) {
  dat$label[dat$g == i] = dat$g[dat$g == i]
  dat$module[dat$g == i] = "MHC-II"
}
for (i in module$gene[module$category == "langerhans"]) {
  dat$label[dat$g == i] = dat$g[dat$g == i]
  dat$module[dat$g == i] = "langerhans"
}


dat = dat[discard_gene2(rownames(dat)),]
dat$label = ""
for (i  in dat$g[dat$group != "non" & dat$module != ""]) {
  dat$label[dat$g == i] = i 
}
dat = dat[dat$g != "SPP1",]
p1 = ggplot(dat, aes(f, q,fill = module,label = label)) +
  geom_point(shape = 21,size=2,color = "grey66")+
  geom_point(data = dat[dat$module != "",],shape = 21,size=4,color = "grey66")+
  geom_text_repel(data = dat[dat$label != "",])+
  scale_fill_manual(values=c( "#DBDB8DFF","#E64B35FF","grey88"),breaks = c("langerhans","MHC-II",""))+
  geom_vline(xintercept = log2(2) , linetype="dashed", color="grey77") + 
  geom_vline(xintercept = -log2(2), linetype="dashed", color="grey77") + 
  geom_hline(yintercept = log(-log10(q_cut) + 1), linetype="dashed", color="grey30") +
  ylab('log((-log10 q) + 1)') + xlab('log2 f')+
  theme_bw()+theme_classic()+
  theme(legend.position = c(0.1,0.15),
        axis.text = element_text(size=14))
ggsave("E:\\adult/langerhans_embryo_adult2.png",p1,width = 5,height = 5,dpi = 300)
write.csv(dat,"E:\\adult/skin_langerhans_embryo_adult.csv")
write.csv(ind_diff_stat,"E:\\adult/skin_langerhans_embryo_vsadult.csv")


#####
table1 = adult.metatable[adult.metatable$subtype == "Kupffer_cell" & adult.metatable$tissue == "Liver",]
table1 = table1[table1$time %in% c("Adult","week12"),]

liverkupffer = the_mat.combine[intersect(use.gene,rownames(the_mat.combine)),table1$Well_ID]
liverkupffer = as.matrix(liverkupffer)
liverkupffer = as.data.frame(liverkupffer)

for (i in unique(module$category)){
  liverkupffer[i,] = colMeans(liverkupffer[module$gene[module$category == i],])
}

liverkupffer["g2m",] = colMeans(liverkupffer[intersect(rownames(liverkupffer),as.character(gene_list[1,])),])
liverkupffer = liverkupffer[c("g2m",unique(module$category)),]


liverkupffer = as.data.frame(t(liverkupffer))
liverkupffer$Well_ID = rownames(liverkupffer)
liverkupffer = melt(liverkupffer)
liverkupffer = merge(liverkupffer,time,by = "Well_ID",all.x = T)
liverkupffer$time[liverkupffer$time != "Adult"] = "Embryo"
tmp = liverkupffer

for(i in c("g2m",unique(module$category))){
  liverkupffer = tmp
  liverkupffer = liverkupffer[liverkupffer$variable == i,]
  
  liverkupffer$time[liverkupffer$time == "Embryo"] = "PCW12"
  liverkupffer$time[liverkupffer$time == "Adult"] = "Adult stage"
  liverkupffer$time = factor(liverkupffer$time, levels=c('PCW12','Adult stage'))
  compaired <- list(c('PCW12','Adult stage'))
  
  liverkupffer$time = factor(liverkupffer$time, levels=c('PCW12','Adult stage'))
  fc = -mean(log2(liverkupffer$value[liverkupffer$time == "Adult stage"]+1))+mean(log2(liverkupffer$value[liverkupffer$time == "PCW12"]+1))
  p1 = ggplot(liverkupffer, aes(x= time, y = log2(value+1),fill = time))+
    geom_boxplot(outlier.color = NA)+
    stat_boxplot(geom ="errorbar",width=0.15)+
    geom_signif(comparisons = compaired,step_increase = 0.1,test =wilcox.test)+
    xlab("")+ylab("")+ylim(0,8)+
    ggtitle(paste0("log2 fc ",round(fc,2)))+
    theme_classic()+
    scale_fill_manual(values = c("#546de5", "#ff4757"))+
    theme(legend.position = 'none')
  ggsave(paste0("E:\\adult/liverkupffer/",i,".png"),p1,width = 2,height = 3.5,dpi = 300)
}

table1$embryo = table1$time
table1$major_type = "macrophage"
table1$tissue = table1$time

pop <- 'macrophage'
tissueA <- "week12"; cutA <- 15
tissueB <- "Adult"; cutB <- 15

n <- 50
ntime <- 70

get_n_cells_per_pop_per_disease(table1, n, pop, tissueA, cutA, diag =F)
get_n_cells_per_pop_per_disease(table1, n, pop, tissueB, cutB, diag =F)

ind_diff <- diff_calc(table1, the_mat.combine[,table1$Well_ID], pop, n, ntime, tissueA, tissueB, cutA, cutB)
ind_diff_stat <- diff_stat(ind_diff, ntime)

f_cut <- 1.5
q_cut <- 0.05
pct_cut <- 0.05
ind_diff_stat_good <- ind_diff_stat[ind_diff_stat$pctA.mean > pct_cut | ind_diff_stat$pctB.mean > pct_cut,]
f <- log2((ind_diff_stat_good$avgB.mean + 0.01)/(ind_diff_stat_good$avgA.mean + 0.01))
q <- log(-(log(ind_diff_stat_good$wilcox.p.mean, 10)) + 1)
dat <- data.frame(g = rownames(ind_diff_stat_good), f = f, q = q)
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

dat$label = ""
dat$module = ""
for (i in module$gene[module$category == "MHC-II"]) {
  dat$label[dat$g == i] = dat$g[dat$g == i]
  dat$module[dat$g == i] = "MHC-II"
}
for (i in module$gene[module$category == "kupffer"]) {
  dat$label[dat$g == i] = dat$g[dat$g == i]
  dat$module[dat$g == i] = "kupffer"
}


dat = dat[discard_gene2(rownames(dat)),]

p1 = ggplot(dat, aes(f, q,fill = module)) +
  geom_point(shape = 21,size=2,color = "grey66")+
  geom_point(data = dat[dat$module != "",],shape = 21,size=4,color = "grey66")+
  scale_fill_manual(values=c( "#FF9896FF","#E64B35FF","grey88"),breaks = c("kupffer","MHC-II",""))+
  geom_vline(xintercept = log2(f_cut) , linetype="dashed", color="grey30") + 
  geom_vline(xintercept = -log2(f_cut), linetype="dashed", color="grey30") + 
  geom_hline(yintercept = log(-log10(q_cut) + 1), linetype="dashed", color="grey30") +
  ylab('log((-log10 q) + 1)') + xlab('log2 f')+
  theme_bw()+theme_classic()+
  theme(legend.position = c(0.2,0.2),
        axis.text = element_text(size=14))
ggsave("E:\\adult/kupffer_embryo_adult.png",p1,width = 5,height = 5,dpi = 300)
write.csv(dat,"E:\\adult/kupffer_embryo_adult.csv")
write.csv(ind_diff_stat,"E:\\adult/kupffer_embryo_adult_ind_diff_stat.csv")


#####
##
table1 = adult.metatable[adult.metatable$subtype == "microglia" & adult.metatable$tissue == "Brain",]
table1 = table1[table1$time %in% c("Adult","week27"),]
brain = the_mat.combine[intersect(use.gene,rownames(the_mat.combine)),table1$Well_ID]
brain = as.matrix(brain)
brain = as.data.frame(brain)

##
table1 = adult.metatable[adult.metatable$subtype == "PraM" & adult.metatable$tissue == "Skin",]
table1 = table1[table1$time %in% c("Adult","week23","week20"),]
skinpram = the_mat.combine[intersect(use.gene,rownames(the_mat.combine)),table1$Well_ID]
skinpram = as.matrix(skinpram)
skinpram = as.data.frame(skinpram)

##
table1 = adult.metatable[adult.metatable$subtype == "langerhans" & adult.metatable$tissue == "Skin",]
table1 = table1[table1$time %in% c("Adult","week12"),]

skinlangerhans = the_mat.combine[intersect(use.gene,rownames(the_mat.combine)),table1$Well_ID]
skinlangerhans = as.matrix(skinlangerhans)
skinlangerhans = as.data.frame(skinlangerhans)

##
table1 = adult.metatable[adult.metatable$subtype == "Kupffer_cell" & adult.metatable$tissue == "Liver",]
table1 = table1[table1$time %in% c("Adult","week12"),]

liverkupffer = the_mat.combine[intersect(use.gene,rownames(the_mat.combine)),table1$Well_ID]
liverkupffer = as.matrix(liverkupffer)
liverkupffer = as.data.frame(liverkupffer)

###
all = cbind(brain[module$gene[module$category == "MHC-II"],],skinpram[module$gene[module$category == "MHC-II"],])
all = cbind(all,skinlangerhans[module$gene[module$category == "MHC-II"],])
all = cbind(all,liverkupffer[module$gene[module$category == "MHC-II"],])

all$gene = rownames(all)
all = reshape2::melt(all)
colnames(all)[2] = "Well_ID"
all = merge(all,time,by = "Well_ID",all.x = T)

all$time[all$time != "Adult"] = "Embryo"
all$tissue = ""
all$tissue[all$Well_ID %in% colnames(brain)] = "brainmicroglia"
all$tissue[all$Well_ID %in% colnames(skinpram)] = "skinpram"
all$tissue[all$Well_ID %in% colnames(skinlangerhans)] = "skinlangerhans"
all$tissue[all$Well_ID %in% colnames(liverkupffer)] = "liverkupffer"

all$new = paste(all$time,all$tissue,sep = "_")
all$value = log2(all$value+1)
all$time = factor(all$time, levels=c('Embryo','Adult'))

all1 = all[all$new %in% c("Adult_skinpram","Embryo_skinpram"),]
all2 = all[all$new %in% c("Adult_brainmicroglia","Embryo_brainmicroglia"),]
all3 = all[all$new %in% c("Adult_skinlangerhans","Embryo_skinlangerhans"),]
all4 = all[all$new %in% c("Adult_liverkupffer","Embryo_liverkupffer"),]

p1 = all1%>%ggplot(aes(time,value))+
  geom_violin(aes(fill=time),scale = "width")+
  facet_grid(all1$gene~.,scales = "free")+
  #scale_fill_manual(values = as.character(big.mac.color1),breaks = names(big.mac.color1))+
  scale_x_discrete("")+
  scale_y_continuous("")+
  scale_fill_manual(values = c("#546de5", "#ff4757"))+
  theme_bw()+
  ggtitle("skinpram")+
  theme(strip.background = element_rect(color = "white", fill = "white"),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(size = 12,angle = 360 ),
        legend.position = "none")
ggsave("E:\\adult/mhc2_skinpram.png",p1,width = 4,height = 10,dpi = 300)

###
all = brain[module$gene[module$category == "microglia"],]
all$gene = rownames(all)
all = reshape2::melt(all)
colnames(all)[2] = "Well_ID"
all = merge(all,time,by = "Well_ID",all.x = T)

all$time[all$time != "Adult"] = "Embryo"
all$tissue = "brainmicroglia"

all$new = paste(all$time,all$tissue,sep = "_")
all$value = log2(all$value+1)
all$time = factor(all$time, levels=c('Embryo','Adult'))

p1 = all%>%ggplot(aes(time,value))+
  geom_violin(aes(fill=time),scale = "width")+
  facet_grid(all$gene~.,scales = "free")+
  #scale_fill_manual(values = as.character(big.mac.color1),breaks = names(big.mac.color1))+
  scale_x_discrete("")+
  scale_y_continuous("")+
  scale_fill_manual(values = c("#546de5", "#ff4757"))+
  theme_bw()+
  ggtitle("brainmicroglia")+
  theme(strip.background = element_rect(color = "white", fill = "white"),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(size = 12,angle = 360 ),
        legend.position = "none")
ggsave("E:\\adult/microglia_brainmicroglia.pdf",p1,width = 4,height = 10,dpi = 300)

###
all = skinpram[module$gene[module$category == "PraM"],]
all$gene = rownames(all)
all = reshape2::melt(all)
colnames(all)[2] = "Well_ID"
all = merge(all,time,by = "Well_ID",all.x = T)

all$time[all$time != "Adult"] = "Embryo"
all$tissue = "skinpram"

all$new = paste(all$time,all$tissue,sep = "_")
all$value = log2(all$value+1)
all$time = factor(all$time, levels=c('Embryo','Adult'))

p1=all%>%ggplot(aes(time,value))+
  geom_violin(aes(fill=time),scale = "width")+
  facet_grid(all$gene~.,scales = "free")+
  #scale_fill_manual(values = as.character(big.mac.color1),breaks = names(big.mac.color1))+
  scale_x_discrete("")+
  scale_y_continuous("")+
  scale_fill_manual(values = c("#546de5", "#ff4757"))+
  theme_bw()+
  ggtitle("skinpram")+
  theme(strip.background = element_rect(color = "white", fill = "white"),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(size = 12,angle = 360 ),
        legend.position = "none")
ggsave("E:\\adult/pram_skinpram.pdf",p1,width = 4,height = 12,dpi = 300)

###

all = skinlangerhans[module$gene[module$category == "langerhans"],]
all$gene = rownames(all)
all = reshape2::melt(all)
colnames(all)[2] = "Well_ID"
all = merge(all,time,by = "Well_ID",all.x = T)

all$time[all$time != "Adult"] = "Embryo"
all$tissue = "skinpram"

all$new = paste(all$time,all$tissue,sep = "_")
all$value = log2(all$value+1)
all$time = factor(all$time, levels=c('Embryo','Adult'))

p1=all%>%ggplot(aes(time,value))+
  geom_violin(aes(fill=time),scale = "width")+
  facet_grid(all$gene~.,scales = "free")+
  #scale_fill_manual(values = as.character(big.mac.color1),breaks = names(big.mac.color1))+
  scale_x_discrete("")+
  scale_y_continuous("")+
  scale_fill_manual(values = c("#546de5", "#ff4757"))+
  theme_bw()+
  ggtitle("skinlangerhans")+
  theme(strip.background = element_rect(color = "white", fill = "white"),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(size = 12,angle = 360 ),
        legend.position = "none")
ggsave("E:\\adult/langerhans_skinlangerhans.pdf",p1,width = 4,height = 6,dpi = 300)

###
all = liverkupffer[module$gene[module$category == "kupffer"],]
all$gene = rownames(all)
all = reshape2::melt(all)
colnames(all)[2] = "Well_ID"
all = merge(all,time,by = "Well_ID",all.x = T)

all$time[all$time != "Adult"] = "Embryo"
all$tissue = "liverkupffer"

all$new = paste(all$time,all$tissue,sep = "_")
all$value = log2(all$value+1)
all$time = factor(all$time, levels=c('Embryo','Adult'))

p1=all%>%ggplot(aes(time,value))+
  geom_violin(aes(fill=time),scale = "width")+
  facet_grid(all$gene~.,scales = "free")+
  #scale_fill_manual(values = as.character(big.mac.color1),breaks = names(big.mac.color1))+
  scale_x_discrete("")+
  scale_y_continuous("")+
  scale_fill_manual(values = c("#546de5", "#ff4757"))+
  theme_bw()+
  ggtitle("liverkupffer")+
  theme(strip.background = element_rect(color = "white", fill = "white"),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(size = 12,angle = 360 ),
        legend.position = "none")
ggsave("E:\\adult/kupffer_liverkupffer.png",p1,width = 4,height = 6,dpi = 300)

###
table1 = adult.metatable[adult.metatable$subtype == "microglia" & adult.metatable$tissue == "Brain",]
table1 = table1[table1$time %in% c("Adult","week27"),]
brain = the_mat.combine[intersect(use.gene,rownames(the_mat.combine)),table1$Well_ID]
brain = as.matrix(brain)
brain = as.data.frame(brain)

##
table1 = adult.metatable[adult.metatable$subtype == "PraM" & adult.metatable$tissue == "Skin",]
table1 = table1[table1$time %in% c("Adult","week23","week20"),]
skinpram = the_mat.combine[intersect(use.gene,rownames(the_mat.combine)),table1$Well_ID]
skinpram = as.matrix(skinpram)
skinpram = as.data.frame(skinpram)

##
table1 = adult.metatable[adult.metatable$subtype == "langerhans" & adult.metatable$tissue == "Skin",]
table1 = table1[table1$time %in% c("Adult","week12"),]

skinlangerhans = the_mat.combine[intersect(use.gene,rownames(the_mat.combine)),table1$Well_ID]
skinlangerhans = as.matrix(skinlangerhans)
skinlangerhans = as.data.frame(skinlangerhans)

##
table1 = adult.metatable[adult.metatable$subtype == "Kupffer_cell" & adult.metatable$tissue == "Liver",]
table1 = table1[table1$time %in% c("Adult","week12"),]

liverkupffer = the_mat.combine[intersect(use.gene,rownames(the_mat.combine)),table1$Well_ID]
liverkupffer = as.matrix(liverkupffer)
liverkupffer = as.data.frame(liverkupffer)

###
use.gene = c("CD74","HLA-DRA","HLA-DRB1","P2RY12","TMEM119","CD83","MRC1","CXCL8",
             "VEGFA","CD207","CD1A","TACSTD2")
all = cbind(brain[use.gene,],skinpram[use.gene,])
all = cbind(all,skinlangerhans[use.gene,])


all$gene = rownames(all)
all = reshape2::melt(all)
colnames(all)[2] = "Well_ID"
all = merge(all,time,by = "Well_ID",all.x = T)

all$time[all$time != "Adult"] = "Embryo"
all$tissue = ""
all$tissue[all$Well_ID %in% colnames(brain)] = "brainmicroglia"
all$tissue[all$Well_ID %in% colnames(skinpram)] = "skinpram"
all$tissue[all$Well_ID %in% colnames(skinlangerhans)] = "skinlangerhans"


all$new = paste(all$time,all$tissue,sep = "_")
all$value = log2(all$value+1)
all$time = factor(all$time, levels=c('Embryo','Adult'))

all$gene = factor(all$gene,levels = use.gene,ordered = T)
all$new = factor(all$new,levels = c("Embryo_brainmicroglia","Adult_brainmicroglia","Embryo_skinpram","Adult_skinpram",
                                    "Embryo_skinlangerhans","Adult_skinlangerhans"),ordered = T)

p1 = all%>%ggplot(aes(new,value))+
  geom_violin(aes(fill=time),scale = "width")+
  facet_grid(all$gene~.,scales = "free")+
  #scale_fill_manual(values = as.character(big.mac.color1),breaks = names(big.mac.color1))+
  scale_x_discrete("")+
  scale_y_continuous("")+
  scale_fill_manual(values = c("#546de5", "#ff4757"))+
  theme_bw()+
  theme(strip.background = element_rect(color = "white", fill = "white"),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(size = 12,angle = 360 ),
        legend.position = "none")
ggsave("E:\\adult/vio.pdf",p1,width = 4,height = 8,dpi = 300)


#####
##compare brain pram

mc_sc3 = adult.metatable[adult.metatable$subtype == "PraM" & adult.metatable$tissue == "Brain",]
mc_sc3 = mc_sc3[mc_sc3$time %in% c("Adult","week23"),]
mc_sc3$embryo = mc_sc3$time
mc_sc3$tissue = mc_sc3$time

mc_sc3$major_type = "macrophage"

n <- 50
ntime <- 100

pop <- 'macrophage'
tissueA <- "Adult"; cutA <- 15
tissueB <- "week23"; cutB <- 15

get_n_cells_per_pop_per_disease(mc_sc3, n, pop, tissueA, cutA, diag =F)
get_n_cells_per_pop_per_disease(mc_sc3, n, pop, tissueB, cutB, diag =F)

ind_diff <- diff_calc(mc_sc3, the_mat.combine, pop, n, ntime, tissueA, tissueB, cutA, cutB)
ind_diff_stat <- diff_stat(ind_diff, ntime)

f_cut <- 1.5
q_cut <- 0.05
pct_cut <- 0.05
ind_diff_stat_good <- ind_diff_stat[ind_diff_stat$pctA.mean > pct_cut | ind_diff_stat$pctB.mean > pct_cut,]
f <- log2((ind_diff_stat_good$avgB.mean + 0.01)/(ind_diff_stat_good$avgA.mean + 0.01))
q <- log(-(log(ind_diff_stat_good$wilcox.p.mean, 10)) + 1)
dat <- data.frame(g = rownames(ind_diff_stat_good), f = f, q = q)
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
deg = discard_gene3(deg)

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
  annotate(geom="text", x=-2, y=3, label=tissueA,  fontface="bold",colour='#2CA02CFF', size=5)+
  annotate(geom="text", x=2, y=3, label=tissueB,  fontface="bold",colour='#FF7F0EFF', size=5)+
  theme_bw()+theme_classic()+xlab("log2(fold change)")+ylab("q value")+
  theme(legend.position = "none",
        axis.text = element_text(size=14))
ggsave("E:\\adult/pram/vo.png",p1,width = 8,height = 8)

#####
table1 = adult.metatable[adult.metatable$subtype == "PraM" & adult.metatable$tissue == "Brain",]
table1 = table1[table1$time %in% c("Adult","week23"),]

skinpram = the_mat.combine[intersect(use.gene,rownames(the_mat.combine)),table1$Well_ID]
skinpram = as.matrix(skinpram)
skinpram = as.data.frame(skinpram)

for (i in unique(module$category)){
  skinpram[i,] = colMeans(skinpram[module$gene[module$category == i],])
}

skinpram["g2m",] = colMeans(skinpram[intersect(rownames(skinpram),as.character(gene_list[1,])),])
skinpram = skinpram[c("g2m",unique(module$category)),]


skinpram = as.data.frame(t(skinpram))
skinpram$Well_ID = rownames(skinpram)
skinpram = melt(skinpram)
skinpram = merge(skinpram,time,by = "Well_ID",all.x = T)
skinpram$time[skinpram$time != "Adult"] = "Embryo"
tmp = skinpram

for(i in c("g2m",unique(module$category))){
  skinpram = tmp
  skinpram = skinpram[skinpram$variable == i,]
  
  skinpram$time[skinpram$time == "Embryo"] = "cs23"
  skinpram$time[skinpram$time == "Adult"] = "Adult stage"
  skinpram$time = factor(skinpram$time, levels=c('cs23','Adult stage'))
  compaired <- list(c('cs23','Adult stage'))
  
  skinpram$time = factor(skinpram$time, levels=c('cs23','Adult stage'))
  fc = mean(skinpram$value[skinpram$time == "Adult stage"])/mean(skinpram$value[skinpram$time == "cs23"])
  p1 = ggplot(skinpram, aes(x= time, y = log2(value+1),fill = time))+
    geom_boxplot(outlier.color = NA)+
    stat_boxplot(geom ="errorbar",width=0.15)+
    geom_signif(comparisons = compaired,step_increase = 0.1,test =wilcox.test)+
    xlab("")+ylab("")+ylim(0,8)+
    ggtitle(paste0("fc ",round((fc),2)))+
    theme_classic()+
    scale_fill_manual(values = c("#546de5", "#ff4757"))+
    theme(legend.position = 'none')
  ggsave(paste0("E:\\adult/brain_mac1/",i,".pdf"),p1,width = 2,height = 3.5,dpi = 300)
}

#####

ab2833 = read.table("E:\\sb40/sb40_umitab/AB2833.txt")
ab2835 = read.table("E:\\sb40/sb40_umitab/AB2835.txt")
adult_microglia = cbind(ab2833,ab2835)

adult_microglia.seurat = adult_microglia
adult_microglia.seurat = CreateSeuratObject(counts = adult_microglia.seurat)

adult_microglia.seurat[["percent.mt"]] <- PercentageFeatureSet(adult_microglia.seurat, pattern = "^MT-")
adult_microglia.seurat[["percent.col"]] <- PercentageFeatureSet(adult_microglia.seurat, pattern = "^COL")
adult_microglia.seurat[["percent.blood"]] <- PercentageFeatureSet(adult_microglia.seurat, 
                                                          features = c("HBA1","HBA2","HBG1","HBG2","HBB","HBE1","HBZ"))
adult_microglia.seurat[["percent.myh"]] <- PercentageFeatureSet(adult_microglia.seurat, pattern = "^MYH")
adult_microglia.seurat[["percent.myl"]] <- PercentageFeatureSet(adult_microglia.seurat, pattern = "^MYL")

VlnPlot(adult_microglia.seurat, features = c("percent.mt","percent.col","percent.blood","percent.myh","percent.myl"), ncol =5)

tmp = adult_microglia.seurat@meta.data
tmp = tmp[tmp$percent.mt < 50,]

qc_cells = rownames(tmp)

adult_microglia.seurat = adult_microglia[setdiff(rownames(adult_microglia),str_subset(rownames(adult_microglia), pattern ='^MT-')),
                                         intersect(qc_cells,colnames(adult_microglia))]
adult_microglia.seurat = CreateSeuratObject(counts = adult_microglia.seurat)
adult_microglia.seurat <- NormalizeData(adult_microglia.seurat, normalization.method =  "LogNormalize")
adult_microglia.seurat <- FindVariableFeatures(adult_microglia.seurat,selection.method = "vst", nfeatures = 2000)
adult_microglia.seurat <- ScaleData(adult_microglia.seurat)
adult_microglia.seurat <- RunPCA(object = adult_microglia.seurat, pc.genes = VariableFeatures(adult_microglia.seurat))
adult_microglia.seurat <- FindNeighbors(adult_microglia.seurat, dims = 1:30)
adult_microglia.seurat <- FindClusters(adult_microglia.seurat, resolution = 0.2)
adult_microglia.seurat <- RunUMAP(object = adult_microglia.seurat, dims = 1:30, do.fast = TRUE)
DimPlot(adult_microglia.seurat)
adult_microglia.seurat.marker = FindAllMarkers(adult_microglia.seurat)

top20 <- adult_microglia.seurat.marker %>% group_by(cluster) %>% top_n(n =10, wt = avg_log2FC)
DoHeatmap(adult_microglia.seurat,features = c(top20$gene,"P2RY12"),slot = "data")
tmp = adult_microglia.seurat@meta.data
tmp = rownames(adult_microglia.seurat@meta.data[adult_microglia.seurat@meta.data$seurat_clusters == 0,])

adult_microglia.seurat = adult_microglia[setdiff(rownames(adult_microglia),str_subset(rownames(adult_microglia), pattern ='^MT-')),
                                   tmp]
adult_microglia.seurat = CreateSeuratObject(counts = adult_microglia.seurat)

the_mat.combine = merge(x = the_mat.combine, y = adult_microglia.seurat)
the_mat.combine <- NormalizeData(the_mat.combine, normalization.method =  "RC")
the_mat.combine = the_mat.combine@assays$RNA@data

head(adult.metatable)
tmp = as.data.frame(tmp)
head(tmp)
colnames(tmp ) = "Well_ID"
tmp$Amp_batch_ID = "Adult"
tmp$time = "Adult"
tmp$tissue = "Brain"
tmp$subtype = "microglia"

adult.metatable = rbind(adult.metatable,tmp)

mc_sc3 = adult.metatable[adult.metatable$subtype == "microglia" & adult.metatable$time == "week27",]
tmp = adult.metatable[adult.metatable$Amp_batch_ID == "Adult",]
mc_sc3 = rbind(mc_sc3,tmp)

mc_sc3$embryo = mc_sc3$time
mc_sc3$tissue = mc_sc3$time

mc_sc3$major_type = "macrophage"

n <- 50
ntime <- 100

pop <- 'macrophage'
tissueA <- "Adult"; cutA <- 15
tissueB <- "week27"; cutB <- 15

get_n_cells_per_pop_per_disease(mc_sc3, n, pop, tissueA, cutA, diag =F)
get_n_cells_per_pop_per_disease(mc_sc3, n, pop, tissueB, cutB, diag =F)

ind_diff <- diff_calc(mc_sc3,the_mat.combine, pop, n, ntime, tissueA, tissueB, cutA, cutB)
ind_diff_stat <- diff_stat(ind_diff, ntime)

f_cut <- 1.5
q_cut <- 0.05
pct_cut <- 0.05
ind_diff_stat_good <- ind_diff_stat[ind_diff_stat$pctA.mean > pct_cut | ind_diff_stat$pctB.mean > pct_cut,]
f <- log2((ind_diff_stat_good$avgB.mean + 0.01)/(ind_diff_stat_good$avgA.mean + 0.01))
q <- log(-(log(ind_diff_stat_good$wilcox.p.mean, 10)) + 1)
dat <- data.frame(g = rownames(ind_diff_stat_good), f = f, q = q)
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
deg = discard_gene3(deg)

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
  annotate(geom="text", x=-2, y=3, label=tissueA,  fontface="bold",colour='#2CA02CFF', size=5)+
  annotate(geom="text", x=2, y=3, label=tissueB,  fontface="bold",colour='#FF7F0EFF', size=5)+
  theme_bw()+theme_classic()+xlab("log2(fold change)")+ylab("q value")+
  theme(legend.position = "none",
        axis.text = element_text(size=14))

table1 = mc_sc3

skinpram = the_mat.combine[intersect(use.gene,rownames(the_mat.combine)),table1$Well_ID]
skinpram = as.matrix(skinpram)
skinpram = as.data.frame(skinpram)

for (i in unique(module$category)){
  skinpram[i,] = colMeans(skinpram[module$gene[module$category == i],])
}

skinpram["g2m",] = colMeans(skinpram[intersect(rownames(skinpram),as.character(gene_list[1,])),])
skinpram = skinpram[c("g2m",unique(module$category)),]


skinpram = as.data.frame(t(skinpram))
skinpram$Well_ID = rownames(skinpram)
skinpram = melt(skinpram)
skinpram = merge(skinpram,time,by = "Well_ID",all.x = T)
skinpram$time[skinpram$time != "Adult"] = "Embryo"
tmp = skinpram

for(i in c("g2m",unique(module$category))){
  skinpram = tmp
  skinpram = skinpram[skinpram$variable == i,]
  
  skinpram$time[skinpram$time == "Embryo"] = "week27"
  skinpram$time[skinpram$time == "Adult"] = "Adult stage"
  skinpram$time = factor(skinpram$time, levels=c('week27','Adult stage'))
  compaired <- list(c('week27','Adult stage'))
  
  skinpram$time = factor(skinpram$time, levels=c('week27','Adult stage'))
  fc = mean(skinpram$value[skinpram$time == "Adult stage"])/mean(skinpram$value[skinpram$time == "week27"])
  p1 = ggplot(skinpram, aes(x= time, y = log2(value+1),fill = time))+
    geom_boxplot(outlier.color = NA)+
    stat_boxplot(geom ="errorbar",width=0.15)+
    geom_signif(comparisons = compaired,step_increase = 0.1,test =wilcox.test)+
    xlab("")+ylab("")+ylim(0,8)+
    ggtitle(paste0("fc ",round((fc),2)))+
    theme_classic()+
    scale_fill_manual(values = c("#546de5", "#ff4757"))+
    theme(legend.position = 'none')
  ggsave(paste0("E:\\adult/brain_mac_lilab/",i,".pdf"),p1,width = 2,height = 3.5,dpi = 300)
}
