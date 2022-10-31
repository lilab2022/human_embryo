### Fig7 module definition
library(clusterProfiler)
library(Rgraphviz)
library(org.Hs.eg.db)
library(tgstat)
library(metacell)
library(pheatmap)
library(readxl)    
library(ggplot2)
library(Seurat)
source("Funcs_module_analysis.R")
setwd("./path")

##############################################################################

# read files required
TF_gene <- read.xlsx("TFs_DatabaseExtract_v_1.01.xlsx")%>%filter(`Is.TF?`=="Yes")%>%pull(HGNC.symbol)
umi <- read.csv("cleanumi.csv",row.names=1)
load("mat.embryo_mat.Rda")
umi_sc <- object; remove(object)
adult_mc <- read.csv("ADULT_50MC.csv")$MC
load("mac.color.rdata")
load("metatable_0926.Rdata")

# read in forbidden gene list 
lat_list <- read_excel_allsheets("./sheets/TableS3v2.xlsx")
lat_genes <- NULL
for(i in names(lat_list)){
  lat_genes <- c(lat_genes,lat_list[[i]][,1])
}

##############################################################################

# preprocess of raw metatable
metadata_mf <- metatable%>%filter(major%in%c("macrophage"))
metadata_mf$time <- factor(metadata_mf$time,levels=c("cs11","cs12","cs13","cs14","cs18","cs19","cs21","cs23","week9","week10","week11","week12","week13","week16","week19","week20","week23","week27","Adult"))
metadata_mf <- metadata_mf %>% filter(time!="Adult")
ind_mf <- metadata_mf$mc %>% unique %>% setdiff(adult_mc)
umi_sc_mf <- umi_sc@mat[,metadata_mf$Well_ID]
anno <- metadata_mf[,c("mc","subtype")] %>% .[!duplicated(.),]

##############################################################################

# compute fp based on python resulted umi per mc 

umi_mf <- umi[,colnames(umi)%in%ind_mf]
min_total_umi=20
f_g_cov = rowSums(umi_mf) > min_total_umi
length(which(f_g_cov))
umi_mf <- umi_mf[f_g_cov,]
mc_meansize = colSums(umi_mf)
ideal_cell_size = pmin(1000, median(colSums(umi_mf)))
g_fp = t(ideal_cell_size*t(umi_mf)/as.vector(mc_meansize))
g_fp_md <- g_fp
fp_reg = 0.5
g_fp_n = (fp_reg+g_fp_md)/apply(fp_reg+g_fp_md, 1, median)
lfp_mf <- log2(g_fp_n)

##############################################################################

### select genes fp>3

kp_term <- c("^CXCL","^TNF")
kp_genes <- foreach(i=kp_term, .combine = c) %do% grep(i, lat_genes, v=T)
lat_genes <- setdiff(lat_genes,kp_genes)
mf_genes = names(which(apply(g_fp_n, 1, max) > 3)) %>% setdiff(lat_genes) %>% setdiff(TF_gene)

### mc level module heatmap

mf_cor_mc = cor(t(g_fp_n[mf_genes, ]))
diag(mf_cor_mc) = NA
blwtrd_cols = colorRampPalette(c('blue', 'white', 'red'))(101)
module_mc <- pheatmap(pmin(pmax(mf_cor_mc, -0.8), 0.8), clustering_method="ward.D2",
                      cutree_rows=20, cutree_cols=20, treeheight_col=0, treeheight_row=0, 
                      cellwidth=10, cellheight=10, fontsize=10, col=blwtrd_cols, show_colnames=T, angle_col = 90)
ggsave(module_mc, filename = "./module_whole.pdf",width = 35,height = 35,units = "in",device = "pdf")


### module big gene list

row_mf_module <- cutree(module_mc$tree_row,k=20)
module_order <- data.frame(gene=module_mc$tree_row$labels[module_mc$tree_row$order])
module_order[,ncol(module_order)+1]=row_mf_module[match(module_order$gene,names(row_mf_module))]
colnames(module_order)[ncol(module_order)]="Module"
rownames(module_order) <- module_order$gene
module_order$Module <- paste("module",module_order$Module,sep="")
module_order$Module <- module_order$Module%>%factor(levels = unique(module_order$Module))%>%as.numeric()
module_order[which(module_order$gene=="CD209"):nrow(module_order),]$Module <- (module_order[which(module_order$gene=="CD209"):nrow(module_order),]$Module)+1
module_order[which(module_order$gene=="FGB"):nrow(module_order),]$Module <- (module_order[which(module_order$gene=="FGB"):nrow(module_order),]$Module)+1

module_order$category <- "un"
module_order[module_order$Module==module_order[module_order$gene=="P2RY12","Module"],"category"] <- "microglia"
module_order[module_order$Module==module_order[module_order$gene=="CD5L","Module"],"category"] <- "kupffer"
module_order[module_order$Module==module_order[module_order$gene=="DNASE1L3","Module"],"category"] <- "gut"
module_order[module_order$Module==module_order[module_order$gene=="CD74","Module"],"category"] <- "MHC-II"
module_order[module_order$Module==module_order[module_order$gene=="CD163","Module"],"category"] <- "core-mf"
module_order[module_order$Module==module_order[module_order$gene=="CD207","Module"],"category"] <- "langerhans"
module_order[module_order$Module==module_order[module_order$gene=="MMP9","Module"],"category"] <- "osteoclast"
module_order[module_order$Module==module_order[module_order$gene=="AFP","Module"],"category"] <- "yolk-sac"
module_order[module_order$Module==module_order[module_order$gene=="VEGFA","Module"],"category"] <- "PraM"

write.xlsx(module_order,"./sheets/module_raw.xlsx")

##############################################################################

# expand module list 

gf = names(which(apply(g_fp_n, 1, max) > 2)) %>% setdiff(lat_genes) %>% setdiff(TF_gene)

# on selected genes
fp_s <- g_fp_n[gf,]%>%t()%>%as.data.frame()
fp_s$mc <- rownames(fp_s)
fp_s <- merge(fp_s,anno[,c("mc","subtype")])%>%as.data.frame()
fp_s <- apply(fp_s[,2:(ncol(fp_s)-1)],2,function(x){tapply(x,fp_s$subtype,median)})%>%t()%>%as.data.frame()

gene_mg <- FfocG(fp_s,module_order,NULL,"microglia")
gene_kf <- FfocG(fp_s,module_order,NULL,"kupffer")
gene_gut <- FfocG(fp_s,module_order,NULL,"gut")
gene_mf <- FfocG(fp_s,module_order,NULL,"core-mf")
gene_MHC <- FfocG(fp_s,module_order,NULL,"MHC-II")
gene_ys <- FfocG(fp_s,module_order,NULL,"yolk-sac")
gene_PraM <- FfocG(fp_s,module_order,NULL,"PraM")
gene_LC <- FfocG(fp_s,module_order,NULL,"langerhans")
gene_OC <- FfocG(fp_s,module_order,NULL,"osteoclast")

md_pro_gut <- Fmodule_sel(g_fp_n,gene_gut,"gut")
md_pro_mf <- Fmodule_sel(g_fp_n,gene_mf,"core-mf")
md_pro_mg <- Fmodule_sel(g_fp_n,gene_mg,"microglia")
md_pro_MHC <- Fmodule_sel(g_fp_n,gene_MHC,"MHC-II")
md_pro_kf <- Fmodule_sel(g_fp_n,gene_kf,"kupffer")
md_pro_ys <- Fmodule_sel(g_fp_n,gene_ys,"yolk-sac")
md_pro_PraM <- Fmodule_sel(g_fp_n,gene_PraM,"PraM")
md_pro_LC <- Fmodule_sel(g_fp_n,gene_LC,"langerhans")
md_pro_OC <- Fmodule_sel(g_fp_n,gene_OC,"osteoclast")

ls_all <- list(md_pro_mg,md_pro_mf,md_pro_MHC,md_pro_gut,md_pro_kf,md_pro_PraM,md_pro_LC,md_pro_OC)

module_order_new <- df_module_big(ls_all,module_order)
module_order_new <- module_order_new %>% filter(category!="un")
table(module_order_new$category)
anyDuplicated(module_order_new$gene)

##############################################################################

### focus on selected module (expanded)
# Fig S3E correlation heatmap on macrophage module
order_sub <- c("HdpM","YsdM_AFP_high","YsdM_AFP_low","pre_microglia","Adrenalgland_macrophage",
               "intestine CD209+ Mφ","intestine CD207+ Mφ","langerhans","gonad_macrophage",
               "osteoclast","pre_PraM","PraM","microglia","Kupffer_cell","red_pulp")
order_mod <- c("core-mf","MHC-II","gut","langerhans","osteoclast","PraM","microglia","kupffer")
module_order_new$category <- factor(module_order_new$category, levels = order_mod)
module_order_new <- module_order_new[order(module_order_new$category),]
module_order_new$Module <- as.factor(module_order_new$category) %>% as.numeric()
gaps_foc <- table(module_order_new$Module)%>%cumsum() %>% as.numeric()
gene_foc <- module_order_new$gene
mf_cor_mc = cor(t(g_fp_n[gene_foc, ]))
diag(mf_cor_mc) = NA
gaps_foc <- table(module_order_new$Module)%>%cumsum() %>% as.numeric()

module_mc_foc <- pheatmap(pmin(pmax(mf_cor_mc, -1), 1),gaps_row = gaps_foc,gaps_col = gaps_foc,
                          cluster_rows = F,cluster_cols = F,  cellwidth=10, cellheight=10, fontsize=10, 
                          col=blwtrd_cols, show_colnames=F, angle_col = 90, legend_labels = c(-1,0,1)) 
ggsave(module_mc_foc, filename = "./figs/S7/module_foc.pdf",width = 27,height = 27,units = "in",device = "pdf")

##############################################################################
# Fig S7A module expression MC vs gene heatmap
fp_module <- g_fp_n[module_order_new$gene,]%>%as.data.frame()
fp_module$module <- module_order_new$category
fp_module <- apply(fp_module[,1:ncol(fp_module)-1],2,function(x){tapply(x,fp_module$module,mean)})%>%as.data.frame()

annotation_col <- metadata_mf[,c("mc","subtype")] %>% .[!duplicated(.),]
annotation_col$subtype <- factor(annotation_col$subtype, levels = order_sub)
annotation_col <- annotation_col[order(annotation_col$subtype),]
annotation_col <- data.frame(subtype=annotation_col$subtype,row.names = annotation_col$mc)
annotation_col <- annotation_col %>% filter(rownames(.)%in%ind_mf)

ann_colors <- list(subtype=mac.color)
lfp_g <- lfp_mf[module_order_new$gene,rownames(annotation_col)]%>%as.data.frame()
gaps_foc <- table(module_order_new$Module)%>%cumsum() %>% as.numeric()

bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
p_mc_g <- pheatmap(lfp_g, clustering_method="ward.D2", treeheight_col=0, treeheight_row=0, 
                   cluster_rows = F, cluster_cols = F,
                   scale = "row",fontsize_row=4, gaps_row = gaps_foc,
                   color =c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
                   show_colnames=F, show_rownames = F,
                   cellwidth = 0.5, cellheight = 3,
                   breaks = bk,legend_breaks = seq(-2,2,1),
                   legend = F, annotation_legend = F,
                   annotation_col = annotation_col, annotation_colors = ann_colors) 
ggsave(p_mc_g,filename = "./figs/S7/module_mc_g.pdf",width = 11,height = 8)
### version2 highlighting genes
p_mc_g_d <- pheatmap(lfp_g, clustering_method="ward.D2", treeheight_col=0, treeheight_row=0, 
                     cluster_rows = F, cluster_cols = F,
                     scale = "row",fontsize_row=6, gaps_row = gaps_foc,
                     color =c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
                     show_colnames=F, show_rownames = T,
                     cellwidth = 0.5, cellheight = 5,
                     breaks = bk,legend_breaks = seq(-2,2,2),
                     legend_labels = c("-2","0","2"),
                     annotation_col = annotation_col, annotation_colors = ann_colors) 
ggsave(p_mc_g_d, filename="./figs/S7/module_mc_g_detail.pdf",width = 15,height = 15)

##############################################################################

# Fig 7A subtype vs module score precise 
sub_mod <- fp_module %>% t() %>% as.data.frame() %>% mutate(mc=rownames(.)) %>% merge(anno)
sub_mod_score <- apply(sub_mod[,c(2:9)],2,function(x){tapply(x,sub_mod$subtype,mean)}) %>% t() %>% as.data.frame()
sub_mod_score <- sub_mod_score[order_mod,order_sub] %>% log2
ann_sub <- data.frame(row.names = colnames(sub_mod_score),subtype=colnames(sub_mod_score))
ann_sub_color <- list(subtype=mac.color)

bk <- c(seq(-2,2,by=0.01))
p_sub_mod <- pheatmap(sub_mod_score, clustering_method="ward.D2", treeheight_col=0, treeheight_row=0, 
                      cluster_rows = F, cluster_cols = F,
                      scale = "row",
                      fontsize_row=15, 
                      color =c(colorRampPalette(colors = c("navy","white"))(length(bk)/2),colorRampPalette(colors = c("white","firebrick3"))(length(bk)/2)),
                      show_colnames=F, show_rownames = T,
                      cellwidth = 20, cellheight = 20,
                      breaks = bk,
                      legend_breaks = seq(-2,2,2),
                      annotation_legend = F,
                      annotation_col = ann_sub, annotation_colors = ann_sub_color) 
ggsave(p_sub_mod, filename="./figs/7A/module_sub_overall.pdf",width = 7,height = 4)

##############################################################################

# Fig 3D GO enrichment 
module_order_new$category <- as.character(module_order_new$category)
ls_term_new <- list()
ls_GO_new <- list()
for(i in unique(module_order_new$category)){
  gene_up <- module_order_new%>%filter(category==i)%>%pull(gene)
  gene_up_id = na.omit(mapIds(x = org.Hs.eg.db,
                              keys = gene_up,
                              keytype = "SYMBOL",
                              column = "ENTREZID"))
  ls_term_new[[unique(module_order_new[module_order_new$category==i,"category"])]] <- gene_up_id
  ls_GO_new[[unique(module_order_new[module_order_new$category==i,"category"])]] = enrichGO(gene = gene_up_id,
                                                                                            OrgDb = org.Hs.eg.db,
                                                                                            keyType = "ENTREZID",
                                                                                            ont = "ALL",
                                                                                            pvalueCutoff = 0.05,
                                                                                            qvalueCutoff = 0.05)
}

GO_term <- list()
for(i in 1:length(ls_GO_new)){
  GO_term[[i]] <- ls_GO_new[[i]]@result[,c("Description","Count","qvalue")] %>% filter(qvalue<0.05)
  GO_term[[i]] <- GO_term[[i]][order(GO_term[[i]]$Count,decreasing = T),] %>% head(3)
  GO_term[[i]]$module <- names(ls_GO_new)[i]
}
GO_term <- do.call(rbind,GO_term) %>% pull(Description)

extra_term <- c("anchored component of plasma membrane","pattern recognition receptor activity","scavenger receptor activity",
                "regulation of cell-cell adhesion","ossification", "bone resorption", "positive regulation of cell adhesion",
                "phagocytic vesicle","regulation of acute inflammatory response","regulation of myeloid leukocyte mediated immunity",
                "vascular endothelial growth factor production","regulation of angiogenesis", "regulation of vasculature development",
                "clathrin-coated vesicle")
ex_term <- c("antigen processing and presentation of exogenous peptide antigen","antigen processing and presentation of exogenous antigen",
             "positive regulation of cell adhesion","response to lipopolysaccharide","external side of plasma membrane","sterol metabolic process",
             "cellular response to lipopolysaccharide","steroid metabolic process")
GO_term <- unique(c(GO_term,extra_term)) %>% setdiff(ex_term)
### GO heatmap
GO_foc <- list()
for(i in 1:length(ls_GO_new)){
  GO_foc[[i]] <- ls_GO_new[[i]]@result[,c("qvalue","Description")]%>%filter(Description%in%GO_term)
  GO_foc[[i]] <- data.frame(qvalue=GO_foc[[i]]$qvalue,row.names = GO_foc[[i]]$Description)%>%t()%>%as.data.frame()
  rownames(GO_foc[[i]]) <- names(ls_GO_new)[i]
}
GO_foc <- do.call(function(...) {
  tmp <- plyr::rbind.fill(...)
  rownames(tmp) <- sapply(GO_foc, function(i){
    rownames(i)
  })
  return(tmp)
}
,GO_foc)%>%t()%>%as.data.frame()
GO_foc[is.na(GO_foc)] <- 1
GO_foc <- -log10(GO_foc)
GO_foc <- GO_foc[,order_mod]

#reorder module and GO term v1: correlation clustering
GO_trail <- GO_foc
GO_trail[GO_trail>0] <- 1

tmp <- c("microglia","core-mf","kupffer","PraM","gut","MHC-II","langerhans","osteoclast")
GO_trail <- GO_trail[,tmp]
order_term <- NULL
tmp <- rownames(GO_trail)
for(i in 1:ncol(GO_trail)){
  n <- length(which(GO_trail[,i]==1))
  order_term <- c(order_term,rownames(GO_trail)[which(GO_trail[,i]==1)])
  tmp <- setdiff(tmp,order_term)
}
order_term <- order_term[!duplicated(order_term)]
order_tmp <- c("microglia","core-mf","kupffer","PraM","gut","MHC-II","langerhans","osteoclast")

# manually order go term
order_term <- c("antigen processing and presentation","blood coagulation",
                "hemostasis","receptor-mediated endocytosis","anchored component of plasma membrane","phagocytic vesicle",
                "pattern recognition receptor activity","cargo receptor activity","scavenger receptor activity",
                "regulation of plasma lipoprotein particle levels",
                "response to molecule of bacterial origin","regulation of angiogenesis","regulation of vasculature development",
                "vascular endothelial growth factor production","myeloid leukocyte differentiation","regulation of acute inflammatory response","regulation of leukocyte cell-cell adhesion","regulation of cell-cell adhesion",
                "regulation of myeloid leukocyte mediated immunity","cell killing",
                "clathrin-coated vesicle","bone resorption","tissue remodeling","cellular response to peptide","ossification")
bk=seq(0,4,0.01)
GO_foc <- GO_foc[order_term, order_tmp]
p_go <- pheatmap(GO_foc, clustering_method="ward.D2", treeheight_col=0, treeheight_row=0, 
                 fontsize_row=18, fontsize_col=18,border_color = "black",
                 color = colorRampPalette(colors = c("white","red"))(length(bk)),
                 show_colnames=T, angle_col = 45,
                 cluster_cols = F,cluster_rows = F,
                 breaks = bk,legend_breaks = seq(0,4,1),
                 cellwidth = 20, cellheight = 20) 
ggsave(p_go,filename = "./figs/7B/heatmap_GO.pdf", width = 10, height = 10)

##############################################################################

### single cell multi score
# cell phase evaluated based on CellCycleScoring function
cellcycle_genes <- read.xlsx("~/metacell/SupplementaryTable/G2M_list.xlsx")
gene_G2M <- cellcycle_genes$Gene_G2M %>% na.omit()
gene_S <- cellcycle_genes$Gene_S %>% na.omit()

rownames(metatable) <- metatable$Well_ID
umi_sc_all <- umi_sc@mat[,metatable$Well_ID]
Seurat <- CreateSeuratObject(counts = umi_sc_all, meta.data = metatable)
Seurat <- NormalizeData(Seurat)

Seurat <- CellCycleScoring(Seurat, s.features = gene_S, g2m.features = gene_G2M, set.ident = TRUE)
head(Seurat[[]])
metadata_cellcycle <- Seurat@meta.data[,c("Well_ID","Phase","S.Score","G2M.Score")]
length(which(metadata_cellcycle$Phase=="G2M"))
metadata_mf <- merge(metadata_mf,metadata_cellcycle)
metatable <- merge(metatable,metadata_cellcycle)

# proliferating cells fraction of each metacell
G2M_count_mc <- tapply(metatable$Phase, metatable$mc, function(x){length(which(x=="G2M"))})
G2M_total_mc <- tapply(metatable$Phase, metatable$mc, length)
identical(names(G2M_count_mc),names(G2M_total_mc))
G2M_pct_mc <- data.frame(mc=names(G2M_count_mc/G2M_total_mc),G2M_pct=G2M_count_mc/G2M_total_mc)
write.table(G2M_pct_mc,"./sheets/G2M_pct_mc.csv")

# single cell module score
load("all_cell_normalized_tab.Rdata")

umi_module <- a_reduced[module_order_new$gene,] %>% as.data.frame()
umi_module$category <- module_order_new$category
cl <-makeCluster(30)
clusterExport(cl,c("umi_module"))
umi_module_mean <- parApply(cl,umi_module[,-ncol(umi_module)],2,function(x){tapply(x,umi_module$category,mean)})
stopCluster(cl)
# umi_module_mean <- log1p(umi_module_mean)

umi_module_pct <- (umi_module_mean/10) %>% t() %>% as.data.frame() %>% mutate(Well_ID=rownames(.))
umi_module_pct <- merge(umi_module_pct, metadata_mf[,c("Well_ID","tissue","time","embryo","metacell","subtype","major","Phase","S.Score","G2M.Score")])
rownames(umi_module_pct) <- umi_module_pct$Well_ID
write.table(umi_module_pct, "./sheets/umi_module_pct.csv")

##############################################################################

# Figure S7 C-D boxplot of MHC-II & core-mf module score of each macrophage subtype

### MHC module score
box_mhc <- ggplot(data=umi_module_pct,aes(x=reorder(subtype, `MHC-II`, FUN = median),y=`MHC-II`))+
  geom_violin(aes(fill=subtype),scale="width")+
  geom_boxplot(aes(fill=subtype),width=0.2,outlier.size = 0.1,color="black",size=0.4)+
  scale_fill_manual(values=mac.color)+
  theme(axis.text.x = element_text(angle=45,vjust=0.8,hjust=0.9,size=10,color="black"))+
  xlab("")+
  ylab("MHC-II score")+
  theme(legend.position = "None",
        axis.text = element_text(size=20,color = "black"),
        axis.title = element_text(size=20),
        axis.text.x = element_blank())

ggsave(box_mhc, file="./figs/S7/MHC_II_boxplot.pdf",width = 8,height = 4)


mhc2_order <- ggplot_build(box_mhc)[["layout"]][["panel_params"]][[1]][["x"]][["breaks"]]

### core-mf module score 
tmp <- umi_module_pct
tmp$subtype <- factor(tmp$subtype,levels = mhc2_order)

box_core.mf <- ggplot(data=tmp,aes(x=subtype,y=`core-mf`))+
  geom_violin(aes(fill=subtype),scale="width")+
  geom_boxplot(aes(fill=subtype),width=0.2,outlier.size = 0.1,color="black",size=0.4)+
  scale_fill_manual(values=mac.color)+
  theme(axis.text.x = element_text(angle=45,vjust=0.8,hjust=0.9,size=10,color="black"))+
  xlab("")+
  ylab("core-mf score")+
  theme(legend.position = "None",
        axis.text = element_text(size=20,color = "black"),
        axis.title = element_text(size=20),
        axis.text.x = element_blank())

ggsave(box_core.mf, file="./figs/S7/core-mf_boxplot.pdf",width = 8,height = 4)

##############################################################################

# Fig 7B quintile plot on module, MHC-II, core-mf, G2M, time 

# compute y scale max
sub_quintile <- c("Kupffer_cell","microglia","intestine CD209+ Mφ","osteoclast","langerhans","PraM","pre_PraM")
max_md(umi_module_pct,sub_quintile)

# plt quintile version1 (except pre-PraM & PraM for Proangiogenic module)
plt_quintile_v1(umi_module_pct,"kupffer","Kupffer_cell")
plt_quintile_v1(umi_module_pct,"microglia","microglia")
plt_quintile_v1(umi_module_pct,"gut","intestine CD209+ Mφ")
plt_quintile_v1(umi_module_pct,"osteoclast","osteoclast")
plt_quintile_v1(umi_module_pct,"langerhans","langerhans")
# plt quintile with significance level for langerhans and CD209+ mf
plt_quintile_signif(umi_module_pct,"gut","intestine CD209+ Mφ")
plt_quintile_signif(umi_module_pct,"langerhans","langerhans")
# plt quintile version1 (pre-PraM & PraM for Proangiogenic module)
umi_PraM <- umi_module_pct %>% mutate(umi_module_pct$subtype)
colnames(umi_PraM)[colnames(umi_PraM)=="subtype"] <- "subtype_ori"
colnames(umi_PraM)[ncol(umi_PraM)] <- "subtype"
umi_PraM[umi_PraM$subtype%in%c("pre_PraM","PraM"),"subtype"] <- "pre_PraM & PraM"
plt_quintile_PraM(umi_PraM,"PraM","pre_PraM & PraM")

##############################################################################

# Fig S7D correlation of module score to MHC-II score
### gut vs MHC-II correlation in gut mf
mc_sel <- anno %>% filter(subtype=="intestine CD209+ Mφ") %>% filter(mc%in%ind_mf) %>% pull(mc)
corr.md <- fp_module[c("MHC-II","gut"),mc_sel] %>% t %>% as.data.frame
r.squared <- lm(corr.md$`MHC-II`~corr.md$`gut`) %>% summary %>% .$r.squared %>% round(4)
lm(corr.md$`gut`~corr.md$`MHC-II`) %>% summary

p <- ggplot(corr.md,aes(x=`MHC-II`,y=`gut`))+
  geom_point(color=mac.color["intestine CD209+ Mφ"],size=3)+
  geom_smooth(method = "glm")+
  ggtitle(paste("R^2=",r.squared,sep=""))+
  ylab("Gut score")+
  xlab("MHC-II score")+
  theme(plot.title =element_text(hjust=0.5,size=25),
        axis.title = element_text(size=25),
        axis.text = element_text(size=20,color="black"))
ggsave(p,filename = "./figs/S7/lm_gut_mhc2.pdf", width = 6, height = 6)

### microglia vs MHC-II module correlation in microglia
mc_sel <- anno %>% filter(subtype=="microglia") %>% filter(mc%in%ind_mf) %>% pull(mc)
corr.md <- fp_module[c("MHC-II","microglia"),intersect(mc_sel,colnames(fp_module))] %>% t %>% as.data.frame
r.squared <- lm(corr.md$`MHC-II`~corr.md$`microglia`) %>% summary %>% .$r.squared %>% round(4)
lm(corr.md$`MHC-II`~corr.md$`microglia`) 

p <- ggplot(corr.md,aes(x=`MHC-II`,y=`microglia`))+
  geom_point(color=mac.color["microglia"],size=3)+
  geom_smooth(method = "glm")+
  ggtitle(paste("R^2=",r.squared,sep=""))+
  ylab("Microglia score")+
  xlab("MHC-II score")+
  theme(plot.title =element_text(hjust=0.5,size=25),
        axis.title = element_text(size=25),
        axis.text = element_text(size=20,color="black"))
ggsave(p,filename = "./figs/S7/lm_mg_mhc2.pdf", width = 6, height = 6)

### kupffer vs MHC-II module correlation in kupffer
mc_sel <- anno %>% filter(subtype=="Kupffer_cell") %>% filter(mc%in%ind_mf) %>% pull(mc)
corr.md <- fp_module[c("MHC-II","kupffer"),intersect(mc_sel,colnames(fp_module))] %>% t %>% as.data.frame
r.squared <- lm(corr.md$`MHC-II`~corr.md$`kupffer`) %>% summary %>% .$r.squared %>% round(4)
lm(corr.md$`MHC-II`~corr.md$`kupffer`) 

p <- ggplot(corr.md,aes(x=`MHC-II`,y=`kupffer`))+
  geom_point(color=mac.color["Kupffer_cell"],size=3)+
  geom_smooth(method = "glm")+
  ggtitle(paste("R^2=",r.squared,sep=""))+
  ylab("Kupffer score")+
  xlab("MHC-II score")+
  theme(plot.title =element_text(hjust=0.5,size=25),
        axis.title = element_text(size=25),
        axis.text = element_text(size=20,color="black"))
ggsave(p,filename = "./figs/S7/lm_kf_mhc2.pdf", width = 6, height = 6)

### osteoclast vs MHC-II module correlation in osteoclast
mc_sel <- anno %>% filter(subtype=="osteoclast") %>% filter(mc%in%ind_mf) %>% pull(mc)
corr.md <- fp_module[c("MHC-II","osteoclast"),intersect(mc_sel,colnames(fp_module))] %>% t %>% as.data.frame
r.squared <- lm(corr.md$`MHC-II`~corr.md$`osteoclast`) %>% summary %>% .$r.squared %>% round(5)
lm(corr.md$`MHC-II`~corr.md$`osteoclast`) 

p <- ggplot(corr.md,aes(x=`MHC-II`,y=`osteoclast`))+
  geom_point(color=mac.color["osteoclast"],size=3)+
  geom_smooth(method = "glm")+
  ggtitle(paste("R^2=",r.squared,sep=""))+
  ylab("Osteoclast score")+
  xlab("MHC-II score")+
  theme(plot.title =element_text(hjust=0.5,size=25),
        axis.title = element_text(size=25),
        axis.text = element_text(size=20,color="black"))
ggsave(p,filename = "./figs/S7/lm_oc_mhc2.pdf", width = 6, height = 6)

### PraM vs MHC-II correlation in PraM & pre_PraM
mc_sel <- anno %>% filter(subtype%in%c("pre_PraM","PraM")) %>% filter(mc%in%ind_mf) %>% pull(mc)
corr.md <- fp_module[c("MHC-II","PraM"),intersect(mc_sel,colnames(fp_module))] %>% t %>% as.data.frame %>% mutate(mc=rownames(.))
corr.md <- merge(corr.md,anno)
r.squared <- lm(corr.md$`MHC-II`~corr.md$`PraM`) %>% summary %>% .$r.squared %>% round(4)
lm(corr.md$`MHC-II`~corr.md$`PraM`) 

p <- ggplot(corr.md,aes(x=`MHC-II`,y=`PraM`))+
  geom_point(aes(color=subtype),size=3)+
  geom_smooth(method = "glm")+
  scale_color_manual(values = mac.color[c("PraM","pre_PraM")])+
  ggtitle(paste("R^2=",r.squared,sep=""))+
  ylab("perivascular score")+
  xlab("MHC-II score")+
  theme(plot.title =element_text(hjust=0.5,size=25),
        axis.title = element_text(size=25),
        axis.text = element_text(size=20,color="black"),
        legend.position = "none")
ggsave(p,filename = "./figs/S7/lm_PraM_mhc2.pdf", width = 6, height = 6)

##############################################################################

# Fig S7E PraM quintile by tissue  ---------------------------------------

# Brain + Head  Gonad+ Testicle+ Ovary
# extract module score df for perivascular macrophage
umi_sc_PraM <- umi_module_pct %>% mutate(umi_module_pct$subtype)
colnames(umi_sc_PraM)[colnames(umi_sc_PraM)=="subtype"] <- "subtype_ori"
colnames(umi_sc_PraM)[ncol(umi_sc_PraM)] <- "subtype"
umi_sc_PraM[umi_sc_PraM$subtype%in%c("pre_PraM","PraM"),"subtype"] <- "pre_PraM & PraM"

# filter tissue with low Freq
umi_sc_PraM <- umi_sc_PraM %>% filter(subtype%in%c("pre_PraM & PraM"))
umi_sc_PraM[umi_sc_PraM$tissue%in%c("Brain","Head"),"tissue"] <- "Brain+Head"
umi_sc_PraM[umi_sc_PraM$tissue%in%c("Female gonad","Male gonad"),"tissue"] <- "Gonad"

tissue_tab <- table(umi_sc_PraM$tissue) %>% sort(decreasing = T)
tissue_sel <- names(tissue_tab)[1:which(tissue_tab==tissue_tab["Brain+Head"])]
umi_sc_PraM <- umi_sc_PraM %>% filter(tissue%in%tissue_sel)
umi_sc_PraM$tissue <- factor(umi_sc_PraM$tissue,levels = tissue_sel)
umi_sc_PraM <- umi_sc_PraM[order(umi_sc_PraM$tissue),]
ggplot(umi_sc_PraM)+geom_bar(aes(tissue))+theme(axis.text.x = element_text(angle=45,vjust = 0.8,hjust = 0.7))

ls_PraM <- list()
for(i in unique(umi_sc_PraM$tissue)){
  ls_PraM[[i]] <- umi_sc_PraM %>% filter(tissue==i)
}

for(i in 1:length(ls_PraM)){
  plt_quintile_tissue(ls_PraM[[i]],"PraM","pre_PraM & PraM")
}
