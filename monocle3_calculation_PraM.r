## monocle3 time

rm(list = ls())


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




require('anndata')
require('chameleon')
require('pheatmap')
require('pracma')
require('stats')

library(dplyr)
library(RcppRoll)
library(matrixStats)
library(monocle3)
library(Seurat)
library(ggplot2)
library(gtools)


data<-read.csv("G://embryo/20220329/cleanumi.csv",header=T,row.names = "X")
seob_data <- CreateSeuratObject(counts = data)
seob_data <- NormalizeData(seob_data, normalization.method =  "LogNormalize")

data<-seob_data@assays$RNA@data
data<-as.matrix(data)
data<-as.data.frame(data)


excluded_genes<-discard_gene2(rownames(data))

data<-data[rownames(data)%in%excluded_genes,]


mcmc1<-read.csv("G://embryo/20220930/subtypes_mf_20220930.csv",header=T) ## loading meta information

#removing adult-associated metacells
adult_mc<-c("X2056","X2350","X282","X871","X891")
mcmc1<-mcmc1[!(mcmc1$MC%in%adult_mc),]

# extracting the metacell IDs used to plot heatmap
all_mac = c(unique(mcmc1$MC[mcmc1$subtype %in% c("YsdM_AFP_low","YsdM_AFP_high","pre_PraM","PraM","microglia","HepM")]))


anno.2 = mcmc1[,c("MC","subtype","major")]
anno.2 = anno.2[!duplicated(anno.2),]
anno.2 = anno.2[anno.2$MC %in% all_mac,]
rownames(anno.2) = anno.2$MC
anno.2 = anno.2[,-1]

tmp<-rownames(anno.2)

################################################
################################################
################################################

### extract umi-tab based on metacell IDs obtained in the last step.

data_sub<-data[,colnames(data)%in%tmp]
data_sub$row_var = rowVars(as.matrix(data_sub[,1:ncol(data_sub)]))


data_sub_ordered<-data_sub[order(data_sub$row_var,decreasing = T),]
data_sub_ordered<-data_sub_ordered[1:2000,] # using the 2000 genes having the largest variance
data_sub_ordered<-data_sub_ordered[,-ncol(data_sub_ordered)] # removing the column of "row_var"

realtime<-read.csv("G://embryo/20220930/subtypes_mf_20220930.csv",header=T)
realtime<-realtime[,c(1,4)] 
realtime<-realtime[realtime$MC%in%tmp,]


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

# restoring intermediate files
write.csv(data.t.t,"G://embryo/20220930/monocle3_data_time_order_PraM.csv") 
write.csv(anno.2,"G://embryo/20220930/monocle3_anno2_PraM.csv")




####################################################

data = read.csv("G://embryo/20220930/monocle3_data_time_order_PraM.csv")
#data = data[-c(1:5),]
rownames(data) = data[,1]

data = data[,-1]
data = as.matrix(data)

anno = read.csv("G://embryo/20220930/monocle3_anno2_PraM.csv",row.names = 1)
anno = anno[colnames(data),]

gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(data,
                         cell_metadata = anno,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")

plot_cells(cds, reduction_method="UMAP") + ggtitle('int.umap')
plot_cells(cds) + ggtitle("label by clusterID")

cds <- cluster_cells(cds)

tmp = names(cds@clusters$UMAP$clusters)
clu = c()
anno$mc = rownames(anno)
anno$subtype[anno$subtype == ""] = "no_name"
for (i in tmp) {
  clu = c(clu,anno$subtype[anno$mc == i])
}

names(clu) = tmp
cds@clusters$UMAP$clusters = clu

png(file="G://embryo/20220930/monocle3_pseudo_time_analysis_PraM.png",width = 1024, height = 800)

plot_cells(cds, show_trajectory_graph = FALSE,group_label_size = 5,
           cell_size = 1) + ggtitle("label by subtype")

dev.off()


cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, group_label_size = 5,cell_size = 0.8,
           label_branch_points = FALSE)
##extrect pseudotime of each mc


cds <- order_cells(cds)

pse = as.data.frame(pseudotime(cds))


write.csv(pse,"G://embryo/20220930/monocle3_pse2_PraM.csv")















