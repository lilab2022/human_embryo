######
###### Fig4B
###### writer: Feng Ruoqing
###### Date : 2022-10-4

microglial.metatable = microglial.metatable[microglial.metatable$tissue %in% c("Skin","Spinalmarrow","Brain","Male gonad"),]
load("E:\\microglial_analysis/microglia_metatable.Rdata")
microglial.metatable$subtype = "microglia"
##skin 
skin1 = metatable[metatable$tissue == "Skin" & metatable$subtype %in% c("PraM","pre_PraM"),]
skin1$tissue = paste0(skin1$tissue,skin1$subtype,sep = "_")
skin1 = skin1[,c("Well_ID","tissue","time","subtype")]

skin2 = metatable[metatable$tissue == "Skin" & metatable$subtype %in% c("langerhans"),]
skin2$tissue = paste0(skin2$tissue,skin2$subtype,sep = "_")
skin2 = skin2[,c("Well_ID","tissue","time","subtype")]

microglial.metatable = microglial.metatable[,c("Well_ID","tissue","time","subtype")]
microglial.metatable = rbind(microglial.metatable,skin1)
microglial.metatable = rbind(microglial.metatable,skin2)

#brain
brain1 = metatable[metatable$tissue == "Brain" & metatable$subtype %in% c("PraM","pre_PraM"),]
brain1$tissue = paste0(brain1$tissue,brain1$subtype,sep = "_")
brain1 = brain1[,c("Well_ID","tissue","time","subtype")]

microglial.metatable = rbind(microglial.metatable,brain1)

##spinalmarrow 
spinalmarrow1 = metatable[metatable$tissue == "Spinalmarrow" & metatable$subtype %in% c("PraM","pre_PraM"),]
spinalmarrow1$tissue = paste0(spinalmarrow1$tissue,spinalmarrow1$subtype,sep = "_")
spinalmarrow1 = spinalmarrow1[,c("Well_ID","tissue","time","subtype")]

microglial.metatable = rbind(microglial.metatable,spinalmarrow1)

##male_gonad
gonad1 = metatable[metatable$tissue == "Male gonad" & metatable$subtype %in% c("osteoclast"),]
gonad1$tissue = "gonad_osteoclast"
gonad1 = gonad1[,c("Well_ID","tissue","time","subtype")]

gonad2 = metatable[metatable$tissue == "Male gonad" & metatable$subtype %in% c("PraM","pre_PraM"),]
gonad2$tissue = "gonad_mac12"
gonad2 = gonad2[,c("Well_ID","tissue","time","subtype")]

gonad3 = metatable[metatable$tissue == "Male gonad" & metatable$subtype %in% c("gonad_macrophage"),]
gonad3$tissue = "gonad_macrophage"
gonad3 = gonad3[,c("Well_ID","tissue","time","subtype")]

microglial.metatable = rbind(microglial.metatable,gonad1)
microglial.metatable = rbind(microglial.metatable,gonad2)
microglial.metatable = rbind(microglial.metatable,gonad3)

use.gene = c("SALL1", "P2RY12", "C3", "TMEM119", "MRC1", "CD163","CD207","DAB2","LYVE1")

tmp_tab = as.data.frame(t(as.matrix(a_reduced[use.gene,])))
tmp_tab$Well_ID = rownames(tmp_tab)

tmp_tab = merge(microglial.metatable,tmp_tab,by = "Well_ID",all.x = T)
#tmp_tab[is.na(tmp_tab)] = 0

tmp_tab = tmp_tab[tmp_tab$time != "Adult",]
tmp_tab = tmp_tab[,c("tissue",use.gene)]
tmp_tab = reshape2::melt(tmp_tab)
tmp_tab = tmp_tab[tmp_tab$tissue %in% c("Brain","Brainpre_PraM_","BrainPraM_",
                                        "Spinalmarrow","Spinalmarrowpre_PraM_","SpinalmarrowPraM_",
                                        "Skin","Skinpre_PraM_","SkinPraM_","Skinlangerhans_"),]
tmp_tab$tissue = factor(tmp_tab$tissue,levels = c("Brain","Brainpre_PraM_","BrainPraM_",
                                                  "Spinalmarrow","Spinalmarrowpre_PraM_","SpinalmarrowPraM_",
                                                  "Skin","Skinpre_PraM_","SkinPraM_","Skinlangerhans_"),ordered = T)

p1 = ggplot(tmp_tab)+
  geom_violin(aes(x = tissue, y = log2(value+1),fill = variable),scale = "width",color = "grey30")+
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
ggsave("E:\\microglial_analysis/microglia_compare_new_normalize.pdf",p1,width = 4,height = 4,dpi = 300)

###
### Fig 4 related to symphony mapping
# import data -------------------------------------------------------------
source("~/metacell/Fig7/script/symphony_utils.R")
library(ggplot2)
library(ggrastr)
library(dplyr)
library(symphony)
setwd("~/metacell/Fig7/symphony/")
load("mat.embryo_mat.Rda")
load("metatable_0926.Rdata")
anno <- metatable[,c("mc","major","subtype")] %>% .[!duplicated(.),]
ref_meta <- metatable

ref_exp_full <- object@mat[,ref_meta$Well_ID]
rownames(ref_meta) <- ref_meta$Well_ID
ref_meta <- ref_meta[colnames(ref_exp_full),]
rm(object)

load("~/metacell/Fig1/data/mac.color.rdata")
load("~/metacell/Fig1/data/major_color.Rdata")
color_con <- c(mac.color,"Ref"="grey")

sub_order <- c("HdpM","YsdM_AFP_high","YsdM_AFP_low","pre_microglia","red_pulp",
               "Kupffer_cell", "pre_PraM","PraM","gonad_macrophage","Adrenalgland_macrophage",
               "intestine CD209+ Mφ","intestine CD207+ Mφ","langerhans","osteoclast","microglia")
color_bar <- c(mac.color[c("pre_PraM","PraM","microglia")], "others"="grey")
mac.color <- mac.color[sub_order]

# plotBasic ---------------------------------------------------------------------

plotBasic2 = function(umap_labels,                # metadata, with UMAP labels in UMAP1 and UMAP2 slots
                      title = 'Query',         # Plot title
                      color.by = 'cell_type',  # metadata column name for coloring
                      facet.by = NULL,         # (optional) metadata column name for faceting
                      color.mapping = NULL,    # custom color mapping
                      legend.position = 'right',
                      pt_size = 2) {  # Show cell type legend
  
  umap_labels_query <- umap_labels[umap_labels[,color.by]!="Ref",]
  umap_labels_ref <- umap_labels[umap_labels[,color.by]=="Ref",]
  
  p = umap_labels %>%
    ggplot(aes(x = UMAP1, y = UMAP2)) + 
    geom_point(data=umap_labels_ref,aes(fill = get(color.by)), size = 2,  shape = 21, color="grey")+
    geom_point(data=umap_labels_query,aes(fill = get(color.by)), size = pt_size,  shape = 21, color="grey30")
  
  if (!is.null(color.mapping)) { p = p + scale_fill_manual(values = color.mapping) }
  
  # Default formatting
  p = p + 
    # theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank())+
    # labs(title = title, color = color.by) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position=legend.position) +
    theme(legend.text = element_text(size=8), legend.title=element_text(size=12)) + 
    guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = 'none')+
    theme(axis.ticks = element_blank())+
    xlab("")+ylab("")
  if(!is.null(facet.by)) {
    p = p + facet_wrap(~get(facet.by)) +
      theme(strip.text.x = element_text(size = 12)) }    
  return(p)
}

###########################################################################
# build reference for major celltype mapping

set.seed(0)
reference = symphony::buildReference(
  ref_exp_full,
  ref_meta,
  vars = NULL,         # variables to integrate over
  theta = c(0.8),
  K = 100,                   # number of Harmony clusters
  verbose = TRUE,            # verbose output
  do_umap = TRUE,            # can set to FALSE if want to run umap separately later
  do_normalize = T,      # set to TRUE if input counts are not normalized yet
  vargenes_method = 'vst',   # method for variable gene selection ('vst' or 'mvp')
  vargenes_groups = NULL, # metadata column specifying groups for variable gene selection 
  topn = 2000,               # number of variable genes to choose per group
  d = 20,                    # number of PCs
  save_uwot_path = './testing_uwot_model_major' # file path to save uwot model
)

save(file="reference_major.rdata",reference)

###########################################################################
# query mapping (Popescu et al. Decoding fetal liver) skin cells 

load("mat_fl_imm.Rda")
load("mat_fs.rda")
load("meta_fl.Rda")
load("meta_fs.rda")

### query exp and metadata

query_exp <- mat_fl_imm
colnames(meta_fl)[colnames(meta_fl)=="cell"] <- "Well_ID"
meta_fl$tissue <- "skin"
meta_fl$sample <- gsub("_.*","",meta_fl$Well_ID)
query_metadata <- meta_fl
rownames(query_metadata) <- query_metadata$Well_ID
identical(rownames(query_metadata),colnames(query_exp))
query_fl = mapQuery(query_exp,             # query gene expression (genes x cells)
                    query_metadata,        # query metadata (cells x attributes)
                    reference,             # Symphony reference object
                    vars = "sample",
                    do_normalize = TRUE,  # perform log(CP10k) normalization on query
                    do_umap = T)        # project query cells into reference UMAP
set.seed(2021)
query_big_fl = knnPredict(query_fl, 
                          reference, 
                          reference$meta_data$major, 
                          k = 30,
                          confidence = TRUE,
                          seed = 1
)
colnames(query_big_fl$meta_data)[(ncol(query_big_fl$meta_data)-1):ncol(query_big_fl$meta_data)] <- c("major_pred_knn","major_pred_knn_prob")

query_anno_fl <- query_big_fl$meta_data
fl_mac_c <- query_anno_fl %>% filter(major_pred_knn=="macrophage") %>% pull(Well_ID)

meta_fl_mac <- meta_fl %>% filter(Well_ID%in%fl_mac_c)
mat_fl_mac <- mat_fl_imm[,meta_fl_mac$Well_ID]
tmp <- read.csv("meta_fl.txt",sep="\t")[,c(1,5,10)]
colnames(tmp) <- c("sample","time","tissue_raw")
meta_fl_mac <- merge(meta_fl_mac, tmp)
meta_fl_mac$time <- gsub(" gestation","",meta_fl_mac$time)

# query mapping (Xu et al.) skin cells 

### query exp and metadata
query_exp <- mat_fs
colnames(meta_fs)[colnames(meta_fs)=="cell"] <- "Well_ID"
meta_fs$tissue <- "skin"
query_metadata <- meta_fs
rownames(query_metadata) <- query_metadata$Well_ID
identical(colnames(query_exp),rownames(query_metadata))
query_fs = mapQuery(query_exp,             # query gene expression (genes x cells)
                    query_metadata,        # query metadata (cells x attributes)
                    reference,             # Symphony reference object
                    vars = "week",        # for fs dataset, each time point represents a sample
                    do_normalize = TRUE,  # perform log(CP10k) normalization on query
                    do_umap = T)        # project query cells into reference UMAP
set.seed(2021)
query_big_fs = knnPredict(query_fs, 
                          reference, 
                          reference$meta_data$major, 
                          k = 30,
                          confidence = TRUE,
                          seed = 1
)

colnames(query_big_fs$meta_data)[(ncol(query_big_fs$meta_data)-1):ncol(query_big_fs$meta_data)] <- c("major_pred_knn","major_pred_knn_prob")

query_anno_fs <- query_big_fs$meta_data
fs_mac_c <- query_anno_fs %>% filter(major_pred_knn=="macrophage") %>% pull(Well_ID)

meta_fs_mac <- meta_fs %>% filter(Well_ID%in%fs_mac_c)
mat_fs_mac <- mat_fs[,meta_fs_mac$Well_ID]

###########################################################################
# decoding buildreference

buildReference_on_UMAP <- function(exp_ref, metadata_ref, vars = NULL, sc2d, K = 100, verbose = FALSE, 
                                   do_umap = TRUE, do_normalize = TRUE, vargenes_method = "vst", 
                                   vargenes_groups = NULL, topn = 2000, tau = 0, theta = 2, 
                                   save_uwot_path = NULL, d = 20, additional_genes = NULL, 
                                   umap_min_dist = 0.1, seed = 111)
{

  set.seed(seed) # for reproducible soft k-means and UMAP
  
  res = list(meta_data = metadata_ref)
  
  if (do_normalize) {
    if (verbose) message('Normalizing')
    exp_ref = normalizeData(exp_ref, 1e4, 'log')
  }
  
  if (verbose) message('Finding variable genes using ', vargenes_method, ' method')
  if (vargenes_method == 'mvp') {
    if (is.null(vargenes_groups)) {
      vargenes_df = findVariableGenes(exp_ref, rep('A', ncol(exp_ref)), num.bin = 20)
    } else { # groups specified
      vargenes_df = findVariableGenes(exp_ref, groups = as.character(metadata_ref[[vargenes_groups]]), 
                                      num.bin = 20)
    }
    var_genes = unique(data.table(vargenes_df)[, head(.SD[order(-.data$gene_dispersion_scaled)], topn), by = .data$group][, .data$symbol])
  } else if (vargenes_method == 'vst') {
    if (is.null(vargenes_groups)) {
      var_genes = vargenes_vst(exp_ref, topn = topn)
    } else { # groups specified
      var_genes = vargenes_vst(exp_ref, groups = as.character(metadata_ref[[vargenes_groups]]), topn = topn)
    }
  } else {
    message("Invalid variable gene selection method. Options are 'vst' or 'mvp'.")
  }
  
  # Add in any additional genes
  if(!is.null(additional_genes)) { 
    if (verbose) message('Adding ', length(additional_genes), ' additional genes')
    var_genes = union(var_genes, additional_genes)   
  }
  if (verbose) message('Total ' , length(var_genes), ' genes for downstream steps')
  
  # Subset gene expression matrix by the desired genes
  exp_ref = exp_ref[var_genes, ]
  
  if (verbose) message('Scaling and PCA')
  vargenes_means_sds = tibble(symbol = var_genes, mean = Matrix::rowMeans(exp_ref))
  vargenes_means_sds$stddev = rowSDs(exp_ref, vargenes_means_sds$mean)
  
  # Scale data
  exp_ref_scaled = scaleDataWithStats(exp_ref, vargenes_means_sds$mean, vargenes_means_sds$stddev, 1)
  
  # PCA
  s = irlba::irlba(exp_ref_scaled, nv = d)
  Z_pca_ref = diag(s$d) %*% t(s$v) # [PCs by cells]
  res$loadings = s$u
  res$vargenes = vargenes_means_sds
  
  # Run Harmony integration
  if (!is.null(vars)) {
    if (verbose) message('Running Harmony integration')
    
    # Run Harmony to harmonize the reference
    ref_harmObj = harmony::HarmonyMatrix(
      data_mat = t(Z_pca_ref),  ## PCA embedding matrix of cells
      meta_data = metadata_ref, ## dataframe with cell labels
      theta = theta,            ## cluster diversity enforcement
      tau = tau,
      vars_use = vars,          ## variable to integrate over
      nclust = K,               ## number of clusters in Harmony model
      max.iter.harmony = 20,
      return_object = TRUE,     ## return the full Harmony model object
      do_pca = FALSE,           ## do not recompute PCs
      verbose = verbose
    )
    
    res$centroids <- t(cosine_normalize_cpp(ref_harmObj$R %*% t(ref_harmObj$Z_corr) , 1))
    res$R <- ref_harmObj$R
    res$betas <- harmony::moe_ridge_get_betas(ref_harmObj)
    res$Z_orig <- Z_pca_ref
    res$Z_corr <- ref_harmObj$Z_corr
    res$K <- K
    res$d <- d
  } else {
    clust_res <- soft_kmeans(Z_pca_ref, K)
    res$centroids <- clust_res$Y
    res$R <- clust_res$R
    res$betas <- NULL
    res$Z_orig <- Z_pca_ref
    res$Z_corr <- Z_pca_ref
  }
  
  # Add row and column names
  colnames(res$Z_orig) = row.names(metadata_ref)
  rownames(res$Z_orig) = paste0("PC_", seq_len(nrow(res$Z_corr)))
  colnames(res$Z_corr) = row.names(metadata_ref)
  rownames(res$Z_corr) = paste0("harmony_", seq_len(nrow(res$Z_corr)))
  
  # Compute reference compression terms
  if (verbose) message('Computing reference compression terms')
  res$cache = compute_ref_cache(res$R, res$Z_corr)
  
  # Compute centroids in harmony PC space
  cluster_sizes = res$cache[[1]] %>% as.matrix()
  centroid_sums = t(res$Z_corr %*% t(res$R)) %>% as.data.frame()
  centroids_pc = sweep(centroid_sums, 1, cluster_sizes, "/")
  colnames(centroids_pc) = paste0("harmony_", seq_len(nrow(res$Z_corr)))
  rownames(centroids_pc) = paste0("centroid_", seq_len(nrow(res$R)))
  res$centroids_pc = centroids_pc
  
  if (do_umap) {
    if (verbose) message('Running UMAP')
    umap <- uwot::umap(
      t(res$Z_corr), n_neighbors = 30, learning_rate = 0.5, init = "laplacian", 
      metric = 'cosine', fast_sgd = FALSE, n_sgd_threads = 1, # for reproducibility
      min_dist = umap_min_dist, n_threads = 4, ret_model = TRUE
    )
    res$umap$embedding = sc2d
    colnames(res$umap$embedding) = c('UMAP1', 'UMAP2')
    
    tmp <- umap
    tmp$embedding <- sc2d
    
    # Since the nn-index component of the uwot model is not able to be saved as an 
    # object, we save the uwot model at a user-defined path.
    if (!is.null(save_uwot_path)) {
      
      # If file already exists, delete it (otherwise will result in an error)
      if (file.exists(save_uwot_path)) {
        if (verbose) message(paste('File already exists at that path... overwriting...'))
        file.remove(save_uwot_path)
      }
      
      model = uwot::save_uwot(tmp, file = save_uwot_path, unload = FALSE, verbose = FALSE)
      res$save_uwot_path = save_uwot_path
      if (verbose) message(paste('Saved uwot model'))
    }
    
    # reference <- res
  }
  return(res)
}

###########################################################################
# symphony based on macrophage 2d on current figure

sc <- read.csv("sc.csv")
rownames(sc) <- sc$Well_ID
sc <- as.matrix(sc[,2:3])
colnames(sc) <- c("V1","V2")
sc_fetal <- intersect(ref_meta[ref_meta$time!="Adult","Well_ID"],rownames(sc))
sc <- sc[sc_fetal,]
dim(sc)

ref_exp_mac <- ref_exp_full[,rownames(sc)]
dim(ref_exp_mac)
ref_meta_mac <- ref_meta[rownames(sc),]
identical(rownames(ref_meta_mac),colnames(ref_exp_mac))

### build reference
set.seed(0)
reference_mac = buildReference_on_UMAP(
  ref_exp_mac,
  ref_meta_mac,
  vars = NULL,         # variables to integrate over
  sc2d = sc,
  theta = c(0.8),
  K = 100,                   # number of Harmony clusters
  verbose = TRUE,            # verbose output
  do_umap = TRUE,            # can set to FALSE if want to run umap separately later
  do_normalize = T,      # set to TRUE if input counts are not normalized yet
  vargenes_method = 'vst',   # method for variable gene selection ('vst' or 'mvp')
  vargenes_groups = NULL, # metadata column specifying groups for variable gene selection 
  topn = 2000,               # number of variable genes to choose per group
  d = 20,                    # number of PCs
  save_uwot_path = './testing_uwot_model_sc2d' # file path to save uwot model
)

save(file="reference_mac.rdata",reference)
# skin mg  ----------------------------------------------------------------

umap_labels_ref = cbind(reference_mac$meta_data, reference_mac$umap$embedding)
### show skin mg location
umap_ref_skin_mg <- umap_labels_ref %>% mutate(sub_skin=.$subtype)
umap_ref_skin_mg[umap_ref_skin_mg$tissue!="Skin","sub_skin"] <- "Ref"
color_skmg <- c("skin microglia"="#1F77B4FF","Ref"="grey")

a1 <- which(umap_ref_skin_mg$tissue=="Skin")
a2 <- setdiff(1:nrow(umap_ref_skin_mg),a1)
umap_ref_skin_mg <- umap_ref_skin_mg[c(a2,a1),]
head(umap_ref_skin_mg)

### 2d on skin microglia
pdf(file = "./reference_mac_skin_mg_sc2d.pdf",width =18, height = 18)
plotBasic2(umap_ref_skin_mg, title = 'Skin microglia', color.by = 'sub_skin',color.mapping = color_con)+theme(legend.position = "None")
dev.off()

# Fig 4C query mapping Popescu et al. on mac ----------------------------------------

identical(colnames(mat_fl_mac),rownames(meta_fl_mac))
query = mapQuery(mat_fl_mac,             # query gene expression (genes x cells)
                 meta_fl_mac,        # query metadata (cells x attributes)
                 reference_mac,             # Symphony reference object
                 vars = "sample",
                 do_normalize = TRUE,  # perform log(CP10k) normalization on query
                 do_umap = T)        # project query cells into reference UMAP
set.seed(2021)
query_sub_fl = knnPredict(query, 
                          reference_mac, 
                          reference_mac$meta_data$subtype, 
                          k = 30,
                          confidence = TRUE,
                          seed = 1
)

colnames(query_sub_fl$meta_data)[(ncol(query_sub_fl$meta_data)-1):ncol(query_sub_fl$meta_data)] <- c("sub_pred_knn","sub_pred_knn_prob")

r_metadata = reference_mac$meta_data
q_metadata = query_sub_fl$meta_data
query_anno_fl_sub <- q_metadata
bar_fl <- data.frame(table(query_anno_fl_sub$sub_pred_knn)[table(query_anno_fl_sub$sub_pred_knn)!=0],Amp="FL")
bar_fl <- bar_fl[order(bar_fl$Freq,decreasing = F),]
bar_fl$Var1 <- factor(bar_fl$Var1,levels = bar_fl$Var1)
ggplot(bar_fl,aes(x=Var1,y=Freq))+
  geom_bar(stat = "identity")+
  coord_flip()+
  theme(axis.title = element_blank())+
  ggtitle("FL predicted subtype distribution")

r_metadata$ref_query = "Ref"
q_metadata$ref_query = 'FL'

r_metadata$sub_pred_knn_prob <- 1

r_metadata <- r_metadata[,c("Well_ID","subtype","sub_pred_knn_prob","ref_query","tissue")]
q_metadata <- q_metadata[,c("Well_ID","sub_pred_knn","sub_pred_knn_prob","ref_query","tissue")]

head(r_metadata)
head(q_metadata)

colnames(r_metadata)[2:4] <- c("cell_type_pred_knn","cell_type_pred_knn_prob","ref_query")
colnames(q_metadata)[2:4] <- c("cell_type_pred_knn","cell_type_pred_knn_prob","ref_query")

meta_data_combined = rbind(r_metadata,q_metadata)
umap_combined = rbind(reference_mac$umap$embedding,query$umap)

umap_combined_labels = cbind(meta_data_combined, umap_combined) 

umap_combined_labels$sort_by <- paste(umap_combined_labels$ref_query,umap_combined_labels$tissue,sep="_")
umap_combined_labels$con <- "Ref"
umap_combined_labels[umap_combined_labels$ref_query=="FL","con"] <- umap_combined_labels[umap_combined_labels$ref_query=="FL","cell_type_pred_knn"]

b1 <- which(umap_combined_labels$ref_query=="Ref")
b2 <- setdiff(1:nrow(umap_combined),b1)
umap_combined_labels <- umap_combined_labels[c(b1,b2),]
umap_combined_labels %>% head(4)
umap_fl <- umap_combined_labels
umap_fl[umap_fl$ref_query=="FL","con"] <- umap_fl[umap_fl$ref_query=="FL","cell_type_pred_knn"]


pdf(file = "./DecodingFL_mapping_sc2d.pdf",width =18, height = 18)
plotBasic2(umap_fl, title = 'DecodingFL_mapping', color.by = 'con',color.mapping = color_con)+theme(legend.position = "None")
dev.off()

# Fig 4C query mapping fetal skin (Xu et al.) on mac ----------------------------------------

identical(colnames(mat_fs_mac),rownames(meta_fs_mac))
query = mapQuery(mat_fs_mac,             # query gene expression (genes x cells)
                 meta_fs_mac,        # query metadata (cells x attributes)
                 reference_mac,             # Symphony reference object
                 vars = "week",
                 do_normalize = TRUE,  # perform log(CP10k) normalization on query
                 do_umap = T)        # project query cells into reference UMAP
set.seed(2021)
query_sub_fs = knnPredict(query, 
                          reference_mac, 
                          reference_mac$meta_data$subtype, 
                          k = 30,
                          confidence = TRUE,
                          seed = 1
)

colnames(query_sub_fs$meta_data)[(ncol(query_sub_fs$meta_data)-1):ncol(query_sub_fs$meta_data)] <- c("sub_pred_knn","sub_pred_knn_prob")

r_metadata = reference_mac$meta_data
q_metadata = query_sub_fs$meta_data
query_anno_fs_sub <- q_metadata
bar_fs <- data.frame(table(query_anno_fs_sub$sub_pred_knn)[table(query_anno_fs_sub$sub_pred_knn)!=0],Amp="FS")
bar_fs <- bar_fs[order(bar_fs$Freq,decreasing = F),]
bar_fs$Var1 <- factor(bar_fs$Var1,levels = bar_fs$Var1)
ggplot(bar_fs,aes(x=Var1,y=Freq))+
  geom_bar(stat = "identity")+
  coord_flip()+
  theme(axis.title = element_blank())+
  ggtitle("fs predicted subtype distribution")

r_metadata$ref_query = "Ref"
q_metadata$ref_query = 'FS'

r_metadata$sub_pred_knn_prob <- 1

r_metadata <- r_metadata[,c("Well_ID","subtype","sub_pred_knn_prob","ref_query","tissue")]
q_metadata <- q_metadata[,c("Well_ID","sub_pred_knn","sub_pred_knn_prob","ref_query","tissue")]

head(r_metadata)
head(q_metadata)

colnames(r_metadata)[2:4] <- c("cell_type_pred_knn","cell_type_pred_knn_prob","ref_query")
colnames(q_metadata)[2:4] <- c("cell_type_pred_knn","cell_type_pred_knn_prob","ref_query")

meta_data_combined = rbind(r_metadata,q_metadata)
umap_combined = rbind(reference_mac$umap$embedding,query$umap)

umap_combined_labels = cbind(meta_data_combined, umap_combined) 

umap_combined_labels$sort_by <- paste(umap_combined_labels$ref_query,umap_combined_labels$tissue,sep="_")
umap_combined_labels$con <- "Ref"
umap_combined_labels[umap_combined_labels$ref_query=="FS","con"] <- umap_combined_labels[umap_combined_labels$ref_query=="FS","cell_type_pred_knn"]

b1 <- which(umap_combined_labels$ref_query=="Ref")
b2 <- setdiff(1:nrow(umap_combined),b1)
umap_combined_labels <- umap_combined_labels[c(b1,b2),]
umap_combined_labels %>% head(4)
umap_fs <- umap_combined_labels
umap_fs[umap_fs$ref_query=="FS","con"] <- umap_fs[umap_fs$ref_query=="FS","cell_type_pred_knn"]


pdf(file = "./DynamicsFS_mapping_sc2d.pdf",width =18, height = 18)
plotBasic2(umap_fs, title = 'DynamicsFS_mapping', color.by = 'con',color.mapping = color_con)+theme(legend.position = "None")
dev.off()

###########################################################################
# experiment data AB2685，AB2686，AB2841，AB2842
# Fig 4C-E &  S4C mapping of MRC1+/-, P2RY12+/- plates in current study

path = "~/metacell/Fig7/data/symphony/"
fileName = dir(path) %>% grep(pattern = ".*txt", ., value = T)
mat_exp <- list()
for(i in unique(fileName)){
  nm <- gsub(".txt","",i)
  mat_exp[[nm]] = read.csv(file = paste(path,i,sep = ""),
                           header = T,sep="\t",stringsAsFactors = F)
}
mat_exp_merge <- do.call(cbind,mat_exp)
mat_exp_sparse <- as.sparse(mat_exp_merge)
meta_exp <- data.frame(row.names=colnames(mat_exp_merge),Well_ID=colnames(mat_exp_merge))
meta_exp$batch <- gsub("\\.W.*","",meta_exp$Well_ID)
# integrating meta info
gating <- data.frame(batch=c("AB2685", "AB2686", "AB2841", "AB2842"), 
                     gating=c("MRC1-P2RY12+","MRC1+P2RY12+-","MRC1-P2RY12+","MRC1+P2RY12+-"),
                     time=c("13 PCW","13 PCW","9 PCW","9 PCW"),
                     tissue=rep("Skin",4))
meta_exp <- merge(meta_exp, gating)
rownames(meta_exp) <- meta_exp$Well_ID
# query mapping
identical(colnames(mat_exp_sparse),rownames(meta_exp))
query = mapQuery(mat_exp_sparse,             # query gene expression (genes x cells)
                 meta_exp,        # query metadata (cells x attributes)
                 reference_mac,             # Symphony reference object
                 do_normalize = TRUE,  # perform log(CP10k) normalization on query
                 do_umap = T)        # project query cells into reference UMAP
set.seed(2021)
query_sub_exp = knnPredict(query, 
                           reference_mac, 
                           reference_mac$meta_data$subtype, 
                           k = 30,
                           confidence = TRUE,
                           seed = 1
)

colnames(query_sub_exp$meta_data)[(ncol(query_sub_exp$meta_data)-1):ncol(query_sub_exp$meta_data)] <- c("sub_pred_knn","sub_pred_knn_prob")

r_metadata = reference_mac$meta_data
q_metadata = query_sub_exp$meta_data
query_anno_exp_sub <- q_metadata
# Fig 4E predicted cell type of 4 plates
bar_exp <- tapply(query_anno_exp_sub$sub_pred_knn,query_anno_exp_sub$batch,table) %>% do.call(cbind,.)
bar_exp <- t(bar_exp) %>% as.data.frame() %>% mutate(sum=rowSums(.)) %>% .[,c("pre_PraM","PraM","microglia","sum")]
bar_exp <- bar_exp %>% mutate(others=bar_exp$sum-(rowSums(bar_exp[,1:3]))) %>% .[c("pre_PraM","PraM","microglia","others")] %>% t()
for(i in 1:ncol(bar_exp)){
  bar_exp[,i] <- bar_exp[,i]/sum(bar_exp[,i])
}
bar_exp <- bar_exp %>% t() %>% as.data.frame() %>% mutate(batch=rownames(.))
bar_exp <- bar_exp[,apply(bar_exp,2,function(x){max(x)!=0})]
bar_exp <- melt(bar_exp, id="batch")
bar_exp <- bar_exp %>% filter(value!=0)
bar_exp <- merge(bar_exp,meta_exp[,c("batch","gating","time")] %>% .[!duplicated(.),])
bar_exp$gating_time <- paste(bar_exp$gating,bar_exp$time,sep="_")
bar_exp$gating_time <- factor(bar_exp$gating_time,levels = c("MRC1+P2RY12+-_9 PCW","MRC1+P2RY12+-_13 PCW","MRC1-P2RY12+_9 PCW","MRC1-P2RY12+_13 PCW"))
colnames(bar_exp)[2:3] <- c("subtype","Fraction")
bar_exp$subtype <- as.character(bar_exp$subtype)
bar_exp$subtype <- factor(bar_exp$subtype,levels = c("pre_PraM","PraM","microglia","others"))
pdf(file = "./exp_composition.pdf",width =5, height = 6)
ggplot(bar_exp,aes(x=gating_time,y=Fraction,fill=subtype))+
  geom_bar(stat = "identity",width=0.8,color = "grey30")+
  ylab("Fraction of cells")+
  scale_fill_manual(values = color_bar)+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle=45, vjust=1, hjust=1,color="black"),
        axis.text.y = element_text(size=12,color="black"))

dev.off()

r_metadata$ref_query = "Ref"
q_metadata$ref_query = paste(q_metadata$gating,q_metadata$time,sep="_")

r_metadata$sub_pred_knn_prob <- 1

r_metadata <- r_metadata[,c("Well_ID","subtype","sub_pred_knn_prob","ref_query","tissue")]
q_metadata <- q_metadata[,c("Well_ID","sub_pred_knn","sub_pred_knn_prob","ref_query","tissue")]

head(r_metadata)
head(q_metadata)

colnames(r_metadata)[2:4] <- c("cell_type_pred_knn","cell_type_pred_knn_prob","ref_query")
colnames(q_metadata)[2:4] <- c("cell_type_pred_knn","cell_type_pred_knn_prob","ref_query")

meta_data_combined = rbind(r_metadata,q_metadata)
umap_combined = rbind(reference_mac$umap$embedding,query$umap)

umap_combined_labels = cbind(meta_data_combined, umap_combined) 

umap_combined_labels$con <- "Ref"


# Fig 4D MRC1-P2RY12+_9 PCW ------------------------------------------------------
umap_tmp <-umap_combined_labels %>% filter(ref_query%in%c("Ref","MRC1-P2RY12+_9 PCW"))

b1 <- which(umap_tmp$ref_query=="Ref")
b2 <- setdiff(1:nrow(umap_tmp),b1)
umap_tmp <- umap_tmp[c(b1,b2),]
umap_tmp %>% head(4)
umap_exp <- umap_tmp
umap_exp[umap_exp$ref_query=="MRC1-P2RY12+_9 PCW","con"] <- umap_exp[umap_exp$ref_query=="MRC1-P2RY12+_9 PCW","cell_type_pred_knn"]


pdf(file = "./MRC1-P2RY12+_9 PCW_mapping_sc2d.pdf",width =18, height = 18)
plotBasic2(umap_exp, title = 'MRC1-P2RY12+_9 PCW_mapping', 
           color.by = 'con',color.mapping = color_con,pt_size = 8)+
  theme(legend.position = "None")
dev.off()

# Fig S4C MRC1-P2RY12+_13 PCW ------------------------------------------------------
umap_tmp <-umap_combined_labels %>% filter(ref_query%in%c("Ref","MRC1-P2RY12+_13 PCW"))

b1 <- which(umap_tmp$ref_query=="Ref")
b2 <- setdiff(1:nrow(umap_tmp),b1)
umap_tmp <- umap_tmp[c(b1,b2),]
umap_tmp %>% head(4)
umap_exp <- umap_tmp
umap_exp[umap_exp$ref_query=="MRC1-P2RY12+_13 PCW","con"] <- umap_exp[umap_exp$ref_query=="MRC1-P2RY12+_13 PCW","cell_type_pred_knn"]


pdf(file = "./MRC1-P2RY12+_13 PCW_mapping_sc2d.pdf",width =18, height = 18)
plotBasic2(umap_exp, title = 'MRC1-P2RY12+_13 PCW_mapping',
           color.by = 'con',color.mapping = color_con,pt_size = 8)+theme(legend.position = "None")
dev.off()

# Fig 4D MRC1+P2RY12+-_9 PCW ------------------------------------------------------
umap_tmp <-umap_combined_labels %>% filter(ref_query%in%c("Ref","MRC1+P2RY12+-_9 PCW"))

b1 <- which(umap_tmp$ref_query=="Ref")
b2 <- setdiff(1:nrow(umap_tmp),b1)
umap_tmp <- umap_tmp[c(b1,b2),]
umap_tmp %>% head(4)
umap_exp <- umap_tmp
umap_exp[umap_exp$ref_query=="MRC1+P2RY12+-_9 PCW","con"] <- umap_exp[umap_exp$ref_query=="MRC1+P2RY12+-_9 PCW","cell_type_pred_knn"]


pdf(file = "./MRC1+P2RY12+-_9 PCW_mapping_sc2d.pdf",width =18, height = 18)
plotBasic2(umap_exp, title = 'MRC1+P2RY12+-_9 PCW_mapping', 
           color.by = 'con',color.mapping = color_con,pt_size = 8)+theme(legend.position = "None")
dev.off()

# Fig S4C MRC1+P2RY12+-_13 PCW ------------------------------------------------------
umap_tmp <-umap_combined_labels %>% filter(ref_query%in%c("Ref","MRC1+P2RY12+-_13 PCW"))

b1 <- which(umap_tmp$ref_query=="Ref")
b2 <- setdiff(1:nrow(umap_tmp),b1)
umap_tmp <- umap_tmp[c(b1,b2),]
umap_tmp %>% head(4)
umap_exp <- umap_tmp
umap_exp[umap_exp$ref_query=="MRC1+P2RY12+-_13 PCW","con"] <- umap_exp[umap_exp$ref_query=="MRC1+P2RY12+-_13 PCW","cell_type_pred_knn"]


pdf(file = "./MRC1+P2RY12+-_13 PCW_mapping_sc2d.pdf",width =18, height = 18)
plotBasic2(umap_exp, title = 'MRC1+P2RY12+-_13 PCW_mapping',
           color.by = 'con',color.mapping = color_con,pt_size = 8)+theme(legend.position = "None")
dev.off()

###########################################################################
# Fig 4C barplot showing fraction of macrophage subsets in 3 datasets including Popescu et al., Xu et al., and current study

Frac_fl <- bar_fl %>% mutate(prop=bar_fl$Freq/sum(bar_fl$Freq))
Frac_fs <- bar_fs %>% mutate(prop=bar_fs$Freq/sum(bar_fs$Freq))
Frac_Li <- data.frame(Var1 = metatable %>% filter(tissue=="Skin"&major=="macrophage") %>% pull(subtype) %>% table() %>% names(),
                      Freq = metatable %>% filter(tissue=="Skin"&major=="macrophage") %>% pull(subtype) %>% table() %>% as.numeric(),
                      Amp = "Li", prop = metatable %>% filter(tissue=="Skin"&major=="macrophage") %>% pull(subtype) %>% table() %>% prop.table() %>% as.numeric())
Frac_all <- rbind(Frac_fs,Frac_fl,Frac_Li)
Frac_all$Var1 <- as.character(Frac_all$Var1)
Frac_all[Frac_all$Var1=="HdpM","Var1"] <- "S100B+ACY3+pMac"
Frac_all$Var1 <- factor(Frac_all$Var1, levels = sub_order)
Frac_all <- Frac_all[order(Frac_all$Var1),]
colnames(Frac_all)[1] <- "subtype"

### show other cell types with grey
Frac_all$subtype <- as.character(Frac_all$subtype)
Frac_all[!Frac_all$subtype%in%c("pre_PraM","PraM","microglia"),"subtype"] <- "others"
tmp <- tapply(Frac_all[Frac_all$subtype=="others",]$prop,Frac_all[Frac_all$subtype=="others",]$Amp,sum) %>% 
  as.table() %>% as.data.frame() %>% mutate(subtype="others")
colnames(tmp)[1:2] <- c("Amp","prop")
Frac_all <- rbind(Frac_all[Frac_all$subtype%in%c("pre_PraM","PraM","microglia"),c("subtype","Amp","prop")], tmp[,c("subtype","Amp","prop")])
Frac_all$subtype <- factor(Frac_all$subtype, levels = c("pre_PraM","PraM","microglia","others"))
### plot bar
p1 <- ggplot(Frac_all)+
  geom_bar(aes(x=Amp, y=prop, fill=subtype),stat = "identity",width = 0.8,color="grey30")+
  scale_fill_manual(values = color_bar)+
  ylab("Cell proportion")+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.title.y = element_text(size=30,color="black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=20,color="black"),
        axis.text.x = element_text(size=20,color="black"))
ggsave(p1, filename = "Fraction_3dataset.pdf",width = 6, height = 6, device = "pdf")
