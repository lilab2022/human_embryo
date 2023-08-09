###Figs3
###writer: Feng Ruoqing & Wang Hao & Li Muxi
###Date: 2022-10-4

## Fig S3A
## The function 'DoMultiBarHeatmap2' (shown as below) was modified from https://rdrr.io/github/RuiyuRayWang/scWidgets/src/R/DoMultiBarHeatmap.R

DoMultiBarHeatmap2 <- function (object, 
                                features = NULL, 
                                cells = NULL, 
                                group.by = "ident", 
                                additional.group.by = NULL, 
                                additional.group.sort.by = NULL, 
                                cols.use = NULL,
                                group.bar = TRUE, 
                                disp.min = -2.5, 
                                disp.max = NULL, 
                                slot = "scale.data", 
                                assay = NULL, 
                                label = TRUE, 
                                size = 5.5, 
                                hjust = 0, 
                                angle = 45, 
                                raster = TRUE, 
                                draw.lines = TRUE, 
                                lines.width = NULL, 
                                group.bar.height = 0.02, 
                                combine = TRUE) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  
  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ", 
                paste(bad.sorts, collapse = ", "))
      }
    }
  }
  
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is_null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]  
      if (!is_null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]  
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }
    
    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]
    
    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }  
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }
      
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }
      
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    
    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])
    
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }
        
        # Default
        cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))  
        
        #Overwrite if better value is provided
        if (!is_null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is_null(names(cols.use[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }
        
        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)  
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(plot + 
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off")) 
        
        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x.major
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                    label.x.pos)
          plot <- plot + geom_text(stat = "identity", 
                                   data = label.x.pos, aes_string(label = "group", 
                                                                  x = "label.x.pos"), y = y.max + y.max * 
                                     0.03 * 0.5, angle = angle, hjust = hjust, 
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) * 
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}



load("G://embryo/20220929/metatable_0929.Rdata")

###
meta_addition<-read.csv("G://embryo/20220930/subtypes_mf_20220930.csv",header=T)
load("G://embryo/20220929/mac.color2.rdata")

data<-read.csv("G://embryo/20220329/cleanumi.csv",header=T,row.names = "X") # loading mc level umi-tab 

data<-data[,colnames(data)%in%meta_addition$MC]

###
excluded.genes<-read.csv("G://embryo/20220616/data/excludeall_gene.csv",header=T) ## loading forbidden genes
data<-data[!rownames(data)%in%excluded.genes$unique.exclude_all.,]

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
                          resolution = 0.3, 
                          random.seed = 1) 

meta<-(seob_data@meta.data)
meta<-cbind(rownames(meta),meta)
colnames(meta)[1]<-"MC"

meta<-merge(meta,meta_addition,by="MC")
rownames(meta)<-meta$MC

seob_data@meta.data<-meta

Idents(seob_data)<-"clustering"

mf.markers2 <- FindAllMarkers(seob_data, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

top10 <- mf.markers2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seob_data, features = top10$gene) + NoLegend()

seob_data$clustering<-factor(x=seob_data$clustering,levels=c("Langerhans","group1","group2","group3","Microglia","Osteoclast"))
seob_data$subtype<-factor(x=seob_data$subtype,levels=c("Langerhans","ACY3+S100B+pMφ","HepM","YsdM_AFP_low","YsdM_AFP_high","Kupffer_cell","red_pulp","pre_PraM","gonad_macrophage","PraM",
                                                       "Adrenalgland_macrophage","intestine CD207+ Mφ","intestine CD209+ Mφ","microglia","osteoclast"))

# Show we can color anything
cols.use <- list(clustering=c('white','white', 'white', 'white', 'white','white'),
                 subtype=c('#DBDB8DFF', '#74350A','#FFBB78FF','#E377C2FF','#9467BDFF', '#FF9896FF','#F7B6D2FF','#2CA02CFF','#17BECFFF','#FF7F0EFF','#BCBD22FF',
                           '#98DF8AFF','#AEC7E8FF','#1F77B4FF','#C49C94FF'))  

####### gene reorder

new_order<-c("Langerhans","group1","group2","group3","Microglia","Osteoclast")

template<-matrix(,nrow=1,ncol=7)
colnames(template)<-colnames(top10)
template<-as.data.frame(template)

for(i in 1:length(unique(top10$cluster))){
  
  
  new_i<-new_order[i]
  
  new_row<-subset(top10,cluster==new_i)
  
  template<-rbind(template,new_row)
  
}

template<-as.data.frame(template)

template<-template[-1,]

write.csv(template,"G://embryo/20220930/gene_reorder_fig3_suppl.csv")


p<-DoMultiBarHeatmap2(seob_data, assay = 'RNA', features = template$gene, group.by='clustering', additional.group.by = 'subtype',additional.group.sort.by = 'subtype',label = FALSE,cols.use=cols.use)+
  scale_fill_gradientn(colors=c("#1B519C","white","#A00021")) 

ggsave(p,filename = "G://embryo/figures/figS3/figS3A.png",width = 20, height = 10, units = "in", device='png')
ggsave(p,filename = "G://embryo/figures/figS3/figS3A.pdf",width = 20, height = 10, units = "in", device='pdf')

#Figure S3B and S3C

load("G://embryo/20220929/metatable_0929.Rdata")
load("G://embryo/20220929/mac.color2.rdata")


simple.meta<-metatable[,c(1,2,5,6,7,10,11,12)]


differentiation<-read.csv("G://embryo/20220808/mf_prog_cytotrace.csv")

colnames(differentiation)<-c("Well_ID","results.CytoTRACE")
differentiation$Cell<-paste0(differentiation$Well_ID,"_human_development")

differentiation<-merge(simple.meta,differentiation,by="Cell")

#remove adult
differentiation<-differentiation[!differentiation$time%in%"Adult",]

#only mf
differentiation<-differentiation[differentiation$major%in%"macrophage",]


p1<-ggplot(differentiation,aes(x = reorder(subtype, -results.CytoTRACE, FUN=median), y = results.CytoTRACE,fill = subtype))+
  scale_fill_manual(values = as.character(mac.color),breaks = names(mac.color))+
  geom_boxplot(color = "grey88",outlier.shape = NA)+
  xlab("Subtype")+ylab("Differentiation capacity")+
  theme_bw()+theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,face="bold"),legend.position = 'none',
        axis.title.y = element_text(face="bold"),
        axis.title.x = element_text(face="bold")
  )
ggsave("G://embryo/figures/figS3/fig3_differentiation_boxplot2.png",p1,width = 12,height = 6,dpi = 300,units="in",device="png")
ggsave("G://embryo/figures/figS3/fig3_differentiation_boxplot2.pdf",p1,width = 12,height = 6,dpi = 300,units="in",device="pdf")

# Figure 3B Upper

meta_addition<-read.csv("G://embryo/20220930/subtypes_mf_20220930.csv",header=T)
meta_addition<-meta_addition[,c(1,3,4,7)]
colnames(meta_addition)[1]<-"mc"
differentiation<-merge(differentiation,meta_addition,by="mc")

group.color<-c("black","#FFD700","darkred","#1F77B4FF","#C49C94FF","#DBDB8DFF")
names(group.color)<-c("group1","group2","group3","Microglia","Osteoclast","Langerhans")

differentiation.x<-differentiation[!differentiation$clustering%in%c("Microglia","Osteoclast","Langerhans"),]

p3<-ggplot(differentiation.x,aes(x = reorder(clustering, -results.CytoTRACE, FUN=median), y = results.CytoTRACE,fill = clustering))+
  scale_fill_manual(values = as.character(group.color),breaks = names(group.color))+
  geom_boxplot(color = "grey88",outlier.shape = NA,width=0.5)+
  xlab("")+ylab("differentiation")+
  theme_bw()+theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,face="bold"),legend.position = 'none',
        axis.title.y = element_text(face="bold"),
        axis.title.x = element_text(face="bold")
  )+
  geom_signif(comparisons = list(c("group1","group3"),
                                 c("group1","group2"),
                                 c("group2","group3")),
              map_signif_level = T,
              step_increase=0.2
  )
ggsave("G://embryo/figures/figS3/figs3_differentiation_boxplot_byGroup.png",p3,width = 12,height = 6,dpi = 300,units="in",device="png")
ggsave("G://embryo/figures/figS3/figs3_differentiation_boxplot_byGroup.pdf",p3,width = 12,height = 6,dpi = 300,units="in",device="pdf")


#### proliferation

proliferation<-read.table("G://embryo/20220814/CellcycleScoring.csv")
proliferation<-proliferation[,c("Well_ID","G2M.Score")]
proliferation$Cell<-paste0(proliferation$Well_ID,"_human_development")

proliferation<-merge(simple.meta,proliferation,by="Cell")
proliferation<-proliferation[,-9]
colnames(proliferation)[3]<-"Well_ID"

#remove adult
proliferation<-proliferation[!proliferation$time%in%"Adult",]

#only mf
proliferation<-proliferation[proliferation$major%in%"macrophage",]


p1<-ggplot(proliferation,aes(x = reorder(subtype, -G2M.Score, FUN=median), y = G2M.Score,fill = subtype))+
  scale_fill_manual(values = as.character(mac.color),breaks = names(mac.color))+
  geom_boxplot(color = "grey88",outlier.shape = NA)+
  xlab("Subtype")+ylab("G2M.score")+
  theme_bw()+theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,face="bold"),legend.position = 'none',
        axis.title.y = element_text(face="bold"),
        axis.title.x = element_text(face="bold")
  )
ggsave("G://embryo/figures/figS3/figS3_proliferation_boxplot.png",p1,width = 12,height = 6,dpi = 300,units="in",device="png")
ggsave("G://embryo/figures/figS3/figS3_proliferation_boxplot.pdf",p1,width = 12,height = 6,dpi = 300,units="in",device="pdf")


# Figure 3B Lower

meta_addition<-read.csv("G://embryo/20220803/subtypes_mf.csv",header=T)
meta_addition<-meta_addition[,c(1,3,4,7)]
colnames(meta_addition)[1]<-"mc"
proliferation<-merge(proliferation,meta_addition,by="mc")

group.color<-c("black","#FFD700","darkred","#1F77B4FF","#C49C94FF","#DBDB8DFF")
names(group.color)<-c("group1","group2","group3","Microglia","Osteoclast","Langerhans")

proliferation.x<-proliferation[!proliferation$clustering%in%c("Microglia","Osteoclast","Langerhans"),]

p3<-ggplot(proliferation.x,aes(x = reorder(clustering, -G2M.Score, FUN=median), y = G2M.Score,fill = clustering))+
  scale_fill_manual(values = as.character(group.color),breaks = names(group.color))+
  geom_boxplot(color = "grey88",outlier.shape = NA,width=0.5)+
  xlab("")+ylab("G2M.score")+
  theme_bw()+theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,face="bold"),legend.position = 'none',
        axis.title.y = element_text(face="bold"),
        axis.title.x = element_text(face="bold")
  )+
  geom_signif(comparisons = list(c("group1","group3"),
                                 c("group1","group2"),
                                 c("group2","group3")),
              map_signif_level = T,
              step_increase=0.2
  )
ggsave("G://embryo/figures/figS3/fig3_proliferation_boxplot_byGroup.png",p3,width = 12,height = 6,dpi = 300,units="in",device="png")
ggsave("G://embryo/figures/figS3/fig3_proliferation_boxplot_byGroup.pdf",p3,width = 12,height = 6,dpi = 300,units="in",device="pdf")

# Figure S3D


rm(list = ls())


####################################################

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

scalar_plot<-function(data_prep2,A,B,A_1,B_1,genelist,A_color,B_color){
  
  
  
  data_prep2$label = ""
  data_prep2$gene = rownames(data_prep2)
  data_prep2$label[data_prep2$gene%in%genelist]=data_prep2$gene[data_prep2$gene%in%genelist]
  
  
  
  corr<-cor(data_prep2[,1],data_prep2[,2])
  corr<-round(corr,digits=4)
  
  colnames(data_prep2)[1:2] = c(A,B)
  p<-ggplot(data_prep2,aes(x = get(A),y = get(B),label = label))+
    geom_point(aes(color = gene_tag), size = 1)+xlim(0,10)+ylim(0,10)+
    geom_text_repel(data = data_prep2[data_prep2 != "",],max.overlaps = 100000)+theme_bw() + #閼冲本娅欑拫鍐╂殻
    
    labs(x = A_1, y = B_1, color = '',subtitle = paste("Correlation:",corr,sep="")) + #閸ф劖鐖ｆ潪瀛樼垼妫版顔曠純?
    
    geom_abline(intercept = 0.585, slope = 1, col = 'black', linetype = 'dashed', size = 0.5) + #鏉??3閸欍儳鏁ゆ禍搴㈠潑閸?? |log2FC|>1 閻ㄥ嫰妲囬崐鑲╁殠
    
    geom_abline(intercept = -0.585, slope = 1, col = 'black', linetype = 'dashed', size = 0.5) +
    
    geom_abline(intercept = 0, slope = 1, col = 'black', linetype = 'dashed', size = 0.8)+
    
    scale_color_manual(values = c(A_color, 'gray', B_color)) +theme(legend.position = 'none')+theme(text = element_text(size = 20))  
  
  ggsave(p,filename = paste("G://embryo/figures/figS3/",A,"_VS_",B,"_scalar_plot.png",sep=""),width = 5, height = 5, units = "in", device='png')
  ggsave(p,filename = paste("G://embryo/figures/figS3/",A,"_VS_",B,"_scalar_plot.pdf",sep=""),width = 5, height = 5, units = "in", device='pdf')
  
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

scalar_plot2<-function(data_prep2,A,B,A_1,B_1,genelist,A_color,B_color){
  
  
  
  data_prep2$label = ""
  data_prep2$gene = rownames(data_prep2)
  data_prep2$label[data_prep2$gene%in%genelist]=data_prep2$gene[data_prep2$gene%in%genelist]
  
  
  
  corr<-cor(data_prep2[,1],data_prep2[,2])
  corr<-round(corr,digits=4)
  
  colnames(data_prep2)[1:2] = c(A,B)
  p<-ggplot(data_prep2,aes(x = get(A),y = get(B),label = label))+
    geom_point(aes(color = gene_tag), size = 1)+xlim(0,10)+ylim(0,10)+
    geom_text_repel(data = data_prep2[data_prep2 != "",],max.overlaps = 100000)+theme_bw() + #閼冲本娅欑拫鍐╂殻
    
    labs(x = A_1, y = B_1, color = '',subtitle = paste("Correlation:",corr,sep="")) + #閸ф劖鐖ｆ潪瀛樼垼妫版顔曠純?
    
    geom_abline(intercept = 0.585, slope = 1, col = 'black', linetype = 'dashed', size = 0.5) + #鏉??3閸欍儳鏁ゆ禍搴㈠潑閸?? |log2FC|>1 閻ㄥ嫰妲囬崐鑲╁殠
    
    geom_abline(intercept = -0.585, slope = 1, col = 'black', linetype = 'dashed', size = 0.5) +
    
    geom_abline(intercept = 0, slope = 1, col = 'black', linetype = 'dashed', size = 0.8)+
    
    scale_color_manual(values = c(A_color, 'gray', B_color)) +theme(legend.position = 'none')+theme(text = element_text(size = 20))  
  
  ggsave(p,filename = paste("G://embryo/figures/figS3/",A,"_VS_",B,"_scalar_plot.png",sep=""),width = 8, height = 8, units = "in", device='png')
  ggsave(p,filename = paste("G://embryo/figures/figS3/",A,"_VS_",B,"_scalar_plot.pdf",sep=""),width = 8, height = 8, units = "in", device='pdf')
  
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

load("G://embryo/20220420/embryo_raw_umi.RData") # loading umi-tab
load("G://embryo/20220929/metatable_0929.Rdata")

#gene_counts<-gene_counts[rownames(gene_counts)%in%excluded_genes,]

mf_id<-read.csv("G://embryo/20220930/subtypes_mf_20220930.csv",header = T) # loading  macrophage meta information

mf_new_mc<-metatable[metatable$mc%in%mf_id$MC,]

excluded.genes<-discard_gene2(rownames(a))
a_reduced<-a[rownames(a)%in%excluded.genes,]


############global normalization

seob_data <- CreateSeuratObject(counts = a_reduced)
seob_data <- NormalizeData(seob_data, normalization.method =  "RC") # without log transformation

a_reduced <- seob_data@assays$RNA@data

### remove adult
adult.id<-subset(metatable,time=="Adult")
a_reduced<-a_reduced[,!colnames(a_reduced)%in%adult.id]

# the name belwo must be consistent with colnames of data
A<-"Kupffer_cell"
B<-"red_pulp"

A_color<-'#F7B6D2FF'
B_color<-'#FF9896FF'

FC_cutoff<-"0.585"
genelist<-c("MARCO","F13A1","FCGR3A","SPP1","RNASE1","CD5L")

data_prep2<-prep(A,B,a_reduced,FC_cutoff,mf_new_mc)

A_1<-"Kupffer"
B_1<-"Red pulp"

scalar_plot2(data_prep2,A,B,A_1,B_1,genelist,A_color,B_color)
#scalar_plot(data_prep2,A,B,A_1,B_1,genelist,A_color,B_color)


A<-"gonad_macrophage"
B<-"Adrenalgland_macrophage"

B_color<-'#17BECFFF'
A_color<-'#BCBD22FF'

FC_cutoff<-"0.585"
genelist<-c("MMP9","SPP1","FCGR3A","LYVE1")

data_prep2<-prep(A,B,a_reduced,FC_cutoff,mf_new_mc)

A_1<-"Male gonad Mφ"
B_1<-"AG Mφ"


scalar_plot2(data_prep2,A,B,A_1,B_1,genelist,A_color,B_color)
#scalar_plot(data_prep2,A,B,A_1,B_1,genelist,A_color,B_color)

###

A<-"gonad_macrophage"
B<-"Adrenalgland_macrophage"

B_color<-'#17BECFFF'
A_color<-'#BCBD22FF'

FC_cutoff<-"0.585"
genelist<-c("MMP9","SPP1","FCGR3A","LYVE1","AXL")

data_prep2<-prep(A,B,a_reduced,FC_cutoff,mf_new_mc)

A_1<-"Male gonad Mφ"
B_1<-"AG Mφ"


scalar_plot2(data_prep2,A,B,A_1,B_1,genelist,A_color,B_color)
#scalar_plot(data_prep2,A,B,A_1,B_1,genelist,A_color,B_color)

### gonad 

osteoclast_all<-subset(metatable,subtype=="osteoclast")

table(osteoclast_all$tissue)

osteoclast_all_bonemarrow_gonad<-osteoclast_all[osteoclast_all$tissue%in%c("Male gonad","Bonemarrow"),]


a_reducedx<-a_reduced[,colnames(a_reduced)%in%osteoclast_all_bonemarrow_gonad$Well_ID]
gonad.id<-osteoclast_all_bonemarrow_gonad[osteoclast_all_bonemarrow_gonad$tissue%in%"Male gonad",]

bonemarrow.id<-osteoclast_all_bonemarrow_gonad[osteoclast_all_bonemarrow_gonad$tissue%in%"Bonemarrow",]

data_A<-a_reducedx[,colnames(a_reducedx)%in%gonad.id$Well_ID]
data_B<-a_reducedx[,colnames(a_reducedx)%in%bonemarrow.id$Well_ID]


A_row_mean = data_calc_row_mean(data_A,"osteoclast_in_gonad")
B_row_mean = data_calc_row_mean(data_B,"osteoclast_in_bonemarrow")

data_prep <-merge(A_row_mean,B_row_mean,by='row.names',all=TRUE)
rownames(data_prep)<-data_prep$Row.names
data_prep<-data_prep[,-1]


data_prep<-log2(data_prep+1)

log2FC = (data_prep[,1])-(data_prep[,2])

data_prep_FC<-cbind(data_prep,log2FC)



gene_tag<-tag(log2FC,as.numeric(0.585))

data_prep2<-cbind(data_prep_FC,gene_tag)

A<-"osteoclast_in_gonad"
B<-"osteoclast_in_bonemarrow"

B_color<-'#FFD700'
A_color<-'#000000'

A_1<-"osteoclast in gonad"
B_1<-"osteoclast in bonemarrow"

genelist<-c("MMP9","SPP1","SELENOP","HLA-DRA","S100A11","ITGB2","IL7R")


scalar_plot2(data_prep2,A,B,A_1,B_1,genelist,A_color,B_color)

## intestine mf 209 vs 207

A<-"intestine CD209+ Mφ"
B<-"intestine CD207+ Mφ"

A_color<-'#AEC7E8FF'
B_color<-'#98DF8AFF'

FC_cutoff<-"0.585"
genelist<-c("CD209","CD207","SELENOP", "ADAMDEC1")

data_prep2<-prep(A,B,a_reduced,FC_cutoff,mf_new_mc)

A_1<-"intest. CD209+ Mφ"
B_1<-"intest. CD207+ Mφ"

scalar_plot2(data_prep2,A,B,A_1,B_1,genelist,A_color,B_color)

#Figure 3s F
### S3F mac composition for each tissue (thymus removed)
meta_S3F <- tapply(metatable$subtype,metatable$tissue,function(x){prop.table(table(x))})
meta_S3F <- meta_S3F[-which(names(meta_S3F)=="Thymus")]
meta_S3F <- do.call(cbind, meta_S3F) %>% as.data.frame() %>% mutate(subtype=rownames(.)) %>% melt()
colnames(meta_S3F)[2:3] <- c("tissue","Freq") 
meta_S3F$subtype <- factor(meta_S3F$subtype, levels = sub_order)
tissue_order <- c("Yolksac","Embryo","Head","AGM","Limb","Heart","Liver","Spleen","Female gonad",
                  "Male gonad", "Bonemarrow","Gastrointestinal tract","Skin","Brain",
                  "Spinalmarrow","Adrenalgland","Kidney","Lung")
tissue.color <- tissue.color[tissue_order]
meta_S3F$tissue <- factor(meta_S3F$tissue, levels = tissue_order)
meta_S3F <- meta_S3F[order(meta_S3F$tissue,meta_S3F$subtype),]

p <- ggplot(meta_S3F)+
  geom_bar(aes(x= tissue,y=Freq,fill=subtype),stat="identity",position = "fill",width = 0.8,color = "grey30")+
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
  scale_fill_manual(values = mac.color,breaks = names(mac.color))

w <- length(unique(meta_S3F$tissue))*0.4
dir.path <- "./S3F/"
ggsave(paste(dir.path ,"mac_tissue_subtype_dynamics.pdf",sep = ""),p,width = w,height = 6)

### S3F tissue composition for each mac subtype 
tmp <- metatable %>% filter(tissue!="Thymus")
meta_S3F_2 <- tapply(tmp$tissue,tmp$subtype,function(x){data.frame(tissue=names(prop.table(table(x))),Freq=as.numeric(prop.table(table(x))))})
for(i in names(meta_S3F_2)){meta_S3F_2[[i]] <- mutate(meta_S3F_2[[i]],subtype=i)}
Freq_complete <- function(Freq){
  for(i in names(Freq)){
    tissue_full <- unique(metatab
