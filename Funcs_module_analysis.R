### Functions for embryonic macrophage module analysis

##############################################################################
### read in forbidden gene list 
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

##############################################################################
### define expanded module geneset
# expand module on anchors
Fmodule_sel <- function(fp,gene,sub){
  ls_all <- list()
  ls_sel <- list()
  fp <- fp[gf,]
  for(i in gene){
    cor <- tail(sort(cor(fp[i, ], t(fp))[1,]), 100)
    ls_sel[[i]] <- data.frame(gene=names(cor),corr=cor)
  }
  df_sel <- do.call(rbind,ls_sel)
  df_sel <- df_sel%>%filter(gene%in%names(which(table(df_sel$gene)==5)))
  top_gene <- tapply(df_sel$corr,df_sel$gene,mean)%>%sort()%>%tail(70)%>%sort(decreasing = T)
  top_gene <- top_gene[which(top_gene>0.5)]
  top_gene <- data.frame(gene=names(top_gene),corr=top_gene)%>%na.omit()%>%t()%>%as.data.frame()
  colnames(top_gene) <- 1:ncol(top_gene)
  # add module information
  top_gene[nrow(top_gene)+1,] <- sub
  rownames(top_gene) <- c(sub,paste(sub,"corr",sep="_"),paste(sub,"module",sep="_"))
  
  return(top_gene)
}

# select anchors
FfocG <- function(fp,module_order,subtype=NULL,module){
  if(is.null(subtype)){
    gene_foc <- module_order%>%filter(category%in%module)%>%pull(gene)
    fp_s <- fp_s[gene_foc,] %>% as.matrix()
    gene_num <- rowMax(fp_s)
    names(gene_num) <- rownames(fp_s)
    gene_num <- sort(gene_num,decreasing = T)
    fp_s <- names(gene_num)[1:5]
  }
  else{
    fp_s <- (fp[,subtype]/apply(fp[,which(!colnames(fp)%in%subtype)],1,max))%>%sort(decreasing = T)%>%names()
    gene_foc <- module_order%>%filter(category%in%module)%>%pull(gene)
    fp_s <- intersect(fp_s, gene_foc)[1:5]
  }
  return(fp_s)
}

# define all modules
df_module_big <- function(ls_module,module_order){
  module_sel <- do.call(function(...) {
    tmp <- plyr::rbind.fill(...)
    rownames(tmp) <- sapply(ls_module, function(i){
      rownames(i)
    })
    return(tmp)
  }
  ,ls_module)%>%t()%>%as.data.frame()
  gene <- NULL
  for(i in 1:(ncol(module_sel)/3)){
    gene <- c(gene,na.omit(as.character(module_sel[,1+(i-1)*3])))
  }
  print(table(gene)[table(gene)>1])
  gene_test <- names(table(gene)[table(gene)>1])
  print(gene_test%in%md_pro_gut[1,])
  print(gene_test%in%md_pro_MHC[1,])
  ### 
  gene_inter <- names(which(table(gene)>1))
  gene_co <- data.frame(gene=module_sel$gut,gut_corr=module_sel$gut_corr)%>%filter(gene%in%gene_inter)
  gene_co <- merge(gene_co, data.frame(gene=module_sel$`MHC-II`,MHCII_corr=module_sel$`MHC-II_corr`))
  gene_co$final <- "gut"
  gene_co[gene_co$gut_corr<gene_co$MHCII_corr,"final"] <- "MHC-II"
  gene_co_gut <- gene_co%>%filter(final=="gut")%>%pull(gene)
  gene_co_ap <- gene_co%>%filter(final=="MHC-II")%>%pull(gene)
  module_sel[module_sel$gut%in%gene_co_ap,"gut_corr"] <- NA
  module_sel[module_sel$`MHC-II`%in%gene_co_gut,"MHC-II_corr"] <- NA
  module_sel[module_sel$gut%in%gene_co_ap,"gut"] <- NA
  module_sel[module_sel$`MHC-II`%in%gene_co_gut,"MHC-II"] <- NA
    ls_sel <- list()
  for(i in 1:(ncol(module_sel)/3)){
    ls_sel[[i]] <- module_sel[,(3*(i-1)+1):(3*(i-1)+3)]
    colnames(ls_sel[[i]]) <- c("gene","corr","category")
  }
  names(ls_sel) <- c("microglia","core-mf","MHC-II","gut","kupffer","PraM","langerhans","osteoclast")
  ls_sel <- lapply(ls_sel,na.omit)
  for( i in names(ls_sel)){
    if(i=="PraM"){
      ls_sel[[i]] <- ls_sel[[i]][1:40,]
    }else{
      ls_sel[[i]] <- ls_sel[[i]][1:20,]
    }
  }
  ls_sel <- do.call(rbind,ls_sel)%>%na.omit()
  ls_sel <- merge(ls_sel[,c("gene","category")],module_order[!duplicated(module_order[,c("Module","category")]),c("Module","category")])
  module_order_new <- rbind(ls_sel[,c("gene","Module","category")],module_order[module_order$category%in%c("un"),c("gene","Module","category")])
  module_order_new <- module_order_new[!duplicated(module_order_new$gene),]
  module_order_new <- module_order_new[order(module_order_new$Module),]
  module_order_new$module_cat <- paste(module_order_new$Module,module_order_new$category,sep="_")
  return(module_order_new)
  return(ls_sel)
}

##############################################################################

# compute max for module score
max_md <- function(umi,sub){
  umi <- umi %>% filter(subtype%in%sub)
  print(paste("MHC-II",ceiling(max(umi$`MHC-II`)),sep=" "))
  print(paste("core-mf",ceiling(max(umi$`core-mf`)),sep=" "))
}

##############################################################################

# save list object for pdf
pdf_save <- function(list,path,width=6,height=6){
  if(!dir.exists(path)) dir.create(path)
  for (i in names(list)){
    pdf(paste(path,i,".pdf",sep=""),width = width,height = height)
    print(list[[i]])
    dev.off()
  }
}

##############################################################################

# Fig 7B plot subtype quintile for multi scores (except pre-PraM & PraM for Proangiogenic module)
plt_quintile_v1 <- function(umi_pct,module,sub){
  
  ls_plt <- list()
  df <- umi_pct%>%filter(subtype==sub)
  df <- df[order(df[,module]),]
  rownames(df) <- 1:nrow(df)
  df$quintile <- factor(sort(rank(row.names(df))%%5)+1)
  
  mhc2_range <- tapply(df$`MHC-II`,df$quintile,summary) %>% do.call(rbind,.) %>% as.data.frame() %>% mutate(quintile=rownames(.))
  core_mf_range <- tapply(df$`core-mf`,df$quintile,summary) %>% do.call(rbind,.) %>% as.data.frame() %>% mutate(quintile=rownames(.))
  md_range <- tapply(df[,module],df$quintile,summary) %>% do.call(rbind,.) %>% as.data.frame() %>% mutate(quintile=rownames(.))
  
  df1 <- tapply(df$Phase, df$quintile, function(x){length(which(x=="G2M"))})
  df2 <- tapply(df$Phase, df$quintile, length)
  phase <- data.frame(G2M_pct=df1/df2*100,`module_score_quintile`=as.factor(1:5),total=df2)
  
  time1 <- do.call(cbind,tapply(df$time, df$quintile, function(x){table(x) %>% prop.table()})) %>% as.data.frame()
  time_info <- t(time1) %>% as.data.frame()
  time_info <- time_info %>% mutate(quintile=rownames(time_info)) %>% melt(id="quintile")
  time_info$quintile <- as.numeric(time_info$quintile) %>% as.factor()
  colnames(time_info) <- c("module_score_quintile","time","time_distribution")
  time_info$time3 <- time_info$time

  time_info$time3 <- gsub("cs.*","CS_stage",time_info$time3)
  time_info[time_info$time3%in%paste("week",9:13,sep=""),"time3"] <- "9-13 PCW"
  time_info[!time_info$time3%in%c("CS_stage","9-13 PCW"),"time3"] <- "After 16 PCW"
  time_info$time3 <- factor(time_info$time3,levels = c("CS_stage","9-13 PCW","After 16 PCW"))
  
  ### MHC-II
  ls_plt[[paste("MHC-II",sub,sep="_")]] <- ggplot(df,aes(x=quintile,y=`MHC-II`))+
    geom_violin(fill=mac.color[sub],color=mac.color[sub],size=0.8,width=0.8)+
    geom_pointrange(data = mhc2_range,
                    mapping = aes(x = quintile, y = Median, ymin=`1st Qu.`,ymax=`3rd Qu.`),
                    color="purple",size=0.2,fill="purple") +
    geom_line(data = mhc2_range,
              mapping = aes(x = quintile, y = Median,group=1),size=2,color="purple")+
    theme(legend.position = "None")+
    xlab(paste(module,"score_quintile",sep="_"))+
    ylim(0,10)+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "None",
          axis.text = element_text(size=30),
          title = element_text(size=30),
          axis.line.y = element_line(color = ifelse(sub=="microglia","black",NA)))

  ## core-mf expression
  ls_plt[[paste("core-mf",sub,sep="_")]] <- ggplot(df,aes(x=quintile,y=`core-mf`))+
    geom_violin(fill=mac.color[sub],color=mac.color[sub],size=0.8,width=0.8)+
    geom_pointrange(data = core_mf_range,color="red",size=0.05,fill="red",
                    mapping = aes(x = quintile, y = Median, ymin=`1st Qu.`,ymax=`3rd Qu.`))+
    geom_line(data = core_mf_range,
              mapping = aes(x = quintile, y = Median,group=1),size=0.3,color="red")+
    xlab(paste(module,"score_quintile",sep="_"))+
    ggtitle(sub)+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "None",
          axis.text = element_text(size=30),
          title = element_text(size=30),
          axis.line.y = element_line(color = ifelse(sub=="microglia","black",NA)),
    )+
    ylim(0,6)
  
  ### module expression
  scale=max(phase$`G2M_pct`)/max(df[,module])
  
  ls_plt[[paste(module,sub,sep="_")]] <- ggplot(df)+
    geom_violin(aes(x=quintile,y=df[,module]),fill="#A6A8A9FF",color="#A6A8A9FF",size=0.8,width=0.8)+
    geom_pointrange(data = md_range,color="red",size=0.3,fill="red",
                    mapping = aes(x = quintile, y = Median, ymin=`1st Qu.`,ymax=`3rd Qu.`)) +
    geom_line(data = md_range,size=2,color="red",mapping = aes(x = quintile, y = Median,group=1))+
    geom_point(data=phase,aes(x=`module_score_quintile`,y=`G2M_pct`/scale),shape=21,size=4,color="black",fill="#219ebc")+
    geom_line(data=phase,aes(x=`module_score_quintile`,y=`G2M_pct`/scale),group=1,color="#219ebc",size=2)+
    scale_y_continuous(sec.axis = sec_axis(~.*scale, name="% proliferating cells")) +
    xlab(paste(module,"score_quintile",sep="_"))+ylab("module score")+
    ggtitle(sub)+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.text.y.left = element_text(size=30,colour = "black"),
          axis.text.y.right = element_text(size=30,colour = "black"),
          axis.line.y.left = element_line(colour = "red",size=2),
          axis.line.y.right = element_line(colour = "#219ebc",size=2),
          legend.position = "None",
          axis.text = element_text(size=30),
          title = element_text(size=30))  
  
  ### developmental time fraction
  ls_plt[[paste("time",sub,sep="_")]] <- ggplot(time_info,aes(x=`module_score_quintile`,y=`time_distribution`))+
    geom_bar(stat="identity",aes(fill=time3),width=0.6)+
    scale_fill_manual(values = c("#2a9d8f","#e9c46a","#e76f51"))+
    xlab(paste(module,"score_quintile",sep="_"))+
    scale_y_continuous(expand=c(0,0))+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "None",
          axis.text = element_text(size=30),
          title = element_text(size=30),
          axis.line.y = element_line(color = ifelse(sub=="microglia","black",NA))
    )
  
  pdf_save(ls_plt,"~/metacell/Fig7/figs/7C/v1/")
  
  return(ls_plt)
}

##############################################################################

# Fig 7B plot subtype quintile for multi scores with significance level for langerhans and CD209+ mf

plt_quintile_signif <- function(umi_pct,module,sub){
  ls_plt <- list()
  df <- umi_pct%>%filter(subtype==sub)
  df <- df[order(df[,module]),]
  rownames(df) <- 1:nrow(df)
  df$quintile <- factor(sort(rank(row.names(df))%%5)+1)
  mhc2_range <- tapply(df$`MHC-II`,df$quintile,summary) %>% do.call(rbind,.) %>% as.data.frame() %>% mutate(quintile=rownames(.))
  ### MHC-II
  ls_plt[[paste("MHC-II_signif",sub,sep="_")]] <- ggplot(df,aes(x=quintile,y=`MHC-II`))+
    geom_violin(fill=mac.color[sub],color=mac.color[sub],size=0.8,width=0.8)+
    geom_pointrange(data = mhc2_range,
                    mapping = aes(x = quintile, y = Median, ymin=`1st Qu.`,ymax=`3rd Qu.`),
                    color="purple",size=0.2,fill="purple") +
    geom_line(data = mhc2_range,
              mapping = aes(x = quintile, y = Median,group=1),size=2,color="purple")+
    theme(legend.position = "None")+
    xlab(paste(module,"score_quintile",sep="_"))+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "None",
          axis.text = element_text(size=30),
          title = element_text(size=30),
          axis.line.y = element_line(color = ifelse(sub=="microglia","black",NA)))
  if(sub=="langerhans"){
    ls_plt[[paste("MHC-II_signif",sub,sep="_")]] <- ls_plt[[paste("MHC-II_signif",sub,sep="_")]]+
      geom_signif(comparisons = list(c(1,2),c(1,3),c(1,4),c(1,5)),test = "wilcox.test",step_increase = 0.15,map_signif_level=T,
                  textsize=10,size=2,tip_length = 0,vjust=0.5)
  }else{
    ls_plt[[paste("MHC-II_signif",sub,sep="_")]] <- ls_plt[[paste("MHC-II_signif",sub,sep="_")]]+
      geom_signif(comparisons = list(c(1,2),c(2,3),c(3,4),c(4,5)),test = "wilcox.test",step_increase = 0.15,map_signif_level=T,
                  textsize=10,size=2,tip_length = 0,vjust=0.5)
  }
  
  pdf_save(ls_plt,"~/metacell/Fig7/figs/S7/",height = 8)
  
  return(ls_plt)
}

##############################################################################

## Fig 7B plot subtype quintile for multi scores (pre-PraM & PraM for Proangiogenic module)
plt_quintile_PraM <- function(umi_pct,module,sub){
  ls_plt <- list()
  mac.color <- c(mac.color,"pre_PraM & PraM"="#A8DADC")
  df <- umi_pct%>%filter(subtype==sub)
  df <- df[order(df[,module]),]
  rownames(df) <- 1:nrow(df)
  df$quintile <- factor(sort(rank(row.names(df))%%5)+1)
  mhc2_range <- tapply(df$`MHC-II`,df$quintile,summary) %>% do.call(rbind,.) %>% as.data.frame() %>% mutate(quintile=rownames(.))
  core_mf_range <- tapply(df$`core-mf`,df$quintile,summary) %>% do.call(rbind,.) %>% as.data.frame() %>% mutate(quintile=rownames(.))
  md_range <- tapply(df[,module],df$quintile,summary) %>% do.call(rbind,.) %>% as.data.frame() %>% mutate(quintile=rownames(.))
  
  
  df1 <- tapply(df$Phase, df$quintile, function(x){length(which(x=="G2M"))})
  df2 <- tapply(df$Phase, df$quintile, length)
  phase <- data.frame(G2M_pct=df1/df2*100,`module_score_quintile`=as.factor(1:5),total=df2)
  
  time1 <- do.call(cbind,tapply(df$time, df$quintile, function(x){table(x) %>% prop.table()})) %>% as.data.frame()
  time_info <- t(time1) %>% as.data.frame()
  time_info <- time_info %>% mutate(quintile=rownames(time_info)) %>% melt(id="quintile")
  time_info$quintile <- as.numeric(time_info$quintile) %>% as.factor()
  colnames(time_info) <- c("module_score_quintile","time","time_distribution")
  time_info$time3 <- time_info$time
  time_info$time3 <- gsub("cs.*","CS_stage",time_info$time3)
  time_info[time_info$time3%in%paste("week",9:13,sep=""),"time3"] <- "9-13 PCW"
  time_info[!time_info$time3%in%c("CS_stage","9-13 PCW"),"time3"] <- "After 16 PCW"
  time_info$time3 <- factor(time_info$time3,levels = c("CS_stage","9-13 PCW","After 16 PCW"))
  
  ori_1 <- do.call(cbind,tapply(df$`subtype_ori`, df$quintile, function(x){table(x) %>% prop.table()})) %>% as.data.frame()
  ori_4 <- t(ori_1) %>% as.data.frame()
  ori_4 <- ori_4 %>% mutate(quintile=rownames(ori_4)) %>% melt(id="quintile")
  ori_4$quintile <- as.numeric(ori_4$quintile) %>% as.factor()
  colnames(ori_4) <- c("module_score_quintile","ori","ori_distribution")
  ori_4$ori2 <- ori_4$ori
  ori_4$ori2 <- factor(ori_4$ori2,levels = c("pre_PraM","PraM"))
  
  ### MHC-II
  ls_plt[[paste("MHC-II",sub,sep="_")]] <- ggplot(df,aes(x=quintile,y=`MHC-II`))+
    geom_violin(fill=mac.color[sub],color=mac.color[sub],size=0.8,width=0.8)+
    geom_pointrange(data = mhc2_range,
                    mapping = aes(x = quintile, y = Median, ymin=`1st Qu.`,ymax=`3rd Qu.`),
                    color="purple",size=0.2,fill="purple") +
    geom_line(data = mhc2_range,
              mapping = aes(x = quintile, y = Median,group=1),size=2,color="purple")+
    theme(legend.position = "None")+
    xlab(paste(module,"score_quintile",sep="_"))+
    ylim(0,10)+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "None",
          axis.text = element_text(size=30),
          title = element_text(size=30),
          axis.line.y = element_line(color = ifelse(sub=="microglia","black",NA)))
  
  ## core-mf expression
  ls_plt[[paste("core-mf",sub,sep="_")]] <- ggplot(df,aes(x=quintile,y=`core-mf`))+
    geom_violin(fill=mac.color[sub],color=mac.color[sub],size=0.8,width=0.8)+
    geom_pointrange(data = core_mf_range,color="red",size=0.05,fill="red",
                    mapping = aes(x = quintile, y = Median, ymin=`1st Qu.`,ymax=`3rd Qu.`))+
    geom_line(data = core_mf_range,
              mapping = aes(x = quintile, y = Median,group=1),size=0.3,color="red")+
    xlab(paste(module,"score_quintile",sep="_"))+
    ggtitle(sub)+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "None",
          axis.text = element_text(size=30),
          title = element_text(size=30),
          axis.line.y = element_line(color = ifelse(sub=="microglia","black",NA)),
    )+
    ylim(0,6)
  
  ### module expression
  scale=max(phase$`G2M_pct`)/max(df[,module])
  
  ls_plt[[paste(module,sub,sep="_")]] <- ggplot(df)+
    geom_violin(aes(x=quintile,y=df[,module]),fill="#A6A8A9FF",color="#A6A8A9FF",size=0.8,width=0.8)+
    geom_pointrange(data = md_range,color="red",size=0.3,fill="red",
                    mapping = aes(x = quintile, y = Median, ymin=`1st Qu.`,ymax=`3rd Qu.`)) +
    geom_line(data = md_range,size=2,color="red",mapping = aes(x = quintile, y = Median,group=1))+
    geom_point(data=phase,aes(x=`module_score_quintile`,y=`G2M_pct`/scale),shape=21,size=4,color="black",fill="#219ebc")+
    geom_line(data=phase,aes(x=`module_score_quintile`,y=`G2M_pct`/scale),group=1,color="#219ebc",size=2)+
    scale_y_continuous(sec.axis = sec_axis(~.*scale, name="% proliferating cells")) +
    xlab(paste(module,"score_quintile",sep="_"))+ylab("module score")+
    ggtitle(sub)+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.text.y.left = element_text(size=30,colour = "black"),
          axis.text.y.right = element_text(size=30,colour = "black"),
          axis.line.y.left = element_line(colour = "red",size=2),
          axis.line.y.right = element_line(colour = "#219ebc",size=2),
          legend.position = "None",
          axis.text = element_text(size=30),
          title = element_text(size=30))  
  
  ### developmental time fraction
  ls_plt[[paste("time",sub,sep="_")]] <- ggplot(time_info,aes(x=`module_score_quintile`,y=`time_distribution`))+
    geom_bar(stat="identity",aes(fill=time3),width=0.6)+
    scale_fill_manual(values = c("#2a9d8f","#e9c46a","#e76f51"))+
    xlab(paste(module,"score_quintile",sep="_"))+
    scale_y_continuous(expand=c(0,0))+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "None",
          axis.text = element_text(size=30),
          title = element_text(size=30),
          axis.line.y = element_line(color = ifelse(sub=="microglia","black",NA))
    )
  
  ### fraction of pre_PraM vs PraM
  ls_plt[[paste("ori",sub,sep="_")]] <- ggplot(ori_4,aes(x=`module_score_quintile`,y=`ori_distribution`))+
    geom_bar(stat="identity",aes(fill=ori2),width=0.6)+
    scale_fill_manual(values = mac.color)+
    xlab(paste(module,"score_quintile",sep="_"))+
    ylab("Fraction of cells")+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=30),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=20,colour = "black"),
          legend.position = "None",
          axis.text = element_text(size=20),
          title = element_text(size=20),
          axis.line.y = element_line(color = "black"),
          axis.ticks.y = element_line(size = 2))+
    scale_y_continuous(expand=c(0,0))
  
  pdf_save(ls_plt,"~/metacell/Fig7/figs/7C/v1/")
  
  return(ls_plt)
}

##############################################################################
### Fig S7E PraM quintile by tissue 
plt_quintile_tissue <- function(umi_pct,module,sub){
  ls_plt <- list()
  mac.color <- c(mac.color,"pre_PraM & PraM"="#A8DADC")
  tissue <- unique(umi_pct$tissue)
  df <- umi_pct%>%filter(subtype==sub)
  df <- df[order(df[,module]),]
  rownames(df) <- 1:nrow(df)
  df$quintile <- factor(sort(rank(row.names(df))%%5)+1)
  mhc2_range <- tapply(df$`MHC-II`,df$quintile,summary) %>% do.call(rbind,.) %>% as.data.frame() %>% mutate(quintile=rownames(.))
  core_mf_range <- tapply(df$`core-mf`,df$quintile,summary) %>% do.call(rbind,.) %>% as.data.frame() %>% mutate(quintile=rownames(.))
  md_range <- tapply(df[,module],df$quintile,summary) %>% do.call(rbind,.) %>% as.data.frame() %>% mutate(quintile=rownames(.))
  
  
  df1 <- tapply(df$Phase, df$quintile, function(x){length(which(x=="G2M"))})
  df2 <- tapply(df$Phase, df$quintile, length)
  phase <- data.frame(G2M_pct=df1/df2*100,`module_score_quintile`=as.factor(1:5),total=df2)
  
  time1 <- do.call(cbind,tapply(df$time, df$quintile, function(x){table(x) %>% prop.table()})) %>% as.data.frame()
  time_info <- t(time1) %>% as.data.frame()
  time_info <- time_info %>% mutate(quintile=rownames(time_info)) %>% melt(id="quintile")
  time_info$quintile <- as.numeric(time_info$quintile) %>% as.factor()
  colnames(time_info) <- c("module_score_quintile","time","time_distribution")
  time_info$time3 <- time_info$time
  time_info$time3 <- gsub("cs.*","CS_stage",time_info$time3)
  time_info[time_info$time3%in%paste("week",9:13,sep=""),"time3"] <- "9-13 PCW"
  time_info[!time_info$time3%in%c("CS_stage","9-13 PCW"),"time3"] <- "After 16 PCW"
  time_info$time3 <- factor(time_info$time3,levels = c("CS_stage","9-13 PCW","After 16 PCW"))
  
  if(sub=="pre_PraM & PraM"){
    ori_1 <- do.call(cbind,tapply(df$`subtype_ori`, df$quintile, function(x){table(x) %>% prop.table()})) %>% as.data.frame()
    ori_4 <- t(ori_1) %>% as.data.frame()
    ori_4 <- ori_4 %>% mutate(quintile=rownames(ori_4)) %>% melt(id="quintile")
    ori_4$quintile <- as.numeric(ori_4$quintile) %>% as.factor()
    colnames(ori_4) <- c("module_score_quintile","ori","ori_distribution")
    ori_4$ori2 <- ori_4$ori
    ori_4$ori2 <- factor(ori_4$ori2,levels = c("pre_PraM","PraM"))
    
  }
  
  ### MHC-II
  ls_plt[[paste(tissue,"MHC-II",sub,sep="_")]] <- ggplot(df,aes(x=quintile,y=`MHC-II`))+
    geom_violin(fill=mac.color[sub],color=mac.color[sub],size=0.8,width=0.8)+
    geom_pointrange(data = mhc2_range,
                    mapping = aes(x = quintile, y = Median, ymin=`1st Qu.`,ymax=`3rd Qu.`),
                    color="purple",size=0.2,fill="purple") +
    geom_line(data = mhc2_range,
              mapping = aes(x = quintile, y = Median,group=1),size=2,color="purple")+
    theme(legend.position = "None")+
    xlab(paste(module,"score_quintile",sep="_"))+
    ylim(0,10)+
    # scale_y_continuous(limits = c(0,8),expand = c(0,6))+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          # axis.text.y = element_blank(),
          axis.text.y = element_text(size=30,colour = "black"),
          legend.position = "None",
          axis.text = element_text(size=30),
          title = element_text(size=30),
          axis.line.y = element_line(color = ifelse(sub=="microglia","black",NA)))

  ### module expression
  scale=max(phase$`G2M_pct`)/max(df[,module])
  
  ls_plt[[paste(tissue,module,sub,sep="_")]] <- ggplot(df)+
    geom_violin(aes(x=quintile,y=df[,module]),fill="#A6A8A9FF",color="#A6A8A9FF",size=0.8,width=0.8)+
    geom_pointrange(data = md_range,color="red",size=0.3,fill="red",
                    mapping = aes(x = quintile, y = Median, ymin=`1st Qu.`,ymax=`3rd Qu.`)) +
    geom_line(data = md_range,size=2,color="red",mapping = aes(x = quintile, y = Median,group=1))+
    geom_point(data=phase,aes(x=`module_score_quintile`,y=`G2M_pct`/scale),shape=21,size=4,color="black",fill="#219ebc")+
    geom_line(data=phase,aes(x=`module_score_quintile`,y=`G2M_pct`/scale),group=1,color="#219ebc",size=2)+
    scale_y_continuous(sec.axis = sec_axis(~.*scale, name="% proliferating cells")) +
    xlab(paste(module,"score_quintile",sep="_"))+ylab("module score")+
    ggtitle(tissue)+
    theme(axis.title.x = element_blank(),
          axis.text.y = element_text(size=30,colour = "black"),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y.left = element_text(size=30,colour = "black"),
          axis.text.y.right = element_text(size=30,colour = "black"),
          axis.line.y.left = element_line(colour = "red",size=2),
          axis.line.y.right = element_line(colour = "#219ebc",size=2),
          legend.position = "None",
          axis.text = element_text(size=30),
          title = element_text(size=30))  
  
  ### developmental time fraction
  ls_plt[[paste(tissue,"time",sub,sep="_")]] <- ggplot(time_info,aes(x=`module_score_quintile`,y=`time_distribution`))+
    geom_bar(stat="identity",aes(fill=time3),width=0.6)+
    scale_fill_manual(values = c("#2a9d8f","#e9c46a","#e76f51"))+
    xlab(paste(module,"score_quintile",sep="_"))+
    scale_y_continuous(expand=c(0,0))+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size=30,colour = "black"),
          legend.position = "None",
          axis.text = element_text(size=30),
          title = element_text(size=30),
          axis.line.y = element_line(color = ifelse(sub=="microglia","black",NA))
    )
  ### proportion of pre-PraM vs PraM
  if(sub=="pre_PraM & PraM"){
    ls_plt[[paste(tissue,"ori",sub,sep="_")]] <- ggplot(ori_4,aes(x=`module_score_quintile`,y=`ori_distribution`))+
      geom_bar(stat="identity",aes(fill=ori2),width=0.6)+
      scale_fill_manual(values = mac.color)+
      xlab(paste(module,"score_quintile",sep="_"))+
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size=30),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=20,colour = "black"),
            legend.position = "None",
            axis.text = element_text(size=20),
            title = element_text(size=20),
            axis.line.y = element_line(color = "black"),
            axis.ticks.y = element_line(size = 2))+
      ylab("Fraction of cells")+
      scale_y_continuous(expand=c(0,0)) 
  }
  pdf_save(ls_plt,paste("~/metacell/Fig7/figs/S7/tissue_quintile/with_Y/",sub,"_tissue/",sep=""))
  return(ls_plt)
}
