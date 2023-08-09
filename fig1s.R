# Figure 1S

#####
###Fig1S A
#####

human_development?

human_development = human_development[, metatable$Cell[metatable$time != "Adult"]]

process_cells_annotation(human_development)

plot_cells_annotation(human_development,type="histogram")
plot_cells_annotation(human_development,type="boxplot")
plot_UMIs_vs_Detected_genes(human_development)


#####
###Fig1S B
#####










#####
###Fig1S C
#####
library(ggsignif)
library(ggplot2)
source("~/Funcs_signif_fig2.R")
setwd("./path")
#load("mac_color.Rdata")
load("~/major_color.Rdata") # loading color scheme for major lineages
load("~/data/metatable_0929.Rdata") # loading single cell metatable 



#######################################################################################

metatable <- metatable %>% filter(time!="Adult")
metatable <- metatable %>% filter(major!="erythrocyte") # remove erythrocyte
metatable <- metatable %>% filter(embryo!="embryo 10") # embryo 10 contain 1 cell only
table(metatable$major) %>% sort
metatable %>% filter(subtype=="discard") %>% pull(major) %>% table 
metatable <- metatable %>% filter(subtype!="discard") # drop discard annotation
metatable$major <- as.factor(metatable$major)

### divide all time into 3 periods
CS_stage = c("cs11","cs12","cs13","cs14","cs18","cs19","cs21","cs23")
`9-13 PCW` = c("week9","week10","week11","week12","week13")
`After 16 PCW` = c("week16","week19","week20","week23","week27")
metatable$period <- metatable$time
metatable$period[metatable$period%in%CS_stage] <- "CS_stage"
metatable$period[metatable$period%in%`9-13 PCW`] <- "9-13 PCW"
metatable$period[metatable$period%in%`After 16 PCW`] <- "After 16 PCW"
### filter out adult cells
meta_fetal <- metatable %>% filter(time!="Adult")
table(meta_fetal$period)
### compute ratio of subtypes in different time periods
ratio <- list()
for(i in unique(meta_fetal$embryo)){
  df <- meta_fetal %>% filter(embryo==i)
  ratio[[i]] <- tapply(df$major,df$period,function(x){prop.table(table(x))}) %>% do.call(cbind,.) %>% as.data.frame()
  ratio[[i]] <- ratio[[i]] %>% mutate(major=rownames(.),period=colnames(.),embryo=i)
  colnames(ratio[[i]])[1] <- "ratio"
}
ratio <- do.call(rbind,ratio)
ratio$period <- factor(ratio$period,levels = c("CS_stage","9-13 PCW","After 16 PCW"))
### draw median horizontal line for each group; control segment length
median <- tapply(ratio$ratio,list(ratio$major,ratio$period),median) %>% melt()
colnames(median) <- c("major","period","median")
median$x_start <- as.numeric(median$period)-0.3
median$x_end <- as.numeric(median$period)+0.3
my_comp <- list(c("CS_stage","9-13 PCW"),c("9-13 PCW","After 16 PCW"),c("CS_stage","After 16 PCW"))
color_grp <- c("#A6A8A9FF","#484A4BFF","#141414FF")
names(color_grp) <- sort(unique(ratio$period))




#######################################################################################

### Fig1S C major cell type proportion significance test in each tissue
### drop small sample of embryo-tissue  pair (e.g., embryo14-adrenalgland 1 cell)
tissue_less <- tapply(meta_fetal$tissue,meta_fetal$embryo, function(x){data.frame(table(x))})
for(i in 1:length(tissue_less)){
  tissue_less[[i]] <- tissue_less[[i]] %>% mutate(embryo=names(tissue_less)[i])
  colnames(tissue_less[[i]])[1] <- "tissue"
  tissue_less[[i]] <- tissue_less[[i]][,c("tissue","embryo","Freq")]
}
tissue_less <- do.call(rbind,tissue_less)
tissue_less <- tissue_less %>% filter(Freq<50)
### some tissue only exist in cs stage, separate and divide them into 3 groups
Ftissue <- meta_fetal[,c("embryo","tissue","time","period")]
Ftissue <- Ftissue[!duplicated(Ftissue),]

tissue_early <- c("AGM","Embryo","Head","Limb","Primitive_gut","female_gonad")
Ftissue_early <- Ftissue %>% filter(tissue%in%tissue_early)
tapply(Ftissue_early$time,Ftissue_early$tissue,table)
meta_fetal_more <- meta_fetal %>% filter(!tissue%in%tissue_early)

### compute ratio for each embryo, each tissue
ratio_tissue_more <- Frac_tissue(meta_fetal_more)

### check embryo number for each period
r_info <- lapply(ratio_tissue_more,function(x){tapply(x$embryo,x$period,function(x){length(unique(x))})})
tissue_auto <- lapply(r_info,function(x)(which(length(which(x>=3))==3))) %>% unlist() %>% names()
tissue_manual <- names(r_info) %>% setdiff(tissue_auto) %>% setdiff("Female gonad")

lapply(ratio_tissue_more[tissue_manual],function(x){tapply(x$embryo,x$time,function(x){length(unique(x))})}) 

#######################################################################################
# tissue with adequate samples for dividing into CS; 9-13 PCW; After 16 PCW
ratio_tissue_auto <- ratio_tissue_more[tissue_auto]

test_tissue_auto <- list()
for(j in names(ratio_tissue_auto)){
  test_all <- data.frame()
  for(i in unique(ratio_tissue_auto[[j]]$major)){
    ratio_sub <- ratio_tissue_auto[[j]] %>% filter(major==i)
    ratio_12 <- ratio_sub %>% filter(period%in%c("CS_stage", "9-13 PCW"))
    ratio_23 <- ratio_sub %>% filter(period%in%c("9-13 PCW", "After 16 PCW"))
    ratio_13 <- ratio_sub %>% filter(period%in%c("CS_stage", "After 16 PCW"))
    
    p.val = c(wilcox.test(ratio~period,data=ratio_12)$p.value,wilcox.test(ratio~period,data=ratio_23)$p.value,
              wilcox.test(ratio~period,data=ratio_13)$p.value)
    test <- data.frame(major=i,p.val_1=p.val[1],p.val_2=p.val[2],p.val_3=p.val[3])
    test_all <- rbind(test_all,test)
  }
  test_tissue_auto[[j]] <- test_all %>% mutate(tissue=j)
}
test_tissue_auto <- do.call(rbind,test_tissue_auto) 
test_tissue_auto <- test_tissue_auto %>% filter(p.val_1<0.05|p.val_2<0.05|p.val_3<0.05)

plt_Frac_tissue_signif(ratio_tissue_auto, test_tissue_auto,"./signif/")

#######################################################################################
# tissue with inadequate samples for dividing into CS; 9-13 PCW; After 16 PCW
# divide time period group for each tissue respectively

ratio_tissue_manual <- ratio_tissue_more[tissue_manual]
for(i in names(ratio_tissue_manual)){
  ratio_tissue_manual[[i]]$period <- as.character(ratio_tissue_manual[[i]]$period)
}

### Spinalmarrow
ratio_tissue_manual$Spinalmarrow[ratio_tissue_manual$Spinalmarrow$time%in%c("cs19","cs21","cs23"),"period"] <- "CS_stage"
ratio_tissue_manual$Spinalmarrow[ratio_tissue_manual$Spinalmarrow$time%in%c("week9","week10","week11"),"period"] <- "9-11 PCW"
ratio_tissue_manual$Spinalmarrow[ratio_tissue_manual$Spinalmarrow$time%in%c("week12","week13","week16"),"period"] <- "12-16 PCW"
ratio_tissue_manual$Spinalmarrow$period <- factor(ratio_tissue_manual$Spinalmarrow$period,levels = c("CS_stage","9-11 PCW","12-16 PCW"))

### Bonemarrow
ratio_tissue_manual$Bonemarrow[ratio_tissue_manual$Bonemarrow$time%in%c("week9","cs23"),"period"] <- "CS23 - 9 PCW"
ratio_tissue_manual$Bonemarrow[ratio_tissue_manual$Bonemarrow$time%in%c("week10","week11","week12"),"period"] <- "10-12 PCW"
ratio_tissue_manual$Bonemarrow[ratio_tissue_manual$Bonemarrow$time%in%c("week13","week16","week23"),"period"] <- "13-23 PCW"
ratio_tissue_manual$Bonemarrow$period <- factor(ratio_tissue_manual$Bonemarrow$period,levels = c("CS23 - 9 PCW","10-12 PCW","13-23 PCW"))

### `Male gonad`
ratio_tissue_manual$`Male gonad`[ratio_tissue_manual$`Male gonad`$time%in%c("cs21"),"period"] <- "CS21"
ratio_tissue_manual$`Male gonad`[ratio_tissue_manual$`Male gonad`$time%in%c("week9","week10","week11"),"period"] <- "9-11 PCW"
ratio_tissue_manual$`Male gonad`[ratio_tissue_manual$`Male gonad`$time%in%c("week12","week13","week16"),"period"] <- "12-16 PCW"
ratio_tissue_manual$`Male gonad`$period <- factor(ratio_tissue_manual$`Male gonad`$period,levels = c("CS21","9-11 PCW","12-16 PCW"))

### Thymus
ratio_tissue_manual$Thymus[ratio_tissue_manual$Thymus$time%in%c("cs21","cs23","week9"),"period"] <- "CS21-week9"
ratio_tissue_manual$Thymus[ratio_tissue_manual$Thymus$time%in%c("week10","week11"),"period"] <- "10-11 PCW"
ratio_tissue_manual$Thymus[ratio_tissue_manual$Thymus$time%in%c("week12","week13","week16"),"period"] <- "12-16 PCW"
ratio_tissue_manual$Thymus$period <- factor(ratio_tissue_manual$Thymus$period,levels = c("CS21-week9","10-11 PCW","12-16 PCW"))

### Spleen
ratio_tissue_manual$Spleen[ratio_tissue_manual$Spleen$time%in%c("week9","week10","week11"),"period"] <- "9-11 PCW"
ratio_tissue_manual$Spleen[ratio_tissue_manual$Spleen$time%in%c("week12","week13"),"period"] <- "12-13 PCW"
ratio_tissue_manual$Spleen[ratio_tissue_manual$Spleen$time%in%c("week16","week20","week23"),"period"] <- "16-23 PCW"
ratio_tissue_manual$Spleen$period <- factor(ratio_tissue_manual$Spleen$period,levels = c("9-11 PCW","12-13 PCW","16-23 PCW"))

### Yolksac
ratio_tissue_manual$Yolksac[ratio_tissue_manual$Yolksac$time%in%c("cs11","cs12"),"period"] <- "CS11-12"
ratio_tissue_manual$Yolksac[ratio_tissue_manual$Yolksac$time%in%c("cs13"),"period"] <- "CS13"
ratio_tissue_manual$Yolksac[ratio_tissue_manual$Yolksac$time%in%c("cs14","cs21"),"period"] <- "CS14-21"
ratio_tissue_manual$Yolksac$period <- factor(ratio_tissue_manual$Yolksac$period,levels = c("CS11-12","CS13","CS14-21"))

### filter tissue
test_tissue_manual <- list()
for(j in names(ratio_tissue_manual)){
  test_all <- data.frame()
  for(i in unique(ratio_tissue_manual[[j]]$major)){
    ratio_sub <- ratio_tissue_manual[[j]] %>% filter(major==i)
    f1 <- sort(unique(ratio_sub$period))[1]
    f2 <- sort(unique(ratio_sub$period))[2]
    f3 <- sort(unique(ratio_sub$period))[3]
    ratio_12 <- ratio_sub %>% filter(period%in%c(f1, f2))
    ratio_23 <- ratio_sub %>% filter(period%in%c(f2, f3))
    ratio_13 <- ratio_sub %>% filter(period%in%c(f1, f3))
    
    p.val = c(wilcox.test(ratio~period,data=ratio_12)$p.value,wilcox.test(ratio~period,data=ratio_23)$p.value,
              wilcox.test(ratio~period,data=ratio_13)$p.value)
    test <- data.frame(major=i,p.val_1=p.val[1],p.val_2=p.val[2],p.val_3=p.val[3])
    test_all <- rbind(test_all,test)
  }
  test_tissue_manual[[j]] <- test_all %>% mutate(tissue=j)
}
test_tissue_manual <- do.call(rbind,test_tissue_manual) 
test_tissue_manual <- test_tissue_manual %>% filter(p.val_1<0.05|p.val_2<0.05|p.val_3<0.05)

plt_Frac_tissue_signif(ratio_tissue_manual, test_tissue_manual,"./signif/")

#######################################################################################
# tissue with few samples for dividing into CS; 9-13 PCW; After 16 PCW
# divide time period group for each tissue respectively

meta_fetal_less <- meta_fetal %>% filter(tissue%in%tissue_early) %>% filter(tissue!="Embryo")
ratio_tissue_less <- Frac_tissue(meta_fetal_less)

for(i in names(ratio_tissue_less)){
  ratio_tissue_less[[i]]$period <- as.character(ratio_tissue_less[[i]]$period)
}

### AGM
ratio_tissue_less$AGM[ratio_tissue_less$AGM$time%in%c("cs13"),"period"] <- "CS13"
ratio_tissue_less$AGM[ratio_tissue_less$AGM$time%in%c("cs14","cs18"),"period"] <- "CS14-18"
ratio_tissue_less$AGM$period <- factor(ratio_tissue_less$AGM$period,levels = c("CS13","CS14-18"))

### Head
ratio_tissue_less$Head[ratio_tissue_less$Head$time%in%c("cs13"),"period"] <- "CS13"
ratio_tissue_less$Head[ratio_tissue_less$Head$time%in%c("cs14"),"period"] <- "CS14"
ratio_tissue_less$Head$period <- factor(ratio_tissue_less$Head$period,levels = c("CS13","CS14"))

### Limb
ratio_tissue_less$Limb[ratio_tissue_less$Limb$time%in%c("cs18","cs19"),"period"] <- "CS18-19"
ratio_tissue_less$Limb[ratio_tissue_less$Limb$time%in%c("cs21","cs23"),"period"] <- "CS21-23"
ratio_tissue_less$Limb$period <- factor(ratio_tissue_less$Limb$period,levels = c("CS18-19","CS21-23"))

### filter tissue
test_tissue_less <- list()
for(j in names(ratio_tissue_less)){
  test_all <- data.frame()
  for(i in unique(ratio_tissue_less[[j]]$major)){
    ratio_sub <- ratio_tissue_less[[j]] %>% filter(major==i)
    p.val = c(wilcox.test(ratio~period,data=ratio_sub)$p.value)
    test <- data.frame(major=i,p.val=p.val)
    test_all <- rbind(test_all,test)
  }
  test_tissue_less[[j]] <- test_all %>% mutate(tissue=j)
}
test_tissue_less <- do.call(rbind,test_tissue_less) 
test_tissue_less <- test_tissue_less %>% filter(p.val<0.05)


# grey palette ------------------------------------------------------------

plt_Frac_tissue_signif_grey_l(ratio_tissue_auto, test_tissue_auto,"./grey_l/")
plt_Frac_tissue_signif_grey_l(ratio_tissue_manual, test_tissue_manual,"./grey_l/")


#####
###Fig1S D
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
### Fig1s E
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
### Fig1s F
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
