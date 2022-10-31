

the_metatable = metatable
the_metatable = the_metatable[the_metatable$major != "erythrocyte",]
the_metatable$new = paste0(the_metatable$tissue,the_metatable$major)
for (new in unique(the_metatable$new)) {
  if (nrow(the_metatable[the_metatable$new == new,])<10){
    the_metatable= the_metatable[-which(the_metatable$new == new),]
  }
}

the_metatable = the_metatable[,c("time","major","Well_ID","tissue")]
the_metatable = the_metatable[-which(the_metatable$time == "Adult"),]
the_metatable = the_metatable[-which(the_metatable$time == "week19"),]
the_metatable$time <- factor(the_metatable$time,levels=time_order,ordered = TRUE)

tmp2 = the_metatable[,c("time","major")]

tmp2 = tmp2[!duplicated(tmp2),]
rownames(tmp2) = NULL
tmp2$rate = 0
for(t in rownames(tmp2)){
  tmp2[t,"rate"] = nrow(the_metatable[the_metatable$time == tmp2[t,"time"] & the_metatable$major == tmp2[t,"major"],])/nrow(the_metatable[the_metatable$time == tmp2[t,"time"],])
}

p1 = ggplot(tmp2,aes(x = time, 
                     stratum = major, 
                     alluvium = major,
                     y = rate,
                     fill = major,
                     label = major)) +
  scale_y_continuous(label = scales::percent_format(),
                     expand=c(0,0)) +
  scale_fill_manual(values = as.character(major_color),breaks = names(major_color)) +
  geom_flow() +
  geom_stratum(width = 0.6,colour = "grey40") +
  theme(axis.title=element_blank(),
        axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.title=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(),
        axis.ticks=element_line(),
        legend.position = "none",
        axis.text = element_text(size = 12))
ggsave("E:\\fig1/ratio_legend.pdf",p1,width =4,height =3,dpi = 300)

###
figre2_major_ratio_per_tissue = function(the_metatable,the_dir,thew , theh){
  tmp = the_metatable
  tmp2 = tmp
  
  for(i in unique(tmp2$tissue)){
    tmp = tmp2[,c("tissue","time","subtype")]
    tmp = tmp[tmp$tissue == i,]
    tmp = tmp[,c("time","subtype")]
    other = setdiff(time_order,c(unique(tmp$time),"Adult","week19"))
    other = as.data.frame(other)
    other$blank = "blank"
    colnames(other) = colnames(tmp)
    tmp = rbind(tmp,other)
    tmp$time <- factor(tmp$time,levels=time_order,ordered = TRUE)
    p2 = ggplot(tmp,aes(x= time,fill = subtype))+
      geom_bar(stat="count",position = "fill",width = 0.8,color = "black")+
      xlab("")+ylab("")+ggtitle(paste(i," total cell number ",nrow(tmp),sep = ""))+
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
    ggsave(paste(the_dir,i,"_tissue.pdf",sep = ""),p2,width = thew,height = theh)
  }
}

### Fig 2 significance test
### Fig 2 B & D
library(ggsignif)
library(ggplot2)
source("Funcs_signif_fig2.R")
setwd("./path")
load("mac_color.Rdata")
load("major_color.Rdata")
load("metatable_0926.Rdata")
#######################################################################################
# Fig 2B major cell type proportion significance test
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

plt_Frac_signif(ratio, "./signif/")


#######################################################################################

### Fig 2D major cell type proportion significance test in each tissue
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

# Figure 2E

load("G://embryo/20220929/metatable_0929.Rdata")


data1<-metatable

##############################
# remove Adult

data1<-data1[!data1$time%in%c("Adult"),]

majorType<-unique(data1$major)

majoryType.cellsum<-vector()

for(i in 1:length(majorType)){
  
  majoryType.cellsum[i]<-nrow(subset(data1,major==majorType[i]))
  
  
}


names(majoryType.cellsum)<-majorType


#####################################3

minorType<-unique(data1$subtype)

minorType.cellsum<-vector()

for(i in 1:length(minorType)){
  
  minorType.cellsum[i]<-nrow(subset(data1,subtype==minorType[i]))
  
  
}

names(minorType.cellsum)<-minorType

tissue.num<-length(unique(data1$tissue))
tissue<-unique(data1$tissue)

tissue_order<-read.csv("G://embryo/20220801/20220720/tissue_order.csv",header=T)

subtype_order<-read.csv("G://embryo/20220802/subtypes_new2.csv",header=T)

heatmap.data<-matrix(data=NA, ncol = tissue.num , nrow =length(minorType), byrow = FALSE, dimnames = NULL)


colnames(heatmap.data)<-as.matrix(tissue_order)
rownames(heatmap.data)<-as.matrix(subtype_order$subtype)
minorType_order<-subtype_order$subtype
tissue_order<-as.character(as.matrix(tissue_order))


for(i in 1:tissue.num){
  
  for(j in 1:length(minorType) ){
    
    
    numerator = nrow(subset(data1,subtype == minorType_order[j] & tissue == tissue_order[i]))
    
    denominator = nrow(subset(data1,subtype == minorType_order[j]))
    
    ratio = numerator/denominator
    
    heatmap.data[j,i]<-ratio
    
  }
}

####

heatmap.data<-t(heatmap.data)

annotation_col<-subtype_order
rownames(annotation_col)<-annotation_col$subtype
annotation_col<-subset(annotation_col,select = -subtype)

ann_colors = list(
  major = major_color
)

bk <- c(seq(0,0.2,by=0.01),seq(0.21,0.6,by=0.01))
p<-pheatmap(heatmap.data,breaks = bk,gaps_col = c(6,9,12,18,23,29,32,35,41),
            color = c(colorRampPalette(colors = c("white","firebrick"))(length(bk)/2),colorRampPalette(colors = c("firebrick","Maroon"))(length(bk)/2)),show_colnames=T,display_numbers = F,cluster_cols =F,cluster_rows = F,fontface="bold",annotation_colors = ann_colors,annotation_col = annotation_col)

ggsave(p,filename = "G://embryo/figures/fig2/fig2E.png",width = 20, height = 8, units = "in", device='png')
ggsave(p,filename = "G://embryo/figures/fig2/fig2E.pdf",width = 20, height = 8, units = "in", device='pdf')
