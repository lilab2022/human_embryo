##Figure 6 
##writer: Dr. Wang Hao
##Date: 2022-10-4

##Figure 6C

load("G://embryo/20220929/metatable_0929.Rdata")
load("G://embryo/20220929/mac.color2.rdata")


load("G://embryo/20220804/FDG_withoutAdult.rdata")

selected.subtype<-c("YsdM_AFP_low","YsdM_AFP_high","pre_PraM","PraM","HepM","microglia","monocyte","DC-MG Prog","ACY3+S100B+pMφ")

time.set<-list(
  
  selected.time0<-c("cs11","cs12"),
  selected.time1<-c("cs11","cs12","cs13","cs14"),
  selected.time2<-c("cs11","cs12","cs13","cs14","cs18","cs19"),
  selected.time3<-c("cs11","cs12","cs13","cs14","cs18","cs19","cs21","cs23"),
  selected.time4<-c("cs11","cs12","cs13","cs14","cs18","cs19","cs21","cs23","week9"),
  selected.time5<-c("cs11","cs12","cs13","cs14","cs18","cs19","cs21","cs23","week9","week11","week12","week13"),
  selected.time6<-c("cs11","cs12","cs13","cs14","cs18","cs19","cs21","cs23","week9","week11","week12","week13","week16","week19","week20","week23","week27")
  
)

fa2.data<-data1[data1$subtype%in%selected.subtype,]

for(i in 1:length(time.set)) {
  
  
  fa2.datax<-fa2.data[fa2.data$time%in%time.set[[i]],]
  
  
  xx<-data1[!data1$Well_ID%in%fa2.datax$Well_ID,]
  xx2<-rbind(fa2.datax,xx)

  xx2$FDG11=NA
  xx2$FDG22=NA
  xx2$FDG11[xx2$Well_ID %in% fa2.datax$Well_ID]=xx2$FDG1[xx2$Well_ID %in% fa2.microgliax$Well_ID]
  xx2$FDG22[xx2$Well_ID %in% fa2.datax$Well_ID]=xx2$FDG2[xx2$Well_ID %in% fa2.microgliax$Well_ID]
  
  p1<-ggplot(xx2)+
    geom_point(aes(x = FDG1, y = FDG2),shape = 19,stroke = 0.2,size = 3,color = "grey88")+
    geom_point(data=xx2[xx2$major=="macrophage",],aes(x=FDG11,y=FDG22,fill=subtype),shape = 21,stroke = 0.2,size = 3)+
    scale_fill_manual(values = as.character(mac.color),breaks = names(mac.color))+
    geom_point(data=xx2[xx2$subtype%in%c("YsdM_AFP_low","YsdM_AFP_high"),],aes(x=FDG11,y=FDG22,fill=subtype),shape = 21,stroke = 0.2,size = 3)+
    theme_bw()+theme(legend.position = 'none')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank()
    )
  
  ggsave(p1,filename = paste("G://embryo/figures/fig6/fig6C/accum_time_bifurcation_",i,".png",sep=""),width = 32, height = 25, dpi = 150, units = "in", device='png')
  ggsave(p1,filename = paste("G://embryo/figures/fig6/fig6C/accum_time_bifurcation_",i,".pdf",sep=""),width = 15, height = 12, units = "in", device='pdf')
  
}


##Figure6 D
load("G://embryo/20220801/20220720/major_color.Rdata")
load("G://embryo/20220929/mac.color2.rdata")


data<-read.csv("G://embryo/20220329/cleanumi.csv",header=T,row.names = "X")

load("G://embryo/20220929/metatable_0929.Rdata")

mcmc1<-metatable

adult.mc<-c("X2056","X2350","X282","X871","X891") # MCs enriched with adult cells

mcmc1 <- mcmc1[!mcmc1$mc%in%adult.mc,]
all_mac = c(unique(mcmc1$mc[mcmc1$subtype %in% c("pre_PraM","PraM","HepM","microglia","YsdM_AFP_high","YsdM_AFP_low")]))

anno.2 = mcmc1[,c("mc","subtype","major")]
anno.2 = anno.2[!duplicated(anno.2),]
anno.2 = anno.2[anno.2$mc %in% all_mac,]
tmp = anno.2$mc 
rownames(anno.2) = tmp
anno.2 = anno.2[,-1]

data_sub<-data[,colnames(data)%in%tmp]
excluded_genes<-discard_gene2(rownames(data_sub))


seob_data <- CreateSeuratObject(counts = data_sub)
seob_data <- NormalizeData(seob_data, normalization.method =  "LogNormalize")



gene_counts<-seob_data@assays[["RNA"]]@counts  #original data
gene_exp <- seob_data@assays[["RNA"]]@data   #normalized data

gene_counts<-gene_counts[rownames(gene_counts)%in%excluded_genes,]
gene_exp<-gene_exp[rownames(gene_exp)%in%excluded_genes,]


gene_counts<-as.matrix(gene_counts)
gene_exp<-as.matrix(gene_exp)

###

cell_info2 <- rownames_to_column(anno.2, 
                                 var = "cell_id")

cell_info <- rownames_to_column(seob_data@meta.data, 
                                var = "cell_id")

cell_info<-merge(cell_info,cell_info2,by="cell_id")



dynob <- wrap_expression(
  expression = t(gene_exp), 
  counts = t(gene_counts), 
  cell_info = cell_info
)


start_cells<-subset(cell_info,subtype=="YsdM_AFP_high")
#start_cells<-subset(cell_info,subtype=="HSC/MPP")


dynob <- add_prior_information(
  dynob,
  start_id = start_cells$cell_id
)

dynob <- add_grouping(
  dynob,
  cell_info$subtype
)
head(dynob$grouping)


guidelines_shiny(dataset = dynob)



# Reproduces the guidelines as created in the shiny app
answers <- dynguidelines::answer_questions(
  multiple_disconnected = FALSE, 
  expect_topology = TRUE, 
  expected_topology = "bifurcation", 
  n_cells = 1003, 
  n_features = 23214, 
  memory = "7GB", 
  prior_information = "start_id", 
  docker = TRUE
)
guidelines <- dynguidelines::guidelines(answers = answers)

methods_selected <- guidelines$methods_selected
methods_selected

#docker is required here for dyno running

model3 <- infer_trajectory(dynob, 
                           method = "slingshot", 
                           give_priors = c("start_id"))


p21<-plot_dimred(model3, 
                 label_milestones = F,
                 plot_milestone_network = F,
                 plot_trajectory = T,
                 size_milestones = 2.5,
                 size_cells = 5,
                 color_density = "grouping",
                 density_cutoff_label = F,
                 border_radius_percentage = 0.4,
                 color_cells = "grouping",
                 groups = tribble(~group_id,~color,"pre_PraM","#2CA02CFF","PraM","#FF7F0EFF","microglia","#1F77B4FF","YsdM_AFP_high","#9467BDFF",
                                  "YsdM_AFP_low","#E377C2FF","HepM","#FFBB78FF"),
                 grouping = dynob$grouping)+theme(legend.position="none")


ggsave("G://embryo/figures/fig6/fig6_bifurcation.png",p21,width = 6,height = 7,dpi = 300,units="in",device="png")
ggsave("G://embryo/figures/fig6/fig6_bifurcation.pdf",p21,width = 6,height = 7,dpi = 300,units="in",device="pdf")

##########

color_matrix3<-read.table("G://embryo/20220930/fp_module(2).csv")
color_matrix3<-t(color_matrix3)


Time<-read.csv("G://embryo/20220903/allcells_realtime_noadult.csv",row.names = "X")
colnames(Time)[1]<-"MC"

rownames(Time)<-Time$MC
source.x<-merge(Time,color_matrix3,by="row.names")
rownames(source.x)<-source.x$Row.names
source.x<-source.x[,c(-1,-2)]
source.x<-source.x[rownames(source.x)%in%cell_info$cell_id,]

source.x<-as.matrix(source.x)

p3<-plot_dimred(model3, 
                label_milestones = F,
                plot_milestone_network = F,
                plot_trajectory = T,
                size_milestones = 2.5,
                size_cells = 5,
                border_radius_percentage = 0.4,
                feature_oi = "the_median",
                color_density = "grouping",
                grouping = dynob$grouping,
                groups = tribble(~group_id,~color,"pre_PraM","#2CA02CFF","PraM","#FF7F0EFF","microglia","#1F77B4FF","YsdM_AFP_high","#9467BDFF","YsdM_AFP_low","#E377C2FF","HepM","#FFBB78FF"),
                color_cells = "feature",
                density_cutoff_label = F,
                expression_source = source.x
)+theme(legend.position="none")


ggsave("G://embryo/figures/fig6/fig6_bifurcation_Time.png",p3,width = 6,height = 7,dpi = 300,units="in",device="png")
ggsave("G://embryo/figures/fig6/fig6_bifurcation_Time.pdf",p3,width = 6,height = 7,dpi = 300,units="in",device="pdf")



p3<-plot_dimred(model3, 
                label_milestones = F,
                plot_milestone_network = F,
                plot_trajectory = T,
                size_milestones = 2.5,
                size_cells = 5,
                border_radius_percentage = 0.4,
                feature_oi = "PraM",
                color_density = "grouping",
                grouping = dynob$grouping,
                groups = tribble(~group_id,~color,"pre_PraM","#2CA02CFF","PraM","#FF7F0EFF","microglia","#1F77B4FF","YsdM_AFP_high","#9467BDFF","YsdM_AFP_low","#E377C2FF","pre_microglia","#FFBB78FF"),
                color_cells = "feature",
                density_cutoff_label = F,
                expression_source = source.x
)+theme(legend.position="none")



ggsave("G://embryo/figures/fig6/fig6_bifurcation_PraM.png",p3,width = 6,height = 7,dpi = 300,units="in",device="png")
ggsave("G://embryo/figures/fig6/fig6_bifurcation_PraM.pdf",p3,width = 6,height = 7,dpi = 300,units="in",device="pdf")


p3<-plot_dimred(model3, 
                label_milestones = F,
                plot_milestone_network = F,
                plot_trajectory = T,
                size_milestones = 2.5,
                size_cells = 5,
                border_radius_percentage = 0.4,
                feature_oi = "microglia",
                color_density = "grouping",
                grouping = dynob$grouping,
                groups = tribble(~group_id,~color,"pre_PraM","#2CA02CFF","PraM","#FF7F0EFF","microglia","#1F77B4FF","YsdM_AFP_high","#9467BDFF","YsdM_AFP_low","#E377C2FF","pre_microglia","#FFBB78FF"),
                color_cells = "feature",
                density_cutoff_label = F,
                expression_source = source.x
)+theme(legend.position="none")



ggsave("G://embryo/figures/fig6/fig6_bifurcation_microglia.png",p3,width = 6,height = 7,dpi = 300,units="in",device="png")
ggsave("G://embryo/figures/fig6/fig6_bifurcation_microglia.pdf",p3,width = 6,height = 7,dpi = 300,units="in",device="pdf")



differentiation<-read.csv("G://embryo/20220808/cytotrace.value.bymc.csv",row.names = "MC")


source.x2<-merge(source.x,differentiation,by=0)

rownames(source.x2)<-source.x2$Row.names
source.x2<-source.x2[,-1]
source.x2<-as.matrix(source.x2)

p3<-plot_dimred(model3, 
                label_milestones = F,
                plot_milestone_network = F,
                plot_trajectory = T,
                size_milestones = 2.5,
                size_cells = 5,
                border_radius_percentage = 0.4,
                feature_oi = "cytotrace",
                color_density = "grouping",
                grouping = dynob$grouping,
                groups = tribble(~group_id,~color,"pre_PraM","#2CA02CFF","PraM","#FF7F0EFF","microglia","#1F77B4FF","YsdM_AFP_high","#9467BDFF","YsdM_AFP_low","#E377C2FF","pre_microglia","#FFBB78FF"),
                color_cells = "feature",
                density_cutoff_label = F,
                expression_source = source.x2
)+theme(legend.position="none")


ggsave("G://embryo/figures/fig6/fig6_bifurcation_differentiation.png",p3,width = 6,height = 7,dpi = 300,units="in",device="png")
ggsave("G://embryo/figures/fig6/fig6_bifurcation_differentiation.pdf",p3,width = 6,height = 7,dpi = 300,units="in",device="pdf")


##Figure 6 E
corr<-cor(meta.info2[,4],meta.info2[,6])
corr<-round(corr,digits=4)



p<-ggplot(meta.info2)+
  geom_point(aes(x=pseudoTime,y=Monocle3,fill=Subtype),shape = 21,stroke = 0.2,size = 3)+
  scale_fill_manual(values = as.character(mac.color),breaks = names(mac.color))+
  theme_bw() + theme(legend.position = 'none')+
  theme(
    
    panel.grid=element_blank()
  )+
  labs(x = "Slingshot", y = "Monocle3", color = '',subtitle = paste("Correlation:",corr,sep="")) + geom_smooth(aes(x=pseudoTime,y=Monocle3),method = "lm")




ggsave(p,filename = paste("G://embryo/20220930/Slingshot_VS_Monocle3__plot_with_dots.png",sep=""),width = 5, height = 5, units = "in", device='png')
ggsave(p,filename = paste("G://embryo/20220930/Slingshot_VS_Monocle3__plot_with_dots.pdf",sep=""),width = 5, height = 5, units = "in", device='pdf')



corr<-cor(meta.info2[,4],meta.info2[,2])
corr<-round(corr,digits=4)



p<-ggplot(meta.info2)+
  geom_point(aes(x=pseudoTime,y=Time,fill=Subtype),shape = 21,stroke = 0.2,size = 3)+
  scale_fill_manual(values = as.character(mac.color),breaks = names(mac.color))+
  theme_bw() + theme(legend.position = 'none')+
  theme(
    
    panel.grid=element_blank()
  )+
  labs(x = "Slingshot", y = "Actual time", color = '',subtitle = paste("Correlation:",corr,sep="")) + geom_smooth(aes(x=pseudoTime,y=Time),method = "lm")




ggsave(p,filename = paste("G://embryo/20220930/Slingshot_VS_ActualTime__plot_with_dots.png",sep=""),width = 5, height = 5, units = "in", device='png')
ggsave(p,filename = paste("G://embryo/20220930/Slingshot_VS_ActualTime__plot_with_dots.pdf",sep=""),width = 5, height = 5, units = "in", device='pdf')




corr<-cor(meta.info2[,4],meta.info2[,5])
corr<-round(corr,digits=4)



p<-ggplot(meta.info2)+
  geom_point(aes(x=pseudoTime,y=gene.module,fill=Subtype),shape = 21,stroke = 0.2,size = 3)+
  scale_fill_manual(values = as.character(mac.color),breaks = names(mac.color))+
  theme_bw() + theme(legend.position = 'none')+
  theme(
    
    panel.grid=element_blank()
  )+
  labs(x = "Slingshot", y = "Proangiogenic score", color = '',subtitle = paste("Correlation:",corr,sep="")) + geom_smooth(aes(x=pseudoTime,y=gene.module),method = "lm")




ggsave(p,filename = paste("G://embryo/20220930/Slingshot_VS_Proangiogenic_score__plot_with_dots.png",sep=""),width = 5, height = 5, units = "in", device='png')
ggsave(p,filename = paste("G://embryo/20220930/Slingshot_VS_Proangiogenic_score__plot_with_dots.pdf",sep=""),width = 5, height = 5, units = "in", device='pdf')




corr<-cor(meta.info2[,5],meta.info2[,3])
corr<-round(corr,digits=4)



p<-ggplot(meta.info2)+
  geom_point(aes(x=pseudoTime,y=cytotrace,fill=Subtype),shape = 21,stroke = 0.2,size = 3)+
  scale_fill_manual(values = as.character(mac.color),breaks = names(mac.color))+
  theme_bw() + theme(legend.position = 'none')+
  theme(
    
    panel.grid=element_blank()
  )+
  labs(x = "Slingshot", y = "Differentiation capacity", color = '',subtitle = paste("Correlation:",corr,sep="")) + geom_smooth(aes(x=pseudoTime,y=cytotrace),method = "lm")




ggsave(p,filename = paste("G://embryo/20220930/Slingshot_VS_Differentiation_capacity__plot_with_dots.png",sep=""),width = 5, height = 5, units = "in", device='png')
ggsave(p,filename = paste("G://embryo/20220930/Slingshot_VS_Differentiation_capacity__plot_with_dots.pdf",sep=""),width = 5, height = 5, units = "in", device='pdf')


corr<-cor(meta.info2[,5],meta.info2[,7])
corr<-round(corr,digits=4)



p<-ggplot(meta.info2)+
  geom_point(aes(x=pseudoTime,y=Monocle3),color="white",shape = 21,stroke = 0.2,size = 3)+
  theme_bw() + theme(legend.position = 'none')+
  theme(
    
    panel.grid=element_blank()
  )+
  labs(x = "Slingshot", y = "Monocle3", color = '',subtitle = paste("Correlation:",corr,sep="")) + geom_smooth(aes(x=pseudoTime,y=Monocle3),method = "lm")

ggsave(p,filename = paste("G://embryo/20220930/Slingshot_VS_Monocle3_plot.png",sep=""),width = 5, height = 5, units = "in", device='png')
ggsave(p,filename = paste("G://embryo/20220930/Slingshot_VS_Monocle3_plot.pdf",sep=""),width = 5, height = 5, units = "in", device='pdf')

corr<-cor(meta.info2[,5],meta.info2[,2])
corr<-round(corr,digits=4)


p<-ggplot(meta.info2)+
  geom_point(aes(x=pseudoTime,y=Time),color="white",shape = 21,stroke = 0.2,size = 3)+
  theme_bw() + theme(legend.position = 'none')+
  theme(
    
    panel.grid=element_blank()
  )+
  labs(x = "Slingshot", y = "Actual time [Day]", color = '',subtitle = paste("Correlation:",corr,sep="")) + geom_smooth(aes(x=pseudoTime,y=Time),method = "lm")

ggsave(p,filename = paste("G://embryo/20220930/Slingshot_VS_ActualTime_plot.png",sep=""),width = 5, height = 5, units = "in", device='png')
ggsave(p,filename = paste("G://embryo/20220930/Slingshot_VS_ActualTime_plot.pdf",sep=""),width = 5, height = 5, units = "in", device='pdf')

corr<-cor(meta.info2[,5],meta.info2[,6])
corr<-round(corr,digits=4)

p<-ggplot(meta.info2)+
  geom_point(aes(x=pseudoTime,y=gene.module),color="white",shape = 21,stroke = 0.2,size = 3)+
  theme_bw() + theme(legend.position = 'none')+
  theme(
    
    panel.grid=element_blank()
  )+
  labs(x = "Slingshot", y = "Proangiogenic score", color = '',subtitle = paste("Correlation:",corr,sep="")) + geom_smooth(aes(x=pseudoTime,y=gene.module),method = "lm")

ggsave(p,filename = paste("G://embryo/20220930/Slingshot_VS_Proangiogenic_score_plot.png",sep=""),width = 5, height = 5, units = "in", device='png')
ggsave(p,filename = paste("G://embryo/20220930/Slingshot_VS_Proangiogenic_score_plot.pdf",sep=""),width = 5, height = 5, units = "in", device='pdf')

corr<-cor(meta.info2[,5],meta.info2[,4])
corr<-round(corr,digits=4)

p<-ggplot(meta.info2)+
  geom_point(aes(x=pseudoTime,y=cytotrace),color="white",shape = 21,stroke = 0.2,size = 3)+
  theme_bw() + theme(legend.position = 'none')+
  theme(
    
    panel.grid=element_blank()
  )+
  labs(x = "Slingshot", y = "Differentiation capacity", color = '',subtitle = paste("Correlation:",corr,sep="")) + geom_smooth(aes(x=pseudoTime,y=cytotrace),method = "lm")

ggsave(p,filename = paste("G://embryo/20220930/Slingshot_VS_Differentiation_capacity_plot.png",sep=""),width = 5, height = 5, units = "in", device='png')
ggsave(p,filename = paste("G://embryo/20220930/Slingshot_VS_Differentiation_capacity_plot.pdf",sep=""),width = 5, height = 5, units = "in", device='pdf')

