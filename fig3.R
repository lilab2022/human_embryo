##Fig3B 
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
