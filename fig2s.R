

for (i in unique(metatable$major)){
  the_metatable = metatable[metatable$major == i,]
  for (t in unique(the_metatable$time)) {
    if (nrow(the_metatable[the_metatable$time == t,])<10){
      the_metatable= the_metatable[-which(the_metatable$time == t),]
    }
  }
  
  tmp = the_metatable
  tmp= tmp[tmp$time != "Adult",]
  
  tmp2 = tmp
  
  
  tmp = tmp2[,c("tissue","time","subtype")]
  tmp = tmp[,c("time","subtype")]
  other = setdiff(time_order,c(unique(tmp$time),"Adult"))
  other = as.data.frame(other)
  other$blank = "blank"
  colnames(other) = colnames(tmp)
  tmp = rbind(tmp,other)
  tmp$time <- factor(tmp$time,levels=time_order,ordered = TRUE)
  p2 = ggplot(tmp,aes(x= time,fill = subtype))+
    geom_bar(stat="count",position = "fill",width = 0.8,color = "black")+
    xlab("")+ylab("")+
    scale_fill_manual(values = c(as.character(subcolor,mac.color),"white"),breaks = c(names(c(subcolor)),"blank"))+
    scale_y_continuous(expand=c(0,0))+
    theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
          axis.text = element_text(size = 12),
          axis.title=element_blank(),
          legend.title=element_blank(),
          panel.background=element_blank(),
          axis.line=element_line(),
          axis.ticks=element_line(),
          legend.position = "none")
  ggsave(paste("E:\\fig1/sub_dynamic/",i,"_sub_dynamic.pdf",sep = ""),p2,width = 8,height = 5,dpi = 300)
}
