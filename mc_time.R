mcmc1 = metatable
mcmc1 = metatable[metatable$time != "Adult",]
anno.2 = mcmc1[,c("mc","subtype")]
anno.2 = anno.2[!duplicated(anno.2),]
anno.2 = anno.2[anno.2$mc %in% mf_id$mc,]
tmp = paste0("X",anno.2$mc)
rownames(anno.2) = tmp
anno.2 = anno.2[,-1]
anno.2$subtype[anno.2$subtype == ""] = "noname"

timeorder = c("cs11","cs12","cs13","cs14","cs18","cs19","cs21","cs23",
              "week9","week10","week11","week12","week13","week16","week19","week20","week23","week27")
the_order = as.data.frame(timeorder)
colnames(the_order) = c("time")
the_order$the_order = rownames(the_order)

realtime = unique(mcmc1$time)
realtime = as.data.frame(realtime)

realtime$tmp = 0
realtime$tmp[realtime$realtime %in% c("cs11","cs12","cs13","cs10")] = 4
realtime$tmp[realtime$realtime %in% c("cs14","cs15")] = 5
realtime$tmp[realtime$realtime %in% c("cs16","cs17")] = 6
realtime$tmp[realtime$realtime %in% c("cs18","cs19")] = 7
realtime$tmp[realtime$realtime %in% c("cs20","cs21","cs22")] = 8
realtime$tmp[realtime$realtime == "cs23"] = 8.5
realtime$tmp[realtime$realtime == "week9"] = 9
realtime$tmp[realtime$realtime == "week10"] = 10
realtime$tmp[realtime$realtime == "week11"] = 11
realtime$tmp[realtime$realtime == "week12"] = 12
realtime$tmp[realtime$realtime == "week13"] = 13
realtime$tmp[realtime$realtime == "week16"] = 16
realtime$tmp[realtime$realtime == "week19"] = 19
realtime$tmp[realtime$realtime == "week20"] = 20
realtime$tmp[realtime$realtime == "week23"] = 23
realtime$tmp[realtime$realtime == "week27"] = 26
#realtime$tmp[realtime$realtime == "Adult"] = 100
#realtime$tmp[realtime$realtime == "Adult"] = 28

colnames(realtime) = c("time","tmp")
realtime$realtime = realtime$tmp*7
realtime = realtime[,c("time","realtime")]

mcmc1 = merge(mcmc1,realtime,by = "time",all.x = T)


the_mc = c()
the_mean = c()
the_median = c()
for (i in unique(mcmc1$mc)) {
  
  the_mc = c(the_mc,i)
  
  tmp1 = as.data.frame(table(mcmc1$time[mcmc1$mc == i]))
  tmp1$percent = tmp1$Freq/sum(tmp1$Freq)
  colnames(tmp1) = c("time","num","freq")
  
  totalcellnum<-sum(tmp1$num)
  totalcellnum.quantile<-round((totalcellnum/4))
  ##totalcellnum.quantile<-round((totalcellnum/10))
  
  tmp1<-subset(mcmc1,mc == i)
  tmp1 = merge(tmp1,the_order,by = "time",all.x=T)
  tmp1 = tmp1[order(as.numeric(tmp1$the_order),decreasing = F),]
  tmp1.cut = tmp1[1:totalcellnum.quantile,]
  
  the_mean = c(the_mean,mean(tmp1.cut$realtime))
  the_median = c(the_median,median(tmp1.cut$realtime))
}

the_mc =data.frame(the_mc,the_mean,the_median)
colnames(the_mc) = c("mc","the_mean","the_median")

write.csv(the_mc,"E://allcells_realtime_noadult.csv")
