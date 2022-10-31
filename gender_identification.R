#####gender

head(metatable)


gtf = read.table("E:\\tmp/gene_intervals_hg38P.ens90.txt",header = T)

ygene = c("KDM5D","TTTY15","EIF1AY","TTTY14","CDD24P4","ZFY","PRKY","GYG2P1","PSMA6P1","USP9Y","TXLNGY",
          "UTY","RPS4Y1","DDX3Y","EIF1AY")
xgene = c("XIST", "PSPHP1")

tmpy = as.matrix(a_reduced[intersect(rownames(a_reduced),ygene),metatable$Well_ID])
tmpy = as.data.frame(colSums(tmpy))
colnames(tmpy) = "y"
tmpy$Well_ID = rownames(tmpy)
 
tmpx = as.matrix(a_reduced[intersect(rownames(a_reduced),xgene),metatable$Well_ID])
tmpx = as.data.frame(colSums(tmpx))
colnames(tmpx) = "x"
tmpx$Well_ID = rownames(tmpx)

tmp = merge(tmpx,tmpy,by = "Well_ID")
tmp = merge(tmp,metatable,by = 'Well_ID')
tmp = tmp[,c("x","y","embryo")]
tmp = tmp[tmp$embryo != "Adult",]
tmp2 = reshape2::melt(tmp)
tmp2$embryo = factor(tmp2$embryo,levels = c("embryo 1","embryo 2","embryo 4","embryo 5","embryo 6",
                                            "embryo 7","embryo 8",
                                            "embryo 9","embryo 10","embryo 11","embryo 12","embryo 13","embryo 14",
                                            "embryo 15","embryo 16","embryo 17","embryo 18","embryo 19","embryo 20",
                                            "embryo 21","embryo 22","embryo 23","embryo 24","embryo 25","embryo 26",
                                            "embryo 27",
                                            "embryo 28","embryo 29","embryo 30","embryo 31","embryo 33",
                                            "embryo 32","embryo 34","embryo 35","embryo 36","embryo 37","embryo 38",
                                            "embryo 39","embryo 40","embryo 41","embryo 43",
                                            "embryo 44","embryo 45","embryo 46","embryo 47","embryo 50","embryo 52",
                                            "embryo 54","embryo 55","embryo 56","embryo57","embryo 58"),ordered = T)
ggplot(tmp2,aes(x = embryo, y = value,fill = variable))+
  geom_violin(scale = "width",trim = T)+
  theme_bw()+theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))


#######
##ferm

tmpy = as.matrix(a_reduced[intersect(rownames(a_reduced),ygene),metatable$Well_ID[metatable$tissue == "Skin" & metatable$subtype == "microglia"]])
tmpy = as.data.frame(colSums(tmpy))
colnames(tmpy) = "y"
tmpy$Well_ID = rownames(tmpy)

tmpx = as.matrix(a_reduced[intersect(rownames(a_reduced),xgene),metatable$Well_ID[metatable$tissue == "Skin" & metatable$subtype == "microglia"]])
tmpx = as.data.frame(colSums(tmpx))
colnames(tmpx) = "x"
tmpx$Well_ID = rownames(tmpx)

tmp = merge(tmpx,tmpy,by = "Well_ID")
tmp = merge(tmp,metatable,by = 'Well_ID')
tmp = tmp[,c("x","y","embryo")]
tmp = tmp[tmp$embryo != "Adult",]
tmp2 = reshape2::melt(tmp)
tmp2$embryo = factor(tmp2$embryo,levels = c("embryo 1","embryo 2","embryo 4","embryo 5","embryo 6",
                                            "embryo 7","embryo 8",
                                            "embryo 9","embryo 10","embryo 11","embryo 12","embryo 13","embryo 14",
                                            "embryo 15","embryo 16","embryo 17","embryo 18","embryo 19","embryo 20",
                                            "embryo 21","embryo 22","embryo 23","embryo 24","embryo 25","embryo 26",
                                            "embryo 27",
                                            "embryo 28","embryo 29","embryo 30","embryo 31","embryo 33",
                                            "embryo 32","embryo 34","embryo 35","embryo 36","embryo 37","embryo 38",
                                            "embryo 39","embryo 40","embryo 41","embryo 43",
                                            "embryo 44","embryo 45","embryo 46","embryo 47","embryo 50","embryo 52",
                                            "embryo 54","embryo 55","embryo 56","embryo57","embryo 58"),ordered = T)
ggplot(tmp2,aes(x = embryo, y = value,fill = variable))+
  geom_violin(scale = "width",trim = T)+
  theme_bw()+theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
