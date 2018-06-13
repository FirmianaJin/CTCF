##################################
###########  part 5  ############# summary of each position 
##################################
summary.bps = function(loc1,loc2,num){
  top.loc = list()
  for(i in 1:20){
    fx = paste0("/Users/firmiana/Desktop/10-mers/GM12878/011618/results.top",i,"/",loc1,"/top",i,loc2,seq(1:num),".alt3.rda")
    top1 = NULL
    for(j in 1:num){
      load(fx[j])
      top1 = rbind(top1,as.vector(summary(as.vector(out))))
    }
    colnames(top1) = c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.")
    top.loc[[i]] = top1
  }
  return(top.loc)
}
# loc1="motif";loc2="seqs";num=15
top.pre.all =summary.bps(loc1 = "pre",loc2 = "pre",num = 10)
top.motif.all = summary.bps(loc1 = "motif",loc2 = "seqs",num = 15)
top.aft.all = summary.bps(loc1 = "pre",loc2 = "pre",num = 10)
top.sum.all = list();for(i in 1:20){x = rbind(top.pre.all[[i]],top.motif.all[[i]],top.aft.all[[i]]);top.sum.all[[i]]=x}
save(top.pre.all,file = "/Users/firmiana/Desktop/10-mers/GM12878/011818/top.pre.all.rda")
save(top.motif.all,file="/Users/firmiana/Desktop/10-mers/GM12878/011818/top.motif.all.rda")
save(top.aft.all,file = "/Users/firmiana/Desktop/10-mers/GM12878/011818/top.aft.all.rda")
save(top.sum.all,file = "/Users/firmiana/Desktop/10-mers/GM12878/011818/top.sum.all.rda")

save.summary1 = function(x,loc2){
  for(i in 1:20){
    dir.create(paste0("/Users/firmiana/Desktop/10-mers/GM12878/011618/results.top",i,"/summary/"))
    sn.pre = paste0("/Users/firmiana/Desktop/10-mers/GM12878/011618/results.top",i,"/summary/top",i,".",loc2,".rda")
    wn.pre = paste0("/Users/firmiana/Desktop/10-mers/GM12878/011618/results.top",i,"/summary/top",i,".",loc2,".csv")
    # sn.pre = paste0("/Users/firmiana/Desktop/10-mers/GM12878/011818/top",i,".",loc2,".rda")
    # wn.pre = paste0("/Users/firmiana/Desktop/10-mers/GM12878/011818/top",i,".",loc2,".csv")
    top.loc = x[[i]];save(top.loc,file = sn.pre)
    top.loc = x[[i]];write.csv(top.loc,file = wn.pre,row.names = F)
  }
}

save.summary1(top.pre.all,loc2 = "pre")
save.summary1(top.aft.all,loc2 = "aft")
save.summary1(top.motif.all,loc2 = "motif")
save.summary1(top.sum.all,loc2 = "sum")



# get deltaSVM value matrix
get.delta.bps = function(loc1,loc2,num){
  top.loc = list()
  for(i in 1:20){
    fx = paste0("/Users/firmiana/Desktop/10-mers/GM12878/011618/results.top",i,"/",loc1,"/top",i,loc2,seq(1:num),".alt3.rda")
    top1.loc = NULL
    for(j in 1:num){
      load(fx[j])
      top1.loc = cbind(top1.loc,out)
    }
    top.loc[[i]] = top1.loc
  }
  return(top.loc)
}

top.loc10pre = get.delta.bps(loc1 = "pre",loc2 = "pre",num = 10)
top.loc15motif = get.delta.bps(loc1 = "motif",loc2 = "seqs",num = 15)
top.loc10aft = get.delta.bps(loc1 = "pre",loc2 = "pre",num = 10)
top.loc35 = list();for(i in 1:20){x = cbind(top.loc10pre[[i]],top.loc15motif[[i]],top.loc10aft[[i]]);top.loc35[[i]]=x}
save(top.loc10pre,file = "/Users/firmiana/Desktop/10-mers/GM12878/011818/top.loc10pre.rda")
save(top.loc15motif,file = "/Users/firmiana/Desktop/10-mers/GM12878/011818/top.loc15motif.rda")
save(top.loc10aft,file = "/Users/firmiana/Desktop/10-mers/GM12878/011818/top.loc10aft.rda")
save(top.loc35,file = "/Users/firmiana/Desktop/10-mers/GM12878/011818/top.loc35.rda")
