##################################
###########  part 2  ############# summary of top sequences
##################################

# Package needed: None
# Input:
# dfn: The path of results from function "seq.score()" (.rda)
# fn: The lacation to store the result

# Output:
# histogram
# brief summary
# return:  correlation between deltaPWM and deltaSVM for top 1000 sequence

summary.var30(dfn,fn){
  ### find top 1000 sequence
  load(dfn)
  de = as.matrix(d)
  # top 1000 sequence based on probability
  num = order(d[,3], decreasing = T)[1:1000]
  sub_de = d[num,]
  save(sub_de,file = "/Users/firmiana/Desktop/10-mers/GM12878/011518/sub_de_prob.rdata")
  
  # top 1000 sequence based on gkmSVM value
  num = order(d[,2],decreasing = T)[1:1000]
  sub_de_gkm = de[num,]
  save(sub_de_gkm,file = "/Users/firmiana/Desktop/10-mers/GM12878/011518/sub_de_gkm.rdata")
  
  sub_de = sub_de_gkm
  seq_list = strsplit(as.character(sub_de[,1]),"")
  
  # create single variant sequence
  seq_ex = function(x){
    data_seq = matrix(rep(x,30),ncol = 10,byrow = T)
    dict = c("A","T","G","C")
    for (i in 1:10) {
      temp <- dict[which(dict!=x[i])]
      for (j in 1:3) {
        data_seq[(i-1)*3+j,i] = temp[j]
      }
    }
    vec=rep(NA,30)
    for(ii in 1:30){
      vec[ii] = paste(data_seq[ii,], collapse = '')
    }
    return(vec)
  }
  m = lapply(seq_list,seq_ex)
  
  # create single variant sequence matrix (30 variant in one row)------ n
  n = matrix(unlist(m),ncol = 30,byrow = T)
  rm(num,seq_ex,m,seq_list)
  
  # create single variant reverse comlementary sequence------ n_c
  n_c1 = chartr("ATGC","TACG",n)
  n_c = n_c1
  for(i in 1:1000){
    for(j in 1:30){
      n_c[i,j]=paste(rev(strsplit(n_c1[i,j],"")[[1]]),collapse = '')
    }
  }
  remove(n_c1,i,j)
  
  # sub_de refers to top 1000 tenmers ranked by gkmSVM weights
  test1 = as.vector(t(n))
  test2 = as.vector(t(n_c))
  m1 = match(test1,d[,1])
  m2 = match(test2,d[,1])
  mm = m1
  mm[which(is.na(mm))]=m2[which(is.na(mm))]
  rm(test1,test2,m1,m2)
  
  # matched sequence for changed sequence
  x.seq = matrix(d[mm,1],ncol = 30,byrow = T)
  # original gkmSVM value for changed sequence
  dist.ori = matrix(d[mm,2],ncol = 30,byrow = T)
  # difference between changed sequence and original seq
  dist = matrix(NA,ncol = 30,nrow = 1000)
  for (i in 1:1000){
    for(j in 1:30){
      dist[i,j]=dist.ori[i,j]-as.numeric(sub_de_gkm[i,2])
    }
  }
  rm(i,j)
  
  # original probability for changed sequence based on probability
  dist.p.ori = matrix(d[mm,3],ncol = 30,byrow = T)
  # difference probability between changed sequence and original seq
  dist.p = matrix(NA,ncol = 30,nrow = 1000)
  for (i in 1:1000){
    for(j in 1:30){
      dist.p[i,j]=exp(dist.p.ori[i,j])-exp(as.numeric(sub_de_gkm[i,3]))
    }
  }
  rm(i,j)
  
  # correlation between deltaPWM and deltaSVM for top 1000 sequence
  g = as.vector(t(dist))
  g.p = as.vector(t(dist.p))
  cr=cor(g,g.p)  # 0.5067374
  
  # histogram for all 30000 changed sequence
  g = unique(as.numeric(as.vector(t(dist))))
  jpeg(file=paste0(fn,"Distribution of single variant.jpeg"))
  hist(g,main = paste0("distribution of ",length(g)," single variant"),xlab = "difference")
  dev.off()
  
  # summary for all 30000 changed sequence
  sink(file = paste0(fn,"summary of ",length(g)," single variant"))
  summary(g)
  sink()
  return(cr)
}

summary.var30(dfn="~/Desktop/gm12878/part1/d.rdata",fn="~/Desktop/gm12878/part2/")
