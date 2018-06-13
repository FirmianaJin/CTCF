##################################
###########  part 6  ############# sequence logo plot
##################################
rm(list = ls())
#### pwm downloaded from Tianlei's folder
pwm = t(as.matrix(read.table("/Users/firmiana/Desktop/10-mers/12878_CTCF/cluster/081617/CTCF.ppm.txt")))
rownames(pwm) = c("A","C","G","T")

vec = strsplit(as.character(count.freq15[1:20,1]),"")
# x=vec[[1]]
seq_ex = function(x){
  data_seq = matrix(rep(x,45),ncol = 15,byrow = T)
  dict = c("A","T","G","C")
  for (i in 1:15) {
    temp <- dict[which(dict!=x[i])]
    for (j in 1:3) {
      data_seq[(i-1)*3+j,i] = temp[j]
    }
  }
  vec=rep(NA,45)
  for(ii in 1:45){
    vec[ii] = paste(data_seq[ii,], collapse = '')
  }
  return(vec)
}
m = lapply(vec,seq_ex)

# create single variant sequence matrix (30 variant in one row)------ n
n = matrix(unlist(m),ncol = 45,byrow = T)
rm(seq_ex,m)

logsum = function(topi){
  m = diag(pwm[topi,1:15])
  log_sum = sum(log(m+0.00001))
  return(log_sum)
}
seq_list = strsplit(as.character(as.vector(n)),"")
value1 = lapply(seq_list,logsum)
value1 = matrix(unlist(value1),byrow = F,ncol = 45)
valueo= lapply(vec,logsum)
valueo= matrix(rep(unlist(valueo),45),byrow = F,ncol = 45)
value = value1-valueo


rm(seq_list,logsum,valueo,value1)

val = matrix(rep(NA,300),ncol = 15,byrow = T)
for(i in 1:15){
  start = 3*i-2
  stop = 3*i
  val[,i] = apply(value[,start:stop],1,mean)
}
rm(i,start,stop)
save(val, file = "/Users/firmiana/Desktop/10-mers/GM12878/022018/val.rda")


m = get.cre()  # threshold got from three cell lines
###   prediction table
tt = matrix(0,ncol = 15, nrow = 20)
for (i in 1:20){
  x = top.sum.all[[i]][11:25,4]
  y = val[i,]
  tt[i,which(x<m)] = 1
};rm(i,x,y)

#### pwm downloaded from Tianlei's folder
pwm = t(as.matrix(read.table("/Users/firmiana/Desktop/10-mers/12878_CTCF/cluster/081617/CTCF.ppm.txt")))
rownames(pwm) = c("A","C","G","T")

load("count.freq15.rda")
z = count.freq15[1:20,1:2]
vec = strsplit(as.character(count.freq15[1:20,1]),"")
logsum = function(topi){
  m = diag(pwm[topi,1:15])
  log_sum = sum(log(m+0.00001))
  return(log_sum)
}
valueo= lapply(vec,logsum)
tt = cbind(z,round(unlist(valueo),digits = 3),tt)
rm(valueo,z,pwm,logsum,vec)
colnames(tt) = c("seq","freq","logsum",1:15)
write.csv(tt,file = "/Users/firmiana/Desktop/10-mers/GM12878/022018/tt.means.csv",row.names = F)

# mean
tt = read.csv("/Users/firmiana/Desktop/10-mers/GM12878/022018/tt.means.csv",header = T)
library(ggseqlogo)
library(ggplot2)
tt.re = tt[,4:18]
tt.re[tt.re==0]<-0.3
seqlist = strsplit(as.character(tt[,1]),"")

pred.logo = list()
for(j in 1: length(seqlist)){  #length(seqlist)
  tcga = matrix(0,ncol = 15,nrow = 4);rownames(tcga) = c("A","C","G","T")
  for(i in 1:15){
    tcga[seqlist[[j]][i],i]=tt.re[j,i]
  }
  pred.logo[[j]]=tcga
}
rm(i,j,tcga,seqlist,tt.re)
# ggseqlogo(pred.logo, ncol=4)


# May 29 revised
pdf("/Users/firmiana/Desktop/10-mers/GM12878/paper/Paper004.pdf",width = 12)
ggplot() +
  geom_logo(pred.logo,method='custom')+
  # annotate('segment', x = 0.5, xend=15.5, y=0.9, yend=0.9, size=1.5,color="black") +
  theme_logo()+
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=10),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  facet_wrap(~seq_group, nrow = 5, scales='free_x')+
  ylab('Prediction')+
  xlab('Top Sequence')+
  theme(
    axis.title.x=element_text(size=20),
    axis.title.y=element_text(size=20) #face='bold', 
  )
dev.off()
