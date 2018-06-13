##################################
#########  part 3.1  ############# prepare datasets for 15-bp ctcf and randomly selected sequence
##################################

# Package needed: Biostrings, BSgenome.Hsapiens.UCSC.hg19
# Input:
# pwmfn: The path of pwm (PWM in "A","C","G","T" row order)
# fn: The lacation to store the result

# Output:
# all matched ctcf location on Genome
get.ctcf(pwmfn,fn){
  library(Biostrings)
  #### pwm downloaded from your folder
  pwm = t(as.matrix(read.table("CTCF.PPM.txt")))
  rownames(pwm) = c("A","C","G","T")
  
  # get reverse-comlementary seq
  pwm2=reverseComplement(pwm)
  
  #### get matched results
  library(BSgenome.Hsapiens.UCSC.hg19)
  motif.score=NULL
  chr.list=paste0("chr", c(1:22,"X"))
  for(chr in chr.list){
    message(chr)
    seq=Hsapiens[[chr]]
    hit1=matchPWM(pwm, unmasked(seq), with.score=T)
    if(length(hit1)>0){
      s1= mcols(hit1)$score
      s1= (s1-minScore(pwm))/(maxScore(pwm)-minScore(pwm))
      gr1=GRanges(seqnames=chr, ranges=hit1@ranges, strand="+", score=s1)
    }else{
      gr1=GRanges()
    }
    
    hit2=matchPWM(pwm2, unmasked(seq), with.score=T)
    if(length(hit2)>0){
      s2=hit2@elementMetadata$score
      s2= (s2-minScore(pwm2))/(maxScore(pwm2)-minScore(pwm2))
      gr2=GRanges(seqnames=chr, ranges=hit2@ranges, strand="-", score=s2)
    }else{
      gr2=GRanges()
    }
    
    gr=sort(c(gr1,gr2))
    motif.score=rbind(motif.score, as.data.frame(gr))
    motif.ctcf=makeGRangesFromDataFrame(motif.score, keep.extra.columns=T)
  }
  rm(chr,first_time,gr,gr1,gr2,hit1,hit2,s1,s2,seq,motif.score,chr.list)
  save(motif.ctcf, file = paste0(fn,"motif.ctcf.rdata"))
}
get.ctcf(pwmfn="~/Desktop/gm12878/part1/CTCF.ppm.txt",fn="~/Desktop/gm12878/part2/")



# Package needed: Biostrings, BSgenome.Hsapiens.UCSC.hg19
# Input:
# number: The number of regions to sample
# span: The nubmer of bp to add to start
# fn: The lacation to store the result
# ctcffn: The path of ctcf results required from function "get.ctcf()"

# Output:
# random selected sequence set.

library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

get.random(ctcffn,fn,number=10000,span=14){
  #chromosomes of interest
  my_chr <- c(1:22,'X')
  my_chr <- gsub(pattern="^", replacement='chr', my_chr)
  
  #initialise list to store chromosome sizes
  my_chr_size <- list()
  for (i in my_chr){
    my_chr_size[[i]] <- length(BSgenome.Hsapiens.UCSC.hg19[[i]])
  }
  
  #checkout my_chr_size
  head(my_chr_size,2)
  save(my_chr_size,file = paste0(fn,"my_chr_size.rda"))
  
  #initialise some vectors for storing random coordinates
  my_random_start  <- vector()
  my_random_end    <- vector()
  my_random_chr    <- vector()
  my_random_strand <- vector()
  
  #loop through number of regions
  for(i in 1:number){
    my_random_chr[i] <- sample(x=my_chr,size=1)
    my_random_strand[i] <- sample(x=c('-','+'),size=1)
    my_max <- my_chr_size[[my_random_chr[i]]]-span
    my_random_start[i] <- floor(runif(n=1, min=20, max=my_max-19))
    my_random_end[i] <- my_random_start[i] + span
  }
  my_random_width <- rep(span+1,number)
  my_random_seq <- data.frame(seqnames=my_random_chr,
                              start=my_random_start,
                              end=my_random_end,
                              width = my_random_width,
                              strand=my_random_strand)
  head(my_random_seq)
  rm(i,my_chr,my_random_chr,my_random_start,my_random_end,my_random_strand,number,span,my_max,my_random_width)
  random_seq = my_random_seq[order(my_random_seq$seqnames,my_random_seq$start),]
  
  load(ctcffn)
  random_seq1 = makeGRangesFromDataFrame(random_seq, keep.extra.columns=T)
  x = random_seq1 %over% motif.ctcf
  randomseq = random_seq[-which(x),]
  save(randomseq,file = paste0(fn,"randomseq.rda"))
}
get.random(num=10000,span=14,fn="~/Desktop/gm12878/part2/",ctcffn="~/Desktop/gm12878/part2/ctcf.rdata")




##################################
#########  part 3.2  ############# create random datasets for cluster matching
##################################

# Package needed: Biostrings, BSgenome.Hsapiens.UCSC.hg19
# Input:
# x
# span: The nubmer of bp to add to start
# fn: The lacation to store the result
# ctcffn: The path of ctcf results required from function "get.ctcf()"

# Output:
# dataframe with 3 columns: sequences, SVM weights, PWM scores

get.datprep(){
  
}
x = randomseq
#get corresponding 33bps for 15-digit motif
y = cbind(x[,1],x[,2]-9,x[,3]+9,x[,4:5])
colnames(y) = colnames(x)[1:5]
y[,4] = 33
y = makeGRangesFromDataFrame(y, keep.extra.columns=T)
save(y, file = "/Users/firmiana/Desktop/10-mers/GM12878/020918/extend.33random.rdata")
#get the corresponding 19bps for 15-digit motif
z = as.data.frame(y)
for(i in 1:15){
  message(i)
  start.p = cbind(z[,1],z[,2]+i-1,z[,2]+i+17,z[,4:5])
  colnames(start.p) = colnames(z)
  psn = makeGRangesFromDataFrame(start.p, keep.extra.columns=T)
  seq = getSeq(BSgenome.Hsapiens.UCSC.hg19,psn)
  seq_list = strsplit(as.character(seq),"")
  fn  = paste0("/Users/firmiana/Desktop/10-mers/GM12878/020918/Seq.ref/list.motif15/motif",i,".rda")
  save(seq_list,file = fn)
}
rm(i,fn,first_time,seq,psn,start.p,z)


pre28 = cbind(x[,1],x[,2]-19,x[,2]+8,x[,4:5])
colnames(pre28) = colnames(x)[1:5]
pre28[,4] = 28
pre28 = makeGRangesFromDataFrame(pre28, keep.extra.columns=T)
save(pre28, file = "/Users/firmiana/Desktop/10-mers/GM12878/020918/extend.pre28random.rdata")
#get the corresponding 19bps for 10-digit pre-motif
z = as.data.frame(pre28)
for(i in 1:10){
  message(i)
  start.p = cbind(z[,1],z[,2]+i-1,z[,2]+i+17,z[,4:5])
  colnames(start.p) = colnames(z)
  psn = makeGRangesFromDataFrame(start.p, keep.extra.columns=T)
  seq = getSeq(BSgenome.Hsapiens.UCSC.hg19,psn)
  seq_list = strsplit(as.character(seq),"")
  fn = paste0("/Users/firmiana/Desktop/10-mers/GM12878/020918/Seq.ref/list.pre10/pre",i,".rda")
  save(seq_list,file = fn)
}
rm(i,fn,first_time,seq,psn,start.p,z)


aft28 = cbind(x[,1],x[,3]-8,x[,3]+19,x[,4:5])
colnames(aft28) = colnames(x)[1:5]
aft28[,4] = 28
aft28 = makeGRangesFromDataFrame(aft28, keep.extra.columns=T)
save(aft28, file = "/Users/firmiana/Desktop/10-mers/GM12878/020918/extend.aft28random.rdata")
#get the corresponding 19bps for 10-digit aft-motif
z = as.data.frame(aft28)
for(i in 1:10){
  message(i)
  start.p = cbind(z[,1],z[,2]+i-1,z[,2]+i+17,z[,4:5])
  colnames(start.p) = colnames(z)
  psn = makeGRangesFromDataFrame(start.p, keep.extra.columns=T)
  seq = getSeq(BSgenome.Hsapiens.UCSC.hg19,psn)
  seq_list = strsplit(as.character(seq),"")
  fn = paste0("/Users/firmiana/Desktop/10-mers/GM12878/020918/Seq.ref/list.aft10/aft",i,".rda")
  save(seq_list,file = fn)
}
rm(i,fn,first_time,seq,psn,start.p,z)
#####


##################################
#########  part 3.3  ############# random cluster
##################################
fn = paste0("/scratch/bioinfo2/yjin85/GM12878/list.motif15/motif",seq(1:15),".rda")
sn = paste0("/scratch/bioinfo2/yjin85/GM12878/motif/rand.motif",seq(1:15),".alt3.rda")
sn2= paste0("/scratch/bioinfo2/yjin85/GM12878/motif/rand.motif",seq(1:15),".alt3.csv")
load("/scratch/bioinfo2/yjin85/GM12878/d.rdata")
de = as.matrix(d)
dir.create("/scratch/bioinfo2/yjin85/GM12878/motif/")


# calculate sum from 19bps 
seq10_deltasum = function(x){
  xlist = c(paste(x[1:10 ], collapse = ''),paste(x[2:11 ], collapse = ''),paste(x[3:12 ], collapse = ''),paste(x[4:13 ], collapse = ''),paste(x[5:14 ], collapse = ''),
            paste(x[6:15 ], collapse = ''),paste(x[7:16 ], collapse = ''),paste(x[8:17 ], collapse = ''),paste(x[9:18 ], collapse = ''),paste(x[10:19 ], collapse = ''))
  xlistr= chartr("ATGC","TACG",xlist)
  
  m1 = d[match(xlist,de[,1]),2]
  m2 = d[match(xlistr,de[,1]),2]
  m1[which(is.na(m1))]=m2[which(is.na(m1))]
  sum = sum(m1)
  return(sum)
}

seq.delta = function(x){
  dict = c("A","T","G","C")
  temp = dict[which(dict!=x[10])]
  
  x1 = x; x2=x; x3=x
  x1[10]=temp[1]; x2[10]=temp[2]; x3[10]=temp[3]
  vec = NA
  vec[1] = seq10_deltasum(x1)-seq10_deltasum(x)
  vec[2] = seq10_deltasum(x2)-seq10_deltasum(x)
  vec[3] = seq10_deltasum(x3)-seq10_deltasum(x)
  return(vec)
}
# rm(x1,x2,x3,xlist1,xlist2,xlist3,vec,temp,dict)

for(i in 1:15){
  print(i)
  load(fn[i])
  sub_seq = seq_list
  output = lapply(sub_seq,seq.delta)
  out = matrix(unlist(output),ncol = 3,byrow = T)
  save(out,file = sn[i])
  write.csv(out,file = sn2[i],row.names = F)
}



##################################
#########  part 3.4  ############# create ctcf datasets for cluster matching
##################################
#get corresponding 33bps for 15-digit motif
x = as.data.frame(motif.ctcf)
y = cbind(x[,1],x[,2]-9,x[,3]+9,x[,4:5])
colnames(y) = colnames(x)[1:5]
y[,4] = 33
y = makeGRangesFromDataFrame(y, keep.extra.columns=T)
save(y, file = "newextend.33ctcf.rdata")

#get the extended corresponding 28bps for 15-digit motif +/- 10-digit
x = as.data.frame(motif.ctcf)
pre28 = cbind(x[,1],x[,2]-19,x[,2]+8,x[,4:5])
colnames(pre28) = colnames(x)[1:5]
pre28[,4] = 28
pre28 = makeGRangesFromDataFrame(pre28, keep.extra.columns=T)
save(pre28, file = "~/Desktop/10-mers/12878_CTCF/cluster/081617/newextend.pre28ctcf.rdata")

aft28 = cbind(x[,1],x[,3]-8,x[,3]+19,x[,4:5])
colnames(aft28) = colnames(x)[1:5]
aft28[,4] = 28
aft28 = makeGRangesFromDataFrame(aft28, keep.extra.columns=T)
save(aft28, file = "~/Desktop/10-mers/12878_CTCF/cluster/081617/newextend.aft28ctcf.rdata")


#get the corresponding 19bps for 15-digit motif
z = as.data.frame(y)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
for(i in 1:15){
  start.p = cbind(z[,1],z[,2]+i-1,z[,2]+i+17,z[,4:5])
  colnames(start.p) = colnames(z)
  psn = makeGRangesFromDataFrame(start.p, keep.extra.columns=T)
  seq = getSeq(BSgenome.Hsapiens.UCSC.hg19,psn)
  seq_list = strsplit(as.character(seq),"")
  fn  = paste0("~/Desktop/10-mers/12878_CTCF/cluster/081617/Seq.ref/list.motif15/seqs",i,".ref.rda")
  save(seq,file = fn)
}
rm(i,fn,first_time,seq,psn,start.p,z)

#get the corresponding 19bps for 10-digit pre-motif
z = as.data.frame(pre28)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
for(i in 1:10){
  message(i)
  start.p = cbind(z[,1],z[,2]+i-1,z[,2]+i+17,z[,4:5])
  colnames(start.p) = colnames(z)
  psn = makeGRangesFromDataFrame(start.p, keep.extra.columns=T)
  seq = getSeq(BSgenome.Hsapiens.UCSC.hg19,psn)
  seq_list = strsplit(as.character(seq),"")
  fn = paste0("~/Desktop/10-mers/12878_CTCF/cluster/081617/Seq.ref/list.pre10/pre",i,".rda")
  save(seq_list,file = fn)
}
rm(i,fn,first_time,seq,psn,start.p,z)

#get the corresponding 19bps for 10-digit aft-motif
z = as.data.frame(aft28)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
for(i in 1:10){
  message(i)
  start.p = cbind(z[,1],z[,2]+i-1,z[,2]+i+17,z[,4:5])
  colnames(start.p) = colnames(z)
  psn = makeGRangesFromDataFrame(start.p, keep.extra.columns=T)
  seq = getSeq(BSgenome.Hsapiens.UCSC.hg19,psn)
  seq_list = strsplit(as.character(seq),"")
  fn = paste0("~/Desktop/10-mers/12878_CTCF/cluster/081617/Seq.ref/list.aft10/aft",i,".rda")
  save(seq_list,file = fn)
}
rm(i,fn,first_time,seq,psn,start.p,z)




##   get the top sequence (the results are based on pwm matrix and hg19 scanning)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

# start here
# get the corresponding 15bps seq in hg19
seq15.hg19= getSeq(BSgenome.Hsapiens.UCSC.hg19,motif.ctcf)
# save(seq15.hg19, file = "/Users/firmiana/Desktop/10-mers/12878_CTCF/cluster/090617/sequence15.hg19.rda")

seq15 = as.vector(seq15.hg19)
length(unique(seq15))  # 24402
seq15 = sort(seq15)
name = unique(seq15)
freq = as.numeric(table(seq15))
count.freq15 = cbind(name,freq)
count.freq15 = count.freq15[order(-freq),]
save(count.freq15,file = "/Users/firmiana/Desktop/10-mers/12878_CTCF/cluster/090617/count.freq15.rda")

load("/Users/firmiana/Desktop/10-mers/12878_CTCF/cluster/090617/count.freq15.rda")
load("/Users/firmiana/Desktop/10-mers/12878_CTCF/cluster/090617/sequence15.hg19.rda")

seq15 = as.vector(seq15.hg19)
dir.create("/Users/firmiana/Desktop/10-mers/GM12878/011618/toprow/")
toprow = paste0("/Users/firmiana/Desktop/10-mers/GM12878/011618/toprow/top.row",seq(1:20),".rdata")
for(i in 1:20){
  x = count.freq15[i,1]
  top.row = which(seq15 == x)
  save(top.row,file = toprow[i])
}


# to check the amount of show-up time for corresponding complementary sequence (always zero)
for(i in 1:20){
  x = count.freq15[i,1]
  y = chartr("ATGC","TACG",x)
  # y = reverse(chartr("ATGC","TACG",x))
  m1 = length(which(seq15 == x))
  m2 = length(which(seq15 == y))
  m = cbind(m1,m2)
  print(m)
}




##################################
#########  part 3.5  ############# ctcf cluster
##################################
dir.create(paste0("/scratch/bioinfo2/yjin85/hg19/results.top",num,"/"))
dir.create(paste0("/scratch/bioinfo2/yjin85/hg19/results.top",num,"/pre/"))
dir.create(paste0("/scratch/bioinfo2/yjin85/hg19/results.top",num,"/aft/"))
dir.create(paste0("/scratch/bioinfo2/yjin85/hg19/results.top",num,"/motif/"))
load("/scratch/bioinfo2/yjin85/GM12878/d.rdata")
load(paste0("/scratch/bioinfo2/yjin85/GM12878/toprow/top.row",num,".rdata"))
de = as.matrix(d)


# calculate sum from 19bps 
seq10_deltasum = function(x){
  xlist = c(paste(x[1:10 ], collapse = ''),paste(x[2:11 ], collapse = ''),paste(x[3:12 ], collapse = ''),paste(x[4:13 ], collapse = ''),paste(x[5:14 ], collapse = ''),
            paste(x[6:15 ], collapse = ''),paste(x[7:16 ], collapse = ''),paste(x[8:17 ], collapse = ''),paste(x[9:18 ], collapse = ''),paste(x[10:19 ], collapse = ''))
  xlistr= chartr("ATGC","TACG",xlist)
  
  m1 = d[match(xlist,de[,1]),2]
  m2 = d[match(xlistr,de[,1]),2]
  m1[which(is.na(m1))]=m2[which(is.na(m1))]
  sum = sum(m1)
  return(sum)
}

seq.delta = function(x){
  dict = c("A","T","G","C")
  temp = dict[which(dict!=x[10])]
  
  x1 = x; x2=x; x3=x
  x1[10]=temp[1]; x2[10]=temp[2]; x3[10]=temp[3]
  vec = NA
  vec[1] = seq10_deltasum(x1)-seq10_deltasum(x)
  vec[2] = seq10_deltasum(x2)-seq10_deltasum(x)
  vec[3] = seq10_deltasum(x3)-seq10_deltasum(x)
  return(vec)
}
# rm(x1,x2,x3,xlist1,xlist2,xlist3,vec,temp,dict)

fn = paste0("/scratch/bioinfo2/yjin85/hg19/data/list.pre10/pre",seq(1:10),".rda")
sn = paste0("/scratch/bioinfo2/yjin85/hg19/results.top",num,"/pre/top",num,"pre",seq(1:10),".alt3.rda")
sn2= paste0("/scratch/bioinfo2/yjin85/hg19/results.top",num,"/pre/top",num,"pre",seq(1:10),".alt3.csv")
for(i in 1:10){
  print(i)
  load(fn[i])
  sub_seq = seq_list[top.row]
  output = lapply(sub_seq,seq.delta)
  out = matrix(unlist(output),ncol = 3,byrow = T)
  save(out,file = sn[i])
  write.csv(out,file = sn2[i],row.names = F)
}


fn = paste0("/scratch/bioinfo2/yjin85/hg19/data/list.motif15/seqs",seq(1:15),".ref.rda")
sn = paste0("/scratch/bioinfo2/yjin85/hg19/results.top",num,"/motif/top",num,"seqs",seq(1:15),".alt3.rda")
sn2= paste0("/scratch/bioinfo2/yjin85/hg19/results.top",num,"/motif/top",num,"seqs",seq(1:15),".alt3.csv")
for(i in 1:15){
  print(i)
  load(fn[i])
  sub_seq = seq_list[top.row]
  output = lapply(sub_seq,seq.delta)
  out = matrix(unlist(output),ncol = 3,byrow = T)
  save(out,file = sn[i])
  write.csv(out,file = sn2[i],row.names = F)
}


fn = paste0("/scratch/bioinfo2/yjin85/hg19/data/list.aft10/aft",seq(1:10),".rda")
sn = paste0("/scratch/bioinfo2/yjin85/hg19/results.top",num,"/aft/top",num,"aft",seq(1:10),".alt3.rda")
sn2= paste0("/scratch/bioinfo2/yjin85/hg19/results.top",num,"/aft/top",num,"aft",seq(1:10),".alt3.csv")
for(i in 1:10){
  print(i)
  load(fn[i])
  sub_seq = seq_list[top.row]
  output = lapply(sub_seq,seq.delta)
  out = matrix(unlist(output),ncol = 3,byrow = T)
  save(out,file = sn[i])
  write.csv(out,file = sn2[i],row.names = F)
}
