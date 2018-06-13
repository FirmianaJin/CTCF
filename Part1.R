##################################
###########  part 1  ############# Dateset preparation
##################################

# Package needed: seqinr, gkmSVM
# Input:
# wd: The path to store the posSet and negSet
# dat: The name of original narrow peak ".bed" file downloaded from ENCODE
# n: The number of posSet/negSet to be used in training steps

# Output:
# ".bed" file
# posSet and negSet
# 10000 randomly selected posSet and negSet (Notice: We use 10000 in our method, it can be other numbers)
get.set(wd,dat,n)={
  setwd(wd)
  bed <- as.data.frame(read.table(paste0(dat,".bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
  bed1 = bed[,1:3]
  bed1 = bed1[order(bed1[,1]),]
  write.table(bed1, file=paste0(dat,"test.bed"), quote=F, sep="\t", row.names=F, col.names=F)
  # get original corresponding number of posSet and negSet
  library(gkmSVM)
  inputBedFN = paste0(wd,dat,"test.bed")
  genNullSeqs(inputBedFN, genomeVersion = 'hg19')
  
  library(seqinr)
  test = read.fasta("posSet.fa")
  # pos
  int = sort(sample.int(length(test),10000,replace = F))
  x = test[int]
  write.fasta(sequences = x, names = names(x), nbchar = 238, file.out = paste0("pos",n,".fa"))
  # neg
  test = read.fasta("negSet.fa")
  int = sort(sample.int(length(test),10000,replace = F))
  y = test[int]
  write.fasta(sequences = y, names = names(y), nbchar = 238, file.out = paste0("neg",n,".fa"))
}

get.set("~/Desktop/gm12878/part1/","ENCFF002DAJ",10000)



# Package needed: gkmSVM
# Input:
# posname: The name of positive set (.fa file)
# nfn: The name of negative set (.fa file)
# fn: The lacation to retrieve the training set and store the result
# testfn: The test sequence file (all 10-mers included)

# Output:
# Kernel file
# 2 CTCF_svmtrain files
# output file with scores for each test sequence
gkm.train(fn,posname,negname,testfn,outname){
  setwd(fn)
  library(gkmSVM)
  # computes kernel
  posfn    = paste0(posname,".fa")
  negfn    = paste0(negname,".fa")
  kernelfn = 'CTCFkernel.txt'
  gkmsvm_kernel(posfn, negfn, kernelfn)
  svmfnprfx = 'ctcf_svmtrain'
  gkmsvm_train(kernelfn,posfn,negfn,svmfnprfx)
  
  svmfnprfx = 'ctcf_svmtrain'
  outfn_nr10 = paste0(outname,".txt")
  gkmsvm_classify(testfn, svmfnprfx, outfn_nr10)
}

gkm.train("~/Desktop/gm12878/part1/","pos10000","neg10000","nr10mers.fa.txt","12878w_nr10")



# Package needed: None
# Input:
# vec: The splited 15-bps sequence
# prob1: The log-transformed PWM for vec
# prob2: The log-transformed PWM for complimentary vec

# Output:
# the largest log-transformed PWM value while matching target sequence to PWM

### 10-mers probability calculation (based on PWM)
tenmers_prob = function(vec,prob1=pr1,prob2=pr3){
  log_sum = rep(NA,12)
  for(ii in 1:6){
    num_s = ii
    num_e = ii+9
    m = diag(prob1[vec,num_s:num_e])
    log_sum[ii] = sum(m)
  }
  
  for(ii in 1:6){
    num_s = ii
    num_e = ii+9
    m = diag(prob2[vec,num_s:num_e])
    log_sum[ii+6] = sum(m)
  }  
  prob_max = max(log_sum)
  return(prob_max)
}

# Package needed: None
# Input:
# pwmfn: The path of pwm (PWM in "A","C","G","T" row order)
# outputfn: The path of output file containing scores for each test sequence
# fn: The lacation to store the result

# Output:
# dataframe with 3 columns: sequences, SVM weights, PWM scores

seq.score(pwmfn,outputfn,fn){
  pwm = t(as.matrix(read.table(pwnfn)))
  rownames(pwm) = c("A","C","G","T")
  prob1 = pwm
  colnames(prob1)=seq(1,15,1)
  rownames(prob1)=c("A","C","G","T")
  pr1=log(prob1)
  # complimentary one
  pr2=pr1[4:1,]
  rownames(pr2)=c("A","C","G","T")
  pr3=pr2[,15:1]
  rm(pr2,pwm,prob1)
  
  ### read in output file --> (tenmer + weight)
  dataset = read.table(outputfn, sep="\t")
  data=strsplit(as.character(dataset[,1]),"")
  
  g = lapply(data,tenmers_prob)
  test1 = unlist(g)
  test1[which(test1==-Inf)]=floor(min(test1[which(test1!=-Inf)]))
  # save(test1,file = "test1.rdata")
  d = cbind(dataset,test1)
  save(d,file = paste0(fn,"d.rda"))
  return(d)
}

seq.score("~/Desktop/gm12878/part1/CTCF.ppm.txt","~/Desktop/gm12878/part1/12878w_nr10.txt",
          "~/Desktop/gm12878/part1/")


### correlation between SVM weights and PWM score (for all 10-mers)
cor(d$V2,d$test1)

### correlation between SVM weights and PWM score (for top 1000 10-mers)
test1 = d[,3];dataset = d[,1:2]
num = test1[order(test1,decreasing = T)][1000]
cor(dataset$V2[test1>num],test1[test1>num])
# plot(dataset$V2[test1>num]~test1[test1>num],pch=20)
