##################################
###########  part 4  ############# threshold & trend plot
##################################

# Package needed: None
# Input:
# fn1: The path of first cell line pre/motif/aft sequences
# fn2: The path of second cell line pre/motif/aft sequences
# fn3: The path of third cell line pre/motif/aft sequences
# Output:
# return:  threshold for prediction based on three cell line

#####  random seq
get.cre(fn1,fn2,fn3){
  sum.bps = function(fn,cell,loc1,loc2,num){
    fx = paste0(fn,loc,"/rand.",loc,seq(1:num),".alt3.rda")
    top1 = NULL
    for(j in 1:num){
      load(fx[j])
      top1 = cbind(top1,as.vector(out))
    }
    return(top1)
  }
  
  cell.all = function(x){
    r1 = sum.bps(fn = "/Users/firmiana/Desktop/10-mers/",cell=x,loc = "pre",num = 10)
    r2 = sum.bps(fn = "/Users/firmiana/Desktop/10-mers/",cell=x,loc = "motif",num = 15)
    r3 = sum.bps(fn = "/Users/firmiana/Desktop/10-mers/",cell=x,loc = "aft",num = 10)
    rand.all = cbind(r1,r2,r3)
    return(rand.all)
  }
  
  rand.gm  = cell.all("GM12878")
  rand.h1  = cell.all("H1-hESC")
  rand.k562= cell.all("K562")
  cre = quantile(c(as.vector(rand.gm),as.vector(rand.h1),as.vector(rand.k562)),probs=0.05,na.rm = T)
  return(cre)
}
cre = get.cre()

pre.plot = function(dist,num=100){
  s = dist[1:num,]
  h = as.numeric(as.vector(t(s)))
  a = rep(1:num,each=30)
  m = as.data.frame(cbind(h,a))
  colnames(m) = c("delta","top_seq")
  m$top_seq = as.factor(m$top_seq)
  mn = apply(s,1,mean)
  se = apply(s,1,sd)
  # pre1 = qnorm(0.975,mean = mn, sd = se)
  # pre2 = qnorm(0.025,mean = mn, sd = se)
  pre1 = mn+se
  pre2 = mn-se
  
  plot(mn~seq(1:num),type="l",pch=20,ylim=c(-7,1.7),cex.lab = 2,
       xlab="Sequence Rank",ylab="Weight Difference") #,main=paste0("The Mean Change of Weights For Top ",num," Motif Sequecnes")
  # points(mn);points(pre1);points(pre2)
  lines(pre1,lty=3,col="grey50");lines(pre2,lty=3,col="grey50")
  abline(h=cre,col="red",lty=4,lwd=2)
}
pre.plot(dist,100)

setwd("~/Desktop/10-mers/GM12878/paper")
pdf(file = "Paper%03d.pdf",onefile = F,width = 12,height = 6)
pre.plot(dist,10)
pre.plot(dist,100)
pre.plot(dist,1000)
dev.off()
