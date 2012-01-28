source('new_points_common.R')
load('preprocess.RData')

# simulated annealing
n <- 25
tmax <- 1000
dcoords.new <- NULL
num.children <- 12
parallelize <- TRUE

# first argument is a period-separated list of indices into the candidates dataframe
e <- commandArgs(TRUE)[1]
e <- as.numeric(unlist(strsplit(e, "\\.")))
de <- candidates[e,]
runid <- commandArgs(TRUE)[2]
e_before <- e

t <- tmax
if(parallelize){
  c1 <- makeForkCluster(num.children)
  kv <- krige.var.par(rbind(dcoords,de[,c("x","y")]),loci,kc,c1)
}else{
  kv <- krige.var(rbind(dcoords,de[,c("x","y")]),loci,kc)
}
vmap <- flipud(matrix(kv,nrow=height,ncol=width,byrow=TRUE))
fitness <- wpe(rmap,vmap)
fitness2 <- mean(sqrt(vmap))
first.fitness <- fitness
first.fitness2 <- fitness2
rm(kv,vmap)

linear.cooling = FALSE

log <- NULL

while(t > 0){
  e2 <- e[sample(seq(1,n),n-1)] # n - 1 sized sample of indices
  de2 <- candidates[e2,]
  while(length(e2) < n){
    p <- sample(seq(1,nrow(candidates)),1)
    if(any(e2 == p)) next
    e2 <- append(e2,p)
    de2 <- rbind(de2,candidates[p,])
  }
  if(parallelize){
    kv <- krige.var.par(rbind(dcoords,de2[,c("x","y")]),loci,kc,c1)
  }else{
    kv <- krige.var(rbind(dcoords,de2[,c("x","y")]),loci,kc)
  }
  vmap <- flipud(matrix(kv,nrow=height,ncol=width,byrow=TRUE))
  new.fitness <- wpe(rmap,vmap)
  new.fitness2 <- mean(sqrt(vmap))
  rm(vmap,kv)
  
  replaced <- TRUE
  deltaf <- new.fitness-fitness
  p <- NA
  if(deltaf < 0){
    fitness <- new.fitness
    fitness2 <- new.fitness2
    e <- e2
    de <- de2
  }else{
    if(linear.cooling){
      p <- t/tmax
    }else{
      temp <- t/tmax
      p <- exp(-100.0*deltaf/temp) # scale up deltaf by 10^2 to get a more meaningful cooling curve
    }
    print(paste(t,"worse :(",fitness,fitness2,deltaf,p))
    if(runif(1) <= p){
      print("accepted badness")
      fitness <- new.fitness
      fitness2 <- new.fitness2
      e <- e2
      de <- de2
    }else{
      replaced <- FALSE
    }
  }
  t <- t - 1
  log <- rbind(log,data.frame(t=t,replaced=replaced,p=p,fitness=fitness,fitness2=fitness2,deltaf=deltaf))
}

if(parallelize){
  stopCluster(c1)
}

wpe.gain <- first.fitness-fitness
akv.gain <- first.fitness2-fitness2

etime <- as.numeric(Sys.time())

print(cat("FITNESS",first.fitness,fitness,first.fitness2,fitness2,etime,""))
print(cat("SAMPLE",e,""))

save(n,first.fitness,first.fitness2,dcoords,wpe.gain,akv.gain,fitness,fitness2,etime,e_before,
     e,log,tmax,candidates,
     file=paste(sep="","sa_slave_",runid,"_",etime,".RData"))

print(cat("DONE",""))
