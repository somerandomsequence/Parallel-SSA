source('new_points_common.R')
load('preprocess.RData')

fname <- commandArgs(TRUE)[1]
load(fname)

parallelize <- TRUE
num.children <- 8

if(!exists("fitness") || !exists("fitness2")){
  de <- candidates[e,]
  if(parallelize){
    c1 <- makeForkCluster(num.children)
    kv <- krige.var.par(rbind(dcoords,de[,c("x","y")]),loci,kc,c1)
  }else{
    kv <- krige.var(rbind(dcoords,de[,c("x","y")]),loci,kc)
  }
  vmap <- flipud(matrix(kv,nrow=height,ncol=width,byrow=TRUE))
  fitness <- wpe(rmap,vmap)
  fitness2 <- mean(sqrt(vmap))
}

wpe.gain <- first.fitness-fitness
akv.gain <- first.fitness2-fitness2

print(cat("FITNESS",first.fitness,fitness,first.fitness2,fitness2,""))
print(cat("SAMPLE",e,""))
