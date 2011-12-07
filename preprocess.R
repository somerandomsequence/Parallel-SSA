source('new_points_common.R')

load('emp_variogram_figs/44:51:DB:00:07:01/44:51:DB:00:07:01_100wimax_0.2_wneg_detrend_trunc_krige.RData')

# filter out erroneous location
d <- d[dcoords$x <= 2018500,]
dcoords <- dcoords[dcoords$x <= 2018500,]

pix.per.meter <- 0.2

# resize sample grid if necessary
#new.pix.per.meter <- 0.02 # 1 per 50m
new.pix.per.meter <- 0.2
if(pix.per.meter != new.pix.per.meter){
  width <- round(diff(range(dcoords$x))*new.pix.per.meter)
  height <- round(diff(range(dcoords$y))*new.pix.per.meter)
  print(paste("Size:",width,"x",height))
  loci <- expand.grid(seq(min(dcoords$x),max(dcoords$x),length.out=width),seq(min(dcoords$y),max(dcoords$y),length.out=height))
  pix.per.meter <- new.pix.per.meter
}

# distance mask
dchull <- dcoords[chull(dcoords),]

# do initial interpolation
k <- krige.conv(coords=dcoords,data=d$sig2,locations=loci,krige=kc)
map <- flipud(matrix(k$predict,nrow=height,ncol=width,byrow=TRUE))
rmap <- roughness(map,height,width)
rm(map,k)

# precompute candidate locations
candidates <- NULL
for(yi in seq(1,height)){
  for(xi in seq(1,width)){
    if(xi < 1 || yi < 1 || xi > width || yi > width) next;
    x <- eastrng[1] + xi/pix.per.meter
    y <- northrng[1] + yi/pix.per.meter
    if(!inhull(dchull,x,y)) next
    #if(!dmask[k,l]) next 
    candidates <- rbind(candidates,data.frame(xi=xi,yi=yi,x=x,y=y,used=FALSE))
    #print(paste(xi,yi,x,y))
  }
}

dcoords$x <- as.integer(dcoords$x)
dcoords$y <- as.integer(dcoords$y)

save(candidates,kc,loci,dcoords,rmap,width,height,file="preprocess.RData")
write.table(candidates,col.names=FALSE,quote=FALSE,file="candidates.txt")
