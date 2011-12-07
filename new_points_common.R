library(geoR)
library(lattice)
library(matlab) # for flipud()
library(splancs) # for pip (point in polygon)
library(parallel) # for krige.var's parallel version

# This is a definition of kriging variance that is based directly 
# on geoR's krige.conv with the nonapplicable stuff stripped out
krige.var.geoR <- function (coords, data, locations, krige){

    if (length(class(krige)) == 0 || class(krige) != "krige.geoR") {
      if (!is.list(krige))
        stop("krige.conv: the argument krige only takes a list or an output of the function krige.control")
      else {
        krige.names <- c("type.krige", "trend.d", "trend.l",
              "obj.model", "beta", "cov.model", "cov.pars",
              "kappa", "nugget", "micro.scale", "dist.epsilon",
              "lambda", "aniso.pars")
        krige.user <- krige
        krige <- list()
        if (length(krige.user) > 0) {
          for (i in 1:length(krige.user)) {
            n.match <- match.arg(names(krige.user)[i],
              krige.names)
            krige[[n.match]] <- krige.user[[i]]
          }
        }
        if (is.null(krige$type.krige))
          krige$type.krige <- "ok"
          if (is.null(krige$trend.d))
            krige$trend.d <- "cte"
          if (is.null(krige$trend.l))
            krige$trend.l <- "cte"
          if (is.null(krige$obj.model))
            krige$obj.model <- NULL
          if (is.null(krige$beta))
            krige$beta <- NULL
          if (is.null(krige$cov.model))
            krige$cov.model <- "matern"
          if (is.null(krige$cov.pars))
            stop("covariance parameters (sigmasq and phi) should be provided in cov.pars")
          if (is.null(krige$kappa))
            krige$kappa <- 0.5
          if (is.null(krige$nugget))
            krige$nugget <- 0
          if (is.null(krige$micro.scale))
            krige$micro.scale <- 0
          if (is.null(krige$dist.epsilon))
            krige$dist.epsilon <- 1e-10
          if (is.null(krige$aniso.pars))
            krige$aniso.pars <- NULL
          if (is.null(krige$lambda))
            krige$lambda <- 1
            krige <- krige.control(type.krige = krige$type.krige,
            trend.d = krige$trend.d, trend.l = krige$trend.l,
              obj.model = krige$obj.model, beta = krige$beta,
              cov.model = krige$cov.model, cov.pars = krige$cov.pars,
              kappa = krige$kappa, nugget = krige$nugget,
              micro.scale = krige$micro.scale, dist.epsilon = krige$dist.epsilon,
              aniso.pars = krige$aniso.pars, lambda = krige$lambda)
        }
    }

    signal <- FALSE
    cov.model <- krige$cov.model
    kappa <- krige$kappa
    lambda <- krige$lambda
    beta <- krige$beta
    cov.pars <- krige$cov.pars
    nugget <- krige$nugget
    micro.scale <- krige$micro.scale
    aniso.pars <- krige$aniso.pars

    coords <- as.matrix(coords)
    dimnames(coords) <- list(NULL, NULL)

    if (is.vector(cov.pars)) {
        sigmasq <- cov.pars[1]
        phi <- cov.pars[2]
        cpars <- c(1, phi)
    }   
    sill.partial <- sum(sigmasq)
    if (sill.partial < 1e-16) {
        tausq.rel <- 0
        tausq.rel.micro <- 0
    }   
    else {
        tausq.rel <- nugget/sum(sigmasq)
        tausq.rel.micro <- micro.scale/sum(sigmasq)
    }

    Vcov <- varcov.spatial(coords = coords, cov.model = cov.model,
        kappa = kappa, nugget = tausq.rel, cov.pars = cpars)$varcov
    trend.d <- matrix(1,nrow=nrow(dcoords))
    trend.l <- unclass(trend.spatial(trend = krige$trend.l, geodata = list(coords = locations)))
    ivtt <- solve(Vcov, trend.d)
    ttivtt <- crossprod(ivtt, trend.d)

    v0 <- loccoords(coords = coords, locations = locations)

    nug.factor <- ifelse(signal, tausq.rel.micro, tausq.rel)
    ind.v0 <- which(v0 < krige$dist.epsilon)
    v0 <- cov.spatial(obj = v0, cov.model = cov.model, kappa = kappa, cov.pars = cpars)
    v0[ind.v0] <- 1 + nug.factor
    ivv0 <- solve(Vcov, v0)
    b <- crossprod(cbind(data, trend.d), ivv0)
    b <- t(trend.l) - b[-1, , drop = FALSE]
    tv0ivv0 <- colSums(v0 * ivv0)
    bitb <- colSums(b * solve(ttivtt, b))
    krige.var <- sill.partial * drop(1 + nug.factor - tv0ivv0 + bitb)
    krige.var[kc$krige.var < 1e-08] <- 0
    return(krige.var)
}

# This definition of kriging variance taken from  equation (1) of Delmelle et al. (2009), and
# is compatable with equation 4.18 in Spatial Statistics by Ripley. It appears to be a constant
# shift off the kriging variance computed by the geoR krige.conv method
krige.var <- function(dcoords,loci,kc){
  cs <- cov.spatial(obj=loccoords(coords=dcoords,locations=loci),cov.model=kc$cov.model,cov.pars=kc$cov.pars,kappa=kc$cov.kappa)
  vcinv <- varcov.spatial(coords=dcoords, cov.model=kc$cov.model,cov.pars=kc$cov.pars,kappa=kc$kappa,nugget=kc$nugget,inv=TRUE)$inverse
  sigmasq <- kc$cov.pars[1]
  sigmak <- lapply(seq(1,nrow(loci)),function(i){ sigmasq - t(cs[,i]) %*% vcinv %*% t(t(cs[,i])) })
  rm(vcinv,cs)
  return(sigmak)
}  

# fourth argument is a cluster made with a command like makeForkCluster(N)
krige.var.par <- function(dcoords,loci,kc,c1){
  cs <- cov.spatial.par(obj=loccoords(coords=dcoords,locations=loci),cov.model=kc$cov.model,cov.pars=kc$cov.pars,kappa=kc$cov.kappa)
  vcinv <- varcov.spatial(coords=dcoords, cov.model=kc$cov.model,cov.pars=kc$cov.pars,kappa=kc$kappa,nugget=kc$nugget,inv=TRUE)$inverse
  sigmasq <- kc$cov.pars[1]
  #c1 <- makeForkCluster(n.par) # doesn't work on windoze, but faster on *nix
  sigmak <- parLapply(c1,seq(1,nrow(loci)),function(i){ sigmasq - t(cs[,i]) %*% vcinv %*% t(t(cs[,i])) })
  #stopCluster(c1)
  rm(vcinv,cs)
  return(sigmak)
} 

roughness <- function(map,height,width,nr=1,pix.per.m=0.2,beta=1.5,alpha=1.0){
  height <- nrow(map)
  width <- ncol(map)
  ret <- map
  neigh <- expand.grid(seq(-nr,nr),seq(-nr,nr))
  dsum <- 0.0
  
  for(k in nrow(neigh)){
    x <- neigh[k,1]
    y <- neigh[k,2]
    if(x == 0 && y == 0) next;
    d <- sqrt(x^2 + y^2)/pix.per.m
    dsum <- dsum + d
  }
  for(i in seq(1,height)){
    for(j in seq(1,width)){
      s <- 0.0
      v <- map[i,j]
      for(k in nrow(neigh)){ 
        x <- neigh[k,1]
        y <- neigh[k,2]
        xi <- x + j
        yi <- y + i
        if(x == 0 && y == 0) next;
        if(xi < 1 || yi < 1 || xi > width || yi > height) next;
        d <- sqrt(x^2 + y^2)/pix.per.m
        v2 <- map[yi,xi]
        s <- s + (d^(-beta) * (v2 - v)^2)/dsum
      }
      ret[i,j] <- s
    }
  }
  ret <- (ret/max(ret))^alpha
  return(ret)
}

wpe <- function(rmap,vmap){
  # note this is not a matrix mult (%*%) so will just multiple rmap[i,i]*vmap[i,i]
  mean(rmap*vmap)
}

rand.x <- function(dcoords){ as.integer(min(dcoords$x) + runif(1)*(max(dcoords$x)-min(dcoords$x))) }
rand.y <- function(dcoords){ as.integer(min(dcoords$y) + runif(1)*(max(dcoords$y)-min(dcoords$y))) }

distance.mask <- function(dmask.thresh=100,width,height,pix.per.meter){
  dmask.thresh <- 100
  dmask <- matrix(FALSE,nrow=height,ncol=width)
  for(i in seq(1,height)){
    for(j in seq(1,width)){
      x0 <- eastrng[1] + j/pix.per.meter
      y0 <- northrng[1] + i/pix.per.meter
      t <- sqrt((dcoords$x-x0)^2 + (dcoords$y-y0)^2)
      if(any(t <= dmask.thresh)) dmask[i,j] <- TRUE
    }
  }
  return(dmask)
}

inhull <- function(hull,x,y){ nrow(pip(data.frame(x=x,y=y),hull)) > 0 }
