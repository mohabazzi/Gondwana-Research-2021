#################################################
### Read TPS file with different point number ###
#################################################

# This function takes a file from tpsDig in which no landmarks are assigned and a different number of points have been generated
# for each specimen. The function will take those points and resample the curve(s) to a set number of points (landmarks)
# that can then be used in 'geomorph'. The 'divide.curve' option can be used to use the number of curves assign in tpsDig or if
# FALSE will combine all points into a single curve prior to resampling. The raw points, resampled, and resample & scaled are
# returned.

read2Dtps.noLMs.old <- function(file,ncurve,divide.curve = TRUE,curve.1.pts = NULL,curve.2.pts = NULL,plot = FALSE,tooth = NULL,makeSliders = TRUE,scale = FALSE) {
  require(sp)
  ignore.case = TRUE
  tpsfile <- scan(file, what = "char", sep = "\n", quiet = TRUE)
  scale.lines <- grep("SCALE=", tpsfile, ignore.case)
  lm.lines <- grep("LM=", tpsfile, ignore.case)
  n <- nspecs <- length(lm.lines)
  k <- 2
  ncurvepts <- as.numeric(sub("POINTS=", "", tpsfile[grep("POINTS=", tpsfile, ignore.case)], ignore.case))
  final.results <- vector('list',5);names(final.results) <- c('rawLMs','points.curve','spacedLMs','coords','sliders')
  labels <- tpsfile[grep("IMAGE=", tpsfile, ignore.case)]
  raw.LM <- vector('list',n); names(raw.LM) <- sub("IMAGE=","",labels)
  spaced.LM <- vector('list',n); names(spaced.LM) <- sub("IMAGE=","",labels)
  num.points <- matrix(,nrow = n, ncol = k) 
  for(i in 1:n) {
    #parse out specimens
    sp <- tpsfile[lm.lines[i]:scale.lines[i]]
    curvepts <- grep("POINTS=", sp, ignore.case)
    ncurvepts <- as.numeric(sub("POINTS=", "", sp[curvepts], ignore.case))
    num.points[i,] <- ncurvepts
    raw.curvepts <- vector('list',2)
    #apex <- ncurvepts[1]
    if(divide.curve == TRUE) {
      spaced.curvepts <- vector('list',2)
      for(j in 1:ncurve) {
        #parse out curves
        curve <- sp[(curvepts[j]+1):(curvepts[j]+ncurvepts[j])]
        curve <- matrix(as.numeric(unlist(strsplit(curve, "\\s+"))), ncol = k, byrow = T)
        raw.curvepts[[j]] <- curve
        #equally space & set number of points for each curve separately
        L.lands <- Line(curve)
        if(j == 1) semi.lands <- spsample(L.lands,curve.1.pts,type="regular")@coords
        if(j == 2) semi.lands <- spsample(L.lands,curve.2.pts,type="regular")@coords
        spaced.curvepts[[j]] <- semi.lands
      }
      spaced.LM[[i]] <- rbind(spaced.curvepts[[1]],spaced.curvepts[[2]])
      raw.LM[[i]] <- rbind(raw.curvepts[[1]],raw.curvepts[[2]][-1,])
      p <- nrow(rbind(spaced.curvepts[[1]],spaced.curvepts[[2]]))
    }
    if(divide.curve == FALSE) {
      for(j in 1:ncurve) {
        #parse out curves
        curve <- sp[(curvepts[j]+1):(curvepts[j]+ncurvepts[j])]
        curve <- matrix(as.numeric(unlist(strsplit(curve, "\\s+"))), ncol = k, byrow = T)
        raw.curvepts[[j]] <- curve
      }
      L.lands <- Line(rbind(raw.curvepts[[1]],raw.curvepts[[2]][-1,]))
      semi.lands <- spsample(L.lands,curve.1.pts,type="regular")@coords
      spaced.LM[[i]] <- semi.lands
      raw.LM[[i]] <- rbind(raw.curvepts[[1]],raw.curvepts[[2]][-1,])
      p <- curve.1.pts
    }
  }
  coords <- array(unlist(spaced.LM),c(p, k, n),dimnames = list(NULL,NULL,names(spaced.LM)))
  if(scale) {
    imscale <- as.numeric(sub("SCALE=", "", tpsfile[grep("SCALE",tpsfile, ignore.case)], ignore.case))
    imscale <- aperm(array(rep(imscale, p * k), c(n, k, p)),c(3, 2, 1))
    coords <- coords * imscale
  }
  final.results$rawLMs <- raw.LM
  final.results$points.curve <- num.points
  final.results$spacedLMs <- spaced.LM
  final.results$coords <- coords
  if(plot == TRUE) {
    par(mfrow = c(1,2))
    plot(raw.LM[[tooth]],xlab = NA,ylab = NA,main = "Digitized Points")
    plot(coords[,,tooth],xlab = NA,ylab = NA,main = "Resampled & Scaled Points")
  }
  if(makeSliders == TRUE) final.results$sliders <- data.frame(before = 1:(p-2),slider = 2:(p-1),after = 3:p)
  return(final.results)
}

read2Dtps.noLMs <- function(file,ncurve,divide.curve = TRUE,curve.1.pts = NULL,curve.2.pts = NULL,plot = FALSE,tooth = NULL,makeSliders = TRUE,scale = FALSE) {
  require(sp)
  ignore.case = TRUE
  tpsfile <- scan(file, what = "char", sep = "\n", quiet = TRUE)
  scale.lines <- grep("SCALE=", tpsfile, ignore.case)
  lm.lines <- grep("LM=", tpsfile, ignore.case)
  n <- nspecs <- length(lm.lines)
  k <- 2
  ncurvepts <- as.numeric(sub("POINTS=", "", tpsfile[grep("POINTS=", tpsfile, ignore.case)], ignore.case))
  final.results <- vector('list',5);names(final.results) <- c('rawLMs','points.curve','spacedLMs','coords','sliders')
  labels <- tpsfile[grep("IMAGE=", tpsfile, ignore.case)]
  raw.LM <- vector('list',n); names(raw.LM) <- sub("IMAGE=","",labels)
  spaced.LM <- vector('list',n); names(spaced.LM) <- sub("IMAGE=","",labels)
  num.points <- matrix(,nrow = n, ncol = k) 
  for(i in 1:n) {
    #parse out specimens
    if(i < n) sp <- tpsfile[lm.lines[i]:lm.lines[i+1]-1]
    else sp <- tpsfile[tail(lm.lines,n = 1):length(tpsfile)]
    curvepts <- grep("POINTS=", sp, ignore.case)
    ncurvepts <- as.numeric(sub("POINTS=", "", sp[curvepts], ignore.case))
    num.points[i,] <- ncurvepts
    raw.curvepts <- vector('list',2)
    #apex <- ncurvepts[1]
    if(divide.curve == TRUE) {
      spaced.curvepts <- vector('list',2)
      for(j in 1:ncurve) {
        #parse out curves
        curve <- sp[(curvepts[j]+1):(curvepts[j]+ncurvepts[j])]
        curve <- matrix(as.numeric(unlist(strsplit(curve, "\\s+"))), ncol = k, byrow = T)
        raw.curvepts[[j]] <- curve
        #equally space & set number of points for each curve separately
        if(j == 1) {
          L.lands <- Line(curve)
          semi.lands <- spsample(L.lands,curve.1.pts,type="regular",offset = c(0,1))@coords
          spaced.curvepts[[j]] <- semi.lands
        }
        if(j == 2) {
          curve <- curve[ncurvepts[j]:1,]
          L.lands <- Line(curve)
          semi.lands <- spsample(L.lands,curve.2.pts-1,type="regular",offset = c(0,1))@coords
          spaced.curvepts[[j]] <- rbind(semi.lands,curve[ncurvepts[j],])
        }
      }
      spaced.LM[[i]] <- rbind(spaced.curvepts[[1]],spaced.curvepts[[2]])
      raw.LM[[i]] <- rbind(raw.curvepts[[1]],raw.curvepts[[2]][-1,])
      p <- nrow(rbind(spaced.curvepts[[1]],spaced.curvepts[[2]]))
    }
    if(divide.curve == FALSE) {
      for(j in 1:ncurve) {
        #parse out curves
        curve <- sp[(curvepts[j]+1):(curvepts[j]+ncurvepts[j])]
        curve <- matrix(as.numeric(unlist(strsplit(curve, "\\s+"))), ncol = k, byrow = T)
        raw.curvepts[[j]] <- curve
      }
      L.lands <- Line(rbind(raw.curvepts[[1]],raw.curvepts[[2]][-1,]))
      semi.lands <- spsample(L.lands,curve.1.pts,type="regular")@coords
      spaced.LM[[i]] <- semi.lands
      raw.LM[[i]] <- rbind(raw.curvepts[[1]],raw.curvepts[[2]][-1,])
      p <- curve.1.pts
    }
  }
  coords <- array(unlist(spaced.LM),c(p, k, n),dimnames = list(NULL,NULL,names(spaced.LM)))
  if(scale) {
    imscale <- as.numeric(sub("SCALE=", "", tpsfile[grep("SCALE",tpsfile, ignore.case)], ignore.case))
    imscale <- aperm(array(rep(imscale, p * k), c(n, k, p)),c(3, 2, 1))
    coords <- coords * imscale
  }
  final.results$rawLMs <- raw.LM
  final.results$points.curve <- num.points
  final.results$spacedLMs <- spaced.LM
  final.results$coords <- coords
  if(plot == TRUE) {
    par(mfrow = c(1,2))
    plot(raw.LM[[tooth]],xlab = NA,ylab = NA,main = "Digitized Points")
    plot(coords[,,tooth],xlab = NA,ylab = NA,main = "Resampled & Scaled Points")
  }
  if(makeSliders == TRUE) final.results$sliders <- data.frame(before = 1:(p-2),slider = 2:(p-1),after = 3:p)
  return(final.results)
}

## function to check teeth

check.lands <- function(coords,curves = NULL) {
  if(.Platform$OS.type == "unix") quartz(width = 18.5,height = 9.5)
  else windows(width = 18.5,height = 9.5)
  par(mfrow = c(5,10))
  if(is.null(curves)) for(i in 1:dim(coords)[3]) plot(coords[,,i],main = i,xlab = NA,ylab = NA, pch = 19)
  else {
    if(class(coords) == "array") {
      for(i in 1:dim(coords)[3]) {
        plot(coords[,,i],main = i,xlab = NA,ylab = NA,type = "n")
        for(j in 1:length(curves)) {
          points(coords[curves[[j]],,i],pch = 19,col = rainbow(length(curves)))
        }
      }
    }
    if(class(coords) == "list") {
      for(i in 1:length(coords)) {
        plot(coords[[i]],main = i,xlab = NA,ylab = NA,type = "n")
        for(j in 1:length(curves)) {
          points(coords[[i]][curves[[j]],],pch = 19,col = rainbow(length(curves)))
        }
      }
    }
  }
}