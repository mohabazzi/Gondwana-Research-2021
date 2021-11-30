# Bootstrapping and rarefaction for landmark data -------------------------

# x - output from gpagen$coords
# groups - factor specifying the group to which each specimen in x belongs
# order - numeric vector specifying the order in which levels(groups) should be plotted
# replicates - single value specifying the number of replicates for the bootstrap
# rarefy - list of the parameters to use during rarefaction (e.g., list(ist(min.N = 5,reps = 999)))

error.plot <- function(gpa.coords, groups, order, blank, replicates = 999, rarefy.par) {
  require(plotrix)
  require(broom)
  options(warn = -1) # Ignore all warnings.
  
  # Bootstrapping
  boot.sampling <- function(data, n) {
    sub.gpa <- data$coords[,,base::sample(1:n,replace = TRUE)]
    sub.data <- geomorph.data.frame(coords = sub.gpa)
    invisible(capture.output(disp <- morphol.disparity(coords~1,data = sub.data,print.progress = FALSE)))
    return(disp)
  }
  grps <- levels(droplevels(groups))
  if(blank == TRUE) grps <- levels(droplevels(groups))[-1]
  res.mat <- matrix(nrow = length(grps),ncol = 3); rownames(res.mat) <- grps; colnames(res.mat) <- c("dispariy","LowerBoot.PI","UpperBoot.PI","LowerBoot.CI","UpperBoot.CI")
  for(i in 1:length(grps)) {
    grp <- which(groups == grps[i])
    # print(gpa.coords[,,grp])
    grp.data <- geomorph.data.frame(coords = gpa.coords[,,grp])
    invisible(capture.output(grp.disp <- morphol.disparity(coords~1,data = grp.data)))
    invisible(capture.output(boot <- base::replicate(n = replicates,expr = boot.sampling(grp.data,dim(grp.data$coords)[3]))))
    if(exists("intervals") != TRUE) stop("Source Intervals.R function first")
    grp.int <- intervals(boot)
    res.mat[i,] <- c(grp.disp,grp.int$conf.low,grp.int$conf.high)
  }
  
  # Rarefaction
  rarefy.sampling <- function(data, n, size) {
    sub.gpa <- data$coords[,,sample(1:n,size = size,replace = FALSE)]
    sub.data <- geomorph.data.frame(coords = sub.gpa)
    invisible(capture.output(rare <- morphol.disparity(coords~1,data = sub.data,print.progress = FALSE)))
    return(rare)
  }
  if(is.null(rarefy.par$min.N)) stop("'min.N' is needed in order to rarefy")
  if(is.null(rarefy.par$reps)) stop("The number of replicates to be run in the rarefaction needs to be specified")
  rare.res <- vector(mode = "list",length = length(grps))
  for(j in order) {
    grp <- which(groups == grps[j])
    Ns <- rarefy.par$min.N:length(grp)
    rare.disp <- matrix(nrow = length(Ns),ncol = 5); colnames(rare.disp) <- c('rare.disp','lower.PI','upper.PI','lower.CI','upper.CI'); rownames(rare.disp) <- as.character(Ns)
    for(i in 1:length(Ns)) {
      grp.data <- geomorph.data.frame(coords = gpa.coords[,,grp])
      invisible(capture.output(raref <- replicate(n = rarefy.par$reps,expr = rarefy.sampling(grp.data,dim(grp.data$coords)[3],size = Ns[i]))))
      stats <- intervals(raref)
      rare.disp[i,] <- c(stats$mean,stats$Lower.PI,stats$Upper.PI,stats$Lower.CI,stats$Upper.CI)
    }
    rare.res[[j]] <- rare.disp; names(rare.res) <- levels(droplevels(groups))
  }

  options(warn = 0)

  return(list(bootstrap.results = res.mat[order,],rarefaction.results = rare.res[order]))
}
