colcat <- function(a, b, pos, x){
  names <- colnames(x)
  if (missing(b)) {
    idx <- grep(a, sapply(names, function(x) substr(x, pos,pos)))
    new_mtx <- x[,idx]
  } else {
    idxa <- grep(a, sapply(names, function(x) substr(x, pos,pos)))
    mat_a <- x[,idxa]
    idxb <- grep(b, sapply(names, function(x) substr(x, pos,pos)))
    mat_b <- x[,idxb]
    new_mtx <- cbind(mat_a,mat_b)
  }
  return(new_mtx)
}



empty_file <- function(mat, ref, equal = FALSE){
  if (equal==FALSE) {
    filt = apply(mat, 1, function(x) ifelse(all(x < ref), 0, 1))
    filtered = mat[-which(filt == 0),]
  } else {
    filt = apply(mat, 1, function(x) ifelse(all(x == ref), 0, 1))
    filtered = mat[-which(filt == 0),]
  }
  return(filtered)
}

combine.mtx <- function(...){
  temp <- list(...)
  for (i in 1:(length(temp)-1)) {
    if (i == 1) {
      new_mtx = merge(temp[[i]], temp[[i+1]], by = 'row.names', all =T)
      rownames(new_mtx) <- new_mtx$Row.names
      new_mtx$Row.names <- NULL
    } else {
      new_mtx = merge(new_mtx, temp[[i+1]], by = 'row.names', all =T)
      rownames(new_mtx) <- new_mtx$Row.names
      new_mtx$Row.names <- NULL
    }
    
  }
  new_mtx[is.na(new_mtx)] <- 0
  return(new_mtx)
}

plot.densities <- function(m, m2, m3, m4, lty=rep(1:6, each=length(palette())), lwd=1, cex=1, legends=TRUE, legend.pos="topleft", legend.cols=1, legend.cex=1, horiz=FALSE, xlim=if(horiz) c(mny,mxy) else c(mnx,mxx),ylim=if(!horiz) c(mny,mxy) else c(mnx,mxx), pch=NULL, type="l", npch=20, col=1:length(d), normalize=FALSE, q=NULL, overall=length(d) > 4, ylab="Density", xlab="", main="", ...) { 
  if (!missing(m2) || !missing(m3) || !missing(m4)) {
    xn <- deparse(substitute(m))
    m <- as.vector(m)
    if (!missing(m2)) { m <- cbind(m,as.vector(m2)); xn <- c(xn,deparse(substitute(m2))) }
    if (!missing(m3)) { m <- cbind(m,as.vector(m3)); xn <- c(xn,deparse(substitute(m3))) }
    if (!missing(m4)) { m <- cbind(m,as.vector(m4)); xn <- c(xn,deparse(substitute(m4))) }
    colnames(m) <- xn
  }
  else if (!is.list(m)){
    m <- as.matrix(m)
  }
  if (is.list(m)) d <- lapply(m, function(x) { suppressWarnings(density(x[is.finite(x)],...)) })
  else d <- apply(m, 2, function(x) { suppressWarnings(density(x[is.finite(x)],...)) })
  lty <- circularize(lty, length(d))
  lwd <- circularize(lwd, length(d))
  cex <- circularize(cex, length(d))
  if (!is.null(q)) {
    if (is.list(m)) qq <- lapply(m, function(x) { suppressWarnings(quantile(x[is.finite(x)], q)) })
    else qq <- apply(m, 2, function(x) { suppressWarnings(quantile(x[is.finite(x)], q)) })
  }
  if (overall) {
    o <- unlist(m)
    n <- names(d)
    d[[length(d)+1]] <- density(o[is.finite(o)], ...)
    names(d) <- c(n, "OVERALL")
    lty[length(d)] <- 1
    lwd[length(d)] <- 3
    col[length(d)] <- 1
    cex[length(d)] <- 1
    if (!is.null(q)) qq[length(d)] <- suppressWarnings(quantile(o[is.finite(o)], q))
  }
  if (normalize) {
    for (i in 1:length(d)) {
      d[[i]]$y <- d[[i]]$y/max(d[[i]]$y)
    }
  }
  mxy <- max(unlist(lapply(d, function(x) max(x$y))))
  mny <- min(unlist(lapply(d, function(x) min(x$y))))
  mxx <- max(unlist(lapply(d, function(x) max(x$x))))
  mnx <- min(unlist(lapply(d, function(x) min(x$x))))
  lty <- rep(lty, length(d), length.out=length(d))
  pch <- ifelse(type=="l",rep(0,length(d)),rep(pch, length(d), length.out=length(d)))
  suppressWarnings(plot(0,0,type="n",ylim=ylim,xlim=xlim, main=main, xlab=xlab, ylab=ylab, ...))
  cols <- if (horiz) c("y","x") else c("x","y")
  for (i in length(d):1) {
    suppressWarnings(lines(d[[i]][cols], col=col[i], lty=lty[i], pch=pch[i], type=type, lwd=lwd[i], cex=cex[i]))
    if (!is.null(q)) {
      for (k in 1:length(q)) {
        if (horiz) abline(h=qq[[i]][k], lty=k, col=col[i])
        else abline(v=qq[[i]][k], lty=k, col=col[i])
      }
    }
  }
  if (type != "l") {
    pch <- rep(pch, length(d), length.out=length(d))
    if (horiz) for (i in 1:length(d)) suppressWarnings(points(d[[i]]$y[seq(1,length(d[[i]]$y), length.out=npch)], d[[i]]$x[seq(1,length(d[[i]]$y), length.out=npch)], col=col[i], pch=pch[i], ...))
    else   for (i in 1:length(d)) suppressWarnings(points(d[[i]]$x[[seq(1,length(d[[i]]$y), length.out=npch)]], d[[i]]$y[seq(1,length(d[[i]]$y), length.out=npch)], col=col[i], pch=pch[i], ...))
  }
  if (length(legends) > 1 || (!is.character(legends) && legends)) {
    if ((length(legends) == length(d)) && (length(legends) > 0)) xl <- legends
    else {
      xl <- if (is.list(m)) names(m) else colnames(m)
      if (overall && !is.null(xl)) xl <- c(xl, "OVERALL")
    }
    if (!is.null(xl)) {
      plot.legend(legend.pos, xl, col=col[1:length(d)], lty=lty, ncol=legend.cols, cex=legend.cex, pch=if (is.null(pch[1])) NULL else pch)
    }
  }
  invisible(d)
}

quantile.normalization <- function(l.mx, verbose=FALSE, averages=FALSE, center=FALSE, q.func=my.quantile, m.func=mean, qmin=min(avg), qmax=max(avg), qmid=1, qshrink=NULL) { 
  #q.func=quantile # if you want to use the embed slower R quantile (slightly different results)
  #convert to list
  retfn <- function(x) x
  if (is.data.frame(l.mx) || is.matrix(l.mx)) { 
    .rn <- rownames(l.mx)
    l.mx <- lapply(as.list(data.frame(l.mx)), function(x) { names(x) <- .rn; x; })
    if (is.data.frame(l.mx))
      retfn <- function(x) { x <- data.frame(x); rownames(x) <- .rn; x; }
    else
      retfn <- function(x) { x <- data.matrix(data.frame(x)); rownames(x) <- .rn; x; }
  }
  
  #scale values to min~max range, or to part of the median/center quantile (qmid)
  q.scale <- function(v, qmin=min(v), qmax=max(v), qmid=1) {
    qs <- quantile(v, p=c((1-qmid)/2,(qmid+1)/2))
    vmin <- qs[1]
    vmax <- qs[2]
    #3-rule eq: (v - vmin) / (vmax - vmin) = (q - qmin) / (qmax - qmin)
    ((v - vmin) / (vmax - vmin)) * (qmax - qmin) + qmin
  }
  
  #median centering
  if (center) {
    if (verbose) cat("Centering...\n")
    # center to 
    md <- unlist(lapply(l.mx, median, na.rm = TRUE))
    mmd <- median(md)
    l.mx <- lapply(l.mx, function(l) l-(median(l,na.rm=TRUE) - mmd))
  }
  
  #shrinking distributions to median/center quantile range
  if (!is.null(qshrink)) {
    if (verbose) cat("Shrinking...\n")
    minq <- median(unlist(lapply(l.mx, quantile, p=(1-qshrink)/2, na.rm = TRUE)))
    maxq <- median(unlist(lapply(l.mx, quantile, p=(qshrink+1)/2, na.rm = TRUE)))
    l.mx <- lapply(l.mx, q.scale, qmin=minq, qmax=maxq, qmid=qshrink)
  }
  
  #maximum number of values
  m <- 0 
  for (i in 1:length(l.mx)) m <- max(m,sum(is.finite(l.mx[[i]])))
  if (verbose) cat("Max-Quantiles:",m,"\n")
  
  #quantiling
  q.mx <- matrix(NA,nrow=m, ncol=length(l.mx))
  for (i in 1:length(l.mx)) { q.mx[,i] <- q.func(l.mx[[i]],probs=(1:m)/m,na.rm=TRUE); if (verbose) cat("Quantiled:",i,"\n"); }
  
  #averaging
  avg <- numeric(length=m)
  if (verbose) cat("Averaging...\n")
  for (i in 1:m) avg[i] <- m.func(q.mx[i,], na.rm=TRUE) 
  
  #scale averages?
  if (min(avg) != qmin  ||  max(avg) != qmax  ||  qmid < 1) avg <- q.scale(avg, qmin=qmin, qmax=qmax, qmid=qmid)
  
  #normalizing
  q.mx <- list()
  if (verbose) cat("Normalizing...\n")
  for (i in 1:length(l.mx)) { 
    if (verbose) cat("Computing ecdf",i,"...\n");
    avg.cdf <- ecdf(l.mx[[i]][is.finite(l.mx[[i]])])
    r <- pmin(pmax(1,avg.cdf(l.mx[[i]])*m),m)
    xq <- avg[r]
    names(xq) <- names(l.mx[[i]])
    q.mx[[i]] <- xq
  }
  names(q.mx) <- names(l.mx)
  
  #return
  if (verbose) cat("Done!\n");
  if (averages) avg else retfn(q.mx)
}

circularize <- function(x, len) {
  if (length(x) < len) {
    x <- rep(x, length.out=len)
  }
  x
}
