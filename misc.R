#' plotROCs
#'
#' Plot ROC curves for given scores.
#'
#' @param scores A data.frame with the different types of scores as columns.
#' @param truth A vector of the true class corresponding to each row of `scores`
#' @param fdr logical; whether to plot FDR instead of FPR on the x axis
#' @param th an optional vector of signifincance thresholds at which to plot
#' points
#' @param printAUC logical; whether to print the AUC
#' @param showLegend logical; whether to print the legend
#' @param prop.wrong.neg Proportion of the negative truths that is wrong
#' @param prop.wrong.pos Proportion of the positive truths that is wrong
#' @param ... passed to `geom_point`
#'
#' @return a ggplot
#' @import ggplot2
#'
#' @examples
#' myscores <- list( test=1:10 )
#' truth <- sample(c(TRUE,FALSE), 10, TRUE)
#' plotROCs(  myscores, truth )
#'
#' @export
plotROCs <- function(scores, truth, fdr=FALSE, th=c(), showLegend=TRUE,
                     prop.wrong.neg=0, prop.wrong.pos=0,
                     printAUC=is.vector(scores) || length(scores)==1,
                     addNull=FALSE, ...){
  truth <- as.integer(as.factor(truth))-1
  scores <- as.data.frame(scores)
  roclist <- lapply(scores, FUN=function(x){
    labels <- truth[order(x, decreasing=TRUE)]
    afpr <- (cumsum(!labels)-(prop.wrong.neg*length(labels)))/((1-prop.wrong.neg)*sum(!labels))
    afpr[afpr<0] <- 0
    data.frame(FPR=cumsum(!labels)/sum(!labels),
               aFPR=afpr,
               FDR=cumsum(!labels)/seq_along(labels),
               TPR=cumsum(labels)/sum(labels))
  })
  methods <- factor( rep(names(roclist),
                         vapply(roclist, FUN.VALUE=integer(1),
                                FUN=function(x) length(x$TPR))
  ),
  levels=names(roclist) )
  d <- data.frame( method=methods,
                   FPR=unlist(lapply(roclist,FUN=function(x) x$FPR)),
                   adjusted.FPR=unlist(lapply(roclist,FUN=function(x) x$aFPR)),
                   TPR=unlist(lapply(roclist,FUN=function(x) x$TPR)),
                   FDR=unlist(lapply(roclist,FUN=function(x) x$FDR)))
  psub <- function(x,y){ x <- x-y; x[x<0] <- 0; x/(1-y) }
  pretop <- function(x,y){ x[x>y] <- y; x/y }
  if(fdr){
    p <- ggplot(d, aes(FDR, TPR)) + geom_path(size=1.2, aes(colour=method))
    auc <- DescTools::AUC(pretop(d$TPR,1-prop.wrong.pos), 1-psub(d$FDR,prop.wrong.neg),
                          from=0, to=0)
    if(addNull) p <- p + geom_vline(xintercept=1-(sum(truth)/length(truth)),
                                    linetype="dotted", colour="red")
  }else{
    if(prop.wrong.neg>0){
      p <- ggplot(d, aes(adjusted.FPR, TPR)) + geom_line(size=1.2, aes(colour=method))
      x <- d$adjusted.FPR
    }else{
      p <- ggplot(d, aes(FPR, TPR)) + geom_line(size=1.2, aes(colour=method))
      x <- psub(d$FPR,prop.wrong.neg)
    }
    auc <- DescTools::AUC(x, pretop(d$TPR,1-prop.wrong.pos), from=0, to=0)
    if(addNull) p <- p + geom_abline(slope=1, linetype="dashed", colour="red")
  }
  dummy <- data.frame(FDR=1,FPR=1,TPR=1,adjusted.FPR=1)
  p <- p + scale_x_continuous(limits=c(0,1), expand=c(-0.01,0.01)) +
    scale_y_continuous(limits=c(0,1), expand=c(-0.01,0.01))
  if(prop.wrong.pos>0) p <- p + geom_hline(yintercept=1-prop.wrong.pos, linetype="dashed") +
    geom_rect(data=dummy,xmin=ifelse(fdr,prop.wrong.neg,0),xmax=1,ymin=1-prop.wrong.pos,ymax=1,alpha=0.2,fill="grey25") +
    #annotate(geom="text", x=0.99, y=1.02-prop.wrong.pos, label=round(1-prop.wrong.pos,3), hjust="right", vjust="bottom") +
    annotate(geom="text", x=0.02+prop.wrong.neg, y=1-prop.wrong.pos/2, label="Homotypic doublets\n(false positives)", vjust="center", hjust="left", size=3.8)
  if(!showLegend) p <- p + theme(legend.position="none")
  if(prop.wrong.neg>0 & fdr) p <- p + geom_vline(xintercept=prop.wrong.neg, linetype="dashed") +
    geom_rect(data=dummy,xmin=0,xmax=prop.wrong.neg,ymin=0,ymax=1-prop.wrong.pos,alpha=0.2,fill="grey25") +
    annotate(geom="text", x=prop.wrong.neg/2, y=(1-prop.wrong.pos)/2, label="Intra-genotype heterotypic\n(false negatives)", vjust="middle", hjust="center", angle=90, size=3.8)
  if(printAUC) p <- p + annotate("text",x=0.99, y=0.01, vjust="bottom", hjust="right",
               label=paste0(ifelse(prop.wrong.neg>0 | prop.wrong.pos>0,"Adjusted",""),
                            ifelse(fdr," AUPRC","AUROC"),"\n", round(auc,3)))
  if(length(th)>0){
    for(th1 in th){
      w <- lapply(scores, FUN=function(x) max(which(sort(x, decreasing=TRUE)>=th1)))
      d2 <- data.frame( method=names(roclist),
                        FPR=mapply(x=roclist,w=w,FUN=function(x,w) x$FPR[w]),
                        adjusted.FPR=mapply(x=roclist,w=w,FUN=function(x,w) x$aFPR[w]),
                        FDR=mapply(x=roclist,w=w,FUN=function(x,w) x$FDR[w]),
                        TPR=unlist(mapply(x=roclist,w=w,FUN=function(x,w) x$TPR[w])) )
      p <- p + geom_point(data=d2, aes(colour=method), ...)
    }
  }
  p
}


.vec2mat <- function(x, size, defval=NA_integer_, tnames=NULL){
  if(length(size)>1){
    if(is.null(tnames)) tnames <- names(size)
    size <- length(size)
  }
  m <- matrix(defval, nrow=size, ncol=size)
  m[lower.tri(m)] <- x
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  if(!is.null(tnames)) colnames(m) <- row.names(m) <- tnames
  m
}

# simulates doublet counts
simDblCounts <- function( ncells=c(1500,200,1000,500,300,800),
                          stickiness=c(1,1,1,1,1,2), size=1, dblrate=0.1,
                          comb=NULL, retMat=TRUE){
  if(is.null(names(ncells))) names(ncells) <- LETTERS[seq_along(ncells)]
  base.e <- ncells/sum(ncells)
  base.e <- base.e %*% t(base.e)
  base.e <- base.e[lower.tri(base.e)]
  s <- stickiness %*% t(stickiness)
  s <- s[lower.tri(s)]
  nd <- sum(ncells)*dblrate
  e <- base.e * s
  truth.comb <- rep(FALSE,length(e))
  if(!is.null(comb)){
    w <- sample.int(length(e),length(comb))
    e[w] <- e[w]*comb
    truth.comb[w] <- comb!=1
  }
  e <- nd * e/sum(e)
  base.e <- nd * base.e/sum(base.e)

  if(is.null(size) || size==0){
    o <- rpois(length(e), e)
  }else{
    o <- rnbinom(length(e),mu=e,size=size)
  }
  if(retMat){
    o <- .vec2mat(o, size=ncells)
    base.e <- .vec2mat(base.e, size=ncells)
    truth.comb <- .vec2mat(truth.comb, size=ncells, defval=FALSE)
  }
  list( observed=o, expected=base.e, log2enrich=log2((o+1)/(base.e+1)),
        truth.comb=truth.comb )
}


#' runParams
#'
#' Runs combinations of scDblFinder parameters on a set of datasets
#'
#' @param ll Named list of datasets (in SCE format, with a 'truth' colData
#' column)
#' @param params Named list of lists, containing alternative parameter values
#' @param eg2 Optional parameter combination table, used to avoid testing all
#' combinations
#' @param seeds integer vector of random seeds
#' @param BPPARAM Optional BiocParallel param for multithreading
#' @param use.precomputed.clusters Logical; whether to use pre-computed
#' clusters (if )
#'
#' @return A data.frame of doublet accuracy metrics across datasets, parameters
#' combinations and seeds
runParams <- function(ll, params=list(
      nrounds=list(NULL,40),
      includePCs=list(1:3, c()),
      max_depth=list(5),
      propRandom=list(0)
    ), eg2=NULL,
    seeds=c(1234, 42), BPPARAM=NULL, use.precomputed.clusters=TRUE,
    FN=scDblFinder){
  library(BiocParallel)
  if(is.null(eg2)){
    eg <- expand.grid( c(list(seeds=seeds), params) )
    eg2 <- as.data.frame(lapply(as.data.frame(eg), FUN=function(x){
      if(is.list(x)) x <- sapply(x, FUN=toString)
      x
    }))
  }else{
    eg <- eg2
  }
  message("Running ", nrow(eg), " combinations on ", length(ll), " datasets.")
  if(is.null(BPPARAM)) BPPARAM <- SerialParam()
  if(is.numeric(BPPARAM)) BPPARAM <- MulticoreParam(min(nrow(eg),BPPARAM))
  #res <- bplapply(seq_len(nrow(eg)), BPPARAM=BPPARAM, FUN=function(i){
  res <- lapply(seq_len(nrow(eg)), FUN=function(i){
    pp <- lapply(eg[i,], FUN=function(x) x[[1]])
    if(is.null(pp$clusters)) pp <- pp[setdiff(names(pp),"clusters")]
    set.seed(pp$seeds)
    pp$seeds <- NULL
    x <- sapply(ll, FUN=function(x){
      if(!is.null(pp$propKnown) && pp$propKnown>0){
        wD <- which(x$truth=="doublet")
        x$known <- FALSE
        x$known[sample(wD, floor(length(wD)*pp$propKnown))] <- TRUE
        pp$knownDoublets <- "known"
        set.seed(pp$seeds)
      }
      pp$propKnown <- NULL
      isTraj <- !is.null(x$is.traj) && x$is.traj
      if( !is.null(pp$clusters) || !use.precomputed.clusters || is.null(ll[[1]]$cluster)){
        wrapper <- function(...) FN(x, ...)
      }else{
        wrapper <- function(...) FN(x, clusters=x$cluster, ...)
      }
      st <- system.time(x <- do.call(wrapper, pp))
      if("known" %in% colnames(colData(x))) x <- x[,-which(x$known)]
      if("scDblFinder.score" %in% colnames(colData(x))){
        s <- split(x$scDblFinder.score, x$truth)
        x$call <- x$scDblFinder.class=="doublet"
      }else{
        s <- split(x$directDoubletScore, x$truth)
        x$call <- x$directDoubletScore>0.5
      }

      pp <- PRROC::pr.curve(s[[1]], s[[2]])
      x$truth <- x$truth=="doublet"
      c(AUPRC=mean(as.numeric(pp[2:3])),
        AUROC=PRROC::roc.curve(s[[1]], s[[2]])[[2]],
        accuracy=sum(x$truth==x$call)/ncol(x),
        TP=sum(x$truth & x$call),
        FP=sum(!x$truth & x$call),
        FN=sum(x$truth & !x$call),
        elapsed=st[[3]])
    })
    data.frame(dataset=colnames(x), t(x))
  })
  res <- cbind(eg2[rep(seq_along(res), sapply(res, nrow)),],
               dplyr::bind_rows(res))
  res$FDR <- res$FP/(res$TP+res$FP)
  res$sensitivity <- res$TP/(res$TP+res$FN)
  rs <- rowsum(res[,c("AUPRC","AUROC")],res$dataset)/nrow(eg)
  tmp <- res[,c("AUPRC","AUROC")]-rs[as.character(res$dataset),]
  colnames(tmp) <- paste0(colnames(tmp),".diff")
  cbind(res, tmp)
}


dblTypesScheme <- function(){
  library(ggplot2)
  set.seed(123)
  d <- data.frame(x=c(rnorm(20,mean=0),rnorm(10,mean=5)),
                  y=c(rnorm(10,mean=5),rnorm(10,mean=0),rnorm(10,mean=0.5)),
                  celltype=rep(LETTERS[1:3],each=10),
                  genotype=sample(rev(letters)[1:3], 30, replace=TRUE))
  wa1 <- which(d$celltype=="A" & d$genotype=="y" & d$x>1)
  wa2 <- which(d$celltype=="C" & d$genotype=="y" & d$y>2)
  wb1 <- which(d$celltype=="A" & d$genotype=="x" & d$y<=min(d$y[which(d$celltype=="A")]))
  wb2 <- which(d$celltype=="B" & d$genotype=="y" & d$y>2)
  wc1 <- which(d$celltype=="C" & d$genotype=="z")
  wc2 <- which(d$celltype=="C" & d$genotype=="x" & d$x>6)

  d2 <- data.frame( x=d$x[c(wa1,wb1,wc1)], xend=d$x[c(wa2,wb2,wc2)],
                    y=d$y[c(wa1,wb1,wc1)], yend=d$y[c(wa2,wb2,wc2)],
                    type=c("within-genotype\nheterotypic", "inter-genotype\nheterotypic",
                           "inter-genotype\nhomotypic") )
  d2$xmean <- (d2$x+d2$xend)/2
  d2$ymean <- (d2$y+d2$yend)/2
  ggplot(d,aes(x,y)) + geom_segment(data=d2, aes(xend=xend,yend=yend), linetype="dashed") +
    geom_point(data=d, aes(shape=celltype, colour=genotype), size=4) +
    labs(x="PC1", y="PC2") + theme_minimal(base_size=10) +
    geom_text(data=d2, aes(x=xmean, y=ymean, label=type), nudge_y=c(0,0,1),
              nudge_x=c(1.1,-1.2,0)) + xlim(range(d$x)+c(-0.2,0.2)) +
    guides(colour=guide_legend(title.position="top"), shape=guide_legend(title.position="top")) +
    theme(axis.ticks=element_blank(), axis.text=element_blank(), legend.position="bottom")
}
