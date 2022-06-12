
library(phateR)
library(RColorBrewer)
library(beeswarm)

#' functions to run PHATE on BayesPrism output
#' default parameters were those used by the spacefold paper
#' @param bp.obj: a (merged) BayesPrism output  
#' @param which.theta indicates whether use first or final theta. Default is final theta.
#' @param anchorCellTypes: a character vector denoting the cell types selected for computing embeddings. 
#'			only cell types show common the spatial pattern across all bp.obj in the bp.list should be selected.
#'			anchorCellTypes should be informative about the physical cordinate of the spot, 
#'			i.e., cell types that are known or suspected to vary with the cordinate.
#'			Spacefold is usually robust to selection of cell types.
#'			#default="all', i.e. uses all cell types
#' @param renorm.to.one, nomalize the seleceted cell type fractions to sum to one. Default TRUE.
#' @param center, scale the seleceted cell type fractions across spatial spots to mean=0. Default TRUE.
#' @param scale, scale the seleceted cell type fractions across spatial spots to sd=1. Default TRUE.
#' @param if.pseudo.axis, return a scaled cordinate with min=zero and max=one.
#' @param if.invert whehter flip the spacefold cordinate. Default is FALSE 
#' @param mds.solver solver used by PHATE
#' @param n.jobs number of threads used by PHATE
#' @param knn the knn paramter used by PHATE
#' @param ... additional paramters passed to PHATE
phate.bp <- function(bp.obj,
					 which.theta='final',
					 anchorCellTypes ="all",
					 renorm.to.one = TRUE,
					 center=TRUE,
					 scale=TRUE,
					 if.pseudo.axis=TRUE,
					 if.invert=FALSE,
					 mds.solver="smacof",
					 n.jobs=20,
					 knn=10,
					 ...){
	
	stopifnot(!is.null(bp.obj$para$meta))

	#other assertion arugments to be written...
	
	if(which.theta =="first") theta <- bp.obj $res$first.gibbs.res$gibbs.theta
	if(which.theta =="final") theta <- bp.obj $res$final.gibbs.theta
	
	if(anchorCellTypes=="all") anchorCellTypes <- colnames(theta) 
	
	#subset using relevant cell types
	theta <- theta[, anchorCellTypes]
	
	if(renorm.to.one) theta <- theta / rowSums(theta)
	
	theta <- scale(theta, center = center, scale = scale)
	
	#can add addtional PCA steps...
	
	phate.res <- phate(data= theta, 
					   ndim= 1, 
					   mds.solver= mds.solver,
					   n.jobs= n.jobs,
					   knn= knn,
					   ...)
	
	cord <- phate.res$embedding
	
	if(if.invert) cord <- -cord
	
	if(if.pseudo.axis) {
		cord <- apply(cord,2,function(cord.i){
								cord.i <- cord.i - min(cord.i)
								cord.i <- cord.i / max(cord.i)
								cord.i
					})
	}
	
	#store spacefold parameters
	bp.obj $para $spacefold <- list(which.theta= which.theta,
									anchorCellTypes= anchorCellTypes,
									renorm.to.one= renorm.to.one,
									center= center,
									scale= scale,
									if.pseudo.axis= if.pseudo.axis,
									ndim=1,
					 				mds.solver= mds.solver,
					 				n.jobs= n.jobs,
					 				knn= knn,
					 				list(...))
		
	#store at the meta entry   
	bp.obj $para$meta$PHATE1 <-  cord
	bp.obj
}


#test:
#WTLI.bp.merged <- phate.bp (WTLI.bp.merged)
#WTLI.bp.merged.top <- subset.bp(WTLI.bp.merged, col.name="region",val="top")

#WTLI.bp.merged.top <- phate.bp (WTLI.bp.merged.top)


#' function to plot beeswarm
#' @param palette a function that generates the palette 
#' @param q.cut upper quantile of spots to plot, default is 0.95 
#'			(=spots containing the top 5% of each cell type) 
#'		if set to "auto", use the $res$background.level
#' @param pdf.prefix the prefix of pdf filename.
#' @param cell.type.order a character vector speficying the order of plotting. 
#'			Default is NULL, which uses the colnames of theta matrix
plot.beeswarm <- function(bp.obj,
						 palette=colorRampPalette(brewer.pal(12, "Paired")),
						 q.cut=0.95,
						 pdf.prefix="output",
						 height=5,
						 width=12,
						 mar=c(8,3,1,1),
						 cell.type.order=NULL){
	
	cord <- as.vector(bp.obj$para$meta$PHATE1)
	stopifnot(!is.null(cord))
		
	which.theta <- bp.obj$para$spacefold$which.theta
	if(which.theta =="first") theta <- bp.obj $res$first.gibbs.res$theta.merged
	if(which.theta =="final") theta <- bp.obj $res$final.gibbs.theta	
	
	Znk <- bp.obj $res$first.gibbs.res$Znk.merged
	
	my.palette <-  palette(ncol(theta))

	cell.type.uniq <-colnames(theta)

	plot.df <- do.call(rbind.data.frame, 
							lapply(1:length(cell.type.uniq), function(idx){
	
							theta.vec <- theta[, cell.type.uniq[idx]]
							Znk.vec <- Znk[, cell.type.uniq[idx]]
							
							if(q.cut=="auto") {
								stopifnot(!is.null(bp.obj$res$background.level))
								theta.cutoff <- bp.obj$res$background.level$theta.cutoffs[cell.type.uniq[idx]]
								Znk.cutoff <- bp.obj$res$background.level$Znk.cutoff[cell.type.uniq[idx]]
							}
							else{
								theta.cutoff <- quantile(theta.vec, prob= q.cut)
								Znk.cutoff <- quantile(Znk.vec,prob= q.cut)
							}
									#browser()					
							data.frame(layout="PHATE",
									   cord= cord, 
									   cell.type = cell.type.uniq[idx],
									   theta = theta.vec,           
									   if.above.background =  theta.vec > theta.cutoff & Znk.vec > Znk.cutoff,
									   col.solid = my.palette[idx],
									   col.theta = sapply(  theta.vec, adjustcolor, col=my.palette[idx]),
									   stringsAsFactors=FALSE)
	}))

	if(!is.null(cell.type.order)) plot.df $cell.type <- factor(plot.df $cell.type,levels= cell.type.order)

	dev.new(height= height, width= width, file=paste(pdf.prefix,"_spacefold_beeswarm.pdf",sep=""))

	par(mar= mar, xpd=TRUE)

	beeswarm(cord ~ cell.type, data = plot.df[plot.df $if.above.background,],
         pch = 16, pwcol = as.character(col.solid),
         xlab = "",
         ylab = "projected coordinate,colored by selected theta",
         cex=0.5,las=2,cex.axis=0.5,
         corral="wrap",
         bty="n",
         main=paste("theta above", q.cut, "quantile"))

	beeswarm(cord ~ cell.type, data = plot.df,
         pch = 16, pwcol = as.character(col.theta),
         xlab = "",
         ylab = "projected coordinate",
         cex=0.5,las=2,cex.axis=0.5,
         corral="wrap",
         bty="n",
         main=paste("color by theta"))
         
	dev.off()

}

#plot.beeswarm (WTLI.bp.merged.top)

