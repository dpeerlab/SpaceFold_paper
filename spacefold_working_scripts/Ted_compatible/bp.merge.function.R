#functions used to merge / split mulitple bayesprism outputs


#' merge the data entry $para$X in the Bayesprism output
#' @param bp.list a named (optional) list of Bayesprism output
merge.X <- function(bp.list){
	
	if(is.null(names(bp.list))) names(bp.list) <- paste("sample_",1:length(bp.list),sep="")

	gene.its <- Reduce(intersect,lapply(bp.list, function(i) colnames(i$para$X) ))
		
	X.merged <- do.call(rbind, lapply(1:length(bp.list), function(idx) {
		idx.X <- bp.list [[idx]]$para$X[, gene.its]
		rownames(idx.X) <- paste(names(bp.list)[idx], rownames(idx.X),sep="-")
		idx.X
	}) )
	
	X.merged
}


#' merge the data entry $para$X$meta in the Bayesprism output
#' @param bp.list a named (optional) list of Bayesprism output
merge.meta <- function(bp.list){
	
	if(is.null(names(bp.list))) names(bp.list) <- paste("sample_",1:length(bp.list),sep="")

	colname.its <- Reduce(intersect,lapply(bp.list, function(i) colnames(i$para$meta) ))
		
	meta.merged <- do.call(rbind.data.frame, 
						   lapply(1:length(bp.list), 
							function(idx) bp.list [[idx]]$para$meta[, colname.its]
					))
	
	meta.merged
}

#' merge the data entry $para$X$feature in the Bayesprism output
#' @param bp.list a named (optional) list of Bayesprism output
merge.feature <- function(bp.list){
	
	if(is.null(names(bp.list))) names(bp.list) <- paste("sample_",1:length(bp.list),sep="")

	rowname.its <- Reduce(intersect,lapply(bp.list, function(i) rownames(i$para$feature) ))
		
	feature.merged <- do.call(rbind.data.frame, 
						   lapply(1:length(bp.list), 
							function(idx) bp.list [[idx]]$para$feature[rowname.its,]
					))
	
	feature.merged
}


merge.2d.matrix <- function(matrix.list){
	cell.type.its <- Reduce(intersect,lapply(matrix.list, colnames))

	do.call(rbind, lapply(1:length(matrix.list), 
							function(idx) {
								matrix.i <- matrix.list[[idx]][, cell.type.its]
								rownames(matrix.i) <- paste(names(matrix.list)[idx], rownames(matrix.i),sep="-")
								matrix.i}) )
}

merge.3d.array <- function(array.list){
	cell.type.its <- Reduce(intersect,lapply(array.list, function(array.i) dimnames(array.i)[[2]]))
	gene.its <- Reduce(intersect,lapply(array.list, function(array.i) dimnames(array.i)[[3]]))
	
	sample.id <- unlist(lapply(1:length(array.list),function(i) paste(names(array.list)[i], dimnames(array.list[[i]])[[1]],sep="-")))
	
	array.merged <- array(0,
						  dim=c(length(sample.id),length(cell.type.its),length(gene.its)),
						  dimnames = list(sample.id, cell.type.its, gene.its))
	 
	sample.size.vec <- unlist(lapply(array.list, function(array.list.i) dim(array.list.i)[1] ))	 
	sample.size.vec <- c(0, sample.size.vec)
	
	for(i in 1:length(array.list)){
		sample.idx <- sample.size.vec[i]+1 : sample.size.vec[i+1]
		array.merged[sample.idx,,] <- array.list [[i]][, cell.type.its, gene.its]
	}
	array.merged

}

#' merge theta in the Bayesprism output
#' @param bp.list a named (optional) list of Bayesprism output
#' @param which.theta indicates whether merge first or final theta.
#' @param state.or.type indicates whether merge cell state or cell type (only works for which.theta="first").
merge.theta <- function(bp.list, 
						which.theta=c("first","final"), 
						state.or.type=c("state","type")){
	
	if(is.null(names(bp.list))) names(bp.list) <- paste("sample_",1:length(bp.list),sep="")

	theta.list <- sapply(bp.list, function(i) {
												if(which.theta =="first") {
													if(state.or.type =="state") theta <- i$res$first.gibbs.res$gibbs.theta
													if(state.or.type =="type") theta <- i$res$first.gibbs.res$theta.merged
												}
												if(which.theta =="final") theta <- i$res$final.gibbs.theta
												return(theta)
											})
	
	merge.2d.matrix(theta.list)
	
}

#' merge Znk in the Bayesprism output
#' @param bp.list a named (optional) list of Bayesprism output
#' @param state.or.type indicates whether merge cell state or cell type (only works for which.theta="first").
merge.Znk <- function(bp.list, 
						state.or.type=c("state","type")){
	
	if(is.null(names(bp.list))) names(bp.list) <- paste("sample_",1:length(bp.list),sep="")

	Znk.list <- sapply(bp.list, function(i) {if(state.or.type =="state") return( i$res$first.gibbs.res$Znk )
											 if(state.or.type =="type") return( i$res$first.gibbs.res$Znk.merged)
											})

	merge.2d.matrix(Znk.list)
}


#' merge Znkg in the Bayesprism output (along the n direction)
#' @param bp.list a named (optional) list of Bayesprism output
#' @param state.or.type indicates whether merge cell state or cell type.
merge.Znkg <- function(bp.list,
						raw.or.norm=c("raw","norm"), 
						state.or.type=c("state","type")){
	
	Znkg.list <- sapply(bp.list, function(i) {if(raw.or.norm=="raw" & state.or.type =="state") return( i$res$first.gibbs.res$Znkg )
											  if(raw.or.norm=="raw" & state.or.type =="type") return( i$res$first.gibbs.res$Znkg.merged)
											  if(raw.or.norm=="norm") return( i$res$first.gibbs.res$Znkg.merged.normed)
											})
	merge.3d.array(Znkg.list)
}


#' merge a list of Bayesprism output to generate a single object in the format of Bayesprism output 
#' @param bp.list a (named) list of BayesPrism object. If name is provided, it will be used as prefix for theta and Znkg. Otherwise a name in the format of sampleN_** is added.
merge.bp <- function(bp.list){
	
	#only merge the important entries in the res
	bp.merged <- list(para=list(X = merge.X(bp.list), 
								meta = merge.meta(bp.list)),
					  res=list(first.gibbs.res=list(gibbs.theta = merge.theta(bp.list, which.theta ="first", state.or.type="state"),
					  								Znkg = merge.Znkg(bp.list, raw.or.norm="raw", state.or.type="state"),
					  								theta.merged = merge.theta(bp.list, which.theta ="first", state.or.type="type"),
					  								Znkg.merged = merge.Znkg(bp.list, raw.or.norm="raw", state.or.type="type"))))	
	
	if(sum(sapply(bp.list,function(i)!is.null(i$para$feature)))==length(bp.list)){
		bp.merged$para$feature <- merge.feature(bp.list)
	}
	
	if(sum(sapply(bp.list,function(i)!is.null(i$res$final.gibbs.theta)))==length(bp.list)){
		bp.merged$res$final.gibbs.theta <-  merge.theta(bp.list, which.theta ="final", state.or.type="type")
	}
	
	if(sum(sapply(bp.list,function(i)!is.null(i$res$first.gibbs.res$Znk.merged)))==length(bp.list)){
		bp.merged$res$first.gibbs.res$Znk <- merge.Znk (bp.list,"state")
		bp.merged$res$first.gibbs.res$Znk.merged <- merge.Znk (bp.list,"type")
	}
	
	if(sum(sapply(bp.list,function(i)!is.null(i$res$first.gibbs.res$Znkg.merged.normed)))==length(bp.list)){
		bp.merged$res$first.gibbs.res$Znkg.merged.normed <- merge.Znkg(bp.list, raw.or.norm="norm", state.or.type="type")
	}
	
	bp.merged
}

#' split a merged BayesPrism ouput by an column of meta; this function is used to generate separate bp.obj for making spatial plot
#' @param by the colname of meta used as key to split; default separate by sample
split.bp <- function(bp.obj,
					 by="sample"){
	
	stopifnot(!is.null(bp.obj$para$meta))
	
	uniq.keys <- unique(bp.obj$para$meta[,by])
	bp.list <- lapply(uniq.keys, function(key.i) subset.bp(bp.obj, 
															col.name=by, 
															val= key.i, 
															operator="=="))
	names(bp.list) <- uniq.keys
	
	bp.list
}





#testing functions:

# load("/lila/data/peer/tinyi/RU_data/visium/dat/bp_dat/from_reprocessed_anndata/LI/merge_finer_states_0930/WTLI1.bp.rdata")
# load("/lila/data/peer/tinyi/RU_data/visium/dat/bp_dat/from_reprocessed_anndata/LI/merge_finer_states_0930/WTLI2.bp.rdata")

# load("/lila/data/peer/tinyi/RU_data/visium/dat/visium_dat/LI/visium.LI.WhiteSketch.rdata")

# WTLI1.bp <- add.meta (bp.obj=WTLI1.bp, meta = bcs_merge, selected.sample="WTLI1")
# WTLI2.bp <- add.meta (bp.obj=WTLI2.bp, meta = bcs_merge, selected.sample="WTLI2")

# label.df <- read.csv("/data/peer/tinyi/RU_data/visium/dat/visium_dat/region_labeled/WTLI1_region.csv")
# WTLI1.bp <- add.region.labels(WTLI1.bp, label.df)

# label.df <- read.csv("/data/peer/tinyi/RU_data/visium/dat/visium_dat/region_labeled/WTLI2_region.csv")
# WTLI2.bp <- add.region.labels(WTLI2.bp, label.df)

# # WTLI.bp <- list(WTLI1.bp, WTLI2.bp)
# # names(WTLI.bp) <- c("WTLI1","WTLI2")


# WTLI.bp.merged <-  merge.bp(WTLI.bp)
# WTLI.bp.split <-  split.bp(WTLI.bp.merged)

# a <-  merge.theta(WTLI.bp, which.theta ="first", state.or.type="state")
# a <-  merge.theta(WTLI.bp, which.theta ="first", state.or.type="type")

# a <- merge.X(WTLI.bp)
# a <- merge.meta(WTLI.bp)
