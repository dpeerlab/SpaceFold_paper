#functions to project gene expression onto the spacefold axis

library("mixtools")
library("mclust")


#set background level using user provided values
#' @param which.theta use first(initial) or final(updated) theta 
#' @param theta.cut cutoff for theta of background spots, a numerical number or a named numerical vector of the same length as ncol(theta).
#' @param Znk.cut cutoff for Znk of background spots, a numerical number or a named numerical vector of the same length as ncol(theta).
set.background.level <-  function(bp.obj, 
								  which.theta=c('first','final'), 
								  theta.cutoffs, 
								  Znk.cutoffs){
		
	Znk <- bp.obj$res$first.gibbs.res$Znk.merged
	if(which.theta=="first") theta <- bp.obj$res$first.gibbs.res$theta.merged
	if(which.theta=="final") theta <- bp.obj$res$final.gibbs.theta
	
	if(is.null(Znk.merged)) {
		bp.obj <- compute.Znk(bp.obj)
		Znk <- bp.obj$res$first.gibbs.res$Znk.merged
	}
	
	if(length(theta.cutoffs)==1) {
		theta.cutoffs <- rep(theta.cutoffs, ncol(theta))
		names(theta.cutoffs) <- colnames(theta)
	}
	if(length(Znk.cutoffs)==1) {
		Znk.cutoffs <- rep(Znk.cutoffs, ncol(theta))
		names(Znk.cutoffs) <- colnames(theta)
	}
	
	stopifnot(sum(names(theta.cutoffs)==colnames(theta))==length(theta.cutoffs))
	stopifnot(sum(names(Znk.cutoffs)==colnames(Znk))==length(Znk.cutoffs))
	
	selected.spot.matrix <- do.call(cbind, 
									lapply(1:ncol(theta), 
											function(col.idx) {
														theta[, col.idx] > theta.cutoffs[col.idx] & 
														Znk[, col.idx] > Znk.cutoffs[col.idx] }))
	colnames(selected.spot.matrix) <- colnames(theta)
	
	bp.obj$res$background.level <- list(theta.cutoffs= theta.cutoffs, 
										Znk.cutoffs= Znk.cutoffs,
										selected.spot.matrix= selected.spot.matrix,
										which.theta= which.theta)
	
	bp.obj
}




#set background level by fitting a mixture model on theta and Znk (for each cell type)
#' @param bp.obj a BayesPrism output object
#' @param which.theta use first(initial) or final(updated) theta 
compute.background.level <-  function(bp.obj, which.theta=c('first','final')){
		
	Znk <- bp.obj$res$first.gibbs.res$Znk.merged
	if(which.theta=="first") theta <- bp.obj$res$first.gibbs.res$theta.merged
	if(which.theta=="final") theta <- bp.obj$res$final.gibbs.theta
	
	if(is.null(Znk)) {
		bp.obj <- compute.Znk(bp.obj)
		Znk <- bp.obj$res$first.gibbs.res$Znk.merged
	}
	
	theta.cutoffs <-c()
	print("fitting mixture models on theta...")
	for(i in 1:ncol(theta)){
		print(colnames(theta)[i])

		fit=gammamixEM(theta[,i], k=2,maxit=10000,maxrestarts=100)
		cls <- apply(fit$posterior,1,which.max)
		if(length(unique(cls))==1) {
			print("refit")
			fit=gammamixEM(theta[,i], k=2,maxit=10000,maxrestarts=100, mom.start=F)
			cls <- apply(fit$posterior,1,which.max)
		}
		print(table(cls))
		theta.cutoffs <- c(theta.cutoffs, median(c(range(theta[cls==1,i]),range(theta[cls==2,i]))))
	}
	names(theta.cutoffs) <- colnames(theta)

	Znk.cutoffs <-c()
	print("fitting mixture models on Znk...")
	for(i in 1:ncol(Znk)){
		print(colnames(Znk)[i])
	
		fit= Mclust(Znk[,i], G=2)
		cls <- fit$classification
		print(table(cls))
		Znk.cutoffs <- c(Znk.cutoffs, median(c(range(Znk[cls==1,i]),range(Znk[cls==2,i]))))
	}
	names(Znk.cutoffs) <- colnames(Znk)
	
	selected.spot.matrix <- do.call(cbind, 
									lapply(1:ncol(theta), 
											function(col.idx) {
														theta[, col.idx] > theta.cutoffs[col.idx] & 
														Znk[, col.idx] > Znk.cutoffs[col.idx] }))
	colnames(selected.spot.matrix) <- colnames(theta)
	
	bp.obj$res$background.level <- list(theta.cutoffs= theta.cutoffs, 
										Znk.cutoffs= Znk.cutoffs,
										selected.spot.matrix= selected.spot.matrix,
										which.theta= which.theta)
	
	bp.obj
}

#test
#WTLI.bp.merged <- compute.background.level(WTLI.bp.merged, "final")



get.gene.idx <- function(bp.obj, selected.genes){
	#determine if selected.genes is symbol or ID, 
	#and return a vector gene.idx denoting the index of selected genes.
	#search by maximum overlapping column
	which.col <- which.max(apply(bp.obj$para$feature,2,
								 function(col.i) length(intersect(as.character(col.i), selected.genes))))
	gene.idx <- match(selected.genes, bp.obj$para$feature[, which.col]) 	
	if(sum(is.na(gene.idx))>0) stop(paste(c("selected.genes", selected.genes[is.na(gene.idx)], "were not found"), collapse=" "))
	gene.idx	
}


#return theta, Znkg array, and the cutoffs 
subset.celltype <- function(bp.obj,
							raw.or.norm=c("raw","norm"),
							selected.cell.types){
	
	#determine if using the original cell types or the regrouped cell types
	#and return a vector celltype.idx denoting the index of selected cell types.
	cell.types.raw <- colnames(bp.obj$res$first.gibbs.res$theta.merged)
	cell.types.rgp <- colnames(bp.obj$res$res.regrouped$theta0)
	
	raw.matched <- length(intersect(cell.types.raw, selected.cell.types))
	rgp.matched <- length(intersect(cell.types.rgp, selected.cell.types))
	
	if(raw.matched >= rgp.matched){
		theta0 <- bp.obj$res$first.gibbs.res$theta.merged
		thetaF <- bp.obj$res$final.gibbs.theta
		Znk <- bp.obj$res$first.gibbs.res$Znk.merged
		Znkg <- bp.obj$res$first.gibbs.res$Znkg.merged
		Znkg.normed <- bp.obj$res$first.gibbs.res$Znkg.merged.normed
		theta.cutoffs <- bp.obj$res$background.level$theta.cutoffs
		background.info <- bp.obj$res$background.level
		cell.type.idx <- match(selected.cell.types, cell.types.raw)
	}
	else{
		theta0 <- bp.obj$res$res.regrouped $theta0
		thetaF <- bp.obj$res$res.regrouped $thetaF
		Znk <- bp.obj$res$res.regrouped$Znk
		Znkg <- bp.obj$res$res.regrouped$Znkg
		Znkg.normed <- bp.obj$res$res.regrouped$Znkg.normed
		background.info <- bp.obj$res$res.regrouped$background.level
		cell.type.idx <- match(selected.cell.types, cell.types.rgp)
	}
	
	stopifnot(sum(is.na(cell.type.idx))==0)

	which.theta <- background.info$which.theta
	if(which.theta=="first") theta <- theta0[, cell.type.idx, drop=F]
	if(which.theta=="final") theta <- thetaF[, cell.type.idx, drop=F]

	Znk <- Znk[, cell.type.idx]
	
	if(raw.or.norm=="raw") Znkg <- Znkg[, cell.type.idx, , drop=F]
	if(raw.or.norm=="norm") Znkg <- Znkg.normed[, cell.type.idx, , drop=F]
	
	selected.spot.matrix <- background.info$selected.spot.matrix
	
	return(list(Znkg= Znkg, 
				theta= theta,
				Znk =  Znk,
				selected.spot.matrix= selected.spot.matrix,
				cord=bp.obj $para$  meta [,'PHATE1']))
								
}


#' function that fits the (binned) expression value against the (binned) cordinate
#' @param exp.vec a numeric vector of expression value in selected spots of one gene in one cell type
#' @param cord.vec a numeric vector of spacefold cordinate in selected spots. same length as exp.vec
#' @param n.bins number of bins, default=20
#' @param bin.by method to bin the data, deviding them into each bin of each sample size (equal.size), 
		# or by the same distance over spacefold axis (equal.step)
#' @param span parameter used by loess curve fitting, default=0.75
fit.exp.curv <- function(exp.vec, 
						 cord.vec, 
						 n.bins=20,
						 bin.by="equal.size",
						 span = 0.75){
	
	quantile.cut.points <- seq(0,1,length.out= n.bins+1 )

	if (bin.by=="equal.size"){
		breaks= quantile(cord.vec,prob= quantile.cut.points)
	}
	if (bin.by=="equal.step"){
		breaks= quantile(cord.vec,prob= quantile.cut.points) 
	}
	
	cord.cuts <- cut(cord.vec,
							breaks= breaks,
							labels=paste("cord",quantile.cut.points[-1],sep="-"),include.lowest=TRUE)
	exp.bin.mean <- by(exp.vec, INDICES= cord.cuts, mean)
	exp.bin.se <- by(exp.vec, INDICES= cord.cuts, function(i) sd(i)/sqrt(length(i)) ) 
	exp.bin.mean.up <- exp.bin.mean +2* exp.bin.se
	exp.bin.mean.low <- exp.bin.mean -2* exp.bin.se

	cord.bin.mid <- by(cord.vec, INDICES= cord.cuts, mean)	
	
	#spline fit
	grid.point <- seq(min(cord.bin.mid,na.rm=T),max(cord.bin.mid,na.rm=T),length.out=1000)

	loess_fit.mean <- predict(loess(exp.bin.mean ~ cord.bin.mid, span= span,na.action="na.omit"), grid.point)
	loess_fit.up <- predict(loess(exp.bin.mean.up ~ cord.bin.mid, span= span,na.action="na.omit"), grid.point)
	loess_fit.low <- predict(loess(exp.bin.mean.low ~ cord.bin.mid, span= span,na.action="na.omit"), grid.point)
	
	return(list(exp.bin.mean = exp.bin.mean,
			    cord.bin.mid = cord.bin.mid, 
			    loess_fit.mean= loess_fit.mean, 
			    loess_fit.up= loess_fit.up, 
			    loess_fit.low= loess_fit.low,
			    grid.point= grid.point,
			    exp.vec= exp.vec,
			    cord.vec= cord.vec))
}




#' function that obtains raw expression data, cordinates, binned data, and fitted curves to plot (after applying the filter)
#' returns a list (over genes) of list (over cell types) of exp.curv object (outputed by fit.exp.curv)
get.plot.dat <- function(bp.obj,
						raw.or.norm=c("raw","norm"),
						selected.genes,
						selected.cell.types,
						bin.by,
						n.bins,
						span){
	
	gene.idx <- get.gene.idx(bp.obj= bp.obj, selected.genes= selected.genes)
	 
	cell.types.dat <- subset.celltype(bp.obj= bp.obj, 
									  raw.or.norm= raw.or.norm,
									  selected.cell.types= selected.cell.types)
	#browser()
	dat.list <- list()
	for(g in 1:length(selected.genes)){
		ct.list<- list() 
		for(k in 1:length(selected.cell.types)){
			spot.idx <- cell.types.dat $selected.spot.matrix[, selected.cell.types[k]]
			exp.vec <- cell.types.dat$Znkg[spot.idx, selected.cell.types[k], gene.idx[g]]
			cord.vec <- cell.types.dat$cord[spot.idx]
			
			ct.list[[k]] <- fit.exp.curv (exp.vec= exp.vec, 
						 								 cord.vec= cord.vec, 
						 								 n.bins= n.bins,
						 								 bin.by= bin.by,
						 								 span = span)
			
		}
		names(ct.list) <- selected.cell.types
		dat.list[[g]] <- ct.list
	} 
	names(dat.list) <- selected.genes
	
	dat.list
}


get.ylim <- function(dat.list,
					 gene.idx,
					 show.raw){
	
	dat.list.gene <- dat.list[[gene.idx]]
	
	fit.dat.vec <- c( unlist(lapply(dat.list.gene, '[[', 'loess_fit.up')),
				  unlist(lapply(dat.list.gene, '[[', 'loess_fit.low'))	 
				 )
					 	
	if(show.raw) range(c(fit.dat.vec, unlist(lapply(dat.list.gene, '[[', 'exp.vec'))))
	else range(fit.dat.vec)				 	
	
}



#' function that plots expression curve
#' @param ct.name character variable of cell type name
#' @param gene.name character variable of gene name
#' @param show.raw logical variable of whether overlay the unbinned Znkg value
#' @param color the color of the fitted curves
#' @param color.raw the color of the raw points
#' @param pch.raw pch of raw points
#' @param ylim limit of y axis, default=NULL (automatically determined)
#' @param xlim limit of x axis, dafault=c(0,1),
#' @param ... other arguments passed to plot 
plot.exp.curv <- function(exp.curv,
						  ct.name=NULL, 
						  gene.name=NULL,
						  show.raw=FALSE,
						  color="black",
						  color.raw= adjustcolor("black",0.75),
						  cex.raw=0.7,
						  pch.raw=16,
						  ylim=NULL,
						  xlim=c(0,1),
						  xlab=NULL,
						  ylab=NULL,
						  axis.pos="left", #c("left","right")
						  ...){
	
	if(is.null(ylim)) ylim <- range(exp.curv$exp.bin.mean)

	if(axis.pos=="left"){
		plot(exp.curv$exp.bin.mean ~ exp.curv$cord.bin.mid,
		 	type="n",
		 	xlim= xlim, 
		 	ylim= ylim,
		 	bty="n",
		 	main=paste(gene.name, "in", ct.name, sep=" "),
		 	xlab=xlab,
		 	ylab= ylab,
		 	...)
	}
	else{
		plot(exp.curv$exp.bin.mean ~ exp.curv$cord.bin.mid,
		 type="n",
		 xlim= xlim, 
		 ylim= ylim,
		 bty="n",
		 main=paste(gene.name, "in", ct.name, sep=" "),
		 xlab=xlab,
		 ylab="",
		 axes=FALSE,
		 ...)
		 mtext(ylab,side=4,line=4,cex=2) 
		axis(4, ylim= ylim)
	}

	
	polygon(c(exp.curv$grid.point,rev(exp.curv$grid.point)),
		c(exp.curv$loess_fit.low,rev(exp.curv$loess_fit.up)),
		border=NA,col=adjustcolor(color,0.25))
	lines(exp.curv$grid.point, exp.curv$loess_fit.mean,lwd=1,col= color)
	points(x= exp.curv$cord.bin.mid, y=exp.curv$exp.bin.mean,pch=16,col= color)
	if(show.raw) points(x= exp.curv$cord.vec, y=exp.curv$exp.vec,pch= pch.raw,col= color.raw,cex= cex.raw)
	NULL
}



get.pdf.name <- function(pdf.idx){
	if(pdf.idx>= 10 & pdf.idx < 100) pdf.name <- paste("tmp.0", pdf.idx,".pdf",sep="")
	if(pdf.idx < 10) pdf.name <- paste("tmp.00", pdf.idx,".pdf",sep="")
	if(pdf.idx >= 100) pdf.name <- paste("tmp.", pdf.idx,".pdf",sep="")
	pdf.name
}




#visualize gene expression level along the projected axis
# use this function only in pdfjam is installed. Otherwise, use plot.cartography.nopdfjam.
#' @param bp.obj a BayesPrism output object
#' @param raw.or.norm plotting raw value or value normalized by size factor
#' @param selected.genes a character vector specifiying the gene name / symbols to be plotted, which is automatically determined by matching columns of  bp.obj$para$feature 
#' @param selected.cell.types a character vector specifiying the cell types / cell type groups to be plotted, which is automatically determined by matching colnames of bp.obj$res$first.gibbs.res$theta.merged and bp.obj$res$res.regrouped$theta0
#' @param bin.by method to bin the data, deviding them into each bin of each sample size (equal.size), 
		# or by the same distance over spacefold axis (equal.step)
#' @param n.bins number of bins, default=20
#' @param span parameter used by loess curve fitting, default=0.75
#' @param show.raw whether plot the unbinned Znkg/Znkg.norm
#' @param overlay.hist whether plot the histogram of the cordiantes containing each cell type separately or on top of the expression curve.
#' @param col.hist color of the histogram
#' @param pdf.prefix the prefix of the pdf output
plot.cartography.pdfjam <- function(bp.obj,
						raw.or.norm=c("raw","norm"),
						selected.genes,
						selected.cell.types,
						bin.by=c("equal.size","equal.step"),
						n.bins=20,
						span=0.75,
						show.raw =FALSE,
						overlay.hist=TRUE,
						col.hist=adjustcolor("grey",0.5),
						pdf.prefix,
						...){
	#assertions...
	plot.dat <- get.plot.dat  (bp.obj, raw.or.norm, selected.genes, selected.cell.types, bin.by, n.bins, span)
	
	color.num <- length(selected.cell.types)
	if(color.num<=8) my.palette <- brewer.pal(color.num,"Dark2")
	if(color.num>8 && color.num<=12) my.palette <- brewer.pal(color.num,"Paired")
	if(color.num > 12)  my.palette <-  colorRampPalette(brewer.pal(12, "Paired"))(color.num)
	

	ylim.all <- lapply(selected.genes, FUN= get.ylim, dat.list= plot.dat, show.raw= show.raw)
	
	
	pdf.idx <- 1
	
	if(!overlay.hist) {
		for (ct.idx in 1:length(selected.cell.types)){
			pdf(get.pdf.name(pdf.idx),useDingbats=F,pointsize=8)
			par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
			par(mar=c(5.1, 5.1, 4.1, 4.1))
			hist(plot.dat[[1]][[ct.idx]]$cord.vec,
				main= selected.cell.types[ct.idx], 
				xlim= c(0,1), 
				breaks=50, 
				col= adjustcolor(my.palette[ct.idx],0.5),
				xlab="SpaceFold cordinate",
				ylab="Cell Type Frequency")
			dev.off()
			pdf.idx <- pdf.idx + 1
			
			for (gene.idx in 1:length(selected.genes)){
				pdf(get.pdf.name(pdf.idx),useDingbats=F,pointsize=8)
				par(mar=c(5.1, 5.1, 4.1, 4.1))
				par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
				plot.exp.curv (plot.dat[[gene.idx]][[ct.idx]],
						  ct.name= selected.cell.types[ct.idx], 
						  gene.name= selected.genes[gene.idx],
						  show.raw= show.raw,
						  ylim= ylim.all[[gene.idx]],
						  color= my.palette[ct.idx],
						  xlab="SpaceFold cordinate",
						  ylab=paste(raw.or.norm, "expression value"),
						  axis.pos= "left")
				dev.off()
				pdf.idx <- pdf.idx + 1
			}
		}
	}
	else{
		for (ct.idx in 1:length(selected.cell.types)){
			for (gene.idx in 1:length(selected.genes)){
				pdf(get.pdf.name(pdf.idx),useDingbats=F,pointsize=8)
				par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
				par(mar=c(6.1, 5.1, 4.1, 6.1))
				hist(plot.dat[[1]][[ct.idx]]$cord.vec,
				 	main= NULL, 
				 	xlim= c(0,1), 
					breaks=50, 
					col= col.hist,
					xlab="SpaceFold cordinate",
					ylab="Cell Type Frequency")
				 par(new=TRUE) 
				 
				 #par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
				 #par(mar=c(5.1, 5.1, 4.1, 2.1))
				 plot.exp.curv (plot.dat[[gene.idx]][[ct.idx]],
						  	ct.name= selected.cell.types[ct.idx], 
						  	gene.name= selected.genes[gene.idx],
						  	show.raw= show.raw,
						  	ylim= ylim.all[[gene.idx]],
						  	color= my.palette[ct.idx],
						  	xlab="SpaceFold cordinate",
						  	ylab=paste(raw.or.norm, "expression value"),
						  	axis.pos= "right")
				 dev.off()
				 pdf.idx <- pdf.idx + 1
			}
		}
	}	

	#merge pdfs
	ncol=length(selected.genes)
	nrow = length(selected.cell.types)
	if(!overlay.hist) ncol <- ncol+1
	system(paste("pdfjam", " ./tmp*.pdf ",  "--nup ", ncol, "x",  nrow , " --landscape  --outfile ", pdf.prefix,".pdf",sep="")) # needs to install pdfjam
	system("rm ./tmp*.pdf")	

	NULL
}


#test

# # selected.genes=c("Reln","Lyve1","Lgr5","Prox1","Pdpn","Itga9","Pdgfra")
# selected.genes=c("ENSMUSG00000025902","ENSMUSG00000025912")

# plot.dat <- get.plot.dat  (bp.obj=WTLI.bp.merged, raw.or.norm="norm", selected.genes, names(grouping.list))
# plot.cartography.pdfjam (bp.obj=WTLI.bp.merged,
							 # raw.or.norm="raw",
							 # selected.genes= selected.genes,
							 # selected.cell.types=names(grouping.list),
							 # bin.by="equal.size",
							 # show.raw =TRUE,
							 # overlay.hist=TRUE,
							 # pdf.prefix="test.raw")


#visualize gene expression level along the projected axis
# make individual plot and output in a single pdf
#' @param bp.obj a BayesPrism output object
#' @param raw.or.norm plotting raw value or value normalized by size factor
#' @param selected.genes a character vector specifiying the gene name / symbols to be plotted, which is automatically determined by matching columns of  bp.obj$para$feature 
#' @param selected.cell.types a character vector specifiying the cell types / cell type groups to be plotted, which is automatically determined by matching colnames of bp.obj$res$first.gibbs.res$theta.merged and bp.obj$res$res.regrouped$theta0
#' @param bin.by method to bin the data, deviding them into each bin of each sample size (equal.size), 
		# or by the same distance over spacefold axis (equal.step)
#' @param n.bins number of bins, default=20
#' @param span parameter used by loess curve fitting, default=0.75
#' @param show.raw whether plot the unbinned Znkg/Znkg.norm
#' @param overlay.hist whether plot the histogram of the cordiantes containing each cell type separately or on top of the expression curve.
#' @param col.hist color of the histogram
#' @param pdf.prefix the prefix of the pdf output
plot.cartography.nopdfjam <- function(bp.obj,
							 raw.or.norm=c("raw","norm"),
							 selected.genes,
							 selected.cell.types,
							 bin.by=c("equal.size","equal.step"),
							 n.bins=20,
							 span=0.75,
							 show.raw =FALSE,
							 overlay.hist=TRUE,
							 col.hist=adjustcolor("grey",0.5),
							 pdf.prefix,...){
	#assertions...
	plot.dat <- get.plot.dat  (bp.obj, raw.or.norm, selected.genes, selected.cell.types, bin.by, n.bins, span)
	
	color.num <- length(selected.cell.types)
	if(color.num<=8) my.palette <- brewer.pal(color.num,"Dark2")
	if(color.num>8 && color.num<=12) my.palette <- brewer.pal(color.num,"Paired")
	if(color.num > 12)  my.palette <-  colorRampPalette(brewer.pal(12, "Paired"))(color.num)
	
	pdf(paste(pdf.prefix,".pdf",sep=""),useDingbats=F,pointsize=8)

	if(!overlay.hist){
		for (ct.idx in 1:length(selected.cell.types)){	
			par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
			par(mar=c(5.1, 5.1, 4.1, 2.1))
			hist(plot.dat[[1]][[ct.idx]]$cord.vec,
				 main= selected.cell.types[ct.idx], 
				 xlim= c(0,1), 
				 breaks=50, 
				 col= adjustcolor(my.palette[ct.idx],0.5),
				 xlab="SpaceFold cordinate",
				 ylab="Cell Type Frequency")
		}
		axis.pos <- "left"
	}
	
	#making plots
	for (gene.idx in 1:length(selected.genes)){
		#get the range
		ylim.gene <- get.ylim (plot.dat, gene.idx, show.raw)
		for (ct.idx in 1:length(selected.cell.types)){
			par(mar=c(5.1, 5.1, 4.1, 5.1))
			par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
				
			if(overlay.hist){
				hist(plot.dat[[1]][[ct.idx]]$cord.vec,
				 main= "", 
				 xlim= c(0,1), 
				 breaks=50, 
				 col= col.hist,
				 xlab="SpaceFold cordinate",
				 ylab="Cell Type Frequency")
				 par(new=TRUE) 
			}
			
			par(mar=c(5.1, 5.1, 4.1, 5.1))
			par(cex.axis=2, cex.lab=2, cex.main=2, cex.sub=2)
			plot.exp.curv (plot.dat[[gene.idx]][[ct.idx]],
						  ct.name= selected.cell.types[ct.idx], 
						  gene.name= selected.genes[gene.idx],
						  show.raw= show.raw,
						  ylim= ylim.gene,
						  color= my.palette[ct.idx],
						  xlab="SpaceFold cordinate",
						  ylab=paste(raw.or.norm, "expression value"),
						  axis.pos= axis.pos)
		}
	}
	
	dev.off()
	NULL
}

#test
# # plot.cartography.nopdfjam (bp.obj=WTLI.bp.merged,
							 # raw.or.norm="norm",
							 # selected.genes= selected.genes,
							 # selected.cell.types=names(grouping.list),
							 # bin.by="equal.size",
							 # show.raw =TRUE,
							 # overlay.hist=FALSE,
							 # pdf.prefix="test")



#visualize gene expression level along the projected axis
# If pdfjam is installed use it. Otherwise,  make individual plot and output in a single pdf.
#' @param bp.obj a BayesPrism output object
#' @param raw.or.norm plotting raw value or value normalized by size factor
#' @param selected.genes a character vector specifiying the gene name / symbols to be plotted, which is automatically determined by matching columns of  bp.obj$para$feature 
#' @param selected.cell.types a character vector specifiying the cell types / cell type groups to be plotted, which is automatically determined by matching colnames of bp.obj$res$first.gibbs.res$theta.merged and bp.obj$res$res.regrouped$theta0
#' @param bin.by method to bin the data, deviding them into each bin of each sample size (equal.size), 
		# or by the same distance over spacefold axis (equal.step)
#' @param n.bins number of bins, default=20
#' @param span parameter used by loess curve fitting, default=0.75
#' @param show.raw whether plot the unbinned Znkg/Znkg.norm
#' @param overlay.hist whether plot the histogram of the cordiantes containing each cell type separately or on top of the expression curve.
#' @param col.hist color of the histogram
#' @param pdf.prefix the prefix of the pdf output
plot.cartography <- function(bp.obj,
						raw.or.norm=c("raw","norm"),
						selected.genes,
						selected.cell.types,
						bin.by=c("equal.size","equal.step"),
						n.bins=20,
						span=0.75,
						show.raw =FALSE,
						overlay.hist=TRUE,
						col.hist=adjustcolor("grey",0.5),
						pdf.prefix,
						...){
	#check if pdfjam is installed
	check.pdfjam <- tryCatch( system("pdfjam -V"), error=function(e) {return (-1)},
											warning=function(w)return(-1))
	if(check.pdfjam != -1){
		plot.cartography.pdfjam (bp.obj,
						raw.or.norm,
						selected.genes,
						selected.cell.types,
						bin.by,
						n.bins,
						span,
						show.raw,
						overlay.hist,
						col.hist,
						pdf.prefix,
						...)
	}
	else{
		plot.cartography.nopdfjam (bp.obj,
						raw.or.norm,
						selected.genes,
						selected.cell.types,
						bin.by,
						n.bins,
						span,
						show.raw,
						overlay.hist,
						col.hist,
						pdf.prefix,
						...)
	}
}


# # #test
# plot.cartography (bp.obj=WTLI.bp.merged,
							 # raw.or.norm="raw",
							 # selected.genes= selected.genes,
							 # selected.cell.types=names(grouping.list),
							 # bin.by="equal.size",
							 # show.raw =TRUE,
							 # overlay.hist=TRUE,
							 # pdf.prefix="test.raw")




