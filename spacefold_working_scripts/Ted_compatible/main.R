#library required:

library(dplyr)  
library(mixtools) 
library(mclust)
library(phateR)
library(RColorBrewer)
library(beeswarm)

#load source files
source.dir <- "/lila/data/peer/tinyi/RU_data/visium/spacefold/"
source (paste(source.dir,"bp.merge.function.R",sep=""))
source (paste(source.dir,"bp.utility.functions.R",sep=""))
source (paste(source.dir,"project.exp.functions.R",sep=""))
source (paste(source.dir,"run.phate.functions.R",sep=""))


#small intestine
#load BayesPrism output
load("/lila/data/peer/tinyi/RU_data/visium/dat/bp_dat/from_reprocessed_anndata/SI/merge_finer_states_0930/WTSI1.bp.rdata")
load("/lila/data/peer/tinyi/RU_data/visium/dat/bp_dat/from_reprocessed_anndata/SI/merge_finer_states_0930/WTSI2.bp.rdata")

#load orignal spaceranger ouptut
load("/lila/data/peer/tinyi/RU_data/visium/dat/visium_dat/SI/visium.SI.White.rdata")
#this file should contain spot annotation file

#spot annotation dataframe should look like this (see tutorial from 10x for details of generating this dataframe from spaceranger output):
#only the first two column is important
head(bcs_merge)
             # barcode sample tissue row col imagerow  imagecol Cluster height
# 1 AAACACCAATAACTGC-1  WTSI2      1  59  19 1437.812  420.0000       1   2000
# 2 AAACAGCTTTCAGAAG-1  WTSI1      1  43   9 1055.000  312.1875       3   2000
# 3 AAACAGGGTCTATATT-1  WTSI1      1  47  13 1140.312  360.6250       1   2000
# 4 AAACAGGGTCTATATT-1  WTSI2      1  47  13 1182.812  348.4375       2   2000
# 5 AAACATTTCCCGGATT-1  WTSI1      1  61  97 1445.312 1385.0000       7   2000
# 6 AAACATTTCCCGGATT-1  WTSI2      1  61  97 1486.562 1370.9375       2   2000
  # width sum_umi sum_gene
# 1  2000   11148     3411
# 2  2000   12317     3864
# 3  2000   44117     5430
# 4  2000   51281     6652
# 5  2000   11016     3265
# 6  2000   12937     4093


#add meta (spot annotation) 
WTSI1.bp <- add.meta (bp.obj= WTSI1.bp, meta = bcs_merge, selected.sample="WTSI1")
WTSI2.bp <- add.meta (bp.obj= WTSI2.bp, meta = bcs_merge, selected.sample="WTSI2")


#merge two output
bp.list <- list(WTSI1.bp, WTSI2.bp)
names(bp.list) <- c("WTSI1","WTSI2")
WTSI.bp <- merge.bp(bp.list)


#and feature (gene annotation)
gene.tab <- read.table("/lila/data/peer/tinyi/RU_data/visium/dat/visium_dat/SI/WTSI2_all_reads_D1_outs/raw_feature_bc_matrix/features.tsv.gz",sep="\t",check.names=F,header=F)[,1:2]
colnames(gene.tab) <- c("id",'symbol')
rownames(gene.tab) <- gene.tab[,"id"]

head(gene.tab)
                                   # id  symbol
# ENSMUSG00000051951 ENSMUSG00000051951    Xkr4
# ENSMUSG00000089699 ENSMUSG00000089699  Gm1992
# ENSMUSG00000102331 ENSMUSG00000102331 Gm19938
# ENSMUSG00000102343 ENSMUSG00000102343 Gm37381
# ENSMUSG00000025900 ENSMUSG00000025900     Rp1
# ENSMUSG00000025902 ENSMUSG00000025902   Sox17

WTSI.bp <- add.feature (WTSI.bp, gene.tab)

#update cell type names (optional)
old.name.vec <- c("Lgr5-_progenitor",
					"enterocyte_proximal_like",
										  "enterocyte_distal_like",
										  "enterocyte_distal_like_low_immature_4",
										  "Jchain+_plasma_paneth",
										  "T_cell_paneth",
										  "CD4/8+_Mzb1+_Trac-")
new.name.vec <- c("TA",
					"enterocyte_bottom_zone",
						   				  "enterocyte_mid_zone",
						   				  "enterocyte_top_zone",
						   				  "plasma",
						   				  "T_cell",
						   				  "plasmacytoid_dendritic")

WTSI.bp <- update.name(WTSI.bp, old.name.vec, new.name.vec)

#compute Znk (fraction of reads in each cell type)
WTSI.bp <- compute.Znk(WTSI.bp)

#normalize by Znk
WTSI.bp <- norm.by.sf(WTSI.bp)

#determine background level, using final theta
WTSI.bp <- compute.background.level (WTSI.bp, which.theta='final')

#run phate
# specify anchorCellTypes to select cell types (if you need to held out a group of cell types) 
WTSI.bp <- phate.bp (WTSI.bp, if.invert=TRUE)

#plot cell type distribution along the spacefold axis
plot.beeswarm(WTSI.bp, pdf.prefix="WTSI")



#regroup the data
#define group of cell types of interest
group.list <- list("Lgr5+_stem"="Lgr5+_stem",
					"Lgr5+_progenitor"="Lgr5+_progenitor",
					"TA"="TA",
					paneth="paneth",
					enterocytes=c("enterocyte_bottom_zone","enterocyte_mid_zone","enterocyte_top_zone"),
					goblet=c("goblet_3","goblet_1","goblet_paneth_2", "goblet_4"))

WTSI.bp <- regroup.bp(WTSI.bp, "type", group.list)


#plot cartography
selected.genes="Sox9"
plot.cartography (bp.obj= WTSI.bp,
							 raw.or.norm="raw",
							 selected.genes= selected.genes,
							 selected.cell.types=names(group.list),
							 bin.by="equal.size",
							 show.raw =TRUE,
							 overlay.hist=TRUE,
							 pdf.prefix="WTSI.Sox9.raw")

plot.cartography (bp.obj= WTSI.bp,
							 raw.or.norm="norm",
							 selected.genes= selected.genes,
							 selected.cell.types=names(group.list),
							 bin.by="equal.size",
							 show.raw =TRUE,
							 overlay.hist=TRUE,
							 pdf.prefix="WTSI.Sox9.norm.withRaw")


plot.cartography (bp.obj= WTSI.bp,
							 raw.or.norm="norm",
							 selected.genes= selected.genes,
							 selected.cell.types=names(group.list),
							 bin.by="equal.size",
							 show.raw =FALSE,
							 overlay.hist=TRUE,
							 pdf.prefix="WTSI.Sox9.norm.woRaw")

plot.cartography (bp.obj= WTSI.bp,
							 raw.or.norm="norm",
							 selected.genes= selected.genes,
							 selected.cell.types=names(group.list),
							 bin.by="equal.size",
							 show.raw =FALSE,
							 overlay.hist=TRUE,
							 n.bins=50,
							 pdf.prefix="WTSI.Sox9.norm.woRaw.bin50")

#
group.list <- list(paneth="paneth",
					"non-paneth_epithelial"=c("Lgr5+_progenitor","Lgr5+_stem", "TA",
												"goblet_3","goblet_1","goblet_paneth_2", "goblet_4",
												"tuft","enterendocrine",
												"enterocyte_bottom_zone","enterocyte_mid_zone", "enterocyte_top_zone"),
					"Lgr5+_epithelial"=c("Lgr5+_progenitor","Lgr5+_stem"),
					"Lgr5-_epithelial"=c("paneth","TA","goblet_3","goblet_1","goblet_paneth_2", "goblet_4",
										"tuft","enterendocrine",
										"enterocyte_bottom_zone","enterocyte_mid_zone", "enterocyte_top_zone"))
WTSI.bp <- regroup.bp(WTSI.bp, "type", group.list)

selected.genes=c("Lgr5","Defa5")
plot.cartography (bp.obj= WTSI.bp,
							 raw.or.norm="raw",
							 selected.genes= selected.genes,
							 selected.cell.types=names(group.list),
							 bin.by="equal.size",
							 show.raw =FALSE,
							 overlay.hist=TRUE,
							 pdf.prefix="WTSI.Lgr5_Defa5.raw")

secreted_factors=c('Ntn1','Il33','Wnt2','Ccl21a','Rspo3','Reln')
plot.cartography (bp.obj= WTSI.bp,
							 raw.or.norm="norm",
							 selected.genes= secreted_factors,
							 selected.cell.types=c('lymphatic','endothelial'),
							 bin.by="equal.size",
							 show.raw =FALSE,
							 overlay.hist=TRUE,
							 pdf.prefix="WTSI.secreted_factors.norm")

# if you need to split the object by visium sample, e.g. to make spatial plots, you may use the split function:
WTSI1.bp.split <- subset.bp(WTSI.bp, col.name="sample","WTSI1")
WTSI2.bp.split <- subset.bp(WTSI.bp, col.name="sample","WTSI2")




#process large intestine dataset

load("/lila/data/peer/tinyi/RU_data/visium/dat/bp_dat/from_reprocessed_anndata/LI/merge_finer_states_0930/WTLI1.bp.rdata")
load("/lila/data/peer/tinyi/RU_data/visium/dat/bp_dat/from_reprocessed_anndata/LI/merge_finer_states_0930/WTLI2.bp.rdata")

load("/lila/data/peer/tinyi/RU_data/visium/dat/visium_dat/LI/visium.LI.WhiteSketch.rdata")

WTLI1.bp <- add.meta (bp.obj=WTLI1.bp, meta = bcs_merge, selected.sample="WTLI1")
WTLI2.bp <- add.meta (bp.obj=WTLI2.bp, meta = bcs_merge, selected.sample="WTLI2")

#large intestine dataset has labels on different regions.
label.df <- read.csv("/data/peer/tinyi/RU_data/visium/dat/visium_dat/region_labeled/WTLI1_region.csv")
WTLI1.bp <- add.region.labels(WTLI1.bp, label.df)

label.df <- read.csv("/data/peer/tinyi/RU_data/visium/dat/visium_dat/region_labeled/WTLI2_region.csv")
WTLI2.bp <- add.region.labels(WTLI2.bp, label.df)

head(label.df)
             # Barcode region
# 1 AAACAAGTATCTCCCA-1 bottom
# 2 AAACAGAGCGACTCCT-1    top
# 3 AAACAGCTTTCAGAAG-1 bottom
# 4 AAACAGTGTTCCTGGG-1 bottom
# 5 AAACATTTCCCGGATT-1 bottom
# 6 AAACCACTACACAGAT-1    top

bp.list <- list(WTLI1.bp, WTLI2.bp)
names(bp.list) <- c("WTLI1","WTLI2")

WTLI.bp <-  merge.bp(bp.list)

gene.tab <- read.table("/lila/data/peer/tinyi/RU_data/visium/dat/visium_dat/SI/WTSI2_all_reads_D1_outs/raw_feature_bc_matrix/features.tsv.gz",sep="\t",check.names=F,header=F)[,1:2]
colnames(gene.tab) <- c("id",'symbol')
rownames(gene.tab) <- gene.tab[,"id"]

WTLI.bp <- add.feature (WTLI.bp, gene.tab)


#update cell type names (optional)
old.name.vec <- c("goblet_2+8")
new.name.vec <- c("goblet_2_8")
WTLI.bp <- update.name(WTLI.bp, old.name.vec, new.name.vec)

#compute Znk (fraction of reads in each cell type)
WTLI.bp <- compute.Znk(WTLI.bp)

#normalize by Znk
WTLI.bp <- norm.by.sf(WTLI.bp)

#determine background level, using final theta
WTLI.bp <- compute.background.level (WTLI.bp, which.theta='final')

#split by region, and then run PHATE


#run phate
# specify anchorCellTypes to select cell types (if you need to held out a group of cell types) 
WTLI.top.bp <- subset.bp(WTLI.bp, "region","top")
WTLI.bottom.bp <- subset.bp(WTLI.bp, "region","bottom")


WTLI.top.bp <- phate.bp (WTLI.top.bp, if.invert=FALSE)
WTLI.bottom.bp <- phate.bp (WTLI.bottom.bp, if.invert=FALSE)

bp.list <- list(WTLI.top.bp, WTLI.bottom.bp)
names(bp.list) <- c("top","bottom")

WTLI.merged.bp <- merge.bp (bp.list)


#plot cell type distribution along the spacefold axis
WTLI.merged.bp$para$spacefold <- WTLI.top.bp$para$spacefold 
plot.beeswarm(WTLI.merged.bp, pdf.prefix="WTLI")




