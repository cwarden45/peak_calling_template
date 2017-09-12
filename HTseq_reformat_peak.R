normalizeTotalExpression = function (geneExpr, totalReads) {
	return(geneExpr / totalReads)
}#end def normalizeTotalExpression

count.defined.values = function(arr, expr.cutoff)
{
	sig.figures = 1
	if (expr.cutoff > 0)
		sig.figures = 0
	expr.cutoff = round(expr.cutoff, digits=sig.figures)
	arr = round(arr, digits=sig.figures)
	return(length(arr[arr > expr.cutoff]))
}#end def count.values

trimmed.counts = function(counts, min.percent, max.percent)
{
	total.counts = sum(counts)
	counts = counts[order(counts)]
	min.index = min.percent * length(counts)
	max.index = max.percent * length(counts)
	counts = counts[min.index:max.index]
	trimmed.counts = sum(counts)
	trimmed.percent = round(100 * trimmed.counts/total.counts, digits=1)
	return(trimmed.percent)
}#end def count.values

param.table = read.table("parameters.txt", header=T, sep="\t")
genome=as.character(param.table$Value[param.table$Parameter == "genome"])
min.expression = as.numeric(as.character(param.table$Value[param.table$Parameter == "fpkm_binding_cutoff"]))
user.folder = as.character(param.table$Value[param.table$Parameter == "Result_Folder"])
sample.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
aligned.stat.file = as.character(param.table$Value[param.table$Parameter == "aligned_stats_file"])
counts.file = as.character(param.table$Value[param.table$Parameter == "counts_file_peak"])
rpkm.file = as.character(param.table$Value[param.table$Parameter == "fpkm_file_peak"])
cpm.file = as.character(param.table$Value[param.table$Parameter == "CPM_file_peak"])
aligned.type = as.character(param.table$Value[param.table$Parameter == "FPKM_norm"])

library(ChIPseeker)
library(GenomicRanges)

sample.table = read.table(sample.file, header=T, sep="\t")
sampleID = as.character(sample.table$sampleID)
sample.label = as.character(sample.table$userID)
dash.flag = grep("-",sample.label)
if(length(dash.flag) > 0){
	print(paste(paste(sample.label[dash.flag],collapse=",")," samples labels have dashes in their labels",sep=""))
}
num.flag = grep("^[0-9]",sample.label)
if(length(num.flag) > 0){
	print(paste(paste(sample.label[num.flag],collapse=",")," samples labels start with numbers",sep=""))
}

if((length(dash.flag) > 0)|(length(num.flag) > 0)){
	stop("Please make sure sample labels do not have dashes and do not start with numbers")
}

count.files = as.character(sample.table$HTseq.peak)

temp.file = count.files[[1]]
temp.table = read.table(temp.file, sep="\t", header=F)
regions = as.character(temp.table[[1]])

count.mat = matrix(nrow=nrow(temp.table), ncol=length(sampleID))
colnames(count.mat) = sample.label

matched.regions = c()

for (i in 1:length(sampleID)){
	temp.file = count.files[[i]]
	temp.table = read.table(temp.file, sep="\t", header=F)
	temp.regions = as.character(temp.table[[1]])
	
	if(length(matched.regions) != 0){
		matched.regions = temp.regions[match(matched.regions, temp.regions, nomatch=0)]
	} else if(!identical(temp.regions, regions)){
		print("Regions not in same order - most likely, different quantification method was used for different samples")
		userAns = readline(prompt="Do you wish to proceed with a subset of matched gene symbols? (y/n): ")
		userAns = tolower(substr(userAns, 1, 1))
		if (userAns != "y"){
			stop("Please re-run htseq-count with all of your samples")
		} else {
			matched.regions = temp.regions[match(regions, temp.regions, nomatch=0)]
		}#end else
	} else {
		count.mat[,i] = temp.table[[2]]
	}#end else
}#end for (i in 1:length(sampleID))

if (length(matched.regions) != 0){
	count.mat = matrix(nrow=length(matched.regions), ncol=length(sampleID))
	colnames(count.mat) = sample.label
	
	for (i in 1:length(sampleID)){
	temp.file = count.files[[i]]
	temp.table = read.table(temp.file, sep="\t", header=F)
	temp.regions = as.character(temp.table[[1]])
	temp.counts = temp.table[[2]]
	
	count.mat[,i] = temp.counts[match(matched.regions, temp.regions, nomatch=0)]
	}#end for (i in 1:length(sampleID))
	
	regions = matched.regions
}#end if (length(matched.regions) != 0)

irrelevant.counts = c("__no_feature", "__ambiguous", "__too_low_aQual","__not_aligned","__alignment_not_unique")
extra.stats = count.mat[match(irrelevant.counts, regions),]
count.mat = count.mat[-match(irrelevant.counts, regions),]
regions = regions[-match(irrelevant.counts, regions)]
rownames(count.mat) = regions

### start code specific to region annotation ###
split.region.info = function(char){
	region.info = unlist(strsplit(char,":"))
	region.chr = region.info[1]
	region.info = unlist(strsplit(region.info[2],"-"))
	region.start = region.info[1]
	region.stop = region.info[2]
	region.list = list(chr=region.chr, start=region.start, stop=region.stop)
	return(region.list)
}#end def split.region.info

grMat = data.frame(t(sapply(as.character(regions), split.region.info)))
gr.chr = as.character(grMat[,"chr"])
gr.start = as.numeric(grMat[,"start"])
gr.stop = as.numeric(grMat[,"stop"])
gr = GRanges(Rle(gr.chr),
			IRanges(start=gr.start, end=gr.stop))

if(genome == "hg38"){
	library(TxDb.Hsapiens.UCSC.hg38.knownGene)
	txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
	
	library(org.Hs.eg.db)
	orgdb = org.Hs.eg.db
} else if(genome == "hg19"){
	library(TxDb.Hsapiens.UCSC.hg19.knownGene)
	txdb = TxDb.Hsapiens.UCSC.hg19.knownGene

	library(org.Hs.eg.db)
	orgdb = org.Hs.eg.db
} else if(genome == "hg18"){
	library(TxDb.Hsapiens.UCSC.hg18.knownGene)
	txdb = TxDb.Hsapiens.UCSC.hg18.knownGene

	library(org.Hs.eg.db)
	orgdb = org.Hs.eg.db
} else if (genome == "mm8"){
	library(TxDb.Mmusculus.UCSC.mm8.knownGene)
	txdb = TxDb.Mmusculus.UCSC.mm8.knownGene

	library(org.Mm.eg.db)
	orgdb = org.Mm.eg.db
} else if (genome == "mm9"){
	library(TxDb.Mmusculus.UCSC.mm9.knownGene)
	txdb = TxDb.Mmusculus.UCSC.mm9.knownGene

	library(org.Mm.eg.db)
	orgdb = org.Mm.eg.db
} else if (genome == "mm10"){
	library(TxDb.Mmusculus.UCSC.mm10.knownGene)
	txdb = TxDb.Mmusculus.UCSC.mm10.knownGene

	library(org.Mm.eg.db)
	orgdb = org.Mm.eg.db
} else {
	stop("Need to add annotations for reference!")
}
			
region.info = data.frame(annotatePeak(gr, TxDb=txdb))
regionIndex = paste(region.info$seqnames,":",region.info$start,"-", region.info$end,sep="")
region.info = region.info[match(regions,regionIndex),]

geneID = keys(orgdb, keytype="ENTREZID")
genecols = c("SYMBOL", "ENTREZID","GENENAME")
genetable = select(orgdb, keys=geneID, columns=genecols, keytype="ENTREZID")

region.length.kb = (gr.stop-gr.start)/1000
region.info = data.frame(regionID=regions, grMat, region.length.kb,
						symbol=genetable$SYMBOL[match(region.info$geneId, genetable$ENTREZID)],
						gene = genetable$GENENAME[match(region.info$geneId, genetable$ENTREZID)],
						region.info[,6:ncol(region.info)])
						
#print(head(region.info))
region.info = data.frame(apply(region.info, 2, as.character))
#print(head(region.info))
						
text.strand = rep("-",nrow(region.info))
text.strand[region.info$geneStrand == 1]="+"
region.info$geneStrand = as.character(region.info$geneStrand)
region.info$geneStrand = text.strand
#print(head(region.info))
### end code specific to region annotation ###

annotated.counts = data.frame(region.info, count.mat)
write.table(annotated.counts, file = counts.file, sep="\t", row.names=F, quote=T)

result.file = paste(user.folder,counts.file,sep="/")
write.table(annotated.counts, file=result.file, row.names=F, quote=F, sep="\t")

intergenic.reads = extra.stats[irrelevant.counts == "__no_feature", ]
region.reads = apply(count.mat, 2, sum)
unique.reads = region.reads + intergenic.reads
percent.region.reads = round(100 * region.reads / unique.reads, digits=1)

if((aligned.type == "aligned")|(aligned.type =="quantified")){
	if(aligned.type == "aligned"){
		aligned.stat.table = read.table(aligned.stat.file, header=T, sep="\t")
		alignedID = as.character(aligned.stat.table$Sample)
		aligned.stat.table = aligned.stat.table[match(sampleID, alignedID),]
	
		aligned.reads = as.numeric(aligned.stat.table$aligned.reads)
	} else if(aligned.type =="quantified"){
		aligned.reads=apply(count.mat, 2, sum)
	}

	total.million.aligned.reads = aligned.reads / 1000000
	print(total.million.aligned.reads)

	rpk = matrix(ncol=ncol(count.mat), nrow=nrow(count.mat))
	for (i in 1:ncol(count.mat)){
		counts = as.numeric(count.mat[,i])
		temp.rpk = counts / as.numeric(region.length.kb)
		rpk[,i] = temp.rpk 
	}
	RPKM = round(log2(t(apply(rpk, 1, normalizeTotalExpression, totalReads = total.million.aligned.reads)) + min.expression), digits=2)
	colnames(RPKM) = sample.label
}else if(aligned.type == "TMM"){
	library(edgeR)
	
	y = DGEList(counts=count.mat, genes=regions)
	y = calcNormFactors(y, method="TMM")
	
	aligned.reads=y$samples$lib.size
	names(aligned.reads)=rownames(y$samples)
	aligned.reads = aligned.reads[match(sample.label, names(aligned.reads))]
	
	RPKM = round(log2(rpkm(y, gene.length = 1000 * as.numeric(region.length.kb))+min.expression), digits=2)
}else{
	stop("Print RPKM_norm must be either 'aligned', 'quantified', or 'TMM'")
}#end else

trimmed.percent = apply(count.mat, 2, trimmed.counts, min.percent=0.3, max.percent=0.95)

coverage.table = data.frame(Sample = sample.label,
							htseq.nofeature.reads =intergenic.reads, percent.unique.region.reads = paste(percent.region.reads,"%",sep=""),
							trimmed.percent=paste(trimmed.percent,"%",sep=""))
write.table(coverage.table, file="Bowtie2_peak_coverage_stats.txt", quote=F, row.names=F, sep="\t")

	
#tables have different file formats for downstream R analysis versus other applications that involve parsing text files
#	--> don't set the Result folder to the working directory, or you may skip regions during DBR analysis
annotated.rpkm = data.frame(region.info, RPKM)
write.table(annotated.rpkm, file = rpkm.file, sep="\t", row.names=F, quote=T)

result.file = paste(user.folder, rpkm.file, sep="/")
write.table(annotated.rpkm, file=result.file, row.names=F, quote=F, sep="\t")

#extra CPM tale

if((aligned.type == "aligned")|(aligned.type =="quantified")){
	if(aligned.type == "aligned"){
		aligned.stat.table = read.table(aligned.stat.file, header=T, sep="\t")
		alignedID = as.character(aligned.stat.table$Sample)
		aligned.stat.table = aligned.stat.table[match(sampleID, alignedID),]
	
		aligned.reads = as.numeric(aligned.stat.table$aligned.reads)
	} else if(aligned.type =="quantified"){
		aligned.reads=apply(count.mat, 2, sum)
	}

	total.million.aligned.reads = aligned.reads / 1000000
	print(total.million.aligned.reads)

	CPM = t(apply(rpk, 1, normalizeTotalExpression, totalReads = total.million.aligned.reads))
	colnames(CPM) = sample.label
}else if(aligned.type == "TMM"){
	library(edgeR)
	
	y = DGEList(counts=count.mat, genes=regions)
	y = calcNormFactors(y, method="TMM")
	
	aligned.reads=y$samples$lib.size
	names(aligned.reads)=rownames(y$samples)
	aligned.reads = aligned.reads[match(sample.label, names(aligned.reads))]
	
	CPM = cpm(y, region.length = 1000 * as.numeric(region.length.kb))
}else{
	stop("Print RPKM_norm must be either 'aligned', 'quantified', or 'TMM'")
}#end else

#tables have different file formats for downstream R analysis versus other applications that involve parsing text files
#	--> don't set the Result folder to the working directory, or you may skip regions during DBR analysis
annotated.cpm = data.frame(region.info, CPM)
write.table(annotated.cpm, file = cpm.file, sep="\t", row.names=F, quote=T)

result.file = paste(user.folder, cpm.file, sep="/")
write.table(annotated.cpm, file=result.file, row.names=F, quote=F, sep="\t")