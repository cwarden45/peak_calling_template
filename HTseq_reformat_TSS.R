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
counts.file = as.character(param.table$Value[param.table$Parameter == "counts_file_TSS"])
rpkm.file = as.character(param.table$Value[param.table$Parameter == "fpkm_file_TSS"])
cpm.file = as.character(param.table$Value[param.table$Parameter == "CPM_file_TSS"])
aligned.type = as.character(param.table$Value[param.table$Parameter == "FPKM_norm"])
tss.dist = as.numeric(as.character(param.table$Value[param.table$Parameter == "flankTSS"]))

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

count.files = as.character(sample.table$HTseq.tss)

temp.file = count.files[[1]]
temp.table = read.table(temp.file, sep="\t", header=F)
promoters = as.character(temp.table[[1]])

count.mat = matrix(nrow=nrow(temp.table), ncol=length(sampleID))
colnames(count.mat) = sample.label

matched.promoters = c()

for (i in 1:length(sampleID)){
	temp.file = count.files[[i]]
	temp.table = read.table(temp.file, sep="\t", header=F)
	temp.promoters = as.character(temp.table[[1]])
	
	if(length(matched.promoters) != 0){
		matched.promoters = temp.promoters[match(matched.promoters, temp.promoters, nomatch=0)]
	} else if(!identical(temp.promoters, promoters)){
		print("Promoters not in same order - most likely, different quantification method was used for different samples")
		userAns = readline(prompt="Do you wish to proceed with a subset of matched gene symbols? (y/n): ")
		userAns = tolower(substr(userAns, 1, 1))
		if (userAns != "y"){
			stop("Please re-run htseq-count with all of your samples")
		} else {
			matched.promoters = temp.promoters[match(promoters, temp.promoters, nomatch=0)]
		}#end else
	} else {
		count.mat[,i] = temp.table[[2]]
	}#end else
}#end for (i in 1:length(sampleID))

if (length(matched.promoters) != 0){
	count.mat = matrix(nrow=length(matched.promoters), ncol=length(sampleID))
	colnames(count.mat) = sample.label
	
	for (i in 1:length(sampleID)){
	temp.file = count.files[[i]]
	temp.table = read.table(temp.file, sep="\t", header=F)
	temp.promoters = as.character(temp.table[[1]])
	temp.counts = temp.table[[2]]
	
	count.mat[,i] = temp.counts[match(matched.promoters, temp.promoters, nomatch=0)]
	}#end for (i in 1:length(sampleID))
	
	promoters = matched.promoters
}#end if (length(matched.promoters) != 0)

irrelevant.counts = c("__no_feature", "__ambiguous", "__too_low_aQual","__not_aligned","__alignment_not_unique")
extra.stats = count.mat[match(irrelevant.counts, promoters),]
count.mat = count.mat[-match(irrelevant.counts, promoters),]
promoters = promoters[-match(irrelevant.counts, promoters)]
rownames(count.mat) = promoters

### start code specific to promoter annotation ###
split.promoter.info = function(char){
	char = gsub("_Promoter","",char)
	promoter.info = unlist(strsplit(char,"_"))
	promoter.gene = promoter.info[1]
	promoter.loc = promoter.info[2]
	
	promoter.info = unlist(strsplit(promoter.loc,":"))
	promoter.chr = promoter.info[1]
	promoter.strand = promoter.info[3]
	promoter.info = unlist(strsplit(promoter.info[2],"-"))
	promoter.start = promoter.info[1]
	promoter.stop = promoter.info[2]
	promoter.list = list(gene=promoter.gene,
						promoter.chr=promoter.chr,
						promoter.start=promoter.start, promoter.stop=promoter.stop,
						promoter.strand=promoter.strand)
	return(promoter.list)
}#end def split.promoter.info

promoter.info = data.frame(t(sapply(as.character(promoters), split.promoter.info)))
promoter.info = apply(promoter.info, 2, as.character)

promoter.length.kb = (2*tss.dist)/1000
promoter.info = data.frame(regionID=promoters, promoter.info, promoter.length.kb)
### end code specific to promoter annotation ###

annotated.counts = data.frame(promoter.info, count.mat)
write.table(annotated.counts, file = counts.file, sep="\t", row.names=F, quote=T)

result.file = paste(user.folder,counts.file,sep="/")
write.table(annotated.counts, file=result.file, row.names=F, quote=F, sep="\t")

intergenic.reads = extra.stats[irrelevant.counts == "__no_feature", ]
promoter.reads = apply(count.mat, 2, sum)
unique.reads = promoter.reads + intergenic.reads
percent.promoter.reads = round(100 * promoter.reads / unique.reads, digits=1)

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
		temp.rpk = counts / as.numeric(promoter.length.kb)
		rpk[,i] = temp.rpk 
	}
	RPKM = round(log2(t(apply(rpk, 1, normalizeTotalExpression, totalReads = total.million.aligned.reads)) + min.expression), digits=2)
	colnames(RPKM) = sample.label
}else if(aligned.type == "TMM"){
	library(edgeR)
	
	y = DGEList(counts=count.mat, genes=promoters)
	y = calcNormFactors(y, method="TMM")
	
	aligned.reads=y$samples$lib.size
	names(aligned.reads)=rownames(y$samples)
	aligned.reads = aligned.reads[match(sample.label, names(aligned.reads))]
	
	RPKM = round(log2(rpkm(y, gene.length = 1000 * as.numeric(promoter.length.kb))+min.expression), digits=2)
}else{
	stop("Print RPKM_norm must be either 'aligned', 'quantified', or 'TMM'")
}#end else

trimmed.percent = apply(count.mat, 2, trimmed.counts, min.percent=0.3, max.percent=0.95)

coverage.table = data.frame(Sample = sample.label,
							htseq.nofeature.reads =intergenic.reads, percent.unique.promoter.reads = paste(percent.promoter.reads,"%",sep=""),
							trimmed.percent=paste(trimmed.percent,"%",sep=""))
write.table(coverage.table, file="Bowtie2_TSS_coverage_stats.txt", quote=F, row.names=F, sep="\t")

	
#tables have different file formats for downstream R analysis versus other applications that involve parsing text files
#	--> don't set the Result folder to the working directory, or you may skip promoters during DBR analysis
annotated.rpkm = data.frame(promoter.info, RPKM)
print(dim(annotated.rpkm))
annotated.rpkm=annotated.rpkm[annotated.rpkm$promoter.length.kb > 0,]
print(dim(annotated.rpkm))
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

	CPM = t(apply(count.mat, 1, normalizeTotalExpression, totalReads = total.million.aligned.reads))
	colnames(CPM) = sample.label
}else if(aligned.type == "TMM"){
	library(edgeR)
	
	y = DGEList(counts=count.mat, genes=promoters)
	y = calcNormFactors(y, method="TMM")
	
	aligned.reads=y$samples$lib.size
	names(aligned.reads)=rownames(y$samples)
	aligned.reads = aligned.reads[match(sample.label, names(aligned.reads))]
	
	CPM = cpm(y, gene.length = 1000 * as.numeric(promoter.length.kb))
}else{
	stop("Print RPKM_norm must be either 'aligned', 'quantified', or 'TMM'")
}#end else

#tables have different file formats for downstream R analysis versus other applications that involve parsing text files
#	--> don't set the Result folder to the working directory, or you may skip promoters during DBR analysis
annotated.cpm = data.frame(promoter.info, CPM)
write.table(annotated.cpm, file = cpm.file, sep="\t", row.names=F, quote=T)

result.file = paste(user.folder, cpm.file, sep="/")
write.table(annotated.cpm, file=result.file, row.names=F, quote=F, sep="\t")