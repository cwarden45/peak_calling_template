param.file = "parameters.txt"
peak.length.dist = "peak_length_dist.png"
peak.length.stats = "peak_length_dist.txt"

max.plot = 1000

param.table = read.table(param.file, header=T, sep="\t")
genome = as.character(param.table$Value[param.table$Parameter == "genome"])
alignment.folder = as.character(param.table$Value[param.table$Parameter == "Alignment_Folder"])
merged.peak.file = as.character(param.table$Value[param.table$Parameter == "mergedPeakGTF"])
peakType = as.character(param.table$Value[param.table$Parameter == "peakType"])
annoType = as.character(param.table$Value[param.table$Parameter == "gtfID"])
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
tss.GTF = as.character(param.table$Value[param.table$Parameter == "promoterGTF_PC"])
merged.GTF = as.character(param.table$Value[param.table$Parameter == "mergedPeakGTF"])

library("GenomicRanges")

if(peakType == "broad"){
	bamSuffix = ".nodup.bam$"
	peakSuffix = "_peaks.broadPeak"
}else if(peakType == "narrow"){
	print("Add narrow peak code")
	stop()
}else{
	print("'peakType' must be 'broad' or 'narrow'")
	stop()	
}

files = list.files(alignment.folder, pattern=bamSuffix)
sampleIDs = gsub(bamSuffix,"",files)

bedFiles = paste(alignment.folder,"/", sampleIDs,"/macs_",sampleIDs,peakSuffix,sep="")
sample.description.table = read.table(sample.description.file, sep="\t", header=T)
longID = sample.description.table$sampleID
sample.label = sample.description.table$userID

matchedIDs = sample.label[match(sampleIDs, longID)]
if(length(matchedIDs) != length(matchedIDs[!is.na(matchedIDs)])){
	print("There are mismatched samples")
	names(matchedIDs) = sampleIDs
	print(matchedIDS)
	stop()
}else{
	bedFiles = bedFiles[match(longID,sampleIDs)]
}

num.peaks = c()

png(peak.length.dist)
for (i in 1:length(bedFiles)){
	print(paste(longID[i]," : ",bedFiles[i],sep=""))
	
	input.table = read.table(bedFiles[i], head=F, sep="\t")
	print(dim(input.table))
	num.peaks[i]=nrow(input.table)
	peak.length = input.table$V3-input.table$V2
	print(max(peak.length))
	peak.length[peak.length > max.plot]=max.plot
	print(length(peak.length[peak.length == max.plot]))
	
			if(i == 1)
				{
					den = density(peak.length, na.rm=T,from=0, to=max.plot)
					plot(den$x, den$y, type="l", xlab = "Peak Length", ylab = "Density",
							xlim=c(0,max.plot), ylim=c(0,0.01), col="gray")
					legend("top",legend=c("sample","merged"),col=c("gray","black"),
							lwd=2, xpd=T, inset =-0.1, ncol=2)
				}#end if(i == 1)
			else
				{
					den = density(peak.length, na.rm=T,from=0, to=max.plot)

					lines(den$x, den$y, type = "l", col="gray")
				}#end else

	gr = reduce(GRanges(Rle(input.table$V1),
    		IRanges(start=input.table$V2, end=input.table$V3)))
			
	if(i == 1){
		peakGr = gr
	}else{
		peakGr = union(peakGr, gr)
	}
}#end for (i in 1:length(bedFiles))
peakGr = reduce(peakGr)

print("Create merged peak table")
merged.peaks = data.frame(peakGr)
peak.length = merged.peaks$width
print(max(peak.length))
	peak.length[peak.length > max.plot]=max.plot
print(length(peak.length[peak.length == max.plot]))
lines(den$x, den$y, type = "l", col="black")
dev.off()

stat.table = data.frame(SampleID=longID, userID=sample.label, bed.file=bedFiles, num.peaks)
write.table(stat.table, peak.length.stats, quote=F, sep="\t", row.names=F)

promoter.table = read.table(tss.GTF, head=F, sep="\t")
print(dim(promoter.table))
if(length(grep("_",promoter.table$V1)) > 0){
	promoter.table = promoter.table[-grep("_",promoter.table$V1),]
	promoter.table$V1 = as.factor(as.character(promoter.table$V1))
	print(dim(promoter.table))
}

extract.gene = function(char){
	char.info = unlist(strsplit(char,split=";"))
	anno.info = as.character(char.info[grep("gene_name", char.info)])
	anno = gsub("\"","",anno.info)
	anno = gsub("gene_name","",anno)
	anno = gsub(" ","",anno)
	return(anno)
}#end def extract.gene

gene = unlist(sapply(as.character(promoter.table$V9), extract.gene))
refGR = GRanges(Rle(promoter.table$V1),
				IRanges(start=promoter.table$V4, end=promoter.table$V5),
				Names=gene, 
				Rle(strand(promoter.table$V7)))

testGR = GRanges(Rle(merged.peaks$seqnames),
					IRanges(start=merged.peaks$start, end=merged.peaks$end))
hits = findOverlaps(refGR, testGR)
overlaps = data.frame(pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)]))
print(head(overlaps))

regionID = paste(merged.peaks$seqnames,":",merged.peaks$start,"-",merged.peaks$end,sep="")
overlapsID = paste(overlaps$seqnames,":",overlaps$start,"-",overlaps$end,sep="")

combine_genes = function(gene.arr){
	gene.arr = unique(gene.arr)
	gene.arr = gene.arr[!is.na(gene.arr)]
	
	if(length(gene.arr) == 0){
		return("other")
	}else{
		return(paste(gene.arr,collapse="_"))
	}
}#end def combine_genes

region.gene = tapply(overlaps$Names, overlapsID, combine_genes)
region.gene=region.gene[match(regionID, names(region.gene))]
region.gene[!is.na(region.gene)]=paste("TSS_",region.gene[!is.na(region.gene)],sep="")
region.gene[is.na(region.gene)]="other"

tabularID = paste(region.gene,regionID,sep="_")

attribute = paste("gene_id \"",regionID,"\"; gene_name \"",region.gene,"\"; transcript_id \"",tabularID,"\"",sep="")
gtf.table = data.frame(chr=merged.peaks$seqnames, source=rep("merged_MACS2_peaks",nrow(merged.peaks)),
						feature=rep("peak",nrow(merged.peaks)),
						start = merged.peaks$start+1,
						end = merged.peaks$end+1,
						score = rep(".",nrow(merged.peaks)),
						strand = rep(".",nrow(merged.peaks)),
						frame = rep(".",nrow(merged.peaks)), attribute)
print(dim(gtf.table))
write.table(gtf.table,merged.GTF, sep="\t", row.names=F, col.names=F, quote=F)