param.table = read.table("parameters.txt", header=T, sep="\t")
genome = as.character(param.table$Value[param.table$Parameter == "genome"])
alignment.folder = as.character(param.table$Value[param.table$Parameter == "Alignment_Folder"])
merged.peak.file = as.character(param.table$Value[param.table$Parameter == "mergedPeakGTF"])
peakType = as.character(param.table$Value[param.table$Parameter == "peakType"])
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
tss.dist = as.numeric(as.character(param.table$Value[param.table$Parameter == "flankTSS"]))

library("ChIPseeker")

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
	names(bedFiles)=matchedIDs
}

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

#use larger window to test validity of smaller window
#promoter= getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
#tagMatrixList = lapply(bedFiles, getTagMatrix, windows=promoter)

#pdf("ChIPseeker_TSS_heatmap_Bowtie2.pdf", width=10, height=5)
#tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
#dev.off()

##double-comment because this function takes even longer than the above function
##pdf("ChIPseeker_TSS_line_coverage.pdf", width=10, height=5)
##plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")
##dev.off()

peakAnnoList = lapply(bedFiles, annotatePeak, TxDb=txdb,
                       tssRegion=c(-tss.dist, tss.dist), verbose=FALSE)

png("ChIPseeker_anno_type_barplot_Bowtie2.png", width=500, height=250)
plot.obj = plotAnnoBar(peakAnnoList)
print(plot.obj)
dev.off()

png("ChIPseeker_TSS_distance_barplot_Bowtie2.png", width=500, height=250)
plot.obj = plotDistToTSS(peakAnnoList, title="Distance of Peak Relative to TSS")
print(plot.obj)
dev.off()
