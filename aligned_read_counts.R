library(GenomicAlignments)

param.table = read.table("parameters.txt", header=T, sep="\t")
alignment.folder = as.character(param.table$Value[param.table$Parameter == "Alignment_Folder"])
length.file = as.character(param.table$Value[param.table$Parameter == "chr_length_file"])
aligned.stats.file = as.character(param.table$Value[param.table$Parameter == "aligned_stats_file"])
nodup.bam = as.character(param.table$Value[param.table$Parameter == "Remove_Duplicates"])
omit.chr = unlist(strsplit(as.character(param.table$Value[param.table$Parameter == "omitChr"]), split=","))

length.table = read.table(length.file, header=T, sep="\t")
chr_length = as.numeric(length.table$Length)
names(chr_length) = as.character(length.table$Chr)

#remove non-canonical chromosomes
nonCanonical = grep("_", names(chr_length))
 if (length(nonCanonical) > 0) {
     chr_length = chr_length[-nonCanonical]
}
#if desired, remove other chromosomes
 if ((length(omit.chr) > 0)&(omit.chr != "none")&(omit.chr != "no")) {
     chr_length = chr_length[-match(omit.chr,names(chr_length))]
}

chromosomes = names(chr_length)
print(chromosomes)

sampleIDs = c()
aligned.reads = c()

bam.files = list.files(alignment.folder, pattern=".bam$")
if(length(grep("sort.bam",bam.files)) > 0){
	bam.files = bam.files[-grep("sort.bam",bam.files)]
}

if(nodup.bam == "no"){
	suffix = ".bam$"
}else if(nodup.bam == "yes"){
	suffix = ".nodup.bam$"
}else{
	print("'Remove_Duplicates' must be 'yes' or 'no'")
	stop()
}

sampleIDs = sub(suffix,"",bam.files)

for(i in 1:length(bam.files)){
	sampleIDs[i]=gsub(suffix,"",bam.files[i])
	inputfile = paste(alignment.folder, bam.files[i], sep="/")
	print(sampleIDs[i])
	print(inputfile)
    
	total_reads = list()
    
	for(chr in chromosomes){
    		print(chr)
        	data = readGAlignments(file = inputfile, use.names = TRUE, 
            					param = ScanBamParam(which = GRanges(chr, IRanges(1, chr_length[chr]))))
		total_reads[[chr]] = unique(as.character(names(data)))
		rm(data)
	}#end  for(chr in chromosomes)
	
	aligned.reads[i] = length(unique(unlist(total_reads)))
	rm(total_reads)
	print(aligned.reads)
}#end for(i in 1:length(bam.files))

stat.table = data.frame(Sample = sampleIDs, aligned.reads = aligned.reads)
write.table(stat.table, aligned.stats.file, row.names=F, sep="\t", quote=F)
print(warnings())
