genome = "hg19"
tss.flanking.length = 1500

promoter.bed = paste("TxDb_",genome,"_promoter_",tss.flanking.length,"bp.gtf",sep="")

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
 
gene.keys = keys(txdb, keytype="GENEID")
txcols = c("TXCHROM", "TXSTART","TXEND","GENEID","TXID","TXNAME","TXSTRAND")
txtable = select(txdb, keys=gene.keys, columns=txcols, keytype="GENEID")
print(dim(txtable))

gene.symbols = keys(orgdb, keytype="SYMBOL")
genecols = c("SYMBOL", "ENTREZID","GENENAME")
genetable = select(orgdb, keys=gene.symbols, columns=genecols, keytype="SYMBOL")
print(dim(genetable))

promoter.chr = c()
promoter.start = c()
promoter.stop = c()
promoter.name = c()
promoter.strand = c()
promoter.transcript = c()
promoter.gene = c()
for (i in 1:nrow(txtable)){
	promoter.chr[i]=as.character(txtable$TXCHROM[i])
	promoter.strand[i] = as.character(txtable$TXSTRAND[i])
	tx.name = as.character(txtable$TXNAME[i])
	tx.geneID = as.character(txtable$GENEID[i])
	if(promoter.strand[i] == "+"){
		TSS = txtable$TXSTART[i]
		promoter.start[i]=TSS-tss.flanking.length
		promoter.stop[i]=TSS+tss.flanking.length
	}else if(promoter.strand[i] == "-"){
		TSS = txtable$TXEND[i]
		promoter.start[i]=TSS-tss.flanking.length
		promoter.stop[i]=TSS+tss.flanking.length
	}else {
		print(paste(promoter.strand[i]," is not a valid strand value",sep=""))
	}
	region = paste(promoter.chr[i],":",promoter.start[i],"-",promoter.stop[i],":",promoter.strand[i],sep="")
	gene.name = as.character(genetable$SYMBOL[match(tx.geneID,genetable$ENTREZID)])
	promoter.name[i]=paste(gene.name,tx.name,region,sep="_")
	promoter.transcript[i]=tx.name
	promoter.gene[i]=gene.name
}#end for (i in 1:nrow(txtable))

attribute = paste("gene_id \"",promoter.name,"\"; gene_name \"",promoter.gene,"\"; transcript_id \"",promoter.transcript,"\"",sep="")
bed.table = data.frame(promoter.chr, source=rep(paste("Bioconductor_",genome,sep=""),length(promoter.name)),
						feature=rep("promoter",length(promoter.name)),
						promoter.start, promoter.stop,
						score = rep(".",length(promoter.name)),promoter.strand,
						frame = rep(".",length(promoter.name)), attribute)
write.table(bed.table,promoter.bed, sep="\t", row.names=F, col.names=F, quote=F)