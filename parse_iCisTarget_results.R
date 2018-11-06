dbaFile = "merged_peaks_[compID]_[criteria].txt"
inputBed = "../Results/i-cisTarget_motifs/input_files/merged_peaks_2kb_[compID]_[direction].bed"
outputResult = "../Results/i-cisTarget_motifs/merged_peaks_2kb_[compID]_[direction]_peaks/statistics.tbl"
outputMapping = "../Results/i-cisTarget_motifs/merged_peaks_2kb_[compID]_[direction]_peaks/input_mapped_to_icistarget_regions.bed"
outputFile = "../Results/i-cisTarget_motifs/merged_peaks_2kb_[compID]_[direction]_peaks.txt"

dba.table = read.delim(dbaFile, head=T, sep="\t")
dbaID = paste(dba.table$chr,":",dba.table$start,"-",dba.table$stop,sep="")
icis.input = read.table(inputBed, head=F, sep="\t")
inputID = paste(icis.input$V1,":",icis.input$V2,"-",icis.input$V3,sep="")
icis.result =  read.table(outputResult, head=F, sep="\t")
icis.map =  read.table(outputMapping, head=F, sep="\t")

print(dim(icis.map))
icis.map=icis.map[icis.map$V8 != "",]
print(dim(icis.map))
map.dba.regionID = paste(icis.map$V1,":",icis.map$V2,"-",icis.map$V3,sep="")
map.regID=paste(icis.map$V5,":",icis.map$V6,"-",icis.map$V7,sep="")

map.symbol = c()
for(i in 1:nrow(icis.map)){
	temp.regionID = map.dba.regionID[i]
	map.symbol[i]=as.character(dba.table$symbol[match(temp.regionID, dbaID)])	
}#end for(i in 1:nrow(icis.map))

icis.map = data.frame(regID=icis.map$V8, map.regID, map.dba.regionID, map.symbol)

icis.regionIDs = c()
dba.regionIDs = c()
dba.symbols = c()

for (i in 1:nrow(icis.result)){
	targets = unlist(strsplit(as.character(icis.result[i,9]),split=";"))
	temp.icisIDs =icis.map$map.regID[match(targets,icis.map$regID)]
	temp.regionIDs = icis.map$map.dba.regionID[match(targets,icis.map$regID)]
	temp.symbols = icis.map$map.symbol[match(targets,icis.map$regID)]
	icis.regionIDs[i]=paste(temp.icisIDs, collapse=";")
	dba.regionIDs[i]=paste(temp.regionIDs, collapse=";")
	dba.symbols[i]=paste(temp.symbols, collapse=";")
}#end for (i in 1:nrow(icis.result))


output.table = data.frame(icis.result[,c(2:9,11)],
							icis.regionIDs, dba.regionIDs, dba.symbols)
names(output.table)=c("Rank","FeatureID","FeatureDescription","FeatureAnnotations","FeatureDatabase","AUC","NES","icis.TargetID","TopTargetRanks","icis.targetPos","dba.regionIDs","dba.symbols")
write.table(output.table, outputFile, row.names=F, sep="\t", quote=F)

