import sys
import re
import os
import subprocess

parameterFile = "parameters.txt"
statFile = "coverage_statistics.txt"

readsFolder = ""
alignmentFolder = ""

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]

	if param == "Reads_Folder":
		readsFolder = value

	if param == "Alignment_Folder":
		alignmentFolder = value
		
if (readsFolder == "") or (readsFolder == "[required]"):
	print "Need to enter a value for 'Reads_Folder'!"
	sys.exit()
	
if (alignmentFolder == "") or (alignmentFolder == "[required]"):
	print "Need to enter a value for 'Alignment_Folder'!"
	sys.exit()


statHandle = open(statFile,"w")
text = "SampleID\tSeqID\tuserID\tTotalReads\tPercent.Aligned\tPercent.Duplicate\tPercent.Proper.chrM\n"
statHandle.write(text)
	
fastqcFolder = readsFolder + "/QC"
fileResults = os.listdir(readsFolder)

for file in fileResults:
	result = re.search("(.*)_S\d+_L\d{3}_R1_001.fastq.gz$",file)
	
	if result:
		sample = result.group(1)
		print sample
		
		r2 = re.search("^(\d+)_.*",sample)
		seqID = r2.group(1)
		
		shortID = re.sub(seqID + "_","",sample)
		
		#get total reads from FastQC
		fastqcPrefix = re.sub(".fastq.gz","",file)
		fastQCtext = fastqcFolder + "/" + fastqcPrefix + "_fastqc/fastqc_data.txt"
		
		inHandle = open(fastQCtext)
		line = inHandle.readline()
		
		lineCount = 0
		
		while line:
			line = re.sub("\n","",line)
			line = re.sub("\r","",line)
			
			lineCount += 1
			
			if lineCount == 7:
				
				totalResult = re.search("Total Sequences\t(\d+)",line)
				if totalResult:
					totalReads = totalResult.group(1)
				else:
					print "Problem parsing FastQC file!\n"
					sys.exit()
			
			line = inHandle.readline()
		
		inHandle.close()
		
		#get alignment stats from samtools flagstat
		sampleSubfolder = alignmentFolder + "/" + sample
		flagstatFile = sampleSubfolder + "/alignment_stats.txt"

		inHandle = open(flagstatFile)
		line = inHandle.readline()
		
		lineCount = 0
		
		while line:
			line = re.sub("\n","",line)
			line = re.sub("\r","",line)
			
			lineCount += 1
			
			if lineCount == 3:
				
				alignedResult = re.search("\((\d+\.\d+%):",line)
				if alignedResult:
					alignedReads = alignedResult.group(1)
				else:
					print "Problem parsing samtools flagstat file!\n"
					sys.exit()
				
			line = inHandle.readline()
		
		inHandle.close()
		
		#get duplicate rate from Picard MarkDuplicates
		duplicateMetrics = sampleSubfolder + "/MarkDuplicates_metrics.txt"

		inHandle = open(duplicateMetrics)
		line = inHandle.readline()
		
		lineCount = 0
		
		while line:
			line = re.sub("\n","",line)
			line = re.sub("\r","",line)
			
			lineCount += 1
			
			if lineCount == 8:
				lineInfo = line.split("\t")
				
				percentDuplicate = 100*float(lineInfo[8])
				percentDuplicate =  '{0:.2f}'.format(percentDuplicate) + "%"
				
			line = inHandle.readline()
		
		inHandle.close()
		
		#get chrM rate from idxstats
		chrMetrics = sampleSubfolder + "/idxstats_no_dup.txt"

		inHandle = open(chrMetrics)
		line = inHandle.readline()
		
		lineCount = 0
		
		chrM_counts = 0
		total_adjusted_counts = 0
		
		while line:
			line = re.sub("\n","",line)
			line = re.sub("\r","",line)
			
			lineInfo = line.split("\t")
			
			chr = lineInfo[0]
			aligned = lineInfo[2]
			unaligned = lineInfo[3]
			
			adj_counts = int(aligned) - int(unaligned)
			if adj_counts < 0:
				adj_counts = 0
				
			if (chr == "chrM")or(chr == "MT"):
				chrM_counts = adj_counts
			
			if chr != "*":
				total_adjusted_counts += adj_counts
				
			line = inHandle.readline()
		
		inHandle.close()

		percentChrM = 100*float(chrM_counts)/float(total_adjusted_counts)
		percentChrM =  '{0:.2f}'.format(percentChrM) + "%"
		
		text = sample + "\t" + seqID + "\t" + shortID + "\t" + totalReads + "\t" + alignedReads + "\t" + percentDuplicate + "\t" + percentChrM + "\n"
		statHandle.write(text)