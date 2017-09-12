import sys
import re
import os

parameterFile = "parameters.txt"
finishedSamples = ()

alignmentFolder = ""
peakGTF = ""
tssGTF = ""
gtfID = ""
email = ""
peStatus = ""
peakType = ""

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]
	
	if param == "Alignment_Folder":
		alignmentFolder = value
		
	if param == "mergedPeakGTF":
		peakGTF = value

	if param == "promoterGTF":
		tssGTF = value		

	if param == "gtfID":
		gtfID = value
		
	if param == "Cluster_Email":
		email = value

	if param == "PE_Reads":
		peStatus = value

	if param == "peakType":
		peakType = value
		
fileResults = os.listdir(alignmentFolder)

submitAll = "master_htseq_queue.sh"
masterHandle = open(submitAll,"w")
text = "#!/bin/bash\n"
masterHandle.write(text)

jobCount = 0

for file in fileResults:
	if peakType == "broad":
		suffix = ".nodup.bam$"
	elif peakType == "narrow":
		suffix = ".bam$"
	else:
		print "'peakType' must be 'broad' or 'narrow'"
		sys.exit()
		
	result = re.search(suffix,file)
	fullPath = os.path.join(alignmentFolder, file)
	
	if result:
		sample = re.sub(suffix,"",file)
		sortResult = re.search(".name.sort.bam",file)
		
		if not sortResult:
			jobCount += 1
		
		if (sample not in finishedSamples) and (not sortResult):
			print sample
	
			shellScript = "htseq_" + sample + ".sh"
			text = "qsub " + shellScript + "\n"
			masterHandle.write(text)
			
			outHandle = open(shellScript, "w")
			text = "#!/bin/bash\n"
			text = text + "#$ -M "+email+"\n"
			text = text + "#$ -m bea\n"
			text = text + "#$ -N cwHT"+str(jobCount)+"\n"
			text = text + "#$ -q short.q\n"
			text = text + "#$ -j yes\n"
			text = text + "#$ -o cwHT"+str(jobCount)+".log\n"
			text = text + "#$ -cwd\n"
			text = text + "#$ -V\n"
			outHandle.write(text)

			#when I've tested single-end RNA-Seq, I've gotten the same results for name and position sorted .bam (and manual says parameter is ignored for single-end data)
			#However, for paired-end RNA-Seq, there can be differences, and the name-sorted .bam seems to be a little more accurate
			if peStatus != "yes":
				print "Add SE code"
				sys.exit()
			else:
				nameSortedBam = sample + ".name.sort.bam"
				sortPrefix = re.sub(".bam$","",nameSortedBam)
				text = "samtools sort -n " + fullPath + " " + sortPrefix + "\n"
				outHandle.write(text)

				countsFile = sample + "_merged_peak_counts.txt"
				text = "htseq-count -f bam -s no -t peak -i " + gtfID + " " + nameSortedBam + " " + peakGTF + " > " + countsFile + "\n"
				outHandle.write(text)

				countsFile = sample + "_promoter_counts.txt"
				text = "htseq-count -f bam -s no -t promoter -i " + gtfID + " " + nameSortedBam + " " + tssGTF + " > " + countsFile + "\n"
				outHandle.write(text)

				text = "rm " + nameSortedBam + "\n"
				outHandle.write(text)