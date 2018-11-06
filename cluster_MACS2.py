import sys
import re
import os

genomeRef = "hs"
#NOTE #1: may need to change species above
#NOTE #2: I've added a variable 'macs_dup_use' to allow MACS to use up to a certain number of duplicates that it defined (but I think it is usually best to leave this alone)

finishedSamples = ()
parameterFile = "parameters.txt"

alignmentFolder = ""
email = ""
peStatus = ""
peakType = ""
dupFlag = ""

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

	if param == "Cluster_Email":
		email = value

	if param == "PE_Reads":
		peStatus = value			

	if param == "peakType":
		peakType = value

	if param == "Remove_Duplicates":
		dupFlag = value
		
fileResults = os.listdir(alignmentFolder)

submitAll = "master_macs2_queue.sh"
masterHandle = open(submitAll,"w")
text = "#!/bin/bash\n"
masterHandle.write(text)

jobCount = 0

for file in fileResults:
	suffix = ".bam$"
	macs_dup_use=""
	if dupFlag == "yes":
		suffix = ".nodup.bam$"
	elif dupFlag == "no":
		macs_dup_use=" --keep-dup 100"
	else:
		print "'Remove_Duplicates' must be 'yes' or 'no'"
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
	
			shellScript = "macs2_" + sample + ".sh"
			text = "qsub " + shellScript + "\n"
			masterHandle.write(text)
			
			outHandle = open(shellScript, "w")
			text = "#!/bin/bash\n"
			text = text + "#$ -M "+email+"\n"
			text = text + "#$ -m bea\n"
			text = text + "#$ -N macs"+str(jobCount)+"\n"
			text = text + "#$ -q short.q\n"
			text = text + "#$ -j yes\n"
			text = text + "#$ -o macs"+str(jobCount)+".log\n"
			text = text + "#$ -cwd\n"
			text = text + "#$ -V\n"
			outHandle.write(text)

			detectPE = " "
			if peStatus == "yes":
				detectPE = " -f BAMPE "
			
			outputPrefix = alignmentFolder + "/" + sample + "/macs_" + sample
			
			if peakType == "broad":
				text = "/net/isi-dcnl/ifs/user_data/Seq/MACS_2.1.1.20160309/bin/macs2 callpeak -t "+ fullPath+" -n "+outputPrefix+" -g "+genomeRef+detectPE+"--broad --broad-cutoff 0.1\n"
				outHandle.write(text)
			elif peakType == "narrow":
				text = "/net/isi-dcnl/ifs/user_data/Seq/MACS_2.1.1.20160309/bin/macs2 callpeak -t "+ fullPath+" -n "+outputPrefix+" -g "+genomeRef+detectPE+"-q 0.1\n"
				outHandle.write(text)
			else:
				print "'peakType' must be 'broad' or 'narrow'"
				sys.exit()