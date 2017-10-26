import sys
import re
import os
import subprocess

parameterFile = "parameters.txt"

finishedSamples = []

java = "/net/isi-dcnl/ifs/user_data/Seq/jre1.8.0_121/bin/java"
jar_path = "/net/isi-dcnl/ifs/user_data/Seq/"

threads = ""
email = ""
readsFolder = ""
bowtieRef = ""
alignmentFolder = ""
pairedStatus = ""
memLimit = ""
dupFlag = ""

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

	if param == "Threads":
		threads = value
		
	if param == "Cluster_Email":
		email = value

	if param == "Bowtie2_Ref":
		bowtieRef = value

	if param == "Alignment_Folder":
		alignmentFolder = value
		
	if param == "PE_Reads":
		pairedStatus = value

	if param == "Remove_Duplicates":
		dupFlag = value
		
	if param == "MEM_Limit":
		memLimit = value

if (memLimit == "") or (memLimit == "[required]"):
	print "Need to enter a value for 'MEM_Limit'!"
	sys.exit()
		
if (readsFolder == "") or (readsFolder == "[required]"):
	print "Need to enter a value for 'Reads_Folder'!"
	sys.exit()

if (threads == "") or (threads == "[required]"):
	print "Need to enter a value for 'Threads'!"
	sys.exit()
	
if (alignmentFolder == "") or (alignmentFolder == "[required]"):
	print "Need to enter a value for 'Alignment_Folder'!"
	sys.exit()
	
if (bowtieRef == "") or (bowtieRef == "[required]"):
	print "Need to enter a value for 'Bowtie2_Ref'!"
	sys.exit()
	
if (email == "") or (email == "[required]"):
	print "Need to enter a value for 'Cluster_Email'!"
	sys.exit()

if (pairedStatus == "") or (pairedStatus == "[required]"):
	print "Need to enter a value for 'PE_Reads'!"
	sys.exit()	

if (dupFlag == "") or (dupFlag == "[required]"):
	print "Need to enter a value for 'Remove_Duplicates'!"
	sys.exit()
else:
	if (dupFlag != "yes") and (dupFlag != "no"):
		print "'Remove_Duplicates' must be 'yes' or 'no'!"
		sys.exit()
	
submitAll = "master_Bowtie2_queue.sh"
masterHandle = open(submitAll,"w")
text = "#!/bin/bash\n"
masterHandle.write(text)
	
fileResults = os.listdir(readsFolder)

jobCount = 0
for file in fileResults:
	result = re.search("(.*)_(S\d+)_L\d{3}_R1_001.fastq$",file)
	
	if result:
		sample = result.group(1)
		barcode = result.group(2)
		
		jobCount += 1
		
		if (sample not in finishedSamples):
			print sample
			shellScript = sample + ".sh"
			text = "qsub " + shellScript + "\n"
			masterHandle.write(text)
			
			outHandle = open(shellScript, "w")
			text = "#!/bin/bash\n"
			text = text + "#$ -M "+email+"\n"
			text = text + "#$ -m bea\n"
			text = text + "#$ -N bow"+str(jobCount)+"\n"
			text = text + "#$ -q all.q\n"
			text = text + "#$ -pe shared "+str(threads)+"\n"
			text = text + "#$ -l vf="+memLimit+"\n"
			text = text + "#$ -j yes\n"
			text = text + "#$ -o bow"+str(jobCount)+".log\n"
			text = text + "#$ -cwd\n"
			text = text + "#$ -V\n"
			outHandle.write(text)
			
			sampleSubfolder = alignmentFolder + "/" + sample
			text = "mkdir " + sampleSubfolder + "\n"
			outHandle.write(text)
									
			if (pairedStatus == "yes"):
				read1 = readsFolder + "/" + file
				read2 = re.sub("_R1_001.fastq$","_R2_001.fastq",read1)
			
				alnSam = sampleSubfolder + "/aligned.sam"
				#add "--local -X 2000" parameters for ATAC-Seq
				text = "bowtie2 --local -X 2000 -p "+ str(threads) + " -x " + bowtieRef + " -1 " + read1 + " -2 " + read2  + " > " + alnSam + "\n"
				outHandle.write(text)			
			elif(pairedStatus == "no"):
				read1 = readsFolder + "/" + file
			
				alnSam = sampleSubfolder + "/aligned.sam"
				#add "--local -X 2000" parameters for ATAC-Seq
				text = "bowtie2 --local -X 2000 -p "+ str(threads) + " " + bowtieRef + " -U " + read1 + " -S " + alnSam + "\n"
				outHandle.write(text)
			else:
				print "'PE_Reads' value must be 'yes' or 'no'"
				sys.exit()

			alnBam = sampleSubfolder + "/aligned.bam"
			text = "samtools view -bS " + alnSam + " > " + alnBam + "\n"
			outHandle.write(text)

			text = "rm " + alnSam + "\n"
			outHandle.write(text)

			tempDir = sampleSubfolder + "/tmp"
			text = "mkdir " + tempDir + "\n"
			outHandle.write(text)
			
			if dupFlag == "yes":
				rgBam = sampleSubfolder + "/rg.bam"
				text = java + " -Xmx" + memLimit + " -Djava.io.tmpdir="+ tempDir + " -jar "+jar_path+"picard-tools-2.5.0/picard.jar AddOrReplaceReadGroups I=" + alnBam + " O=" + rgBam + " SO=coordinate RGID=1 RGLB=ChIP-Seq RGPL=Illumina RGPU="+barcode+" RGCN=COH RGSM=" + sample + " CREATE_INDEX=true\n"
				outHandle.write(text)

				text = "rm " + alnBam + "\n"
				outHandle.write(text)
				
				statsFile = sampleSubfolder + "/alignment_stats.txt"
				text = "samtools flagstat " + rgBam + " > " + statsFile + "\n"
				outHandle.write(text)

				statFile = sampleSubfolder + "/idxstats.txt"
				text = "samtools idxstats " + rgBam + " > " + statFile + "\n"
				outHandle.write(text)
				
				duplicateMetrics = sampleSubfolder + "/MarkDuplicates_metrics.txt"
				filteredBam = alignmentFolder + "/" + sample + ".nodup.bam"
				text = java + " -Xmx" + memLimit + " -Djava.io.tmpdir="+ tempDir + " -jar "+jar_path+"picard-tools-2.5.0/picard.jar MarkDuplicates I=" + rgBam + " O=" + filteredBam + " M=" + duplicateMetrics+" REMOVE_DUPLICATES=true CREATE_INDEX=true\n"
				outHandle.write(text)

				statFile = sampleSubfolder + "/idxstats_no_dup.txt"
				text = "samtools idxstats " + filteredBam + " > " + statFile + "\n"
				outHandle.write(text)
				
				text = "rm " + rgBam + "\n"
				outHandle.write(text)
				
				text = "rm " + re.sub(".bam$",".bai",rgBam) + "\n"
				outHandle.write(text)	
			else:
				rgBam = alignmentFolder + "/" + sample + ".bam"
				text = java + " -Xmx" + memLimit + " -Djava.io.tmpdir="+ tempDir + " -jar "+jar_path+"picard-tools-2.5.0/picard.jar AddOrReplaceReadGroups I=" + alnBam + " O=" + rgBam + " SO=coordinate RGID=1 RGLB=ChIP-Seq RGPL=Illumina RGPU="+barcode+" RGCN=COH RGSM=" + sample + " CREATE_INDEX=true\n"
				outHandle.write(text)

				text = "rm " + alnBam + "\n"
				outHandle.write(text)
				
				statsFile = sampleSubfolder + "/alignment_stats.txt"
				text = "samtools flagstat " + rgBam + " > " + statsFile + "\n"
				outHandle.write(text)
				
				statFile = sampleSubfolder + "/idxstats.txt"
				text = "samtools idxstats " + rgBam + " > " + statFile + "\n"
				outHandle.write(text)

			text = "rm -R " + tempDir + "\n"
			outHandle.write(text)
				
			if (pairedStatus == "yes"):
				text = "gzip "+ read1 +"\n"
				outHandle.write(text)
			
				text = "gzip "+ read2 +"\n"
				outHandle.write(text)			
			elif(pairedStatus == "no"):
				text = "gzip "+ read1 +"\n"
				outHandle.write(text)