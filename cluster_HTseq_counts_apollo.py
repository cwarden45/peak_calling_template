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

	if param == "Remove_Duplicates":
		dupFlag = value
		
fileResults = os.listdir(alignmentFolder)

submitAll = "master_htseq_queue.sh"
masterHandle = open(submitAll,"w")
text = "#!/bin/bash\n"
masterHandle.write(text)

jobCount = 0

for file in fileResults:
	if dupFlag == "yes":
		suffix = ".nodup.bam$"
	elif dupFlag == "no":
		suffix = ".bam$"
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
	
			shellScript = "htseq_" + sample + ".sh"
			text = "sbatch " + shellScript + "\n"
			masterHandle.write(text)
			
			outHandle = open(shellScript, "w")
			text = "#!/bin/bash\n"
			text = text + "#SBATCH -J cwHT"+str(jobCount)+"\n"
			text = text + "#SBATCH --mail-type=ALL\n"
			text = text + "#SBATCH --mail-user="+email+"\n"
			text = text + "#SBATCH -n 1\n"#one thread
			text = text + "#SBATCH -N 1\n"
			text = text + "#SBATCH --mem=4g\n"
			text = text + "#SBATCH --time=12:00:00\n"
			text = text + "#SBATCH --output=cwHT"+str(jobCount)+".log\n\n"
			
			#text = text + "module load Python/2.7.14-foss-2017a\n\n"
			#code above is supposed to be sufficient for htseq-count, but I need to use alternative set-up below
			text = text + "module load Python/3.6.1-foss-2017a\n"
			#installed with "python setup.py install --user", with git code
			#with executable /net/isi-dcnl/ifs/user_data/Seq/cwarden/More_Programs/htseq/build/scripts-3.6/htseq-count
			text = text + "module load samtools/1.6\n\n"
			outHandle.write(text)

			#when I've tested single-end RNA-Seq, I've gotten the same results for name and position sorted .bam (and manual says parameter is ignored for single-end data)
			#However, for paired-end RNA-Seq, there can be differences, and the name-sorted .bam seems to be a little more accurate
			if peStatus != "yes":
				countsFile = sample + "_merged_peak_counts.txt"
				text = "/net/isi-dcnl/ifs/user_data/Seq/cwarden/More_Programs/htseq/build/scripts-3.6/htseq-count -f bam -r pos -s no -t peak -i " + gtfID + " " + fullPath + " " + peakGTF + " > " + countsFile + "\n"
				outHandle.write(text)

				countsFile = sample + "_promoter_counts.txt"
				text = "/net/isi-dcnl/ifs/user_data/Seq/cwarden/More_Programs/htseq/build/scripts-3.6/htseq-count -f bam -r pos -s no -t promoter -i " + gtfID + " " + fullPath + " " + tssGTF + " > " + countsFile + "\n"
				outHandle.write(text)
			else:
				nameSortedBam = sample + ".name.sort.bam"
				text = "/opt/SAMtools/1.6/bin/samtools sort -n -o "+nameSortedBam+ " " + fullPath + "\n"
				outHandle.write(text)

				countsFile = sample + "_merged_peak_counts.txt"
				text = "/net/isi-dcnl/ifs/user_data/Seq/cwarden/More_Programs/htseq/build/scripts-3.6/htseq-count -f bam -s no -t peak -i " + gtfID + " " + nameSortedBam + " " + peakGTF + " > " + countsFile + "\n"
				outHandle.write(text)

				countsFile = sample + "_promoter_counts.txt"
				text = "/net/isi-dcnl/ifs/user_data/Seq/cwarden/More_Programs/htseq/build/scripts-3.6/htseq-count -f bam -s no -t promoter -i " + gtfID + " " + nameSortedBam + " " + tssGTF + " > " + countsFile + "\n"
				outHandle.write(text)

				text = "rm " + nameSortedBam + "\n"
				outHandle.write(text)