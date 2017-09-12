import sys
import re
import os

igvtools = "/opt/igvtools_2.3.91/igvtools"
finishedSamples = ()

alignmentFolder = "../hg19_Bowtie2_Alignment"
genome = "hg19"

fileResults = os.listdir(alignmentFolder)

for file in fileResults:
	bamResult = re.search("(.*).bam$",file)
	if bamResult:
		sample = bamResult.group(1)
		if sample not in finishedSamples:
			print sample
			fullPath = alignmentFolder + "/" + file

			tdf = fullPath + ".tdf"
			command = igvtools + " count -z 7 " + fullPath + " " + tdf + " " + genome
			os.system(command)