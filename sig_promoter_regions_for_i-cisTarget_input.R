
dba.file = "[comp id]_[criteria].txt"
upID = "Trt Increased Binding"
downID = "Trt Decreased Binding"

upBed = "../Results/i-cisTarget_motifs/input_files/merged_peaks_2kb_[comp id]_UP.bed"
downBed = "../Results/i-cisTarget_motifs/input_files/merged_peaks_2kb_[comp id]_DOWN.bed"

input.table = read.delim(dba.file, head=T, sep="\t")

print(dim(input.table))
input.table = input.table[((input.table$annotation == "Promoter (1-2kb)")|(input.table$annotation == "Promoter (<=1kb)"))&!is.na(input.table$annotation),]
print(dim(input.table))

up.table = input.table[(input.table$status == upID),]
up.table = up.table[,2:4]
write.table(up.table, upBed, col.names=F, row.names=F, quote=F, sep="\t")

down.table = input.table[(input.table$status == downID),]
down.table = down.table[,2:4]
write.table(down.table, downBed, col.names=F, row.names=F, quote=F, sep="\t")
