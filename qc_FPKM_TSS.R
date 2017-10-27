colLab = function(n, labelColors, clusMember) { 
   if(is.leaf(n)) { 
       a <- attributes(n) 
	   #print(a)
       # clusMember - vector of sample names (ordered to match label color.palette)
       # labelColors - a vector of color.palette for the above grouping 
       labCol <- labelColors[clusMember == a$label]
	   #print(labCol)
       attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol) 
   } 
   n 
}

count.defined.values <- function(arr)
{
	return(length(arr[!is.na(arr)]))
}#end def count.values

param.table = read.table("parameters.txt", header=T, sep="\t")
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
cpm.file = as.character(param.table$Value[param.table$Parameter == "fpkm_file_TSS"])
cluster.distance = as.character(param.table$Value[param.table$Parameter == "cluster_distance"])
plot.groups = unlist(strsplit(as.character(param.table$Value[param.table$Parameter == "plot_groups"]), split=","))
plot.types = unlist(strsplit(as.character(param.table$Value[param.table$Parameter == "plot_type"]), split=","))
min.expression = as.numeric(as.character(param.table$Value[param.table$Parameter == "fpkm_binding_cutoff"]))

fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black","darkgreen","thistle1","tan","orchid1",colors())

sample.table = read.table(sample.description.file, sep="\t", header=T)
userID = as.character(sample.table$userID)

normalized.table = read.table(cpm.file, sep="\t", header=T)
normalized.mat = normalized.table[,match(userID, names(normalized.table))]
expr.max = ceiling(max(normalized.mat, na.rm = T))
expr.min = ceiling(min(normalized.mat, na.rm = T))

quantiles = round(apply(normalized.mat, 2, quantile, probs=c(0.01,0.05,0.25,0.5,0.75,0.95,0.99)), digits=2)
quantiles = data.frame(Percentage = rownames(quantiles), quantiles)
write.table(quantiles, "TSS_log2_FPKM_quantiles.txt",sep="\t", quote=F, row.names=F)


for (i in 1:length(plot.groups)){
	group = plot.groups[i]
	group.type = plot.types[i]
	print(group)
	
	
	temp.mat = normalized.mat[,!is.na(sample.table[,group])]
	print(dim(temp.mat))
	qc.grp = sample.table[,group]
	qc.grp = qc.grp[!is.na(sample.table[,group])]
	clusterID = userID[!is.na(sample.table[,group])]
	
	#you can also try adding center=F, if you are concerned zero-centering is skewing results
	pca.values = prcomp(na.omit(data.matrix(temp.mat)))
	pc.values = data.frame(pca.values$rotation)
	variance.explained = (pca.values $sdev)^2 / sum(pca.values $sdev^2)
	pca.table = data.frame(PC = 1:length(variance.explained), percent.variation = variance.explained, t(pc.values))

	pca.text.file = paste("TSS_FPKM_",group,"_pca_values.txt",sep="")
	write.table(pca.table, pca.text.file, quote=F, row.names=F, sep="\t")
	
	groups = levels(as.factor(as.character(qc.grp)))
	num.sample.types = length(groups)
	color.palette <- fixed.color.palatte[1:length(groups)]

	labelColors = rep("black",times=ncol(temp.mat))
	if(group.type == "continuous"){
		library("RColorBrewer")
		plot.var = as.numeric(qc.grp)
		plot.var.min = min(plot.var, na.rm=T)
		plot.var.max = max(plot.var, na.rm=T)
		
		plot.var.range = plot.var.max - plot.var.min
		plot.var.interval = plot.var.range / continuous.color.breaks
		
		color.range = colorRampPalette(c("green","black","orange"))(n = continuous.color.breaks)
		plot.var.breaks = plot.var.min + plot.var.interval*(0:continuous.color.breaks)
		for (j in 1:continuous.color.breaks){
			#print(paste(plot.var.breaks[j],"to",plot.var.breaks[j+1]))
			labelColors[(plot.var >= plot.var.breaks[j]) &(plot.var <= plot.var.breaks[j+1])] = color.range[j]
		}#end for (j in 1:continuous.color.breaks)
	}else{
		for (j in 1:length(groups)){
			labelColors[qc.grp == as.character(groups[j])] = color.palette[j]
		}#end for (j in 1:length(groups))
	}

	pca.file = paste("TSS_FPKM_pca_by_",group,".png",sep="")
	png(file=pca.file)
	if(group.type == "continuous"){
		par(mar = par("mar") + c(0,0,0,5))
		plot(pc.values$PC1, pc.values$PC2, col = labelColors, xlab = paste("PC1 (",round(100* variance.explained[1] , digits = 2),"%)", sep = ""),
				ylab = paste("PC2 (",round(100* variance.explained[2] , digits = 2),"%)", sep = ""),
				pch=19, main=paste("Color by ",group,sep=""))
		legend("right",legend=c(round(plot.var.max,digits=1),rep("",length(color.range)-2),round(plot.var.min,digits=1)),
				col=rev(color.range),  pch=15, inset=-0.2, xpd=T, y.intersp = 0.4, cex=0.8, pt.cex=1.5)
	}else{
		par(mar = par("mar") + c(0,0,0,7))
		plot(pc.values$PC1, pc.values$PC2, col = labelColors, xlab = paste("PC1 (",round(100* variance.explained[1] , digits = 2),"%)", sep = ""),
				ylab = paste("PC2 (",round(100* variance.explained[2] , digits = 2),"%)", sep = ""),
				pch=19, main=paste("Color by ",group,sep=""))
		legend("right",legend=groups,col=color.palette,
				xpd=T, inset=-0.3, pch=19)
	}
	dev.off()

	box.file = paste("TSS_FPKM_box_plot_by_",group,".png",sep="")
	png(file=box.file)
	boxplot(temp.mat, col=labelColors, xaxt='n', main = group)
	legend("topright",legend=groups,col=color.palette,  pch=19)
	dev.off()

	cluster.file = paste("TSS_FPKM_cluster_by_",group,".png",sep="")
	if(cluster.distance == "Euclidean"){
		dist1 <- dist(as.matrix(t(temp.mat)))
	}else if (cluster.distance == "Pearson_Dissimilarity"){
		cor.mat = cor(as.matrix(temp.mat))
		dis.mat = 1 - cor.mat
		dist1 <- as.dist(dis.mat)
	}else{
		stop("cluster_distance must be 'Euclidean' or 'Pearson_Dissimilarity'")
	}
	clusMember <- groups
	hc <- hclust(dist1)
	dend1 <- as.dendrogram(hc)
	png(file = cluster.file)
	dend1 <- dendrapply(dend1, colLab, labelColors=labelColors, clusMember=clusterID) 
	a <- attributes(dend1) 
	attr(dend1, "nodePar") <- c(a$nodePar, lab.col = labelColors) 
	op <- par(mar = par("mar") + c(0,0,0,10)) 
	plot(dend1, horiz=T, main = group)
	par(op) 
	dev.off()

	hist.file = paste("TSS_FPKM_density_by_",group,".png",sep="")
	png(file = hist.file)
	for (i in 1:ncol(temp.mat))
		{		
			data <- as.numeric(t(temp.mat[,i]))
			
			if(i == 1)
				{
					den <- density(data, na.rm=T,from=expr.min, to=expr.max, adjust=2)
					expr <- den$x
					freq <- den$y
					plot(expr, freq, type="l", xlab = paste("Log2(FPKM + ",min.expression,")",sep=""), ylab = "Density",
							xlim=c(expr.min,expr.max), ylim=c(0,0.2), col=labelColors[i],
					    		main = group)
					if(group.type == "continuous"){
						legend("right",legend=c(round(plot.var.max,digits=1),rep("",length(color.range)-2),round(plot.var.min,digits=1)),
								col=rev(color.range),  pch=15, inset=-0.4, xpd=T, y.intersp = 0.4, cex=0.8, pt.cex=1.5)
					}else{
						legend("right",legend=groups,col=color.palette,  pch=19, inset=-0.6, xpd=T)
					}
				}#end if(i == 1)
			else
				{
					den <- density(data, na.rm=T,from=expr.min, to=expr.max, adjust=2)
					expr <- den$x
					freq <- den$y
					lines(expr, freq, type = "l", col=labelColors[i])
				}#end else
		}#end for (i in 1:length(ncol(temp.mat)))
	dev.off()	
}#end for (group in plot.groups)
