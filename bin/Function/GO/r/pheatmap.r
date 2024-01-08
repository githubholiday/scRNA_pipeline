pheatmap2<-function(Input = Null,
		Output = Null,
		main = "main",
		width = 8
		)
{
	library(pheatmap)
	d<-read.table(Input,header=T,stringsAsFactors=FALSE, row.names = 1,sep="\t",check.names=F)
	max_str_len <- max(nchar(row.names(d)))
	tt = row.names(d)
	row.names(d)[1] = paste( tt[1] , paste(rep(" ",  max_str_len - nchar(tt[1]) + 8 ), collapse='') , collapse='')
	row.names(d)[2] = paste( tt[2] , paste(rep(" ",  max_str_len - nchar(tt[2]) + 8 ), collapse='') , collapse='')

	srtbottom=max(strwidth(colnames(d),units="inches",cex=1.5,font=2))+0.1
	strright=max(strwidth(rownames(d),units="inches",cex=1.5,font=2))+0.1
	figheight=nrow(d)*(strheight("P",units="inches",cex=1.5,font=2)+0.1)
	print (main)
	titlelen=strwidth(main,units="inches",cex=1.5,font=2)
	#graphics.off()
	unlink("Rplots.pdf")
	height=srtbottom+figheight+1
	if (length(args) !=4) {
		width=strright * 1 +0.8*ncol(d)+titlelen+0.5
	}
	if (width > 20 && length(args) == 3){
		width = 20
	}
#	width=strright * 1 +0.8*ncol(d)+1.8
#	pdf(Output,w=width,h=height)
#	par(font=2,font.axis=2,font.lab=2,cex.axis=1.5,cex.main=2,mai=c(srtbottom,1.5,1,strright))
#	par(font=2,font.axis=2,font.lab=2,cex.axis=1.5,cex.main=1.5)
	d[d>0.05]=0.06
	pheatmap(as.matrix(d),display_numbers = F,color=colorRampPalette(rev(c("gray","yellow","#F58020","#F99F1B","red")))(10),scale="none",main=main,legend=T,cellheight=20,cellwidth=20,cluster_rows = F,cluster_cols=F,legend_breaks=c(min(d),0.05,0.06),legend_labels=c(0,0.05,1),filename=Output,width=width,height=height)

		dev.off()
}

args=commandArgs(T)
	if (length(args) !=3 && length(args) !=4){
		print ("Rscript FunctionEnrichment.r  <Input> <Output> <main> <width>(alternative)")
		print ("Example : /annoroad/bioinfo/PMO/zhaohm/bin/R-3.0.1/bin/Rscript FunctionEnrichment.r GODescription.xls GODescription.pdf GO 12")
		q()
	}
#source("/annoroad/bioinfo/PMO/zhaohm/Rexample/GO/pheatmap.r")
pheatmap2(Input=args[1],Output=args[2],main=args[3],width=args[4])

