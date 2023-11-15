pheatmap2<-function(Input = Null,
		Output = Null,
		main = "main"
		)
{
	library(pheatmap)
	d<-read.table(Input,header=T,stringsAsFactors=FALSE, row.names = 1,sep="\t",check.names=F)
	srtbottom=max(strwidth(colnames(d),units="inches",cex=1.5,font=2))+0.1
	strright=max(strwidth(rownames(d),units="inches",cex=1.5,font=2))+0.1
	figheight=nrow(d)*(strheight("P",units="inches",cex=1.5,font=2)+0.15)
	graphics.off()
	unlink("Rplots.pdf")
	height=srtbottom+figheight+1
	width=strright+0.8*ncol(d)+8
	pdf(Output,w=width,h=height)
#par(font=2,font.axis=2,font.lab=2,cex.axis=1.5,cex.main=2,mai=c(srtbottom,1.5,1,strright))
	par(font=2,font.axis=2,font.lab=2,cex.axis=1.5,cex.main=2)
		d[d>0.05]=0.06
		max_d = max(d)*0.99
		pheatmap(as.matrix(d),display_numbers = F,color=colorRampPalette(rev(c("gray","yellow","#F58020","#F99F1B","red")))(10),scale="none",main=main,fontsize=20,legend=T,cellheight=20,cellwidth=30,fontsize_row=20,fontsize_col=20,cluster_rows = F,cluster_cols=F,legend_breaks=c(min(d),max_d,0.06),legend_labels=c(0,max_d,1))

		dev.off()
}

args=commandArgs(T)
	if (length(args) !=3){
		print ("Rscript FunctionEnrichment.r  <Input> <Output> <main>")
		print ("Example : /annoroad/bioinfo/PMO/zhaohm/bin/R-3.0.1/bin/Rscript FunctionEnrichment.r GODescription.xls GODescription.pdf GO")
		q()
	}
#source("/annoroad/bioinfo/PMO/zhaohm/Rexample/GO/pheatmap.r")
pheatmap2(Input=args[1],Output=args[2],main=args[3])

