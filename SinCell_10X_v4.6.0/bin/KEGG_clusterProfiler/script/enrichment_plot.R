# renxue
## draw enrichment picture
library(getopt)

para <- matrix(c(
    "help",   "h",	0,	"logical",
    "input",	"i",	1,	"character",
    "title",	"t",	1,	"character",
    "number",	"n",	1,	"integer",
    "facet",	"f",	1,	"logical",
    "ylab",	"y",	1,	"character",
    "prefix",	"p",	1,	"character"), byrow = TRUE, ncol = 4)
opt <- getopt(para, debug = FALSE)
print_usage <- function(para = NULL) {
    cat(getopt(para, usage = TRUE))
    cat("
    Options:
    --help		h	NULL		get this help
    --prefix	p	character		outfile prefix for enrichment report[forced]
    --input 	i	file	input file that the first column is genelist [forced]
    --title		t	character	title of the picture [optional,default is blank]
    --facet		t	NULL	whether facet based on ONTOLOGY [optional,default is T]
    --ylab		y	character	y-axis lable
    --number	n	int	 number of term to draw for picture [optional,default is all]
    \n")
    q(status = 1)
}
if (!is.null(opt$help))	{
    print_usage(para)
}
if (is.null(opt$input)) {
    cat("Please input the go enrichment result ...\n\n")
     print_usage(para)
}
if (is.null(opt$prefix)) {
    cat("Please speicfify the output prefix ...\n\n")
    print_usage(para)
}
if (is.null(opt$title)) {
	title = ""
}else{
	title = opt$title
}
if (is.null(opt$ylab)) {
	ylab = "Function"
}else{
	ylab = opt$ylab
}
if (is.null(opt$facet)) {
	facet = T
}else{
	facet = opt$facet
}

if (is.null(opt$number)) {
	number = -1
}else{
	number = opt$number
}
library(ggplot2)
library(dplyr)
data=read.table(opt$input,header=T,sep='\t',quote="")
if (facet == T ){
	data %>% arrange(ONTOLOGY,p.adjust) -> data1
}else{
	data %>% arrange(p.adjust) -> data1
}

if(number != -1 && length(data[,1]) > number ){
	print("total row number is greater than -n,so cut data")
	if (facet == T ){
		data2 = data1 %>% group_by(ONTOLOGY) %>% filter(row_number() <= number)
	}else{
		data2 = data1 %>% filter(row_number() <= number)
	}
}else{
	data2 =data1
}

data_final = data2

data_len = length(data_final[,1])
if (data_len <1){
	print("the input file has 0 line，so do not draw picture")
	q(status = 0)
}

max_description_length = max(sapply(as.vector(data_final[,1]) , nchar))
width_len = (data_len)*0.002+9 + max_description_length*0.5 
height_len = data_len*0.03+5
#height_len =9
#library(clusterProfiler)
#library(enrichplot)

## 计算
gr = apply(data_final,1,function(x){eval(parse(text=x["GeneRatio"]))})
br = apply(data_final,1,function(x){eval(parse(text=x["BgRatio"]))})
fc = gr/br
data2=data.frame("gr"=gr,"br"=br,"fc"=fc)
data_new = cbind(data_final,data2)
if (facet == T){
	data_new %>% arrange(ONTOLOGY,Count,Description)%>% mutate(order=row_number(),order_factor=factor(order,labels=Description)) -> data_sort1
}else{
	data_new %>% arrange(Count,Description)%>% mutate(order=row_number(),order_factor=factor(order,labels=Description)) -> data_sort1
}
p=ggplot(data = data_sort1) + geom_bar(aes(x=order_factor,y=Count,fill=p.adjust,color=p.adjust),stat = 'identity',position="dodge",width=0.8)
#p=p+facet_grid(.~ONTOLOGY,scales="free_x", space="free_x" )
if (facet == T){ p=p+facet_grid(ONTOLOGY~.,scales="free", space="free" )}
p=p+coord_flip()
p=p+labs(title=title,y="Gene Count", x=ylab)
p=p+scale_color_gradient(high = "blue",low = "red")
p=p+scale_fill_gradient(high = "blue",low = "red")
p=p+theme(plot.title=element_text(hjust=0.5))

pdf(paste(opt$prefix, ".barplot.pdf", sep = ""), width=width_len, height=height_len)
p
dev.off()

if (facet == T){
	data_new %>% arrange(ONTOLOGY,gr,Description)%>% mutate(order=row_number(),order_factor=factor(order,labels=Description)) -> data_sort2
}else{
	data_new %>% arrange(gr,Description)%>% mutate(order=row_number(),order_factor=factor(order,labels=Description)) -> data_sort2
}
d=ggplot(data = data_sort2) + geom_point(aes(x=order_factor,y=gr,size=Count,fill=p.adjust,color=p.adjust,))
if (facet == T){d=d+facet_grid(ONTOLOGY~.,scales="free", space="free" )}
d=d+coord_flip()
d=d+labs(title=title,y="GeneRatio", x=ylab)
d=d+scale_fill_gradient(high = "blue",low = "red")
d=d+scale_color_gradient(high = "blue",low = "red")
d=d+theme(plot.title=element_text(hjust=0.5))

pdf(paste(opt$prefix, ".dotplot.pdf", sep = ""), width=width_len, height=height_len)
d
dev.off()


##将生成png放在后边
print("开始生成png ......")
png(paste(opt$prefix, ".barplot.png", sep = ""), width=width_len, height=height_len, units="in",res=300,type ="cairo")
p
dev.off()

png(paste(opt$prefix, ".dotplot.png", sep = ""), width=width_len, height=height_len, units="in",res=300,type ="cairo")
d
dev.off()
