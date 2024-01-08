library(getopt)



para <- matrix(c(
    "help",   "h",	0,	"logical",
    "prefix",	"p",	1,	"character",
    "species",   "s",   1, "character",
    "input",	"i",	1,	"character"), byrow = TRUE, ncol = 4)

opt <- getopt(para, debug = FALSE)
print_usage <- function(para = NULL) {
    cat(getopt(para, usage = TRUE))
    cat("
    Options:
    --help		h	NULL		get this help
    --input 	i	character	indir for expression file[forced]
    --species   s   character   species, homo_sapiens or mus_musculus [forced]
    --prefix	p	character	the prefix for outputfiles [forced]
    \n")
    q(status = 1)
}
if (!is.null(opt$help))	{
    print_usage(para)
}
if (is.null(opt$input)) {
    cat("Please input the data file1 ...\n\n")
     print_usage(para)
}
if (is.null(opt$prefix)) {
    cat("Please speicfify the output prefix ...\n\n")
    print_usage(para)
}

if (is.null(opt$species)) {
    cat("Please input species ...\n\n")
     print_usage(para)
}

library(readr)
library(dplyr)
library(clusterProfiler)
library(enrichplot)

if (grep("homo_sapiens", opt$species, ignore.case = TRUE)) {
    library(org.Hs.eg.db)
    species_database <- org.Hs.eg.db
}else if (grep("mus_musculus", opt$species, ignore.case = TRUE)) {
    library(org.MM.eg.db)
    species_database <- org.MM.eg.db
}else {
     cat("Please right input species ...\n\n")
     print_usage(para)
}



#/annoroad/data1/bioinfo/PMO/suyalei/software/Anaconda/minconda3/envs/r_4/bin/R
input_table <- read_tsv(file = opt$input)
de_table <- filter(input_table, Significant == "yes")

ego2 <- enrichGO(gene         = de_table$Gene,
                OrgDb         = species_database,
                keyType       = "ENSEMBL",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1)

pdf(paste(opt$prefix, ".barplot.pdf", sep = ""), width = 10, height = 6.18)
barplot(ego2, showCategory = 20)
dev.off()

png(paste(opt$prefix, ".barplot.png", sep = ""), width = 10, height = 6.18, units="in",res=400,type ="cairo")
barplot(ego2, showCategory = 20)
dev.off()

pdf(paste(opt$prefix, ".dotplot.pdf", sep = ""), width = 10, height = 6.18)
dotplot(ego2, showCategory = 30)
dev.off()


png(paste(opt$prefix, ".dotplot.png", sep = ""), width = 10, height = 6.18, units="in",res=400,type ="cairo")
dotplot(ego2, showCategory = 30)
dev.off()

write.table(ego2,
            paste(opt$prefix, "go_enrichment.xls", sep = "_"),
            quote = FALSE , sep="\t" , row.names=FALSE)

print("go enrichment finished")

select_table <- input_table  %>% filter( pval<0.1 )
genelist <- select_table$Log2FoldChange
names(genelist) <- as.character(select_table$Gene)
genelist <- sort(genelist, decreasing = TRUE)

#ego3 <- gseGO(geneList     = genelist,
#              OrgDb        = species_database,
#              keyType      = "ENSEMBL",
#              ont          = "ALL",
#              minGSSize    = 10,
#              maxGSSize    = 500,
#              pvalueCutoff = 1,
#              verbose      = FALSE)
#write.table(ego3,
#            paste(opt$prefix, "gse_enrichment.xls", sep = "_"),
#            quote = FALSE,sep="\t" , row.names=FALSE)

#for (i in row(ego3)) {
#    GO_id = gsub( ":", "", ego3$ID[i])
#    #print(paste(opt$prefix, GO_id, "gse.pdf", sep = "."))
#    pdf(paste(opt$prefix, GO_id, "gse.pdf", sep = "."),  width = 10, height = 6.18)
#    p1 <- gseaplot(ego3, geneSetID = i, title = paste(ego3$ID[i], ego3$Description[i], sep = " "))
#    print(p1)
#    dev.off()
    
#}
