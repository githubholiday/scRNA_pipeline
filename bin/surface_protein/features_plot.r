args = commandArgs(T)
library(Seurat)
library(cowplot)
library(ggplot2)
rds1 = readRDS(args[1])
rds = subset(rds1, orig.ident == args[2])
prefix = args[2]
outdir = args[3]
ref = args[4]
gene = args[5]
ck = strsplit(args[5], ',')[[1]]

genes <- rownames(rds)
#head(genes)

a = read.csv(ref, header = T)

gene_exit = c()
for (i in a$name) {
  t <- gsub('_', '-', i)
  if (t %in% genes){
    gene_exit = c(gene_exit, i)
  }
}

gene_name = read.csv(gene, header = T)


for (i in gene_exit) {
  name = gene_name[gene_name$PRO_ID==i,]$Genename
  
  if (name %in% genes) {
    t <- gsub('_', '-', i)
    print(c(t, name))
    graph = paste(outdir, paste('FeaturePlot',i,'antibody.pdf', sep='_'), sep='/')
    p1 = FeaturePlot(object = rds, features = c(t, name, ck), cols = c("grey", "blue"), reduction = "umap")
    ggsave(p1, filename=graph, width=12, height=10)
  }
}

