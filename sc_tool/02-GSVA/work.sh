#第一步进行GSVA分析


GSVA:
    mkdir -p $(outdir)
    $(Rscript) $(script_dir)/gsva_3rna.r -r $(rda) -d $(database) -o $(outdir)/$(prefix).gsva.csv

