#/annoroad/data1/software/install/anaconda3/Anaconda3-2023.03-1/bin/python
"""
QC CellRanger Count
Input dir: /path/CellRangerCount/sample/
the basename sample as the sample-name
the tree follow filess should be in the dir:
    - ./outs/filtered_gene_bc_matrices/barcodes.tsv.gz
    - ./outs/filtered_gene_bc_matrices/features.tsv.gz
    - ./outs/filtered_gene_bc_matrices/matrix.mtx.gz

conf.yml contents like the followings

metadata:
  species: human #human, mouse, other

QC:
  signatures:
    HB_score:
      - Hbb-bt
      - Hbb-bs
      - Hbb-bh3
    MT_score: 
      - mt-1
      - mt-2
      - mt-3
  gene_starts:
    ribo_score: 
      - Rpl
      - Rps
  min_genes: 200
  min_cells: 3
"""

import argparse
import os, sys
import matplotlib.pyplot as plt
import scanpy as sc 
import scrublet as scr
import pandas as pd
Bin = os.path.abspath(os.path.dirname(__file__))


def human_qc_mt_ribo_Hb(andata, Vlnkeys, min_genes=200, min_cells=3):
    # 如果1个gene 在少于3个细胞中表达，则过滤掉这个gene
    # 如果一个细胞表达少于200个基因，则过滤掉这个细胞
    sc.pp.filter_cells(andata, min_genes=min_genes, inplace=True)
    sc.pp.filter_genes(andata, min_cells=min_cells, inplace=True) 
    
    andata.var['MT'] = andata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    andata.var['ribo'] = andata.var_names.str.startswith(('RPL', 'RPS'))
    andata.var['HB'] = andata.var_names.str.contains(("^HB[^(P)]"))
    
    sc.pp.calculate_qc_metrics(andata, qc_vars=['MT', 'ribo', 'HB'], percent_top=None, log1p=False, inplace=True)
    
    
    
    Vlnkeys.append('pct_counts_MT')
    Vlnkeys.append('pct_counts_ribo')
    Vlnkeys.append('pct_counts_HB')
    
    return andata, Vlnkeys

def mouse_qc_mt_ribo_Hb(andata, Vlnkeys, min_genes=200, min_cells=3):
    # 如果1个gene 在少于3个细胞中表达，则过滤掉这个gene
    # 如果一个细胞表达少于200个基因，则过滤掉这个细胞
    sc.pp.filter_cells(andata, min_genes=min_genes, inplace=True)
    sc.pp.filter_genes(andata, min_cells=min_cells, inplace=True) 
    
    andata.var['MT'] = andata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(andata, qc_vars=['MT'], percent_top=None, log1p=False, inplace=True)
    
    andata.var['ribo'] = andata.var_names.str.startswith(('Rpl', 'Rps'))
    # annotate the group of ribo genes as 'ribo'
    sc.pp.calculate_qc_metrics(andata, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True)
    
    #Hb_genes = ["Hba1","Hba2","Hbb","Hbd","Hbe1","Hbg1","Hbg2","Hbm","Hbq1","Hbz"]
    Hb_genes = ["Hbb-bt","Hbb-bs","Hbb-bh3","Hbb-bh2","Hbb-bh1","Hbb-bh0","Hbb-y","Hbs1l","Hba-x","Hba-a1","Hbq1b","Hba-a2","Hbq1a","Hbp1","Hba-ps4","Hbegf"]
    andata.var['HB'] = andata.var_names.isin(Hb_genes)  # annotate the group of mitochondrial genes as 'mt'
    #andata.var['HB'] = andata.var_names.str.startswith(('HB'))
    sc.pp.calculate_qc_metrics(andata, qc_vars=['HB'], percent_top=None, log1p=False, inplace=True)
    Vlnkeys.append('pct_counts_MT')
    Vlnkeys.append('pct_counts_ribo')
    Vlnkeys.append('pct_counts_HB')
    
    return andata, Vlnkeys

def calculate_by_conf(andata, Vlnkeys, confdic):
    if 'signatures' in confdic['QC'].keys():
        for k, v in confdic['QC']['signatures'].items():
            andata.var[k] = andata.var_names.isin(v)
            sc.pp.calculate_qc_metrics(andata, qc_vars=[k], percent_top=None, log1p=False, inplace=True)
            Vlnkeys.append('pct_counts_{0}'.format(k))

    if 'gene_starts' in confdic['QC'].keys():
        for k, v in confdic['QC']['gene_starts'].items():
            andata.var[k] = andata.var_names.str.startswith(tuple(v))
            sc.pp.calculate_qc_metrics(andata, qc_vars=[k], percent_top=None, log1p=False, inplace=True)
            Vlnkeys.append('pct_counts_{0}'.format(k))
    
    return andata, Vlnkeys

def expected_doublet_rate(cell_num):
    exp_d_rate = 0.023
    
    if cell_num > 18000:
        exp_d_rate = 0.14
    elif cell_num > 15000:
        exp_d_rate = 0.1
    elif cell_num > 10000:
        exp_d_rate = 0.076
    elif cell_num > 9000:
        exp_d_rate = 0.069
    elif cell_num > 8000:
        exp_d_rate = 0.061
    elif cell_num > 7000:
        exp_d_rate = 0.054
    elif cell_num > 6000:
        exp_d_rate = 0.046
    elif cell_num > 5000:
        exp_d_rate = 0.039
    elif cell_num > 4000:
        exp_d_rate = 0.031
    elif cell_num > 3000:
        exp_d_rate = 0.023
    elif cell_num <= 3000:
        exp_d_rate = 0.016
    
    return(exp_d_rate)


def do_doublets_analysis(sample, sample_adata, outdir):
    print(sample , " doublets analysis start...")
   
    #sample_adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
    
    #counts_matrix = sample_adata.X
    #genes = sample_adata.var

    #print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
    #print('Number of genes in gene list: {}'.format(len(genes)))
    #print('Expected doublet rate: {0}'.format(expected_doublet_rate(counts_matrix.shape[0])))

    #scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=expected_doublet_rate(counts_matrix.shape[0]))
    scrub = scr.Scrublet(sample_adata.X,  expected_doublet_rate=0.06)
    #doublet_scores, predicted_doublets = scrub.scrub_doublets()
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    #scrub.call_doublets()
    
    fig = plt.figure(figsize=(8,4),constrained_layout=False, facecolor='white')
    scrub.plot_histogram()
    plt.savefig(outdir + "/{0}_scrublets_histogram.pdf".format(sample))
    
    #print('Running UMAP...')
    #scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    #print('UMAP Done.')
    
    #sample_adata.obs["predicted_doublets_"] = scrub.predicted_doublets_
    #sample_adata.obs["doublet_scores_obs_"] = scrub.doublet_scores_obs_
    
    #fig = plt.figure(figsize=(8,4),constrained_layout=False, facecolor='white')
    #scrub.plot_embedding('UMAP', order_points=True)
    #plt.savefig(outdir + "/{0}_scrublets_UMAP.pdf".format(sample))
    
    #sample_adata.write(outdir + "/{0}.scrublet.h5ad".format(sample))
    print('{0} finished doublets'.format(sample))
    b=predicted_doublets
    ratio = list(b).count(True)/len(sample_adata.obs.index)
    print("The {0} predicted doublets Ratio:{1}".format(sample, ratio))
    barcode=list(sample_adata.obs.index) 
    df2=pd.DataFrame({'Barcode':barcode,'Doublet_scores':doublet_scores,'Predicted_doublets':predicted_doublets})
    df2.to_csv("{}/{}_Predicted_doublet_scores.csv".format(outdir, sample),index=False)
    
    with open('{}/{}_doublet_ratio.xls'.format(outdir, sample), 'w') as fi:
        fi.write('Sample\t{}\n'.format(sample))
        fi.write('DoubletCell\t{:.3f}%\n'.format(ratio*100))
    
    #return scrub.predicted_doublets_, scrub.doublet_scores_obs_
    return predicted_doublets, doublet_scores
  
def do_cell_cycle(andata, species='human'):
    cell_cycle_genes = [x.strip() for x in open(os.path.join(Bin,'./regev_lab_cell_cycle_genes.txt'))]
    if species=='mouse':
        cell_cycle_genes = [x.capitalize() for x in cell_cycle_genes]

    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in andata.var_names]

    sc.pp.normalize_per_cell(andata, counts_per_cell_after=1e4)
    sc.pp.log1p(andata)
    sc.pp.scale(andata)
    
    sc.tl.score_genes_cell_cycle(andata, s_genes=s_genes, g2m_genes=g2m_genes)
    return andata


    
    
def Vlnplot_Ks(Andata, Vlnkeys, sample, outdir):
    fig, axs = plt.subplots(1, len(Vlnkeys), figsize=(14,4))
    
    for i, v in enumerate(Vlnkeys):
        sc.pl.violin(Andata, v,
                     jitter=0, multi_panel=False, show=False, ax=axs[i])
        if i == 0:
            axs[i].set_ylabel(sample)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, '{0}_QC.pdf'.format(sample)))
    
    
    
    
def main():
    parser=argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='')
    parser.add_argument('-i','--indir',help='sample CellRanger_Count dir, basename as sample name', dest='indir', type=str, required=True)
    parser.add_argument('-o','--outdir',help='out put dir',dest='outdir', type=str, required=True)
    parser.add_argument('-s','--sample',help='sample name', dest='sample', type=str, required=True)
    parser.add_argument('-p','--species',help='species', dest='species', type=str, required=True)
    parser.add_argument('--min_genes', help='min_genes of qc', dest='min_genes', default=200)
    parser.add_argument('--min_cells', help='min_cells of qc', dest='min_cells', default=3)
    
    args=parser.parse_args()
    sample = args.sample
    outdir = args.outdir
    species = args.species
    min_genes = int(args.min_genes)
    min_cells = int(args.min_cells)
    
    print('read single cell data...')
    adata = sc.read_10x_mtx(args.indir)
                            #var_names='gene_symbols',
                            # use gene symbols for the variable names (variables-axis index)
                            #cache=True

    vlnkeys = ['n_genes_by_counts', 'total_counts']

    if species == '9606':
        adata, vlnkeys = human_qc_mt_ribo_Hb(adata, vlnkeys, min_genes=min_genes, min_cells=min_cells)
        andata = do_cell_cycle(adata.copy(), species='human')
        adata.obs['S_score'] = andata.obs['S_score']
        adata.obs['G2M_score'] = andata.obs['G2M_score']
        adata.obs['phase'] = andata.obs['phase']
        #vlnkeys.append('S_score')
        #vlnkeys.append('G2M_score')
    elif species == '10090':
        adata, vlnkeys = mouse_qc_mt_ribo_Hb(adata, vlnkeys, min_genes=min_genes, min_cells=min_cells)
        andata = do_cell_cycle(adata.copy(), species='mouse')
        adata.obs['S_score'] = andata.obs['S_score']
        adata.obs['G2M_score'] = andata.obs['G2M_score']
        adata.obs['phase'] = andata.obs['phase']
        #vlnkeys.append('S_score')
        #vlnkeys.append('G2M_score')
        
    else:
        pass

    #adata, vlnkeys = calculate_by_conf(adata, vlnkeys, conf)

    predicted_doublets_, doublet_scores_obs_ = do_doublets_analysis(sample, adata, outdir)
    #print(pd.DataFrame(predicted_doublets_).head(2))
    #print(pd.DataFrame(doublet_scores_obs_).head(2))
    
    adata.obs["predicted_doublets"] = predicted_doublets_
    adata.obs["doublet_scores"] = doublet_scores_obs_
    vlnkeys.append('doublet_scores')
    
    adata.write(outdir + "/{0}.QC.h5ad".format(sample))
    Vlnplot_Ks(adata, vlnkeys, sample, outdir)



if __name__ == '__main__':
    main()


    
