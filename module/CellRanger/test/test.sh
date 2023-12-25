## 测试 Count_Matrix 

module load singularity/3.7.3

## 测试fq
singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/tuchengfang/01_Pipeline/Sif/scRNA//cellranger_7.1.0.sif make -f /work/share/acuhtwkcu9/liutao/seqwisdom/4_scRNA/module/CellRanger/makefile indir=/work/share/acuhtwkcu9/tuchengfang/04_Project/20231219SC01ZD005/Filter/Filter_Result/ outdir=/work/share/acuhtwkcu9/liutao/seqwisdom/4_scRNA/module/CellRanger/test/output log_file=/work/share/acuhtwkcu9/liutao/seqwisdom/4_scRNA/module/CellRanger/test/output/log.txt lib=clean fq_id=3D sample_id=U3d ref_dir=/work/share/acuhtwkcu9/tuchengfang/database//refdata-cellranger-Rattus_norvegicus.Rnor_6.0.89.chr data_type=fastq force_cell= Count_Matrix 

## 测试h5
singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/tuchengfang/01_Pipeline/Sif/scRNA//cellranger_7.1.0.sif make -f /work/share/acuhtwkcu9/liutao/seqwisdom/4_scRNA/module/CellRanger/makefile indir=/work/share/acuhtwkcu9/liutao/seqwisdom/4_scRNA/module/CellRanger/test/input outdir=/work/share/acuhtwkcu9/liutao/seqwisdom/4_scRNA/module/CellRanger/test/output log_file=/work/share/acuhtwkcu9/liutao/seqwisdom/4_scRNA/module/CellRanger/test/output/log.txt lib=clean fq_id=3D sample_id=U3d ref_dir=/work/share/acuhtwkcu9/tuchengfang/database//refdata-cellranger-Rattus_norvegicus.Rnor_6.0.89.chr data_type=h5 force_cell= Count_Matrix 


## 测试Stat
singularity exec --bind /work/share/acuhtwkcu9/:/work/share/acuhtwkcu9/ /work/share/acuhtwkcu9/tuchengfang/01_Pipeline/Sif/scRNA//seurat_py.sif make -f /work/share/acuhtwkcu9/liutao/seqwisdom/4_scRNA/module/CellRanger/makefile indir=/work/share/acuhtwkcu9/tuchengfang/04_Project/20231219SC01ZD005/Filter/Filter_Result/ outdir=/work/share/acuhtwkcu9/liutao/seqwisdom/4_scRNA/module/CellRanger/test/output log_file=/work/share/acuhtwkcu9/liutao/seqwisdom/4_scRNA/module/CellRanger/test/output/log.txt lib=clean fq_id=3D sample_id=U3d ref_dir=/work/share/acuhtwkcu9/tuchengfang/database//refdata-cellranger-Rattus_norvegicus.Rnor_6.0.89.chr species=11071   data_type=fastq Stat