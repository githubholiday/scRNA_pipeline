s1 <- readRDS("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/03_Spatial_tr/03_Epi_Analysis/result/GAQ.epi.rds")
s2 <- readRDS("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/03_Spatial_tr/03_Epi_Analysis/result/HRE.epi.rds")
s3 <- readRDS("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/03_Spatial_tr/03_Epi_Analysis/result/LFL.epi.rds")
s4 <- readRDS("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/03_Spatial_tr/03_Epi_Analysis/result/LGJ.epi.rds")
s5 <- readRDS("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/03_Spatial_tr/03_Epi_Analysis/result/LH.epi.rds")
s6 <- readRDS("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/03_Spatial_tr/03_Epi_Analysis/result/STQ.epi.rds")
s7 <- readRDS("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/03_Spatial_tr/03_Epi_Analysis/result/YCD.epi.rds")
s8 <- readRDS("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/03_Spatial_tr/03_Epi_Analysis/result/ZFY.epi.rds")
combine <- merge(s1, y=c(s2,s3,s4,s5,s6,s7,s8))

outfile <- "/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/03_Spatial_tr/03_Epi_Analysis/result/merge.rds"
saveRDS( combine, outfile )