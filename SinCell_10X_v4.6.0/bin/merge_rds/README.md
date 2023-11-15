# 10X filter_cell
* 模块功能：使用Seurat软件按比较组样本合并多样本rds
* 模块版本： v0.0.1
* 作者：leiguo
* 邮箱：leiguo@genome.cn

### 软件环境
* R: v4.0.3
* Seurat: v4.1.1

### 资源消耗
随合并样本个数变化(9个样本)：
* 25G,  1cpu,  ~60m

### 使用方法
make -f this_make rds_dir= cmp_name= integrating_dir= config_file= Merge
target：
* Merge：按比较组样本合并多样本rds，并整合分析成一个rds文件保存

参数：
* rds_dir: 样本过滤后rds存放目录，sample/sample_filter_cell.rds 【必选|路径】；
* cmp_name: 比较组名称，示例A_VS_B 【必选|字符串】；
* integrating_dir：合并rds输出目录 【必选|路径】；
* config_file：比较组及相关参数配置文件 【必选|路径】。


### 输入文件示例
input/
config.ini
H1/H1_filter_cell.rds  #过滤后rds
P1/P1_filter_cell.rds

### 输出文件示例
output/
1_QC/P_VS_H.rds  #按比较组合并，且整合分析后rds
    
### 输出文件
1. cmp.rds
按比较组合并，且整合分析后rds，用于后续分析

### 注意事项
无

