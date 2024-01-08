# single_cell updata_lims
* 模块功能: 单细胞分析结果指标按样本插入sci_bioinfo 10xgenomics和c4表
* 模块版本： v0.0.1
* 作者：leiguo
* 邮箱：leiguo@genome.cn

### 软件环境
* python3: v3.11.0

### 资源消耗
* 1G,  1cpu,  ~1m

### 使用方法
make -f mk_updata config_ini= sample_id= type= cellranger_stat= filter_stat= double_stat= software= Lims
target：
* Lims：插入lims

参数：
config_ini：流程配置config用于识别样本信息【必选|路径】；
sample_id: 样本名称【必选|字符串】；
type：流程类型，10X 或 C4【必选|字符串】;
cellranger_stat：cellranger的结果统计文件 【必选|路径】
filter_stat：细胞过滤的统计文件 【必选|路径】；
double_stat：双细胞的统计文件，双细胞未加到过滤中是单独的文件 【必选|路径】；
software：分析软件路径，用于记录软件版本信息 【必选|路径】


### 输入文件示例
input/
├── analysis_summary.xls
├── Cell_filter_stat.xls
├── config.ini
└── hCAF1_doublet_ratio.xls

输入统计文件要求：
    首列为统计指标
    后续列为样本和样本对应指标结果，按照列名和样本名对应查找数据。

### 输出文件示例
output/
    tb_analysis_sincell_rna  按项目号+样本插入lims数据库，已存在的进行覆盖。

### 输出文件
无

### 注意事项
1. 由于lims数据库字段同流程统计字段名称不统一，故script/updata_lims.py脚本增加配置文件-j/--json_config, 按照流程类型对字段名称进行一一对应。
   后续增加流程类型，需要增加该文件的配置。


