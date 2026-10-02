# 执行与验收 SOP

## 输入与环境

原始 FASTQ 保持只读。`inventory` 全量检查 gzip/FASTQ 四行结构、两端 read ID/读数，并记录压缩文件 SHA256；需要读取全部文件，耗时与测序量有关。`validate` 仅检查配置和文件存在，不能代替完整 inventory。

nf-core 固定 3.2.2，流程参数以本地已使用版本为基线；模板没有声明使用最新版本。默认仅标记 H3K27ac target duplicates，保持真实 paired fragments，不做 ATAC shift。需记录官方 test 的报告，再审阅真实数据。

每个 run 使用新 ID，保存参数、样本和配置快照，输出到独立目录。不能复用已有 run ID。Nextflow work 目录为 `work/<run_id>`。当前没有 resume 命令；重启旧 run 需要按其 `command.json` 和原 work 路径显式加入 `-resume`，保留实际执行记录。不要在没有完成验收时把运行产物当作正式结果。

## QC 决策

至少审核有效目标片段数、spike-in 比例、比对率、重复/复杂度、片段长度分布、背景对照有效深度、峰数/峰宽/FRiP、组内相关性与 PCA、代表基因浏览器信号。阈值依靶标与实验预先规定，不固定用一个数字判所有 CUT&Tag。样本排除必须有单独决策记录。

若 spike-in 默认计数与序列证据不符，使用 `config/profiles/spikein_overlap_diagnostic.config` 在固定 read 子集上与默认参数比较；该配置不自动启用。第一批 overlap 修正是案例证据，不能未经验证套到所有文库。保存全量计数、比对参数和实验加入阶段，再接受比例因子。

IgG 极浅或默认峰近空时，比较不同背景、重复处理与无对照策略。无对照结果保留诊断标签。按独立生物样本评估重复支持与 FRiP、黑名单、轨迹后决定峰集。当前模板计数入口接收已验收 master BED，不自动创建峰集。

## 原始计数

`master.bed` 必须为 BED4：chrom/start/end/unique_peak_id；0-based 半开区间。`fragments.tsv` 含 `sample_id` 和 `fragments_bed`，一行一个独立生物样本。BED 每行必须代表一个完整 paired fragment，不能把 R1/R2 各当一个 fragment。重叠多个峰的 fragment 将计入多个峰；建议 master 区间合并至不重叠并记录合并规则。

```bash
pixi run -e analysis python scripts/count_fragments.py --peaks master.bed --samples fragments.tsv --out results/counts/raw.tsv
```

## 差异分析

metadata.tsv 含 `sample_id biological_sample_id condition replicates_confirmed` 四列，制表符分隔；confirmed 写 TRUE。sample_id 必须匹配 counts 列，每个 biological_sample_id 仅出现一次。只支持 `~ condition` 两组对比；复杂 batch、配对或重复测量设计需扩展实现，不能直接套用。

```bash
pixi run -e analysis Rscript scripts/differential.R results/counts/raw.tsv metadata.tsv KD Control results/differential/model1 conventional
pixi run -e analysis Rscript scripts/differential.R results/counts/raw.tsv metadata.tsv KD Control results/differential/model2 spikein spikein.tsv
```

spikein.tsv 含 `sample_id spikein_fragments calibration_accepted`；accepted 写 TRUE。DESeq2 的 size factor 为 spike-in counts / 其几何平均数，归一化 counts 因此与该值相除。默认过滤 total count ≥10、FDR<0.05、|log2FC|≥1，目前 R 入口固定这些参数；project.json 的同名参数是计划元数据，修改后需同步修改 R 实现并记录。CPM BigWig 与 median-ratio DESeq2 是不同归一化。

## 注释与 ORA

```bash
pixi run -e analysis Rscript scripts/annotate.R master.bed reference.sqlite org.Hs.eg.db results/annotation.tsv
pixi run -e analysis Rscript scripts/enrichment.R gain_gene_ids.txt tested_universe.txt term2gene.tsv results/ORA.tsv
```

物种 OrgDb/TxDb 需显式安装并登记来源/版本/哈希，不随模板下载。注释启动子默认为 TSS±2kb。ORA 输入与背景是去重后的基因 ID，必须来自同一已测试峰宇宙、同一区域类别；TERM2GENE 为两列带表头的 term/gene 表，采用相同 ID 体系。gain/loss 分开，可有同一基因同时关联两种方向。保留全量结果和资源哈希。通路富集是峰关联基因分析，不是基因表达差异。

Motif 当前按规划执行外部 HOMER，尚未包装通用入口：目标与背景需匹配长度、GC 和同一测试峰宇宙，排除窗口重叠，保留实际 knownResults 文件，不能只凭退出码验收。

## 交付和存储

交付清单含实验设计、比较方向、结果/模型版本、QC 决策、完整结果、限制和版本。检查所有软链接及 IGV 路径可迁移。缓存/work 清理前保存任务命令、日志和清单，并保护 FASTQ、发布结果、引用的参考与 IGV 文件；当前没有自动删除入口。
