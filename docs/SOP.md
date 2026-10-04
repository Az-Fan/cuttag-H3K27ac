# 模块 SOP

当前初始化、环境、自动计划、续跑和交付说明以 [OPERATIONS.md](OPERATIONS.md) 为准；可选方案证据见 [OPTION_VALIDATION.md](OPTION_VALIDATION.md)。

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
pixi run -e analysis Rscript scripts/differential.R results/counts/raw.tsv metadata.tsv KD Control results/differential/model1 conventional model_review.json
pixi run -e analysis Rscript scripts/differential.R results/counts/raw.tsv metadata.tsv KD Control results/differential/model2 spikein spikein.tsv model_review_spikein.json
```

spikein.tsv 含 `sample_id spikein_fragments calibration_accepted calibration_group`；accepted 写 TRUE。当前正式 `~ condition` 模型要求本次对比的 target 样本处于同一校准组；多个组分别居中后的系数不能进入忽略组别的单一模型。DESeq2 的 size factor 为 spike-in counts / 其几何平均数，归一化 counts 因此与该值相除。默认过滤 total count ≥10、FDR<0.05、|log2FC|≥1，R 入口支持文末列出的环境变量，project.json 的同名参数是计划元数据，执行独立 R 命令时须同步。CPM BigWig 与 median-ratio DESeq2 是不同归一化。

## 注释与 ORA

```bash
pixi run -e analysis Rscript scripts/annotate.R master.bed reference.sqlite org.Hs.eg.db results/annotation.tsv
pixi run -e analysis Rscript scripts/enrichment.R gain_gene_ids.txt tested_universe.txt term2gene.tsv results/ORA.tsv
```

物种 OrgDb/TxDb 需显式安装并登记来源/版本/哈希，不随模板下载。注释启动子默认为 TSS±2kb。ORA 输入与背景是去重后的基因 ID，必须来自同一已测试峰宇宙、同一区域类别；TERM2GENE 为两列带表头的 term/gene 表，采用相同 ID 体系。gain/loss 分开，可有同一基因同时关联两种方向。保留全量结果和资源哈希。通路富集是峰关联基因分析，不是基因表达差异。

Motif 当前按规划执行外部 HOMER，尚未包装通用入口：目标与背景需匹配长度、GC 和同一测试峰宇宙，排除窗口重叠，保留实际 knownResults 文件，不能只凭退出码验收。

## 交付和存储

交付清单含实验设计、比较方向、结果/模型版本、QC 决策、完整结果、限制和版本。检查所有软链接及 IGV 路径可迁移。缓存/work 清理前保存任务命令、日志和清单，并保护 FASTQ、发布结果、引用的参考与 IGV 文件；当前没有自动删除入口。

## v0.2 专项模块

`qc.official_test_run_id` 必须指向当前项目已成功的官方 test run，并匹配版本/profile，方可执行真实 alignment/production。

### 共识峰

峰 manifest 是 TSV：`group biological_sample_id peaks_bed`。peaks_bed 相对 manifest 目录解析。同一生物样本多文件先在样本内合并，因此不能增加支持度。`--fraction` 定义每个组需要的独立样本比例，向上取整；默认 2/3。只保留达到支持度的碱基区间，随后组间 union/merge。默认逐碱基扣除 blacklist 重叠部分，再按 `--min-width` 丢弃过短片段；`--blacklist-mode drop` 仅用于比较“任意重叠即丢弃整个区间”的敏感性结果。provenance 会记录被扣除和因长度丢弃的区间。

```bash
python scripts/consensus.py --manifest peaks.tsv --blacklist blacklist.bed --out results/peaks/consensus1
```

### spike-in 审计与正式 fragment count

manifest 为 TSV：`sample_id mode bowtie2_log`，可加 `calibration_group`，不同组分别计算相对系数。同一模式同一样本只能一行；log 相对于 manifest 目录解析。统计 concordantly exactly 1 与 >1 的 paired fragments；多重比对仍在计数中，报告要求审核 MAPQ/多重比对规则。不要把不同 read 子集、技术拆分单位与合并后生物样本混在同一校准集合。

```bash
python scripts/spikein_audit.py --manifest spikein_logs.tsv --out results/qc/spikein1 --equal-amount-confirmed --added-at after_tagmentation
```

`spikein_audit.py` 仅审计 Bowtie2 日志中的 exactly-1 与 multi pair 统计，不是正式校准计数。正式计数使用 `spikein_fragments.py`，要求输入 spike-in-only BAM；manifest 每个测序 unit 一行，必须显式填写 `calibration_group`，脚本校验 BAM 参考序列字典、proper primary pair、两端 MAPQ，并只计一次 paired fragment。技术 unit 汇总到 biological sample 后输出各 MAPQ 结果和组内系数。只有 `role=target` 的生物样本参与归一化系数计算；control/IgG 计数保留作审核，不参与 target DE 的几何均值。

```text
sample_id biological_sample_id unit_id condition group role calibration_group bam
KD1 KD1 KD1_lane1 KD KD_H3K27ac target prep_A data/lambda/KD1.bam
```

```bash
python scripts/spikein_fragments.py --manifest spikein_bams.tsv --reference data/external/lambda.fa --mapq 20 30 --out results/qc/spikein_fragments_01
```

正式差异输入应使用对应策略目录下的 `target_biological_sample_counts.tsv`，不要把 IgG/control 行并入 target DE size factors。`calibration_group` 是必填字段，不能省略后让不同制备/加入批次默认合并。

所有方案都保持 `calibration_accepted=False`。正式接受前需审阅相同输入读段、唯一比对/MAPQ 规则、target+spike-in 竞争比对、等量加样、加入时点、校准范围和完整数据计数。竞争比对诊断可用 `spikein_crossmap.py` 在固定配对 FASTQ 子集上比较唯一 target、唯一 spike-in、跨参考和 ambiguous pairs；子集结果不可代替全量 BAM 计数。

### Peak caller 诊断

manifest TSV 列为 `sample_id target_bam control_bam bam_policy`，目标 BAM 必须是 duplicate-retained；BAM 路径相对 manifest 所在目录。只跑 no-IgG 时 `control_bam` 可留空。默认比较 MACS2 narrow/broad × IgG 默认缩放/scale-to-large/no-IgG，共六种候选；`--shapes`、`--backgrounds`、`--duplicate-modes all auto` 可缩小或扩展网格。不会自动选择生产结果。完成的每个候选会导出 `artifact_sets/<strategy>.tsv`，供 `strategy_matrix.py` 继续下游比较。除峰数外，还要比较区间交叠、峰宽、FRiP、黑名单交叠、重复支持和下游结论敏感性。

```bash
python scripts/peak_diagnostics.py --manifest bams.tsv --gsize 2700000000 --out results/qc/caller1
# 添加 --execute 才执行；诊断参数需在执行前固定并记录。
```

### Motif 与可视化

目标/背景应预先限定为同一测试宇宙的远端候选区域，脚本不替代区域分类。生成峰中心等长窗口，样本集内贪心排除重叠，背景额外去除目标重叠窗口。默认为 200bp、目标/背景各至少 50 窗口，不足则跳过并记录。

```bash
python scripts/motif.py --target distal_gain.bed --background distal_tested_other.bed --fasta genome.fa --out results/motif/gain1
# 添加 --execute 使用 HOMER，并检查 knownResults.txt。
pixi run -e analysis Rscript scripts/visualize.R results/differential/model1 results/figures/model1
```

Motif 入口当前固定 vertebrates 已知 motif；其他物种需修改选择。FASTA 同路径 .fai 是必需输入，越界窗口会被排除并记录。HOMER 命令返回 0 也必须有非空 knownResults。可视化生成描述性 PCA、Spearman、火山图和最多 100 峰的行 Z-score 热图；不将该 PCA 称为 VST PCA。

差异阈值现在可通过 `CUTTAG_ALPHA`、`CUTTAG_LFC`、`CUTTAG_MIN_COUNT` 环境变量设置，默认 .05/1/10；完整模型方法文件记录实用值。独立 R 命令仍用环境变量；推荐 run_differential.py，它会从项目配置传参并核对审核记录中的同名阈值。

### 最终报告

release.json 含 `project`、`limitations` 字符串数组、`artifacts` 数组。每个 artifact 包含 `path`、可选 `label`、`required`（默认 true）；相对路径基于 release.json 所在目录。生成 HTML 审阅页与 SHA256 清单，缺少必需文件/悬空文件链接即失败。现在会实际复制所选文件/目录并修复包内 IGV 路径，可用 verify_release.py 在搬迁后重新校验；科学状态仍为 REVIEW_REQUIRED。

```bash
python scripts/finalize.py --manifest release.json --out results/release/review1
```

## paired fragment 和信号轨迹

```bash
pixi run -e analysis python scripts/fragments.py --bam sample.bam --out results/fragments/sample
```

脚本保留重复标记的 primary proper pairs，不做 Tn5 shift；MAPQ 阈值同时用于两端，移除 duplicate 需显式传参。BED 依据 proper pair 的 TLEN 还原片段跨度；对异常/零 TLEN 丢弃。需用 samtools flagstat/片段分布抽查和原 pipeline fragments 核对。

信号输入 manifest TSV: `sample_id fragments_bed spikein_scale`。CPM 以保留的 target paired fragments 为分母。可选 spikein_scale 是经审核的 track multiplier（通常为组内几何均数/样本 spike-in fragments）；不能直接把 DESeq2 sizeFactor 当轨迹 multiplier。IGV XML 使用相对 BigWig 路径，--genome 必须显式指定组装。spike-in 轨迹还须传 --review，绑定 track_manifest_sha256、sizes_sha256 及校准证据。

```bash
pixi run -e analysis python scripts/tracks.py --manifest tracks.tsv --sizes genome.fa.fai --genome hg38 --out results/signal/run1
```

## QC atlas 与阶段编排

QC 输入 TSV 列为 `sample_id target_fragments spikein_fraction duplication_rate frip peak_count median_fragment_length`；未知值留空并明确画作 unavailable。该 atlas 是描述性摘要，没有内置通用 PASS/FAIL 阈值。

```bash
pixi run -e analysis python scripts/qc_atlas.py --metrics qc_metrics.tsv --out results/qc/atlas1
```

可用 `workflow.py` 顺序运行声明好的命令数组，确保依赖步骤通过后才运行下游；不得以 shell 字符串传命令。每步骤有独立命令/日志/状态，适用于下游整理，不替代 nf-core 的 Nextflow 调度或人工 QC 签核。

如主机的 Pixi 包缓存不可写，使用项目内独立缓存：

```bash
PIXI_CACHE_DIR="$PWD/shared_cache/pixi_cache" pixi install --locked -e default
PIXI_CACHE_DIR="$PWD/shared_cache/pixi_cache" pixi install --locked -e analysis
PIXI_CACHE_DIR="$PWD/shared_cache/pixi_cache" pixi run -e analysis selftest
```

机器级镜像若持续 TLS 失败，需修正 Pixi mirror 配置或网络链路；不要绕过 TLS 校验。容器和 Nextflow 缓存同样保存在项目的 `shared_cache/`，该目录已忽略且不得提交。


## 核查后新增的必要判断（以本节为准）

`qc.production_review` 指向 JSON 审核记录。请先固定 config 和 samples，运行 inventory 保存完整输入清单，再生成草稿：

```bash
python scripts/prepare_review.py --config config/project.json --kind production --evidence path/to/qc_report.md --out reviews/production.json
python scripts/prepare_review.py --config config/project.json --kind model --counts counts.tsv --metadata metadata.tsv --evidence path/to/peak_and_design_review.md --out reviews/model.json
python scripts/run_differential.py --config config/project.json --counts counts.tsv --metadata metadata.tsv --review reviews/model.json --out results/differential/model1
```

生成器只绑定当前内容哈希，所有接受字段初始为 false。审核者填写 reviewer、reviewed_at、reason，并根据证据填写 decisions；记录存在不等于接受。先在 config 填好 production_review 路径和已经作出的实验/QC 决策，再生成其审核文件，避免生成后改配置使哈希失效。改变 inputs、证据文件或模型/阈值后需要新的审核，不能复用旧记录。

生产审核还绑定 results/qc/fastq_inventory.json 和 results/qc/reference_audit.json；执行前新生成的 FASTQ 清单必须与已审核清单哈希一致，因此相同路径下的数据替换也会失效。生产决定含 upstream_accepted、control_strategy_accepted；模型决定含 upstream_accepted、peaks_accepted、replicates_confirmed、simple_design_accepted。校准分析另需 spikein_counting_accepted、spikein_calibration_accepted、等量加入前提、加入时点和 calibration_scope；模型记录绑定 spikein_sha256。当前只实现等量 spike-in 和 ~condition，不支持时拒绝执行。

直接 differential.R 调用的最后一个参数现在必须是 review.json；它核对 counts/metadata/spikein 哈希、对比/归一化/阈值、审核字段和证据文件。配置化 Python 入口另外核对 config 哈希和配置中的生物样本映射。原始输入检查现在由真实 run 自动执行，不只检查文件存在。

workflow 可给每步指定 outputs 文件列表；任一缺失即算失败并拦截依赖。未声明 outputs 的通用命令只具有退出码检查，不应把它说成通过了产物完整性验收。


## v0.4.0 行为更新

`finalize.py` 现在实际复制产物并生成可点击导航，不再只是绝对路径索引。`workflow.py --resume` 会校验已声明输入和产物哈希。`peak_diagnostics.py` 默认使用 tools.json 固定容器，先运行 fetch_tools.py。`collect_qc.py` 从 fragments/peaks 产生 QC 表，`prepare_universe.py` 从完整统计结果和峰注释准备 ORA/motif 背景；操作细节见 OPERATIONS。
