# 创建、运行、验收与交付

## 支持范围

当前默认是 Linux、paired-end CUT&Tag、H3K27ac、hg38、MACS2 narrow peaks、两组独立生物样本的 `~ condition`。技术拆分必须先合并为同一统计样本。批次/配对设计不自动退化为简单模型；需另行实现并验收。每组只有一个独立生物样本时仍可生成 QC、轨迹和描述性峰集，但正式 DESeq2 入口被阻止。

## 干净初始化

```bash
python3 scripts/init_project.py --out /path/to/new-study --project-id new_study
cd /path/to/new-study
```

只复制模板代码、配置、文档、测试、锁文件；不复制原始数据、结果、缓存和 Git 历史，也不继承已批准的实验决策。拒绝覆盖已有目录。`template_origin.json` 保留模板源文件哈希。新目录中的示例样本表需替换为真实登记。

## 环境准备与验收

```bash
pixi install --locked -e analysis
python3 scripts/fetch_tools.py --tool macs2
python3 scripts/install_bioc_data.py --out logs/bioc_data_setup
pixi run -e analysis python scripts/acceptance.py --out results/acceptance_01
pixi run -e analysis python scripts/cuttag.py test --run-id official_test_01
# 同一官方微型输入，另测模板 MACS2 narrow / 保留 target 重复 / CPM 组合
pixi run -e analysis python scripts/cuttag.py test --test-template-params --run-id template_params_01
```

MACS2 诊断默认使用 `config/tools.json` 中固定 SHA256 的 2.2.7.1 容器，与 nf-core 3.2.2 中的版本一致。本机发现原生 Bioconda 2.2.9.1 build 5 在真正调用 callpeak 时缺少 `__log_finite`；因此已移出默认依赖。需要 Apptainer 容器权限。`--native-macs2` 仅用于显式测试其他已验证的原生安装，不作为默认回退。

Bioconda 部分数据包需要安装后额外下载。`install_bioc_data.py` 只读取当前环境已登记的 post-link 包和来源，按上游 MD5 校验，保存 SHA256、安装日志及加载检查。不能仅凭 `pixi install` 返回或 `macs2 --version` 判断运行环境正常。系统 R 的通过不能代替 Pixi R 的通过。

遇到代理证书问题，可使用系统可信证书：`pixi install --locked -e analysis --tls-root-certs system`。不要关闭 TLS 校验。若本机全局镜像配置不可用，可显式 `--no-config`；项目缓存可设为 `PIXI_CACHE_DIR="$PWD/shared_cache/pixi_cache"`。

`acceptance.py` 是本地软件验收，`--quick` 仅跑契约测试。全套涵盖真实 samtools/bedtools/BigWig、干净初始化的下游工作流、4 组片段过滤策略、3 种 MACS2 控制方案、2 种 DESeq2 归一化、ORA、显式 TxDb/OrgDb 注释与 HOMER。官方 nf-core 测试单独运行；任何计算 PASS 都不等同于科学审核通过。

## 上游与参考

编辑 `config/project.json`、`config/samples.tsv` 后，按 SOP 执行 validate → inventory → alignment → review → production。inventory 现在同时产生 `results/qc/reference_audit.json`：校验 FASTA 序列、可选 FAI 的序列名/长度、GTF 与 blacklist 坐标范围、线粒体序列名，并对参考文件计算哈希。默认要求未压缩 FASTA/GTF/BED。坐标兼容不能证明 assembly/release 身份，应在审核中记录来源。

生产审核必须绑定 reference audit 和 FASTQ inventory；执行前会重新生成并比对，换参考或换输入后旧审核失效。先完成配置，再生成审核草稿；修改配置哈希后必须生成新的审核记录。

## 下游计划与审核停点

模板自身 production run 可尝试严格适配：

```bash
python3 scripts/import_nfcore.py --config config/project.json \
  --run results/runs/production_01 --out results/import_01
```

适配器只接受 nf-core/cutandrun 3.2.2、Bowtie2、MACS2 narrow 或 broad、target duplicate-retained 的已完成生产输出；逐项核对配置/样本快照，缺失产物不猜测替代文件。旧项目或不同布局通过 `examples/artifacts.tsv` 显式登记：每生物样本一个已经合并的 BAM 和一个选定峰文件。`bam_policy=duplicates_retained` 是需审核的声明，不能恢复已经去掉的重复片段。内容相同的 BAM 不允许登记为两个统计样本。

```bash
pixi run -e analysis python scripts/plan_downstream.py --config config/project.json \
  --artifacts results/import_01/artifacts.tsv --out results/downstream_01
pixi run -e analysis python scripts/workflow.py --workflow results/downstream_01/workflow.json \
  --out results/execution_01
# 中断后校验并跳过内容完全相同的成功步骤
pixi run -e analysis python scripts/workflow.py --workflow results/downstream_01/workflow.json \
  --out results/execution_01 --resume
```

自动计划包含参考审核、每样本 paired fragments、共识峰、原始计数、同一 master peak universe 的 QC、CPM BigWig/IGV。单次计划支持 `--mapq`、`--fraction` 和 `--remove-duplicates`。不会自动进入正式差异分析：先审阅 QC/峰/重复独立性，再执行 prepare_review.py 和 run_differential.py。Spike-in 轨迹仍需独立的系数审核，不从试跑结果自动套用。

### 生成并比较下游方案矩阵

先为每种已运行的 peak calling 方案准备一个 artifact TSV；每份表必须包含完全相同的生物样本 ID，且目标 BAM 保留 duplicate flags。`examples/strategy_matrix.json` 演示如何声明峰集及下游候选轴。候选数量是各轴组合数乘以 artifact set 数，默认最多 48 项，避免意外展开过大的网格。

2/3 支持度请写 `0.6666666666666666`，不要四舍五入为 `0.6667`：支持样本数按向上取整计算，后者在三个重复时会要求 3/3。自动计划会核对 BAM 序列字典和 peak 坐标与 target FASTA 索引一致。已有 v0.5.0/v0.5.1 矩阵若使用过 `0.6667`，请修正配置并用新目录重跑共识及下游；旧产物保留为历史结果。

```bash
pixi run -e analysis python scripts/strategy_matrix.py \
  --config config/project.json --matrix examples/strategy_matrix.json \
  --out results/strategy_matrix_plan_01
# 审核 strategy_matrix.json 和各候选 workflow.json 后再执行
pixi run -e analysis python scripts/strategy_matrix.py \
  --config config/project.json --matrix examples/strategy_matrix.json \
  --out results/strategy_matrix_run_01 --execute
pixi run -e analysis python scripts/compare_strategy_matrix.py \
  --matrix results/strategy_matrix_run_01/strategy_matrix.json \
  --out results/strategy_comparison_01
```

矩阵执行中断且配置、矩阵、artifact 表均未变化时，在同一输出目录加 `--resume --execute` 继续；成功的候选按哈希跳过，已失败的候选仍保留失败证据，修复原因后请用新输出目录重跑。

比较报告在完成候选峰的 pooled union 上重新计算 FRiP，并在相同区间上重计数以计算 common-universe replicate correlation。它另列每种方案自己的 peak-universe correlation（仅用于本方案内部 QC，不用于候选排序）、峰数、宽度、覆盖度和 master intervals 的两两覆盖一致性。失败候选会保留状态；未完成的候选不会伪装成零值。

正式差异模型仍需逐个完成审核和运行。将审核过的 DE complete tables 登记到 `candidate_id / normalization / result_tsv / review_json` 表，再运行 `compare_models.py`；它要求模型使用相同对比、设计和过滤阈值，并报告坐标匹配率、log2FC 相关性及方向/显著性变化。互相重叠的 peak universe 使用一对一双向覆盖匹配。此报告不选择归一化或 peak caller；spike-in 是否可解释仍需结合加样证据和校准范围。

续跑会校验源代码、配置、锁文件、声明输入、依赖、解释器和成功产物的哈希。变化时要求新建计划与产物目录。失败步骤若留下任何声明产物，不覆盖重跑；保留失败证据，用新目录运行。没有自动删除原始数据或清理工作目录的功能。

## 自动 QC 与可选方案

`collect_qc.py --manifest samples.tsv --out results/qc_metrics_01` 的输入列为 `sample_id, fragments_bed, peaks_bed`（TSV）。FRiP 分母是所有保留的 paired fragments；与多个峰重叠的同一片段只计一次；重复坐标行按原始片段数保留。不同样本比较应使用同一个 peak universe 和过滤规则。不能从 BED 推断 PCR duplication；不能从 target fragments 推断 spike-in fraction，这两项明确为 NA。

可选方案的命令、真实数据试验与推荐见 [OPTION_VALIDATION.md](OPTION_VALIDATION.md)。

## 富集和 motif 宇宙

```bash
python3 scripts/prepare_universe.py --results results/model/complete.tsv \
  --annotation results/annotation.tsv --out results/universe_01 --alpha 0.05 --lfc 1
```

必须提供完整 DE 表和带稳定 peak_id、BED 坐标、geneId、annotation 的峰注释。只将有限 padj 和 log2FC 的峰纳入可判定显著性的宇宙；过滤/未校正的峰不会混入 ORA 背景。输出 all/promoter/distal 三套去重 Entrez gene lists，以及 gain/loss 的 motif target/background BED。一个基因的不同峰可以同时 gain 和 loss，单独记录，不强制归成单一方向。远端注释是位置关联，不是已验证 enhancer–gene 链接。每个数据库仍需查看 `.mapping.tsv` 和 `.audit.txt` 中的有效背景。

## 可搬迁交付

release JSON 中每项可指定 `path`（文件或目录）、`destination`（包内相对路径）、`required`。整体选择 QC/轨迹目录可保留 HTML 图片和 IGV 资源。目录中有软链接则要求显式选取已解析文件，防止无意递归收集外部数据。

```bash
python3 scripts/finalize.py --manifest release.json --out results/delivery_01
python3 scripts/verify_release.py --package results/delivery_01
# 搬到其他路径后重新运行 verify_release；不依赖原始 source 路径
```

打包会复制文件、重写已选 IGV 资源相对路径、生成可点击入口和逐文件哈希；缺必需文件、缺 IGV 资源、目标路径冲突会失败。`state=REVIEW_REQUIRED` 保留科学交付审核，不因文件齐全而自动批准结论。
