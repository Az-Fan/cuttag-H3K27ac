# CUT&Tag 项目模板 v0.3.1

**核查结论：仍是模块化开发模板，尚非完整验收的稳定全流程。** 详细问题、修复和剩余缺口见 [核查报告](docs/AUDIT.md)。

从 ACLY-CUTTAG 两批项目沉淀的独立模板。当前版本提供样本校验、完整 FASTQ 检查、nf-core 运行计划/执行、共识峰、spike-in 计数审计、MACS2 诊断矩阵、片段计数、DESeq2、峰注释、离线 ORA、HOMER motif、统计可视化和交付审阅入口。完整建设规划见 [docs/BUILD_PLAN.md](docs/BUILD_PLAN.md)。

## 快速开始

复制整个模板目录作为新项目，编辑 `config/project.json` 与 `config/samples.tsv`，登记参考和原始数据。路径相对于项目根目录，也支持绝对路径。示例是每组一个样本，仅供结构演示；真实数据和参考需自行配置。

```bash
python3 scripts/cuttag.py validate --allow-missing
python3 scripts/cuttag.py plan --stage alignment
# 安装依赖并生成项目自己的锁文件
pixi install -e default
# 分析、轨迹和 QC 图表模块
pixi install -e analysis
pixi run python scripts/cuttag.py test
# 将成功的 test run ID 写入 qc.official_test_run_id
pixi run python scripts/cuttag.py validate
pixi run python scripts/cuttag.py inventory
pixi run python scripts/cuttag.py run --stage alignment
# 审阅 MultiQC、片段长度、重复、参考/spike-in 后更新 qc 配置
pixi run python scripts/cuttag.py run --stage production
```

`plan` 只保存计划；允许未准备好的输入。`run` 检查文件存在和磁盘资源，生产阶段要求结构化 QC/对照策略审核记录并绑定当前配置与样本表哈希；spike-in 还要求已知加入阶段、等量前提、解释范围与接受的计数方法。成功运行只标记 `COMPUTATIONAL_PASS`，数据验收须另行记录。官方 test 的成功需要在决策记录中确认。

## 样本登记

每行是一对测序 FASTQ。一个生物样本拆成多份测序文件时，沿用相同 `sample_id`、`biological_sample_id`、`group` 和 `replicate`，改变 `unit_id` 和 FASTQ 路径。不同文库使用不同 `library_id`。nf-core 的相同 group/replicate 输入会归到同一重复；最终产物名称以实际运行输出为准。

`sample_id` 是统计样本单位；独立生物样本不能共享 group/replicate。文件名、相同 index 或高相关性都不能代替建库记录确认。每组少于两个独立生物样本时，正式差异分析入口不可用。

## 分析模块

参见 [docs/SOP.md](docs/SOP.md) 中的文件契约与命令。各模块可独立运行，须用明确的已验收输入。默认实例为 H3K27ac narrow peaks，其他靶标须修改策略。

- `scripts/cuttag.py`：配置与设计校验、FASTQ 完整性/配对/哈希、nf-core 官方 test、运行计划和生产日志。
- `scripts/fragments.py`、`scripts/tracks.py`：从 BAM 构建 paired fragment 与 CPM/Spike-in BigWig。`scripts/count_fragments.py`：已验收 BED4 峰集 × 每生物样本 fragment BED 的原始计数。
- `scripts/run_differential.py`：从项目配置传递模型参数、检查峰集接受和样本映射；`scripts/differential.R` 也独立核对模型审核记录，避免直接调用时丢失必要判断。
- `scripts/consensus.py` 与 `scripts/peak_diagnostics.py`：以独立样本为单位构建共识峰并诊断 MACS2 背景策略。`scripts/qc_atlas.py` 汇总 QC，`scripts/workflow.py` 按显式依赖顺序执行命令。`scripts/annotate.R`：显式 TxDb SQLite、OrgDb 和 BED 坐标转换。
- `scripts/enrichment.R`：使用明确背景和离线 TERM2GENE 的 ORA。

## 当前验证范围

26 项 Python 行为测试已通过，所有 R 脚本通过语法检查。Pixi 依赖已成功解析并生成 `pixi.lock`；系统 R 环境中的合成 DESeq2 和两个源项目的设计登记回放见 [PROJECT_STATUS.md](PROJECT_STATUS.md)。

本版本仍是开发版。真实 FASTQ 上的 nf-core 官方/生产流程、HOMER 和完整 ORA 工作流尚待验证；本轮已对本地 TxDb 的峰 ID/坐标和微型 BAM→BigWig 做实际验证，范围详见核查报告。

新增模块用法详见 [SOP](docs/SOP.md)：`consensus.py`、`spikein_audit.py`、`peak_diagnostics.py`、`motif.py`、`visualize.R`、`finalize.py`。命令式脚本相互独立，由明确输入契约衔接；目前没有统一自动执行所有下游步骤的编排器。spike-in 全量比对需用审计选定参数运行 Bowtie2，再交日志给审计入口。

paired fragment 生成、CPM/可选 spike-in BigWig 与 QC atlas 模块也已加入；可用 `scripts/workflow.py` 将明确的命令按依赖顺序串联。
