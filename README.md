# CUT&Tag 项目模板 v0.1.0

从 ACLY-CUTTAG 两批项目沉淀的独立模板。当前版本提供样本校验、完整 FASTQ 检查、nf-core 运行计划/执行、片段计数、DESeq2、峰注释和离线 ORA 入口。完整建设规划见 [docs/BUILD_PLAN.md](docs/BUILD_PLAN.md)。

## 快速开始

复制整个模板目录作为新项目，编辑 `config/project.json` 与 `config/samples.tsv`，登记参考和原始数据。路径相对于项目根目录，也支持绝对路径。示例是每组一个样本，仅供结构演示；真实数据和参考需自行配置。

```bash
python3 scripts/cuttag.py validate --allow-missing
python3 scripts/cuttag.py plan --stage alignment
# 安装依赖并生成项目自己的锁文件
pixi install
pixi run python scripts/cuttag.py test
pixi run python scripts/cuttag.py validate
pixi run python scripts/cuttag.py inventory
pixi run python scripts/cuttag.py run --stage alignment
# 审阅 MultiQC、片段长度、重复、参考/spike-in 后更新 qc 配置
pixi run python scripts/cuttag.py run --stage production
```

`plan` 只保存计划；允许未准备好的输入。`run` 检查文件存在和磁盘资源，生产阶段要求 `qc.upstream_accepted=true`；spike-in 生产还要求 `spikein.calibration_accepted=true`。成功运行只标记 `COMPUTATIONAL_PASS`，数据验收须另行记录。官方 test 的成功需要在决策记录中确认。

## 样本登记

每行是一对测序 FASTQ。一个生物样本拆成多份测序文件时，沿用相同 `sample_id`、`biological_sample_id`、`group` 和 `replicate`，改变 `unit_id` 和 FASTQ 路径。不同文库使用不同 `library_id`。nf-core 的相同 group/replicate 输入会归到同一重复；最终产物名称以实际运行输出为准。

`sample_id` 是统计样本单位；独立生物样本不能共享 group/replicate。文件名、相同 index 或高相关性都不能代替建库记录确认。每组少于两个独立生物样本时，正式差异分析入口不可用。

## 分析模块

参见 [docs/SOP.md](docs/SOP.md) 中的文件契约与命令。各模块可独立运行，须用明确的已验收输入。默认实例为 H3K27ac narrow peaks，其他靶标须修改策略。

- `scripts/cuttag.py`：配置与设计校验、FASTQ 完整性/配对/哈希、nf-core 运行计划和日志。
- `scripts/count_fragments.py`：已验收 BED4 峰集 × 每生物样本 fragment BED 的原始计数。
- `scripts/differential.R`：独立重复核验、常规或 spike-in size factors、完整结果/模型/MA 图。
- `scripts/annotate.R`：显式 TxDb SQLite、OrgDb 和 BED 坐标转换。
- `scripts/enrichment.R`：使用明确背景和离线 TERM2GENE 的 ORA。

## 当前验证范围

Python 契约测试覆盖技术拆分、重复映射冲突、FASTQ 重用、单样本阻止正式推断。未执行真实测序数据、官方 nf-core test、R 分析或环境安装；此版本仍是开发版，尚不能称稳定版。`pixi.lock` 将在首次安装/解析后生成并需纳入新项目版本管理。

待后续实现：自动 consensus/重复支持峰构建、spike-in 参数审计与全量重计数、caller 诊断矩阵、通用 motif/轨迹图、自动最终报告与清理入口、迁移回放与端到端验证。相关手工方法和验收规则见 SOP 与建设规划。
