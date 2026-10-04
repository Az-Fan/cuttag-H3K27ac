# 分析日志

模板实例化后，每个阶段追加：做了什么、目的、输入/输出、软件与参数、QC 与异常、决策证据、下一步。实际运行配置快照与日志存 results/runs/<run_id>。

## 模板 v0.2.0 构建

新增共识峰、spike-in 配对计数审计、MACS2 背景敏感性矩阵、HOMER motif、统计可视化和交付索引。生成 Pixi 锁文件并修正 Python 3.12/MACS2 兼容性，采用 Python 3.11。官方 test 成功 run ID 现在是实际执行的门禁。合成 DESeq2/可视化与两批设计回放通过；环境安装因镜像网络中断尚待完成。

## 模板 v0.3.0 补全

加入 BAM paired fragment 生成、fragment CPM/Spike-in BigWig 与相对 IGV 会话、QC atlas 及显式依赖 workflow runner。对片段流程改为按名称排序后流式解析，不将全 BAM SAM 文本载入内存。11 个 Python 契约/模块测试通过，合成 DESeq2 和可视化、R/Python 语法检查及源设计只读回放通过。Pixi 环境安装遇主机 mirror 与官方频道 TLS EOF，故 nf-core 官方 test 尚未运行。

## 完整性与判断核查

核查 Git 46e2d24 后发现并修复多处执行门禁和数据契约缺陷，详细证据、13 项问题及未完成清单写入 docs/AUDIT.md。新增绑定输入/证据哈希的 review、配置化差异入口和 reviewer 草稿生成器；直接 R 也核对 review。26 项 Python 测试、微型真实工具流程、合成模型和本地 ChIPseeker 注释实测通过。现状明确为开发模板，尚未达到完整稳定全流程标准。

## 2026-10-04 — v0.4.1 spike-in/peak sensitivity

- Replaced whole-consensus-interval blacklist removal with configurable basewise subtraction (default) and retained drop mode only for sensitivity.
- Added paired lambda fragment counting from spike-in-only BAMs with pair-level MAPQ20/30, primary/proper-pair checks, reference dictionary verification and biological-sample aggregation. Bowtie2 summary audit now reports unique and multi pairs separately and labels all size factors diagnostic-only.
- Added competitive target+spike-in mapping audit for matched FASTQ subsets; added synthetic integration fixture.
- Expanded MACS2 comparison to narrow/broad crossed with default IgG scaling, scale-to-large and no-IgG diagnostic; supports narrow/broad import and interval-width summaries.
- Enforced consistent spike-in scaling across all tracks in a manifest.
- Validation: 44 unittest contracts pass; spike-in fixture pass; six MACS2 synthetic options pass with planted loci recovered; template `validate --allow-missing` has no errors and correctly blocks formal inference for example 1-vs-1 design.
- Real source evidence: Batch2 existing comparison supports narrow primary + broad sensitivity and shows IgG processing materially changes calls. Full paired-fragment MAPQ re-count from all source BAMs was not completed in this session; Batch1 has only Bowtie2 summary logs, not lambda BAMs. No scientific spike-in calibration is accepted.

## 2026-10-03 — v0.4.0 项目闭环与可选方案验收

补齐初始化、参考检查、严格上游适配、下游 DAG、内容一致续跑、自动片段 QC、tested universe、可搬迁交付及环境验收。保留原实验生物重复、IgG、spike-in 与模型审核门禁。新增参数对照执行器与峰覆盖比较，记录四个真实技术单元的 spike-in 子集结果；未执行两批全量重分析，未自动接受校准。

官方 nf-core 3.2.2 test 已执行，退出码 0，但 trace 有两个被忽略的 Preseq 失败。新增 task audit 公开这些失败。原生 Bioconda MACS2 存在 __log_finite 运行时错误，已从默认依赖移除，以固定哈希的 MACS2 2.2.7.1 容器替代，并通过三种背景策略的已知信号测试。补装并验证 Bioconductor 数据依赖，避免系统 R 的成功掩盖 Pixi R 问题。

小型证据保存在 docs/validation/2026-10-03，完整本机日志和产物保存在 results/；当前实验方案推荐及范围见 docs/OPTION_VALIDATION.md。
