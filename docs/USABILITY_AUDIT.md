# v0.5.2 使用前核查（2026-10-05）

结论：在模板声明范围内，可以用于新项目。范围为 Linux/Pixi、paired-end H3K27ac、MACS2 narrow/broad，以及具有独立生物重复的两组 `~ condition` 模型。新项目仍需填写真实输入、准备兼容参考，并按 SOP 审核 QC、对照、重复独立性和校准证据。

## 本次发现并修复

| 问题 | 实际影响 | 修复与验证 |
|---|---|---|
| 示例矩阵将 2/3 写成 `0.6667` | `ceil(3 × 0.6667)` 为 3，三重复中仅两重复支持的区域被错误排除于所称的 2/3 基线 | 改为 `0.6666666666666666`；用三个输入 BED 实跑，验证支持数为 2、预期区间保留 |
| 新增 common-universe 计数测试依赖 bedtools，CI 未安装 | 系统 Python 下测试失败，Pixi 本机通过不能保证 CI 可运行 | CI 显式安装 bedtools；在系统 Python 加单独 bedtools 的精简环境下跑完 58 项测试 |
| 比较器允许缺失执行记录时使用矩阵中的 PASS 标签 | 不完整执行证据可能被当成已完成结果 | 缺少 workflow 状态、步骤或产物哈希时拒绝；校验候选间样本/生物身份/条件映射一致；测试删除记录和改变条件 |
| 外部 target BAM 和 peak 缺少与参考的强制绑定 | 错误参考的 BAM/peak 可能进入下游计数 | 自动计划强制传入 FASTA 索引；核对 BAM 序列名/长度字典和 peak 坐标边界；不匹配在计数前失败 |
| DE 区域匹配反复扫描同染色体整个前缀 | 大峰集比较可能很慢 | 用起点和累计最大终点缩小扫描；验证长包围区间不会遗漏；50,000 对相同区间全部匹配，本机约 0.18 秒（非性能保证） |

同时修正 common-universe 测试数据，使候选 count 矩阵与 fragment 文件一致；对固定区间重新计数后，检查每样本准确计数及相关性，而非只检查输出文件存在。汇总 TSV 补充 `min_width`，版本信息同步为 0.5.2。

2/3 修复影响从旧示例复制的策略矩阵。`consensus.py` 原本的默认 `2/3` 不受此问题影响。已按旧矩阵分析的项目应新建目录重跑受影响的共识和下游，不能凭此直接认定此前两个项目的结果全部错误。

## 核查保留的关键行为

- 同一生物样本的技术 unit 不增加统计重复数；1 vs 1 无正式 DE 入口。
- 目标计数为 proper primary paired fragments，双端 MAPQ 过滤、重复策略显式声明、没有 ATAC shift。
- 共识采用明确的 support-core / reproducible-union 定义，blacklist 默认 subtract，最小宽度 50 bp，并保留 1/50/100 bp 敏感性。
- 固定 pooled-union FRiP 与相关性用于跨候选比较；候选自身相关性保留为内部 QC。过滤规则相同且只改变 peak universe 时，common-universe 相关性相同是预期行为，它不能独立选择 peak 边界。
- DESeq2 接收整数原始计数；校准按样本 spike-in 计数设置 size factors，CPM 轨迹与 DE 分开。不同校准组不能直接混入当前简单模型。
- ORA/motif 背景使用明确的检验宇宙；最近基因保留为位置关联；旧审核和有变化的输入不能直接续跑复用。

计数语义对照 [bedtools coverage 文档](https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html)，原始计数与 size factors 对照 [DESeq2 官方说明](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)。实现正确性主要由本地真实工具 fixture 和边界测试检验。

## 使用边界

复杂批次、配对或交互设计需要另行实现；模板会拒绝静默替换成简单模型。全量真实 HPAEC 数据的 lambda calibration、IgG 有效深度和方案优劣仍需实际数据验证。当前核查与合成测试不能替代这些实验判断，也没有将历史 nf-core 官方测试重新标记为本次运行。

操作入口见 [OPERATIONS](OPERATIONS.md)，本版测试证据见 [验收记录](validation/2026-10-05/local_acceptance_v052.json)。
