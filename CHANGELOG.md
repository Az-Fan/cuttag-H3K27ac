# v0.5.2 — 新数据入口与方案矩阵核查

- 修正示例矩阵支持度 `0.6667` 导致三重复实际要求 3/3 的问题，使用 `0.6666666666666666` 表示 2/3；加入真实三样本共识行为测试。
- 自动下游计划强制校验 target BAM 序列字典和输入 peak 坐标与目标参考一致。
- 方案比较拒绝缺失执行记录、缺失产物哈希或样本身份/条件不一致；汇总表显式报告 min_width。
- common-universe 回归测试使用与片段文件一致的原始计数并检验准确数值；CI 安装所需 bedtools。
- DE 区域匹配跳过已结束区间，避免常见不重叠峰集每次扫描整条染色体前缀。

# v0.5.1 — 宽度敏感性与统一评价区间

- 共识峰最小宽度默认 50 bp，示例比较 1/50/100 bp；统一评价区间重新计数并计算重复相关性。
- crossmap 增加 biological_sample_id/unit_id；清除重复的源代码哈希函数与字典 key。

# v0.5.0 — 显式候选矩阵

- 多个 peak artifact sets 与下游过滤/共识候选组成独立计划，可执行、续跑、比较，并保留科学审核停点。

# v0.4.1 — spike-in 与峰策略审计完善

- 默认逐碱基扣除 blacklist；支持对照用的 whole-interval drop，并显式提供 support-core / reproducible-union 共识 universe。
- 新增 spike-in-only paired-fragment 计数器：proper primary pairs、双端 MAPQ20/30、参考字典校验、技术 unit 汇总到 biological sample；加入 target+spike-in 竞争比对审计。
- Bowtie2 summary 审计拆分唯一/多重比对，并禁止将日志估算误作为正式 DE/轨迹因子。
- MACS2 诊断扩展为 narrow/broad × 3 种 IgG 背景；nf-core adapter 支持对应 broad/narrow 产物。
- 轨迹拒绝同一 manifest 混用 spike-in 与 CPM；补充合成端到端 fixture、SOP 与策略验证说明。

# v0.4.0 — 可重复验收与项目操作闭环

- 干净初始化，nf-core 3.2.2 严格结果适配，显式 BAM/peak 登记和下游 DAG 计划。
- 内容哈希续跑、锁定执行、缺失产物阻断；参考字典/坐标校验及生产审核绑定。
- 从 paired fragments 自动收集 QC；tested gene/motif universe 和 ORA 映射损失报告。
- 固定输入的 4 模式 spike-in 执行器、按 calibration_group 分别计算系数、峰方案覆盖比较。
- 固定哈希 MACS2 容器替换本机失败的原生依赖；补全 Bioconductor data packages，增加运行时 doctor 和集成 acceptance。
- 可复制交付包、IGV 相对路径修复、搬迁后文件校验；GitHub 契约 CI。
- 本机官方 nf-core test 退出码 0，133 个任务完成、2 个 Preseq 失败被上游忽略；明确记录任务级警告。
- 必要科学判断保持人工审核，不自动增加生物重复、接受 spike-in 或挑选显著峰最多的模型。

# 变更记录

## 0.3.1 — 必要判断与数据契约核查修复

- 修复 spike-in filtering 分支与审核后比对配置未连接的问题。
- 正式生产与差异分析需要内容绑定的审核记录；直接调用 differential.R 的最后一个参数现在必须是 review.json。新增 prepare_review.py 和配置化 run_differential.py。旧调用方式需按 SOP 更新。
- 加强样本/模型/数值/产物校验，修正峰 ID 传递、片段路径解析及 motif 边界处理。
- 26 项 Python 测试及多项真实工具小型集成验证通过；完整性结论与未完成项见 docs/AUDIT.md，仍为开发版。
