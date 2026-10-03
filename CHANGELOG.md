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
