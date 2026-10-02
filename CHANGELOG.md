# 变更记录

## 0.3.1 — 必要判断与数据契约核查修复

- 修复 spike-in filtering 分支与审核后比对配置未连接的问题。
- 正式生产与差异分析需要内容绑定的审核记录；直接调用 differential.R 的最后一个参数现在必须是 review.json。新增 prepare_review.py 和配置化 run_differential.py。旧调用方式需按 SOP 更新。
- 加强样本/模型/数值/产物校验，修正峰 ID 传递、片段路径解析及 motif 边界处理。
- 26 项 Python 测试及多项真实工具小型集成验证通过；完整性结论与未完成项见 docs/AUDIT.md，仍为开发版。
