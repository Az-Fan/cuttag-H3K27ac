# 分析日志

模板实例化后，每个阶段追加：做了什么、目的、输入/输出、软件与参数、QC 与异常、决策证据、下一步。实际运行配置快照与日志存 results/runs/<run_id>。

## 模板 v0.2.0 构建

新增共识峰、spike-in 配对计数审计、MACS2 背景敏感性矩阵、HOMER motif、统计可视化和交付索引。生成 Pixi 锁文件并修正 Python 3.12/MACS2 兼容性，采用 Python 3.11。官方 test 成功 run ID 现在是实际执行的门禁。合成 DESeq2/可视化与两批设计回放通过；环境安装因镜像网络中断尚待完成。

## 模板 v0.3.0 补全

加入 BAM paired fragment 生成、fragment CPM/Spike-in BigWig 与相对 IGV 会话、QC atlas 及显式依赖 workflow runner。对片段流程改为按名称排序后流式解析，不将全 BAM SAM 文本载入内存。11 个 Python 契约/模块测试通过，合成 DESeq2 和可视化、R/Python 语法检查及源设计只读回放通过。Pixi 环境安装遇主机 mirror 与官方频道 TLS EOF，故 nf-core 官方 test 尚未运行。
