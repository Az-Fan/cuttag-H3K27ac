# 模板开发状态

版本：0.3.0，开发版。

已实现：按 proper paired BAM 生成片段 BED/片段长度；fragment CPM 与可选已审 Spike-in BigWig/相对 IGV session；QC PDF/HTML atlas；显式依赖步骤的顺序编排器。

已实现：配置与样本层级校验、FASTQ 完整性/哈希、nf-core 官方 test/生产入口、按生物样本的共识峰、原始片段计数、spike-in 日志审计、MACS2 三策略诊断、两组 DESeq2、ChIPseeker 注释、离线 ORA、HOMER 已知 motif、PCA/相关性/火山图/热图、交付 HTML 审阅与校验和清单。

已验证：11 项 Python 行为测试；所有 Python/R 入口语法检查；使用系统已有 DESeq2 在 300 峰、6 个独立合成样本上拟合和可视化，差异方向通过；两批源项目设计登记只读回放通过。第一批 3 vs 3 仍未核实重复来源，第二批 1 vs 1 正确阻止正式推断。

环境：pixi.lock 成功生成。安装失败于主机 Pixi 配置的 Westlake mirror TLS/隧道连接；绕过镜像后官方频道也遇到 TLS 中断。锁定环境中的自测与官方 nf-core test 尚未运行。

尚待验证：真实 FASTQ 上的官方/生产流程、外部 HOMER、MACS2 诊断、参考注释和 ORA 端到端、BigWig/IGV 可迁移性。

尚待完善：spike-in 比对参数矩阵自动执行/全量计数；完整端到端真实项目回放；参考/数据库自动校验与可搬迁交付包；受保护存储清理。统一步骤编排器已支持显式依赖，但仍需项目级 workflow manifest。当前仍不能标为已验收稳定版。
