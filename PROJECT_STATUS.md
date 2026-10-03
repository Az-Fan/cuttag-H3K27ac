# v0.4.0 验收状态

本版已形成可复用的项目操作闭环：初始化 → 环境与参考检查 → nf-core 上游 → 审核 → 自动下游计划 → 模型/注释/富集/motif → 可搬迁交付。采用分阶段审核，不会为追求一键运行而自动接受生物重复、IgG、spike-in 或统计模型。

## 已补齐

- 干净项目初始化，模板来源文件哈希；科学批准不继承。
- FASTA/FAI 字典、GTF/BED 坐标范围和参考哈希审核；生产执行重新验证。
- nf-core 3.2.2 特定输出布局的严格适配；不同布局使用显式产物表；拒绝样本映射变化和不同统计样本复用同一 BAM 内容。
- fragments → 共识峰 → counts → 自动 QC → CPM tracks 的可执行计划；明确的模型审核停点。
- 哈希一致才续跑、进程锁、失败日志、缺失产物阻断；Nextflow 忽略失败单独暴露。
- FRiP 按 paired fragment 计数；同一片段跨多个峰不重复计数；不可推断的 QC 项记 NA。
- 同输入 spike-in 4 模式比较和按校准组分别计算相对因子；峰方案区间覆盖比较。
- tested gene/motif universe、promoter/distal 分层、双向基因和 ORA 映射损失报告。
- MACS2 固定容器、Bioconductor 数据依赖安装/校验、运行时 doctor、可重复 acceptance 与契约 CI。
- 实际复制交付文件、相对 IGV 路径、可点击入口、搬迁后的哈希验证。

## 已有验收证据

- 44 项 Python 契约/回归测试通过，包括技术拆分、审核失效、参考不兼容、错误坐标、缺失文件、续跑污染、适配器样本映射、校准组隔离和共识支持度网格。
- 最终本机集成套件 8/8 通过（含 44 项契约测试），在 Pixi 分析环境运行：真实微型 BAM/fragment/count/BigWig；干净新项目 DAG 与搬迁交付；4 组 MAPQ/重复设置；3 种 MACS2 控制策略；2 种 DESeq2 归一化、ORA、显式 TxDb/OrgDb 注释、HOMER；无变化信号/无有效 padj 的绘图降级也通过。
- 两批原项目各 2 个已剪切技术单元，4 × 20,000 pairs × 4 种 Bowtie2 设置实际测试；定量结果与条件性推荐见 [OPTION_VALIDATION](docs/OPTION_VALIDATION.md)。未把第二批拆分文件当生物重复。
- 官方 nf-core/cutandrun 3.2.2 `test,apptainer` 运行结束：退出码 0；133 个任务完成，2 个 Preseq 失败被上游忽略。任务级状态单独保存，不能表述成“所有任务无错误”。
- 模板默认 MACS2 narrow、保留 target 重复、CPM 参数组合也已在官方微型输入上独立跑完，退出码 0；具体任务警告及峰产物数见 [模板参数测试](docs/validation/2026-10-03/template_params_test.json)。
- 原生 MACS2 二进制故障已实际复现；默认改为已验证且固定 SHA256 的 2.2.7.1 容器。Pixi 环境缺失的 GenomeInfoDbData/GO.db/TxDb 数据包已按注册表校验补齐；人类 OrgDb 已加入锁文件和安装步骤。

原始集成验收记录和后续补充验证见 [docs/validation/2026-10-03](docs/validation/2026-10-03)。只保留小型验证摘要进入 Git；原始 reads、BAM、数据库、缓存和完整运行结果不入库。

## 适用边界与仍需完成的实验验收

1. 当前默认支持 paired-end H3K27ac、MACS2 narrow、两组独立样本 `~ condition`。批次、配对、复杂交互设计及其他靶标/调用器需单独实现或配置并验收；不会静默简化模型。
2. 第一批独立生物重复仍须实验记录确认；第二批已确认 1 vs 1，只能描述性分析。模板不替用户批准这些事实。
3. spike-in 子集结果支持 overlap 候选，但仍需合并技术单元后的全量计数、双端质量/多重比对策略，以及等量、加入阶段和校准范围证据。
4. IgG 策略、去重、MAPQ、共识支持度和归一化在具体项目上的优劣，需要真实数据的位点检查、重复一致性与敏感性比较。软件验收不能代替该判断。
5. 显式输入适配器的布局契约已有测试；跨所有 nf-core 参数/版本的自动识别、原两批全量重新分析、de novo motif/GSEA、多个峰调用器与全部组合的性能比较，不在本版已验收声明内。
6. 交付包仅验证文件及链接资源，仍标记 REVIEW_REQUIRED。没有执行破坏性清理，也不会自动删除原始数据。

运行步骤：[OPERATIONS](docs/OPERATIONS.md)；方法与审核：[SOP](docs/SOP.md)；历史缺陷审计：[AUDIT](docs/AUDIT.md)。
