# v0.5.2 开发状态

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
- blacklist 默认逐碱基扣除，另有 whole-interval-drop 敏感性开关；consensus 支持 reproducible-union 与 support-core 两种可比较定义。
- spike-in-only BAM 正式 fragment counter：proper primary pair、双端 MAPQ20/30、参考字典约束、技术 unit 到 biological sample 汇总和按校准组系数；另有竞争参考跨映射审计。Bowtie2 summary 明确只属诊断。
- peak caller sensitivity 覆盖 narrow/broad × 默认 IgG / scale-to-large / no-IgG × keep-dup all/auto，共 12 组合；轨迹输入禁止部分样本 spike-in、部分 CPM。
- 新增方案矩阵：多个 caller/background artifact sets 可与 MAPQ、duplicate、consensus fraction、peak universe 组合；每个候选保存独立 DAG、状态和输入哈希，矩阵可在配置/输入不变时续跑。
- 新增固定 pooled candidate peak union 的 FRiP、重复样本 raw-count correlation、master interval overlap 对照报告；新增经过审核的 DE complete tables 的 reciprocal-overlap/LFC/direction 对照入口。
- 共识峰 blacklist subtraction 后默认移除短于 50 bp 的残片；策略矩阵纳入 1/50/100 bp 的完整敏感性轴，50 bp 仅为可复核的模板基线而非普适生物学阈值。策略报告对固定 pooled-union 区间重计数并给出可跨策略比较的 common-universe replicate correlation，同时保留候选自身峰集相关性作为内部 QC。
- spike-in competitive cross-map manifest 与正式计数统一 `sample_id / biological_sample_id / unit_id` 身份，可分别报告技术 unit 并追溯到统计样本；删除策略矩阵中重复的源码哈希函数、校验和 JSON key。
- Spike-in formal manifest 强制填写 calibration_group；目标组内仅 target 生物样本参与 size factor，单独输出 target biological count 表；QNAME 分组使用 samtools collate 流式处理，未比对 pair 纳入竞争比对分类。
- MACS2 诊断允许仅 no-IgG 运行，并支持显式选择 shape/background/`keep-dup all|auto` 网格。

## v0.4.1 历史验收证据

- 44 项 Python 契约/回归测试通过，包括技术拆分、审核失效、参考不兼容、错误坐标、缺失文件、续跑污染、适配器样本映射、校准组隔离和共识支持度网格。
- 最终本机集成套件 8/8 通过（含 44 项契约测试），在 Pixi 分析环境运行：真实微型 BAM/fragment/count/BigWig；干净新项目 DAG 与搬迁交付；4 组 MAPQ/重复设置；3 种 MACS2 控制策略；2 种 DESeq2 归一化、ORA、显式 TxDb/OrgDb 注释、HOMER；无变化信号/无有效 padj 的绘图降级也通过。
- 两批原项目各 2 个已剪切技术单元，4 × 20,000 pairs × 4 种 Bowtie2 设置实际测试；定量结果与条件性推荐见 [OPTION_VALIDATION](docs/OPTION_VALIDATION.md)。未把第二批拆分文件当生物重复。
- 官方 nf-core/cutandrun 3.2.2 `test,apptainer` 运行结束：退出码 0；133 个任务完成，2 个 Preseq 失败被上游忽略。任务级状态单独保存，不能表述成“所有任务无错误”。
- 模板默认 MACS2 narrow、保留 target 重复、CPM 参数组合也已在官方微型输入上独立跑完，退出码 0；具体任务警告及峰产物数见 [模板参数测试](docs/validation/2026-10-03/template_params_test.json)。
- 原生 MACS2 二进制故障已实际复现；默认改为已验证且固定 SHA256 的 2.2.7.1 容器。Pixi 环境缺失的 GenomeInfoDbData/GO.db/TxDb 数据包已按注册表校验补齐；人类 OrgDb 已加入锁文件和安装步骤。
- v0.4.1：46 项 Python 契约/回归测试、spike-in BAM/竞争映射合成 fixture、Python 编译和 `cuttag.py validate --allow-missing` 通过；六组 MACS2 narrow/broad × 背景合成比较均找回植入富集位点。`acceptance.py` 已纳入 spike-in fixture。

原始集成验收记录和后续补充验证见 [docs/validation/2026-10-03](docs/validation/2026-10-03)。只保留小型验证摘要进入 Git；原始 reads、BAM、数据库、缓存和完整运行结果不入库。

## 适用边界与仍需完成的实验验收

1. 当前默认支持 paired-end H3K27ac、MACS2 narrow、两组独立样本 `~ condition`。批次、配对、复杂交互设计及其他靶标/调用器需单独实现或配置并验收；不会静默简化模型。
2. 第一批独立生物重复仍须实验记录确认；第二批已确认 1 vs 1，只能描述性分析。模板不替用户批准这些事实。
3. spike-in 子集结果支持 overlap 候选；新计数器已通过合成验证，但真实数据的全量 MAPQ fragment 复核尚未完成。仍需 cross-map 全量/代表性审查，以及等量、加入阶段和校准范围证据。当前两批数据都未由本版自动接受为正式归一化。
4. IgG 策略、去重、MAPQ、共识支持度和归一化在具体项目上的优劣，需要真实数据的位点检查、重复一致性与敏感性比较。软件验收不能代替该判断。
5. 显式输入适配器的布局契约已有测试；跨所有 nf-core 参数/版本的自动识别、原两批全量重新分析、de novo motif/GSEA、多个峰调用器与全部组合的性能比较，不在本版已验收声明内。
6. 交付包仅验证文件及链接资源，仍标记 REVIEW_REQUIRED。没有执行破坏性清理，也不会自动删除原始数据。

## v0.5.2 使用前核查与本机验收

- 完整 acceptance 9/9、58 项 Python 单元/契约测试通过；精简系统 Python＋bedtools 环境也通过全部 58 项测试，Python 编译及 Pixi 锁文件检查通过。[本版验收记录](docs/validation/2026-10-05/local_acceptance_v052.json) 保存状态与当前源码/示例/CI 配置哈希。
- 修正示例 2/3 支持度被四舍五入成 0.6667 后实际要求 3/3 的错误；补齐 target BAM 字典、peak 坐标、比较执行证据及样本映射校验；CI 显式安装 bedtools。详见 [使用前核查](docs/USABILITY_AUDIT.md)。
- 结论：可用于声明范围内的新项目，必要科学审核仍保留；没有将合成测试当作新实验的全量数据验收，也未独立确认本次提交的远端 Actions 状态。

## v0.5.1 本机验收（历史版本）

- Pixi 分析环境下完整 `acceptance.py` 9/9 通过；54 项 Python 单元/契约测试通过。新增验证覆盖 common-universe replicate correlation、blacklist subtraction 后的 50 bp 默认最小宽度和 crossmap 的 unit/biological sample 身份追溯。完整命令、状态和 58 个源码哈希见 [v0.5.1 验收记录](docs/validation/2026-10-04/local_acceptance_v051.json)。
- 合成项目 fixture 实际运行两个 strategy candidates，并在 pooled union 上统一重计数；项目初始化、片段/共识/count/QC/track DAG、续跑、交付搬迁验证均通过。依然只验证软件行为，不代表具体真实实验的候选方案已完成或被接受。

## v0.5.0 本机验收（历史版本）

- Pixi 分析环境下 `acceptance.py` 9/9 通过，其中包括 52 项 Python 单元/契约测试、spike-in BAM/竞争比对、真实微型下游方案矩阵 DAG、固定 pooled-union FRiP 比较、MACS2 12 种 shape/background/duplicate 组合、两种 DESeq2 归一化、注释和 motif fixture。完整摘要及脚本哈希见 [本版验收记录](docs/validation/2026-10-04/local_acceptance.json)。
- 新 matrix 支持候选计划登记、逐项运行、单候选失败保留和输入/模板未变化后的续跑。模型比较 fixture 核实了审核标志、相同设计/阈值/样本元数据哈希以及区域效应与方向变化汇总。
- 上述均为软件/合成数据验收。peak caller 候选需由用户在真实数据上执行并审核；DE 稳定性需要人工审核并正式运行每个模型。复杂/批次/配对设计仍未纳入正式模型矩阵。

运行步骤：[OPERATIONS](docs/OPERATIONS.md)；方法与审核：[SOP](docs/SOP.md)；历史缺陷审计：[AUDIT](docs/AUDIT.md)。
