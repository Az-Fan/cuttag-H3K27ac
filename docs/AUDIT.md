# 完整性、正确性与必要判断核查

核查起点：Git `46e2d24`（v0.3.0）。范围：所有模板脚本、项目配置、文档与既有测试；对照两个来源项目的设计记录，以及本地 nf-core/cutandrun 3.2.2 源码的 spike-in 分支。

## 结论

此前“模板主体补齐”的表述过于乐观：已覆盖主要阶段，但尚不是完成科学验收和端到端回放的稳定全流程。部分脚本是独立入口，部分要求只有文字说明，且发现了会影响数据处理和统计前提的实际缺陷。本轮已修复下表的问题，并保留未完成清单。不能把已有 11 个测试通过等同于全部工具或方法正确。

## 已修复问题

| 编号 | 原问题与影响 | 本轮处理 |
|---|---|---|
| F01 | filtering-only 固定 CPM；本地 nf-core 的 prepare_genome 仅在 Spikein 模式建立 spike-in 分支，导致启用 spike-in 的前期 QC 缺少计数 | 启用时在 filtering-only 阶段也设置 Spikein；此时仍不进入下游峰/归一化阶段 |
| F02 | overlap 诊断配置存在但控制器不能传入；审核修正后生产仍可能用默认 no-overlap | 配置增加 alignment_profile；计划和执行均传入明确配置；生产审核绑定该文件哈希；显式 FASTA 禁止继承默认 spike-in 索引 |
| F03 | upstream/calibration 仅有布尔开关，无 IgG 背景策略依据；未知加入阶段、未确认等量仍可放行 | 生产需 QC/control-strategy 审核记录、证据文件及 config/sample 哈希；spike-in 要求当前实现支持的等量前提、加入阶段、解释范围与计数/校准接受 |
| F04 | peaks_accepted 无人读取；直接调用 R 不检查峰集/QC 决策 | 新增配置化差异入口；直接 R 入口同样要求模型审核记录；绑定 counts、metadata、spike-in 表及模型/阈值 |
| F05 | config 中阈值与设计不传给 R，可能悄悄运行不同模型 | run_differential.py 显式传递阈值、方向、归一化；不支持的 batch/配对公式明确拒绝，不自动简化成 ~condition |
| F06 | 字符串 "false" 在 Python 中为真；缺失比较组/损坏样本表检查不充分 | 布尔字段严格类型检查；比较组及方向、样本列值、阈值与归一化检查；禁止同一 pipeline group 混合条件、禁止前导零 replicate 混淆整数编号 |
| F07 | R 仅按前四列拿坐标；spike-in 重复 sample ID 可被 match 静默取首行 | 按字段名拿坐标；拒绝非法/缺失坐标、重复列名、额外 count 样本、重复 spike-in ID |
| F08 | 注释后按输入行顺序重赋 peak ID，依赖第三方保持排序 | GRanges metadata 携带 peak_id；按 ID 重排并断言 chromosome/start/end 一致；真实 TxDb 实测通过。原实现存在排序假设风险，并非已证明每次都发生错配 |
| F09 | tracks 不核查已接受的 spike-in；NaN 可越过 <=0；重复 sample ID 覆盖输出；IGV genome=unknown；大 BED 全量读入哈希 | spike-in 轨迹强制审核记录；拒绝非有限/非正系数和重复 ID；显式 --genome；流式哈希 |
| F10 | motif 窗口未检查染色体右边界；background 与 target 冲突删除未留痕 | 强制 FASTA/.fai；越界、目标背景重叠写排除清单；最少窗口数必须为正 |
| F11 | workflow 只看退出码；命令没生成预期结果也可成功 | steps.outputs 声明后检查文件存在，缺失则失败并拦截依赖。未声明 outputs 的通用命令仍仅有退出码契约 |
| F12 | 片段路径解析依赖 cwd；QC 的 nan 被当可绘制数值；片段 filter/sort 命令未留档 | fragment manifest 路径相对 manifest；严格 BED4；QC 非有限值列为缺失、比例越界拒绝；保存 filter/sort 命令与不成对/TLEN 异常丢弃计数 |
| F13 | run 不要求 FASTQ 完整检查，只检查文件存在 | 真实 run 在 Nextflow 前执行完整配对/结构/gzip/哈希检查；生产时还比对已审核输入清单哈希，失败写入 input_integrity 状态。大文件因此有额外全量读取开销 |

审核记录是科研决策的可追溯载体，不能自动证明人为填写的样本独立性或证据真实性。模板不会替使用者发明接受理由。为非等量 spike-in、复杂 batch/配对设计建模仍需单独实现，不能把布尔字段改为 true 作为替代。

## 必要判断的实际位置

| 判断 | 当前状态 |
|---|---|
| 技术拆分不能变成独立重复 | validate 检查同一 biological_sample_id 的组/重复映射；共识峰按声明的生物样本去重；R 拒绝同一生物 ID 重复出现。独立性事实仍需建库证据 |
| 1 vs 1 不做正式条件推断 | 校验和 R 均阻止；描述性上游仍可运行。每组 >=2 是程序最低条件，不代表实验设计质量已充分 |
| IgG 是否过浅/背景策略是否可用 | 诊断矩阵保留，无自动“选峰最多方案”；生产要求 control_strategy_accepted 及证据。尚无从 BAM 自动判定浅对照的统一数值阈值 |
| target duplicates 与 ATAC shift | H3K27ac 默认保留 duplicate-marked fragments；无 ATAC shift；去重必须显式选择。上游参数仍须按 assay/靶标审核 |
| spike-in 身份、加入阶段、计数几何、归一化解释 | 生产/差异/轨迹的接受条件加强；计数审计输出默认未接受；能自动读取日志，不能自动证明实验前提或完成所有比对模式 |
| 峰集筛选与 raw counts | 基于声明的生物样本支持度；peak universe 的接受由审核证据记录；差异输入为整数 raw counts。整数检查无法鉴别被人为四舍五入的归一化矩阵来源 |
| 富集背景与远端基因关联 | ORA 检查输入是显式 universe 子集；最近基因只表示位置关联的解释仍保留。tested universe、promoter/distal 分层与 ID 映射仍需准备，未自动接全 |
| motif 目标/背景来源 | 检查等长、重叠、边界和最低数量；尚未自动证明均来自同一实际检验峰集或均为远端候选区域 |
| 计算成功 ≠ 科学接受 | 状态仍为 COMPUTATIONAL_PASS；release 仍 REVIEW_REQUIRED；人工决策不会由退出码自动改变 |

## 实际验证范围

- 26 项 Python 测试通过，其中新增 15 项覆盖错误布尔类型、无效设计/对比、缺失实验前提、过期审核输入、未接受峰集、spike-in filtering 计划与结果文件缺失。
- 真实 samtools/bedtools/pyBigWig：微型合成 BAM → paired fragments → peak counts → CPM BigWig → QC PDF/HTML，通过。检查双端 MAPQ、重复保留/去重、CPM 分母、manifest 相对路径、NaN 处理。
- 使用系统 R 和已有 DESeq2：300 峰、6 个合成样本，经配置化差异入口+审核记录拟合并绘图，通过，方向正确。
- 使用本地 hg38 TxDb SQLite 与系统 ChIPseeker：3 个跨染色体乱序峰，peak_id 与 BED 坐标完整保留。
- `pixi lock --check --offline --no-config` 通过：当前 manifest/lock 一致。局部可执行文件可用，不据此声称所有安装脚本、R 数据包或容器已完整验收。
- 全部 Python/R 脚本语法检查通过。真实 nf-core 官方 test、HOMER 完整流程、实际文库的 end-to-end 验收未执行。

## 仍不完整的部分

1. nf-core 输出到生物样本 fragment/peak/统计 metadata 的自动适配，统一带输入输出契约的项目 workflow 示例与断点续跑。
2. spike-in 参数矩阵自动比对与全数据重计数；当前只有配置与日志审计，生产审核需引用外部完成的证据。
3. 从 BAM/peaks 自动汇总 FRiP、有效背景、复杂度、相关性等 QC 指标；现有 atlas 读取准备好的表。
4. FASTA/GTF/blacklist/TxDb 版本、染色体兼容性和资源哈希统一校验；目前部分只检查存在，模板没有下载并鉴定所有参考。
5. ORA 基因宇宙/映射与 promoter/distal 分层、motif 宇宙验证、不同峰方案×不同归一化的完整敏感性汇总。
6. 完整移交包（文件复制、路径重写、IGV 迁移验证）、受保护清理和完整 provenance；目前 release 只是带 SHA256 的审阅索引。
7. 清洁锁定环境下所有外部模块和官方测试、真实项目回放的验收。

这些缺口包含实现缺口，不能全部归因于网络。当前结论：有可用且经针对性验证的模块，必要判断已经补强，但尚未达到原规划的完整稳定项目模板标准。
