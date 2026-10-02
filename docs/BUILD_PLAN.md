# CUT&Tag 项目模板（规划稿）

本模板从 `ACLY-CUTTAG` 和 `ACLY-CUTTAG-batch2` 两个项目沉淀而来，目标是建立可复制、可审计、不会把技术重复误当生物学重复的 CUT&Tag/CUT&RUN 项目骨架。本文先定义模板范围与实施顺序；模板应经过小型合成数据和真实项目迁移验证后再标为稳定版。

## 设计原则

1. **实验设计是配置数据，不从文件名猜测。** 每个生物样本、文库、测序 lane、FASTQ 文件组分别登记；技术拆分可以合并或作为同一生物样本输入，但不能增加生物学重复数。
2. **计算完成与科学验收分开。** 每个阶段分别记录软件状态、数据 QC 状态、实验设计状态和分析解释限制；上游跑完不自动允许下游推断。
3. **所有输入与结果可追溯。** 记录 SHA256、来源、参考版本、参数、命令、软件/容器版本、时间、退出码和日志；不得依赖隐式绝对路径或漂移的 latest 标签。
4. **不覆盖已发布结果。** 每次运行写入带时间戳或 run ID 的目录；正式结果经过验收后再通过明确的 current 指针或报告索引发布。
5. **大文件与源码分开。** 项目仓库放配置、脚本、文档、小型测试数据；FASTQ、BAM、参考索引、容器和缓存存外置数据/缓存目录，并以清单与校验和登记。
6. **探索性敏感性分析明示其身份。** IgG/no-IgG、CPM/median-ratio/spike-in、peak caller 等比较要分别产出并说明适用边界，不自动选取“最好看”的方案。

## 建议目录

```text
cuttag-project/
├── README.md
├── PROJECT_STATUS.md
├── analysis_log.md
├── .gitignore
├── pixi.toml / pixi.lock
├── config/
│   ├── project.yaml                 # 物种、基因组、靶标、比较、方向定义
│   ├── experiment_design.yaml       # 生物样本/文库/技术拆分/分组/重复关系
│   ├── samplesheet.csv              # nf-core 输入样本表，由校验脚本生成或核验
│   ├── params/                      # test、alignment、production、diagnostic 参数
│   ├── references.yaml              # FASTA/GTF/blacklist/spike-in 来源与哈希
│   └── profiles/                    # 本地执行资源配置，不含机器专属路径
├── data/
│   ├── raw/                         # 可链接至外置只读原始数据
│   ├── external/                    # 外部资源清单；大文件本身可外置
│   └── manifest/                    # 文件、校验和、来源、交付记录
├── results/
│   ├── runs/<run_id>/               # 不可变运行结果与 pipeline_info
│   ├── qc/                          # 汇总 QC 与验收表
│   ├── peaks/                       # caller/策略分开的峰集
│   ├── signal/                      # CPM、spike-in 等分开
│   ├── differential/<model_id>/
│   ├── annotation/
│   ├── downstream/
│   └── release/                     # 经审核的交付索引/报告
├── scripts/
│   ├── 00_validate_project.py
│   ├── 01_inventory_fastq.py
│   ├── 02_run_nfcore.sh
│   ├── 03_qc_report.py
│   ├── 04_spikein_audit.py
│   ├── 05_build_peak_matrix.R
│   ├── 06_differential.R
│   ├── 07_annotate_peaks.R
│   ├── 08_downstream.R
│   └── 90_finalize_release.py
├── docs/
│   ├── SOP.md
│   ├── QC_ACCEPTANCE.md
│   ├── METHODS.md
│   ├── INTERPRETATION_LIMITS.md
│   ├── TROUBLESHOOTING.md
│   └── decisions/
└── tests/
    ├── fixtures/                    # 极小合成文件，不放真实样本
    └── expected_contracts/
```

数据根目录、Nextflow work/cache、Apptainer 缓存和临时目录由 `config/paths.env` 或环境变量配置；路径文件不提交，提供 `.example`。结果清理脚本只能清除由模板登记的缓存/work 路径，必须先显示待删路径、估计空间与保留产物，并输出删除清单及校验记录。

## 配置与样本契约

`experiment_design.yaml` 至少区分：`biological_sample_id`、`library_id`、`sequencing_unit_id`、`fastq_r1/r2`、`condition`、`assay_target`、`control_type`、`biological_replicate_id`、`technical_split_of`、`batch`、`lane`、`input_amount`、`spikein_amount`、`spikein_added_at`。未知字段明确写 `unknown`，并带来源/确认日期；不能把未知静默当作默认值。

校验器应检查：文件存在/可读、R1/R2 配对和 read 数、样本 ID 唯一性、角色合法性、对照归属、分组与重复映射、参考/黑名单/Spike-in 版本和哈希、设计矩阵可估计性、每组独立生物样本数。样本表中的技术拆分标识不得生成额外的生物重复。单生物样本/组项目允许上游描述性分析，但差异推断应标注“不具备生物学重复推断条件”，并关闭正式 DEG/差异峰推断入口；用户显式要求探索性分析时，报告必须保留此限制。

Spike-in 配置同时记录物种/序列 accession、FASTA 与索引来源、加入量、加入时点、对照样本角色、比对几何/过滤规则。不得只依赖 nf-core 默认计数；如计数与 FastQC 序列证据、单端筛查或预期明显不符，先执行预注册的参数审计，再审批生产比例因子。

## 全流程和阶段门禁

| 阶段 | 主要工作 | 关键产物 | 进入下一阶段的条件 |
|---|---|---|---|
| 0 项目建档 | 记录实验设计、来源、靶标、对照、参考和问题 | 项目配置、设计确认单、待确认项 | 样本身份和分组有来源；未知显式标出 |
| 1 输入检查 | FASTQ inventory、哈希、配对、读数、污染/接头快速检查 | 输入 manifest、初始 QC | 所有输入可追溯，损坏/缺失有处理记录 |
| 2 测试运行 | 锁定 nf-core/cutandrun 版本、容器与本地参考，运行官方 test 与小规模真实子集 | versions、test report、命令日志 | 环境与输出契约通过；空间/资源预算可行 |
| 3 上游生产 | trimming、比对、过滤、重复标记、片段与 MultiQC | run 目录、BAM、pipeline_info、QC 表 | 任务成功且独立 QC 审阅；不因成功码自动验收 |
| 4 QC 审查 | 比对率、复杂度/重复、片段长度、FRiP、TSS/背景、重复相关性、样本身份 | QC atlas、样本决策表 | 每样本 PASS/REVIEW/FAIL 有证据和签核；排除有理由 |
| 5 Spike-in/信号 | 参考确认、计数参数审计、全量计数、比例因子与 BigWig | spike-in 审计、CPM 与校准轨迹 | 计数方法能解释序列证据；实验加入时点足以支持目标解释 |
| 6 峰与矩阵 | 比较合理 peak caller/control 方案、重复支持、黑名单处理、计数矩阵 | 方案分开的 master peaks 与 raw counts | 峰集有浏览器/FRiP/重复支持审查，保留全部模型定义 |
| 7 差异分析 | 根据设计采用可支持的统计单位与模型；归一化敏感性分析 | 完整统计表、设计矩阵和模型诊断 | 有足够独立样本及可估计设计；否则仅描述，不报告推断结论 |
| 8 注释与下游 | 峰注释、区域分类、ORA/motif、轨迹/热图 | 全量表、图、数据库/版本背景 | 说明最近基因不等于靶基因，motif 不等于 TF 结合，富集为关联 |
| 9 交付与归档 | 汇总方法、限制、结果导航、运行清单和存储清理 | release report、provenance、manifest | 独立复核文件完整、路径可迁移、声明与证据一致 |

门禁状态采用 `NOT_STARTED / RUNNING / COMPUTATIONAL_PASS / QC_REVIEW / ACCEPTED / BLOCKED`；每个状态记录负责人、时间、证据文件和理由。`COMPUTATIONAL_PASS` 不等于 `ACCEPTED`。

## 分析选择需要保持为显式模型

- **Peak calling：** 有可靠匹配 IgG 时按预设方案使用；极浅或异常 IgG 时，单独报告其可用性，比较背景策略/无对照敏感性结果并保留诊断标签。不能将无对照诊断峰直接改称正式峰集。
- **Spike-in 与 CPM：** 同时生成时输出独立目录和单位，禁止跨归一化直接比较轨迹高度。记录 λ/E. coli 等 spike-in 以及加入实验步骤；外源 DNA 只能校准其加入之后经历的步骤。
- **重复与峰支持：** 配置按组定义支持比例，分母为独立生物样本数；技术拆分须先按样本合并或合适处理。报告中用“文件/文库支持”还是“生物重复支持”必须与设计一致。
- **差异分析：** 输入原始片段 counts，不用 BigWig。保存完整测试峰宇宙、设计公式、对比方向、过滤门槛、size factors、软件版本及模型诊断。比较峰集不同时基于基因组区间重叠，不直接按 peak ID 比较。
- **功能解释：** 位置注释、ORA 背景、motif 区域与背景集、数据库版本和失败状态全部留档；无显著结果也作为有效结果交付。

## 环境与复现

建议采用 Pixi 管理 Nextflow、Java、Apptainer 和小型审计工具；锁定 `pixi.lock`，nf-core 流程使用明确 release/tag，容器使用固定 digest 或经验证的版本。以 `profiles/test.config`、`profiles/local.config`、`profiles/production.config` 分层；生产配置只放通用资源，不写开发机绝对路径。脚本从项目根目录解析相对路径，所有命令写入 run manifest。离线参考资源登记来源、版本、下载日期、SHA256、索引构建命令和索引工具版本。

至少设置 `validate`、`versions`、`nf-test`/官方 test、`dry-run` 和 `run` 入口；先验证再生产。环境安装、参考索引和容器缓存写到独立缓存区，不提交 `.pixi/`、Nextflow work、容器层和中间 BAM。

## 日志、存储和交付

每次运行保存配置快照、输入 manifest、git commit、软件版本、命令行、Nextflow trace/report/timeline/DAG、日志、输出文件校验和及状态。任何手动覆盖、修复和过滤另写决策记录。清理前列出文件和大小；原始输入、已发布结果、当前报告和关键 provenance 受保护；删除后对受保护产物做大小/哈希复核。软链接目标必须在 manifest 中解析，交付前检查悬空链接及绝对路径。

交付首页概括样本设计、实验靶标、比较方向、完成阶段、主要 QC、可用结果、未完成项、归一化/峰方案与解释限制。明确“文件数不等于独立生物重复数”；明确峰信号不等于 RNA 表达，最近基因不等于调控靶点，spike-in 不能单独证明全局修饰变化。

## 实施路线

1. **整理规范和契约：** 把本规划、配置 schema、目录、阶段状态表与样本映射规则定稿。
2. **建立最小骨架：** README、模板配置、样本表、Pixi/Nextflow 配置、路径示例、日志/provenance 契约和 `.gitignore`。
3. **优先完成验证入口：** 项目/FASTQ/样本表/参考校验，输出机器可读报告；验证重复层级和统计可估计性。
4. **完成上游可复跑：** test、nf-core 生产入口、资源配置、QC 汇总和断点续跑说明。
5. **完成 CUT&Tag 专项审计：** spike-in 参数矩阵/全量计数、IgG 深度及 caller 敏感性、重复支持峰、片段 counts、CPM/Spike-in 轨迹。
6. **加入统计与下游可选模块：** 差异、注释、ORA、motif、图表；各自有独立输入契约和报告限制。
7. **以两批项目回放验证：** 仅把配置和必要脚本映射到模板，分别验证多重复但来源需审查、单生物样本多文件拆分两种设计，确保模板正确拦截不适当推断。
8. **发布版本：** 补齐 SOP、故障排查、示例输出与迁移说明，记录模板版本/变更日志；之后再将此目录标记为稳定模板。

## 当前已知来源经验

- 第一批完成上游和差异分析，但重复来源曾需核查；λ 计数受到 Bowtie2 overlap 参数影响，IgG 极浅也使峰调用需诊断。修正计数、背景与重复支持规则必须成为有记录的决策，而不能固化成隐性默认值。
- 第二批最终确认每组只有 1 个生物样本、三份文件属于测序拆分；因此模板要允许上游 QC、信号轨迹和描述性峰诊断，同时阻止将三份文件写成 3 个生物重复并据此进行正式差异推断。
- 两项目配置曾含机器绝对路径、缓存和软链接；模板改用项目相对路径/环境变量，并在交付与清理时解析链接和校验受保护文件。
- 两批有历史配置、诊断和生产参数并存的情况；模板需区分 `test / diagnostic / production` 参数文件，并将每次实际运行配置快照进对应 run 目录。
