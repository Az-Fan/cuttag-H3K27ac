# CUT&Tag 项目模板 v0.4.1

从 ACLY-CUTTAG 两批项目沉淀的独立模板，面向 paired-end H3K27ac：提供上游运行、科学审核停点、可续跑下游、可选方案对照和可搬迁交付。支持范围及操作步骤见 [操作指南](docs/OPERATIONS.md)，方法细节见 [SOP](docs/SOP.md)。

当前已通过本机合成集成验收和官方 nf-core 测试运行，并对两批真实数据的 spike-in 子集做了参数比较。**这不代表所有参数在新实验上都已得到生物学验证。** 实际结果、条件性推荐和限制见 [可选方案验证](docs/OPTION_VALIDATION.md) 与 [验收状态](PROJECT_STATUS.md)。

## 创建新项目

```bash
python3 scripts/init_project.py --out /path/to/new-study --project-id new_study
cd /path/to/new-study
pixi install --locked -e analysis
python3 scripts/fetch_tools.py --tool macs2
python3 scripts/install_bioc_data.py --out logs/bioc_data_setup
pixi run -e analysis python scripts/acceptance.py --out results/acceptance_01
```

初始化不复制数据、缓存、结果或 Git 历史，也不继承科学审核结论。编辑 `config/project.json`、`config/samples.tsv`，登记真实参考和原始 FASTQ；示例是每组一个生物样本，不满足正式差异推断要求。

## 工作流程

```text
样本/实验登记 → 环境和官方 test → FASTQ/参考完整性 → alignment
  → QC、对照及 spike-in 参数比较 → 人工科学审核 → production
  → 结果适配/显式登记 → fragments、共识峰、counts、QC、轨迹
  → 峰和统计模型审核 → 差异/注释 → 同宇宙 ORA、motif
  → 图表/方法/限制 → 复制交付包 → 搬迁后的哈希校验
```

```bash
pixi run -e analysis python scripts/cuttag.py validate --allow-missing
pixi run -e analysis python scripts/cuttag.py test --run-id official_test_01
# 将成功的 test run ID 登记到 qc.official_test_run_id，并审阅 task_audit.json
pixi run -e analysis python scripts/cuttag.py inventory
pixi run -e analysis python scripts/cuttag.py run --stage alignment
# 按 SOP 完成生产审核；通过后运行 production
pixi run -e analysis python scripts/cuttag.py run --stage production
```

生产审核绑定配置、样本表、FASTQ 清单和参考哈希；换输入后旧审核失效。`COMPUTATIONAL_PASS` 只表示计算入口完成；Nextflow 中被忽略的失败由 `task_audit.json` 单独暴露。

## 下游自动计划

用严格限定版本/布局的 `import_nfcore.py` 适配模板生产结果，或按 `examples/artifacts.tsv` 显式登记每个生物样本的已合并 BAM 和选定峰文件。

```bash
pixi run -e analysis python scripts/plan_downstream.py --config config/project.json \
  --artifacts examples/artifacts.tsv --out results/downstream_01
pixi run -e analysis python scripts/workflow.py --workflow results/downstream_01/workflow.json \
  --out results/execution_01
```

生成原始片段计数、同一峰宇宙的 FRiP/片段长度 QC、CPM BigWig 与 IGV。`--resume` 只复用哈希一致的成功步骤；变化或已有失败产物要求新目录。正式差异分析仍经 `prepare_review.py` / `run_differential.py` 审核入口，不会自动接受样本、峰策略或校准。

## 必须保留的判断

- 技术拆分不是生物重复；1 vs 1 只允许描述性分析。只支持 `~ condition`，不隐式忽略批次/配对因素。
- 重复片段、MAPQ、IgG 缩放和峰集规则都有明确基线和敏感性测试；不按峰数最多选方案。
- Spike-in 正式计数使用 spike-in-only proper paired fragments、两端 MAPQ 门槛及技术 unit 到 biological sample 汇总；竞争比对、等量、加入阶段和解释范围仍需审核，子集测试不替代全量计数。
- ORA/motif 使用明确的 tested universe，记录基因映射损失；邻近基因不等于已验证调控靶点。
- 缺失 QC 明确为 NA；缺必需产物会失败，计算通过不自动改为科学接受。

可复用命令、当前推荐及验证证据均在 [OPTION_VALIDATION.md](docs/OPTION_VALIDATION.md)。交付请使用 `finalize.py` 和 `verify_release.py`，详见 [操作指南](docs/OPERATIONS.md)。
