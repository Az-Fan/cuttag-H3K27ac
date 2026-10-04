# 可选方案：测试证据与推荐

本页将软件实现通过与实验方案选择分开记录。测试时间：2026-10-03。每个新项目需要在固定输入上重做适用的敏感性比较；不能把本页某个样本的参数直接当作所有 CUT&Tag 文库的结论。

## 已实际执行的比较

| 决策 | 实际测试 | 当前推荐及边界 |
|---|---|---|
| 重复片段 | 合成 proper-pair BAM：保留/去除 duplicate 标记，与 MAPQ 20/30 交叉，共 4 组；准确得到 3、2、2、1 个片段 | H3K27ac 模板以保留重复作为明确基线，去重作敏感性分析；需同时考虑文库复杂度、深度、IgG 和位点分布。软件测试不能决定真实重复是否为 PCR 产物 |
| MAPQ | 同上，含 MAPQ 10、25、60 的配对片段；另有两端 MAPQ 不一致的测试 | 默认两端 MAPQ≥20；30 是更严格的敏感性候选。不能仅因总片段数较多选择阈值 |
| Spike-in 比对/计数 | 两批项目各取 Control/KD 的一个已剪切技术单元，各前 20,000 对，同一子集运行 4 种模式；模板另实现 spike-in-only BAM proper-primary-pair、双端 MAPQ20/30 fragment 计数和 target+spike-in 竞争比对 | 当前数据优先验证 end-to-end + overlap；保留 no-overlap 诊断，dovetail/local 不作无条件默认。子集计数不能进入正式校准；完整 paired-fragment 结果需按新流程生成并审阅 |
| IgG 策略 | 固定 MACS2 2.2.7.1 容器，BAMPE、keep-dup all：IgG 默认缩放 / scale-to-large / 无 IgG，共 3 组；均找回合成富集位点 | 优先评估匹配 IgG 的深度和结构；默认缩放是基线，另两组作敏感性证据。不能以峰数最大选策略。浅 IgG 放大可能放大噪声 |
| 共识峰支持度/黑名单 | 独立生物样本支持度、技术文件重复不增加支持、相接区间不算重叠、blacklist 逐碱基扣除都有行为测试 | 默认组内 2/3，向上取整；support-core 与 reproducible-union 两个 universe 并列评估。whole-interval drop 仅保留为敏感性模式；比较共识覆盖范围、重复一致性、下游方向稳定性 |
| 归一化 | 300 个合成峰 × 6 个样本，常规 DESeq2 和已知 spike-in 系数分别实跑；验证对比方向、size factor、normalized count 除法 | 两种模型回答的尺度不同，不按显著峰更多来选。确认等量、加入时点、计数方法和解释范围后才可接受 spike-in |
| ORA / motif 背景 | tested universe 的过滤/坐标/双向基因测试；clusterProfiler 实跑植入通路；HOMER 已知 motif 在合成 DNA 上实跑；不足窗口明确跳过 | 使用同一区域类别、同一统计宇宙背景；报告 ID 映射损失。远端峰最近基因只作关联。当前不声称 de novo motif、GSEA、不同峰调用器全面比较已验收 |

## 真实数据 spike-in 子集结果

每行均为同一份技术单元的前 20,000 对已经剪切的 reads。第二批的 unit1 是深度拆分文件之一，**不是独立生物重复**。跨批次不混合计算校准系数；输出显式 `calibration_group`。本次没有运行全量重计数。

| 技术单元 | end-to-end no-overlap | end-to-end overlap | end-to-end dovetail | local overlap |
|---|---:|---:|---:|---:|
| B1 Control unit1 | 0 | 10,782 | 10,782 | 10,795 |
| B1 KD unit1 | 0 | 2,806 | 2,806 | 2,810 |
| B2 Control unit1 | 0 | 10,764 | 10,764 | 10,781 |
| B2 KD unit1 | 0 | 2,760 | 2,760 | 2,762 |

这里的单位是 Bowtie2 concordant pairs，包含 concordantly exactly 1 和 >1，并不等于已接受的唯一比对高质量校准片段。原始表及子集 FASTQ 哈希见 [counts](validation/2026-10-03/spikein_subset_counts.tsv) 和 [provenance](validation/2026-10-03/spikein_subset_provenance.json)。

结果支持排查短片段与 `--no-overlap` 冲突：允许 overlap 恢复大量计数，dovetail 无额外收益，local 只略增。推荐先把 end-to-end + overlap 带入全量、MAPQ/多重比对、片段长度及实验加入阶段审核。不能由此宣称 lambda 校准已成立，也不能认定第一批 3 vs 3 的独立性已确认。

Bowtie2 对 overlapping mates 和 dovetail 的定义见[官方手册](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)。MACS 对小数据集放大缩放的解释见[官方 callpeak 文档](https://macs3-project.github.io/MACS/docs/callpeak.html)；本次实际执行的是上述固定 MACS2 容器。

## 可复用命令

```bash
# 输入 TSV: sample_id calibration_group fastq_1 fastq_2；请用已剪切、已正确合并的 paired reads
pixi run -e analysis python scripts/spikein_compare.py \
  --manifest examples/spikein_reads.tsv --fasta data/external/lambda.fa \
  --max-pairs 20000 --out results/spikein_matrix_01
# 选择后的全量计数；0 表示全量，仍不会自动接受校准
pixi run -e analysis python scripts/spikein_compare.py \
  --manifest examples/spikein_reads.tsv --fasta data/external/lambda.fa \
  --max-pairs 0 --modes end_to_end_overlap --out results/spikein_full_01

# 三种 IgG 策略，固定版本容器已通过 fetch_tools.py 安装后使用
pixi run -e analysis python scripts/peak_diagnostics.py \
  --manifest bams.tsv --gsize 2700000000 --out results/peak_matrix_01 --execute
# 输入 TSV: sample_id strategy peaks_bed
python3 scripts/compare_peaks.py --manifest peak_comparison.tsv --out results/peak_comparison_01
```

`spikein_compare.py` 保存同一输入子集、参数、日志和哈希；默认只取前 N 对，不宣称是随机代表性抽样。`read1_mapq20_pairs` 只是 read1 的描述性计数；正式计数策略还要审核两端质量。输入的不同技术文件不能冒充不同统计样本。每个 calibration_group 必须是可比较的实验加入/制备集合。

对比方法输出峰覆盖 bp 交集、Jaccard、双方被保留比例；QC 的 FRiP 应另外使用同一 master universe。任何方法都不会将 `accepted` 自动改成 true，也不自动将所有选项套在真实全数据上消耗资源。

## 官方流程与环境验收

nf-core/cutandrun 3.2.2 官方 `test,apptainer` 于本机实际跑完：133 个任务完成，2 个 Preseq 任务失败且被上游策略忽略，Nextflow 退出码 0。详见 [trace 摘要](validation/2026-10-03/official_test.json)。这验证了官方测试配置，不能替代模板 H3K27ac 配置在新实验上的全量验收。`audit_run.py` 将失败/重试/缺失 trace 单独报告，防止把退出码 0 描述为所有任务均成功。

原生 Bioconda MACS2 2.2.9.1 build 5 的 `__log_finite` 问题在本机真实复现，和[上游缺陷记录](https://github.com/bioconda/bioconda-recipes/issues/59362)一致。当前默认诊断使用已验证且固定哈希的 MACS2 2.2.7.1 容器；不会静默切换到未经测试的原生程序。


同一官方微型输入还使用 `--test-template-params` 跑通模板默认 MACS2 narrow / target 保留重复 / CPM 组合。退出码、逐任务摘要、参数和实际 narrowPeak 产物计数见 [模板参数测试](validation/2026-10-03/template_params_test.json)。这检验参数连接，未将官方测试数据解释为真实 H3K27ac 生物实验。
