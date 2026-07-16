# DIT-HAP Pipeline 性能优化审计与设计

日期: 2026-07-16
项目基准: `Spore2YES6_1328`(1 sample × 8 timepoint × 1 condition,无生物学重复分支)
状态: 设计 / 待评审

---

## 1. 背景与方法

本文档对 DIT-HAP snakemake 上游 pipeline 逐步进行性能审计,并给出分层的优化方案。
所有瓶颈判断基于 **实测**,而非猜测:

- 从 `projects/Spore2YES6_1328/logs/**` 的日志时间戳计算每个 rule 的墙钟耗时;
- 从 `results/**` 的中间文件体积评估 I/O 压力;
- 通读各步骤脚本,区分"已向量化/多核"与"逐行 Python 循环/单核"。

优化目标:在**不改变任何数值结果**的前提下缩短总墙钟时间;结构性重构(会改数值实现或 CLI 契约的)单独标注并放到最后一层。

---

## 2. 数据规模与并行度约束(关键前提)

基准项目是 **单样本 × 8 时间点 × 单条件**。这带来一个决定性约束:

> read-processing 阶段按 `{sample}_{timepoint}_{condition}` 展开,可并行的 job 数 **只有 ~8 个**(每 timepoint 一个),不是几十上百。

推论:
1. **"放开调度并行度"收益有限** —— 因为 job 本来就少。真正的杠杆在**单个 rule 内部**(PBL/PBR 拆分、脚本向量化、脚本内多核)。
2. **例外:内存声明 bug 会把这 8 路也压成 1 路** —— 见 Tier 1.2,这是当前最廉价的墙钟收益。
3. **depletion 阶段是全局单文件**(112,987 条 insertion 汇总成一张表),没有任何 job 级并行度,curve_fitting 的 29 分钟**只能靠脚本内并行解决**。

---

## 3. 实测耗时基线(单样本,当前串行执行)

| 步骤 | 脚本/工具 | 单文件实测 | 8 tp 累计(当前串行) | 特征 |
|---|---|---|---|---|
| 6 bam_to_tsv | parse_bam_to_tsv.py | 300–751s | ~53 min | 单核解析循环;200GB 内存声明压死并行 |
| 7 filter_aligned_reads | filter_aligned_reads.py | 206–529s | ~37 min | 单核 pandas;PBL/PBR 同 rule 串行 |
| 8 extract_insertion_sites | extract_insertion_sites.py | 129–328s | ~20 min | `.apply(axis=1)` 逐行;PBL/PBR 串行 |
| 10 annotate | annotate_genomic_features.py | ~330s | — | 3× `.apply(axis=1)` + groupby Python 循环 |
| **15 insertion curve_fitting** | curve_fitting.py | **1729s** | **~29 min** | **单核 scipy.minimize 逐行,112,987 条** |
| 17 gene curve_fitting | curve_fitting.py | 60s | 1 min | 同上,仅 4,498 条 |

其余(fastp / cutadapt / bwa / samtools)已多线程,不是主要瓶颈,但存在**关键路径**与 **PBL/PBR 串行**的调度问题(见 Tier 2)。

### 中间产物磁盘占用(意外发现)
| 目录 | 体积 | 是否 `temp()` |
|---|---|---|
| `6_filtered/` | **18 GB** | ❌ 否 |
| `4_sorted/` | 8.3 GB | ❌ 否 |
| `5_tabulated/` | (已清) | ✅ 是 |

`6_filtered`、`4_sorted` 是纯中间产物却未声明 `temp()`,长期占用几十 GB 磁盘。

---

## 4. 优化方案(按收益/风险分层)

### Tier 1 — 高收益、低风险、不改数值(优先做)

**1.1 curve_fitting 脚本内并行化 —— 最大单点收益 ✅ 已实现(2026-07-16)**
- 现状: [curve_fitting.py:529-542](../../workflow/scripts/depletion_scoring/curve_fitting.py#L529-L542) 是逐行 `scipy.minimize`(maxiter=3000)的 Python for 循环,单核,112,987 条 → 1729s。
- 每条 insertion 完全独立(embarrassingly parallel)。
- 实现: 新增模块级 picklable worker `fit_and_augment`,串行循环换成 `joblib.Parallel(n_jobs=jobs)`(默认保序 → 输出与串行逐字节等价)。新增 `-j/--jobs` 参数,两个 curve_fitting rule 传 `-j {threads}`(`threads: 16`)。numpy import 前 `os.environ.setdefault("OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS","1")` 防 worker 内 BLAS 过订阅。joblib 加进 env。移除了 tqdm 进度条(与并行不兼容)。
- **实测: 1729s → 136s(-j 16),12.7× 加速**,RSS ~675MB。

**端到端验证(2026-07-16,真实 snakemake --use-conda 路径):**
- 强制重跑 step-15 rule(mtime 触发,隔离单 job),stats conda env 因新增 joblib 重建成功 → 验证了 joblib 在正式环境可解析、`-j {threads}` 正确传参(日志确认 "Fitting sigmoid curves with 16 worker(s)")。真实拟合耗时 **138s**,与手动实测一致。
- **并行等价性(核心命题,已确认):** 固定线程环境下,全量 112,987 行 serial(-j 1)vs parallel(-j 16)输出**逐字节 `diff` 一致**(body md5 `92ef79…`),-j16 两次运行也一致。**joblib 并行本身精确、可复现、无损。**

- **根因排查:一个我最初误判、最终查明的问题(诚实记录):**
  - 现象: 通过 snakemake 跑的 -j16(body md5 `b1b322…`)与直接调脚本的 -j16(`92ef79…`)在 **1,588 行**上有真实 R² 差异(非舍入),其中数百行是 R²=0.83↔-2.3 的灾难性翻转。
  - **我最初错误地归因为**"跨 conda 环境的 BLAS 构建不同"——错的,两者用的是同一个 env(`3fae…`)。我还一度误用了**编辑脚本前产生的陈旧输出**做基线,导致虚假的"13,188 行差异"警报。
  - **实际根因(已证):** snakemake 的 `shell.py` 会把 `OMP_NUM_THREADS/OPENBLAS_NUM_THREADS` 等**设为该 rule 的 `threads` 值**(这里 16)并注入 shell。脚本原来的 `os.environ.setdefault(...)` **无法覆盖已设值** → 在 snakemake 下 BLAS 实际以 16 线程运行。多线程 BLAS 让 `scipy.minimize` 在病态拟合上收敛到不同局部最优。
  - **验证:** 直接以 `OMP_NUM_THREADS=16` 跑脚本 → 精确复现 snakemake 的 `b1b322…`;`OMP=1` → `92ef79…`。
  - **非本次并行化引入**(原串行版也在 snakemake 的 OMP=threads 下运行),但这是真实的**不可复现性 bug**。

- **修复(已应用):** 把线程 env 变量从 `setdefault` 改为**无条件赋值 `=1`**(6 个变量含 VECLIB/GOTO),保证无论如何启动都单线程 BLAS、输出线程数无关且可复现。
- **e2e 最终验证(已通过):** 修复后经真实 snakemake 路径,`--cores 16` 与 `--cores 1` 输出 body md5 **均为 `92ef79…`**,彼此一致且等于直接调用的 OMP=1 参考。**输出线程数无关、可复现、serial==parallel。** 买回的 12.7× 加速在正式环境成立且数值正确。

**1.2 bam_to_tsv 内存声明 → 4000 —— ⚠️ 稳健性修正,非加速(更正:2026-07-16)**
- **原判断有误。** 我最初声称 `mem_mb=200000`(200GB)导致 bam_to_tsv 被迫串行、修复后 53min→10min。**经实测推翻:** snakemake 的 `mem_mb` 只有在命令行显式传 `--resources mem_mb=<总量>` 时才约束调度;本项目实际用 `snakemake --use-conda --cores 16`(无 profile、无 `--resources`),该声明**从未**影响调度。用玩具 workflow 验证: `--cores 4` 下 4 个 200GB job 照样并行(span 3.0s),仅 `--cores 4 --resources mem_mb=200000` 时才串行(span 12s)。
- **真正让 bam_to_tsv 看似串行的是 `threads: 8`**(见 §2.1):`--cores 16` 下并发度 = `floor(16/8)=2`。上一次真实运行的日志时间戳证实 8 个 timepoint 基本是 2 路并发跑完。
- 保留 `mem_mb=4000` 改动: 无害且更稳健(将来若有人用 `--resources` 调度不会被 200GB 卡死),但**不计入任何加速**。

**1.3 ~~给 `6_filtered` / `4_sorted` 加 `temp()`~~ —— 已决定不做(DP2)**
- 保持现状,`6_filtered`(18GB)/`4_sorted`(8.3GB)保留供 debug / QC 复查。此项不纳入实施范围。

### Tier 2 — read-processing 并发重设计(批 2 焦点)

#### 调度模型(实测确立的正确前提)
- 本项目用 `snakemake --use-conda --cores 16` 调用,**无 profile、无 `--resources`**。
- 在此模式下,并发度**完全由 `rule.threads` 决定**: 同时运行的某 rule 实例数 = `floor(cores / rule.threads)`。`mem_mb` 不参与(见 §1.2 更正)。
- 主机有 **64 核**,但 pipeline 只用 `--cores 16` —— 提高 `--cores` 本身就是一条免费杠杆(独立于代码改动)。

#### 当前各 read-processing rule 的 threads 与真实多核利用
| rule | threads | 并发@16 | 工具/脚本真的用满 threads 吗? |
|---|---|---|---|
| fastp | 6 | 2 | ✅ fastp 多线程 |
| junction_classification (cutadapt) | 6 | 2 | ✅ cutadapt `--cores` |
| fastqc | 4 | 4 | ⚠️ threads=并发文件数,但 4 文件被拆成 **4 条串行命令**(见 2.3) |
| bwa_mem_mapping | 8 | 2 | ✅ bwa-mem2;但 PBL/PBR 同 rule 串行 |
| samtools sort/index | 2 | 8 | 部分;PBL/PBR 串行 |
| **bam_to_tsv** | **8** | **2** | ❌ **解析循环单核 Python**,threads 只喂 pysam 解压 → 声明 8 但基本用不到,白白把并发压到 2 |
| **filter_aligned_reads** | **(无)=1** | **16** | ❌ 单核 pandas |
| **extract_insertion_sites** | **(无)=1** | **16** | ❌ 单核 pandas + `.apply` |

关键错配: **bam_to_tsv 用 threads 换不到吞吐,却牺牲了并发**;而 filter/extract 是真单核却已能 16 路并发(受限于上游产物就绪节奏 + IO)。

#### 2.1 bam_to_tsv `threads: 8 → 1` ✅ 已实现(2026-07-16,基于实测)
- 实测 `parse_bam_to_tsv.py` 在同一 BAM(8.37M alignments)下:

  | threads | CPU% | 墙钟 | 输出 md5 |
  |---|---|---|---|
  | 1 | 99% | 177s | 6a4a97… |
  | 4 | 102% | 179s | 6a4a97…(相同) |
  | 8 | 102% | 174s | 6a4a97…(相同) |

- 结论: pysam 解压线程**零收益**(CPU% 始终 ~100% = 单核),瓶颈是单核 Python 解析循环。`threads: 8` 换不到吞吐却把并发压到 `floor(16/8)=2`。
- 实现: 降到 `threads: 1` → 并发 `floor(16/1)=16`,8 个 timepoint 一次铺满。这是 bam_to_tsv 阶段的真实加速(约 2→8 路,受 IO 与上游就绪节奏限制)。输出三档 md5 一致,零数值风险。
- 后续: 若要真正加速单个 bam_to_tsv 的解析吞吐(而非靠并发),属 Tier 4.2 解析循环重写。

#### 2.2 拆分 PBL/PBR 为独立 job
- 现状: bam_to_tsv、filter_aligned_reads、extract_insertion_sites 三个 rule 都在**同一 shell 块里串行跑 PBL 再跑 PBR**。
- 方案: 引入 `fragment ∈ {PBL, PBR}` wildcard,拆成独立 job → job 数翻倍,配合合理的 threads 让 snakemake 铺满核。
- 注意: bwa mapping 已在同 rule 串行 PBL/PBR;拆分需同步改下游 input 引用。风险低但触及 wildcard 结构,须 `-n` 验证 DAG。

**2.3 FastQC 从 mapping 关键路径摘除 + 合并调用**
- 现状 A: [read_processing.smk:136-139](../../workflow/rules/read_processing.smk#L136-L139) 让 `bwa_mem_mapping` 把 4 个 FastQC html 当输入 —— 人为依赖,mapping 根本不需要 QC 结果。
- 方案 A: 从 bwa mapping 的 input 里删掉 4 个 fastqc 依赖,让 FastQC 与 mapping 并行(FastQC 只需汇入 MultiQC)。
- 现状 B: [read_processing.smk:113-122](../../workflow/rules/read_processing.smk#L113-L122) 把 4 个文件拆成 4 条串行 `fastqc` 命令,`--threads 4` 形同虚设(FastQC 的 threads = 同时处理几个文件)。
- 方案 B: 合并成一条 `fastqc --threads 4 f1 f2 f3 f4` → ~4x。
- 预期: 把 FastQC 移出关键路径,mapping 提前开跑。

**2.4 samtools sort/index PBL/PBR 串行 → 拆分**
- [read_processing.smk:181-190](../../workflow/rules/read_processing.smk#L181-L190) PBL/PBR 串行且仅 2 线程。可拆分或提高线程。收益小(只喂 QC,不在 insertion 关键路径)。

#### 免费杠杆: 提高 `--cores`
- 主机 64 核,当前只用 16。把 `--cores` 提到 32/48 无需任何代码改动即可放大所有 rule 的并发度。文档记录,由用户按机器负载决定。

### Tier 3 — 中收益、中风险、脚本内向量化(改实现,需回归验证)

**3.1 extract_insertion_sites 去 `.apply(axis=1)`**
- [extract_insertion_sites.py:173](../../workflow/scripts/read_processing/extract_insertion_sites.py#L173) `calculate_insertion_coordinate` 逐行 apply。插入坐标是 `ref_start`/`ref_end` + strand 的向量运算,可用 `np.where` 向量化。
- 风险: 需逐字节比对新旧输出;有 chunk 逻辑要保持一致。

**3.2 annotate_genomic_features 去 3× apply**
- [annotate_genomic_features.py:306-325](../../workflow/scripts/read_processing/annotate_genomic_features.py#L306-L325) 三处 `.apply(axis=1)` + 一个 groupby Python for 循环。可部分向量化。
- 数据量 40MB 级,单文件几分钟,收益中等。

**3.3 filter_aligned_reads 多核 chunk**
- chunk 内已向量化(好),但 chunk 之间串行单核。可用进程池并行处理 chunk,或直接依赖 2.1 的 PBL/PBR 拆分获得 2x。
- 权衡: 引入进程池会增加复杂度;若 2.1 已拆分,单核 chunk 可能够用。

### Tier 4 — 高收益、高风险、结构性重构(单独立项)

**4.1 中间产物 TSV → parquet**
- 现状: 步骤间全部用 TSV 文本落盘再读回。`6_filtered` 单文件达 2.3GB 文本。pandas 读写文本 CSV 比 parquet/feather 慢约一个量级,且体积大。
- 方案: read-processing 链的中间产物改 parquet(列式、压缩、带 dtype)。
- 影响: 触及多个脚本的 CLI 契约(`-i/-o` 的格式假设)与 rule 的 output 后缀,以及 `chunk_size` 分块读的方式(parquet 用 row group)。属大改。
- 预期: read-processing 链 I/O 显著下降 + 磁盘占用下降。**建议作为独立设计单独评审。**

**4.2 bam_to_tsv 解析循环重写**
- [parse_bam_to_tsv.py:347-391](../../workflow/scripts/read_processing/parse_bam_to_tsv.py#L347-L391) 每 read 构造多个 frozen dataclass、逐字段格式化。`threads=8` 只作用于解压不作用于解析。
- 方案: 减少每 read 的对象分配(直接拼字段、只取下游真正用到的列);或评估是否能用 pysam 更底层 API。
- 权衡: 收益需 profile 确认(可能解压才是瓶颈);属大改,建议先做 1.2 放开并行,再决定是否值得重写。

---

## 5. 建议实施顺序与决策点

推荐分批推进,每批独立验证后再进入下一批:

1. **批 1(先做,当天可完成,零数值风险):** Tier 1.1 curve_fitting 并行化 + Tier 1.2 内存声明修正。
   - 这两项拿下最大的两块墙钟(29 min + 53 min),改动小、可逐字节验证输出一致。
2. **批 2(调度层):** Tier 2.1 PBL/PBR 拆分 + Tier 2.2 FastQC 摘路径&合并。需改 wildcard 结构,回归 `-n` dry-run 验证 DAG。
3. **批 3(脚本向量化):** Tier 3,逐脚本改 + 输出比对。
4. **批 4(重构,单独立项):** Tier 4 parquet / 解析重写。

### 决策点(已确认 2026-07-16)
- **DP1 — curve_fitting 并行后端: `joblib`。** 用 `joblib.Parallel(n_jobs=threads)` + `delayed`,依赖加进 `workflow/envs/statistics_and_figure_plotting.yml`。worker 内设 `OMP_NUM_THREADS=1` 防 numpy/scipy 过订阅。
- **DP2 — `temp()` 策略: 都不标,保持现状。** `6_filtered`(18GB)/`4_sorted`(8.3GB)继续保留,供 debug / QC 复查。§4 Tier 1.3 从实施范围中移除。
- **DP3 — Tier 4 不本轮做。** parquet 重构与 bam 解析重写单独立项。本轮范围 = 批 1 + 批 2 + 批 3(Tier 1.1/1.2 + Tier 2 + Tier 3)。

### 验证方法(每批通用)
- 数值一致性: 对改动前后的产物做 `sort` 后 `diff`(或按 key 排序后逐列比对,浮点用容差)。
- DAG 正确性: `snakemake -n --use-conda` 确认依赖图与目标不变。
- 墙钟对比: 复用本文档 §3 的日志时间戳方法,改动前后同项目对比。

---

## 6. 预期总收益(批 1 + 批 2,粗估)

| 阶段 | 当前 | 优化后(粗估) |
|---|---|---|
| bam_to_tsv(8 tp) | ~53 min 串行 | ~10 min 并行 |
| filter + extract | ~57 min | ~30 min(PBL/PBR 拆分) |
| insertion curve_fitting | ~29 min | ~2–3 min(16 核) |
| **关键路径合计** | 数量级 ~2.5 h | **有望压到 <1 h** |

数值精确加速取决于机器核数与 IO;以上为方向性估计,须以批 1 落地后的实测为准。


---

## 批 2 实施记录(2026-07-16,commit f731e42)

**2.2 PBL/PBR 拆分 —— 已实现(用户决策:引入 fragment wildcard 拆成独立 job)**
- `bam_to_tsv` / `filter_aligned_reads` / `extract_insertion_sites` 三个 rule 从"单 job 内串行跑 PBL 再 PBR"改为带 `{fragment}` 通配符(Snakefile 约束 `PBL|PBR`)的独立 job。job 数 8→16(每个 rule),两 fragment 并发。
- `merge_strand_insertions` 改为直接引用 `7_insertions/*.PBL.tsv` 与 `*.PBR.tsv` 显式路径(原 named output 已不存在)。
- **QC 契约保护(用户决策:三个都拆+日志合并规则):** filter 日志现在按 fragment 分开。新增 `merge_filter_logs` rule 把 PBL+PBR 两个日志 `cat` 回一个合并日志,保持 `extract_mapping_filtering_statistics.py` 的契约(用日志文件名 stem 做样本名、期望单文件含两块 FILTERING SUMMARY)。`mapping_filtering_statistics` 输入改指向合并日志。

**2.3 FastQC 摘除关键路径 —— 已实现(用户决策:摘依赖+合并命令+走 QC 聚合)**
- 从 `bwa_mem_mapping` 的 input 删掉 4 个 FastQC html 人为依赖 → mapping 不再被 QC 阻塞,二者并行。
- FastQC 4 条串行命令合并成 1 条 `fastqc --threads 4 file1 file2 file3 file4`(FastQC 的 threads = 同时处理几个文件,原来形同虚设)→ ~4x。
- FastQC 不孤儿:经 `multiqc_preprocessing` 可达(验证:对 multiqc target 用 `-F` 会调度 `fastqc_junction_classification` = 8 job)。

**验证:**
- `snakemake -n` 全量 clean parse;强制 QC target 时 `bam_to_tsv=16 / filter_aligned_reads=16 / merge_filter_logs=8`,符合预期。
- 单 timepoint 经真实 snakemake 重跑拆分后的 extract + merge_strand:`7_insertions` PBL/PBR 与 `8_merged` 输出**逐字节一致**于拆分前。拆分只改调度,不改数值。
- QC 解析器在合并日志上产出正确的单样本表(PBL/PBR 两列齐全、样本名无 fragment 后缀)。

**待补(下一步):** 真实全流程端到端墙钟对比(拆分前 vs 后),确认 read-processing 阶段实际压缩幅度。2.4(samtools sort/index 拆分)收益小、未做。
