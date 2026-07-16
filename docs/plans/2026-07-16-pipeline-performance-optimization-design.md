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
- **实测: 1729s → 136s(-j 16),12.7× 加速**,RSS ~675MB。子集 3000 行 serial(-j 1)vs parallel(-j 8)三个输出文件全部 `diff` 一致,success rate 99.975% 与原版一致。

**1.2 修正 bam_to_tsv 的内存声明 —— 最廉价的墙钟收益 ✅ 已实现(2026-07-16)**
- 现状: [read_processing.smk:259](../../workflow/rules/read_processing.smk#L259) 声明 `mem_mb=200000`(200GB)。但该脚本是 flat-memory 流式解析(逐 read,不累积),实际只吃几百 MB–2GB。
- 后果: Snakemake 以为每个 job 要 200GB,在 `--cores 16` 下**几乎无法同时跑两个** → 8 个 timepoint 被迫串行。
- 实现: 改成 `mem_mb=4000`(留余量),附注释说明理由。
- 预期: 8 路 bam_to_tsv 可并行 → 该阶段墙钟从 ~53 min 压到 ~10 min(受核数/IO限制)。

**1.3 ~~给 `6_filtered` / `4_sorted` 加 `temp()`~~ —— 已决定不做(DP2)**
- 保持现状,`6_filtered`(18GB)/`4_sorted`(8.3GB)保留供 debug / QC 复查。此项不纳入实施范围。

### Tier 2 — 中收益、低风险、调度/依赖层(不改算法)

**2.1 拆分 PBL/PBR 为独立 job**
- 现状: bam_to_tsv、filter_aligned_reads、extract_insertion_sites 三个 rule 都在**同一个 shell 块里串行跑 PBL 再跑 PBR**。
- 方案: 引入 `fragment ∈ {PBL, PBR}` wildcard,拆成两个独立 job,让 Snakemake 并行调度。
- 预期: 每个 rule ~2x(在有空闲核时)。风险低,但要改 rule 的 wildcard 结构与下游 input 引用。

**2.2 FastQC 从 mapping 关键路径摘除 + 合并调用**
- 现状 A: [read_processing.smk:136-139](../../workflow/rules/read_processing.smk#L136-L139) 让 `bwa_mem_mapping` 把 4 个 FastQC html 当输入 —— 人为依赖,mapping 根本不需要 QC 结果。
- 方案 A: 从 bwa mapping 的 input 里删掉 4 个 fastqc 依赖,让 FastQC 与 mapping 并行(FastQC 只需汇入 MultiQC)。
- 现状 B: [read_processing.smk:113-122](../../workflow/rules/read_processing.smk#L113-L122) 把 4 个文件拆成 4 条串行 `fastqc` 命令,`--threads 4` 形同虚设(FastQC 的 threads = 同时处理几个文件)。
- 方案 B: 合并成一条 `fastqc --threads 4 f1 f2 f3 f4` → ~4x。
- 预期: 把 FastQC 移出关键路径,mapping 提前开跑。

**2.3 samtools sort/index PBL/PBR 串行 → 拆分**
- [read_processing.smk:181-190](../../workflow/rules/read_processing.smk#L181-L190) PBL/PBR 串行且仅 2 线程。可拆分或提高线程。收益小(只喂 QC,不在 insertion 关键路径)。

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

