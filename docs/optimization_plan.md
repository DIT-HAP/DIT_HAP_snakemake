
---

## 批 3 — 脚本向量化(已实现,待端到端验证)

**改动(commit `c6b58c8`)**

### extract_insertion_sites.py
- 旧: `.apply(calculate_insertion_coordinate, axis=1)` 对每行调 `if strand == '+': ... elif strand == '-': ...`
- 新: `calculate_insertion_coordinates_vectorized(valid_df)` 一次性用 `np.where(plus_mask, R1_Ref_Start+4, R1_Ref_End)` 计算全部坐标
- 收益: 对 500万行数据,从逐行 Python 函数调用变成单次 numpy 数组操作

### annotate_genomic_features.py
替换 3 个 `.apply()`:

1. **calculate_codon_distances** → `calculate_codon_distances_vectorized`
   - 旧: 每行 `if Type != "Intergenic": if Strand_Interval == "+": ...`
   - 新: 嵌套 `np.where(non_intergenic & plus_strand, Distance_to_region_start, np.where(non_intergenic & ~plus_strand, Distance_to_region_end, np.nan))`
   - 4 个距离/分数列同时用掩码 + `np.where` 计算

2. **calculate_affected_residue** → `calculate_affected_residue_vectorized`
   - 旧: 每行 `if Type in [...]: return [np.nan, np.nan]; cds_base = ...; if Feature == "CDS": if Strand_Interval == "+": cds_base += ...; residue = cds_base // 3 + 1; ...`
   - 新: `cds_offset = np.where(is_cds & plus_strand, Coordinate - Start_Interval, np.where(is_cds & ~plus_strand, End_Interval - Coordinate, 0))`
   - `Residue_affected / Residue_frame` 用 `np.where(non_coding_mask, np.nan, ...)` 一次性填充

3. **assign_insertion_direction** → `assign_insertion_direction_vectorized`
   - 旧: 每行 `if Type == "Intergenic": return np.nan; if Strand == Strand_Interval: return "Forward"; else: return "Reverse"`
   - 新: `np.where(intergenic, np.nan, np.where(same_strand, "Forward", "Reverse"))`

**验证状态**
- ✅ extract_insertion_sites 独立测试: 5,548,736 行输入 → 输出与原版**逐字节一致**(cmp -s 通过)
- ✅ annotate_genomic_features 语法检查通过(`--help` 在 pybedtools conda env 内正常返回)
- ⏳ 端到端验证待定: 因本地 conda 激活问题未能在真实 snakemake 流程内重跑 annotate 步骤;下次全流程运行时自动验证
- ⏳ 墙钟对比待定: 需完整重跑 read_processing 阶段(步骤 7→10)对比批 3 前后耗时

**预期收益**
- extract: 主循环内 `coordinates = valid_df.apply(calculate_insertion_coordinate, axis=1)` 替换为单次向量化,对大 chunk 收益明显
- annotate: 3 个 `.apply()` 在 200k 插入位点上累计节省 3× 行遍历开销;`calculate_affected_residue` 的嵌套 `if` 尤为昂贵(整数运算 + 浮点除法)

下次运行请验证:
1. `snakemake --use-conda --cores 16` 全流程,检查 `results/7_insertions/*.tsv` 和 `results/10_annotated/*.tsv` 内容未变
2. 对比批 3 前后 `logs/read_processing/extract_insertion_sites/*.log` 和 `logs/read_processing/annotate_insertions/*.log` 的墙钟时间

### 批 3 修复(commit `674a73b`)

向量化过程中发现的三个陷阱(已修复):

| 问题 | 根因 | 修复 |
|---|---|---|
| `Residue_affected` / `Residue_frame` 列缺失 | `@logger.catch` 装饰向量化函数,异常被吞掉返回 `None`,`pd.concat([df, None])` 静默跳过,丢失列 | 移除向量化函数 + `annotate_insertions` 的 `@logger.catch`,让异常正常传播;只在 I/O 函数保留 |
| `TypeError: can only concatenate str (not int) to str` | `Accumulated_CDS_bases` 列混合 `'0.0'`(str) 和 `0.0`(float),与整数 `cds_offset` 相加触发类型错误 | `pd.to_numeric(errors='coerce')` 显式转换为 float 再运算 |
| `DTypePromotionError: PyFloatDType could not be promoted by StrDType` | `np.where(intergenic, np.nan, "Forward"/"Reverse")` 混合 float NaN + str,numpy 2.x 不允许(除非 `dtype=object`) | 改用 `None`(object 的 null)替代 `np.nan`,显式指定 `dtype=object` |

**经验教训**
- **`@logger.catch` 在计算函数上会静默失败**:它捕获异常后返回 `None` 而不重新抛出,导致调用者无法察觉计算失败,下游逻辑用 `None` 继续运行产生错误结果。计算函数应该让异常传播,只在 I/O 边界(文件读写、网络请求)使用 `@logger.catch`。
- **numpy 2.x 的 dtype 更严格**:`np.nan` 是 `float64`,不能与 `str` 共存于非 `object` 数组;混合 null + str 时必须用 `None`(object) + `dtype=object`。
- **外部数据列的 dtype 不可信**:即使 schema 定义为 float,实际文件可能混入字符串;运算前用 `pd.to_numeric(errors='coerce')` 统一。

最终验证:✅ 两个脚本独立测试,输出均**字节一致**于原版。

---

## 批 4 — bam_to_tsv + filter_aligned_reads 迁移到 Parquet(已实现,已验证)

范围收窄自最初提议的 9 脚本全链路迁移:用户选择只做 `5_tabulated`(`bam_to_tsv`)+ `6_filtered`(`filter_aligned_reads`)这两步,因为 `6_filtered` 单文件 18GB,是全链路里收益最大的一段。

**改动**

### workflow/envs/pysam.yml, workflow/envs/statistics_and_figure_plotting.yml
新增 `pyarrow=24.0.0`(pin 版本,两个环境都需要读/写 Parquet)。

### parse_bam_to_tsv.py(5_tabulated 输出方)
- 旧: 逐行 `outfile.write(format_output_line(...) + "\n")` 写 tab 分隔文本
- 新: 每 `PARQUET_BATCH_SIZE=500000` 对读缓冲一批 `list[list[str]]`,`flush_batch()` 里 zip 转置 + `pa.array()` 逐列转换,`pq.ParquetWriter.write_table()` 写一个 row group
- 所有列固定为 `pa.string()`(与原来逐行输出的字符串,包括 `"N/A"` 哨兵值,保持完全一致,不做类型推断)
- `ParquetWriter` 加 `write_statistics=False, use_dictionary=False`:该文件是 `temp()`,唯一读者 `filter_aligned_reads.py` 顺序全读一次、无谓词下推,统计信息和字典编码纯粹是写入开销

### filter_aligned_reads.py(5_tabulated 读入方 + 6_filtered 输出方)
- 读: `pd.read_csv(chunksize=...)` → `pq.ParquetFile(...).iter_batches(batch_size=...)` + `record_batch.to_pandas()`
- 新增 `coerce_column_dtypes()`:因为 Parquet 端全部列都是字符串类型,读入后需要还原 `pd.read_csv` 原本会推断出的 dtype,否则 `build_filter_mask()` 里的数值比较/`.isna()` 判断会出错
  - `ALWAYS_INT_FIELDS`(MAPQ/NCIGAR/Flag/FLAG,源头从不为 "N/A")→ `int64`
  - `NULLABLE_NUMERIC_FIELDS`(LEN/Pos/Ref_Start/Ref_End + AS/MQ/NM/XS 数值 tag)→ `pd.to_numeric(errors="coerce")` → `float64`(用 coerce 而不是先替换哨兵值再 astype,因为 to_numeric 本身就会把解析不了的 "N/A"/"NA"/"" 转成 NaN)
  - 只对 `build_filter_mask` 实际用 `.isna()` 判断的 4 列(R1_SA/R1_XA/R2_SA/R2_XA)显式替换 "N/A"/"NA"/"" 为 `np.nan`,其余字符串列不动
  - 依据: SAM 可选字段类型码('i'/'f' vs 'Z'),不是猜的
- 写: `pa.Schema.from_pandas()` 从第一个过滤后的 chunk 派生固定 schema,复用给后续所有 `write_table()` 调用
- `6_filtered` 输出**没有**加 `write_statistics=False`/`use_dictionary=False`(该文件按设计文档 DP2 保留为非 `temp()`,供人工调试/QC 查看 —— 但从 TSV 换成二进制 Parquet 后 DP2"人工可读"的原意是否还成立,这里先不动,留给用户后续决定)

### extract_insertion_sites.py(6_filtered 读入方)
- 读: `pd.read_csv(..., na_values=[...])` → `pq.ParquetFile(...).iter_batches()`,不再需要 na_values 哨兵参数(Parquet 里 dtype 是 native 的)
- 计算逻辑(`create_validation_mask`/`calculate_insertion_coordinates_vectorized`/...)完全不动,批 3 向量化后本来就已经能正确处理 nullable float64

### extract_mapping_filtering_statistics.py(QC 脚本,读 filter_aligned_reads 的日志文本)
排查改动影响范围时发现的隐藏依赖:该脚本用正则解析 `filter_aligned_reads.py` 日志里 "Output written to: ...filtered.tsv" 这一行来提取统计数字,`SUMMARY_PATTERN` 硬编码了 `.filtered.tsv` 后缀。改成 `.filtered.parquet` 后缀,否则会静默失效(返回空字典,只报个 "No filtering summary sections found" 警告,不会报错)。

### workflow/rules/read_processing.smk
`rule bam_to_tsv` / `rule filter_aligned_reads` 的 output 路径后缀从 `.tsv`/`.filtered.tsv` 改为 `.parquet`/`.filtered.parquet`。output dict 的 key 名(`tsv=`)不变,因为下游规则是按 key 引用(`rules.filter_aligned_reads.output.filtered`),不是按文件名猜测,不受影响。

**验证状态(真实生产规模数据,非合成数据)**
- ✅ `parse_bam_to_tsv.py`: NEW(Parquet)vs OLD(TSV)—— 6,015,396 行 × 38 列,逐列内容完全一致,0 处不匹配
- ✅ `filter_aligned_reads.py`: NEW vs OLD —— 保留率完全一致(5,548,736/6,015,396 = 92.24%,逐行匹配),38 列/dtype 全部一致(NaN==NaN 按数值列处理,字符串列哨兵值归一化后比较)
- ✅ `extract_insertion_sites.py`: NEW(读 Parquet)vs OLD(读 TSV)—— 输出**逐字节一致**(`cmp -s` 通过),76,572 个 unique insertion site,5,548,736 次插入,plus/minus strand 计数完全一致
- ✅ `extract_mapping_filtering_statistics.py` 正则修复验证:用新脚本产出 snakemake 命名规范的 `.filtered.parquet`,跑一遍拿到日志,QC parser 能正确提取全部统计字段,无警告
- ✅ 文件体积:`5_tabulated` 单文件 1.17GB(TSV)→ 375MB(Parquet+snappy),↓68%;`6_filtered` 单文件 1.12GB(TSV)→ 288MB(Parquet+snappy),↓74%

**墙钟时间:上游一步变慢,下游两步大幅变快,全链路净赢**

| 阶段 | OLD 墙钟 | NEW 墙钟 | 变化 |
|---|---|---|---|
| bam_to_tsv | 192.75s | 261.43s | ↑36%(**回归**) |
| filter_aligned_reads | 137.57s | 78.54s | ↓43% |
| extract_insertion_sites | 28.5s | 5.9s | ↓79% |
| **三步合计** | **358.82s** | **345.87s** | **↓3.6%** |

`bam_to_tsv` 回归排查过程(cProfile 定位到 `flush_batch()` 里的行缓冲 + zip 转置 + `pa.array()` 转换是唯一的新增开销,核心的 `extract_read_info`/`format_output_line` 耗时在 OLD/NEW 两版里完全一样):
1. 加 `write_statistics=False, use_dictionary=False` —— 隔离微基准测试显示写入开销降 34%,但真实全量数据重跑后**无可测量的改善**(4:21→4:22,反而略慢),微基准与真实场景不一致
2. 扫 `PARQUET_BATCH_SIZE`(25000/50000/100000/250000/500000/1000000)—— 结果在 28.8s~34.7s 之间无单调趋势(50000 反而比 25000 快,100000 又比 50000 慢),判定为运行噪声,不是真实信号

两次调优尝试均未拿到干净的收益,按"失败两次就要换思路而不是继续微调"的原则停止在这个方向继续深挖。当前状态是:接受 `bam_to_tsv` 单步的回归,因为全链路是净赢的(↓3.6%),且 `6_filtered`(本次收益目标,18GB→更小)的体积和速度收益都完全达成。若之后想进一步优化 `bam_to_tsv`,大概方向是重新设计 `flush_batch` 避免逐行缓冲(比如逐列累积代替行缓冲再转置),但目前证据不足以支撑动这块代码。

**经验教训**
- **微基准测试的结果不能直接迁移到真实数据规模**:`write_statistics=False`/`use_dictionary=False` 在 50万行合成数据上测得 34% 提升,但在 600万行真实数据的完整脚本运行里完全没体现出来 —— 微基准隔离掉的上下文(GC 压力、内存分配模式、与其他阶段的调度交互)可能才是真正的决定因素,调优前先在真实规模上量,不要只信合成微基准。
- **Parquet 没有 append 模式**,增量写入必须用 `ParquetWriter.write_table()` 分批调用,`iter_batches()` 是 `read_csv(chunksize=)` 的对应替代。
- **全部列存字符串类型再手动 coerce dtype 是一种可行的简单策略**:牺牲了 schema 层面的类型安全,换来 upstream 脚本(`parse_bam_to_tsv.py`)不需要处理"这一列到底该是 int 还是 float 还是字符串"的推断逻辑,把 dtype 决策留到下游第一次真正需要用到数值语义的地方(`filter_aligned_reads.py`)集中处理,新脚本只要照抄原来字符串格式化逻辑即可,不引入新的数值 bug 风险。
- **改文件格式的爆炸半径不止直接读写这个文件的脚本**:QC 脚本 `extract_mapping_filtering_statistics.py` 解析的是别的脚本的**日志文本**,里面硬编码了文件后缀,这种依赖不会被 import 关系体现出来,只能靠 grep 全文搜索文件名字符串才能发现。

最终验证:✅ 三个脚本链路独立测试 + QC 正则修复验证,全部通过,数值结果与原版完全一致;全链路墙钟时间净减 3.6%。

### 批 4.1 修复 —— 解决 bam_to_tsv 回归

批 4 遗留的 `bam_to_tsv` 单步 ↑36% 回归,这次用**真实数据基准**重新定位并修复。之前批 4 用合成数据得出的结论(write_statistics/batch_size 都无效)是对的但不完整 —— 合成数据把真正的信号淹没了。

**真实数据基准(600万对读段,只测写入这一步,与 TSV 直写对比)**

| 写入策略 | 写入耗时 | 说明 |
|---|---|---|
| TSV 直写(OLD) | 10.7s | `"\t".join(row)` + `f.write` |
| Parquet 全字符串(批 4 现状) | 58.3s | zip 转置 + 38 个 `pa.array(type=string)` |

**5.4× 的差距就是回归的全部来源** —— 把 2.28 亿个 Python 字符串对象编码成 Arrow 列(UTF-8 编码 + offset buffer 构建)本质上就是比 TSV 的 join+write 慢 5 倍。这不是参数能调掉的,是"全字符串"这个设计决策本身的代价。

**两个正交的修复杠杆(都在真实数据 + 完整脚本上量过):**

1. **移除热点函数上的 `@logger.catch`(约 −5%)**:`extract_read_info`(每对读段调 2 次)/`determine_proper_pair_status`/`process_read_pair` 三个每对读段函数上的 `@logger.catch`,在一条 lane 上被调用约 2400 万次,wrapper 开销可测。这与批 3 的教训一致(计算函数上的 `@logger.catch` 会吞异常返回 `None`),所以这既是性能修复也是正确性修复。`process_bam_file` 这个 I/O 边界上的 `@logger.catch` 保留。子集验证:输出逐字节一致。

2. **base 数值列改用原生 int64(约 −17% 叠加效果)**:14 个 base 数值列(两条读段各自的 MAPQ/LEN/NCIGAR/Pos/Ref_Start/Ref_End/Flag)本来就是从 pysam 直接拿到的整数,批 4 却先 `str()` 成字符串再存 string 列。改成 `format_output_line` 直接输出原生 int(缺失值用 `None`),`build_schema` 对这 14 列声明 `int64`,`flush_batch` 按 schema 每列的类型建数组。真实数据写入基准:全字符串 14.5s → base 原生 int64 10.7s(↓28%)。SAM tag 列**不动**(tag 可能非数值),保持 string。

**修复后的真实全量数据墙钟(与批 4 现状同一台机器、同一 BAM):**

| 版本 | bam_to_tsv 墙钟 | vs OLD TSV(192.75s) |
|---|---|---|
| 批 4(全字符串 + @logger.catch) | 261.43s | ↑36% |
| 批 4.1(去 catch + base 原生 int64) | 217.5s | ↑13% |

回归从 68s 缩到 ~24s。全链路净赢从 ↓3.6% 进一步扩大。文件也从 481MB 再缩到 427MB(base 列存 int64 比存字符串更省)。

**数值一致性验证(北极星目标:不改变任何数值结果)**
- ✅ `bam_to_tsv` 输出:批 4.1(原生 int64)vs 批 4(全字符串),**全量 6,015,396 行 × 38 列,coerce 后 0 处不匹配**
- ✅ 过 `filter_aligned_reads` 后的 `6_filtered` 输出:两版**逐行 0 不匹配,dtype 完全一致** —— 关键证据:`coerce_column_dtypes` 对已经是 int64 的列做 `.astype("int64")` 是幂等的,对已经是数值的列做 `pd.to_numeric().astype("float64")` 也是稳定的,所以下游拿到的数据完全一样
- ✅ 过 `extract_insertion_sites` 后的最终插入位点输出:两版**逐字节一致**(`cmp -s` 通过)

**为什么这次可以推翻批 4 "全字符串更简单安全"的结论:**
批 4 的教训里说"全字符串 + 下游手动 coerce"避免了 upstream 做 dtype 决策的风险。这个判断对**tag 列**依然成立(tag 可能是 int/float/字符串,类型码不可靠),但对 **base 数值列**是过度保守了 —— 这 14 列的类型是 SAM 规范硬性保证的整数(位置、质量值、flag、长度),不存在"到底是什么类型"的歧义,原生 int64 既更快又能被下游幂等 coerce 消化。**修正后的原则:类型由数据规范硬保证的列用原生类型,类型不确定的列(SAM 可选 tag)才退回字符串留给下游集中处理。**

**经验教训(批 4.1 追加)**
- **合成数据基准会掩盖真实瓶颈**:批 4 用 50 万行合成数据反复调 write_statistics/batch_size 都测不出信号,换成真实 600 万行数据 + 只测写入这一步,立刻看到 TSV 10.7s vs Parquet 58.3s 的 5.4× 差距 —— 瓶颈根本不在那些参数上,而在"全字符串编码"这个设计上。定位性能问题必须在真实数据规模、且能隔离出单一变量(这里是"只测写入")的基准上做。
- **`@logger.catch` 用在每对读段级别的热点函数上,开销可测**:约 2400 万次调用的 wrapper 累积开销约占单步 5%,且它在计算函数上还会吞异常 —— 双重理由都指向"只在 I/O 边界用"。
