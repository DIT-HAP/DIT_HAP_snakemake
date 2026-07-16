
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
