
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
