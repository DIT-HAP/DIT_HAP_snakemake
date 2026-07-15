# Config & Sample Sheet Schema 规范化设计

**日期**: 2026-07-09
**范围**: 仅 `/data/c/yangyusheng_optimized/DIT_HAP`
**状态**: 已确认，待实施
**前置**: 上游重构（分支 refactor/upstream-cleanup）已完成

---

## 1. 目标

给 config YAML 和 sample sheet TSV 引入 Snakemake 官方标准的 JSON Schema 校验（`snakemake.utils.validate`），实现：
- 加载时立即校验（而非跑到一半才因缺键/类型错炸）
- schema 内嵌每个参数的 `description`（含单位/取值范围/示例）作为权威文档
- sample sheet 的列结构、类型、必填得到约束

## 2. 已确认的设计决策

- **方式**：JSON Schema + `snakemake.utils.validate()`（官方标准）。
- **条件必填（决策 A）**：用 JSON Schema `if/then` —— 当 `merge_similar_timepoints: true` 时才 require `similar_timepoints / merged_timepoint / drop_columns`。
- **默认值（决策 B）**：schema 给可选参数标 `default` + `description`（起文档作用，`validate` 不回填），代码保持 `config.get(...)` 兜底。仅 `use_DEseq2_for_biological_replicates` 属此类（default: false）。
- **约束严格度**：严格 —— 类型 + 数值范围（minimum/maximum）+ 字符串 pattern（adapter 序列 `^[ACGTNacgtn]+$`、release_version 日期 `^\d{4}-\d{2}-\d{2}$`、fastq 路径后缀 `\.(fq|fastq)\.gz$`）。

## 3. 目录结构

```
workflow/schemas/
├── config.schema.yaml    # 校验 config YAML（顶层 object）
└── samples.schema.yaml   # 校验 sample sheet 每一行
```

## 4. config.schema.yaml

JSON Schema draft-07（YAML 写法），`type: object`。

**必填核心键（17 个，9/9 全有）**：workdir, snakemake_wrapper_version, sample_sheet, project_name, multiqc_config, merge_similar_timepoints, initial_time_point, hard_filtering_cutoff, chunk_size, aligned_read_filtering, adapter_sequence, adapter_sequence_r2, Pombase_release_version, PBL_adapter, PBL_reverseComplement_adapter, PBR_adapter, PBR_reverseComplement_adapter。

**可选带 default**：`use_DEseq2_for_biological_replicates` (boolean, default false)。

**曲线拟合**：`time_points`（array of number）—— 必填（8/9，spikein 缺；spikein 不跑曲线拟合，但为一致性仍建议要求或单独说明。实施时确认：若 spikein 确实不需要，改为非 required 并在 description 注明）。

**条件必填（if/then）**：
```yaml
if:
  properties: {merge_similar_timepoints: {const: true}}
then:
  required: [similar_timepoints, merged_timepoint, drop_columns]
```

**嵌套 object**：`aligned_read_filtering` 定义为嵌套 schema —— `read_1_filtering` / `read_2_filtering` 各含 mapq_threshold(0-60)、nm_threshold(≥0)、ncigar_value(≥0)、no_sa(bool)、no_xa(bool)；顶层 require_proper_pair(bool)。

**约束示例**：
- `hard_filtering_cutoff`: integer, minimum 0
- `chunk_size`: integer, minimum 1
- `mapq_threshold`: integer, 0–60
- adapter 系列: string, pattern `^[ACGTNacgtn]+$`
- `Pombase_release_version`: string, pattern `^\d{4}-\d{2}-\d{2}$`
- `time_points`: array, items number
- `similar_timepoints` / `drop_columns`: array of string
- 每个参数都带 `description`（单位/取值范围/示例/在哪个规则用到）

## 5. samples.schema.yaml

逐行校验（`validate(df, schema)` 按 properties 校验每行）：
```yaml
$schema: "http://json-schema.org/draft-07/schema#"
description: "Row schema for the DIT-HAP sample sheet TSV"
type: object
properties:
  Sample:     {type: string, description: "Sample/library ID; groups timepoints of one library"}
  Timepoint:  {type: string, description: "Timepoint label (e.g. YES0, YES1); a column in the count matrix"}
  Condition:  {type: string, description: "Experimental condition (e.g. YES)"}
  read1:      {type: string, pattern: "\\.(fq|fastq)\\.gz$", description: "Absolute path to R1 fastq.gz"}
  read2:      {type: string, pattern: "\\.(fq|fastq)\\.gz$", description: "Absolute path to R2 fastq.gz"}
required: [Sample, Timepoint, Condition, read1, read2]
```

**关键**：Snakefile 用 `pd.read_csv(config["sample_sheet"], sep="\t", dtype=str)` 读，保证列都是 string（避免 Timepoint 形如 `0` 被读成 int 导致校验失败）。不校验路径存在（那由 Snakemake DAG 处理）。

## 6. Snakefile 集成

```python
from snakemake.utils import validate

configfile: snakemake_config_file
validate(config, "workflow/schemas/config.schema.yaml")   # 加载即校验

sample_sheet = pd.read_csv(config["sample_sheet"], sep="\t", dtype=str)
validate(sample_sheet, "workflow/schemas/samples.schema.yaml")  # 逐行校验
```

## 7. 验证

- 对全部 9 个 config + 6 个 sample sheet 各跑 `snakemake -n`，确认通过校验、DAG 正常。
- 负向测试：删一个必填键 / 改错一个类型，确认 schema 拦下并给出清晰报错。
- 确认 `dtype=str` 读入不破坏既有 rule（wildcard 仍是 string，一致）。

## 8. 提交

- 一个 commit 加两个 schema + Snakefile 集成 + `dtype=str` 改动。
- 分支：继续用 refactor/upstream-cleanup（或按用户指示新开分支）。
