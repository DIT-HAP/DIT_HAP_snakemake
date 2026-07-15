# DIT-HAP 上游流程重构设计

**日期**: 2026-07-09
**范围**: 仅 `/data/c/yangyusheng_optimized/DIT_HAP`（`DIT_HAP_pipeline` 不动）
**状态**: 已确认，待实施

---

## 1. 背景与目标

原项目 `DIT_HAP_pipeline` 将上游 Snakemake 流程与下游分析（enrichment / clustering / 机器学习 / thesis figures 等 notebooks 及配套 `src` 模块）混在一个仓库里。本次重构把**上游**抽到干净的新仓库 `DIT_HAP`，并对其代码结构与风格做现代化优化。

### 分离边界

- **DIT_HAP（上游，本仓库）**：从原始测序 reads 一路跑到 **gene-level DR/DL 表格**，外加 QC 报告。
- **DIT_HAP_pipeline（下游，保持不动）**：其 notebooks 消费上游产出的 gene-level 表格，做后续生物学分析。
- **交接接口**：gene-level DR/DL 表格（`results/{project_name}/16_gene_level_depletion_analysis/` 与 `17_gene_level_curve_fitting/`）。

分离本身在 DIT_HAP 中已基本完成（下游 notebooks 与 `src` 模块未带入）。本次聚焦于清理「重组做到一半」留下的连接点，并把上游代码提升到统一的现代标准。

### 四项工作

1. 命名与目录结构统一
2. `run:` 块归位 + 孤儿脚本接线
3. 脚本逐个彻底重构到 `python-script-conventions` skill 新标准
4. 文档 / 配置清理与收尾

---

## 2. 第 1 段：命名与目录结构统一

确立四个规范模块名，全线贯穿（目录名 / 规则 include / log 路径 / message）：

| 模块 | 规则文件 | 脚本目录 | log 目录 |
|------|---------|---------|---------|
| `reference_data` | reference_data.smk | scripts/reference_data/ | logs/{project_name}/reference_data/ |
| `read_processing` | read_processing.smk | scripts/read_processing/ | logs/{project_name}/read_processing/ |
| `depletion_scoring` | depletion_scoring.smk | scripts/depletion_scoring/ | logs/{project_name}/depletion_scoring/ |
| `quality_control` | quality_control.smk | scripts/quality_control/ | logs/{project_name}/quality_control/ |

**改动清单**：

1. **合并 reference 目录**：`scripts/reference/fetch_pombase_datasets.sh` → `scripts/reference_data/`；删空目录 `scripts/reference/`；reference_data.smk 中 shell 路径同步更新。
2. **修脚本引用路径**：depletion_scoring.smk 中所有 `workflow/scripts/depletion_analysis/` → `workflow/scripts/depletion_scoring/`（当前指向不存在的目录，会导致规则失败）。
3. **统一 log 路径**：
   - `logs/{project_name}/preprocessing/` → `logs/{project_name}/read_processing/`
   - `logs/{project_name}/depletion_analysis/` → `logs/{project_name}/depletion_scoring/`
   - reference_data.smk 中 `logs/preparation/...` → `logs/{project_name}/reference_data/...`（补上缺失的 `{project_name}` 层级）
4. **results 步骤编号保留不动**（01–17）。

---

## 3. 第 2 段：`run:` 块归位 + 孤儿脚本接线

三个内联 `run:` 块替换为已存在但尚未接线的独立脚本。

| 规则位置 | 当前 `run:` 块 | 目标脚本（已存在） |
|---------|--------------|------------------|
| read_processing.smk:426 | merge_similar_timepoints | `read_processing/merge_similar_timepoints.py` |
| read_processing.smk:459 | concat_counts_and_annotations | `read_processing/concat_counts_and_annotations.py` |
| depletion_scoring.smk:185 | r_square_as_weights | `depletion_scoring/compute_r2_weights.py` |

**改动清单**：

1. 三个 `run:` 块 → `shell:` 调用对应脚本，参数经命令行（`-i/-o` 等）传入，与其余规则风格一致。
2. **接线前逐个核对**脚本的输入 / 输出 / 参数签名与 `run:` 块逻辑是否一致；`concat_counts_and_annotations` 需确认 argparse 能接收多文件列表输入（`input.counts` / `input.annotations` 各为多文件）。
3. **逻辑冲突时以脚本（改进版）为准**。
4. 接线后规则中不再有任何 `run:` 块，全部走 `shell:` + 脚本。

**收益**：`run:` 块无法用 `conda:` 隔离、难以单独测试；抽成脚本后可独立运行与验证。

---

## 4. 第 3 段：脚本逐个彻底重构（24 个脚本 ≈ 8000 行）

**标准来源**：`python-script-conventions` skill（Modern Python 3.12+）。
原 `PYTHON_STYLE_GUIDE.md` 作废。

### 目标标准要点

- **7 段分区**：IMPORTS → DECORATORS → CONSTANTS → CONFIG → LOGGING → CORE → MAIN
- **类型系统**：`X | None`（非 `Optional`）、native `list`/`dict`、`type` 别名、`def func[T]`
- **配置对象**：`@dataclass(kw_only=True, slots=True, frozen=True)`
- **枚举**：`StrEnum` 用于成组字符串标记（列名 / 模式 / 状态），不用于单例常量
- **控制流**：`match...case` 替代长 if-elif；guard clauses 处理边界
- **日志**：loguru（禁 print）、`@logger.catch` 于核心函数、`--verbose` flag
- **路径**：`pathlib.Path`（禁 `os.path`）
- **入口**：`parse_args()` + `if __name__ == "__main__": sys.exit(main())`，失败返回 1
- **docstring**：模块级按 §3.2 格式（title → description → Input → Output → Usage → metadata）；函数级单行
- **库模块 vs standalone**：被 import 的库模块省略 LOGGING/CONFIG/MAIN 段
- **文档语言**：公开仓库用英文

### 重构铁律（风险控制）

- **只改结构与风格，不改算法逻辑。**
- 每个脚本重构后必须**验证输出不变**（用已有中间产物或真实数据对照）。
- 每个脚本（或一小批）独立 commit，便于回溯。

### 已知明确缺口（优先项）

- `quality_control/gene_coverage_analysis.py` — 5 处 `print()` → loguru
- `depletion_scoring/gene_level_depletion_analysis.py` — 缺 section 分区
- `read_processing/reads_hard_filtering.py` — 缺 section 分区

### 待重构脚本清单（按模块）

- **reference_data/**：extract_genome_region.py
- **read_processing/**：parse_bam_to_tsv.py, filter_aligned_reads.py, extract_insertion_sites.py, merge_strand_insertions.py, concatenate_timepoint_data.py, annotate_genomic_features.py, merge_similar_timepoints.py, concat_counts_and_annotations.py, reads_hard_filtering.py
- **depletion_scoring/**：def_ctr_insertions.py, impute_missing_values_using_FR.py, insertion_level_depletion_analysis_has_replicates.py, insertion_level_depletion_analysis_no_replicates.py, curve_fitting.py, compute_r2_weights.py, gene_level_depletion_analysis.py
- **quality_control/**：extract_mapping_filtering_statistics.py, PBL_PBR_correlation_analysis.py, read_count_distribution_analysis.py, insertion_orientation_analysis.py, insertion_density_analysis.py, gene_coverage_analysis.py, distribution_of_curve_fitting_results.py

---

## 5. 第 4 段：文档 / 配置清理与收尾

### 配置

1. **修 config workdir**（重要）：所有 `config/*.yaml` 的 `workdir:` 由 `DIT_HAP_pipeline` → `DIT_HAP`（Snakefile 的 `workdir:` 已正确）。
2. **删下游 config 字段**：经扫描确认，**无下游专用残留字段**，此项无需执行。

### 文档（DIT_HAP 当前无 README/CLAUDE.md，按需新建）

3. **README.md**（新建）：只讲上游到 gene-level 表格 + QC；注明「下游分析见 DIT_HAP_pipeline」。英文。
4. **CLAUDE.md**（新建）：新模块名、正确的 output structure、无下游 notebooks 段落；Python 风格指向 `python-script-conventions` skill。
5. **openspec**：不引入。
6. **resources/**：上游仅用 PomBase 下载数据与 `Hayles_2013` 表；下游用的一堆 xlsx/csv 本就未带入，确认即可。

### 验证关卡

- 全程 `snakemake -n`（dry-run）+ `snakemake --lint`。
- 每阶段独立 commit。

---

## 6. 实施节奏

先落盘本设计文档并 commit，再用 writing-plans 生成详细实施步骤清单，然后逐步执行。

四段按序推进：第 1+2 段（让上游 dry-run 跑通）→ 第 3 段（脚本重构）→ 第 4 段（文档）。
