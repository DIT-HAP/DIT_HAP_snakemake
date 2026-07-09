# DIT-HAP 上游流程重构 Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 把 DIT_HAP 上游流程（跑到 gene-level DR/DL 表格 + QC）清理成命名统一、无内联 `run:` 块、脚本全部符合 python-script-conventions skill 新标准的干净仓库。

**Architecture:** 分四段推进——(1) 命名/目录结构统一，(2) `run:` 块归位+孤儿脚本接线，(3) 24 个脚本逐个重构到 Modern Python 3.12+ 标准，(4) 文档/配置清理。第 1+2 段让上游 `snakemake -n` 跑通；第 3 段每个脚本重构后用"行为不变"验证；第 4 段收尾。

**Tech Stack:** Snakemake 8+, Python 3.12, loguru, pandas, pydeseq2, pybedtools；风格标准来自 `python-script-conventions` skill。

**关键参考路径:**
- 风格标准：`python-script-conventions` skill（用 Skill 工具加载）+ 其 `templates/agent_script_template.py`
- 黄金对照数据：`/data/c/yangyusheng_optimized/DIT_HAP_pipeline/results/HD_DIT_HAP_generationRAW/`（旧脚本产物，仅供人工核对数值合理性；编号命名与新版不同）
- 设计文档：`docs/plans/2026-07-09-upstream-refactor-design.md`

**全局铁律:**
- DIT_HAP_pipeline 目录绝对不动。
- 脚本重构只改结构/风格，不改算法逻辑。
- 每段（或每脚本）独立 commit。
- 每次改完规则跑 `snakemake -n` 确认 DAG 不报错。

---

## 阶段 A：命名与目录结构统一

### Task A1: 合并 reference 目录

**Files:**
- Move: `workflow/scripts/reference/fetch_pombase_datasets.sh` → `workflow/scripts/reference_data/fetch_pombase_datasets.sh`
- Modify: `workflow/rules/reference_data.smk:23`（shell 路径）

**Step 1:** 用 `git mv` 移动脚本（保留可执行位与历史）：
```bash
git mv workflow/scripts/reference/fetch_pombase_datasets.sh workflow/scripts/reference_data/fetch_pombase_datasets.sh
rmdir workflow/scripts/reference 2>/dev/null || true
```

**Step 2:** 修 reference_data.smk 中 `download_pombase_data` 规则的 shell 行：
`workflow/scripts/reference/fetch_pombase_datasets.sh` → `workflow/scripts/reference_data/fetch_pombase_datasets.sh`

**Step 3:** 验证目录已空且引用已更新：
```bash
ls workflow/scripts/reference 2>&1   # 期望: No such file
grep -rn "scripts/reference/" workflow/   # 期望: 无输出
```

**Step 4:** 暂不 commit（与 A2/A3 合并为一个"命名统一"commit）。

---

### Task A2: 修 depletion_scoring 脚本引用路径

**Files:**
- Modify: `workflow/rules/depletion_scoring.smk`（所有 shell 中的脚本目录）

**Step 1:** 全量替换该文件中 `workflow/scripts/depletion_analysis/` → `workflow/scripts/depletion_scoring/`（6 处：def_ctr_insertions, impute_missing_values_using_FR, insertion_level_depletion_analysis_has_replicates, insertion_level_depletion_analysis_no_replicates, curve_fitting, gene_level_depletion_analysis）。

**Step 2:** 验证：
```bash
grep -rn "scripts/depletion_analysis" workflow/   # 期望: 无输出
grep -rn "scripts/depletion_scoring" workflow/rules/depletion_scoring.smk | wc -l   # 期望 >=6
```

**Step 3:** 逐个确认引用的脚本真实存在：
```bash
for s in def_ctr_insertions impute_missing_values_using_FR insertion_level_depletion_analysis_has_replicates insertion_level_depletion_analysis_no_replicates curve_fitting gene_level_depletion_analysis; do test -f "workflow/scripts/depletion_scoring/$s.py" && echo "OK $s" || echo "MISSING $s"; done
```
期望：全部 OK。

---

### Task A3: 统一 log 目录命名

**Files:**
- Modify: `workflow/rules/read_processing.smk`
- Modify: `workflow/rules/depletion_scoring.smk`
- Modify: `workflow/rules/reference_data.smk`

**Step 1:** read_processing.smk 中 `logs/{project_name}/preprocessing/` → `logs/{project_name}/read_processing/`（全量替换）。

**Step 2:** depletion_scoring.smk 中 `logs/{project_name}/depletion_analysis/` → `logs/{project_name}/depletion_scoring/`（全量替换）。

**Step 3:** reference_data.smk 中 `logs/preparation/` → `logs/{project_name}/reference_data/`。注意该文件的规则当前 log 无 `{project_name}` 层级也无 wildcard——需确认这些规则的 wildcard 上下文（`{release_version}`）。因这些规则不带 `project_name` wildcard，log 路径改为 `logs/reference_data/...{release_version}.log`（保持无 project_name，与规则 wildcard 一致），或统一加 project_name（需评估是否引入未定义 wildcard）。**实施时先确认 reference_data 规则是否能访问 project_name 变量**（Snakefile 中 `project_name` 是全局 Python 变量，可用 f-string 注入）。采用 f-string 注入：`f"logs/{project_name}/reference_data/download_pombase_data_{{release_version}}.log"`。

**Step 4:** 验证 dry-run 不报错：
```bash
snakemake -n 2>&1 | tail -20   # 期望: 无 error，DAG 构建成功
grep -rn "preprocessing\|depletion_analysis\|preparation" workflow/rules/*.smk   # 期望: 仅剩注释或 message 文案，无 log 路径
```

**Step 5: Commit 阶段 A**
```bash
git add workflow/
git commit -m "refactor: unify module naming across dirs, script paths, and log paths"
```

---

## 阶段 B：`run:` 块归位 + 孤儿脚本接线

> 每个 Task 先核对脚本签名 vs run 块逻辑，冲突以脚本为准。

### Task B1: merge_similar_timepoints 接线

**Files:**
- Read: `workflow/scripts/read_processing/merge_similar_timepoints.py`（确认 argparse 签名）
- Modify: `workflow/rules/read_processing.smk`（rule merge_similar_timepoints 的 `run:` → `shell:`）

**Step 1:** 读脚本，记录它的 CLI 参数（input/output/similar_timepoints/merged_timepoint/drop_columns 如何传入）。

**Step 2:** 对照当前 `run:` 块逻辑（read_csv index_col=[0,1,2,3] → 求和 → drop → sort → to_csv），确认脚本等价或为改进版。若脚本缺少某参数入口，补脚本的 argparse（属"接线"必要改动，不算逻辑改动）。

**Step 3:** 把 `run:` 块替换为 `shell:` 调用，params 经 CLI 传入。列表型参数（similar_timepoints、drop_columns）用空格拼接或多值传入，与脚本 argparse 的 `nargs` 匹配。

**Step 4:** dry-run 验证：
```bash
snakemake -n results/HD_DIT_HAP_generationRAW/11_merged/... 2>&1 | tail   # 用一个具体 merged 目标
```

**Step 5:** Commit：
```bash
git add workflow/rules/read_processing.smk workflow/scripts/read_processing/merge_similar_timepoints.py
git commit -m "refactor: wire merge_similar_timepoints rule to script"
```

---

### Task B2: concat_counts_and_annotations 接线

**Files:**
- Read: `workflow/scripts/read_processing/concat_counts_and_annotations.py`
- Modify: `workflow/rules/read_processing.smk`（rule concat_counts_and_annotations）

**Step 1:** 读脚本，重点确认 argparse 能否接收**多文件列表**输入（counts 多文件 + annotations 多文件），以及两个输出（counts + annotations）如何传参。

**Step 2:** 对照 `run:` 块：counts 按 `Path(f).name.split(".")[0]` 做 key 后 concat 成双层列 MultiIndex；annotations concat 后 reset→drop_duplicates→set_index。确认脚本逻辑等价（冲突以脚本为准）。

**Step 3:** `run:` → `shell:`。多文件用 `{input.counts}`（Snakemake 展开为空格分隔）传给脚本 `nargs="+"` 参数。

**Step 4:** dry-run：
```bash
snakemake -n results/HD_DIT_HAP_generationRAW/12_concatenated/raw_reads.tsv 2>&1 | tail
```

**Step 5:** Commit：
```bash
git add workflow/rules/read_processing.smk workflow/scripts/read_processing/concat_counts_and_annotations.py
git commit -m "refactor: wire concat_counts_and_annotations rule to script"
```

---

### Task B3: r_square_as_weights 接线

**Files:**
- Read: `workflow/scripts/depletion_scoring/compute_r2_weights.py`
- Modify: `workflow/rules/depletion_scoring.smk`（rule r_square_as_weights）

**Step 1:** 读脚本，确认 CLI 签名（input 拟合统计 → output weights）。

**Step 2:** 对照 `run:` 块：R2 clip → confidence → 找 `_fitted` 列 → weights = 1-confidence 广播到各 timepoint 列。**注意 run 块用了 `col.rstrip("_fitted")`（rstrip 是字符集删除，有潜在 bug）**——若脚本已修正为 `removesuffix`，以脚本为准。

**Step 3:** `run:` → `shell:`。

**Step 4:** dry-run：
```bash
snakemake -n results/HD_DIT_HAP_generationRAW/15_insertion_level_curve_fitting/insertions_LFC_fitted_with_r_square_as_weights.tsv 2>&1 | tail
```

**Step 5:** Commit：
```bash
git add workflow/rules/depletion_scoring.smk workflow/scripts/depletion_scoring/compute_r2_weights.py
git commit -m "refactor: wire r_square_as_weights rule to compute_r2_weights.py script"
```

---

### Task B4: 确认无残留 run 块 + 全局 dry-run

**Step 1:**
```bash
grep -rn "run:" workflow/rules/   # 期望: 无输出
```

**Step 2:** 逐个 config 做 dry-run（覆盖 replicate 与 no-replicate 两个分支）：
```bash
snakemake -n --configfile config/config_HD_generationRAW.yaml 2>&1 | tail -5
snakemake -n --configfile config/config_LD_generationRAW.yaml 2>&1 | tail -5
```
期望：两者 DAG 均构建成功。

**Step 3:** `snakemake --lint 2>&1 | tail -30`，记录 lint 警告（不强制清零，但记录）。

---

## 阶段 C：脚本逐个重构到 python-script-conventions

> **执行前**：用 Skill 工具加载 `python-script-conventions`，并 skill_view 其 `templates/agent_script_template.py`。
> **每个脚本流程固定**（见下方 Task C-template），按依赖顺序从数据流上游到下游。
> **验证策略**：由于无测试框架也无 DIT_HAP 自产基线，采用双重验证——
>   (a) **静态**：`python -c "import ast; ast.parse(open('FILE').read())"` 语法通过 + `ruff check FILE`（若可用）+ 人工核对 checklist；
>   (b) **行为**：若该脚本能在合理时间内用小样本或已有中间产物跑通，则重构前后各跑一次，`diff` 输出（数值型用 `python` 读入对比 `.equals()` 或数值容差）。无法跑通的（依赖全流程上游产物）则依赖静态验证 + 逐行 diff 审查逻辑段未变。

### Task C-template（对每个脚本套用）

**Step 1:** 读整脚本，记录其**逻辑段落**（哪些函数做什么），标出 I/O 边界与算法核心。

**Step 2:** 判定 standalone script（有 main/CLI）还是 library module（被 import）。上游脚本几乎都是 standalone。

**Step 3:** 按 7 段布局重排：IMPORTS(三组排序) → DECORATORS → CONSTANTS/StrEnum → CONFIG(frozen dataclass) → LOGGING(setup + --verbose) → CORE(@logger.catch, guard clauses, match/case) → MAIN(parse_args + sys.exit(main()))。

**Step 4:** 现代化改写（**不动算法**）：
- `Optional[X]`→`X | None`；`typing.List/Dict`→`list/dict`
- 配置类→`@dataclass(kw_only=True, slots=True, frozen=True)`
- 成组字符串标记→`StrEnum`
- 长 if-elif→`match/case`；边界→guard clauses
- `print`→`logger`；`os.path`→`Path`；`.format()/%`→f-string
- 模块 docstring 按 §3.2（title/description/Input/Output/Usage/metadata，英文）
- 函数 docstring 单行；删无用 import

**Step 5:** 静态验证：
```bash
python -c "import ast,sys; ast.parse(open(sys.argv[1]).read())" <FILE>
ruff check <FILE> 2>/dev/null || echo "ruff 不可用，跳过"
```

**Step 6:** 行为验证（能跑则跑，见上方策略）；不能跑则逐段 diff 审查算法未变。

**Step 7:** Commit 单脚本：
```bash
git add <FILE>
git commit -m "refactor(<module>): modernize <script>.py to python-script-conventions"
```

### 重构顺序（按数据流 + 先易后难）

**C1–C7 read_processing:** parse_bam_to_tsv → filter_aligned_reads → extract_insertion_sites → merge_strand_insertions → concatenate_timepoint_data → annotate_genomic_features → (merge_similar_timepoints, concat_counts_and_annotations, reads_hard_filtering 已在 B 阶段接线，此处补齐重构；reads_hard_filtering 缺 section 分区优先)

**C8 reference_data:** extract_genome_region

**C9–C15 depletion_scoring:** def_ctr_insertions → impute_missing_values_using_FR → insertion_level_depletion_analysis_no_replicates → insertion_level_depletion_analysis_has_replicates → curve_fitting → compute_r2_weights → gene_level_depletion_analysis(缺 section 分区优先)

**C16–C22 quality_control:** extract_mapping_filtering_statistics → PBL_PBR_correlation_analysis → read_count_distribution_analysis → insertion_orientation_analysis → gene_coverage_analysis(有 5 处 print 优先) → insertion_density_analysis(908 行，最大) → distribution_of_curve_fitting_results

**Step 末:** 全部重构完成后再跑一次全局 dry-run + lint 确认规则与脚本 CLI 仍匹配（重构可能改了 argparse）：
```bash
snakemake -n --configfile config/config_HD_generationRAW.yaml 2>&1 | tail
snakemake -n --configfile config/config_LD_generationRAW.yaml 2>&1 | tail
```

---

## 阶段 D：文档 / 配置清理

### Task D1: 修 config workdir

**Files:** `config/*.yaml`（8+ 个）

**Step 1:** 扫描当前指向：
```bash
grep -rn "workdir:" config/
```

**Step 2:** 全部 `workdir: /data/c/yangyusheng_optimized/DIT_HAP_pipeline` → `.../DIT_HAP`。

**Step 3:** 验证：
```bash
grep -rn "DIT_HAP_pipeline" config/   # 期望: 无输出
```

**Step 4:** Commit：
```bash
git add config/
git commit -m "fix: point config workdir to DIT_HAP"
```

---

### Task D2: 新建 README.md（英文，仅上游）

**Files:** Create `README.md`

**Step 1:** 写 README：概述（上游到 gene-level DR/DL 表格 + QC）、安装、快速开始、四模块架构（reference_data/read_processing/depletion_scoring/quality_control）、output structure（对齐真实 results 编号）、"下游分析见 DIT_HAP_pipeline"指引。英文。

**Step 2:** Commit：
```bash
git add README.md && git commit -m "docs: add upstream-only README"
```

---

### Task D3: 新建 CLAUDE.md

**Files:** Create `CLAUDE.md`

**Step 1:** 写 CLAUDE.md：项目概述、常用命令、Snakemake 配置说明、四模块架构（新命名）、output structure、Python 脚本约定**指向 `python-script-conventions` skill**（不再内嵌旧 STYLE_GUIDE）。

**Step 2:** Commit：
```bash
git add CLAUDE.md && git commit -m "docs: add CLAUDE.md pointing to python-script-conventions skill"
```

---

### Task D4: 最终验证

**Step 1:** 全模块 dry-run 两分支：
```bash
snakemake -n --configfile config/config_HD_generationRAW.yaml 2>&1 | tail
snakemake -n --configfile config/config_LD_generationRAW.yaml 2>&1 | tail
```

**Step 2:** `snakemake --lint`。

**Step 3:** 全仓 grep 确认无旧命名残留：
```bash
grep -rn "depletion_analysis\|preprocessing\|preparation" workflow/rules/ | grep -v "message\|#"
grep -rn "run:" workflow/rules/
grep -rn "scripts/reference/" workflow/
```
期望：均无输出（或仅 message 文案）。

**Step 4:** Commit 收尾（如有）。

---

## 验证清单（贯穿全程）

- [ ] 每改规则跑 `snakemake -n` 不报错
- [ ] 脚本重构后 `ast.parse` 通过 + checklist 核对
- [ ] 无 `run:` 块残留
- [ ] 无旧命名（preprocessing/depletion_analysis/preparation 作为路径）残留
- [ ] config workdir 指向 DIT_HAP
- [ ] 每段独立 commit
