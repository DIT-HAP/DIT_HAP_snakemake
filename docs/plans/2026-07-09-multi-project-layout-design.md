# 多 Project 布局重构设计

**日期**: 2026-07-09
**范围**: 仅 `/data/c/yangyusheng_optimized/DIT_HAP`
**状态**: 已确认，待实施
**前置**: 上游重构 + config/sample-sheet schema 校验已完成（分支 refactor/upstream-cleanup）

---

## 1. 目标

把当前"config 集中在根 config/、输出散在 results/reports/logs/{project_name}/"的布局，改为**每个 project 自包含**的 `projects/{project_name}/` 布局，同时保持 reference 数据、代码、全局配置的中心共享。既支持日常共享代码/reference 跑多 project，也支持完成的 project 整体打包归档。

## 2. 已确认的决策

- **布局**：`projects/{project_name}/` 收纳该 project 的 config + sample_sheet + results + reports + logs。
- **身份键**：`project_name`（如 `HD_DIT_HAP_generationRAW`），目录名 = project_name。
- **共享资源**：完全自包含 —— 被多 config 共用的 sample sheet 复制多份到各 project。
- **中心共享（不进 project）**：`resources/pombase_data/{release_version}/`（reference，按版本共享）、`workflow/`（代码/schema）、`config/DIT_HAP.mplstyle` + `config/multiqc_config.yml`（全局配置）。
- **project 内 config 位置**：放 `projects/{project_name}/config/` 子目录（与 results 等平级）。
- **模板**：根 `config/` 下加 `config.template.yaml` + `sample_sheet.template.tsv` 作为新建起点。

## 3. 目标布局

```
DIT_HAP/
├── Snakefile
├── config/
│   ├── DIT_HAP.mplstyle              # 全局共享
│   ├── multiqc_config.yml            # 全局共享
│   ├── config.template.yaml          # 新建 project 的 config 起点（全字段+注释）
│   └── sample_sheet.template.tsv     # 新建 sample sheet 起点（表头+示例）
├── projects/
│   ├── HD_DIT_HAP_generationRAW/
│   │   ├── config/
│   │   │   ├── config.yaml
│   │   │   └── sample_sheet.tsv
│   │   ├── results/
│   │   ├── reports/
│   │   └── logs/
│   └── ... (9 个 project)
├── resources/pombase_data/{release_version}/   # 中心共享（按版本）
└── workflow/  (schemas / rules / scripts / envs 中心共享)
```

## 4. config → project 映射（9 个）

| project_name (目录) | 原 config | sample sheet（复制进 project/config/） |
|---|---|---|
| HD_DIT_HAP_generationRAW | config_HD_generationRAW.yaml | sample_sheet_HD.tsv |
| HD_DIT_HAP_generationPLUS1 | config_HD_generationPLUS1.yaml | sample_sheet_HD.tsv (共享,复制) |
| HD_DIT_HAP | config_HD.yaml | sample_sheet_HD.tsv (共享,复制) |
| LD_DIT_HAP_generationRAW | config_LD_generationRAW.yaml | sample_sheet_LD.tsv |
| LD_DIT_HAP_generationPLUS1 | config_LD_generationPLUS1.yaml | sample_sheet_LD.tsv (共享,复制) |
| HD_diploid | config_HD_diploid.yaml | sample_sheet_diploid.tsv |
| LD_haploid | config_LD_haploid.yaml | sample_sheet_LD_haploid.tsv |
| Spikein | config_spikein.yaml | sample_sheet_spikein.tsv |
| Spore2YES6_1328 | config_1328_spore2YES6.yaml | sample_sheet_LD1328_spore2YES6.tsv |

## 5. Snakefile 选择机制

```python
# 选择当前 project（改这一行切换实验）
project = "HD_DIT_HAP_generationRAW"
# project = "HD_DIT_HAP_generationPLUS1"
# ... 其余 project 注释列出

configfile: f"projects/{project}/config/config.yaml"
validate(config, "workflow/schemas/config.schema.yaml")

# 一致性校验：目录名 == config 里的 project_name
assert config["project_name"] == project, (
    f"project dir '{project}' != config project_name '{config['project_name']}'"
)

sample_sheet = pd.read_csv(config["sample_sheet"], sep="\t", dtype=str)
validate(sample_sheet, "workflow/schemas/samples.schema.yaml")
```

- `project` 变量是唯一真相，用于拼路径。
- 断言 `project_name == project`，杜绝目录名与 config 字段错配。

## 6. 路径改写（改动量最大）

4 个 `.smk` 文件的输出/日志前缀：
- `results/{project_name}/...` → `projects/{project_name}/results/...`
- `reports/{project_name}/...` → `projects/{project_name}/reports/...`
- `logs/{project_name}/...` → `projects/{project_name}/logs/...`
- `resources/pombase_data/...` **不变**（中心共享）

规则里是 f-string 与 `{project_name}` wildcard 混用，需逐文件核对替换，**不误伤 resources 路径**。

每个 project 的 `config.yaml`：
- `sample_sheet:` → `projects/{project_name}/config/sample_sheet.tsv`
- `workdir` 保持指向仓库根（DIT_HAP）

## 7. 迁移步骤

1. 建根 `config/config.template.yaml` + `config/sample_sheet.template.tsv`。
2. 对 9 个 project：建 `projects/{project_name}/config/`，`git mv` config 进去改名 `config.yaml`；sample sheet 复制进去改名 `sample_sheet.tsv`（共享的复制多份）。
3. 改各 config.yaml 的 `sample_sheet:` 路径。
4. 改 4 个 .smk 的路径前缀。
5. 改 Snakefile（选 project + 校验一致）。

## 8. 验证

- 9 个 project 各跑 `snakemake -n`：schema 通过、DAG 正常、输出落 `projects/{name}/`、resources 指中心。
- 负向：目录名 ≠ config project_name → 断言拦下。
- 归档：`tar tf projects/{name}/...` 确认自包含（含 config+sheet+输出，不含 reference）。

## 9. 提交拆分

- commit 1：模板 + config/sample_sheet 搬迁 + config 内 sample_sheet 路径改写。
- commit 2：规则路径前缀改写 + Snakefile 选择机制。
