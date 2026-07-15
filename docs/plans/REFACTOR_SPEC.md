# 单脚本重构规范（阶段 C 通用指令）

> 每个 implementer 重构 **一个** 脚本时遵循本规范。目标：把脚本改写到 `python-script-conventions` skill 的 Modern Python 3.12+ 标准，**只改结构/风格，绝不改算法或行为**。

## STEP 0：加载标准
1. 用 Skill 工具加载 `python-script-conventions`（打印完整标准）。
2. 读模板：`/data/a/yangyusheng/.claude/skills/python-script-conventions/templates/agent_script_template.py`

## 绝对铁律
1. **不改算法 / 行为。** 重构后对相同输入必须产生相同输出：相同的 pandas 操作、相同的数值逻辑、相同的 I/O 格式、**完全相同的 CLI flag（名称/类型/required/nargs 语义）**、相同日志信息的含义。
2. **保留这些易被"顺手改坏"的行为**（若脚本里有）：
   - `df.groupby(..., axis=1)` / 任何 pandas `axis=1` 弃用写法——**保持原样**，不要"修复"。
   - 多索引读写：`index_col=[...]`, `header=[...]`。
   - 分块处理 `chunksize`、逐行/逐块循环逻辑。
   - pysam / pybedtools / subprocess 调用的参数与顺序。
   - 输出写法：`to_csv` 的 sep/header/index 参数。
3. 只动目标文件一个。**不碰** `/data/c/yangyusheng_optimized/DIT_HAP_pipeline`。
4. 不引入 `snakemake.script` 的 `snakemake` 对象——这些是命令行脚本，保持 argparse。

## 要改的（对齐 skill）
- **7 段布局**：IMPORTS → DECORATORS(无则省) → GLOBAL CONSTANTS & ENUMS → CONFIGURATION & DATACLASSES → LOGGING SETUP → CORE LOGIC → MAIN EXECUTION。用模板的 `# ===` banner。
- **shebang + PEP723**：加 `#!/usr/bin/env python3`；PEP 723 依赖块注释掉（照模板，列出该脚本真实第三方依赖如 pandas/loguru/pysam/pybedtools/numpy/scipy 等）。
- **模块 docstring**：§3.2 格式（Title === / 描述 / Input / Output / Usage / Author-Date-Version），**英文**，描述该脚本真实功能。Author 用 `Yusheng Yang (guidance) + Claude (implementation)`，Date `2026-07-09`，Version `1.0.0`。
- **imports**：三组（stdlib / 数据处理如 pandas,numpy / 第三方如 loguru,pysam），组内字母序，全部置顶，**删未用**。
- **类型**：`Optional[X]`→`X | None`；`Tuple/List/Dict[...]`→`tuple/list/dict[...]`；相应清理 `from typing import ...`。
- **配置类**：用 `@dataclass(kw_only=True, slots=True, frozen=True)`。若原 `__post_init__` 里给 self 赋值（frozen 禁止），把该规范化移到构造前（main 里）或用 `object.__setattr__`；`mkdir` 之类不赋值属性的副作用可留在 `__post_init__`。用关键字实例化。
- **日志**：`setup_logger(log_level)`（loguru，禁 print），保留 `--verbose`→DEBUG。禁 `print()`——所有输出走 logger。
- **错误处理**：核心函数用 `@logger.catch`。删掉"只是 log 后再 raise"的冗余 try/except（decorator 已带 traceback）。**保留**用于特定控制流的 try/except（如 KeyError 回退、逐块跳过、特定异常不同处理）——这类**不要删**。保留真正的 `raise ValueError(...)` 校验。
- **控制流**：长 if-elif 链且分支是固定类别时用 `match/case`；边界用 guard clause 提前返回。**不要为改而改**——若原 if-elif 不是明显的类别分发，保持原样。
- **StrEnum**：仅当脚本里有一组**成组、重复使用**的字符串标记（列名/模式/状态）才引入；单例常量用普通 UPPERCASE 常量。**不要硬造 StrEnum**（YAGNI）。
- **docstring**：函数/类单行。
- **MAIN**：`def parse_args()` + `def main() -> int:`，成功 return 0 / 预期失败 return 1，`if __name__ == "__main__": sys.exit(main())`。按模板：main() 内用 try/except 返回 1，**不要**给 main 加 `@logger.catch`（否则吞掉退出码）。

## 验证
1. `python -c "import ast; ast.parse(open('<FILE>').read())"` → 无错。
2. `ruff check <FILE>`（若可用）→ 报告输出；与 skill 冲突的 nit 可忽略但要报告。
3. 逐行自检算法未变（尤其 pandas 操作、循环、I/O）。
4. 核对 parse_args 的 flag 与调用它的 Snakemake 规则完全一致（prompt 里会给出该规则的调用行）。

## 提交
```
git add <FILE>
git commit -m "refactor(<module>): modernize <script>.py to python-script-conventions

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

## 汇报
- 验证输出（ast、ruff）。
- 一句话确认 "Algorithm unchanged" + 你如何保留关键行为。
- CLI flag 未变的确认。
- commit SHA。
- 任何判断取舍；若 skill 与"保留行为"冲突，**保留行为**并说明。
