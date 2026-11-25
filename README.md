# CMIP Future Projection

用于 CMIP6/CMIP5 未来情景分析的可复用 skill 与 Python 工具。通过一份 YAML 配置，完成原始数据质量检查、季节或区域指数、可选偏差订正、多模式汇总、结果图表，以及明确选择的日尺度极值分析。

*A reusable skill and configuration-driven Python workflow for CMIP projections, regional indices, explicit bias adjustment, ensemble summaries, and optional daily extremes.*

默认不进行偏差订正，也不自动将年度均值送入 GEV 拟合。先明确输入与研究问题，再选择需要的计算分支。这里的“通用”指可配置、可检查、可扩展；不同变量和事件仍需要相应的科学定义。

## 内容

| 位置 | 用途 |
|---|---|
| [cmip-future-projection/SKILL.md](cmip-future-projection/SKILL.md) | AI 助手的任务入口、工作流程与能力边界 |
| [scripts/run.py](cmip-future-projection/scripts/run.py) | 命令行入口 |
| [scripts/cmipkit/](cmip-future-projection/scripts/cmipkit/) | 数据读取、指数、订正、集合、极值、绘图与流程代码 |
| [templates/](cmip-future-projection/templates/) | 五类可修改的分析配置 |
| [references/](cmip-future-projection/references/) | 方法、配置、变量配方、绘图及排错说明 |
| [tests/](cmip-future-projection/tests/) | 合成数据和数值回归检查 |
| [agents/openai.yaml](cmip-future-projection/agents/openai.yaml) | 支持该格式的助手界面元数据 |

整个 `cmip-future-projection/` 目录构成 skill；命令行工具也可以独立使用。接入支持本地 skill 的助手时，按该助手的配置方式加载这个完整目录，并保留相对资源路径。不要只复制 SKILL.md 而丢掉脚本与模板。

## 安装依赖

建议使用独立的 Python 3.11 或 3.12 环境。这两个版本也是仓库 CI 的目标环境。

~~~bash
python -m venv .venv
~~~

Linux/macOS 激活：

~~~bash
source .venv/bin/activate
~~~

Windows PowerShell 激活：

~~~powershell
.venv\Scripts\Activate.ps1
~~~

然后在仓库根目录执行：

~~~bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python cmip-future-projection/scripts/run.py doctor
~~~

核心依赖支持 NetCDF、非标准日历、CSV/NPZ、统计和普通绘图。仅在需要时安装可选依赖：

~~~bash
python -m pip install -r requirements-optional.txt
~~~

Cartopy 用于地图背景；GRIB 需要 cfgrib 和可用的 ecCodes；Dask 供适用的数据处理扩展使用。安装 Dask 不会自动把当前全部计算变成惰性或分布式计算。Cartopy 的海岸线数据也可能需要事先准备；离线环境应使用已缓存的数据或清楚标注的普通经纬度图。

requirements 文件给出兼容范围，并非精确锁定的环境。正式研究应另行记录实际安装版本、数据版本和运行配置。

## 快速开始

在仓库根目录创建一个独立分析目录：

~~~bash
python cmip-future-projection/scripts/run.py init myrun --template grid
~~~

修改生成的 `myrun/config.yaml`，至少核实路径、变量名、单位、频率、区域和时间段。模板中的示例数据并不随仓库提供；模板通过结构校验不代表这些文件存在。

~~~bash
python cmip-future-projection/scripts/run.py validate myrun/config.yaml
python cmip-future-projection/scripts/run.py validate myrun/config.yaml --check-inputs
python cmip-future-projection/scripts/run.py run myrun/config.yaml --dry-run
python cmip-future-projection/scripts/run.py run myrun/config.yaml
~~~

`validate` 先检查配置；`--check-inputs` 检查所有选中文件的元数据、时间轴和请求时段覆盖。`run --dry-run` 执行配置要求的质量预检，不计算正式预估。报告会区分元数据检查与完整数值扫描，列出缺段及其他警告。

**相对数据和输出路径以配置文件所在目录为基准。** 例如，`myrun/config.yaml` 内的 `data/models` 指向 `myrun/data/models`，与启动命令时的工作目录无关。路径支持 `~` 和环境变量，文件模式支持 `{model}`、`{experiment}`、`{member}` 等占位符。

## 选择模板

| `--template` | 适合的任务 | 关键定义 |
|---|---|---|
| `grid` | 月尺度区域温度、情景时间序列与空间变化 | 日历月长加权的季节平均；模板显式采用均值订正 |
| `series` | 已计算好的年度指数 CSV/NPZ | 保留上游指数定义，不重建格点；默认无订正和 GEV |
| `precip` | 月尺度降水与 QDM | 输出季节平均率 mm/day；不是把月平均率相加后当总量 |
| `regional-extremes` | 区域日序列的季节或年度块最大值 | 先每天做空间平均，再取季节最大值并拟合 GEV |
| `daily` | 日尺度 POT/GPD 极值 | 固定历史参考阈值、按真实日期去簇、按有效观测时长估计事件率 |

这些是起点，不是必须照用的区域、模式、情景和校正方法。默认配置的 `correction.method` 为 `none`，模板则可能为对应示例明确选择 `mean` 或 QDM。修改方法前，请阅读 [方法说明](cmip-future-projection/references/methods.md)。

已有结果只需要重绘时，可以使用明确的绘图适配脚本。CLI 的 `--stages` 会补齐所需上游计算，但不是从所有已保存结果恢复任意计算阶段的完整断点续跑系统。

## 原始 CMIP 数据质量检查（v1.2.0）

首次接入一批原始数据，先运行分块数值检查并保存报告：

~~~bash
python cmip-future-projection/scripts/run.py validate myrun/config.yaml --check-data --report-json myrun/input-quality.json
~~~

该命令检查配置选中的时间、层次与区域，不启动预估。报告列出问题文件、时间缺段、重叠记录、单位/日历/成员/网格冲突，以及 NaN、无穷值、全空时次和数值范围。普通 `--check-inputs` 只检查元数据和坐标，不声称已扫描全部数值。

正式运行可在 YAML 中设置：

~~~yaml
quality:
  scan: full
  coverage: error
  require_units: true
  max_missing_fraction: null
  valid_min: null
  valid_max: null
  chunk_size: 32
  max_chunk_values: 2000000
runtime:
  max_input_values: 30000000
~~~

默认 `scan: headers`、`coverage: warn`、`require_units: false`，以兼容已有资料；有明确的缺测率或值域阈值时，正式运行自动启用完整数值扫描。值域按配置转换后、时间聚合前的单位解释。缺测比例分母包含所选时空中的静态掩膜和重复文件记录，不能把适用于温度的阈值直接套给海洋变量。

读取器先裁剪再加载，并逐季节年处理年度分析；日尺度读取在加载前检查累计数组规模。`max_input_values` 是数组值数量保护，不是精确内存上限。损坏或失配的缓存会重建，NetCDF、NPZ 和运行状态采用临时文件完成后原子替换，避免把未写完的文件当作成品。

正常分析写入 `quality_report.json`；`manifest.json` 和 `report.md` 记录检查级别。若预检失败，需要留存报告时使用 `--report-json PATH`，命令以非零状态退出。原始数据不会被修改，缺测和异常不会被自动填补或裁剪。完整操作与边界见[原始数据质量检查](cmip-future-projection/references/quality-control.md)。

## 输入与能力范围

| 输入或需求 | 当前支持方式 |
|---|---|
| 一维经纬度坐标的规则格点 NetCDF | 原生格点流程；按配置选变量、层次、区域与季节 |
| 兼容的规则格点 GRIB | 使用 cfgrib/ecCodes；先确认消息分组、单位和时间语义 |
| 年度 CSV 或 NPZ | 原生 series 流程；CSV 每年一行，NPZ 含 years/values |
| 每日格点区域极值 | daily 分支；只支持单成员策略，并保持历史—情景成员一致 |
| WorldClim GeoTIFF、曲线或旋转网格、非结构网格 | 需要经验证的输入适配或预处理 |
| 每日站点 CSV、热浪持续时间、复合事件 | 需要额外的数据适配或事件算法 |
| ESGF 下载、守恒重网格、统计降尺度、非平稳极值 | 不作为当前 CLI 已实现的自动功能 |

十二个月的多年气候态不能被当作逐年或逐日样本。插值到更细网格也不会自动获得更细尺度的气候信息。日尺度 anomaly 尚未实现，配置会明确拒绝该选项；series 的异常、标准化和去趋势应在上游采用固定历史参数处理。

## 结果与复现

结果写入配置的 `project.outdir`，按启用的阶段产生：

- `config.resolved.yaml`：解析后的实际配置。
- `manifest.json`：输入、模型/成员、阶段状态、失败记录和产物信息。
- `quality_report.json`：本次输入的检查范围、逐文件发现和质量摘要。
- `report.md` 与 `tables/`：结果说明和数值表格。
- `corrected/`、`params/`、`series/`、`ensemble/`、`extremes/`：对应计算产物。
- `figures/` 和 `logs/`：图件和日志。

manifest 区分预检发现的 `preflight_inputs` 与实际消费的 `actual_inputs`，也区分本次生成/更新的产物和本次没有生成的旧文件。不要仅凭目录中某张图片存在，就声称本次运行生成了它。

默认启用严格配置和处理检查，避免未知参数或模型错误被静默忽略。更改科学配置后，建议使用新的输出目录；若明确使用 `--overwrite`，程序会把原运行目录移到带 `.previous-<时间戳>` 后缀的位置，再建立新的运行目录。

可用 `periods.validation: [起年, 止年]` 指定年度历史诊断的评估窗口。它筛选参考与模式共同存在的年度样本，不自动下载额外数据，也不自动构成独立验证；若与校准窗口重叠，会给出提示。正式比较仍需判断样本独立性与科学问题是否匹配。

## 方法上的重要约定

- 月平均量按各自日历的月份长度加权；缺失时间与不完整季节受到显式检查。`sum` 不会自动把率积分为总量。
- 均值与方差订正会改变变化幅度；QDM 保留的是设计中的分位数变化，并不保证时间趋势完全不变。当前 QDM 使用整个配置未来窗口估计未来 CDF。
- 不自动将季节均温、季节总量或任意年度指数解释为块最大值。年度 GEV 需要显式启用并说明样本类型。
- SciPy 的 GEV 形状参数采用 `c = −ξ`；GPD 则采用 `c = ξ`。同样本拟合后的普通 KS p 值仅为未校准诊断，不作为“全部拟合必须通过”的验收门槛。
- 模型间分位区间表示模型离散度，不是置信区间或校准后的概率。成员平均会抑制变率，不能作为极值分析的默认处理。
- 所有当前极值拟合均为平稳模型。百年变化中的固定窗口拟合需要明确解释，不自动代表每一个未来年份的风险。

完整公式、边界条件和来源见 [方法说明](cmip-future-projection/references/methods.md) 与 [方法来源](cmip-future-projection/references/sources.md)。

## 本地检查与 CI

在仓库根目录执行：

~~~bash
python cmip-future-projection/tests/test_methods.py
python cmip-future-projection/tests/test_io_contracts.py
python cmip-future-projection/tests/test_raw_io.py
python cmip-future-projection/tests/test_quality.py
python cmip-future-projection/tests/test_quality_workflow.py
python cmip-future-projection/tests/test_statistical_contracts.py
python cmip-future-projection/tests/test_workflow.py
python cmip-future-projection/tests/test_daily.py
python cmip-future-projection/scripts/run.py selftest
~~~

GitHub Actions 配置对 Python 3.11 和 3.12 运行同样的检查，只使用合成数据，不需要用户服务器目录或研究数据。首次远端运行成功前，不将 CI 配置本身描述为“已通过 CI”。

这些检查覆盖已知的数据处理、统计计算和流程回归问题，**不替代真实数据验证、方法适用性判断、独立验证期评估或论文结论审查**。本仓库不附带无法复现其输入来源的旧项目结果图。

## 进一步阅读

- [操作流程](cmip-future-projection/references/workflow.md)
- [配置参考](cmip-future-projection/references/config-reference.md)
- [变量与数据配方](cmip-future-projection/references/variables.md)
- [绘图说明](cmip-future-projection/references/plotting.md)
- [故障排查](cmip-future-projection/references/troubleshooting.md)
