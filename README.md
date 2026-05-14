# 长江渔业生态模型 README

本项目是 `2026A-长江流域人类活动对长江渔业的影响` 建模文档的 Python 实现，覆盖 Q1-Q5 的常微分方程仿真、稳定性分析、最大 Lyapunov 指数判别、EI 综合生态指数计算，以及关键缺失参数的局部灵敏度分析。

## 1. 项目结构

```text
config.py
models_ode.py
stability_analysis.py
plot_utils.py
main.py
README.md
附件1 长江鱼类洄游路径.pdf
附件2 长江鱼类的食物链.pdf
附件3 2022年长江流域水生生物资源及生境状况公报.pdf
outputs/
```

各模块职责如下：

- `config.py`：集中管理初始值、参数默认值、参数来源、区间赋值、情景配置和输出目录。
- `models_ode.py`：实现 Q1-Q5 的 ODE 方程与公共辅助函数。
- `stability_analysis.py`：实现平衡点求解、Jacobian、特征值分析、最大 Lyapunov 指数和参数扫描。
- `plot_utils.py`：统一输出时间序列图、相图、分岔图、EI 对比图和灵敏度曲线。
- `main.py`：按 Q1-Q5 顺序执行完整流程，并补充参数表导出与灵敏度分析。

## 2. 环境与运行方式

建议在 `LLM_Classification` conda 环境中运行。

```powershell
conda activate LLM_Classification
python main.py
```

如果不切换环境，也可以直接指定解释器：

```powershell
D:\anaconda\envs\LLM_Classification\python.exe main.py
```

运行完成后会自动生成：

- `outputs/figures/`：所有图表
- `outputs/data/`：时序数据、EI 结果、灵敏度分析结果
- `outputs/reports/`：摘要报告、参数来源表、参数区间表、参考来源索引表

## 3. 程序主要运行思路

整体流程是：

`公报/附件/文献 -> config.py 参数组织 -> models_ode.py 数值仿真 -> stability_analysis.py 指标分析 -> plot_utils.py 图表输出`

各题执行逻辑如下：

1. Q1：分别模拟 `P-Z-L-B` 和 `V-R-C-M` 两个子系统，对比禁渔前后 10 年演化。
2. Q2：用 `p_list` 表示坝群阻断概率，并按 `H = 1 - prod(1 - p_i)` 聚合阻断压力，比较低/中/高阻断和有无放流的鲟类恢复曲线，再单独模拟江豚。
3. Q3：构建功能群层面的广义 Lotka-Volterra 食物网，先仿真，再求平衡点和 Jacobian 特征值。
4. Q4：构造三维资源-中层-顶层系统，生成稳定极限环和混沌吸引子示例，并用最大 Lyapunov 指数判别。
5. Q5：在污染和入侵压力下模拟本土种与外来种竞争，计算 EI 及其各分量。
6. 灵敏度分析：对区间赋值的关键参数做单因素扫描，输出参数-响应曲线和结果表。

## 4. 参数赋值原则

当前参数体系分为三层：

1. 公报/附件驱动的初始值和基线量级
2. 文献支撑的生态学参数
3. 缺失参数的区间赋值与灵敏度分析

实现时采用以下约定：

- 模型状态量统一解释为“相对生物量/相对数量/标准化丰度”，不直接把公报统计量原样抄成状态变量。
- 当公报只能提供趋势或物种恢复信号时，用于确定量级和相对大小，不伪造高精度。
- `H` 不作为独立输入，而由 `p_list` 统一推导。
- Q4 的混沌参数是展示性参数，用于说明动力学现象，不解释为实测标定值。

## 5. 初始值表

下表只列最关键的状态量，完整版本会在运行后导出到 `outputs/reports/initial_conditions_table.csv`。

| 变量 | 所属模型 | 默认值 | 来源 | 赋值逻辑 |
|---|---|---:|---|---|
| `P` | Q1 浮游植物 | 82.0 | 2022 公报 | 根据 2020-2022 水生生物恢复信号映射为较高初级生产者量级 |
| `Z` | Q1 浮游动物 | 30.0 | 附件2 食物链 | 位于浮游植物与滤食鱼之间，设为中层资源量级 |
| `L` | Q1 鲢鱼 | 19.3 | 2022 公报 | 结合鱼类恢复趋势映射到中等恢复量级 |
| `B` | Q1 鳙鱼 | 15.2 | 附件2 食物链 | 作为较高营养级滤食鱼，初值略低于鲢鱼 |
| `C` | Q1 草鱼 | 18.0 | 附件2 + 文献 | 结合四大家鱼栖息地受坝体影响的研究做中等量级设定 |
| `S` | Q2 鲟类 | 1.2 | 2022 公报 | 公报表明天然繁殖仍弱，因此设置为极低基线 |
| `D` | Q2 江豚 | 12.49 | 官方通报 | 依据 2022 年江豚数量 1249 头，按 100:1 缩放 |
| `R,F,T,I` | Q3 功能群 | 120, 46, 12, 9 | 公报 + 附件2 | 分别表示资源层、中层鱼类、顶层捕食者、入侵群体 |
| `N,I` | Q5 本土/入侵 | 63.8, 8.5 | 公报 + 附件2 | 反映“本土占优、入侵已存在但未主导”的背景 |

## 6. 文献参数表

运行后完整表格保存在 `outputs/reports/parameter_sources_table.csv`。下面列出最关键的一批。

| 参数 | 含义 | 默认值 | 来源类型 | 参考依据 |
|---|---|---:|---|---|
| `r_P, r_Z, r_L, r_B` | Q1 各群体增长率 | 0.36-1.10 | 文献 | [Tian et al., 2006](https://doi.org/10.1016/j.ecolmodel.2005.09.003) 的生态模型参数化量级 |
| `alpha_PZ, alpha_PL, alpha_ZB` | Q1 摄食强度 | 0.017-0.020 | 文献 | [Holling, 1959](https://doi.org/10.4039/Ent91293-5) 支撑捕食项和摄食函数解释 |
| `mu_C, mu_M, mu_L, mu_B` | Q1 鱼类基础死亡率 | 0.10-0.12 | 文献 | [Tian et al., 2006](https://doi.org/10.1016/j.ecolmodel.2005.09.003) 的生物量模型参数量级 |
| `rho` | Q2 阻断压力放大系数 | 0.90 | 区间+文献背景 | 受洄游受阻和繁殖损失共同影响，量级参考坝体生态影响研究 |
| `p_list` | 坝群阻断概率 | `[0.18,0.24,0.16]` | 附件1 | 由洄游路径和坝体阻隔结构转换为概率输入 |
| `A[1,0], A[1,2]` | Q3 关键 GLV 相互作用 | 0.010, -0.010 | 文献+区间 | 附件2 给方向，数值由 GLV 量级估计并进入灵敏度分析 |
| `eps, omega` | Q4 周期扰动幅值与频率 | 0.05/0.28, 2.0/5.6 | 展示性参数 | [May, 1976](https://doi.org/10.1038/261459a0) 支撑复杂动力学展示 |
| `theta_N, theta_I` | Q5 污染对容纳量压缩 | 0.95, 0.38 | 文献+区间 | [Rubach et al., 2006](https://pubmed.ncbi.nlm.nih.gov/16193755/) 支撑污染降低增长潜力/容纳量 |
| `gamma_N, gamma_I` | Q5 污染诱导额外死亡率 | 0.13, 0.05 | 文献+区间 | 污染应激增加死亡率的经典生态毒理方向 |
| `alpha_NI, beta_IN` | Q5 本土-入侵竞争强度 | 0.72, 0.58 | 文献+区间 | [Palmer et al., 2017](https://doi.org/10.1016/j.jenvman.2017.03.060) 与入侵鱼类综述支撑 |

## 7. 区间参数表与灵敏度分析

下列参数缺少长江专属精确实测值，因此采用区间赋值，并强制纳入单因素灵敏度分析。完整区间表见 `outputs/reports/parameter_ranges_table.csv`。

| 参数 | 默认值 | 区间 | 原因 | 响应量 |
|---|---:|---:|---|---|
| `rho` | 0.90 | 0.60-1.20 | 阻断对繁殖成功率损失的放大作用存在明显不确定性 | 鲟类 10 年末恢复量 |
| `stocking_width` | 0.12 | 0.06-0.20 | 放流脉冲持续时间受数值近似和管理方式影响 | 鲟类 10 年末恢复量 |
| `stocking_scale` | 1.00 | 0.60-1.40 | 放流规模年度波动较大 | 鲟类 10 年末恢复量 |
| `p_list_scale` | 1.00 | 0.70-1.30 | 通道可通过性和坝群联合作用存在不确定性 | 鲟类 10 年末恢复量 |
| `A_10` | 0.010 | 0.006-0.016 | 资源层对中层鱼类补给强度缺少专属标定 | 最大实部特征值 |
| `A_12` | -0.010 | -0.016 至 -0.006 | 顶层捕食对中层压制强度不确定 | 最大实部特征值 |
| `theta_N` | 0.95 | 0.70-1.20 | 本土鱼类污染敏感性区域差异大 | `EI`、`N_final`、`I_final` |
| `theta_I` | 0.38 | 0.20-0.60 | 入侵种通常更耐扰动，但缺乏统一定量值 | `EI`、`N_final`、`I_final` |
| `gamma_N` | 0.13 | 0.08-0.18 | 污染诱导额外死亡率受污染类型影响 | `EI`、`N_final`、`I_final` |
| `gamma_I` | 0.05 | 0.02-0.10 | 入侵种额外死亡率通常低于本土种 | `EI`、`N_final`、`I_final` |
| `alpha_NI` | 0.72 | 0.50-0.95 | 入侵种对本土种竞争压力难直接观测 | `EI`、`N_final`、`I_final` |
| `beta_IN` | 0.58 | 0.40-0.80 | 本土种对入侵种反制能力不确定 | `EI`、`N_final`、`I_final` |

灵敏度分析的实现方式固定为：

- 只做关键参数的单因素局部分析
- 每次只改变一个参数
- 其余参数固定在默认值
- 输出 `outputs/data/sensitivity_analysis.csv`
- 同时输出对应的参数-响应曲线图

## 8. H 参数的统一实现

`H` 表示坝群阻断导致的综合破碎化压力。程序中不允许在情景层手写一个常数 `H`，统一采用：

```python
H = 1 - prod(1 - p_i)
```

这里的 `p_i` 来自附件1 对洄游路径和坝体阻隔结构的描述，经建模转换成阻断概率。这样 Q2 鲟类恢复和 Q5 EI 里的阻断压力使用同一套来源，避免口径不一致。

## 9. 主要输出

### 图表

- Q1 时间序列图
- Q2 鲟类有/无放流与不同阻断率恢复曲线
- Q2 江豚恢复图
- Q3 GLV 时间序列图
- Q4 极限环相图、混沌相图、分岔图
- Q5 EI 横向条形图和雷达图
- 各关键参数灵敏度曲线图

### 数据与报告

- `outputs/data/q5_ei_metrics.csv`
- `outputs/data/sensitivity_analysis.csv`
- `outputs/reports/summary.csv`
- `outputs/reports/source_library.csv`
- `outputs/reports/initial_conditions_table.csv`
- `outputs/reports/parameter_sources_table.csv`
- `outputs/reports/parameter_ranges_table.csv`

## 10. 如何解读灵敏度结果

- 如果某条灵敏度曲线斜率很大，说明模型结果对该参数更敏感，后续应优先寻找更可靠的实测或文献约束。
- 如果某条曲线近似平坦，说明该参数在当前区间内不是主要不确定性来源。
- Q2 重点看鲟类 10 年末恢复量是否随阻断或放流参数显著改变。
- Q3 重点看最大实部特征值是否跨过 0，因为这代表局部稳定性发生变化。
- Q5 重点看 `EI`、`N_final` 和 `I_final` 的方向变化，判断污染敏感参数和竞争参数对生态评价的影响。

## 11. 引用来源

### 长江背景、公报与附件

1. 长江流域水生生物资源及生境状况公报（2022年）官方发布页：[gov.cn](https://www.gov.cn/lianbo/bumen/202310/content_6908750.htm)
2. 官方通报：2022 年长江江豚数量达到 1249 头：[english.www.gov.cn](https://english.www.gov.cn/statecouncil/ministries/202302/28/content_WS63fdfe52c6d0a757729e761d.html)
3. 附件1《长江鱼类洄游路径》
4. 附件2《长江鱼类的食物链》

### 种群增长与生态动力学

1. Tian, R. C. 2006. Toward standard parameterizations in marine biological modeling. *Ecological Modelling*. DOI: [10.1016/j.ecolmodel.2005.09.003](https://doi.org/10.1016/j.ecolmodel.2005.09.003)
2. May, R. M. 1976. Simple mathematical models with very complicated dynamics. *Nature*. DOI: [10.1038/261459a0](https://doi.org/10.1038/261459a0)

### 捕食、竞争与食物链结构

1. Holling, C. S. 1959. The Components of Predation as Revealed by a Study of Small-Mammal Predation of the European Pine Sawfly. DOI: [10.4039/Ent91293-5](https://doi.org/10.4039/Ent91293-5)
2. Xu, L. et al. 2012. Spatiotemporal patterns of the fish assemblages downstream of the Gezhouba Dam on the Yangtze River. DOI: [10.1007/s11427-012-4349-0](https://doi.org/10.1007/s11427-012-4349-0)
3. Yi, Y. et al. 2010. Impact of the Gezhouba and Three Gorges Dams on habitat suitability of carps in the Yangtze River. DOI: [10.1016/j.jhydrol.2010.04.018](https://doi.org/10.1016/j.jhydrol.2010.04.018)
4. Chen, D. et al. 2019. Regime shift in fish assemblage structure in the Yangtze River following construction of the Three Gorges Dam. *Scientific Reports*. [文章链接](https://www.nature.com/articles/s41598-019-38993-x)

### 污染效应与入侵建模

1. Rubach, M. N. et al. 2006. Meta-analysis of intrinsic rates of increase and carrying capacity of populations affected by toxic and other stressors. [PubMed](https://pubmed.ncbi.nlm.nih.gov/16193755/)
2. Palmer, J. et al. 2017. Dynamic models in research and management of biological invasions. DOI: [10.1016/j.jenvman.2017.03.060](https://doi.org/10.1016/j.jenvman.2017.03.060)
3. Freshwater Fish Invasions: A Comprehensive Review. 2022. DOI: [10.1146/annurev-ecolsys-032522-015551](https://doi.org/10.1146/annurev-ecolsys-032522-015551)

## 12. 后续最值得继续完善的方向

- 用更细的监测数据替代部分相对量级初值
- 对 `p_list` 建立更明确的坝体-河段映射表
- 将灵敏度分析扩展到多参数联合不确定性分析
- 如果要投稿或答辩展示，可把 `outputs/reports/*.csv` 进一步整理成论文附录表格
