# 存储所有模型的参数、初始条件和情景定义，提供构建参数集的函数
from __future__ import annotations

from copy import deepcopy
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "outputs"
FIGURE_DIR = OUTPUT_DIR / "figures"
DATA_DIR = OUTPUT_DIR / "data"
REPORT_DIR = OUTPUT_DIR / "reports"

# 统一管理全局模拟时间设置。
# 当前主流程默认模拟 10 年，并把时间离散为 2001 个点，兼顾曲线平滑度和运行效率。
SIMULATION_YEARS = 10.0
TIME_POINTS = 2001
TIME_SPAN = (0.0, SIMULATION_YEARS)
TIME_GRID = None

# 绘图风格集中放在这里，方便统一调中文字体、分辨率和整体观感。
PLOT_STYLE = {
    "font_candidates": ["SimHei", "Microsoft YaHei", "Arial Unicode MS", "sans-serif"],
    "dpi": 180,
    "save_dpi": 300,
    "figsize": (10, 6),
    "phase_figsize": (7.5, 6.5),
    "palette": "Set2",
    "context": "talk",
}

# 这是项目级的参数说明摘要，用于快速解释“模型里哪些设定来自文档，哪些来自生态学估计”。
PARAMETER_METADATA = {
    "document_based": [
        "禁渔通过降低捕捞强度 F 进入系统",
        "污染通过降低环境容纳量 K_i、提高死亡率 mu_i 进入系统",
        "鲟类阻断率 H 为派生参数，按 H = 1 - prod(1 - p_i) 计算",
        "EI 由生物多样性、水质、保护物种恢复、食物链连通度与污染/入侵/阻断压力综合构成",
    ],
    "ecology_estimates": [
        "所有内禀增长率、摄食系数、转化效率、容纳量和污染敏感系数为生态学合理估计值",
        "Q2 放流默认为分段注入，不直接离散 delta 函数",
        "Q4 混沌示例参数为展示性选择，服务于稳定极限环/混沌吸引子可视化",
        "情景层不直接定义 H，所有阻断压力均通过 p_list 聚合得到",
    ],
}

# EI 指数的默认权重。
# 这里只控制各指标在综合评分中的相对重要性，不参与种群动力学方程本身。
DEFAULT_EI_WEIGHTS = {
    # 这里的 H 是 EI 指数里“阻断压力项”的权重，不是阻断参数本身。
    "B": 0.22,
    "Q": 0.18,
    "Rp": 0.18,
    "Cf": 0.14,
    "E": 0.10,
    "I": 0.10,
    "H": 0.08,
}

# Q1 子模型 1：浮游植物-浮游动物-鲢鱼-鳙鱼系统参数。
# 重点体现初级生产、滤食关系、污染影响和禁渔前后的捕捞差异。
Q1_PLANKTON_PARAMS = {
    "r_P": 1.10,
    "K_P": 135.0,
    "alpha_PZ": 0.020,
    "alpha_PL": 0.017,
    "theta_P": 0.16,
    "r_Z": 0.82,
    "K_Z": 92.0,
    "e_Z": 0.28,
    "alpha_ZB": 0.018,
    "theta_Z": 0.13,
    "r_L": 0.41,
    "K_L": 58.0,
    "e_L": 0.18,
    "mu_L": 0.12,
    "F_L": 0.24,
    "r_B": 0.36,
    "K_B": 52.0,
    "e_B": 0.20,
    "mu_B": 0.11,
    "F_B": 0.22,
    "E": 0.18,
}

# Q1 子模型 2：水草-底栖资源-草鱼-青鱼系统参数。
# 这一组更偏向草食链和底栖资源链，对应四大家鱼中的另外两类典型功能群。
Q1_MACROPHYTE_PARAMS = {
    "r_V": 0.92,
    "K_V": 128.0,
    "alpha_VC": 0.021,
    "theta_V": 0.10,
    "r_R": 0.72,
    "K_R": 104.0,
    "theta_R": 0.08,
    "alpha_RM": 0.015,
    "r_C": 0.34,
    "K_C": 54.0,
    "e_C": 0.22,
    "mu_C": 0.10,
    "F_C": 0.25,
    "r_M": 0.31,
    "K_M": 46.0,
    "e_M": 0.20,
    "mu_M": 0.10,
    "F_M": 0.23,
    "E": 0.14,
}

# Q2 鲟类恢复模型参数。
# 这里同时包含自然增长、死亡、阻断、捕捞和人工放流五类作用。
Q2_STURGEON_PARAMS = {
    "r_S": 0.18,
    "K_S": 42.0,
    "mu_S": 0.09,
    "rho": 0.90,
    "F_S": 0.06,
    # p_list 是唯一的原始阻断输入，后续所有 H 都由它聚合计算得到。
    "p_list": [0.18, 0.24, 0.16],
    "stocking_schedule": [
        {"time": 1.0, "magnitude": 1.8},
        {"time": 3.0, "magnitude": 2.1},
        {"time": 5.0, "magnitude": 2.4},
        {"time": 7.0, "magnitude": 2.4},
    ],
    "stocking_width": 0.12,
    "stocking_mode": "piecewise",
}

# Q2 江豚模型参数。
# 江豚响应主要通过鱼类资源供给、污染和人为干扰三个因素体现。
Q2_PORPOISE_PARAMS = {
    "r_D": 0.11,
    "K_D": 24.0,
    "alpha_D": 0.010,
    "fish_half_sat": 32.0,
    "fish_resource": 58.0,
    "E": 0.18,
    "theta_D": 0.16,
    "human_disturbance": 0.04,
}

# Q3 广义 Lotka-Volterra 食物网参数。
# labels/r/K/A/F/theta 共同定义功能群层面的增长、相互作用、捕捞和污染压力。
Q3_GLV_PARAMS = {
    "labels": ["R", "F", "T", "I"],
    "r": [0.88, 0.24, 0.15, 0.28],
    "K": [160.0, 80.0, 28.0, 36.0],
    "A": [
        [-0.010, -0.008, 0.000, -0.005],
        [0.010, -0.012, -0.010, -0.007],
        [0.000, 0.012, -0.016, -0.005],
        [0.006, -0.012, 0.000, -0.011],
    ],
    "F": [0.00, 0.11, 0.03, 0.02],
    "E": 0.18,
    "theta": [0.05, 0.08, 0.10, 0.06],
}

# Q4 稳定极限环示例参数。
# 这组参数的目标是生成“MLE < 0”的持续振荡轨迹，用于展示规则周期行为。
Q4_LIMIT_CYCLE_PARAMS = {
    "r0": 1.05,
    "K": 1.15,
    "a": 0.95,
    "b": 0.62,
    "c": 0.55,
    "d": 0.48,
    "e": 0.72,
    "f": 0.32,
    "theta": 0.05,
    "E0": 0.18,
    "eps": 0.05,
    "omega": 2.0,
}

# Q4 混沌示例参数。
# 这组参数经过展示性选择，使系统更容易出现“MLE > 0”的复杂吸引子。
Q4_CHAOS_PARAMS = {
    "r0": 1.42,
    "K": 1.0,
    "a": 1.16,
    "b": 0.78,
    "c": 0.72,
    "d": 0.78,
    "e": 0.90,
    "f": 0.18,
    "theta": 0.07,
    "E0": 0.26,
    "eps": 0.28,
    "omega": 5.6,
}

# Q5 本土种-入侵种竞争模型参数。
# 污染通过压缩容纳量、提高死亡率进入模型，竞争项决定两类群体的长期占优关系。
Q5_NATIVE_INVASIVE_PARAMS = {
    "r_N": 0.46,
    "K_N0": 95.0,
    "theta_N": 0.95,
    "mu_N0": 0.10,
    "gamma_N": 0.13,
    "alpha_NI": 0.72,
    "r_I": 0.55,
    "K_I0": 88.0,
    "theta_I": 0.38,
    "mu_I0": 0.08,
    "gamma_I": 0.05,
    "beta_IN": 0.58,
    "E": 0.18,
}

# 所有子模型的默认初始状态。
# 这些值统一解释为相对生物量或相对数量，其来源和缩放逻辑在下方元信息表中详细说明。
INITIAL_CONDITIONS = {
    "q1_plankton": [82.0, 30.0, 19.3, 15.2],
    "q1_macrophyte": [76.0, 55.0, 18.0, 13.5],
    "q2_sturgeon": [1.2],
    "q2_porpoise": [12.49],
    "q3_glv": [120.0, 46.0, 12.0, 9.0],
    "q4_cycle": [0.82, 0.36, 0.24],
    "q4_chaos": [0.83, 0.34, 0.25],
    "q5_native_invasive": [63.8, 8.5],
}

# 统一管理公报、附件和论文的索引信息。
# 后续 README 和导出的参数表都只引用 source_id，真正的题名、年份和链接都从这里查。
SOURCE_LIBRARY = {
    "bulletin_2022_gov": {
        "type": "bulletin",
        "title": "《长江流域水生生物资源及生境状况公报（2022年）》发布：长江流域水生生物资源量呈恢复态势",
        "year": 2023,
        "link": "https://www.gov.cn/lianbo/bumen/202310/content_6908750.htm",
        "note": "官方通报中提到 2022 年重点水域监测到鱼类 193 种，较 2020 年增加 25 种，完整性指数总体回升。",
    },
    "attachment_1_migration": {
        "type": "attachment",
        "title": "附件1 长江鱼类洄游路径",
        "year": 2026,
        "link": "附件1 长江鱼类洄游路径.pdf",
        "note": "用于确定鲟类和四大家鱼的洄游/阻断结构，不直接提供连续动力学参数。",
    },
    "attachment_2_foodweb": {
        "type": "attachment",
        "title": "附件2 长江鱼类的食物链",
        "year": 2026,
        "link": "附件2 长江鱼类的食物链.pdf",
        "note": "用于确定食物链方向、营养级关系和功能群映射。",
    },
    "gov_porpoise_2022": {
        "type": "bulletin",
        "title": "Number of finless porpoises in China exceeds 1,200",
        "year": 2023,
        "link": "https://english.www.gov.cn/statecouncil/ministries/202302/28/content_WS63fdfe52c6d0a757729e761d.html",
        "note": "2022 年科学考察结果显示长江江豚数量达到 1,249 头，比 2017 年增长 23.42%。",
    },
    "yi2010_carps": {
        "type": "literature",
        "title": "Impact of the Gezhouba and Three Gorges Dams on habitat suitability of carps in the Yangtze River",
        "year": 2010,
        "link": "https://doi.org/10.1016/j.jhydrol.2010.04.018",
        "note": "支持四大家鱼与大坝、产卵栖息地、水动力条件之间的关系。",
    },
    "xu2012_assemblage": {
        "type": "literature",
        "title": "Spatiotemporal patterns of the fish assemblages downstream of the Gezhouba Dam on the Yangtze River",
        "year": 2012,
        "link": "https://doi.org/10.1007/s11427-012-4349-0",
        "note": "支持阻断和水工扰动对鱼类群落结构的影响。",
    },
    "chen2019_shift": {
        "type": "literature",
        "title": "Regime shift in fish assemblage structure in the Yangtze River following construction of the Three Gorges Dam",
        "year": 2019,
        "link": "https://www.nature.com/articles/s41598-019-38993-x",
        "note": "支持大坝影响鱼类群落组成、自然繁殖和生态结构。",
    },
    "holling1959": {
        "type": "literature",
        "title": "The Components of Predation as Revealed by a Study of Small-Mammal Predation of the European Pine Sawfly",
        "year": 1959,
        "link": "https://doi.org/10.4039/Ent91293-5",
        "note": "支持捕食/摄食响应项的生态学解释。",
    },
    "may1976": {
        "type": "literature",
        "title": "Simple mathematical models with very complicated dynamics",
        "year": 1976,
        "link": "https://doi.org/10.1038/261459a0",
        "note": "支持三维食物链系统出现复杂动力学和混沌现象。",
    },
    "tian2006_parameterization": {
        "type": "literature",
        "title": "Toward standard parameterizations in marine biological modeling",
        "year": 2006,
        "link": "https://doi.org/10.1016/j.ecolmodel.2005.09.003",
        "note": "支持增长率、摄食、死亡率、转化效率等生态学参数量级和参数化原则。",
    },
    "rubach2006_toxic": {
        "type": "literature",
        "title": "Meta-analysis of intrinsic rates of increase and carrying capacity of populations affected by toxic and other stressors",
        "year": 2006,
        "link": "https://pubmed.ncbi.nlm.nih.gov/16193755/",
        "note": "支持污染通过降低增长能力/容纳量并提高额外压力进入模型。",
    },
    "palmer2017_invasion": {
        "type": "literature",
        "title": "Dynamic models in research and management of biological invasions",
        "year": 2017,
        "link": "https://doi.org/10.1016/j.jenvman.2017.03.060",
        "note": "支持入侵种竞争模型和动态参数扫描的建模思路。",
    },
    "freshwater_invasion_review": {
        "type": "literature",
        "title": "Freshwater Fish Invasions: A Comprehensive Review",
        "year": 2022,
        "link": "https://doi.org/10.1146/annurev-ecolsys-032522-015551",
        "note": "支持入侵鱼类对本土群落的竞争压力和风险设定。",
    },
}

# 为每个状态变量补充“默认值 + 来源 + 缩放逻辑”说明。
# 这里强调的是“为什么初值取这个量级”，而不是把公报原始统计量直接当成模型状态量。
INITIAL_CONDITION_METADATA = {
    "q1_plankton": {
        "P": {
            "default": 82.0,
            "unit": "相对生物量",
            "source_type": "bulletin",
            "source_id": "bulletin_2022_gov",
            "basis": "依据公报中 2020-2022 年水生生物资源恢复趋势，将初级生产者设为较高基线量级。",
            "scaling_logic": "采用 0-150 相对生物量尺度，资源恢复态势映射为 80 左右的初值。",
        },
        "Z": {
            "default": 30.0,
            "unit": "相对生物量",
            "source_type": "attachment",
            "source_id": "attachment_2_foodweb",
            "basis": "依据浮游动物位于浮游植物之上、四大家鱼之下的营养级位置设定相对量级。",
            "scaling_logic": "按食物链中层丰度低于初级生产者但高于鱼类的原则赋值。",
        },
        "L": {
            "default": 19.3,
            "unit": "相对生物量",
            "source_type": "bulletin",
            "source_id": "bulletin_2022_gov",
            "basis": "参考公报中 2022 年鱼类监测到 193 种的恢复信号，并结合四大家鱼生态地位进行缩放。",
            "scaling_logic": "用鱼类恢复趋势构造 10-25 的中等鱼类初值区间，默认取 19.3。",
        },
        "B": {
            "default": 15.2,
            "unit": "相对生物量",
            "source_type": "attachment",
            "source_id": "attachment_2_foodweb",
            "basis": "依据鳙鱼位于滤食鱼中较高营养级位置，初值略低于鲢鱼。",
            "scaling_logic": "相对于鲢鱼取 0.75-0.85 比例，默认 15.2。",
        },
    },
    "q1_macrophyte": {
        "V": {
            "default": 76.0,
            "unit": "相对生物量",
            "source_type": "bulletin",
            "source_id": "bulletin_2022_gov",
            "basis": "水生植被按资源恢复背景设为较高值。",
            "scaling_logic": "在相对尺度上略低于浮游植物，反映主河道与通江湖泊水生植被恢复不完全。",
        },
        "R": {
            "default": 55.0,
            "unit": "相对生物量",
            "source_type": "attachment",
            "source_id": "attachment_2_foodweb",
            "basis": "底栖资源为青鱼摄食基础，按中高资源水平设定。",
            "scaling_logic": "设置在 40-70 区间内，保持高于草鱼和青鱼初值。",
        },
        "C": {
            "default": 18.0,
            "unit": "相对生物量",
            "source_type": "literature",
            "source_id": "yi2010_carps",
            "basis": "结合四大家鱼产卵栖息地受水坝影响的文献背景，草鱼初值设为中等恢复状态。",
            "scaling_logic": "保持与鲢鱼同量级，但受草食资源约束略低于部分滤食群。",
        },
        "M": {
            "default": 13.5,
            "unit": "相对生物量",
            "source_type": "attachment",
            "source_id": "attachment_2_foodweb",
            "basis": "青鱼作为更高营养级鱼类，初值低于草鱼。",
            "scaling_logic": "采用草鱼初值的 70%-80%。",
        },
    },
    "q2_sturgeon": {
        "S": {
            "default": 1.2,
            "unit": "相对生物量",
            "source_type": "bulletin",
            "source_id": "bulletin_2022_gov",
            "basis": "公报指出中华鲟自然繁殖形势仍严峻，因此设为极低基线。",
            "scaling_logic": "在相对尺度下设置为接近濒危低值，不与放流数量直接等同。",
        },
    },
    "q2_porpoise": {
        "D": {
            "default": 12.49,
            "unit": "相对数量",
            "source_type": "bulletin",
            "source_id": "gov_porpoise_2022",
            "basis": "依据 2022 年江豚数量 1249 头。",
            "scaling_logic": "按 100:1 缩放为模型相对数量，得到 12.49。",
        },
    },
    "q3_glv": {
        "R": {
            "default": 120.0,
            "unit": "相对生物量",
            "source_type": "bulletin",
            "source_id": "bulletin_2022_gov",
            "basis": "资源层按公报恢复态势设定为较高基线。",
            "scaling_logic": "作为 GLV 资源层，取略高于 Q1 初级生产者的综合水平。",
        },
        "F": {
            "default": 46.0,
            "unit": "相对生物量",
            "source_type": "attachment",
            "source_id": "attachment_2_foodweb",
            "basis": "中层鱼类综合了四大家鱼和中小型鱼资源。",
            "scaling_logic": "中层鱼类设为资源层的 35%-40%。",
        },
        "T": {
            "default": 12.0,
            "unit": "相对数量",
            "source_type": "bulletin",
            "source_id": "gov_porpoise_2022",
            "basis": "顶层捕食者参考江豚恢复规模。",
            "scaling_logic": "采用与江豚相近的相对数量水平。",
        },
        "I": {
            "default": 9.0,
            "unit": "相对生物量",
            "source_type": "bulletin",
            "source_id": "bulletin_2022_gov",
            "basis": "公报提示外来物种需要警惕，因此设为中低起始量级。",
            "scaling_logic": "保证低于本土中层鱼类，但不是接近零。",
        },
    },
    "q5_native_invasive": {
        "N": {
            "default": 63.8,
            "unit": "相对生物量",
            "source_type": "bulletin",
            "source_id": "bulletin_2022_gov",
            "basis": "以长江鱼类资源恢复为基线，本土种初值高于入侵种。",
            "scaling_logic": "按鱼类群落恢复态势设为中高值。",
        },
        "I": {
            "default": 8.5,
            "unit": "相对生物量",
            "source_type": "literature",
            "source_id": "freshwater_invasion_review",
            "basis": "外来鱼类在多数河流中相对占比仍低于本土优势群体，但已形成持续压力。",
            "scaling_logic": "保持为本土种的 10%-15%。",
        },
    },
}

# 记录关键生态学参数的含义、量纲、默认值和来源类型。
# 这部分主要服务于三件事：README 说明、论文附录追溯、后续参数校准。
PARAMETER_SOURCES = {
    "Q1_PLANKTON_PARAMS": {
        "r_P": {"meaning": "浮游植物内禀增长率", "unit": "1/年", "default": 1.10, "source_type": "literature", "source_id": "tian2006_parameterization", "sensitivity": False},
        "r_Z": {"meaning": "浮游动物内禀增长率", "unit": "1/年", "default": 0.82, "source_type": "literature", "source_id": "tian2006_parameterization", "sensitivity": False},
        "alpha_PZ": {"meaning": "浮游动物摄食浮游植物强度", "unit": "1/(生物量·年)", "default": 0.020, "source_type": "literature", "source_id": "holling1959", "sensitivity": True},
        "alpha_PL": {"meaning": "鲢鱼摄食浮游植物强度", "unit": "1/(生物量·年)", "default": 0.017, "source_type": "estimated_range", "source_id": "attachment_2_foodweb", "sensitivity": True},
        "alpha_ZB": {"meaning": "鳙鱼摄食浮游动物强度", "unit": "1/(生物量·年)", "default": 0.018, "source_type": "estimated_range", "source_id": "attachment_2_foodweb", "sensitivity": True},
        "mu_L": {"meaning": "鲢鱼自然死亡率", "unit": "1/年", "default": 0.12, "source_type": "literature", "source_id": "tian2006_parameterization", "sensitivity": False},
        "mu_B": {"meaning": "鳙鱼自然死亡率", "unit": "1/年", "default": 0.11, "source_type": "literature", "source_id": "tian2006_parameterization", "sensitivity": False},
        "theta_P": {"meaning": "污染对浮游植物的敏感系数", "unit": "无量纲", "default": 0.16, "source_type": "literature", "source_id": "rubach2006_toxic", "sensitivity": False},
        "theta_Z": {"meaning": "污染对浮游动物的敏感系数", "unit": "无量纲", "default": 0.13, "source_type": "literature", "source_id": "rubach2006_toxic", "sensitivity": False},
    },
    "Q1_MACROPHYTE_PARAMS": {
        "r_V": {"meaning": "水草内禀增长率", "unit": "1/年", "default": 0.92, "source_type": "literature", "source_id": "tian2006_parameterization", "sensitivity": False},
        "r_R": {"meaning": "底栖资源恢复率", "unit": "1/年", "default": 0.72, "source_type": "estimated_range", "source_id": "attachment_2_foodweb", "sensitivity": False},
        "alpha_VC": {"meaning": "草鱼利用水草强度", "unit": "1/(生物量·年)", "default": 0.021, "source_type": "attachment", "source_id": "attachment_2_foodweb", "sensitivity": True},
        "alpha_RM": {"meaning": "青鱼利用底栖资源强度", "unit": "1/(生物量·年)", "default": 0.015, "source_type": "attachment", "source_id": "attachment_2_foodweb", "sensitivity": True},
        "theta_V": {"meaning": "污染对水草敏感系数", "unit": "无量纲", "default": 0.10, "source_type": "literature", "source_id": "rubach2006_toxic", "sensitivity": False},
    },
    "Q2_STURGEON_PARAMS": {
        "r_S": {"meaning": "鲟类内禀增长率", "unit": "1/年", "default": 0.18, "source_type": "estimated_range", "source_id": "xu2012_assemblage", "sensitivity": False},
        "rho": {"meaning": "阻断压力放大系数", "unit": "无量纲", "default": 0.90, "source_type": "estimated_range", "source_id": "attachment_1_migration", "sensitivity": True},
        "F_S": {"meaning": "鲟类捕捞强度", "unit": "1/年", "default": 0.06, "source_type": "bulletin", "source_id": "bulletin_2022_gov", "sensitivity": False},
        "p_list": {"meaning": "各关键大坝阻断概率", "unit": "概率", "default": [0.18, 0.24, 0.16], "source_type": "attachment", "source_id": "attachment_1_migration", "sensitivity": True},
        "stocking_width": {"meaning": "单次放流作用时间窗", "unit": "年", "default": 0.12, "source_type": "estimated_range", "source_id": "gov_porpoise_2022", "sensitivity": True},
        "stocking_schedule": {"meaning": "放流时刻与放流强度", "unit": "相对数量/年", "default": "4次分段注入", "source_type": "estimated_range", "source_id": "gov_porpoise_2022", "sensitivity": True},
    },
    "Q2_PORPOISE_PARAMS": {
        "r_D": {"meaning": "江豚恢复率", "unit": "1/年", "default": 0.11, "source_type": "bulletin", "source_id": "gov_porpoise_2022", "sensitivity": False},
        "theta_D": {"meaning": "污染对江豚敏感系数", "unit": "无量纲", "default": 0.16, "source_type": "literature", "source_id": "rubach2006_toxic", "sensitivity": False},
        "human_disturbance": {"meaning": "航运和人类活动扰动强度", "unit": "1/年", "default": 0.04, "source_type": "estimated_range", "source_id": "bulletin_2022_gov", "sensitivity": False},
    },
    "Q3_GLV_PARAMS": {
        "r": {"meaning": "功能群增长率向量", "unit": "1/年", "default": [0.88, 0.24, 0.15, 0.28], "source_type": "literature", "source_id": "tian2006_parameterization", "sensitivity": False},
        "A": {"meaning": "功能群相互作用矩阵", "unit": "1/(生物量·年)", "default": "4x4 矩阵", "source_type": "estimated_range", "source_id": "attachment_2_foodweb", "sensitivity": True},
        "F": {"meaning": "功能群外部损失项", "unit": "1/年", "default": [0.00, 0.11, 0.03, 0.02], "source_type": "bulletin", "source_id": "bulletin_2022_gov", "sensitivity": False},
        "theta": {"meaning": "功能群污染敏感系数", "unit": "无量纲", "default": [0.05, 0.08, 0.10, 0.06], "source_type": "literature", "source_id": "rubach2006_toxic", "sensitivity": False},
    },
    "Q4_LIMIT_CYCLE_PARAMS": {
        "eps": {"meaning": "周期扰动幅值", "unit": "无量纲", "default": 0.05, "source_type": "literature", "source_id": "may1976", "sensitivity": False},
        "omega": {"meaning": "周期扰动频率", "unit": "1/年", "default": 2.0, "source_type": "literature", "source_id": "may1976", "sensitivity": False},
    },
    "Q4_CHAOS_PARAMS": {
        "eps": {"meaning": "周期扰动幅值", "unit": "无量纲", "default": 0.28, "source_type": "literature", "source_id": "may1976", "sensitivity": False},
        "omega": {"meaning": "周期扰动频率", "unit": "1/年", "default": 5.6, "source_type": "literature", "source_id": "may1976", "sensitivity": False},
    },
    "Q5_NATIVE_INVASIVE_PARAMS": {
        "r_N": {"meaning": "本土鱼类增长率", "unit": "1/年", "default": 0.46, "source_type": "literature", "source_id": "freshwater_invasion_review", "sensitivity": False},
        "r_I": {"meaning": "入侵鱼类增长率", "unit": "1/年", "default": 0.55, "source_type": "literature", "source_id": "freshwater_invasion_review", "sensitivity": False},
        "theta_N": {"meaning": "污染对本土种容纳量压缩强度", "unit": "无量纲", "default": 0.95, "source_type": "literature", "source_id": "rubach2006_toxic", "sensitivity": True},
        "theta_I": {"meaning": "污染对入侵种容纳量压缩强度", "unit": "无量纲", "default": 0.38, "source_type": "estimated_range", "source_id": "freshwater_invasion_review", "sensitivity": True},
        "gamma_N": {"meaning": "污染提高本土种死亡率的强度", "unit": "1/年", "default": 0.13, "source_type": "literature", "source_id": "rubach2006_toxic", "sensitivity": True},
        "gamma_I": {"meaning": "污染提高入侵种死亡率的强度", "unit": "1/年", "default": 0.05, "source_type": "estimated_range", "source_id": "freshwater_invasion_review", "sensitivity": True},
        "alpha_NI": {"meaning": "入侵种对本土种竞争系数", "unit": "无量纲", "default": 0.72, "source_type": "literature", "source_id": "palmer2017_invasion", "sensitivity": True},
        "beta_IN": {"meaning": "本土种对入侵种竞争系数", "unit": "无量纲", "default": 0.58, "source_type": "literature", "source_id": "palmer2017_invasion", "sensitivity": True},
    },
}

# 对缺乏长江专属精确标定的关键参数给出区间。
# 这些区间既是“不确定性声明”，也是灵敏度分析扫描的直接输入。
PARAMETER_RANGES = {
    "Q2_STURGEON_PARAMS": {
        "rho": {"default": 0.90, "lower": 0.60, "upper": 1.20, "reason": "阻断对繁殖成功率和洄游损失的放大程度缺乏统一实测值", "method": "single_factor_scan"},
        "stocking_width": {"default": 0.12, "lower": 0.06, "upper": 0.20, "reason": "放流脉冲作用时间窗取决于数值近似和管理节律", "method": "single_factor_scan"},
        "stocking_scale": {"default": 1.00, "lower": 0.60, "upper": 1.40, "reason": "放流强度可因年度投放规模差异而变化", "method": "single_factor_scan"},
        "p_list_scale": {"default": 1.00, "lower": 0.70, "upper": 1.30, "reason": "阻断概率受坝群联合作用和通道可通过性不确定性影响", "method": "single_factor_scan"},
    },
    "Q3_GLV_PARAMS": {
        "A_10": {"default": 0.010, "lower": 0.006, "upper": 0.016, "reason": "资源层对中层鱼类补给强度缺乏长江专属标定", "method": "single_factor_scan"},
        "A_12": {"default": -0.010, "lower": -0.016, "upper": -0.006, "reason": "顶层捕食对中层鱼类压制强度具有明显不确定性", "method": "single_factor_scan"},
    },
    "Q5_NATIVE_INVASIVE_PARAMS": {
        "theta_N": {"default": 0.95, "lower": 0.70, "upper": 1.20, "reason": "本土鱼类对污染导致容纳量下降的敏感程度具有区域差异", "method": "single_factor_scan"},
        "theta_I": {"default": 0.38, "lower": 0.20, "upper": 0.60, "reason": "入侵种通常对环境扰动更耐受，但缺乏统一定量值", "method": "single_factor_scan"},
        "gamma_N": {"default": 0.13, "lower": 0.08, "upper": 0.18, "reason": "污染诱导额外死亡率受污染类型和鱼类耐受性影响", "method": "single_factor_scan"},
        "gamma_I": {"default": 0.05, "lower": 0.02, "upper": 0.10, "reason": "入侵种额外死亡率预期小于本土种，但存在不确定性", "method": "single_factor_scan"},
        "alpha_NI": {"default": 0.72, "lower": 0.50, "upper": 0.95, "reason": "入侵竞争压力与入侵种功能性状和栖位重叠有关", "method": "single_factor_scan"},
        "beta_IN": {"default": 0.58, "lower": 0.40, "upper": 0.80, "reason": "本土种反制入侵种能力缺乏长江统一量化结果", "method": "single_factor_scan"},
    },
}

# 明确规定“哪些参数需要做灵敏度分析、怎么扫、看什么响应量”。
# 主程序会逐条读取这个清单，自动生成 CSV 和曲线图，因此这里相当于灵敏度分析的任务配置表。
SENSITIVITY_SPECS = [
    {"name": "Q2_rho", "module": "Q2", "kind": "scalar", "param_group": "q2_sturgeon", "param_name": "rho", "values": [0.60, 0.75, 0.90, 1.05, 1.20], "response": "sturgeon_final"},
    {"name": "Q2_stocking_width", "module": "Q2", "kind": "scalar", "param_group": "q2_sturgeon", "param_name": "stocking_width", "values": [0.06, 0.10, 0.12, 0.16, 0.20], "response": "sturgeon_final"},
    {"name": "Q2_stocking_scale", "module": "Q2", "kind": "stocking_scale", "param_group": "q2_sturgeon", "param_name": "stocking_scale", "values": [0.60, 0.80, 1.00, 1.20, 1.40], "response": "sturgeon_final"},
    {"name": "Q2_p_list_scale", "module": "Q2", "kind": "list_scale", "param_group": "q2_sturgeon", "param_name": "p_list", "values": [0.70, 0.85, 1.00, 1.15, 1.30], "response": "sturgeon_final"},
    {"name": "Q3_A_10", "module": "Q3", "kind": "matrix_entry", "param_group": "q3_glv", "param_name": "A", "index": [1, 0], "values": [0.006, 0.008, 0.010, 0.013, 0.016], "response": "max_real_eigen"},
    {"name": "Q3_A_12", "module": "Q3", "kind": "matrix_entry", "param_group": "q3_glv", "param_name": "A", "index": [1, 2], "values": [-0.016, -0.013, -0.010, -0.008, -0.006], "response": "max_real_eigen"},
    {"name": "Q5_theta_N", "module": "Q5", "kind": "scalar", "param_group": "q5_native_invasive", "param_name": "theta_N", "values": [0.70, 0.82, 0.95, 1.08, 1.20], "response": "EI"},
    {"name": "Q5_theta_I", "module": "Q5", "kind": "scalar", "param_group": "q5_native_invasive", "param_name": "theta_I", "values": [0.20, 0.30, 0.38, 0.48, 0.60], "response": "EI"},
    {"name": "Q5_gamma_N", "module": "Q5", "kind": "scalar", "param_group": "q5_native_invasive", "param_name": "gamma_N", "values": [0.08, 0.10, 0.13, 0.15, 0.18], "response": "EI"},
    {"name": "Q5_gamma_I", "module": "Q5", "kind": "scalar", "param_group": "q5_native_invasive", "param_name": "gamma_I", "values": [0.02, 0.04, 0.05, 0.08, 0.10], "response": "EI"},
    {"name": "Q5_alpha_NI", "module": "Q5", "kind": "scalar", "param_group": "q5_native_invasive", "param_name": "alpha_NI", "values": [0.50, 0.61, 0.72, 0.84, 0.95], "response": "EI"},
    {"name": "Q5_beta_IN", "module": "Q5", "kind": "scalar", "param_group": "q5_native_invasive", "param_name": "beta_IN", "values": [0.40, 0.49, 0.58, 0.69, 0.80], "response": "EI"},
]

# 情景层只负责描述“在基线参数上做哪些管理或压力扰动”。
# 所有场景都尽量保持含义单一，例如污染、入侵、放流、阻断背景分别管理，避免语义混杂。
SCENARIOS = {
    "baseline": {
        "description": "禁渔基线恢复情景",
        "q1_plankton": {"F_L": 0.05, "F_B": 0.05, "E": 0.16},
        "q1_macrophyte": {"F_C": 0.05, "F_M": 0.05, "E": 0.12},
        "q2_sturgeon": {"F_S": 0.02},
        "q2_porpoise": {"E": 0.16, "human_disturbance": 0.03},
        "q3_glv": {"F": [0.00, 0.04, 0.02, 0.02], "E": 0.15},
        "q5_native_invasive": {"E": 0.16},
    },
    "ban_no_pollution": {
        "description": "禁渔+无污染",
        "q5_native_invasive": {"E": 0.02},
        # EI 场景不再手写 H，而是引用一个阻断背景情景来统一计算。
        "ei_context": {"E": 0.02, "I_pressure": 0.22, "fragmentation_scenario": "sturgeon_medium_H"},
    },
    "ban_high_pollution": {
        "description": "禁渔+高污染",
        "q5_native_invasive": {"E": 0.58},
        "ei_context": {"E": 0.58, "I_pressure": 0.22, "fragmentation_scenario": "sturgeon_medium_H"},
    },
    "ban_invasion": {
        "description": "禁渔+入侵物种压力",
        "q5_native_invasive": {"E": 0.18, "r_I": 0.74, "alpha_NI": 0.96},
        "ei_context": {"E": 0.18, "I_pressure": 0.68, "fragmentation_scenario": "sturgeon_medium_H"},
    },
    "sturgeon_stocking_off": {
        "description": "鲟类无放流",
        "q2_sturgeon": {"stocking_schedule": []},
    },
    "sturgeon_stocking_on": {
        "description": "鲟类有放流",
        "q2_sturgeon": {"stocking_schedule": deepcopy(Q2_STURGEON_PARAMS["stocking_schedule"])},
    },
    "sturgeon_low_H": {
        "description": "鲟类低阻断率",
        "q2_sturgeon": {"p_list": [0.08, 0.10, 0.06]},
    },
    "sturgeon_medium_H": {
        "description": "鲟类中阻断率",
        "q2_sturgeon": {"p_list": [0.18, 0.24, 0.16]},
    },
    "sturgeon_high_H": {
        "description": "鲟类高阻断率",
        "q2_sturgeon": {"p_list": [0.38, 0.42, 0.30]},
    },
}


def build_base_parameter_sets() -> dict[str, dict]:
    # 返回一份可安全修改的基础参数副本，避免后续情景覆盖时污染原始默认参数。
    return {
        "q1_plankton": deepcopy(Q1_PLANKTON_PARAMS),
        "q1_macrophyte": deepcopy(Q1_MACROPHYTE_PARAMS),
        "q2_sturgeon": deepcopy(Q2_STURGEON_PARAMS),
        "q2_porpoise": deepcopy(Q2_PORPOISE_PARAMS),
        "q3_glv": deepcopy(Q3_GLV_PARAMS),
        "q4_limit_cycle": deepcopy(Q4_LIMIT_CYCLE_PARAMS),
        "q4_chaos": deepcopy(Q4_CHAOS_PARAMS),
        "q5_native_invasive": deepcopy(Q5_NATIVE_INVASIVE_PARAMS),
    }


def merge_params(base: dict, updates: dict | None) -> dict:
    # 用浅层键覆盖的方式组合“基础参数 + 情景修正项”，满足本项目当前配置结构。
    merged = deepcopy(base)
    if updates:
        merged.update(deepcopy(updates))
    return merged
