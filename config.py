# 存储所有模型的参数、初始条件和情景定义，提供构建参数集的函数
from __future__ import annotations

from copy import deepcopy
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "outputs"
FIGURE_DIR = OUTPUT_DIR / "figures"
DATA_DIR = OUTPUT_DIR / "data"
REPORT_DIR = OUTPUT_DIR / "reports"

SIMULATION_YEARS = 10.0
TIME_POINTS = 2001
TIME_SPAN = (0.0, SIMULATION_YEARS)
TIME_GRID = None

PLOT_STYLE = {
    "font_candidates": ["SimHei", "Microsoft YaHei", "Arial Unicode MS", "sans-serif"],
    "dpi": 180,
    "save_dpi": 300,
    "figsize": (10, 6),
    "phase_figsize": (7.5, 6.5),
    "palette": "Set2",
    "context": "talk",
}

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

INITIAL_CONDITIONS = {
    "q1_plankton": [72.0, 28.0, 16.0, 12.0],
    "q1_macrophyte": [68.0, 52.0, 15.0, 11.0],
    "q2_sturgeon": [4.2],
    "q2_porpoise": [6.2],
    "q3_glv": [115.0, 42.0, 10.0, 8.0],
    "q4_cycle": [0.82, 0.36, 0.24],
    "q4_chaos": [0.83, 0.34, 0.25],
    "q5_native_invasive": [58.0, 10.0],
}

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
    merged = deepcopy(base)
    if updates:
        merged.update(deepcopy(updates))
    return merged
