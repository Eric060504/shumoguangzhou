# ODE模型定义文件，包含所有题目的ODE函数实现，以及一些辅助函数
from __future__ import annotations

import math
from typing import Iterable

import numpy as np


def safe_state(y: Iterable[float]) -> np.ndarray:
    # 把状态变量统一截断到非负区间，避免数值误差产生负生物量。
    return np.maximum(np.asarray(y, dtype=float), 0.0)


def compute_fragmentation(p_list: Iterable[float]) -> float:
    # 按“至少被一个大坝阻断”的累计概率聚合，得到文档中的 H。
    values = np.clip(np.asarray(list(p_list), dtype=float), 0.0, 0.999)
    return 1.0 - float(np.prod(1.0 - values))


def stocking_input(t: float, schedule: list[dict], mode: str = "piecewise", width: float = 0.1) -> float:
    # 根据放流计划返回 t 时刻的放流输入强度，供鲟类方程直接调用。
    if not schedule:
        return 0.0
    if mode == "pulse":
        # 脉冲近似模式：用窄高斯代替理想 delta 脉冲，便于数值积分。
        total = 0.0
        sigma = max(width / 2.5, 1e-3)
        for item in schedule:
            total += item["magnitude"] * math.exp(-0.5 * ((t - item["time"]) / sigma) ** 2)
        return total
    # 分段注入模式：在短时间窗内按常量速率补入，更稳定也更直观。
    total = 0.0
    half_width = max(width / 2.0, 1e-6)
    for item in schedule:
        if abs(t - item["time"]) <= half_width:
            total += item["magnitude"] / width
    return total


def pollution_adjusted_params(
    K0: float,
    mu0: float,
    theta: float,
    gamma: float,
    E: float,
) -> tuple[float, float]:
    # 污染同时通过压缩环境容纳量和提高死亡率影响种群。
    K_eff = K0 * math.exp(-theta * E)
    mu_eff = mu0 + gamma * E
    return K_eff, mu_eff


def q1_plankton_filterfish_ode(t: float, y: np.ndarray, params: dict) -> np.ndarray:
    # Q1 子模型 1：描述浮游植物、浮游动物、鲢鱼、鳙鱼之间的连续时间演化。
    P, Z, L, B = safe_state(y)
    dP = (
        params["r_P"] * P * (1.0 - P / params["K_P"])
        - params["alpha_PZ"] * P * Z
        - params["alpha_PL"] * P * L
        - params["theta_P"] * params["E"] * P
    )
    dZ = (
        params["r_Z"] * Z * (1.0 - Z / params["K_Z"])
        + params["e_Z"] * params["alpha_PZ"] * P * Z
        - params["alpha_ZB"] * Z * B
        - params["theta_Z"] * params["E"] * Z
    )
    dL = (
        params["r_L"] * L * (1.0 - L / params["K_L"])
        + params["e_L"] * params["alpha_PL"] * P * L
        - params["mu_L"] * L
        - params["F_L"] * L
    )
    dB = (
        params["r_B"] * B * (1.0 - B / params["K_B"])
        + params["e_B"] * params["alpha_ZB"] * Z * B
        - params["mu_B"] * B
        - params["F_B"] * B
    )
    return np.array([dP, dZ, dL, dB], dtype=float)


def q1_macrophyte_benthic_ode(t: float, y: np.ndarray, params: dict) -> np.ndarray:
    # Q1 子模型 2：描述水草、底栖资源、草鱼、青鱼之间的营养级作用。
    V, R, C, M = safe_state(y)
    dV = (
        params["r_V"] * V * (1.0 - V / params["K_V"])
        - params["alpha_VC"] * V * C
        - params["theta_V"] * params["E"] * V
    )
    dR = (
        params["r_R"] * R * (1.0 - R / params["K_R"])
        - params["alpha_RM"] * R * M
        - params["theta_R"] * params["E"] * R
    )
    dC = (
        params["r_C"] * C * (1.0 - C / params["K_C"])
        + params["e_C"] * params["alpha_VC"] * V * C
        - params["mu_C"] * C
        - params["F_C"] * C
    )
    dM = (
        params["r_M"] * M * (1.0 - M / params["K_M"])
        + params["e_M"] * params["alpha_RM"] * R * M
        - params["mu_M"] * M
        - params["F_M"] * M
    )
    return np.array([dV, dR, dC, dM], dtype=float)


def q2_sturgeon_ode(t: float, y: np.ndarray, params: dict) -> np.ndarray:
    # Q2 鲟类模型：把繁殖恢复、阻断压力、捕捞压力和人工放流放进同一条方程。
    S = safe_state(y)[0]
    # 鲟类模型中的 H 始终由 p_list 计算，不接受外部直接给定常数 H。
    H = compute_fragmentation(params["p_list"])
    stocking = stocking_input(
        t,
        schedule=params.get("stocking_schedule", []),
        mode=params.get("stocking_mode", "piecewise"),
        width=params.get("stocking_width", 0.1),
    )
    dS = (
        params["r_S"] * S * (1.0 - S / params["K_S"])
        - params["mu_S"] * S
        - params["rho"] * H * S
        - params["F_S"] * S
        + stocking
    )
    return np.array([dS], dtype=float)


def q2_porpoise_ode(t: float, y: np.ndarray, params: dict) -> np.ndarray:
    # Q2 江豚模型：用鱼类资源供给、污染和人类扰动刻画江豚恢复能力。
    D = safe_state(y)[0]
    fish_resource = max(params["fish_resource"], 0.0)
    intake = params["alpha_D"] * D * fish_resource / (params["fish_half_sat"] + fish_resource)
    dD = (
        params["r_D"] * D * (1.0 - D / params["K_D"])
        + intake
        - params["theta_D"] * params["E"] * D
        - params["human_disturbance"] * D
    )
    return np.array([dD], dtype=float)


def q3_glv_foodweb_ode(t: float, y: np.ndarray, params: dict) -> np.ndarray:
    # Q3 广义 Lotka-Volterra 模型：面向功能群而不是单物种展开食物网动力学。
    x = safe_state(y)
    r = np.asarray(params["r"], dtype=float)
    K = np.asarray(params["K"], dtype=float)
    A = np.asarray(params["A"], dtype=float)
    F = np.asarray(params["F"], dtype=float)
    theta = np.asarray(params["theta"], dtype=float)
    E = float(params["E"])
    # A @ x 给出功能群之间的净相互作用，保留 GLV 写法便于后续 Jacobian 分析。
    interaction = A @ x
    dx = x * (r * (1.0 - x / K) + interaction - F - theta * E)
    return dx.astype(float)


def q4_cycle_chaos_ode(t: float, y: np.ndarray, params: dict) -> np.ndarray:
    # Q4 三维食物链模型：用于展示稳定极限环与混沌吸引子的相变行为。
    X, Y, Z = safe_state(y)
    # Q4 通过周期扰动 r(t) 来触发从稳定振荡到复杂动力学的过渡。
    r_t = params["r0"] * (1.0 - params["theta"] * params["E0"]) + params["eps"] * math.sin(params["omega"] * t)
    dX = r_t * X * (1.0 - X / params["K"]) - params["a"] * X * Y
    dY = params["b"] * params["a"] * X * Y - params["c"] * Y * Z - params["d"] * Y
    dZ = params["e"] * params["c"] * Y * Z - params["f"] * Z
    return np.array([dX, dY, dZ], dtype=float)


def q5_native_invasive_ode(t: float, y: np.ndarray, params: dict) -> np.ndarray:
    # Q5 竞争模型：比较本土鱼类与外来入侵种在污染背景下的竞争结果。
    N, I = safe_state(y)
    K_N, mu_N = pollution_adjusted_params(
        K0=params["K_N0"],
        mu0=params["mu_N0"],
        theta=params["theta_N"],
        gamma=params["gamma_N"],
        E=params["E"],
    )
    K_I, mu_I = pollution_adjusted_params(
        K0=params["K_I0"],
        mu0=params["mu_I0"],
        theta=params["theta_I"],
        gamma=params["gamma_I"],
        E=params["E"],
    )
    dN = params["r_N"] * N * (1.0 - (N + params["alpha_NI"] * I) / K_N) - mu_N * N
    dI = params["r_I"] * I * (1.0 - (I + params["beta_IN"] * N) / K_I) - mu_I * I
    return np.array([dN, dI], dtype=float)
