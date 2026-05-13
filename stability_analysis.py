# 包含 Jacobian 矩阵计算、平衡点求解、最大 Lyapunov 指数计算、参数扫描等功能的稳定性分析工具函数，支持符号计算和数值分析两种方式，适用于各种 ODE 模型的局部稳定性和混沌行为分析
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve


@dataclass
class EigenResult:
    # 统一封装特征值分析结果，便于主程序直接写入摘要报告。
    eigenvalues: np.ndarray
    max_real_part: float
    is_locally_stable: bool


def solve_equilibrium(func, y0: np.ndarray, params: dict) -> np.ndarray:
    # 用数值根求解器寻找平衡点，并截断负值以保持生态量含义。
    y0 = np.asarray(y0, dtype=float)
    solution = fsolve(lambda y: func(0.0, y, params), y0)
    return np.maximum(solution, 0.0)


def symbolic_jacobian(state_symbols: list[sp.Symbol], rhs_exprs: list[sp.Expr]) -> sp.Matrix:
    # 根据符号形式的状态方程自动构造 Jacobian，避免手工求偏导出错。
    return sp.Matrix(rhs_exprs).jacobian(sp.Matrix(state_symbols))


def evaluate_jacobian(
    jacobian_expr: sp.Matrix,
    state_symbols: list[sp.Symbol],
    equilibrium: np.ndarray,
    substitutions: dict | None = None,
) -> np.ndarray:
    # 先代入参数，再代入平衡点，得到可直接做特征值分析的数值 Jacobian。
    eq_map = {symbol: float(value) for symbol, value in zip(state_symbols, equilibrium)}
    full_map = {}
    if substitutions:
        full_map.update(substitutions)
    full_map.update(eq_map)
    evaluated = np.array(jacobian_expr.evalf(subs=full_map), dtype=float)
    return evaluated


def eigen_analysis(matrix: np.ndarray) -> EigenResult:
    # 提取特征值并用最大实部判断局部稳定性，是稳定性分析的直接输出。
    eigenvalues = np.linalg.eigvals(matrix)
    max_real = float(np.max(np.real(eigenvalues)))
    return EigenResult(
        eigenvalues=eigenvalues,
        max_real_part=max_real,
        is_locally_stable=max_real < 0.0,
    )


def max_lyapunov_exponent(
    func,
    y0: np.ndarray,
    params: dict,
    t_span: tuple[float, float],
    dt: float = 0.05,
    delta0: float = 1e-8,
    transient_fraction: float = 0.2,
) -> float:
    # 双轨道重标定法：持续比较参考轨道与微扰轨道的发散速度。
    y = np.asarray(y0, dtype=float)
    direction = np.ones_like(y, dtype=float)
    direction /= np.linalg.norm(direction)
    y_perturbed = y + delta0 * direction
    times = np.arange(t_span[0], t_span[1], dt)
    transient_time = t_span[0] + transient_fraction * (t_span[1] - t_span[0])
    sum_logs = 0.0
    count = 0

    for start, end in zip(times[:-1], times[1:]):
        sol_ref = solve_ivp(func, (start, end), y, args=(params,), t_eval=[end], method="RK45")
        sol_pert = solve_ivp(func, (start, end), y_perturbed, args=(params,), t_eval=[end], method="RK45")
        y = sol_ref.y[:, -1]
        y_perturbed = sol_pert.y[:, -1]
        diff = y_perturbed - y
        dist = np.linalg.norm(diff)
        if dist <= 1e-16:
            dist = 1e-16
        if end >= transient_time:
            sum_logs += np.log(dist / delta0)
            count += 1
        # 每一步把微扰长度重新拉回 delta0，避免误差被无限放大或衰减。
        y_perturbed = y + delta0 * diff / dist

    if count == 0:
        return float("nan")
    return sum_logs / (count * dt)


def parameter_sweep(
    func,
    y0: np.ndarray,
    base_params: dict,
    sweep_key: str,
    values: list[float],
    t_span: tuple[float, float],
    summary_fn=None,
    t_eval: np.ndarray | None = None,
) -> list[dict]:
    # 对单个关键参数做逐点扫描，支撑恢复曲线、分岔图和灵敏度比较。
    results = []
    for value in values:
        params = dict(base_params)
        params[sweep_key] = value
        sol = solve_ivp(func, t_span, y0, args=(params,), t_eval=t_eval, method="RK45")
        record = {"parameter": value, "solution": sol}
        if summary_fn is not None:
            record.update(summary_fn(sol, params))
        results.append(record)
    return results
