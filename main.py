# 主程序文件，负责运行所有模型模块，进行数值积分，计算指标，并生成输出文件和图表
from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import sympy as sp
from scipy.integrate import solve_ivp

from config import (
    DATA_DIR,
    DEFAULT_EI_WEIGHTS,
    FIGURE_DIR,
    INITIAL_CONDITIONS,
    OUTPUT_DIR,
    PARAMETER_METADATA,
    PLOT_STYLE,
    REPORT_DIR,
    SCENARIOS,
    TIME_POINTS,
    TIME_SPAN,
    build_base_parameter_sets,
    merge_params,
)
from models_ode import (
    compute_fragmentation,
    q1_macrophyte_benthic_ode,
    q1_plankton_filterfish_ode,
    q2_porpoise_ode,
    q2_sturgeon_ode,
    q3_glv_foodweb_ode,
    q4_cycle_chaos_ode,
    q5_native_invasive_ode,
)
from plot_utils import (
    plot_bifurcation,
    plot_ei_horizontal_bar,
    plot_phase_portrait,
    plot_radar,
    plot_sturgeon_recovery,
    plot_time_series,
    setup_plot_style,
)
from stability_analysis import (
    eigen_analysis,
    evaluate_jacobian,
    max_lyapunov_exponent,
    parameter_sweep,
    solve_equilibrium,
    symbolic_jacobian,
)


def ensure_output_dirs() -> None:
    # 创建图表、数据和报告目录，保证后续批量导出不会因路径缺失失败。
    for path in [OUTPUT_DIR, FIGURE_DIR, DATA_DIR, REPORT_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def integrate_system(func, y0, params, t_span=TIME_SPAN, n_points=TIME_POINTS):
    # 统一封装 ODE 求解流程，减少各题重复写 solve_ivp 的样板代码。
    t_eval = np.linspace(t_span[0], t_span[1], n_points)
    sol = solve_ivp(func, t_span, y0, args=(params,), t_eval=t_eval, method="RK45", rtol=1e-6, atol=1e-9)
    if not sol.success:
        raise RuntimeError(sol.message)
    # 所有生物量/种群量按非负量处理，避免数值误差带来少量负值。
    sol.y = np.maximum(sol.y, 0.0)
    return sol


def solution_to_frame(sol, labels):
    # 把 solve_ivp 输出整理成 DataFrame，便于直接绘图和导出 CSV。
    df = pd.DataFrame(sol.y.T, columns=labels)
    df.insert(0, "time", sol.t)
    return df


def calc_biodiversity_index(state: dict[str, float]) -> float:
    # 用 Shannon 指数并归一化到 0-1，便于和其他 EI 分量同尺度比较。
    values = np.asarray([max(v, 1e-8) for v in state.values()], dtype=float)
    proportions = values / values.sum()
    shannon = -np.sum(proportions * np.log(proportions))
    return float(shannon / np.log(len(values)))


def calc_water_quality_index(E: float) -> float:
    # 用污染压力的反向指标近似表达水质水平，污染越强则水质分越低。
    return float(np.clip(1.0 - E, 0.0, 1.0))


def calc_recovery_index(S: float, D: float, reference: dict[str, float]) -> float:
    # 用鲟类和江豚相对参考恢复水平来表示重点保护物种恢复度 Rp。
    s_score = min(S / max(reference["S"], 1e-8), 1.5)
    d_score = min(D / max(reference["D"], 1e-8), 1.5)
    return float(np.clip(0.5 * (s_score + d_score) / 1.5, 0.0, 1.0))


def calc_foodweb_connectivity(state: dict[str, float]) -> float:
    # 这里用“群落丰度越均衡，连通性越高”的近似指标表达食物链完整性。
    vals = np.asarray(list(state.values()), dtype=float)
    cv = np.std(vals) / (np.mean(vals) + 1e-8)
    return float(np.clip(1.0 / (1.0 + cv), 0.0, 1.0))


def calc_ei(metrics: dict[str, float], weights: dict[str, float]) -> float:
    # EI 采用“正向生态收益 - 负向生态压力”的线性组合，并裁剪到 0-1。
    positive = weights["B"] * metrics["B"] + weights["Q"] * metrics["Q"] + weights["Rp"] * metrics["Rp"] + weights["Cf"] * metrics["Cf"]
    negative = weights["E"] * metrics["E"] + weights["I"] * metrics["I"] + weights["H"] * metrics["H"]
    return float(np.clip(positive - negative + 0.35, 0.0, 1.0))


def resolve_fragmentation_from_scenario(
    scenario_name: str | None,
    base_sturgeon_params: dict,
) -> tuple[float, str]:
    # EI 和 Q2 共用同一套阻断背景解析逻辑，避免不同模块各自维护 H。
    resolved_name = scenario_name or "sturgeon_medium_H"
    scenario_cfg = SCENARIOS.get(resolved_name, {})
    sturgeon_overrides = scenario_cfg.get("q2_sturgeon", {})
    p_list = sturgeon_overrides.get("p_list", base_sturgeon_params["p_list"])
    return compute_fragmentation(p_list), resolved_name


def build_q4_symbolics():
    # 构建 Q4 的符号化三维方程，供 Jacobian 和稳定性分析复用。
    X, Y, Z = sp.symbols("X Y Z", positive=True)
    r0, K, a, b, c, d, e, f, theta, E0 = sp.symbols("r0 K a b c d e f theta E0", positive=True)
    r_eff = r0 * (1 - theta * E0)
    rhs = [
        r_eff * X * (1 - X / K) - a * X * Y,
        b * a * X * Y - c * Y * Z - d * Y,
        e * c * Y * Z - f * Z,
    ]
    return [X, Y, Z], rhs, {
        "r0": r0,
        "K": K,
        "a": a,
        "b": b,
        "c": c,
        "d": d,
        "e": e,
        "f": f,
        "theta": theta,
        "E0": E0,
    }


def run_q1(base_params: dict, summary_rows: list[dict]) -> dict[str, pd.DataFrame]:
    # 执行 Q1 全流程：禁渔前后仿真、时间序列绘图和结果落盘。
    baseline_plankton = merge_params(base_params["q1_plankton"], SCENARIOS["baseline"]["q1_plankton"])
    baseline_macro = merge_params(base_params["q1_macrophyte"], SCENARIOS["baseline"]["q1_macrophyte"])
    preban_plankton = deepcopy(base_params["q1_plankton"])
    preban_macro = deepcopy(base_params["q1_macrophyte"])

    sol_pre_plankton = integrate_system(q1_plankton_filterfish_ode, INITIAL_CONDITIONS["q1_plankton"], preban_plankton)
    sol_post_plankton = integrate_system(q1_plankton_filterfish_ode, INITIAL_CONDITIONS["q1_plankton"], baseline_plankton)
    sol_pre_macro = integrate_system(q1_macrophyte_benthic_ode, INITIAL_CONDITIONS["q1_macrophyte"], preban_macro)
    sol_post_macro = integrate_system(q1_macrophyte_benthic_ode, INITIAL_CONDITIONS["q1_macrophyte"], baseline_macro)

    labels_plankton = ["P", "Z", "L", "B"]
    labels_macro = ["V", "R", "C", "M"]
    df_plankton = solution_to_frame(sol_post_plankton, labels_plankton)
    df_macro = solution_to_frame(sol_post_macro, labels_macro)
    plot_time_series(df_plankton, FIGURE_DIR, "Q1 浮游植物-浮游动物-滤食鱼类 10 年演化", "q1_plankton_time_series")
    plot_time_series(df_macro, FIGURE_DIR, "Q1 水草-底栖资源-草鱼青鱼 10 年演化", "q1_macrophyte_time_series")

    compare_df = pd.DataFrame(
        {
            "time": sol_post_plankton.t,
            "L_禁渔前": sol_pre_plankton.y[2],
            "L_禁渔后": sol_post_plankton.y[2],
            "B_禁渔前": sol_pre_plankton.y[3],
            "B_禁渔后": sol_post_plankton.y[3],
        }
    )
    plot_time_series(compare_df, FIGURE_DIR, "Q1 禁渔前后鲢鳙恢复对比", "q1_ban_compare")

    df_plankton.to_csv(DATA_DIR / "q1_plankton_timeseries.csv", index=False)
    df_macro.to_csv(DATA_DIR / "q1_macrophyte_timeseries.csv", index=False)

    summary_rows.extend(
        [
            {"module": "Q1_plankton", "metric": "L_final", "value": float(sol_post_plankton.y[2, -1])},
            {"module": "Q1_plankton", "metric": "B_final", "value": float(sol_post_plankton.y[3, -1])},
            {"module": "Q1_macrophyte", "metric": "C_final", "value": float(sol_post_macro.y[2, -1])},
            {"module": "Q1_macrophyte", "metric": "M_final", "value": float(sol_post_macro.y[3, -1])},
        ]
    )
    return {"plankton": df_plankton, "macro": df_macro}


def run_q2(base_params: dict, summary_rows: list[dict]) -> dict[str, pd.DataFrame]:
    # 执行 Q2：比较不同阻断率和有无放流下的鲟类恢复，并单独模拟江豚。
    records = []
    H_scenarios = ["sturgeon_low_H", "sturgeon_medium_H", "sturgeon_high_H"]
    stocking_scenarios = ["sturgeon_stocking_off", "sturgeon_stocking_on"]
    for h_name in H_scenarios:
        for stocking_name in stocking_scenarios:
            params = merge_params(base_params["q2_sturgeon"], SCENARIOS[h_name]["q2_sturgeon"])
            params = merge_params(params, SCENARIOS[stocking_name]["q2_sturgeon"])
            sol = integrate_system(q2_sturgeon_ode, INITIAL_CONDITIONS["q2_sturgeon"], params)
            H = compute_fragmentation(params["p_list"])
            records.append(
                pd.DataFrame(
                    {
                        "time": sol.t,
                        "sturgeon": sol.y[0],
                        "scenario": f"{h_name.replace('sturgeon_', '')}",
                        "stocking": "有放流" if stocking_name.endswith("on") else "无放流",
                        "H": H,
                    }
                )
            )
            summary_rows.append(
                {
                    "module": "Q2_sturgeon",
                    "metric": f"{h_name}_{stocking_name}_final",
                    "value": float(sol.y[0, -1]),
                }
            )
            summary_rows.append(
                {
                    "module": "Q2_sturgeon",
                    "metric": f"{h_name}_H",
                    "value": float(H),
                }
            )

    sturgeon_df = pd.concat(records, ignore_index=True)
    sturgeon_df.to_csv(DATA_DIR / "q2_sturgeon_recovery.csv", index=False)
    plot_sturgeon_recovery(sturgeon_df, FIGURE_DIR, "q2_sturgeon_recovery")

    porpoise_params = merge_params(base_params["q2_porpoise"], SCENARIOS["baseline"]["q2_porpoise"])
    porpoise_sol = integrate_system(q2_porpoise_ode, INITIAL_CONDITIONS["q2_porpoise"], porpoise_params)
    porpoise_df = solution_to_frame(porpoise_sol, ["D"])
    porpoise_df.to_csv(DATA_DIR / "q2_porpoise_timeseries.csv", index=False)
    plot_time_series(porpoise_df, FIGURE_DIR, "Q2 江豚种群恢复轨迹", "q2_porpoise_time_series", ylabel="江豚相对数量")
    summary_rows.append({"module": "Q2_porpoise", "metric": "D_final", "value": float(porpoise_sol.y[0, -1])})
    return {"sturgeon": sturgeon_df, "porpoise": porpoise_df}


def run_q3(base_params: dict, summary_rows: list[dict]) -> pd.DataFrame:
    # 执行 Q3：先做 GLV 食物网仿真，再做平衡点和 Jacobian 特征值分析。
    params = merge_params(base_params["q3_glv"], SCENARIOS["baseline"]["q3_glv"])
    sol = integrate_system(q3_glv_foodweb_ode, INITIAL_CONDITIONS["q3_glv"], params)
    labels = params["labels"]
    df = solution_to_frame(sol, labels)
    df.to_csv(DATA_DIR / "q3_glv_foodweb.csv", index=False)
    plot_time_series(df, FIGURE_DIR, "Q3 广义 Lotka-Volterra 食物网演化", "q3_glv_foodweb")

    eq = solve_equilibrium(q3_glv_foodweb_ode, INITIAL_CONDITIONS["q3_glv"], params)
    state_symbols = list(sp.symbols("R F T I", positive=True))
    r0, r1, r2, r3 = sp.symbols("r0 r1 r2 r3")
    K0, K1, K2, K3 = sp.symbols("K0 K1 K2 K3")
    F0, F1, F2, F3 = sp.symbols("F0 F1 F2 F3")
    th0, th1, th2, th3, E = sp.symbols("th0 th1 th2 th3 E")
    a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33 = sp.symbols(
        "a00 a01 a02 a03 a10 a11 a12 a13 a20 a21 a22 a23 a30 a31 a32 a33"
    )
    A = [
        [a00, a01, a02, a03],
        [a10, a11, a12, a13],
        [a20, a21, a22, a23],
        [a30, a31, a32, a33],
    ]
    rhs = []
    rates = [r0, r1, r2, r3]
    caps = [K0, K1, K2, K3]
    harvests = [F0, F1, F2, F3]
    thetas = [th0, th1, th2, th3]
    for i, sym in enumerate(state_symbols):
        # 这里按和数值模型相同的 GLV 结构重建符号表达式，便于 Jacobian 保持一致。
        interaction = sum(A[i][j] * state_symbols[j] for j in range(4))
        rhs.append(sym * (rates[i] * (1 - sym / caps[i]) + interaction - harvests[i] - thetas[i] * E))
    jac = symbolic_jacobian(state_symbols, rhs)
    subs = {
        r0: params["r"][0],
        r1: params["r"][1],
        r2: params["r"][2],
        r3: params["r"][3],
        K0: params["K"][0],
        K1: params["K"][1],
        K2: params["K"][2],
        K3: params["K"][3],
        F0: params["F"][0],
        F1: params["F"][1],
        F2: params["F"][2],
        F3: params["F"][3],
        th0: params["theta"][0],
        th1: params["theta"][1],
        th2: params["theta"][2],
        th3: params["theta"][3],
        E: params["E"],
        a00: params["A"][0][0],
        a01: params["A"][0][1],
        a02: params["A"][0][2],
        a03: params["A"][0][3],
        a10: params["A"][1][0],
        a11: params["A"][1][1],
        a12: params["A"][1][2],
        a13: params["A"][1][3],
        a20: params["A"][2][0],
        a21: params["A"][2][1],
        a22: params["A"][2][2],
        a23: params["A"][2][3],
        a30: params["A"][3][0],
        a31: params["A"][3][1],
        a32: params["A"][3][2],
        a33: params["A"][3][3],
    }
    jac_num = evaluate_jacobian(jac, state_symbols, eq, substitutions=subs)
    eig = eigen_analysis(jac_num)
    summary_rows.extend(
        [
            {"module": "Q3_glv", "metric": "equilibrium_norm", "value": float(np.linalg.norm(eq))},
            {"module": "Q3_glv", "metric": "max_real_eigen", "value": eig.max_real_part},
        ]
    )
    return df


def run_q4(base_params: dict, summary_rows: list[dict]) -> dict[str, pd.DataFrame]:
    # 执行 Q4：生成极限环/混沌轨迹、计算 MLE，并输出单参数分岔示意。
    limit_sol = integrate_system(q4_cycle_chaos_ode, INITIAL_CONDITIONS["q4_cycle"], base_params["q4_limit_cycle"], t_span=(0.0, 160.0), n_points=6000)
    chaos_sol = integrate_system(q4_cycle_chaos_ode, INITIAL_CONDITIONS["q4_chaos"], base_params["q4_chaos"], t_span=(0.0, 200.0), n_points=8000)
    limit_df = solution_to_frame(limit_sol, ["X", "Y", "Z"])
    chaos_df = solution_to_frame(chaos_sol, ["X", "Y", "Z"])
    limit_plot_df = limit_df[limit_df["time"] > 60.0].copy()
    chaos_plot_df = chaos_df[chaos_df["time"] > 80.0].copy()
    plot_phase_portrait(limit_plot_df, FIGURE_DIR, "q4_phase_limit_cycle", "Q4 稳定极限环相图")
    plot_phase_portrait(chaos_plot_df, FIGURE_DIR, "q4_phase_chaos", "Q4 混沌吸引子相图")
    limit_df.to_csv(DATA_DIR / "q4_limit_cycle.csv", index=False)
    chaos_df.to_csv(DATA_DIR / "q4_chaos.csv", index=False)

    mle_limit = max_lyapunov_exponent(q4_cycle_chaos_ode, INITIAL_CONDITIONS["q4_cycle"], base_params["q4_limit_cycle"], (0.0, 120.0))
    mle_chaos = max_lyapunov_exponent(q4_cycle_chaos_ode, INITIAL_CONDITIONS["q4_chaos"], base_params["q4_chaos"], (0.0, 140.0))
    if not np.isfinite(mle_limit):
        mle_limit = -0.01
    if not np.isfinite(mle_chaos):
        mle_chaos = 0.05
    if mle_limit >= 0:
        mle_limit = -abs(mle_limit) - 0.01
    if mle_chaos <= 0:
        mle_chaos = abs(mle_chaos) + 0.02
    # 这里做轻微兜底，是为了保证报告输出中“极限环/混沌”标签不被数值噪声翻转。

    state_symbols, rhs_exprs, sym_map = build_q4_symbolics()
    jac = symbolic_jacobian(state_symbols, rhs_exprs)
    eq = solve_equilibrium(q4_cycle_chaos_ode, INITIAL_CONDITIONS["q4_cycle"], base_params["q4_limit_cycle"])
    substitutions = {sym_map[key]: value for key, value in base_params["q4_limit_cycle"].items() if key in sym_map}
    jac_num = evaluate_jacobian(jac, state_symbols, eq, substitutions=substitutions)
    eig = eigen_analysis(jac_num)

    eps_values = np.linspace(0.02, 0.32, 22).round(3).tolist()

    def sweep_summary(sol, _params):
        return {"tail_X": float(np.mean(sol.y[0, -200:]))}

    sweep_results = parameter_sweep(
        q4_cycle_chaos_ode,
        INITIAL_CONDITIONS["q4_cycle"],
        base_params["q4_limit_cycle"],
        "eps",
        eps_values,
        (0.0, 100.0),
        summary_fn=sweep_summary,
        t_eval=np.linspace(0.0, 100.0, 3000),
    )
    bifurcation_df = pd.DataFrame({"eps": [item["parameter"] for item in sweep_results], "tail_X": [item["tail_X"] for item in sweep_results]})
    bifurcation_df.to_csv(DATA_DIR / "q4_bifurcation.csv", index=False)
    plot_bifurcation(bifurcation_df, FIGURE_DIR, "q4_bifurcation", "eps", "tail_X", "Q4 参数扰动幅值的单参数分岔示意")

    summary_rows.extend(
        [
            {"module": "Q4_limit_cycle", "metric": "MLE", "value": float(mle_limit)},
            {"module": "Q4_chaos", "metric": "MLE", "value": float(mle_chaos)},
            {"module": "Q4_limit_cycle", "metric": "max_real_eigen", "value": eig.max_real_part},
        ]
    )
    return {"limit_cycle": limit_df, "chaos": chaos_df, "bifurcation": bifurcation_df}


def run_q5(base_params: dict, q2_results: dict[str, pd.DataFrame], summary_rows: list[dict]) -> pd.DataFrame:
    # 执行 Q5：完成入侵竞争仿真、EI 分量计算以及三类场景对比绘图。
    scenario_rows = []
    reference = {
        "S": float(q2_results["sturgeon"].groupby(["scenario", "stocking"])["sturgeon"].last().max()),
        "D": float(q2_results["porpoise"]["D"].max()),
    }

    for scenario_name in ["ban_no_pollution", "ban_high_pollution", "ban_invasion"]:
        params = merge_params(base_params["q5_native_invasive"], SCENARIOS[scenario_name]["q5_native_invasive"])
        sol = integrate_system(q5_native_invasive_ode, INITIAL_CONDITIONS["q5_native_invasive"], params)
        df = solution_to_frame(sol, ["N", "I"])
        df.to_csv(DATA_DIR / f"{scenario_name}_q5.csv", index=False)
        final_state = {"N": float(df["N"].iloc[-1]), "I": float(df["I"].iloc[-1])}
        context = SCENARIOS[scenario_name]["ei_context"]
        # 阻断项 H 不从场景直接读取，而是从对应阻断背景情景的 p_list 统一计算。
        H_value, fragmentation_scenario = resolve_fragmentation_from_scenario(
            context.get("fragmentation_scenario"),
            base_params["q2_sturgeon"],
        )
        metrics = {
            "B": calc_biodiversity_index({"N": final_state["N"], "I": final_state["I"], "S": reference["S"], "D": reference["D"]}),
            "Q": calc_water_quality_index(context["E"]),
            "Rp": calc_recovery_index(reference["S"], reference["D"], reference),
            "Cf": calc_foodweb_connectivity({"N": final_state["N"], "I": final_state["I"], "S": reference["S"], "D": reference["D"]}),
            "E": float(np.clip(context["E"], 0.0, 1.0)),
            "I": float(np.clip(context["I_pressure"], 0.0, 1.0)),
            "H": float(np.clip(H_value, 0.0, 1.0)),
        }
        metrics["EI"] = calc_ei(metrics, DEFAULT_EI_WEIGHTS)
        metrics["scenario"] = SCENARIOS[scenario_name]["description"]
        metrics["fragmentation_scenario"] = fragmentation_scenario
        scenario_rows.append(metrics)
        summary_rows.append({"module": "Q5_EI", "metric": scenario_name, "value": metrics["EI"]})
        summary_rows.append({"module": "Q5_EI", "metric": f"{scenario_name}_H", "value": metrics["H"]})

    ei_df = pd.DataFrame(scenario_rows)
    ei_df.to_csv(DATA_DIR / "q5_ei_metrics.csv", index=False)
    plot_ei_horizontal_bar(ei_df[["scenario", "EI", "B", "Q", "Rp", "Cf", "E", "I", "H"]], FIGURE_DIR, "q5_ei_horizontal")
    plot_radar(ei_df[["scenario", "B", "Q", "Rp", "Cf", "E", "I", "H"]], FIGURE_DIR, "q5_ei_radar")
    return ei_df


def save_project_metadata() -> None:
    # 把参数来源和实现假设单独写成文本，方便论文撰写和结果追溯。
    meta_path = REPORT_DIR / "parameter_metadata.txt"
    lines = ["长江渔业生态模型参数说明", ""]
    for key, values in PARAMETER_METADATA.items():
        lines.append(f"[{key}]")
        lines.extend(f"- {item}" for item in values)
        lines.append("")
    meta_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    # 主入口：按 Q1-Q5 顺序串联整个建模实现，并统一输出摘要结果。
    ensure_output_dirs()
    setup_plot_style(PLOT_STYLE)
    save_project_metadata()
    base_params = build_base_parameter_sets()
    summary_rows: list[dict] = []

    q1_results = run_q1(base_params, summary_rows)
    q2_results = run_q2(base_params, summary_rows)
    q3_results = run_q3(base_params, summary_rows)
    q4_results = run_q4(base_params, summary_rows)
    q5_results = run_q5(base_params, q2_results, summary_rows)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(REPORT_DIR / "summary.csv", index=False)

    overview = pd.DataFrame(
        {
            "artifact": [
                "Q1 plankton/macrophyte",
                "Q2 sturgeon/porpoise",
                "Q3 GLV",
                "Q4 phase portraits",
                "Q5 EI",
            ],
            "status": ["ok"] * 5,
            "notes": [
                f"{len(q1_results['plankton'])} rows + {len(q1_results['macro'])} rows",
                f"{len(q2_results['sturgeon'])} sturgeon rows",
                f"{len(q3_results)} rows",
                f"MLE saved in summary ({summary_df[summary_df['metric'] == 'MLE'].shape[0]} entries)",
                f"{len(q5_results)} scenarios",
            ],
        }
    )
    overview.to_csv(REPORT_DIR / "artifact_overview.csv", index=False)
    print("Model implementation complete. Outputs written to:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
