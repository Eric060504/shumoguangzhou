# 画图工具函数，封装了常用的绘图操作，使用Matplotlib和Seaborn库，支持时间序列、相位图、分岔图、EI指标对比等多种类型的图表生成，并提供统一的样式设置和保存功能
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def setup_plot_style(style_cfg: dict) -> None:
    sns.set_theme(style="whitegrid", context=style_cfg.get("context", "talk"), palette=style_cfg.get("palette", "Set2"))
    mpl.rcParams["figure.dpi"] = style_cfg.get("dpi", 180)
    mpl.rcParams["savefig.dpi"] = style_cfg.get("save_dpi", 300)
    mpl.rcParams["axes.unicode_minus"] = False
    mpl.rcParams["font.sans-serif"] = style_cfg.get("font_candidates", ["SimHei", "sans-serif"])


def save_figure(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_dir / f"{stem}.png", bbox_inches="tight")
    fig.savefig(out_dir / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def plot_time_series(df: pd.DataFrame, out_dir: Path, title: str, stem: str, ylabel: str = "生物量 / 相对数量") -> None:
    fig, ax = plt.subplots(figsize=(11, 6.5))
    value_cols = [col for col in df.columns if col != "time"]
    melted = df.melt(id_vars="time", value_vars=value_cols, var_name="变量", value_name="数值")
    sns.lineplot(data=melted, x="time", y="数值", hue="变量", linewidth=2.5, ax=ax)
    ax.set_title(title)
    ax.set_xlabel("时间 / 年")
    ax.set_ylabel(ylabel)
    ax.legend(title="")
    save_figure(fig, out_dir, stem)


def plot_sturgeon_recovery(df: pd.DataFrame, out_dir: Path, stem: str) -> None:
    fig, ax = plt.subplots(figsize=(11, 6.5))
    sns.lineplot(data=df, x="time", y="sturgeon", hue="scenario", style="stocking", linewidth=2.6, ax=ax)
    ax.set_title("鲟类在不同阻断率与放流情景下的恢复曲线")
    ax.set_xlabel("时间 / 年")
    ax.set_ylabel("鲟类生物量")
    ax.legend(title="")
    save_figure(fig, out_dir, stem)


def plot_phase_portrait(df: pd.DataFrame, out_dir: Path, stem: str, title: str) -> None:
    fig = plt.figure(figsize=(8, 6.8))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(df["X"], df["Y"], df["Z"], lw=1.5, color="#1f77b4")
    ax.set_title(title)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    save_figure(fig, out_dir, stem)


def plot_bifurcation(df: pd.DataFrame, out_dir: Path, stem: str, x_col: str, y_col: str, title: str) -> None:
    fig, ax = plt.subplots(figsize=(10, 6.3))
    sns.scatterplot(data=df, x=x_col, y=y_col, s=40, color="#c44e52", ax=ax)
    ax.set_title(title)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    save_figure(fig, out_dir, stem)


def plot_ei_horizontal_bar(df: pd.DataFrame, out_dir: Path, stem: str) -> None:
    plot_df = df.melt(id_vars="scenario", var_name="metric", value_name="value")
    fig, ax = plt.subplots(figsize=(12.5, 7.0))
    sns.barplot(data=plot_df, y="metric", x="value", hue="scenario", orient="h", ax=ax)
    ax.set_title("不同情境下 EI 及其指标分量对比")
    ax.set_xlabel("归一化指标值")
    ax.set_ylabel("")
    ax.legend(title="")
    save_figure(fig, out_dir, stem)


def plot_radar(df: pd.DataFrame, out_dir: Path, stem: str) -> None:
    metrics = [col for col in df.columns if col != "scenario"]
    angles = np.linspace(0, 2 * np.pi, len(metrics), endpoint=False).tolist()
    angles += angles[:1]
    fig, ax = plt.subplots(figsize=(8.5, 8.5), subplot_kw={"polar": True})
    for _, row in df.iterrows():
        values = row[metrics].tolist()
        values += values[:1]
        ax.plot(angles, values, linewidth=2.2, label=row["scenario"])
        ax.fill(angles, values, alpha=0.12)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(metrics)
    ax.set_title("EI 指标雷达图")
    ax.legend(loc="upper right", bbox_to_anchor=(1.25, 1.1))
    save_figure(fig, out_dir, stem)


def plot_sensitivity_curve(
    df: pd.DataFrame,
    out_dir: Path,
    stem: str,
    title: str,
    x_col: str = "value",
    y_col: str = "response_value",
    hue_col: str | None = None,
) -> None:
    # 用统一风格输出单因素灵敏度曲线，便于比较参数扫描对响应量的影响方向和强弱。
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    if hue_col and hue_col in df.columns:
        sns.lineplot(data=df, x=x_col, y=y_col, hue=hue_col, marker="o", linewidth=2.4, ax=ax)
    else:
        sns.lineplot(data=df, x=x_col, y=y_col, marker="o", linewidth=2.4, ax=ax, color="#4c72b0")
    ax.set_title(title)
    ax.set_xlabel("参数取值")
    ax.set_ylabel("响应量")
    if hue_col and hue_col in df.columns:
        ax.legend(title="")
    save_figure(fig, out_dir, stem)

