"""
可视化模块
============
生成学术期刊级别的图表，用于论文插图和交互式展示。

图表类型:
  1. 纤维化进程曲线 (F, E, T, I vs 时间)
  2. PV曲线对比 (正常 vs 不同纤维化程度)
  3. 呼吸周期模拟 (压力-容积-流量)
  4. 药物干预对比 (无干预 vs 各药物)
  5. 肺功能指标随纤维化变化
  6. 联合用药效果对比

输出格式: PNG (300dpi, 学术出版质量)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")  # 非交互式后端，避免GUI依赖
import matplotlib.pyplot as plt
from matplotlib import rcParams

# 中文字体配置 (Windows环境)
rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "DejaVu Sans"]
rcParams["axes.unicode_minus"] = False
rcParams["figure.dpi"] = 150
rcParams["savefig.dpi"] = 300
rcParams["savefig.bbox"] = "tight"

# 学术论文配色方案 (Nature风格)
COLORS = {
    "normal": "#2196F3",       # 蓝色 - 正常
    "mild": "#FF9800",         # 橙色 - 轻度
    "moderate": "#F44336",     # 红色 - 中度
    "severe": "#9C27B0",       # 紫色 - 重度
    "fibroblast": "#E91E63",   # 粉色 - 成纤维细胞
    "ecm": "#FF5722",          # 深橙 - ECM
    "tgfbeta": "#4CAF50",      # 绿色 - TGF-β
    "inflammation": "#00BCD4", # 青色 - 炎症
    "drug": "#3F51B5",         # 靛蓝 - 药物干预
    "combo": "#009688",        # 青绿 - 联合用药
}


def plot_fibrosis_progression(result, save_path=None, title=None):
    """
    绘制纤维化进程曲线

    Parameters
    ----------
    result : dict
        fibrosis_model.simulate()的输出
    save_path : str, optional
        保存路径
    title : str, optional
        图表标题
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    t = result["t"]

    # 肌成纤维细胞
    axes[0, 0].plot(t, result["F"], color=COLORS["fibroblast"], linewidth=2)
    axes[0, 0].set_ylabel("肌成纤维细胞密度 F(t)", fontsize=11)
    axes[0, 0].set_title("A. 肌成纤维细胞动力学", fontsize=12, fontweight="bold")
    axes[0, 0].set_ylim(-0.05, 1.05)
    axes[0, 0].axhline(y=0.5, color="gray", linestyle="--", alpha=0.5)
    axes[0, 0].fill_between(t, 0, result["F"], alpha=0.15, color=COLORS["fibroblast"])

    # ECM
    axes[0, 1].plot(t, result["E"], color=COLORS["ecm"], linewidth=2)
    axes[0, 1].set_ylabel("ECM/胶原密度 E(t)", fontsize=11)
    axes[0, 1].set_title("B. ECM沉积动力学", fontsize=12, fontweight="bold")
    axes[0, 1].set_ylim(-0.05, 1.05)
    axes[0, 1].axhline(y=0.5, color="gray", linestyle="--", alpha=0.5)
    axes[0, 1].fill_between(t, 0, result["E"], alpha=0.15, color=COLORS["ecm"])

    # TGF-β
    axes[1, 0].plot(t, result["T"], color=COLORS["tgfbeta"], linewidth=2)
    axes[1, 0].set_ylabel("TGF-β1浓度 T(t)", fontsize=11)
    axes[1, 0].set_title("C. TGF-β1动力学", fontsize=12, fontweight="bold")
    axes[1, 0].set_ylim(-0.05, 1.05)
    axes[1, 0].fill_between(t, 0, result["T"], alpha=0.15, color=COLORS["tgfbeta"])

    # 炎症
    axes[1, 1].plot(t, result["I"], color=COLORS["inflammation"], linewidth=2)
    axes[1, 1].set_ylabel("炎症因子水平 I(t)", fontsize=11)
    axes[1, 1].set_title("D. 炎症因子动力学", fontsize=12, fontweight="bold")
    axes[1, 1].set_ylim(-0.05, 1.05)
    axes[1, 1].fill_between(t, 0, result["I"], alpha=0.15, color=COLORS["inflammation"])

    for ax in axes.flat:
        ax.set_xlabel("时间 (年)", fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    if title:
        fig.suptitle(title, fontsize=14, fontweight="bold", y=1.02)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        plt.close()
    else:
        plt.savefig("fibrosis_progression.png")
        plt.close()

    return save_path or "fibrosis_progression.png"


def plot_pv_curves(respiratory_model, E_values=None, save_path=None):
    """
    绘制不同纤维化程度的PV曲线对比

    Parameters
    ----------
    respiratory_model : RespiratoryModel
        呼吸力学模型实例
    E_values : list of float, optional
        要对比的ECM密度值列表
    save_path : str, optional
        保存路径
    """
    if E_values is None:
        E_values = [0.05, 0.2, 0.4, 0.6, 0.8]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # 左图: 线性PV曲线
    for E in E_values:
        P, V = respiratory_model.pv_curve(E)
        C = respiratory_model.compute_compliance(E)
        label = f"E={E:.2f} (C={C:.3f} L/cmH2O)"
        color = _E_to_color(E)
        ax1.plot(V, P, color=color, linewidth=2, label=label)

    ax1.set_xlabel("肺容积 (L)", fontsize=12)
    ax1.set_ylabel("跨肺压 (cmH2O)", fontsize=12)
    ax1.set_title("A. 静态PV曲线", fontsize=13, fontweight="bold")
    ax1.legend(fontsize=9, loc="upper left")
    ax1.grid(True, alpha=0.3)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # 右图: Sigmoid PV曲线
    for E in E_values:
        P, V = respiratory_model.pv_curve_sigmoid(E)
        C = respiratory_model.compute_compliance(E)
        label = f"E={E:.2f} (C={C:.3f} L/cmH2O)"
        color = _E_to_color(E)
        ax2.plot(V, P, color=color, linewidth=2, label=label)

    ax2.set_xlabel("肺容积 (L)", fontsize=12)
    ax2.set_ylabel("跨肺压 (cmH2O)", fontsize=12)
    ax2.set_title("B. S形PV曲线 (更符合真实肺力学)", fontsize=13, fontweight="bold")
    ax2.legend(fontsize=9, loc="upper left")
    ax2.grid(True, alpha=0.3)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        plt.close()
    else:
        plt.savefig("pv_curves_comparison.png")
        plt.close()

    return save_path or "pv_curves_comparison.png"


def plot_drug_comparison(results_dict, save_path=None):
    """
    绘制药物干预对比图

    Parameters
    ----------
    results_dict : dict
        {标签: simulate_result} 字典
    save_path : str, optional
        保存路径
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    color_list = list(COLORS.values())[:len(results_dict)]

    for i, (label, result) in enumerate(results_dict.items()):
        color = color_list[i] if i < len(color_list) else "gray"
        ls = "--" if "干预" in label or label != "无干预" else "-"

        axes[0].plot(result["t"], result["F"], color=color, linewidth=2,
                     linestyle=ls, label=label)
        axes[1].plot(result["t"], result["E"], color=color, linewidth=2,
                     linestyle=ls, label=label)
        axes[2].plot(result["t"], result["T"], color=color, linewidth=2,
                     linestyle=ls, label=label)

    axes[0].set_title("A. 肌成纤维细胞密度", fontsize=12, fontweight="bold")
    axes[0].set_ylabel("F(t)", fontsize=11)

    axes[1].set_title("B. ECM密度", fontsize=12, fontweight="bold")
    axes[1].set_ylabel("E(t)", fontsize=11)

    axes[2].set_title("C. TGF-β1浓度", fontsize=12, fontweight="bold")
    axes[2].set_ylabel("T(t)", fontsize=11)

    for ax in axes:
        ax.set_xlabel("时间 (年)", fontsize=11)
        ax.set_ylim(-0.05, 1.05)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        plt.close()
    else:
        plt.savefig("drug_comparison.png")
        plt.close()

    return save_path or "drug_comparison.png"


def plot_lung_function_metrics(respiratory_model, save_path=None):
    """
    绘制肺功能指标随纤维化程度的变化

    Parameters
    ----------
    respiratory_model : RespiratoryModel
        呼吸力学模型实例
    save_path : str, optional
        保存路径
    """
    E_range = np.linspace(0.05, 0.9, 100)

    C_list = []
    FVC_list = []
    DLco_list = []
    W_list = []
    stiffness_list = []

    for E in E_range:
        m = respiratory_model.compute_lung_function_metrics(E)
        C_list.append(m["C_lung_mL_cmH2O"])
        FVC_list.append(m["FVC_percent_predicted"])
        DLco_list.append(m["DLco_percent_predicted"])
        W_list.append(m["P_for_500mL_cmH2O"])
        stiffness_list.append(m["stiffness_ratio"])

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 肺顺应性
    axes[0, 0].plot(E_range, C_list, color=COLORS["normal"], linewidth=2)
    axes[0, 0].axhline(y=200, color="gray", linestyle="--", alpha=0.5, label="正常值 (~200)")
    axes[0, 0].axhspan(50, 100, alpha=0.1, color="red", label="IPF范围 (50-100)")
    axes[0, 0].set_ylabel("肺顺应性 (mL/cmH2O)", fontsize=11)
    axes[0, 0].set_title("A. 肺顺应性 vs 纤维化", fontsize=12, fontweight="bold")
    axes[0, 0].legend(fontsize=9)

    # FVC
    axes[0, 1].plot(E_range, FVC_list, color=COLORS["ecm"], linewidth=2)
    axes[0, 1].axhline(y=80, color="red", linestyle="--", alpha=0.5, label="诊断阈值 (80%)")
    axes[0, 1].set_ylabel("FVC占预计值 (%)", fontsize=11)
    axes[0, 1].set_title("B. FVC vs 纤维化", fontsize=12, fontweight="bold")
    axes[0, 1].legend(fontsize=9)

    # DLco
    axes[1, 0].plot(E_range, DLco_list, color=COLORS["tgfbeta"], linewidth=2)
    axes[1, 0].axhline(y=80, color="red", linestyle="--", alpha=0.5, label="诊断阈值 (80%)")
    axes[1, 0].set_ylabel("DLco占预计值 (%)", fontsize=11)
    axes[1, 0].set_title("C. DLco vs 纤维化", fontsize=12, fontweight="bold")
    axes[1, 0].legend(fontsize=9)

    # 呼吸驱动压
    axes[1, 1].plot(E_range, W_list, color=COLORS["severe"], linewidth=2)
    axes[1, 1].set_ylabel("产生500mL潮气量所需压力 (cmH2O)", fontsize=11)
    axes[1, 1].set_title("D. 呼吸做功 vs 纤维化", fontsize=12, fontweight="bold")

    for ax in axes.flat:
        ax.set_xlabel("ECM密度 E", fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        plt.close()
    else:
        plt.savefig("lung_function_metrics.png")
        plt.close()

    return save_path or "lung_function_metrics.png"


def plot_breathing_cycle(respiratory_model, E_values=None, save_path=None):
    """
    绘制不同纤维化程度的呼吸周期对比

    Parameters
    ----------
    respiratory_model : RespiratoryModel
        呼吸力学模型实例
    E_values : list of float, optional
        ECM密度值列表
    save_path : str, optional
        保存路径
    """
    if E_values is None:
        E_values = [0.05, 0.4, 0.8]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for E in E_values:
        t, P, V, flow = respiratory_model.breathing_cycle(E)
        C = respiratory_model.compute_compliance(E)
        label = f"E={E:.2f} (C={C:.3f} L/cmH2O)"
        color = _E_to_color(E)

        axes[0].plot(t, V, color=color, linewidth=2, label=label)
        axes[1].plot(V, P, color=color, linewidth=2, label=label)

    axes[0].set_xlabel("时间 (秒)", fontsize=12)
    axes[0].set_ylabel("肺容积 (L)", fontsize=12)
    axes[0].set_title("A. 呼吸周期 — 容积变化", fontsize=13, fontweight="bold")

    axes[1].set_xlabel("肺容积 (L)", fontsize=12)
    axes[1].set_ylabel("跨肺压 (cmH2O)", fontsize=12)
    axes[1].set_title("B. 动态PV环", fontsize=13, fontweight="bold")

    for ax in axes:
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        plt.close()
    else:
        plt.savefig("breathing_cycle.png")
        plt.close()

    return save_path or "breathing_cycle.png"


def _E_to_color(E):
    """将ECM密度映射为颜色"""
    if E < 0.15:
        return COLORS["normal"]
    elif E < 0.35:
        return COLORS["mild"]
    elif E < 0.6:
        return COLORS["moderate"]
    else:
        return COLORS["severe"]
