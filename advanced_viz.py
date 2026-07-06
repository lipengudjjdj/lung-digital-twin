"""
高级可视化模块
====================
包含：
1. 蒙特卡洛置信区间图 (所有药物组误差棒)
2. 肺呼吸周期2D动画 (不同纤维化阶段)
3. FVC-ECM关系散点图
4. 参数热力图 (参数-参数-输出)

AI Generated: created with DuMate assistance
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch
from matplotlib import animation
from fibrosis_model import FibrosisODEModel
from respiratory_model import RespiratoryModel
from drug_intervention import DrugIntervention
from fvc_simulator import FVCSimulator
from config import FIBROSIS_ODE, DRUG_INTERVENTION


class AdvancedVisualization:
    """高级可视化"""

    def __init__(self):
        self.drug_module = DrugIntervention()
        self.fvc_sim = FVCSimulator()
        # Nature配色
        self.colors = {
            "base": "#37474F",
            "normal": "#42A5F5",
            "ipf": "#E91E63",
            "nintedanib": "#9C27B0",
            "pirfenidone": "#FF9800",
            "huangqi": "#4CAF50",
            "danshen": "#00BCD4",
            "gusuibu": "#795548",
            "baizhu": "#8BC34A",
            "danggui": "#F44336",
            "chuanxiong": "#3F51B5",
            "combo": "#FF5722",
        }

    def plot_drug_comparison_with_ci(self, n_mc=50, perturbation=0.12, save_path=None):
        """
        绘制药物对比+蒙特卡洛置信区间(误差棒)

        Parameters
        ----------
        n_mc : int
            蒙特卡洛运行次数
        perturbation : float
            参数扰动幅度
        save_path : str, optional
            保存路径
        """
        drug_ids = ["nintedanib", "pirfenidone", "huangqi", "danshen", "gusuibu", "baizhu", "danggui", "chuanxiong"]
        
        fig, axes = plt.subplots(1, 2, figsize=(16, 7))
        
        # 左图: ECM 3.5yr密度(误差棒)
        ax1 = axes[0]
        base_model = FibrosisODEModel()
        base_result = base_model.simulate()
        t = base_result["t"]
        idx_35 = np.argmin(np.abs(t - 3.5))
        E_base_35 = base_result["E"][idx_35]
        
        # 蒙特卡洛计算各药物组
        drug_means = []
        drug_stds = []
        drug_labels = ["安慰剂"]
        drug_colors_list = [self.colors["ipf"]]
        
        # 安慰剂MC
        # 只对ODE数值参数加噪声, 排除t_span等非数值字段
        ode_numeric_keys = ["alpha", "beta", "gamma", "delta", "epsilon",
                            "sigma", "zeta", "eta", "theta", "F0", "E0", "T0", "I0"]
        placebo_E35 = []
        for _ in range(n_mc):
            params = FIBROSIS_ODE.copy()
            for key in ode_numeric_keys:
                noise = np.random.normal(0, perturbation)
                params[key] = max(params[key] * (1 + noise), 0.01)
            try:
                m = FibrosisODEModel(params=params)
                r = m.simulate()
                placebo_E35.append(r["E"][idx_35])
            except Exception:
                pass
        
        drug_means.append(np.mean(placebo_E35))
        drug_stds.append(np.std(placebo_E35))
        
        for did in drug_ids:
            drug_info = DRUG_INTERVENTION.get(did, {})
            label = drug_info.get("name_cn", did)
            drug_labels.append(label)
            drug_colors_list.append(self.colors.get(did, "#607D8B"))
            
            mc_E35 = []
            for _ in range(n_mc):
                params = self.drug_module.apply_drug_to_ode_params(did)
                for key in ode_numeric_keys:
                    if key in params:
                        noise = np.random.normal(0, perturbation)
                        params[key] = max(params[key] * (1 + noise), 0.01)
                try:
                    D = params.get("D_intervention", 0)
                    m = FibrosisODEModel(params=params)
                    m.set_drug_intervention(D)
                    r = m.simulate()
                    mc_E35.append(r["E"][idx_35])
                except Exception:
                    pass
            
            drug_means.append(np.mean(mc_E35) if mc_E35 else 0)
            drug_stds.append(np.std(mc_E35) if len(mc_E35) > 1 else 0)
        
        # 联合用药
        combo_mc = []
        for _ in range(n_mc):
            params = self.drug_module.apply_combination(["huangqi", "danshen"])
            for key in ode_numeric_keys:
                if key in params:
                    noise = np.random.normal(0, perturbation)
                    params[key] = max(params[key] * (1 + noise), 0.01)
            try:
                D = params.get("D_intervention", 0)
                m = FibrosisODEModel(params=params)
                m.set_drug_intervention(D)
                r = m.simulate()
                combo_mc.append(r["E"][idx_35])
            except Exception:
                pass
        
        drug_labels.append("黄芪+丹参")
        drug_colors_list.append(self.colors["combo"])
        drug_means.append(np.mean(combo_mc) if combo_mc else 0)
        drug_stds.append(np.std(combo_mc) if len(combo_mc) > 1 else 0)
        
        x_pos = np.arange(len(drug_labels))
        bars = ax1.bar(x_pos, drug_means, yerr=[1.96*s for s in drug_stds],
                       color=drug_colors_list, edgecolor="black", linewidth=0.8,
                       capsize=5, error_kw={"linewidth": 1.5})
        
        ax1.set_xticks(x_pos)
        ax1.set_xticklabels(drug_labels, rotation=45, ha="right", fontsize=10)
        ax1.set_ylabel("ECM密度 E(3.5yr)", fontsize=12)
        ax1.set_title(f"A. 各干预组ECM密度 (MC n={n_mc}, 95%CI)", fontsize=13, fontweight="bold")
        ax1.axhline(y=drug_means[0], color="gray", linestyle="--", alpha=0.5)
        ax1.grid(axis="y", alpha=0.3)
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        
        # 右图: FVC下降率对比(误差棒)
        ax2 = axes[1]
        fvc_means = []
        fvc_stds = []
        
        # 安慰剂FVC
        placebo_fvc = []
        for _ in range(n_mc):
            params = FIBROSIS_ODE.copy()
            for key in ode_numeric_keys:
                noise = np.random.normal(0, perturbation)
                params[key] = max(params[key] * (1 + noise), 0.01)
            try:
                m = FibrosisODEModel(params=params)
                r = m.simulate(t_span=(0, 1), t_eval_points=100)
                E_prog = r["E"][-1] - r["E"][0]
                FVC_decline = self.fvc_sim.k_FVC_ECM * E_prog
                placebo_fvc.append(FVC_decline)
            except Exception:
                pass
        
        fvc_means.append(np.mean(placebo_fvc))
        fvc_stds.append(np.std(placebo_fvc))
        
        for did in drug_ids:
            mc_fvc = []
            for _ in range(n_mc):
                params = self.drug_module.apply_drug_to_ode_params(did)
                for key in ode_numeric_keys:
                    if key in params:
                        noise = np.random.normal(0, perturbation)
                        params[key] = max(params[key] * (1 + noise), 0.01)
                try:
                    D = params.get("D_intervention", 0)
                    m = FibrosisODEModel(params=params)
                    m.set_drug_intervention(D)
                    r = m.simulate(t_span=(0, 1), t_eval_points=100)
                    E_prog = r["E"][-1] - r["E"][0]
                    FVC_decline = self.fvc_sim.k_FVC_ECM * E_prog
                    mc_fvc.append(FVC_decline)
                except Exception:
                    pass
            
            fvc_means.append(np.mean(mc_fvc) if mc_fvc else 0)
            fvc_stds.append(np.std(mc_fvc) if len(mc_fvc) > 1 else 0)
        
        # 联合FVC
        combo_fvc = []
        for _ in range(n_mc):
            params = self.drug_module.apply_combination(["huangqi", "danshen"])
            for key in ode_numeric_keys:
                if key in params:
                    noise = np.random.normal(0, perturbation)
                    params[key] = max(params[key] * (1 + noise), 0.01)
            try:
                D = params.get("D_intervention", 0)
                m = FibrosisODEModel(params=params)
                m.set_drug_intervention(D)
                r = m.simulate(t_span=(0, 1), t_eval_points=100)
                E_prog = r["E"][-1] - r["E"][0]
                FVC_decline = self.fvc_sim.k_FVC_ECM * E_prog
                combo_fvc.append(FVC_decline)
            except Exception:
                pass
        
        fvc_means.append(np.mean(combo_fvc) if combo_fvc else 0)
        fvc_stds.append(np.std(combo_fvc) if len(combo_fvc) > 1 else 0)
        
        # 临床参考线
        ax2.bar(x_pos, fvc_means, yerr=[1.96*s for s in fvc_stds],
                color=drug_colors_list, edgecolor="black", linewidth=0.8,
                capsize=5, error_kw={"linewidth": 1.5})
        
        # INPULSIS参考
        ax2.axhline(y=239.9, color="#9E9E9E", linestyle="--", alpha=0.7, label="INPULSIS安慰剂")
        ax2.axhline(y=114.1, color="#9C27B0", linestyle="--", alpha=0.7, label="INPULSIS尼达尼布")
        ax2.axhline(y=131.2, color="#FF9800", linestyle="--", alpha=0.7, label="ASCEND吡非尼酮")
        
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(drug_labels, rotation=45, ha="right", fontsize=10)
        ax2.set_ylabel("FVC年下降率 (mL/yr)", fontsize=12)
        ax2.set_title(f"B. FVC年下降率 (MC n={n_mc}, 95%CI)", fontsize=13, fontweight="bold")
        ax2.legend(fontsize=8, loc="upper right")
        ax2.grid(axis="y", alpha=0.3)
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path)
            plt.close()
        else:
            plt.savefig("drug_comparison_ci.png")
            plt.close()
        
        return save_path or "drug_comparison_ci.png"

    def plot_fvc_ecm_relationship(self, save_path=None):
        """
        绘制FVC与ECM密度的非线性关系
        """
        model = FibrosisODEModel()
        result = model.simulate()
        t = result["t"]
        E = result["E"]
        
        # FVC映射
        k_FVC_ECM = 450.0
        FVC_baseline = 3500.0
        FVC_predicted = 4000.0
        E_progress = E - E[0]
        FVC = FVC_baseline - k_FVC_ECM * E_progress
        FVC = np.maximum(FVC, 1500.0)
        FVC_pct = (FVC / FVC_predicted) * 100
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        
        # 左图: FVC% vs ECM
        ax1.plot(E, FVC_pct, color="#2196F3", linewidth=2.5)
        ax1.axhline(y=80, color="red", linestyle="--", alpha=0.5, label="诊断阈值 80%")
        ax1.axhline(y=50, color="darkred", linestyle="--", alpha=0.5, label="重度 50%")
        ax1.fill_between(E, 80, 100, alpha=0.1, color="green", label="轻度")
        ax1.fill_between(E, 50, 80, alpha=0.1, color="orange", label="中度")
        ax1.fill_between(E, 0, 50, alpha=0.1, color="red", label="重度")
        ax1.set_xlabel("ECM密度 E(t)", fontsize=12)
        ax1.set_ylabel("FVC占预计值 (%)", fontsize=12)
        ax1.set_title("A. FVC与ECM密度关系", fontsize=13, fontweight="bold")
        ax1.legend(fontsize=9)
        ax1.grid(True, alpha=0.3)
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        
        # 右图: 时间轴上的FVC%和ECM双Y轴
        ax2_twin = ax2.twinx()
        l1 = ax2.plot(t, FVC_pct, color="#2196F3", linewidth=2.5, label="FVC%")
        l2 = ax2_twin.plot(t, E, color="#E91E63", linewidth=2.5, label="ECM密度")
        
        ax2.set_xlabel("时间 (年)", fontsize=12)
        ax2.set_ylabel("FVC占预计值 (%)", fontsize=12, color="#2196F3")
        ax2_twin.set_ylabel("ECM密度 E(t)", fontsize=12, color="#E91E63")
        ax2.axhline(y=80, color="red", linestyle="--", alpha=0.3)
        
        lines = l1 + l2
        labels = [l.get_label() for l in lines]
        ax2.legend(lines, labels, fontsize=10, loc="center right")
        ax2.set_title("B. FVC%与ECM时间序列", fontsize=13, fontweight="bold")
        ax2.grid(True, alpha=0.3)
        ax2.spines["top"].set_visible(False)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path)
            plt.close()
        else:
            plt.savefig("fvc_ecm_relationship.png")
            plt.close()
        
        return save_path or "fvc_ecm_relationship.png"

    def generate_breathing_animation(self, n_stages=4, save_path=None):
        """
        生成肺呼吸周期2D动画(不同纤维化阶段)

        Parameters
        ----------
        n_stages : int
            纤维化阶段数(0=正常, 1-3=轻中重度)
        save_path : str, optional
            保存路径(GIF)
        """
        fig, axes = plt.subplots(1, n_stages, figsize=(16, 4.5))
        if n_stages == 1:
            axes = [axes]
        
        stage_labels = ["正常肺", "轻度纤维化", "中度纤维化", "重度纤维化"]
        fibrosis_ratios = [1.0, 2.0, 3.5, 6.0]
        compliance_ratios = [1.0, 0.6, 0.3, 0.12]
        
        # 肺轮廓参数
        theta = np.linspace(0, 2 * np.pi, 100)
        
        # 动画帧
        n_frames = 60
        breath_cycle = np.sin(np.linspace(0, 2 * np.pi, n_frames))  # 呼吸周期
        
        def draw_lung(ax, stage_idx, frame):
            ax.clear()
            
            ratio = fibrosis_ratios[stage_idx]
            C_ratio = compliance_ratios[stage_idx]
            phase = breath_cycle[frame]
            
            # 肺体积变化(潮气量)
            V_base = 2.5  # L (FRC)
            V_tidal = 0.5 * C_ratio  # 潮气量随顺应性降低
            V_current = V_base + V_tidal * (phase + 1) / 2
            
            # 2D肺轮廓(椭圆)
            a = 1.0 + 0.15 * (phase + 1) / 2 * C_ratio  # 横轴
            b = 1.2 + 0.2 * (phase + 1) / 2 * C_ratio   # 纵轴
            
            # 正常肺轮廓
            x_normal = a * np.cos(theta)
            y_normal = b * np.sin(theta)
            
            # 纤维化肺(略小, 不规则)
            shrink = 1.0 / np.sqrt(ratio / fibrosis_ratios[0])
            irregularity = 0.03 * (ratio - 1) * np.sin(7 * theta + frame * 0.1)
            
            x_lung = a * shrink * np.cos(theta) * (1 + irregularity)
            y_lung = b * shrink * np.sin(theta) * (1 + irregularity)
            
            # 绘制
            color_intensity = min(0.3 + 0.15 * (ratio - 1), 0.85)
            lung_color = (0.9, 1 - color_intensity, 1 - color_intensity)
            
            ax.fill(x_lung, y_lung, color=lung_color, alpha=0.7, edgecolor="black", linewidth=2)
            
            # 气道
            ax.plot([0, 0], [b * shrink + 0.1, b * shrink + 0.5], color="black", linewidth=2)
            ax.plot([-0.3, 0, 0.3], [b * shrink + 0.5, b * shrink + 0.3, b * shrink + 0.5],
                    color="black", linewidth=2)
            
            # 纤维化斑点(越严重越多)
            n_spots = int((ratio - 1) * 8)
            np.random.seed(42 + stage_idx)
            for _ in range(n_spots):
                sx = np.random.uniform(-0.7, 0.7) * shrink
                sy = np.random.uniform(-0.8, 0.8) * shrink
                sr = np.random.uniform(0.02, 0.06)
                spot = plt.Circle((sx, sy), sr, color="#B71C1C", alpha=0.6)
                ax.add_patch(spot)
            
            # 标注
            C_0 = 0.2
            C_current = C_0 * C_ratio
            ax.set_title(f"{stage_labels[stage_idx]}\nC={C_current:.3f} L/cmH2O", fontsize=11, fontweight="bold")
            ax.set_xlim(-1.8, 1.8)
            ax.set_ylim(-1.8, 2.2)
            ax.set_aspect("equal")
            ax.axis("off")
        
        def animate(frame):
            for i in range(n_stages):
                draw_lung(axes[i], i, frame)
            fig.suptitle(f"肺呼吸周期 - 不同纤维化阶段 (帧 {frame+1}/{n_frames})",
                        fontsize=13, fontweight="bold")
            return axes
        
        anim = animation.FuncAnimation(fig, animate, frames=n_frames, interval=100, blit=False)
        
        if save_path:
            anim.save(save_path, writer="pillow", fps=10)
        else:
            anim.save("breathing_animation.gif", writer="pillow", fps=10)
        
        plt.close()
        return save_path or "breathing_animation.gif"

    def plot_parameter_heatmap(self, param1="gamma", param2="delta",
                                n_grid=20, save_path=None):
        """
        参数-参数-输出热力图
        固定其他参数,扫描param1和param2的二维空间,计算ECM 3.5yr输出

        Parameters
        ----------
        param1 : str
            第一个参数名(横轴)
        param2 : str
            第二个参数名(纵轴)
        n_grid : int
            网格分辨率
        save_path : str, optional
            保存路径
        """
        # 参数扫描范围
        p1_range = np.linspace(FIBROSIS_ODE[param1] * 0.3, FIBROSIS_ODE[param1] * 2.0, n_grid)
        p2_range = np.linspace(FIBROSIS_ODE[param2] * 0.3, FIBROSIS_ODE[param2] * 2.0, n_grid)
        
        output_matrix = np.zeros((n_grid, n_grid))
        
        for i, p2_val in enumerate(p2_range):
            for j, p1_val in enumerate(p1_range):
                params = FIBROSIS_ODE.copy()
                params[param1] = p1_val
                params[param2] = p2_val
                try:
                    model = FibrosisODEModel(params=params)
                    result = model.simulate()
                    t = result["t"]
                    idx_35 = np.argmin(np.abs(t - 3.5))
                    output_matrix[i, j] = result["E"][idx_35]
                except Exception:
                    output_matrix[i, j] = np.nan
        
        fig, ax = plt.subplots(figsize=(9, 7))
        
        im = ax.imshow(output_matrix, extent=[p1_range[0], p1_range[-1], p2_range[0], p2_range[-1]],
                       origin="lower", aspect="auto", cmap="YlOrRd")
        
        # 标注基础参数位置
        ax.scatter([FIBROSIS_ODE[param1]], [FIBROSIS_ODE[param2]], color="blue",
                   s=150, marker="*", zorder=5, edgecolors="black", linewidth=1.5,
                   label="校准参数")
        
        # 临床参考区域
        ax.contour(p1_range, p2_range, output_matrix, levels=[0.15, 0.20, 0.30, 0.50],
                   colors="black", linewidths=0.8, alpha=0.5)
        
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("ECM密度 E(3.5yr)", fontsize=12)
        
        ax.set_xlabel(f"{param1} ({param1}={FIBROSIS_ODE[param1]:.3f})", fontsize=12)
        ax.set_ylabel(f"{param2} ({param2}={FIBROSIS_ODE[param2]:.3f})", fontsize=12)
        ax.set_title(f"参数空间扫描: {param1} vs {param2}\nECM 3.5yr输出热力图", fontsize=13, fontweight="bold")
        ax.legend(fontsize=10)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path)
            plt.close()
        else:
            plt.savefig("parameter_heatmap.png")
            plt.close()
        
        return save_path or "parameter_heatmap.png"


if __name__ == "__main__":
    viz = AdvancedVisualization()
    
    print("[1/4] 药物对比+置信区间...")
    viz.plot_drug_comparison_with_ci(n_mc=30, perturbation=0.12)
    
    print("[2/4] FVC-ECM关系...")
    viz.plot_fvc_ecm_relationship()
    
    print("[3/4] 呼吸动画...")
    viz.generate_breathing_animation()
    
    print("[4/4] 参数热力图...")
    viz.plot_parameter_heatmap(param1="gamma", param2="delta", n_grid=15)
    
    print("所有高级可视化完成！")
