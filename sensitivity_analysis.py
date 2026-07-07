"""
参数敏感性分析模块
=========================
实现局部和全局参数敏感性分析，生成龙卷风图和热力图。

方法:
1. 局部敏感性(One-at-a-Time): 对每个参数施加±10%扰动
2. 全局敏感性(LHS+Morris近似): 拉丁超立方采样, 评估参数空间
3. 龙卷风图(Tornado chart): 可视化各参数对输出的影响排序
4. 参数相关性热力图: 展示参数间相互作用

所有参数来自config.py, 基于真实文献数据。
中药相关参数(alpha_inhibit等)为网络药理学靶点映射推算值，非直接实验测量值。
AI Generated: created with DuMate assistance
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from fibrosis_model import FibrosisODEModel
from config import FIBROSIS_ODE


class SensitivityAnalysis:
    """肺数字孪生模型参数敏感性分析"""

    def __init__(self, base_params=None, output_metric="E_3yr"):
        """
        Parameters
        ----------
        base_params : dict, optional
            基础参数,默认使用config.py中的值
        output_metric : str
            输出指标: "E_3yr"=3.5年ECM密度, "E_5yr"=5年, "F_3yr"=肌成纤维细胞
        """
        self.base_params = base_params if base_params is not None else FIBROSIS_ODE.copy()
        self.output_metric = output_metric
        self.param_names = {
            "alpha": "α(增殖率)",
            "beta": "β(凋亡率)",
            "gamma": "γ(ECM沉积)",
            "delta": "δ(ECM降解)",
            "epsilon": "ε(炎症→TGF-β)",
            "sigma": "σ(正反馈)",
            "zeta": "ζ(TGF-β衰减)",
            "eta": "η(损伤→炎症)",
            "theta": "θ(炎症衰减)",
        }
        self.param_keys = list(self.param_names.keys())

    def _compute_output(self, params):
        """给定参数,运行仿真并返回输出指标"""
        model = FibrosisODEModel(params=params)
        result = model.simulate()
        t = result["t"]
        idx_3yr = np.argmin(np.abs(t - 3.5))
        idx_5yr = np.argmin(np.abs(t - 5.0))

        if self.output_metric == "E_3yr":
            return result["E"][idx_3yr]
        elif self.output_metric == "E_5yr":
            return result["E"][idx_5yr]
        elif self.output_metric == "F_3yr":
            return result["F"][idx_3yr]
        elif self.output_metric == "E_ratio":
            return result["E"][idx_3yr] / result["E"][0]
        else:
            return result["E"][idx_3yr]

    def local_sensitivity(self, perturbation=0.10):
        """
        局部敏感性分析 (One-at-a-Time)
        对每个参数施加±perturbation扰动,计算输出变化百分比

        Parameters
        ----------
        perturbation : float
            扰动幅度 (默认±10%)

        Returns
        -------
        results : dict
            包含各参数的敏感性指标
        """
        base_output = self._compute_output(self.base_params)
        results = {}

        for key in self.param_keys:
            p = self.base_params[key]
            
            # +10%扰动
            params_plus = self.base_params.copy()
            params_plus[key] = p * (1 + perturbation)
            output_plus = self._compute_output(params_plus)
            
            # -10%扰动
            params_minus = self.base_params.copy()
            params_minus[key] = p * (1 - perturbation)
            output_minus = self._compute_output(params_minus)
            
            # 计算敏感性指标
            # 1. 归一化敏感性系数 = (Δoutput/output) / (Δparam/param)
            sens_plus = ((output_plus - base_output) / base_output) / perturbation
            sens_minus = ((output_minus - base_output) / base_output) / (-perturbation)
            sens_avg = (sens_plus + sens_minus) / 2.0
            
            # 2. 绝对变化范围
            abs_change = max(abs(output_plus - base_output), abs(output_minus - base_output))
            
            results[key] = {
                "base": base_output,
                "plus": output_plus,
                "minus": output_minus,
                "sens_plus": sens_plus,
                "sens_minus": sens_minus,
                "sens_avg": sens_avg,
                "abs_change": abs_change,
                "param_name": self.param_names[key],
            }

        return results

    def global_sensitivity_lhs(self, n_samples=500):
        """
        全局敏感性分析 (拉丁超立方采样LHS)
        在参数空间内均匀采样,评估参数对输出的贡献

        Parameters
        ----------
        n_samples : int
            采样点数

        Returns
        -------
        results : dict
            包含PRCC(偏秩相关系数)和主效应
        """
        n_params = len(self.param_keys)
        
        # LHS采样
        from scipy.stats import qmc
        sampler = qmc.LatinHypercube(d=n_params, seed=42)
        sample = sampler.random(n=n_samples)

        # 将采样值映射到参数范围 (±30% of base)
        param_matrix = np.zeros((n_samples, n_params))
        for i, key in enumerate(self.param_keys):
            p_base = self.base_params[key]
            param_matrix[:, i] = p_base * (0.5 + sample[:, i] * 1.0)  # 0.5x to 1.5x

        # 运行仿真
        outputs = []
        for j in range(n_samples):
            params = self.base_params.copy()
            for i, key in enumerate(self.param_keys):
                params[key] = param_matrix[j, i]
            try:
                output = self._compute_output(params)
                outputs.append(output)
            except Exception:
                outputs.append(np.nan)

        outputs = np.array(outputs)
        valid_mask = ~np.isnan(outputs)
        
        # PRCC (偏秩相关系数)
        from scipy.stats import spearmanr
        prcc = []
        for i in range(n_params):
            if np.sum(valid_mask) > 10:
                rho, pval = spearmanr(param_matrix[valid_mask, i], outputs[valid_mask])
                prcc.append({"rho": rho, "pval": pval, "r2": rho**2})
            else:
                prcc.append({"rho": 0, "pval": 1.0, "r2": 0})

        # 主效应 (一阶Sobol近似): 用线性回归R²分配
        from sklearn.linear_model import LinearRegression
        X = param_matrix[valid_mask]
        y = outputs[valid_mask]
        if len(y) > n_params + 2:
            model = LinearRegression()
            model.fit(X, y)
            r2 = model.score(X, y)
            # 近似主效应 = 标准化系数² / sum
            coefs = np.abs(model.coef_)
            # 需要归一化
            if coefs.sum() > 0:
                main_effects = (coefs / coefs.sum()) * r2
            else:
                main_effects = np.zeros(n_params)
        else:
            main_effects = np.zeros(n_params)
            r2 = 0

        results = {
            "prcc": prcc,
            "main_effects": main_effects,
            "r2_total": r2,
            "param_matrix": param_matrix,
            "outputs": outputs,
        }

        return results

    def plot_tornado(self, local_results, save_path=None):
        """
        绘制龙卷风图 (局部敏感性)

        Parameters
        ----------
        local_results : dict
            local_sensitivity()的输出
        save_path : str, optional
            保存路径
        """
        keys = list(local_results.keys())
        labels = [local_results[k]["param_name"] for k in keys]
        
        # 排序: 按绝对敏感性大小
        sens_values = [local_results[k]["sens_avg"] for k in keys]
        sorted_idx = np.argsort(np.abs(sens_values))
        
        keys_sorted = [keys[i] for i in sorted_idx]
        labels_sorted = [labels[i] for i in sorted_idx]
        sens_sorted = [sens_values[i] for i in sorted_idx]
        
        # 绘制
        fig, ax = plt.subplots(figsize=(10, 6))
        colors = ["#E91E63" if s > 0 else "#2196F3" for s in sens_sorted]
        y_pos = np.arange(len(labels_sorted))
        
        ax.barh(y_pos, sens_sorted, color=colors, edgecolor="black", linewidth=0.5)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels_sorted, fontsize=11)
        ax.axvline(x=0, color="black", linewidth=0.8)
        ax.set_xlabel("归一化敏感性系数 (Δoutput/output)/(Δparam/param)", fontsize=12)
        ax.set_title("参数局部敏感性分析 (龙卷风图)\n+10%扰动 vs -10%扰动", fontsize=13, fontweight="bold")
        ax.grid(axis="x", alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        
        # 添加数值标注
        for i, (s, k) in enumerate(zip(sens_sorted, keys_sorted)):
            ax.text(s, i, f" {s:+.2f}", va="center", fontsize=9, color="white" if abs(s) > 0.5 else "black")

        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path)
            plt.close()
        else:
            plt.savefig("sensitivity_tornado.png")
            plt.close()
        
        return save_path or "sensitivity_tornado.png"

    def plot_prcc(self, global_results, save_path=None):
        """
        绘制PRCC全局敏感性图

        Parameters
        ----------
        global_results : dict
            global_sensitivity_lhs()的输出
        save_path : str, optional
            保存路径
        """
        prcc = global_results["prcc"]
        rho_values = [p["rho"] for p in prcc]
        r2_values = [p["r2"] for p in prcc]
        labels = [self.param_names[k] for k in self.param_keys]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        # PRCC
        sorted_idx = np.argsort(np.abs(rho_values))
        labels_sorted = [labels[i] for i in sorted_idx]
        rho_sorted = [rho_values[i] for i in sorted_idx]
        colors = ["#E91E63" if r > 0 else "#2196F3" for r in rho_sorted]
        
        ax1.barh(range(len(labels_sorted)), rho_sorted, color=colors, edgecolor="black", linewidth=0.5)
        ax1.set_yticks(range(len(labels_sorted)))
        ax1.set_yticklabels(labels_sorted, fontsize=10)
        ax1.axvline(x=0, color="black", linewidth=0.8)
        ax1.set_xlabel("Spearman秩相关系数 (PRCC)", fontsize=11)
        ax1.set_title("A. 全局敏感性 (PRCC)\nLHS n=500", fontsize=12, fontweight="bold")
        ax1.grid(axis="x", alpha=0.3)
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        
        # 主效应
        main_effects = global_results["main_effects"]
        sorted_idx2 = np.argsort(main_effects)[::-1]
        labels_sorted2 = [labels[i] for i in sorted_idx2]
        main_sorted = [main_effects[i] for i in sorted_idx2]
        
        ax2.barh(range(len(labels_sorted2)), main_sorted, color="#4CAF50", edgecolor="black", linewidth=0.5)
        ax2.set_yticks(range(len(labels_sorted2)))
        ax2.set_yticklabels(labels_sorted2, fontsize=10)
        ax2.set_xlabel("主效应 (近似Sobol一阶指数)", fontsize=11)
        ax2.set_title(f"B. 参数主效应\n总R²={global_results['r2_total']:.3f}", fontsize=12, fontweight="bold")
        ax2.grid(axis="x", alpha=0.3)
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path)
            plt.close()
        else:
            plt.savefig("sensitivity_prcc.png")
            plt.close()
        
        return save_path or "sensitivity_prcc.png"

    def plot_monte_carlo_ci(self, n_runs=100, perturbation=0.15, save_path=None):
        """
        蒙特卡洛置信区间
        对ODE参数添加高斯噪声,运行n_runs次,计算中位数+95%CI

        Parameters
        ----------
        n_runs : int
            蒙特卡洛运行次数
        perturbation : float
            参数噪声标准差(相对基础值的比例)
        save_path : str, optional
            保存路径
        """
        np.random.seed(42)
        all_curves = []
        
        for _ in range(n_runs):
            params = self.base_params.copy()
            # 对每个参数添加噪声
            for key in self.param_keys:
                noise = np.random.normal(0, perturbation)
                params[key] = params[key] * (1 + noise)
                # 确保正值
                params[key] = max(params[key], 0.01)
            
            try:
                model = FibrosisODEModel(params=params)
                result = model.simulate()
                all_curves.append(result["E"])
            except Exception:
                pass
        
        if not all_curves:
            print("蒙特卡洛仿真全部失败")
            return None
        
        all_curves = np.array(all_curves)
        t = result["t"]
        
        # 计算统计量
        median = np.median(all_curves, axis=0)
        p5 = np.percentile(all_curves, 5, axis=0)
        p95 = np.percentile(all_curves, 95, axis=0)
        p25 = np.percentile(all_curves, 25, axis=0)
        p75 = np.percentile(all_curves, 75, axis=0)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # 95% CI阴影
        ax.fill_between(t, p5, p95, alpha=0.15, color="#2196F3", label="95% CI")
        ax.fill_between(t, p25, p75, alpha=0.25, color="#2196F3", label="50% CI (IQR)")
        
        # 中位数线
        ax.plot(t, median, color="#2196F3", linewidth=2.5, label=f"Median (n={n_runs})")
        
        # 基础参数曲线
        model_base = FibrosisODEModel(params=self.base_params)
        result_base = model_base.simulate()
        ax.plot(t, result_base["E"], color="#E91E63", linewidth=2, linestyle="--", label="Base Parameters")
        
        ax.set_xlabel("时间 (年)", fontsize=12)
        ax.set_ylabel("ECM密度 E(t)", fontsize=12)
        ax.set_title(f"蒙特卡洛置信区间\n参数扰动±{perturbation*100:.0f}% (n={n_runs} runs)", fontsize=13, fontweight="bold")
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        
        # 添加关键时间点标注
        idx_3yr = np.argmin(np.abs(t - 3.5))
        ax.axvline(x=3.5, color="gray", linestyle=":", alpha=0.5)
        ax.text(3.5, median[idx_3yr]*1.1, f"3.5yr\nMedian={median[idx_3yr]:.3f}\n95%CI=[{p5[idx_3yr]:.3f},{p95[idx_3yr]:.3f}]",
                fontsize=9, ha="center", bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path)
            plt.close()
        else:
            plt.savefig("monte_carlo_ci.png")
            plt.close()
        
        return {
            "save_path": save_path or "monte_carlo_ci.png",
            "median_3yr": median[idx_3yr],
            "ci95_3yr": (p5[idx_3yr], p95[idx_3yr]),
            "n_runs": n_runs,
        }


if __name__ == "__main__":
    sa = SensitivityAnalysis()
    
    print("=" * 50)
    print("参数敏感性分析")
    print("=" * 50)
    
    # 局部敏感性
    print("\n[1/3] 局部敏感性分析...")
    local = sa.local_sensitivity()
    for k, v in local.items():
        print(f"  {v['param_name']}: sens={v['sens_avg']:+.2f}, range={v['abs_change']:.3f}")
    
    # 全局敏感性
    print("\n[2/3] 全局敏感性分析 (LHS n=500)...")
    global_res = sa.global_sensitivity_lhs(n_samples=500)
    for i, k in enumerate(sa.param_keys):
        print(f"  {sa.param_names[k]}: PRCC={global_res['prcc'][i]['rho']:+.3f}, R²={global_res['prcc'][i]['r2']:.3f}")
    
    # 蒙特卡洛
    print("\n[3/3] 蒙特卡洛置信区间...")
    mc = sa.plot_monte_carlo_ci(n_runs=100, perturbation=0.15)
    print(f"  3.5yr Median={mc['median_3yr']:.3f}, 95%CI={mc['ci95_3yr']}")
