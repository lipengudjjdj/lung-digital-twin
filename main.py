"""
肺数字孪生模型 — 主程序 (v3.0)
=================================
v3.0新增（基于最新文献）:
  1. YAP/TAZ机械转导恶性循环 (Advanced Science 2025)
  2. 成纤维细胞亚型动力学 Csmd1+/Cd248+ (Nature空间转录组学)
  3. TGF-β梯度双角色 (Science Advances 2026, 左为团队)
  4. 药代动力学PK-PD联合建模 (Eur J Clin Pharmacol 2018等)
  5. Rentosertib (TNIK抑制剂) / Nerandomilast (PDE4B抑制剂) 新药数据
  6. Zhou et al. Nat Biotechnol 2026 数字孪生架构参考

用法:
  python main.py               # 运行完整仿真并生成所有图表
  python main.py --gui         # 启动Gradio交互界面
  python main.py --sensitivity # 仅运行敏感性分析
  python main.py --paper       # 仅生成论文草稿
  python main.py --extended    # 运行7变量扩展模型

所有参数来源见 config.py 和《肺数字孪生_参数数据库与文献来源.md》

AI Generated: created with DuMate assistance
"""

import os
import sys
import numpy as np

from config import (
    FIBROSIS_ODE, RESPIRATORY_MODEL, DRUG_INTERVENTION,
    NORMAL_LUNG, IPF_PATHOLOGY, COUPLING, REFERENCES
)
from fibrosis_model import FibrosisODEModel
from respiratory_model import RespiratoryModel
from drug_intervention import DrugIntervention
from visualization import (
    plot_fibrosis_progression,
    plot_pv_curves,
    plot_drug_comparison,
    plot_lung_function_metrics,
    plot_breathing_cycle,
)
from sensitivity_analysis import SensitivityAnalysis
from fvc_simulator import FVCSimulator
from paper_generator import PaperGenerator
from advanced_viz import AdvancedVisualization
# v3.0新增模块
from mechanotransduction_model import ExtendedFibrosisModel
from pk_pd_model import PharmacokineticModel, PKPDCoupledModel


class LungDigitalTwin:
    """肺数字孪生模型 — 顶层整合类 (v3.0)"""

    def __init__(self):
        self.fibrosis_model = FibrosisODEModel()
        self.respiratory_model = RespiratoryModel()
        self.drug_intervention = DrugIntervention()
        self.sensitivity = SensitivityAnalysis()
        self.fvc_sim = FVCSimulator()
        self.paper_gen = PaperGenerator()
        self.adv_viz = AdvancedVisualization()
        # v3.0新增
        self.extended_model = ExtendedFibrosisModel()
        self.pk_pd = {}
        for drug_id in ["nintedanib", "pirfenidone", "rentosertib", "nerandomilast"]:
            self.pk_pd[drug_id] = PKPDCoupledModel(drug_id)
        self.results = {}

    def run_full_simulation(self, t_span=(0, 10), drug_start_time=None):
        """
        运行完整的数字孪生仿真

        Parameters
        ----------
        t_span : tuple
            仿真时间范围 (年)
        drug_start_time : float, optional
            药物介入时间 (年)

        Returns
        -------
        results : dict
            所有仿真结果
        """
        print("=" * 60)
        print("  肺数字孪生模型 (Lung Digital Twin) v3.0 — 仿真启动")
        print("  所有参数来自已发表文献，详见参数数据库文档")
        print("  v3.0新增: YAP/TAZ机械转导+成纤维细胞亚型+PK-PD")
        print("=" * 60)

        # 1. 无干预的纤维化进程
        print("\n[1/9] 模拟无干预纤维化进程...")
        self.fibrosis_model.D = 0.0
        result_no_drug = self.fibrosis_model.simulate(t_span=t_span)
        self.results["no_drug"] = result_no_drug

        # 2. 各药物干预仿真 (含v3.0新增Rentosertib/Nerandomilast)
        print("[2/9] 模拟药物干预效果...")
        drug_ids = list(self.drug_intervention.drugs.keys())
        for drug_id in drug_ids:
            drug = self.drug_intervention.drugs[drug_id]
            modified_params = self.drug_intervention.apply_drug_to_ode_params(drug_id)
            model = FibrosisODEModel(params=modified_params)
            D_factor = modified_params.get("D_intervention", 0)
            model.set_drug_intervention(D_factor)
            result = model.simulate(t_span=t_span, drug_start_time=drug_start_time)
            self.results[drug_id] = result
            print(f"  [OK] {drug['name_cn']} 仿真完成")

        # 3. 联合用药仿真
        print("[3/9] 模拟联合用药...")
        combo_groups = [
            (["huangqi", "danshen"], "黄芪+丹参"),
            (["huangqi", "gusuibu"], "黄芪+骨碎补"),
            (["chuanxiong", "danshen"], "川芎+丹参"),
            (["kushen", "huangqi"], "苦参+黄芪"),
        ]
        for drug_ids_combo, combo_name in combo_groups:
            combo_params = self.drug_intervention.apply_combination(drug_ids_combo)
            combo_model = FibrosisODEModel(params=combo_params)
            D_combo = combo_params.get("D_intervention", 0)
            combo_model.set_drug_intervention(D_combo)
            result_combo = combo_model.simulate(t_span=t_span, drug_start_time=drug_start_time)
            key = f"combo_{'_'.join(drug_ids_combo)}"
            self.results[key] = result_combo
            print(f"  [OK] {combo_name} 联合仿真完成")

        # 4. 呼吸力学验证
        print("[4/9] 呼吸力学模型验证...")
        resp_validation = self.respiratory_model.validate_model()
        for k, v in resp_validation.items():
            print(f"  {k}: {v}")

        # 5. FVC仿真
        print("[5/9] FVC仿真与临床数据对比...")
        for drug_id in [None, "nintedanib", "pirfenidone", "huangqi", "danshen"]:
            # FVC仿真固定使用4年窗(INPULSIS试验为52周≈1年，用4年更合理)
            fvc_result = self.fvc_sim.simulate_fvc(drug_id=drug_id, t_span=(0, 4))
            key = f"fvc_{drug_id or 'placebo'}"
            self.results[key] = fvc_result
            name = drug_id or "安慰剂"
            print(f"  {name}: FVC下降率 = {fvc_result['fvc_decline_rate']:.1f} mL/yr")

        # 6. ODE模型验证
        print("[6/9] ODE模型验证...")
        ode_validation = self.fibrosis_model.validate_against_clinical(result_no_drug)
        print(f"  3.5年ECM密度: {ode_validation['E_at_3yr']:.3f}")
        print(f"  ECM增长倍数: {ode_validation['E_ratio_to_baseline']:.1f}x (临床参考: {IPF_PATHOLOGY['collagen_increase_factor']}x)")
        print(f"  稳态ECM: {ode_validation['E_steady_state']:.3f}")
        if ode_validation["time_to_50pct_fibrosis"]:
            print(f"  达50%纤维化: {ode_validation['time_to_50pct_fibrosis']:.2f}年")

        # === v3.0新增 ===
        # 7. 7变量扩展模型 (YAP/TAZ + 成纤维细胞亚型)
        print("[7/9] 扩展模型仿真 (YAP/TAZ + 成纤维细胞亚型)...")
        ext_result = self.extended_model.simulate(t_span=t_span)
        self.results["extended_no_drug"] = ext_result
        t_ext = ext_result["t"]
        idx_3yr = np.argmin(np.abs(t_ext - 3.5))
        print(f"  ECM增长: {ext_result['E'][idx_3yr]/ext_result['E'][0]:.1f}x")
        print(f"  YAP活性: {ext_result['Y'][idx_3yr]:.3f}")
        print(f"  分泌型/修复型比: {ext_result['F_s'][idx_3yr]/max(ext_result['F_r'][idx_3yr],0.001):.2f}")

        # 8. PK-PD耦合仿真 (新药)
        print("[8/9] PK-PD耦合仿真...")
        for drug_id in ["nintedanib", "pirfenidone", "rentosertib", "nerandomilast"]:
            base_model = FibrosisODEModel()
            pkpd_result = self.pk_pd[drug_id].simulate_pkpd_fibrosis(
                base_model, t_years=t_span[1], drug_start_year=drug_start_time or 2.0
            )
            self.results[f"pkpd_{drug_id}"] = pkpd_result
            print(f"  [OK] {drug_id} PK-PD仿真完成 (avg PD effect: {pkpd_result['avg_pd_effect']:.3f})")

        # 9. 扩展模型+药物干预
        print("[9/9] 扩展模型+药物干预仿真...")
        for drug_id in ["nintedanib", "rentosertib", "nerandomilast"]:
            drug_changes = self.drug_intervention.get_drug_params(drug_id)
            ext_model_drug = ExtendedFibrosisModel()
            # 简化: 使用PK-PD稳态效应修改参数
            D_drug = drug_changes.get("gamma_inhibit", 0) * 0.6
            ext_model_drug.set_drug_intervention(D_drug)
            ext_drug_result = ext_model_drug.simulate(
                t_span=t_span, drug_start_time=drug_start_time or 2.0
            )
            self.results[f"extended_{drug_id}"] = ext_drug_result
            print(f"  [OK] 扩展模型+{drug_id} 完成")

        return self.results

    def generate_all_figures(self, output_dir=None):
        """
        生成所有图表 (包括新增的敏感性分析、FVC对比等)

        Parameters
        ----------
        output_dir : str, optional
            输出目录
        """
        if output_dir is None:
            output_dir = os.path.dirname(os.path.abspath(__file__))
            output_dir = os.path.join(output_dir, "output_figures")

        os.makedirs(output_dir, exist_ok=True)

        print(f"\n生成图表到: {output_dir}")

        # === 原始5张图 ===
        # 图1: 纤维化进程曲线
        print("  [1/11] 纤维化进程曲线...")
        plot_fibrosis_progression(
            self.results["no_drug"],
            save_path=os.path.join(output_dir, "fig1_fibrosis_progression.png"),
            title="肺纤维化进程仿真 (无干预)"
        )

        # 图2: PV曲线对比
        print("  [2/11] PV曲线对比...")
        plot_pv_curves(
            self.respiratory_model,
            save_path=os.path.join(output_dir, "fig2_pv_curves.png")
        )

        # 图3: 药物干预对比 (含新增中药)
        print("  [3/11] 药物干预对比...")
        drug_compare = {
            "无干预": self.results["no_drug"],
            "黄芪": self.results["huangqi"],
            "丹参": self.results["danshen"],
            "骨碎补": self.results["gusuibu"],
            "川芎": self.results["chuanxiong"],
            "苦参": self.results["kushen"],
            "尼达尼布": self.results["nintedanib"],
            "吡非尼酮": self.results["pirfenidone"],
        }
        plot_drug_comparison(
            drug_compare,
            save_path=os.path.join(output_dir, "fig3_drug_comparison.png")
        )

        # 图4: 肺功能指标
        print("  [4/11] 肺功能指标...")
        plot_lung_function_metrics(
            self.respiratory_model,
            save_path=os.path.join(output_dir, "fig4_lung_function_metrics.png")
        )

        # 图5: 呼吸周期
        print("  [5/11] 呼吸周期...")
        plot_breathing_cycle(
            self.respiratory_model,
            save_path=os.path.join(output_dir, "fig5_breathing_cycle.png")
        )

        # === 新增6张图 ===
        # 图6: FVC仿真 vs 临床数据
        print("  [6/11] FVC仿真 vs 临床数据对比...")
        self.fvc_sim.plot_fvc_comparison(
            save_path=os.path.join(output_dir, "fig6_fvc_clinical_comparison.png")
        )

        # 图7: 参数局部敏感性(龙卷风图)
        print("  [7/11] 参数局部敏感性分析(龙卷风图)...")
        local_sens = self.sensitivity.local_sensitivity()
        self.sensitivity.plot_tornado(
            local_sens,
            save_path=os.path.join(output_dir, "fig7_tornado_sensitivity.png")
        )

        # 图8: 全局敏感性(PRCC)
        print("  [8/11] 全局敏感性分析(PRCC)...")
        global_sens = self.sensitivity.global_sensitivity_lhs(n_samples=300)
        self.sensitivity.plot_prcc(
            global_sens,
            save_path=os.path.join(output_dir, "fig8_prcc_sensitivity.png")
        )

        # 图9: 蒙特卡洛置信区间
        print("  [9/11] 蒙特卡洛置信区间...")
        self.sensitivity.plot_monte_carlo_ci(
            n_runs=100, perturbation=0.15,
            save_path=os.path.join(output_dir, "fig9_monte_carlo_ci.png")
        )

        # 图10: 药物对比+误差棒
        print("  [10/11] 药物对比(含95%CI误差棒)...")
        self.adv_viz.plot_drug_comparison_with_ci(
            n_mc=30, perturbation=0.12,
            save_path=os.path.join(output_dir, "fig10_drug_comparison_ci.png")
        )

        # 图11: FVC-ECM关系 + 参数热力图
        print("  [11/11] FVC-ECM关系 + 参数热力图...")
        self.adv_viz.plot_fvc_ecm_relationship(
            save_path=os.path.join(output_dir, "fig11_fvc_ecm_relationship.png")
        )
        self.adv_viz.plot_parameter_heatmap(
            param1="gamma", param2="delta", n_grid=15,
            save_path=os.path.join(output_dir, "fig11b_parameter_heatmap.png")
        )

        # 呼吸动画(GIF)
        print("  [BONUS] 生成呼吸周期动画...")
        try:
            self.adv_viz.generate_breathing_animation(
                save_path=os.path.join(output_dir, "breathing_animation.gif")
            )
        except Exception as e:
            print(f"    动画生成跳过: {e}")

        print(f"\n所有图表已生成到: {output_dir}")

    def run_sensitivity_analysis(self, output_dir=None):
        """
        单独运行敏感性分析

        Parameters
        ----------
        output_dir : str, optional
            输出目录
        """
        if output_dir is None:
            output_dir = os.path.dirname(os.path.abspath(__file__))
            output_dir = os.path.join(output_dir, "output_figures")

        os.makedirs(output_dir, exist_ok=True)

        print("=" * 60)
        print("  参数敏感性分析")
        print("=" * 60)

        # 局部敏感性
        print("\n[1/3] 局部敏感性分析 (OAT ±10%)...")
        local = self.sensitivity.local_sensitivity()
        for k, v in local.items():
            print(f"  {v['param_name']}: sens={v['sens_avg']:+.2f}")
        self.sensitivity.plot_tornado(
            local,
            save_path=os.path.join(output_dir, "sensitivity_tornado.png")
        )

        # 全局敏感性
        print("\n[2/3] 全局敏感性分析 (LHS n=500)...")
        global_res = self.sensitivity.global_sensitivity_lhs(n_samples=500)
        for i, k in enumerate(self.sensitivity.param_keys):
            print(f"  {self.sensitivity.param_names[k]}: PRCC={global_res['prcc'][i]['rho']:+.3f}")
        self.sensitivity.plot_prcc(
            global_res,
            save_path=os.path.join(output_dir, "sensitivity_prcc.png")
        )

        # 蒙特卡洛
        print("\n[3/3] 蒙特卡洛置信区间...")
        mc = self.sensitivity.plot_monte_carlo_ci(
            n_runs=100, perturbation=0.15,
            save_path=os.path.join(output_dir, "monte_carlo_ci.png")
        )
        if mc:
            print(f"  3.5yr Median={mc['median_3yr']:.3f}, 95%CI={mc['ci95_3yr']}")

        # 参数热力图
        print("\n[BONUS] 参数空间热力图...")
        self.adv_viz.plot_parameter_heatmap(
            param1="gamma", param2="delta", n_grid=20,
            save_path=os.path.join(output_dir, "parameter_heatmap.png")
        )

        print("\n敏感性分析完成！")

    def generate_paper(self):
        """生成论文草稿"""
        print("=" * 60)
        print("  论文草稿生成")
        print("=" * 60)
        paper, path = self.paper_gen.generate_full_paper()
        print(f"论文草稿已保存至: {path}")
        return path

    def generate_report(self, output_path=None):
        """
        生成仿真报告

        Parameters
        ----------
        output_path : str, optional
            报告输出路径
        """
        if output_path is None:
            output_dir = os.path.dirname(os.path.abspath(__file__))
            output_path = os.path.join(output_dir, "simulation_report.md")

        result = self.results.get("no_drug")
        if result is None:
            print("请先运行 run_full_simulation()")
            return

        # 获取关键指标
        t = result["t"]
        E = result["E"]
        F = result["F"]
        T = result["T"]

        idx_3yr = np.argmin(np.abs(t - 3.5))
        idx_5yr = np.argmin(np.abs(t - 5.0))

        report = f"""# 肺数字孪生模型 — 仿真报告 (v2.0)

> 生成时间: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M')}
> 本项目代码由AI辅助生成，已标注

## 一、模型概述

本数字孪生模型整合了：
1. **纤维化进程ODE模型** — 基于Suki & Bates (2024)正反馈环路理论
2. **呼吸力学模型** — 肺顺应性随纤维化动态变化
3. **药物干预模块** — 基于中药网络药理学靶点数据（10种中药/西药）
4. **参数敏感性分析** — 局部OAT + 全局LHS(PRCC) + 蒙特卡洛CI
5. **FVC仿真** — 与INPULSIS/ASCEND临床试验数据对比验证
6. **论文自动生成** — 方法学和结果章节草稿

## 二、模型验证

### ODE模型 vs 临床数据

| 指标 | 模型输出 | 临床参考 | 来源 |
|------|----------|----------|------|
| 3.5年ECM密度 | {E[idx_3yr]:.3f} | — | — |
| ECM增长倍数 | {E[idx_3yr]/E[0]:.1f}x | {IPF_PATHOLOGY['collagen_increase_factor']}x | Suki & Bates 2024 |
| 稳态ECM | {E[-1]:.3f} | — | — |

### 呼吸力学验证

| 指标 | 模型值 | 临床参考 |
|------|--------|----------|
| 正常肺顺应性 | {self.respiratory_model.compute_compliance(0.05):.3f} L/cmH2O | 0.2 L/cmH2O |
| IPF肺顺应性 (E=0.8) | {self.respiratory_model.compute_compliance(0.8):.3f} L/cmH2O | 0.05-0.10 L/cmH2O |
| 刚度增加比 (E=0.8) | {0.2/self.respiratory_model.compute_compliance(0.8):.1f}x | 6.4x |

## 三、药物干预效果

| 药物 | E(3.5年) | E(5年) | ECM降低率 |
|------|----------|--------|-----------|
"""

        E_no_drug = self.results["no_drug"]["E"]
        all_drug_ids = list(self.drug_intervention.drugs.keys()) + [
            k for k in self.results if k.startswith("combo_")
        ]
        for drug_id in all_drug_ids:
            if drug_id == "no_drug":
                continue
            if drug_id in self.results:
                if drug_id.startswith("combo_"):
                    drug_name = drug_id.replace("combo_", "").replace("_", "+")
                else:
                    drug_name = self.drug_intervention.drugs.get(drug_id, {}).get("name_cn", drug_id)
                E_drug = self.results[drug_id]["E"]
                reduction = (1 - E_drug[idx_3yr] / E_no_drug[idx_3yr]) * 100 if E_no_drug[idx_3yr] > 0 else 0
                report += f"| {drug_name} | {E_drug[idx_3yr]:.3f} | {E_drug[idx_5yr]:.3f} | {reduction:.1f}% |\n"

        # FVC验证
        report += """
## 四、FVC仿真验证

| 组别 | 模型FVC下降率(mL/yr) | 临床FVC下降率(mL/yr) | 来源 |
|------|---------------------|---------------------|------|
"""
        fvc_placebo = self.results.get("fvc_placebo")
        fvc_nin = self.results.get("fvc_nintedanib")
        fvc_pir = self.results.get("fvc_pirfenidone")

        if fvc_placebo:
            report += f"| 安慰剂 | {fvc_placebo['fvc_decline_rate']:.1f} | 239.9 | INPULSIS-1[15] |\n"
        if fvc_nin:
            report += f"| 尼达尼布 | {fvc_nin['fvc_decline_rate']:.1f} | 114.7 | INPULSIS-1[15] |\n"
        if fvc_pir:
            report += f"| 吡非尼酮 | {fvc_pir['fvc_decline_rate']:.1f} | N/A (ASCEND以%pred为终点) | ASCEND[15b] |\n"

        report += f"""
## 五、参考文献

"""
        for ref_id, ref_text in REFERENCES.items():
            report += f"[{ref_id}] {ref_text}\n"

        report += """
## 六、声明

本项目代码由AI (DuMate) 辅助生成，所有模型参数均来自已发表的学术文献。
模型仅用于学术研究，不作为临床诊断或治疗依据。
"""

        with open(output_path, "w", encoding="utf-8") as f:
            f.write(report)

        print(f"仿真报告已生成: {output_path}")
        return output_path


def launch_gradio():
    """启动Gradio交互界面"""
    try:
        import gradio as gr
    except ImportError:
        print("Gradio未安装，正在安装...")
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install", "gradio"])
        import gradio as gr

    twin = LungDigitalTwin()

    def simulate_fibrosis(t_max, damage_level, alpha, gamma, sigma):
        """Gradio回调: 纤维化进程仿真"""
        params = FIBROSIS_ODE.copy()
        params["alpha"] = alpha
        params["gamma"] = gamma
        params["sigma"] = sigma
        params["t_span"] = (0, t_max)

        model = FibrosisODEModel(params=params)
        model.set_damage_level(damage_level)
        result = model.simulate()

        plot_fibrosis_progression(result, save_path="temp_fibrosis.png")
        return "temp_fibrosis.png"

    def simulate_drug(drug_name, dose, t_max, drug_start):
        """Gradio回调: 药物干预仿真"""
        drug_id_map = {
            "黄芪": "huangqi", "丹参": "danshen", "甘草": "gancao",
            "当归": "danggui", "白术": "baizhu", "骨碎补": "gusuibu",
            "川芎": "chuanxiong", "麦冬": "maitong", "苦参": "kushen",
            "尼达尼布": "nintedanib", "吡非尼酮": "pirfenidone",
            "Rentosertib": "rentosertib", "Nerandomilast": "nerandomilast",
        }

        drug_id = drug_id_map.get(drug_name, "huangqi")

        model_no = FibrosisODEModel()
        result_no = model_no.simulate(t_span=(0, t_max))

        modified = twin.drug_intervention.apply_drug_to_ode_params(drug_id, dose)
        model_drug = FibrosisODEModel(params=modified)
        D = modified.get("D_intervention", 0)
        model_drug.set_drug_intervention(D)
        result_drug = model_drug.simulate(
            t_span=(0, t_max),
            drug_start_time=drug_start if drug_start > 0 else None
        )

        results_dict = {
            "无干预": result_no,
            f"{drug_name} (剂量={dose:.0%})": result_drug,
        }
        plot_drug_comparison(results_dict, save_path="temp_drug.png")
        return "temp_drug.png"

    def simulate_pv(E_value):
        """Gradio回调: PV曲线"""
        plot_pv_curves(
            twin.respiratory_model,
            E_values=[0.05, E_value, max(E_value + 0.2, 0.8)],
            save_path="temp_pv.png"
        )
        return "temp_pv.png"

    def run_sensitivity():
        """Gradio回调: 敏感性分析"""
        sa = SensitivityAnalysis()
        local = sa.local_sensitivity()
        sa.plot_tornado(local, save_path="temp_sensitivity.png")
        return "temp_sensitivity.png"

    def run_fvc():
        """Gradio回调: FVC仿真"""
        twin.fvc_sim.plot_fvc_comparison(save_path="temp_fvc.png")
        return "temp_fvc.png"

    # 构建界面
    with gr.Blocks(title="肺数字孪生模型 v3.0", theme=gr.themes.Soft()) as demo:
        gr.Markdown("# 肺数字孪生模型 (Lung Digital Twin) v3.0")
        gr.Markdown("基于真实文献参数的肺纤维化进程仿真与药物干预评估")
        gr.Markdown("> v3.0: YAP/TAZ机械转导 + PK-PD联合建模 + 新药 | AI辅助生成")

        with gr.Tab("纤维化进程"):
            with gr.Row():
                with gr.Column():
                    t_max = gr.Slider(1, 20, value=10, step=1, label="仿真时长 (年)")
                    damage = gr.Slider(0, 1, value=0.3, step=0.05, label="损伤水平")
                    alpha = gr.Slider(0.1, 2.0, value=FIBROSIS_ODE["alpha"], step=0.1, label="alpha (增殖率)")
                    gamma = gr.Slider(0.05, 1.0, value=FIBROSIS_ODE["gamma"], step=0.02, label="gamma (ECM沉积率)")
                    sigma = gr.Slider(0.05, 1.0, value=FIBROSIS_ODE["sigma"], step=0.05, label="sigma (正反馈率)")
                    btn1 = gr.Button("运行仿真", variant="primary")
                with gr.Column():
                    out1 = gr.Image(label="纤维化进程曲线")

            btn1.click(simulate_fibrosis, [t_max, damage, alpha, gamma, sigma], out1)

        with gr.Tab("药物干预"):
            with gr.Row():
                with gr.Column():
                    drug_name = gr.Dropdown(
                        choices=["黄芪", "丹参", "甘草", "当归", "白术", "骨碎补",
                                 "川芎", "麦冬", "苦参", "尼达尼布", "吡非尼酮",
                                 "Rentosertib", "Nerandomilast"],
                        value="黄芪", label="选择药物"
                    )
                    dose = gr.Slider(0.1, 1.0, value=1.0, step=0.1, label="剂量水平")
                    drug_t = gr.Slider(1, 20, value=10, step=1, label="仿真时长 (年)")
                    drug_start = gr.Slider(0, 5, value=2, step=0.5, label="药物介入时间 (年)")
                    btn2 = gr.Button("运行仿真", variant="primary")
                with gr.Column():
                    out2 = gr.Image(label="药物干预对比")

            btn2.click(simulate_drug, [drug_name, dose, drug_t, drug_start], out2)

        with gr.Tab("PV曲线"):
            with gr.Row():
                with gr.Column():
                    E_val = gr.Slider(0.05, 0.9, value=0.4, step=0.05, label="ECM密度 (纤维化程度)")
                    btn3 = gr.Button("生成PV曲线", variant="primary")
                with gr.Column():
                    out3 = gr.Image(label="PV曲线对比")

            btn3.click(simulate_pv, [E_val], out3)

        with gr.Tab("敏感性分析"):
            with gr.Row():
                btn4 = gr.Button("运行敏感性分析", variant="primary")
            with gr.Row():
                out4 = gr.Image(label="龙卷风图")

            btn4.click(run_sensitivity, [], out4)

        with gr.Tab("FVC仿真"):
            with gr.Row():
                btn5 = gr.Button("运行FVC仿真", variant="primary")
            with gr.Row():
                out5 = gr.Image(label="FVC vs 临床数据")

            btn5.click(run_fvc, [], out5)

    demo.launch(share=False, server_name="0.0.0.0", server_port=7860)


if __name__ == "__main__":
    if "--gui" in sys.argv:
        launch_gradio()
    elif "--sensitivity" in sys.argv:
        twin = LungDigitalTwin()
        twin.run_sensitivity_analysis()
    elif "--paper" in sys.argv:
        twin = LungDigitalTwin()
        twin.generate_paper()
    elif "--extended" in sys.argv:
        # 仅运行7变量扩展模型
        ext = ExtendedFibrosisModel()
        result = ext.simulate()
        print("\n7变量扩展模型结果:")
        t = result["t"]
        idx_3yr = np.argmin(np.abs(t - 3.5))
        for key in ["F", "E", "T", "I", "Y", "F_s", "F_r"]:
            print(f"  {key}(3.5yr) = {result[key][idx_3yr]:.3f}")
    else:
        # 运行完整仿真
        twin = LungDigitalTwin()
        twin.run_full_simulation()
        twin.generate_all_figures()
        twin.generate_report()
