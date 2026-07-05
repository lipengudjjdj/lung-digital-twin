"""
肺数字孪生模型 — 主程序
=========================
整合纤维化进程ODE模型、呼吸力学模型、药物干预模块，
提供完整的数字孪生仿真和交互式界面。

用法:
  python main.py           # 运行完整仿真并生成所有图表
  python main.py --gui     # 启动Gradio交互界面

所有参数来源见 config.py 和《肺数字孪生_参数数据库与文献来源.md》
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


class LungDigitalTwin:
    """肺数字孪生模型 — 顶层整合类"""

    def __init__(self):
        self.fibrosis_model = FibrosisODEModel()
        self.respiratory_model = RespiratoryModel()
        self.drug_intervention = DrugIntervention()
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
        print("  肺数字孪生模型 (Lung Digital Twin) — 仿真启动")
        print("  所有参数来自已发表文献，详见参数数据库文档")
        print("=" * 60)

        # 1. 无干预的纤维化进程
        print("\n[1/4] 模拟无干预纤维化进程...")
        self.fibrosis_model.D = 0.0
        result_no_drug = self.fibrosis_model.simulate(t_span=t_span)
        self.results["no_drug"] = result_no_drug

        # 2. 各药物干预仿真
        print("[2/4] 模拟药物干预效果...")
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
        print("[3/4] 模拟联合用药 (黄芪+丹参)...")
        combo_params = self.drug_intervention.apply_combination(["huangqi", "danshen"])
        combo_model = FibrosisODEModel(params=combo_params)
        D_combo = combo_params.get("D_intervention", 0)
        combo_model.set_drug_intervention(D_combo)
        result_combo = combo_model.simulate(t_span=t_span, drug_start_time=drug_start_time)
        self.results["combo_huangqi_danshen"] = result_combo
        print("  [OK] 黄芪+丹参联合仿真完成")

        # 4. 呼吸力学验证
        print("[4/4] 呼吸力学模型验证...")
        resp_validation = self.respiratory_model.validate_model()
        for k, v in resp_validation.items():
            print(f"  {k}: {v}")

        # ODE模型验证
        ode_validation = self.fibrosis_model.validate_against_clinical(result_no_drug)
        print(f"\nODE模型验证:")
        print(f"  3.5年ECM密度: {ode_validation['E_at_3yr']:.3f}")
        print(f"  ECM增长倍数: {ode_validation['E_ratio_to_baseline']:.1f}x (临床参考: {IPF_PATHOLOGY['collagen_increase_factor']}x)")
        print(f"  稳态ECM: {ode_validation['E_steady_state']:.3f}")
        if ode_validation["time_to_50pct_fibrosis"]:
            print(f"  达50%纤维化: {ode_validation['time_to_50pct_fibrosis']:.2f}年")

        return self.results

    def generate_all_figures(self, output_dir=None):
        """
        生成所有图表

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

        # 图1: 纤维化进程曲线
        print("  [1/5] 纤维化进程曲线...")
        plot_fibrosis_progression(
            self.results["no_drug"],
            save_path=os.path.join(output_dir, "fig1_fibrosis_progression.png"),
            title="肺纤维化进程仿真 (无干预)"
        )

        # 图2: PV曲线对比
        print("  [2/5] PV曲线对比...")
        plot_pv_curves(
            self.respiratory_model,
            save_path=os.path.join(output_dir, "fig2_pv_curves.png")
        )

        # 图3: 药物干预对比
        print("  [3/5] 药物干预对比...")
        drug_compare = {
            "无干预": self.results["no_drug"],
            "黄芪": self.results["huangqi"],
            "丹参": self.results["danshen"],
            "尼达尼布": self.results["nintedanib"],
            "吡非尼酮": self.results["pirfenidone"],
        }
        plot_drug_comparison(
            drug_compare,
            save_path=os.path.join(output_dir, "fig3_drug_comparison.png")
        )

        # 图4: 肺功能指标
        print("  [4/5] 肺功能指标...")
        plot_lung_function_metrics(
            self.respiratory_model,
            save_path=os.path.join(output_dir, "fig4_lung_function_metrics.png")
        )

        # 图5: 呼吸周期
        print("  [5/5] 呼吸周期...")
        plot_breathing_cycle(
            self.respiratory_model,
            save_path=os.path.join(output_dir, "fig5_breathing_cycle.png")
        )

        print(f"\n所有图表已生成到: {output_dir}")

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

        report = f"""# 肺数字孪生模型 — 仿真报告

> 生成时间: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M')}
> ⚠️ 本项目代码由AI辅助生成，已标注

## 一、模型概述

本数字孪生模型整合了：
1. **纤维化进程ODE模型** — 基于Suki & Bates (2024)正反馈环路理论
2. **呼吸力学模型** — 肺顺应性随纤维化动态变化
3. **药物干预模块** — 基于中药网络药理学靶点数据

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
| 正常肺顺应性 | {self.respiratory_model.compute_compliance(0.05):.3f} L/cmH₂O | 0.2 L/cmH₂O |
| IPF肺顺应性 (E=0.8) | {self.respiratory_model.compute_compliance(0.8):.3f} L/cmH₂O | 0.05-0.10 L/cmH₂O |
| 刚度增加比 (E=0.8) | {0.2/self.respiratory_model.compute_compliance(0.8):.1f}x | 6.4x |

## 三、药物干预效果

| 药物 | E(3.5年) | E(5年) | ECM降低率 |
|------|----------|--------|-----------|
"""

        E_no_drug = self.results["no_drug"]["E"]
        for drug_id in ["huangqi", "danshen", "gancao", "nintedanib", "pirfenidone", "combo_huangqi_danshen"]:
            drug_name = self.drug_intervention.drugs.get(drug_id, {}).get("name_cn", drug_id)
            if drug_id in self.results:
                E_drug = self.results[drug_id]["E"]
                reduction = (1 - E_drug[idx_3yr] / E_no_drug[idx_3yr]) * 100 if E_no_drug[idx_3yr] > 0 else 0
                report += f"| {drug_name} | {E_drug[idx_3yr]:.3f} | {E_drug[idx_5yr]:.3f} | {reduction:.1f}% |\n"

        report += f"""
## 四、参考文献

"""
        for ref_id, ref_text in REFERENCES.items():
            report += f"[{ref_id}] {ref_text}\n"

        report += """
## 五、声明

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

        # 生成图表
        plot_fibrosis_progression(result, save_path="temp_fibrosis.png")
        return "temp_fibrosis.png"

    def simulate_drug(drug_name, dose, t_max, drug_start):
        """Gradio回调: 药物干预仿真"""
        drug_id_map = {
            "黄芪": "huangqi", "丹参": "danshen", "甘草": "gancao",
            "当归": "danggui", "白术": "baizhu",
            "尼达尼布": "nintedanib", "吡非尼酮": "pirfenidone",
        }

        drug_id = drug_id_map.get(drug_name, "huangqi")

        # 无干预
        model_no = FibrosisODEModel()
        result_no = model_no.simulate(t_span=(0, t_max))

        # 有干预
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

    # 构建界面
    with gr.Blocks(title="肺数字孪生模型", theme=gr.themes.Soft()) as demo:
        gr.Markdown("# 🫁 肺数字孪生模型 (Lung Digital Twin)")
        gr.Markdown("基于真实文献参数的肺纤维化进程仿真与药物干预评估")
        gr.Markdown("> ⚠️ 本项目由AI辅助生成 | 所有参数来自已发表文献")

        with gr.Tab("纤维化进程"):
            with gr.Row():
                with gr.Column():
                    t_max = gr.Slider(1, 20, value=10, step=1, label="仿真时长 (年)")
                    damage = gr.Slider(0, 1, value=0.3, step=0.05, label="损伤水平")
                    alpha = gr.Slider(0.1, 2.0, value=0.8, step=0.1, label="α (增殖率)")
                    gamma = gr.Slider(0.1, 2.0, value=0.6, step=0.1, label="γ (ECM沉积率)")
                    sigma = gr.Slider(0.1, 2.0, value=0.3, step=0.1, label="σ (正反馈率)")
                    btn1 = gr.Button("运行仿真", variant="primary")
                with gr.Column():
                    out1 = gr.Image(label="纤维化进程曲线")

            btn1.click(simulate_fibrosis, [t_max, damage, alpha, gamma, sigma], out1)

        with gr.Tab("药物干预"):
            with gr.Row():
                with gr.Column():
                    drug_name = gr.Dropdown(
                        choices=["黄芪", "丹参", "甘草", "当归", "白术", "尼达尼布", "吡非尼酮"],
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

    demo.launch(share=False, server_name="0.0.0.0", server_port=7860)


if __name__ == "__main__":
    if "--gui" in sys.argv:
        launch_gradio()
    else:
        # 运行完整仿真
        twin = LungDigitalTwin()
        twin.run_full_simulation()
        twin.generate_all_figures()
        twin.generate_report()
