"""
肺纤维化进程 ODE 模型
========================
基于 Suki & Bates 2024 [Ref.1] 的正反馈环路理论：
  胶原沉积 → 组织变硬 → 异常力学信号(mechanotransduction) → 更多胶原沉积

模型方程:
  dF/dt = α·T·(1-F) - β·F            肌成纤维细胞动力学
  dE/dt = γ·F·(1-E) - δ·E·D          ECM/胶原沉积动力学
  dT/dt = ε·I + σ·E - ζ·T             TGF-β1动力学 (含ECM正反馈)
  dI/dt = η·damage - θ·I              炎症因子动力学

所有参数来源见 config.py
"""

import numpy as np
from scipy.integrate import solve_ivp
from config import FIBROSIS_ODE, IPF_PATHOLOGY


class FibrosisODEModel:
    """肺纤维化进程的常微分方程模型"""

    def __init__(self, params=None):
        """
        初始化模型参数

        Parameters
        ----------
        params : dict, optional
            自定义参数字典，未提供则使用config.py中的默认值
        """
        p = params if params is not None else FIBROSIS_ODE

        # 肌成纤维细胞参数
        self.alpha = p["alpha"]    # TGF-β驱动的增殖率
        self.beta = p["beta"]      # 凋亡率

        # ECM参数
        self.gamma = p["gamma"]    # ECM沉积率
        self.delta = p["delta"]    # ECM降解率

        # TGF-β参数
        self.epsilon = p["epsilon"]  # 炎症→TGF-β转化率
        self.sigma = p["sigma"]      # ECM→TGF-β正反馈率
        self.zeta = p["zeta"]        # TGF-β衰减率

        # 炎症参数
        self.eta = p["eta"]        # 损伤→炎症转化率
        self.theta = p["theta"]    # 炎症衰减率

        # 初始条件
        self.F0 = p["F0"]
        self.E0 = p["E0"]
        self.T0 = p["T0"]
        self.I0 = p["I0"]

        # 仿真参数
        self.t_span = p["t_span"]
        self.t_eval_points = p["t_eval_points"]

        # 药物干预因子 (0=无干预)
        self.D = 0.0

        # 损伤水平 (0=无损伤, 1=持续强损伤)
        self.damage = 0.3  # 默认中等损伤水平，对应IPF慢性损伤

    def set_drug_intervention(self, D):
        """
        设置药物干预因子

        Parameters
        ----------
        D : float
            药物干预强度 (0=无干预, 1=完全抑制ECM沉积)
        """
        self.D = np.clip(D, 0.0, 1.0)

    def set_damage_level(self, damage):
        """
        设置初始损伤水平

        Parameters
        ----------
        damage : float
            损伤水平 (0=无, 1=持续强损伤)
        """
        self.damage = np.clip(damage, 0.0, 1.0)

    def _ode_system(self, t, y):
        """
        ODE方程组

        Parameters
        ----------
        t : float
            时间 (年)
        y : array-like
            状态向量 [F, E, T, I]

        Returns
        -------
        dydt : list
            导数向量
        """
        F, E, T, I = y

        # 确保状态变量在合理范围内
        F = np.clip(F, 0.0, 1.0)
        E = np.clip(E, 0.0, 1.0)
        T = np.clip(T, 0.0, 1.0)
        I = np.clip(I, 0.0, 1.0)

        # dF/dt: 肌成纤维细胞动力学
        # TGF-β驱动的增殖 - 凋亡
        dFdt = self.alpha * T * (1 - F) - self.beta * F

        # dE/dt: ECM/胶原沉积动力学
        # 成纤维细胞驱动的沉积 - 降解(受药物抑制)
        dEdt = self.gamma * F * (1 - E) - self.delta * E * (1 + self.D)

        # dT/dt: TGF-β1动力学
        # 炎症来源 + ECM正反馈(关键!) - 自然衰减
        # 注: σ·E项是Suki & Bates [Ref.1] 所描述的正反馈环路的核心:
        #   胶原沉积(ECM增加) → 组织变硬 → 异常力学信号 → 更多TGF-β产生
        dTdt = self.epsilon * I + self.sigma * E - self.zeta * T

        # dI/dt: 炎症因子动力学
        # 持续损伤驱动 - 自然衰减
        dIdt = self.eta * self.damage - self.theta * I

        return [dFdt, dEdt, dTdt, dIdt]

    def simulate(self, t_span=None, t_eval_points=None, drug_start_time=None):
        """
        运行仿真

        Parameters
        ----------
        t_span : tuple, optional
            仿真时间范围 (年)
        t_eval_points : int, optional
            输出时间点数
        drug_start_time : float, optional
            药物开始介入时间 (年)。若提供，在该时间点启用药物干预

        Returns
        -------
        result : dict
            包含时间序列和状态变量的字典
        """
        if t_span is None:
            t_span = self.t_span
        if t_eval_points is None:
            t_eval_points = self.t_eval_points

        y0 = [self.F0, self.E0, self.T0, self.I0]
        t_eval = np.linspace(t_span[0], t_span[1], t_eval_points)

        if drug_start_time is not None and drug_start_time > t_span[0]:
            # 分两阶段仿真: 药物介入前 和 药物介入后
            # 阶段1: 无药物
            D_backup = self.D
            self.D = 0.0
            t_eval_1 = t_eval[t_eval <= drug_start_time]
            if len(t_eval_1) == 0 or t_eval_1[-1] < drug_start_time:
                t_eval_1 = np.append(t_eval_1, drug_start_time)

            sol1 = solve_ivp(
                self._ode_system,
                [t_span[0], drug_start_time],
                y0,
                t_eval=t_eval_1,
                method="RK45",
                rtol=1e-8,
                atol=1e-10,
            )

            # 阶段2: 启用药物
            self.D = D_backup
            y0_drug = sol1.y[:, -1]
            t_eval_2 = t_eval[t_eval > drug_start_time]

            if len(t_eval_2) > 0:
                sol2 = solve_ivp(
                    self._ode_system,
                    [drug_start_time, t_span[1]],
                    y0_drug,
                    t_eval=t_eval_2,
                    method="RK45",
                    rtol=1e-8,
                    atol=1e-10,
                )

                # 合并两阶段结果
                t_combined = np.concatenate([sol1.t, sol2.t])
                y_combined = np.concatenate([sol1.y, sol2.y], axis=1)
            else:
                t_combined = sol1.t
                y_combined = sol1.y
        else:
            # 单阶段仿真
            sol = solve_ivp(
                self._ode_system,
                t_span,
                y0,
                t_eval=t_eval,
                method="RK45",
                rtol=1e-8,
                atol=1e-10,
            )
            t_combined = sol.t
            y_combined = sol.y

        result = {
            "t": t_combined,                    # 时间 (年)
            "F": np.clip(y_combined[0], 0, 1),  # 肌成纤维细胞密度
            "E": np.clip(y_combined[1], 0, 1),  # ECM密度
            "T": np.clip(y_combined[2], 0, 1),  # TGF-β1浓度
            "I": np.clip(y_combined[3], 0, 1),  # 炎症因子水平
        }

        return result

    def compute_fibrosis_degree(self, E):
        """
        从ECM密度计算纤维化程度

        Parameters
        ----------
        E : float or array
            ECM密度 (归一化 0-1)

        Returns
        -------
        fibrosis_degree : float or array
            纤维化程度 (0=正常, 1=最严重)
        """
        from config import COUPLING
        return np.clip(E * COUPLING["E_to_fibrosis_degree"], 0.0, 1.0)

    def validate_against_clinical(self, result):
        """
        将模型输出与临床数据进行对比验证

        Parameters
        ----------
        result : dict
            simulate()返回的结果

        Returns
        -------
        validation : dict
            验证指标字典
        """
        # 提取关键时间点的数据
        t = result["t"]
        E = result["E"]
        F = result["F"]

        # 1. 检查纤维化进展时间尺度
        # IPF中位生存期3.5年 [Ref.7]，纤维化应在2-5年内达到严重水平
        idx_3yr = np.argmin(np.abs(t - 3.5))
        E_at_3yr = E[idx_3yr]

        # 2. 检查ECM增长倍数
        # IPF胶原含量增加2.5倍 [Ref.1]
        E_ratio = E_at_3yr / result["E"][0] if result["E"][0] > 0 else 0

        # 3. 检查稳态ECM水平
        E_steady = E[-1]

        # 4. 计算达到50%纤维化的时间
        idx_50 = np.where(E >= 0.5)[0]
        t_50 = t[idx_50[0]] if len(idx_50) > 0 else None

        validation = {
            "E_at_3yr": E_at_3yr,
            "E_ratio_to_baseline": E_ratio,
            "E_steady_state": E_steady,
            "time_to_50pct_fibrosis": t_50,
            "clinical_ref": {
                "median_survival_years": IPF_PATHOLOGY["median_survival_years"],
                "collagen_increase_factor": IPF_PATHOLOGY["collagen_increase_factor"],
            },
        }

        return validation

    def get_model_description(self):
        """返回模型的数学描述（用于论文方法部分）"""
        desc = """
肺纤维化进程数学模型
=====================

基于Suki & Bates (2024)的正反馈环路理论，我们构建了包含4个状态变量的
耦合常微分方程组，描述肺纤维化的动态进程：

状态变量：
  F(t): 肌成纤维细胞密度 (归一化，0-1)
  E(t): 细胞外基质(ECM)/胶原密度 (归一化，0-1)
  T(t): TGF-β1浓度 (归一化，0-1)
  I(t): 炎症因子水平 (归一化，0-1)

方程组：
  dF/dt = α·T·(1-F) - β·F
  dE/dt = γ·F·(1-E) - δ·E·(1+D)
  dT/dt = ε·I + σ·E - ζ·T
  dI/dt = η·damage - θ·I

关键机制：
  1. σ·E项：ECM→TGF-β正反馈（异常力学转导 mechanotransduction）
     胶原沉积→组织变硬→力学信号→更多TGF-β产生 [Suki & Bates 2024]
  2. D项：药物干预因子，增强ECM降解
  3. (1-F)和(1-E)项：容量限制（组织空间有限）

参数值及来源见config.py
"""
        return desc


if __name__ == "__main__":
    # 快速测试
    model = FibrosisODEModel()
    result = model.simulate()
    validation = model.validate_against_clinical(result)

    print("=" * 50)
    print("肺纤维化ODE模型 — 快速验证")
    print("=" * 50)
    print(f"3.5年时ECM密度: {validation['E_at_3yr']:.3f}")
    print(f"ECM增长倍数(相对基线): {validation['E_ratio_to_baseline']:.1f}倍")
    print(f"稳态ECM密度: {validation['E_steady_state']:.3f}")
    print(f"达到50%纤维化时间: {validation['time_to_50pct_fibrosis']:.2f}年" if validation['time_to_50pct_fibrosis'] else "未达到50%")
    print(f"\n临床参考:")
    print(f"  中位生存期: {validation['clinical_ref']['median_survival_years']}年")
    print(f"  胶原增加倍数: {validation['clinical_ref']['collagen_increase_factor']}倍")
