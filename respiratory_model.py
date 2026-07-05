"""
呼吸力学模型
================
描述肺的呼吸力学行为，将纤维化进程模型的输出（ECM密度）
映射为肺顺应性变化，生成PV曲线和呼吸周期模拟。

核心方程:
  C_lung = C_normal * max(1 - k_stiff * E, C_min_ratio)
  V(P) = C_lung * P + V_FRC
  P(t) = P_baseline + P_amplitude * sin(2π·f·t)

参数来源见 config.py
"""

import numpy as np
from config import RESPIRATORY_MODEL, NORMAL_LUNG, IPF_PATHOLOGY, COUPLING


class RespiratoryModel:
    """呼吸力学仿真模型"""

    def __init__(self, params=None):
        """
        初始化呼吸力学参数

        Parameters
        ----------
        params : dict, optional
            自定义参数字典
        """
        p = params if params is not None else RESPIRATORY_MODEL

        # 顺应性参数
        self.C_normal = p["C_normal"]      # 正常肺顺应性 (L/cmH2O)
        self.C_IPF = p["C_IPF"]            # IPF肺顺应性 (L/cmH2O)
        self.k_stiff = p["k_stiff"]        # 纤维化→顺应性耦合系数

        # 容积参数
        self.V_FRC = p["V_FRC"]            # FRC (L)
        self.V_TLC = p["V_TLC"]            # TLC (L)

        # 呼吸周期参数
        self.RR = p["RR"]                  # 呼吸频率 (次/分)
        self.R_aw = p["R_aw"]              # 气道阻力 (cmH2O·L^-1·s)

        # 压力参数
        self.P_baseline = p["P_baseline"]
        self.P_amplitude = p["P_amplitude"]

        # 耦合约束
        self.C_min_ratio = COUPLING["C_min_ratio"]

    def compute_compliance(self, E):
        """
        根据ECM密度计算肺顺应性

        这是"数字孪生"的核心耦合: 病理状态 → 生理功能

        Parameters
        ----------
        E : float or array-like
            ECM密度 (归一化 0-1)，来自纤维化ODE模型

        Returns
        -------
        C_lung : float or array
            肺顺应性 (L/cmH2O)

        Notes
        -----
        使用指数衰减模型:
          C_lung = C_normal * C_min_ratio + (C_normal - C_normal * C_min_ratio) * exp(-k * E)

        校准依据:
          - E≈0.05 (正常) → C≈0.186 L/cmH2O (接近静态0.2)
          - E≈0.40 (中度纤维化) → C≈0.087 L/cmH2O (2.3x刚度)
          - E≈0.80 (严重纤维化) → C≈0.031 L/cmH2O (6.4x刚度)
          - 对应IPF肺组织刚度升高6.4倍 [Ref.10]
          - IPF肺顺应性实测0.05-0.1 L/cmH2O [Ref.1]
        """
        # 指数衰减模型: 顺应性随ECM指数下降
        # C = C_min + (C_normal - C_min) * exp(-k * E)
        C_min = self.C_normal * self.C_min_ratio
        C_lung = C_min + (self.C_normal - C_min) * np.exp(-self.k_stiff * E)
        return C_lung

    def pv_curve(self, E, P_range=None, n_points=200):
        """
        生成静态PV曲线

        Parameters
        ----------
        E : float
            ECM密度 (归一化 0-1)
        P_range : tuple, optional
            跨肺压范围 (cmH2O)
        n_points : int
            采样点数

        Returns
        -------
        P : array
            跨肺压 (cmH2O)
        V : array
            肺容积 (L)

        Notes
        -----
        正常PV曲线呈S形，这里用简化的线性模型：
          V = C_lung * P + V_FRC
        更精确的版本可在后续迭代中加入Sigmodal拟合
        """
        if P_range is None:
            P_range = (0, 30)  # 跨肺压0-30 cmH2O覆盖临床范围

        C_lung = self.compute_compliance(E)
        P = np.linspace(P_range[0], P_range[1], n_points)

        # 静态PV关系: V = C * P + V_FRC
        V = C_lung * P + self.V_FRC

        # 物理约束: 容积不能超过TLC或低于0
        V = np.clip(V, 0, self.V_TLC)

        return P, V

    def pv_curve_sigmoid(self, E, P_range=None, n_points=200):
        """
        生成S形静态PV曲线（更符合真实肺力学）

        使用Sigmodal拟合，反映肺泡渐进复张的特征

        Parameters
        ----------
        E : float
            ECM密度 (归一化 0-1)
        P_range : tuple, optional
            跨肺压范围 (cmH2O)
        n_points : int
            采样点数

        Returns
        -------
        P : array
            跨肺压 (cmH2O)
        V : array
            肺容积 (L)

        Notes
        -----
        S形PV曲线方程 (Venegas et al. 1998):
          V = V_min + (V_max - V_min) / (1 + exp(-(P - P_inflection) / k))

        纤维化影响:
          - 降低V_max (限制性通气)
          - 增大P_inflection (需要更大压力才能膨胀)
          - 减小k (曲线变陡峭→顺应性降低范围集中)
        """
        if P_range is None:
            P_range = (0, 30)

        P = np.linspace(P_range[0], P_range[1], n_points)
        C_lung = self.compute_compliance(E)

        # 正常参数 (基于Berne & Levy [Ref.13] 和 West [Ref.14])
        V_min = 0.5    # 最小容积 (L), 对应残气量附近
        V_max_normal = self.V_TLC  # 正常最大容积

        # 纤维化对参数的影响
        fibrosis_factor = np.clip(E * COUPLING["E_to_fibrosis_degree"], 0, 1)

        # V_max随纤维化降低 (限制性通气, TLC < 80%预计值 [IPF诊断标准])
        V_max = V_max_normal * (1 - 0.35 * fibrosis_factor)

        # P_inflection随纤维化增大 (肺变硬，需要更大压力)
        P_inflection_normal = 5.0   # 正常拐点压力 (cmH2O)
        P_inflection = P_inflection_normal * (1 + 2.0 * fibrosis_factor)

        # k随纤维化减小 (曲线变陡)
        k_normal = 3.0
        k = k_normal * C_lung / self.C_normal

        # Sigmoid PV曲线
        V = V_min + (V_max - V_min) / (1 + np.exp(-(P - P_inflection) / k))
        V = np.clip(V, 0, self.V_TLC)

        return P, V

    def breathing_cycle(self, E, n_cycles=3, dt=0.02):
        """
        模拟呼吸周期

        Parameters
        ----------
        E : float
            ECM密度 (归一化 0-1)
        n_cycles : int
            模拟呼吸周期数
        dt : float
            时间步长 (秒)

        Returns
        -------
        t : array
            时间 (秒)
        P : array
            跨肺压 (cmH2O)
        V : array
            肺容积 (L)
        flow : array
            气流速度 (L/s)
        """
        C_lung = self.compute_compliance(E)
        f = self.RR / 60.0  # 呼吸频率 (Hz)

        # 总时间
        T_total = n_cycles / f
        t = np.arange(0, T_total, dt)

        # 跨肺压: 正弦波模拟呼吸肌驱动
        P = self.P_baseline + self.P_amplitude * np.sin(2 * np.pi * f * t)

        # 肺容积 (含阻力效应)
        # 动态方程: P = V/C + R * dV/dt
        # 解: V(t) = C * P_amplitude / sqrt(1 + (w*R*C)^2) * sin(wt - phi) + C * P_baseline + V_FRC
        w = 2 * np.pi * f
        impedance = np.sqrt(1 + (w * self.R_aw * C_lung) ** 2)
        phase = np.arctan(w * self.R_aw * C_lung)

        V = (C_lung * self.P_amplitude / impedance
             * np.sin(w * t - phase)
             + C_lung * self.P_baseline
             + self.V_FRC)
        V = np.clip(V, 0, self.V_TLC)

        # 气流速度
        flow = np.gradient(V, dt)

        return t, P, V, flow

    def compute_lung_function_metrics(self, E):
        """
        计算关键肺功能指标

        Parameters
        ----------
        E : float
            ECM密度 (归一化 0-1)

        Returns
        -------
        metrics : dict
            肺功能指标字典
        """
        C_lung = self.compute_compliance(E)
        fibrosis_degree = np.clip(E * COUPLING["E_to_fibrosis_degree"], 0, 1)

        # 肺顺应性
        # 正常: 0.2 L/cmH2O, IPF: 0.05-0.1 L/cmH2O

        # 估算FVC (用力肺活量)
        # IPF: TLC < 80%预计值, FVC下降
        # 参考INPULSIS: FVC年下降150-250 mL/年
        TLC_predicted = self.V_TLC
        TLC_actual = TLC_predicted * (1 - 0.35 * fibrosis_degree)
        FVC = TLC_actual - 0.5  # 粗略: FVC ≈ TLC - RV

        # FVC占预计值百分比
        FVC_predicted = TLC_predicted - 0.5
        FVC_percent = (FVC / FVC_predicted) * 100 if FVC_predicted > 0 else 0

        # 估算DLco (一氧化碳弥散量)
        # IPF早期即出现DLco下降, < 80%预计值
        # 纤维化影响气体交换面积和弥散距离
        DLco_percent = 100 * (1 - 0.5 * fibrosis_degree)

        # 呼吸功
        # W = ∫P dV, 纤维化时顺应性↓, 同样潮气量需要更大压力
        P_required = NORMAL_LUNG["V_T"] / 1000.0 / C_lung  # 产生500mL潮气量所需压力
        W_breathing = 0.5 * P_required * NORMAL_LUNG["V_T"] / 1000.0  # 简化估算

        metrics = {
            "C_lung_L_cmH2O": C_lung,
            "C_lung_mL_cmH2O": C_lung * 1000,
            "compliance_ratio": C_lung / self.C_normal,
            "stiffness_ratio": self.C_normal / C_lung,  # 刚度比 = 1/顺应性比
            "TLC_L": TLC_actual,
            "TLC_percent_predicted": (TLC_actual / TLC_predicted) * 100,
            "FVC_L": FVC,
            "FVC_percent_predicted": FVC_percent,
            "DLco_percent_predicted": DLco_percent,
            "P_for_500mL_cmH2O": P_required,
            "work_of_breathing_J": W_breathing,
            "fibrosis_degree": fibrosis_degree,
            "clinical_ref": {
                "C_normal": f"{self.C_normal} L/cmH2O",
                "C_IPF_range": f"{IPF_PATHOLOGY['C_L_IPF']-0.025}-{IPF_PATHOLOGY['C_L_IPF']+0.025} L/cmH2O",
                "stiffness_increase": f"{IPF_PATHOLOGY['stiffness_increase_factor']}x [Ref.10]",
            },
        }

        return metrics

    def validate_model(self):
        """
        验证呼吸力学模型是否符合已知临床数据

        Returns
        -------
        validation : dict
            验证结果
        """
        # 正常肺验证
        C_normal = self.compute_compliance(0.05)  # 正常ECM密度

        # IPF肺验证
        C_IPF_mild = self.compute_compliance(0.4)   # 轻度纤维化
        C_IPF_severe = self.compute_compliance(0.8)  # 重度纤维化

        validation = {
            "C_normal_model": f"{C_normal:.3f} L/cmH2O",
            "C_normal_ref": f"{self.C_normal} L/cmH2O [Ref.13,14]",
            "C_IPF_mild_model": f"{C_IPF_mild:.3f} L/cmH2O",
            "C_IPF_severe_model": f"{C_IPF_severe:.3f} L/cmH2O",
            "C_IPF_ref_range": "0.05-0.10 L/cmH2O [Ref.1]",
            "stiffness_ratio_mild": f"{self.C_normal/C_IPF_mild:.1f}x",
            "stiffness_ratio_severe": f"{self.C_normal/C_IPF_severe:.1f}x",
            "stiffness_ratio_ref": f"{IPF_PATHOLOGY['stiffness_increase_factor']}x [Ref.10]",
        }

        return validation


if __name__ == "__main__":
    model = RespiratoryModel()
    validation = model.validate_model()

    print("=" * 50)
    print("呼吸力学模型 — 验证")
    print("=" * 50)
    for k, v in validation.items():
        print(f"  {k}: {v}")
