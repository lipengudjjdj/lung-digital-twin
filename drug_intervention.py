"""
药物干预模块
==============
基于中药网络药理学靶点数据，将中药/西药的干预效果映射到ODE模型参数上。

数据来源:
  - 万方数据2025, 136首肺纤维化专利复方 (频次统计)
  - 中药网络药理学分析 (2026年3月, 基于TCMSP/BatMan-TCM平台)
  - INPULSIS试验 (尼达尼布) [Ref.15]
  - ASCEND试验 (吡非尼酮) [Ref.15b]
  - GENESIS-IPF (rentosertib) [Ref.16]
  - FIBRONEER-IPF (nerandomilast) [Ref.17]

重要声明:
  中药(黄芪、丹参等9味)的抑制系数(alpha_inhibit等)为基于网络药理学靶点映射的
  模型参数，非直接实验测量值。具体而言，通过TCMSP等平台获取中药-靶点-通路映射，
  再根据靶点与ODE参数的对应关系推算抑制系数。这些参数反映的是通路层面的相对
  抑制强度，用于模型仿真中的定性/半定量比较，不应被视为精确药效学参数。

原理:
  药物→靶点→信号通路→ODE参数调整
  例如: 黄芪→TGFB1/SMAD3→TGF-β/Smad通路→降低α(增殖率)和σ(正反馈率)
"""

import numpy as np
from config import DRUG_INTERVENTION, FIBROSIS_ODE


class DrugIntervention:
    """药物干预模型：将药物作用映射到ODE参数修改"""

    def __init__(self):
        """初始化药物数据库"""
        self.drugs = DRUG_INTERVENTION.copy()

    def get_drug_list(self):
        """
        获取可用药物列表

        Returns
        -------
        list of dict
            药物信息列表
        """
        result = []
        for key, drug in self.drugs.items():
            result.append({
                "id": key,
                "name_cn": drug["name_cn"],
                "targets": drug.get("targets", []),
                "pathways": drug.get("pathways", []),
                "frequency": drug.get("frequency", "N/A"),
            })
        return result

    def get_drug_params(self, drug_id, dose_level=1.0):
        """
        获取药物对ODE参数的修改量

        Parameters
        ----------
        drug_id : str
            药物ID (如 "huangqi", "danshen", "nintedanib"等)
        dose_level : float
            剂量水平 (0-1, 1=标准剂量, 0.5=半剂量)

        Returns
        -------
        param_changes : dict
            ODE参数修改量字典
        """
        if drug_id not in self.drugs:
            raise ValueError(f"未知药物ID: {drug_id}。可用: {list(self.drugs.keys())}")

        drug = self.drugs[drug_id]
        dose = np.clip(dose_level, 0.0, 1.0)

        changes = {
            "drug_name": drug["name_cn"],
            "drug_id": drug_id,
            "dose_level": dose,
        }

        # 获取所有抑制系数
        for param_key in ["alpha_inhibit", "beta_inhibit", "gamma_inhibit",
                          "delta_inhibit", "epsilon_inhibit", "sigma_inhibit",
                          "zeta_inhibit", "eta_inhibit", "theta_inhibit"]:
            inhibit = drug.get(param_key, 0.0) * dose
            changes[param_key] = inhibit

        # 获取所有增强系数 (如delta_enhance促进ECM降解, beta_enhance促进凋亡)
        for param_key in ["delta_enhance", "beta_enhance", "epsilon_enhance"]:
            enhance = drug.get(param_key, 0.0) * dose
            changes[param_key] = enhance

        return changes

    def apply_drug_to_ode_params(self, drug_id, dose_level=1.0, base_params=None):
        """
        将药物干预应用到ODE参数上，返回修改后的参数

        Parameters
        ----------
        drug_id : str
            药物ID
        dose_level : float
            剂量水平 (0-1)
        base_params : dict, optional
            基础ODE参数，默认使用config.py中的值

        Returns
        -------
        modified_params : dict
            修改后的ODE参数
        """
        if base_params is None:
            base_params = FIBROSIS_ODE.copy()

        changes = self.get_drug_params(drug_id, dose_level)
        modified = base_params.copy()

        # 映射: alpha_inhibit → 降低 alpha
        # 修改后的 alpha = alpha * (1 - inhibit)
        for param_name in ["alpha", "beta", "gamma", "delta",
                           "epsilon", "sigma", "zeta", "eta", "theta"]:
            inhibit_key = f"{param_name}_inhibit"
            if inhibit_key in changes and changes[inhibit_key] > 0:
                original = modified[param_name]
                modified[param_name] = original * (1.0 - changes[inhibit_key])
                modified[f"{param_name}_original"] = original
                modified[f"{param_name}_change"] = -changes[inhibit_key] * original

        # 映射: delta_enhance → 增加 delta (促进ECM降解)
        #        beta_enhance → 增加 beta (促进肌成纤维细胞凋亡)
        for param_name in ["delta", "beta", "epsilon"]:
            enhance_key = f"{param_name}_enhance"
            if enhance_key in changes and changes[enhance_key] > 0:
                original = modified.get(param_name, FIBROSIS_ODE.get(param_name, 0))
                modified[param_name] = original * (1.0 + changes[enhance_key])
                if f"{param_name}_original" not in modified:
                    modified[f"{param_name}_original"] = original

        # 计算综合药物干预因子D (用于ECM降解增强项)
        # D = max(gamma_inhibit, alpha_inhibit) * dose_level
        D = max(changes.get("gamma_inhibit", 0), changes.get("alpha_inhibit", 0))
        # 如果有delta_enhance，D因子也增加
        D += changes.get("delta_enhance", 0) * 0.5
        modified["D_intervention"] = D

        return modified

    def apply_combination(self, drug_ids, dose_levels=None):
        """
        联合用药: 多种药物组合应用

        Parameters
        ----------
        drug_ids : list of str
            药物ID列表
        dose_levels : list of float, optional
            对应剂量水平列表，默认均为1.0

        Returns
        -------
        modified_params : dict
            修改后的ODE参数
        """
        if dose_levels is None:
            dose_levels = [1.0] * len(drug_ids)

        if len(drug_ids) != len(dose_levels):
            raise ValueError("药物ID列表和剂量水平列表长度必须一致")

        # 从基础参数开始
        modified = FIBROSIS_ODE.copy()

        # 累计各药物的抑制效果
        total_inhibit = {}
        for param_name in ["alpha", "beta", "gamma", "delta",
                           "epsilon", "sigma", "zeta", "eta", "theta"]:
            total_inhibit[param_name] = 0.0

        combo_info = []
        for drug_id, dose in zip(drug_ids, dose_levels):
            changes = self.get_drug_params(drug_id, dose)
            combo_info.append({
                "drug": changes["drug_name"],
                "dose": dose,
            })

            for param_name in total_inhibit:
                inhibit_key = f"{param_name}_inhibit"
                total_inhibit[param_name] += changes.get(inhibit_key, 0.0)

            # 增强效果
            for param_name in ["delta", "beta", "epsilon"]:
                enhance_key = f"{param_name}_enhance"
                if enhance_key not in modified:
                    modified[enhance_key] = 0.0
                modified[enhance_key] += changes.get(enhance_key, 0.0)

        # 限制总抑制不超过90%（防止过度抑制导致不合理结果）
        for param_name in total_inhibit:
            total_inhibit[param_name] = min(total_inhibit[param_name], 0.9)

        # 应用到参数
        for param_name in total_inhibit:
            original = modified[param_name]
            modified[param_name] = original * (1.0 - total_inhibit[param_name])
            modified[f"{param_name}_original"] = original
            modified[f"{param_name}_inhibit_total"] = total_inhibit[param_name]

        # 应用增强效果
        for param_name in ["delta", "beta", "epsilon"]:
            enhance_key = f"{param_name}_enhance"
            if enhance_key in modified and modified[enhance_key] > 0:
                original = modified[param_name]
                modified[param_name] = original * (1.0 + modified[enhance_key])
                if f"{param_name}_original" not in modified:
                    modified[f"{param_name}_original"] = original

        # 综合D因子
        D = max(total_inhibit.get("gamma", 0), total_inhibit.get("alpha", 0))
        D += modified.get("delta_enhance", 0) * 0.5
        modified["D_intervention"] = D
        modified["combination"] = combo_info

        return modified

    def get_therapeutic_effect_summary(self, drug_id, dose_level=1.0):
        """
        获取药物干预效果的文字总结（用于界面展示）

        Parameters
        ----------
        drug_id : str
            药物ID
        dose_level : float
            剂量水平

        Returns
        -------
        summary : str
            效果总结文本
        """
        changes = self.get_drug_params(drug_id, dose_level)
        drug = self.drugs[drug_id]

        lines = [f"【{drug['name_cn']}】干预效果 (剂量水平: {dose_level:.0%})"]
        lines.append(f"  靶点: {', '.join(drug.get('targets', []))}")
        lines.append(f"  通路: {', '.join(drug.get('pathways', []))}")

        if drug.get("frequency"):
            lines.append(f"  临床频次: {drug['frequency']}/136首方剂")

        # 主要抑制效果
        effects = []
        for param_name, label in [
            ("alpha", "肌成纤维细胞增殖"),
            ("gamma", "ECM沉积"),
            ("sigma", "ECM→TGF-β正反馈"),
            ("eta", "炎症驱动"),
        ]:
            inhibit_key = f"{param_name}_inhibit"
            val = changes.get(inhibit_key, 0)
            if val > 0:
                effects.append(f"  ↓{label}: 抑制{val:.0%}")

        lines.extend(effects)

        # 临床参考数据
        if "FVC_decline_treated" in drug:
            lines.append(f"  临床参考: FVC年下降率 {drug['FVC_decline_treated']} mL/年")
            lines.append(f"  数据来源: {drug.get('ref', '')}")

        return "\n".join(lines)


if __name__ == "__main__":
    # 快速测试
    di = DrugIntervention()

    print("=" * 50)
    print("药物干预模块 — 测试")
    print("=" * 50)

    # 测试单药
    print("\n--- 黄芪干预效果 ---")
    print(di.get_therapeutic_effect_summary("huangqi"))

    print("\n--- 尼达尼布干预效果 ---")
    print(di.get_therapeutic_effect_summary("nintedanib"))

    # 测试联合用药
    print("\n--- 黄芪+丹参联合 ---")
    combo = di.apply_combination(["huangqi", "danshen"])
    print(f"  alpha: {FIBROSIS_ODE['alpha']:.3f} → {combo['alpha']:.3f}")
    print(f"  gamma: {FIBROSIS_ODE['gamma']:.3f} → {combo['gamma']:.3f}")
    print(f"  sigma: {FIBROSIS_ODE['sigma']:.3f} → {combo['sigma']:.3f}")
