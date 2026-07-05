# 肺数字孪生模型 — 仿真报告

> 生成时间: 2026-07-05 23:49
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
| 3.5年ECM密度 | 0.131 | — | — |
| ECM增长倍数 | 2.6x | 2.5x | Suki & Bates 2024 |
| 稳态ECM | 0.454 | — | — |

### 呼吸力学验证

| 指标 | 模型值 | 临床参考 |
|------|--------|----------|
| 正常肺顺应性 | 0.156 L/cmH2O | 0.2 L/cmH2O |
| IPF肺顺应性 (E=0.8) | 0.032 L/cmH2O | 0.05-0.10 L/cmH2O |
| 刚度增加比 (E=0.8) | 6.2x | 6.4x |

## 三、药物干预效果

| 药物 | E(3.5年) | E(5年) | ECM降低率 |
|------|----------|--------|-----------|
| 黄芪 | 0.092 | 0.129 | 30.0% |
| 丹参 | 0.089 | 0.127 | 31.6% |
| 甘草 | 0.106 | 0.152 | 19.3% |
| 尼达尼布 | 0.079 | 0.107 | 39.8% |
| 吡非尼酮 | 0.085 | 0.118 | 34.9% |
| combo_huangqi_danshen | 0.062 | 0.076 | 52.8% |

## 四、参考文献

[1] Suki B, Bates JHT. Mathematical Modeling of the Healthy and Diseased Lung, Springer 2024, Ch.8
[2] Zhou X, et al. Digital twins of ex vivo human lungs. Nat Biotechnol 2026
[7] Martinez FJ, et al. Idiopathic pulmonary fibrosis. Nat Rev Dis Primers 2017;3:17074
[8] Ye Z & Hu Y. TGF-β1 in idiopathic pulmonary fibrosis. Int J Mol Med 2021;48:132
[9] Zhao R, et al. Sustained amphiregulin expression drives progressive fibrosis. Cell Stem Cell 2024
[10] Pulmonary-Targeted NPs for IPF Therapy. Advanced Science 2025
[13] Berne & Levy Physiology, 8th Edition (Koeppen & Stanton)
[14] West's Respiratory Physiology: The Essentials, 10th Edition
[15] Richeldi L, et al. INPULSIS Trial. NEJM 2014;370:2071-2082

## 五、声明

本项目代码由AI (DuMate) 辅助生成，所有模型参数均来自已发表的学术文献。
模型仅用于学术研究，不作为临床诊断或治疗依据。
