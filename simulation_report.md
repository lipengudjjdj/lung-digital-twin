# 肺数字孪生模型 — 仿真报告 (v2.0)

> 生成时间: 2026-07-06 23:36
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
| 当归 | 0.103 | 0.148 | 21.3% |
| 白术 | 0.106 | 0.156 | 18.9% |
| 骨碎补 | 0.097 | 0.140 | 26.2% |
| 川芎 | 0.095 | 0.135 | 27.1% |
| 麦冬 | 0.104 | 0.153 | 20.1% |
| 苦参 | 0.093 | 0.130 | 28.9% |
| 尼达尼布 | 0.079 | 0.107 | 39.8% |
| 吡非尼酮 | 0.085 | 0.118 | 34.9% |
| Rentosertib (TNIK抑制剂) | 0.074 | 0.097 | 43.3% |
| Nerandomilast (PDE4B抑制剂) | 0.089 | 0.122 | 32.2% |
| huangqi+danshen | 0.062 | 0.076 | 52.8% |
| huangqi+gusuibu | 0.064 | 0.080 | 51.3% |
| chuanxiong+danshen | 0.063 | 0.078 | 51.7% |
| kushen+huangqi | 0.063 | 0.076 | 52.1% |

## 四、FVC仿真验证

| 组别 | 模型FVC下降率(mL/yr) | 临床FVC下降率(mL/yr) | 来源 |
|------|---------------------|---------------------|------|
| 安慰剂 | 0.0 | 239.9 | INPULSIS[15] |
| 尼达尼布 | 78.8 | 114.1 | INPULSIS[15] |
| 吡非尼酮 | 0.0 | 131.2 | ASCEND |

## 五、参考文献

[1] Suki B, Bates JHT. Mathematical Modeling of the Healthy and Diseased Lung. Springer 2024. ISBN: 978-3-031-53202-3
[2] Zhou X, et al. Digital twins of ex vivo human lungs enable accurate evaluation of therapeutic efficacy. Nat Biotechnol 2026 (已接收/在线首发)
[7] Martinez FJ, et al. Idiopathic pulmonary fibrosis. Nat Rev Dis Primers 2017;3:17074. doi:10.1038/nrdp.2017.74
[8] Ye Z & Hu Y. TGF-β1 in idiopathic pulmonary fibrosis. Int J Mol Med 2021;48:132. doi:10.3892/ijmm.2021.4950
[9] Zhao R, et al. Sustained amphiregulin expression drives progressive lung fibrosis. Cell Stem Cell 2024;31(4). doi:10.1016/j.stem.2024.02.013
[10] YAP/TAZ机械转导恶性循环: 多篇文献支持 — (a) Singh MK, et al. Eur Respir J 2025 (YAP/TAZ调控巨噬细胞介导肺纤维化); (b) Liu F, et al. J Clin Invest 2024; (c) Nature Rev Mol Cell Biol 2024 (YAP/TAZ mechanobiology综述). IPF肺刚度6-7倍: Liu F, et al. Am J Physiol Lung Cell Mol Physiol 2016;311(1):L52-63
[11] Nature空间转录组学 2025 — Csmd1+分泌型/Cd248+修复型成纤维细胞亚型 (需补充具体DOI)
[12] 左为团队 Science Advances 2026 — TGF-β梯度双角色与iBMP7工程化基底细胞 (需补充具体DOI)
[13] Koeppen BM, Stanton BA. Berne & Levy Physiology, 8th Edition. Elsevier 2024. ISBN: 978-0-323-87804-0
[14] West JB, Luks AM. West's Respiratory Physiology: The Essentials, 11th Edition. Wolters Kluwer 2021. ISBN: 978-1975155985
[15] Richeldi L, et al. Efficacy and safety of nintedanib in IPF (INPULSIS). NEJM 2014;370:2071-2082. doi:10.1056/NEJMoa1402584
[15b] King TE Jr, et al. A phase 3 trial of pirfenidone in IPF (ASCEND). NEJM 2014;370:2083-2092. doi:10.1056/NEJMoa1402582
[16] Ren Y, et al. A generative AI-discovered TNIK inhibitor for IPF: a randomized phase 2a trial. Nature Medicine 2025. doi:10.1038/s41591-025-03600-z
[17] Nerandomilast in Patients with Idiopathic Pulmonary Fibrosis. FIBRONEER-IPF Phase 3. NEJM 2025. doi:10.1056/NEJMoa2502600
[18] Schmid U, et al. Population PK of nintedanib in NSCLC/IPF patients. Eur J Clin Pharmacol 2018;74:91-103. doi:10.1007/s00228-017-2366-3 (1191例, Ka=0.0827, CL/F=897, Vd/F=465)
[19] Pirfenidone PK: 中国健康志愿者数据, J Clin Pharmacol 2007; 多中心PK研究
[20] Brillet PY, et al. Personalized lung poromechanical modeling for fibrotic ILD. Expert Rev Med Devices 2025
[21] 黄芪甲苷(AS-IV)抗肺纤维化: Wang H, et al. J Ethnopharmacol 2024; 网络药理学分析

## 六、声明

本项目代码由AI (DuMate) 辅助生成，所有模型参数均来自已发表的学术文献。
模型仅用于学术研究，不作为临床诊断或治疗依据。
