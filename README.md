# 🫁 肺数字孪生模型 (Lung Digital Twin)

> **AI生成声明**: 本项目代码由 AI (DuMate) 辅助生成，已在每个文件中标注。所有模型参数均来自已发表的学术文献，每项参数均标注出处。

## 项目简介

基于真实文献参数的肺纤维化数字孪生模型，用于模拟和研究特发性肺纤维化 (IPF) 的病理进程、呼吸力学变化及药物干预效果。

### 核心特性

- **纤维化进程ODE模型**: 基于 Suki & Bates (2024) 正反馈环路理论，用4个耦合微分方程描述 TGF-β → 成纤维细胞增殖 → ECM沉积 的动态过程
- **呼吸力学模型**: 肺顺应性随纤维化程度动态变化，生成PV曲线和呼吸周期仿真
- **药物干预模块**: 基于中药网络药理学靶点数据，支持单药和联合用药评估
- **交互式界面**: Gradio Web界面，支持参数调节和实时仿真

### 创新点

1. **首次将纤维化进程ODE模型与呼吸力学模型耦合**，实现"病理-功能"双向映射
2. **融合中药网络药理学靶点**，实现中药干预的数字孪生模拟（文献空白领域）
3. **交互式可视化平台**，直观展示纤维化发展和药物干预效果

## 项目结构

```
lung_digital_twin/
├── config.py               # 全局参数配置（所有真实参数+文献来源）
├── fibrosis_model.py       # 纤维化进程ODE模型
├── respiratory_model.py    # 呼吸力学模型
├── drug_intervention.py    # 药物干预模块
├── visualization.py        # 可视化模块
├── main.py                 # 主程序入口
├── README.md               # 本文件
├── AI_GENERATED.md         # AI生成声明
└── output_figures/         # 输出图表（运行后生成）
```

## 快速开始

### 环境要求

- Python 3.9+
- 依赖包: numpy, scipy, matplotlib

### 安装

```bash
pip install numpy scipy matplotlib
# 可选（交互界面）:
pip install gradio
```

### 运行

```bash
# 运行完整仿真并生成所有图表
python main.py

# 启动Gradio交互界面
python main.py --gui
```

### 快速测试

```python
from fibrosis_model import FibrosisODEModel

model = FibrosisODEModel()
result = model.simulate()

# 查看关键指标
validation = model.validate_against_clinical(result)
print(f"3.5年ECM密度: {validation['E_at_3yr']:.3f}")
print(f"ECM增长倍数: {validation['E_ratio_to_baseline']:.1f}x")
```

## 模型说明

### 纤维化进程ODE模型

```
dF/dt = α·T·(1-F) - β·F            肌成纤维细胞动力学
dE/dt = γ·F·(1-E) - δ·E·(1+D)      ECM沉积动力学
dT/dt = ε·I + σ·E - ζ·T             TGF-β1动力学 (含正反馈)
dI/dt = η·damage - θ·I              炎症因子动力学
```

其中 σ·E 项是核心创新：ECM沉积→组织变硬→异常力学信号(mechanotransduction)→更多TGF-β产生

### 呼吸力学模型

```
C_lung(E) = C_normal × max(1 - k_stiff × E, C_min_ratio)
V(P) = C_lung × P + V_FRC
```

### 药物干预映射

药物 → 靶点 → 信号通路 → ODE参数修改

| 药物 | 主要靶点 | 主要通路 | 核心抑制 |
|------|----------|----------|----------|
| 黄芪 | TGFB1, AKT1 | TGF-β/Smad | ↓增殖率α, ↓正反馈σ |
| 丹参 | PI3K, AKT1 | PI3K/AKT | ↓ECM沉积率γ |
| 甘草 | RELA, NFKB1 | NF-κB | ↓炎症η |
| 尼达尼布 | PDGFR, FGFR | RTK通路 | ↓增殖率α, ↓沉积率γ |
| 吡非尼酮 | TGFB1 | TGF-β/Smad | ↓增殖率α |

## 参数来源

所有参数均来自已发表文献，详见《肺数字孪生_参数数据库与文献来源.md》

核心参考文献：
1. Suki B, Bates JHT. Mathematical Modeling of the Healthy and Diseased Lung, Springer 2024
2. Zhou X, et al. Digital twins of ex vivo human lungs. Nat Biotechnol 2026
3. Richeldi L, et al. INPULSIS Trial. NEJM 2014

## 免责声明

本项目仅用于学术研究目的，不作为临床诊断或治疗依据。模型结果需结合临床实际情况进行解读。

## 许可证

MIT License
