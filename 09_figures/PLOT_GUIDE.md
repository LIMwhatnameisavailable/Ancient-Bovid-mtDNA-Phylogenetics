# PLOT_GUIDE.md
# 温江多遗址牛族分析 · 可视化规范 v1.0
# 位置：~/LIM/IVPP/analysis/09_figures/PLOT_GUIDE.md

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 0. 当前已生成图一览（最终版）

所有图片输出于 09_figures/，同时保存 300dpi PNG + SVG。

| 编号 | 文件名 | Python函数 | 内容 |
|------|--------|-----------|------|
| Fig1 | fig1_qc_coverage | fig1_qc_coverage() | 覆盖度分层柱状图，HIGH/MID/LOW 三色，20×和100×阈值线 |
| Fig2 | fig2_species_composition | fig2_species_composition() | 左：mtDNA谱系归属饼图；右：BLAST匹配度水平条形图 |
| Fig3 | fig3_id_tree | fig3_id_tree() | ML鉴定树（cladogram），BS≥70标注，不一致候选红字标记 |
| Fig4 | fig4_quality_scatter | fig4_coverage_vs_n() | 覆盖度 vs N% 散点图，6个重点样本引导线标注 |
| Fig5 | fig5_site_comparison | fig5_site_comparison() | 西南塔遗址×主遗址 mtDNA谱系柱状对比 |
| Fig6 | fig6_blast_identity | fig6_blast_identity() | BLAST匹配度直方图 + 逐样本点图 |
| Fig7 | fig7_morphology_vs_mtdna | fig7_morphology_vs_mtdna() | 形态鉴定×mtDNA谱系 4×2矩阵热图 |
| Fig8 | fig8_haplotype_network | step4_plot_network.py | 单倍型网络图（MST+spring layout, Python实现）|

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 1A. 各图详细说明

### Fig1：覆盖度分层柱状图（fig1_qc_coverage）

- **数据源**：logs/sample_metadata.csv → coverage 列
- **布局**：单图，(10, 5) 英寸
- **内容**：
  - 所有样本按覆盖度降序排列的柱状图
  - 三色分层：HIGH（≥100×，天蓝#4A90D9）/ MID（20–99×，青绿#2A9D8F）/ LOW（<20×，灰色#AAAAAA）
  - 两条阈值线：100×（琥珀黄虚线--）、20×（灰色虚线:）
  - 不一致候选（Z9863/Z9867）在x轴标签后加星号"*"，柱子上方红色加粗ID标注
  - LOW覆盖样本（Z9853/Z9864）柱子上方红色标注覆盖度值+"×"
- **图例**：位于右上角，5项（HIGH/MID/LOW色块 + 2条阈值线）

### Fig2：物种鉴定组成（fig2_species_composition）

- **数据源**：species_summary.csv + sample_metadata.csv
- **布局**：1×2 左右排列，(12, 5) 英寸
- **左图**：饼图
  - 仅显示牦牛（紫色#6A4C93）/ 普通牛（琥珀黄#E9C46A）占比
  - 饼边无标签，图例位于饼正下方（lower center, bbox_to_anchor）
- **右图**：BLAST匹配度水平条形图
  - x轴范围 94–100.55%，99% 阈值线（琥珀黄虚线）
  - 样本按匹配度升序排列
  - 不一致候选样本柱体用红框（edgecolor=#E64B35）
  - QC caution 样本（Z9853/Z9859/Z9864/Z9865）用灰色柱体
  - 正常样本按mtDNA谱系着色（紫/琥珀黄）
- **图例**：右下角，4项（牦牛紫/牛琥珀黄/不一致候选红框/QC caution灰）

### Fig3：ML鉴定树（fig3_id_tree）

- **数据源**：04_tree_identification/id_tree.treefile
- **布局**：(12, 9) 英寸，cladogram（所有枝长设为1）
- **内容**：
  - IQ-TREE2 ML树，去除水牛外群后展示
  - Bootstrap支持值≥70显示（灰色斜体，12pt）
  - 样本标签着色：牦牛支系=紫色#6A4C93，普通牛支系=琥珀黄#8A6A16，不一致候选=朱红#E64B35
  - 后缀标记：[*] = 不一致候选（红色加粗），[+] = QC caution（灰色加粗），[?] = 形态不确定
- **图例**：底部单行5列（牦牛紫/牛黄/*=不一致候选/?=形态不确定/+=QC caution）
- **脚注**："Cladogram topology shown for mtDNA clade assignment; branch lengths are not proportional."

### Fig4：覆盖度 vs N% 散点图（fig4_quality_scatter）

- **数据源**：01_qc/qc_report.csv × logs/sample_metadata.csv
- **布局**：(8, 6) 英寸，x轴对数尺度
- **内容**：
  - 所有样本散点，按mtDNA谱系着色（紫=牦牛/琥珀黄=普通牛）
  - 三条参考线：20×（灰色虚线）、100×（琥珀黄虚线）、5% N%（红色虚线）
  - 6个重点样本（Z9853/Z9859/Z9863/Z9864/Z9865/Z9867）用直箭头引导线标注
  - 标注位置手动固定，避免重叠
  - 不一致候选用红色箭头+红色文字，其余用灰色
  - 标签带白色半透明背景框
- **图例**：右上角，5项（牦牛/普通牛/100×/20×/5% N）

### Fig5：子遗址物种对比（fig5_site_comparison）

- **数据源**：sample_grouping.csv × sample_metadata.csv
- **布局**：1×2 左右排列，(10, 5.5) 英寸
- **内容**：
  - 左：西南塔遗址（n=2），右：主遗址（n=16）
  - 每图显示 Yak / Cattle 两个柱体（紫/琥珀黄）
  - 柱顶标注具体数量
  - 西南塔遗址标题附加 "(descriptive comparison only)"
- **图例**：底部居中（Yak紫 / Cattle琥珀黄）
- **脚注**："Discordant candidates are listed in Fig7"（斜体灰色）

### Fig6：BLAST匹配度分布（fig6_blast_identity）

- **数据源**：03_blast/species_summary.csv
- **布局**：1×2 左右排列，(12, 5) 英寸
- **左图**：匹配度直方图（bin=0.5%）
  - 牦牛（紫）与普通牛（琥珀黄）分色叠加
  - 99% 阈值线
- **右图**：逐样本匹配度点图
  - 99% 和 95% 两条参考线
  - 不一致候选加 [*] 红色标记，QC caution 加 [+] 灰色标记
- **图例**：左图含（Yak/牛/99%阈值）

### Fig7：形态鉴定×mtDNA谱系矩阵图（fig7_morphology_vs_mtdna）

- **数据源**：sample_grouping.csv（含 morph_cat 和 mtdna_cat 交叉表）
- **布局**：(7.6, 5.2) 英寸
- **内容**：
  - 4列（Bovine / Cattle / Yak / Uncertain）× 2行（Yak mtDNA / Cattle mtDNA）
  - 每格显示样本计数 n + ID列表
  - 配色：不一致候选=浅红底+红框，形态不确定=浅黄底，concordant=浅绿底，无数据=近白#FAFAFA
  - Z9863/Z9867 在格中红色加粗显示 + 星号
- **图例**：底部中心 2×2（Discordant candidate / Bovine morphology / Uncertain morphology / No sample）

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 1B. 字体系统

优先级链：Times New Roman → Georgia → DejaVu Serif → serif

| 场景 | 规格 |
|---|---|
| 图题 | Bold, 14pt |
| 坐标轴标签 | 变量斜体/单位正体, 12pt |
| 刻度标签 | 正体, 10pt |
| 图例 | 正体, 11pt |
| 标注文字 | 正体, 9pt |
| 子图标签 (a)(b) | Bold, 12pt，左上角 |

流程图/架构图节点文字：Inter 或思源黑体（中文场景）

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 2. 配色系统

### 2.1 项目语义配色（不得更改颜色-含义对应关系）

```python
COLORS = {
    # ── 样本分组（系统发育树 / 网络图核心配色）──
    "ancient_wjd_main":   "#4A90D9",  # 温江多主遗址（天蓝）
    "ancient_wjd_tower":  "#2A9D8F",  # 西南塔遗址（青绿）
    "modern_tibet_yak":   "#1A3A5C",  # 现代西藏牦牛（深蓝）
    "wild_yak":           "#6A4C93",  # 野牦牛（深紫）
    "modern_cattle":      "#E9C46A",  # 现代普通牛（琥珀黄）
    "ancient_other":      "#8B5E3C",  # 其他遗址古代样本（棕）
    "modern_other":       "#D3D3D3",  # 其他现代参考（浅灰）

    # ── 特殊标注 ──
    "prior_yak":          "#6A4C93",  # Z9867先验牦牛（深紫星形）
    "prior_cattle":       "#E9C46A",  # Z9863先验普通牛（琥珀黄星形）
    "hybrid_candidate":   "#E64B35",  # 杂交候选（朱红）
    "low_coverage":       "#AAAAAA",  # 低覆盖度样本（中灰，虚线边框）

    # ── 统计图 ──
    "bar_ancient":        "#4A90D9",  # 古代样本柱（天蓝）
    "bar_modern":         "#D3D3D3",  # 现代对照柱（浅灰）
    "threshold_line":     "#E9C46A",  # 参考线（琥珀黄虚线）
    "bsp_median":         "#4A90D9",  # BSP中位线（天蓝）
    "bsp_hpd":            "#4A90D9",  # BSP 95%HPD阴影（天蓝，alpha=0.2）

    # ── 通用 ──
    "bg_panel":           "#FAFAFA",
    "text_primary":       "#333333",
    "text_secondary":     "#666666",
    "axis_line":          "#444444",
    "grid_line":          "#E0E0E0",
    "success":            "#2A9D8F",
    "warning":            "#E9C46A",
    "error":              "#E64B35",
}
```

### 2.2 连续型色图

| 场景 | 色图 |
|---|---|
| 遗传距离热图 | `Blues` |
| 有符号误差 | `RdBu_r` |
| 多样性分布 | `YlOrRd` |

禁止使用：`jet` / `rainbow` / `hsv`

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 3. 图幅与布局

| 图表类型 | 尺寸（英寸） |
|---|---|
| 单图 | (6, 5) |
| 1×2 对比图 | (12, 5) |
| 1×4 消融对比 | (18, 5) |
| 2×2 | (10, 8) |
| 系统发育树 | (10, 14) |
| 单倍型网络图 | (10, 10) |
| BSP曲线 | (8, 5) |
| 多样性对比柱状图 | (8, 5) |

间距：`layout='constrained'` 优先；
DPI：草稿100，正式PNG 300，出版PNG 600，SVG矢量。

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 4. 坐标轴、线条、图例（速查）

坐标轴：隐藏top/right spine，bottom/left linewidth=1.0，color=#444444
刻度：direction=out，major length=4，minor length=2
网格：默认关闭，折线/柱状图按需启用（--，0.5pt，#E0E0E0，alpha=0.7）
主数据线：linewidth=2.0
参考线：linewidth=1.2，--，#E9C46A
图例：frameon=True，framealpha=0.9，edgecolor=#CCCCCC，fancybox=False

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 5. 系统发育树专项（FigTree导出后AI精修）

节点着色：见 §2.1 样本分组配色
Z9863（先验普通牛）：琥珀黄星形标注
Z9867（先验牦牛）：深紫星形标注
杂交候选：朱红 #E64B35 标注
低覆盖度样本（Z9853/Z9864）：灰色虚线边框
Bootstrap支持值：≥70显示，<70不显示
导出：SVG → 09_figures/Fig3_Tree_raw.svg

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 6. 单倍型网络图专项（PopART导出后AI精修）

分组着色：见 §2.1
节点大小：与单倍型频次成比例
连线长度：代表突变步数，1步=1个突变
导出：SVG → 09_figures/Fig4_Network_raw.svg

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 7. rcParams 全局模板

```python
import matplotlib.pyplot as plt

def apply_global_style():
    plt.rcParams.update({
        'font.family':        'serif',
        'font.serif':         ['Times New Roman','Georgia','DejaVu Serif'],
        'mathtext.fontset':   'stix',
        'axes.unicode_minus': False,
        'font.size':          10,
        'axes.titlesize':     14,
        'axes.labelsize':     12,
        'xtick.labelsize':    10,
        'ytick.labelsize':    10,
        'legend.fontsize':    11,
        'figure.titlesize':   14,
        'axes.linewidth':     1.0,
        'axes.edgecolor':     '#444444',
        'axes.facecolor':     'white',
        'axes.labelcolor':    '#333333',
        'axes.spines.top':    False,
        'axes.spines.right':  False,
        'xtick.direction':    'out',
        'ytick.direction':    'out',
        'xtick.major.size':   4.0,
        'ytick.major.size':   4.0,
        'xtick.color':        '#333333',
        'ytick.color':        '#333333',
        'lines.linewidth':    2.0,
        'legend.frameon':     True,
        'legend.framealpha':  0.9,
        'legend.edgecolor':   '#CCCCCC',
        'legend.fancybox':    False,
        'figure.facecolor':   'white',
        'figure.dpi':         100,
        'savefig.dpi':        300,
        'savefig.bbox':       'tight',
        'savefig.facecolor':  'white',
        'axes.grid':          False,
        'axes.prop_cycle': plt.cycler(color=[
            '#4A90D9','#E64B35','#2A9D8F',
            '#E9C46A','#6A4C93','#8B5E3C',
        ]),
    })
```

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 8. 输出规范

```python
def save_figure(fig, output_dir, stem, dpi=300):
    from pathlib import Path
    p = Path(output_dir); p.mkdir(parents=True, exist_ok=True)
    fig.savefig(p/f"{stem}.png", dpi=dpi, bbox_inches='tight', facecolor='white')
    fig.savefig(p/f"{stem}.svg", format='svg', bbox_inches='tight', facecolor='white')
```

命名规范：`fig{N}_{内容描述}.png/svg`，全小写，下划线分隔，无空格无中文
例：fig1_qc_coverage.png / fig3_ml_tree.svg

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 9. 禁止事项

- 使用 jet / rainbow / hsv 色图
- 纯黑 #000000 作为文字或轴线色（用 #333333 / #444444）
- 仅保存PNG不保存SVG
- PNG DPI低于300
- 文件名含空格、中文或特殊符号
- PPT中使用截图替代矢量图
