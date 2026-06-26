# ANALYSIS_STEPS.md
# 温江多遗址牛族分析 · 流程步骤 v1.2
# 使用前：检查 logs/ 中各步 DONE 时间戳确认最新运行

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## ✅ 完成状态总览

| Step | 内容 | 状态 | 关键输出 |
|------|------|:----:|---------|
| 0 | 数据准备 | ✅ DONE | 18序列合并+QC报告 |
| 1 | 参考面板构建 | ✅ DONE | 51条统一命名参考面板 |
| 2a–2c | 物种鉴定（BLAST+ML鉴定树） | ✅ DONE | 7 yak / 11 cattle; Z9863/Z9867不一致候选确认 |
| 3 | 系统发育分析（精细树） | ✅ DONE | TIM2+F+I+G4; K0914=YakX谱系(BS=99) |
| 3b | 去除lowQ重建精细树 | ✅ DONE | 不一致候选鲁棒性验证通过 |
| 4 | 单倍型分析 | ✅ DONE | 13个单倍型(16样本); 按谱系分组多样性 |
| 4b | PopART NEXUS输出 | ✅ DONE | 06_haplo/popart_input.nex |
| 4c | **单倍型网络图(Fig8)** | ✅ DONE | Python生成,09_figures/fig8_haplotype_network |
| 5 | BEAST2 BSP种群历史 | ❌ **未执行** | 可选探索性分析; 07_beast/为空 |
| 6 | k-mer ML分类器(RF+SVM) | ✅ DONE | RF-BLAST 100%一致; SVM不可信 |

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 环境配置（首次运行）

```bash
conda create -n bovini python=3.10 -y && conda activate bovini
conda install -c bioconda mafft iqtree modeltest-ng beast2 -y
pip install biopython pandas matplotlib seaborn scikit-learn networkx
```
GUI工具（本地安装）：FigTree / Tracer / BEAUti / DnaSP

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## Step 0 · 数据准备  ✅

操作：
  1. 合并18条FASTA → 00_raw/wenjiangduo_18seqs.fasta
  2. 生成 logs/sample_metadata.csv
  3. QC报告 → 01_qc/qc_report.csv

**结果速查**（logs/step0_prepare.log）：
  - 温江多18个古代牛族线粒体基因组
  - 长度范围：16,322–16,341 bp
  - 覆盖度分层：HIGH(≥100×) 8条 / MID(20–99×) 8条 / LOW(<20×) 2条
  - LOW样本：Z9853(2.1×) / Z9864(10.9×) — 排除于群体统计
  - 3个样本N% ≥ 5%：Z9853(30.0% WARN) / Z9859(8.7% WRN) / Z9865(7.4%)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## Step 1 · 参考序列下载与面板构建  ✅

**参考面板结构**（v2.0，详见 02_ref/README.md）：

| 分组 | 数量 | 内容 |
|------|:----:|------|
| **REF** | 1 | 主外群：非洲水牛 NC_020617 |
| **CORE** | 2 | 家水牛(第二外群) + 印度牛 I1 |
| **LIT1** | 1 | 古代家牦牛 YakA (Chen et al. 2023) |
| **LIT2** | 9 | 古代 YakX 谱系 (Gilardet & Oppenheimer et al. 2025) |
| **LIT3** | 13 | 东亚 P1a 家牛 ×11 + 欧洲原牛 P ×2 |
| **LIT4** | 15 | 中国古代家牛/原牛 C (Cai & Kim et al. 2025) |
| **MOD** | 10 | 现代牦牛 ×5 + 柴达木牛T3 ×5 |
| **合计** | **51** | 统一命名规范; 存档于 `02_ref/renamed/` |

⚠️ 外群变更：从 Bubalus_bubalis_NC006295 **改为** 非洲水牛 NC_020617
  （被 Chen 2023, Cai 2025, Gilardet 2025 一致采用）
⚠️ 禁止直接从 `panel_modern/` / `panel_ancient/` / `core/` 读取建树

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## Step 2 · 物种鉴定（线粒体谱系推断） ✅

### 2a. BLAST（禁止 -remote）
见 03_blast/results.csv → 解析 → species_summary.csv

### 2b. BLAST解析（scripts/step2_parse_blast.py）

**结果速查**（logs/step2_parse_blast.log）：
  - **牦牛 7 条**：Z9853, Z9854, Z9855, Z9856, Z9858, Z9863, K0914
  - **普通牛 11 条**：K0909, Z9850, Z9851, Z9852, Z9857, Z9859, Z9860, Z9861, Z9864, Z9865, Z9867
  - 不一致候选（形态×mtDNA）：
    - Z9863（形态=普通牛）→ BLAST 牦牛 99.83% / 100%cov → DISCORDANT
    - Z9867（形态=牦牛）→ BLAST 普通牛 98.92% / 100%cov → DISCORDANT
  - 低置信度样本：Z9853(95.10%id/72%cov/2.1×), Z9859(97.09%/94%), Z9865(96.27%)

**先验检查报告**（纳入PPT）：
  - Z9863 ✅ 先验cattle → mtDNA=yak → **DISCORDANT**
  - Z9867 ✅ 先验yak → mtDNA=cattle → **DISCORDANT**
  → 两个不一致候选均经BLAST确认，为下一步系统发育验证提供输入

### 2c. 鉴定树（IQ-TREE2）

**结果速查**（logs/step2c_id_tree.log）：
  - 22序列(18WJD+4ref), 15,790 bp, 最佳模型 TIM2+F+G4
  - 两大支系均获 BS=100 支持
  - Z9863 → **牦牛支系** BS=100 ✅ 不一致候选确认
  - Z9867 → **普通牛支系** BS=100 ✅ 不一致候选确认
  - Z9853 → 牦牛支系 BS=85（但覆盖度仅2.1×, 审慎解释）
  - **两个子遗址均发现牦牛+普通牛共存**

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## Step 3 · 系统发育分析（精细树） ✅

脚本：scripts/step3_filter_alignment.py + IQ-TREE2 → scripts/step3_postprocess.py（后处理）

**结果速查**（logs/step3_phylo.log）：
  - 69序列(18WJD+51ref), 16,338 bp, 最佳模型 **TIM2+F+I+G4**
  - Bovini 根 → Bos 支系(BS=100) → 普通牛(BS=100, 45末端) + 牦牛(BS=92, 22末端)
  - Z9863 → **牦牛支系** BS=92 ✅ 不一致候选再次确认
  - Z9867 → **普通牛支系** BS=100 ✅ 不一致候选再次确认
  - K0914 → **LIT2_YakX 支系** BS=99（详见 Step 3b 敏感性测试）

**PPT报告要点**：
  1. 精细树在有51条参考面板语境下确认了两个不一致候选
  2. 西南塔(K0909/K0914)与主遗址同谱系样本聚在同一支系，非异常值
  3. 双外群策略（非洲水牛+家水牛）有效锚定根部

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## Step 3b · 去除lowQ序列重建精细树 ✅

**目的**：排除低质量序列检验拓扑稳定性
**排除序列**（6条）：
  - WJD：Z9853(2.1×), Z9859(21.2×), Z9864(10.9×), Z9865(21.9×)
  - LIT2：AG018_lowQ(9.68×), AG028_lowQ(3.92×)
**重建**：63序列; 模型仍为 TIM2+F+I+G4

**结果速查**（logs/step3_remove_lowQ.log）：
  - Z9863 → **牦牛支系** BS=93 ✅ 鲁棒
  - Z9867 → **普通牛支系** BS=100, Z9850姐妹末梢(BS=99) ✅ 鲁棒
  - **两个不一致候选排除lowQ序列后位置不变 → 非人工假象**

### K0914 YakX 敏感性测试

**观测**（logs/phylo_sensitivity_test.log）：
  - K0914(西南塔, 149× HIGH) 在 Step 3 精细树中与 LIT2_YakX 聚(BS=99)
  - 担忧：两条低覆盖度YakX参考(AG018/AG028_lowQ)可能造成拉偏

**测试**：排除所有6条lowQ序列后重建
**结果**：K0914 仍与 LIT2_YakX 聚(BS=99, sister to UP174+UP175)
**结论**：✅ **敏感性测试通过** — K0914 的 YakX 谱系归属是真实的

**⚠️ PPT 表述注意事项**（n=1, 2025年新发表谱系）：
  - "K0914的线粒体谱系与 Gilardet & Oppenheimer et al. (2025) 报道的 YakX 分支一致"
  - BLAST 无法区分 YakA/YakX（NCBI 数据库未完整收录 YakX 序列）
  - 仅系统发育树含 LIT2 参考面板时才能检出
  - 此发现将温江多个体与青藏高原-新疆古代牦牛谱系联系

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## Step 4 · 单倍型分析  ✅

### 4a. 单倍型分配（scripts/step4_haplotype.py）
  - 排除 LOW: Z9853, Z9864
  - 去除 gap/N/ambig 位点，仅比较干净序列
  - 输出：06_haplo/haplotype_assignment.csv

**结果速查**（logs/step4_haplotype.log）：
  - 16 样本, 清理后 10,668 bp, 变异位点 575
  - **13 个单倍型**, Hd=0.9667
  - **⚠️ 全样本 π=0.0256 不可直接报告**（混合谱系差异造成的假高值）

**按mtDNA谱系分组重复多样性才是有效指标**：

| 指标 | **牦牛组(N=6)** | **普通牛组(N=10)** |
|------|:--------------:|:----------------:|
| 清理位点 | 15,879 bp | 14,385 bp |
| 变异位点 | 98 | 31 |
| 单倍型数 | 6 (全部unique) | 9 (Z9850+Z9867共享H02) |
| **Hd** | **1.0000** | **0.9778** |
| **π** | **0.002120 (≈0.21%)** | **0.000582 (≈0.058%)** |

**解读**：
  - 牦牛 π≈0.21%：谱系内中等多样性，符合预期
  - 普通牛 π≈0.058%：T3谱系极低多样性，与东亚黄牛瓶颈效应一致
  - 两个子遗址间 **0 个共享单倍型**：西南塔独有 2 个(H03=K0909, H04=K0914), 主遗址独有 11 个

**PPT报告要点**：
  - 必须按谱系分组报告 Hd/π，不能报告混合值（否则审稿人会质疑）
  - H02 共享：Z9850(形态?) + Z9867(牦牛形态→cattle mtDNA) 共有同一H02单倍型
  - 清晰说明Z9859/Z9865已由BLAST确认为普通牛，纳入牛组

### 4b. PopART NEXUS 生成（scripts/step4_to_nexus.py）
  - 输出：06_haplo/popart_input.nex
  - TRAITS 分组：Ancient_WJD_Main / Ancient_WJD_Tower
  - LOW 覆盖度样本保留于 NEXUS，但标注排除

### 4c. 单倍型网络图 Fig8（scripts/step4_plot_network.py）

基于最小生成树(MST) + spring layout，Python 实现。
**详见**：PLOT_GUIDE.md §6（单倍型网络图专项）
**输出**：09_figures/fig8_haplotype_network.png / .svg

**特征**：
  - 节点大小 ∝ 频次；色泽：Main=天蓝#4A90D9, Tower=青绿#2A9D8F
  - 连线标注突变步数；LOW样本灰色虚线边框
  - 两个分离聚类：牦牛支系(H01/H04/H07/H08/H13) ↔ 普通牛支系(其他)
  - 支系间距离：~530–545 bp（深层次分化），支系内距离：1–8 bp

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## Step 5 · [可选] BEAST2 BSP 种群历史  ❌ 未执行

⚠️ **本步为可选探索性分析**，因样本量小(≤16条)，置信区间必然较宽。
   **当前状态：07_beast/ 为空，未执行**。

如后续需要执行：
  - 输入：仅牦牛谱系的温江多样本（排除 LOW 和 cattle 谱系样本）
  - 写入 logs/beast_input_ids.txt 后再运行 BEAUti/BEAST2
  - 突变率：1.26×10⁻⁸ sub/site/year

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## Step 6 · [辅助方法] k-mer ML 分类器  ✅

脚本：scripts/step6_kmer_classifier.py
参数：k=4（256维）, 训练集10 MOD参考(5yak+5cattle T3)
分类器：RF(200棵) + SVM(RBF核, probability=True)

**结果速查**（logs/step6_kmer.log）：
  - **RF vs BLAST：18/18 (100%) 一致 ✅** — 独立方法验证通过
  - SVM vs BLAST：12/18 (67%) — ⚠️ **不可信**

**⚠️ SVM 不可信原因**：
  - 10 训练样本 × 256 维特征(p >> n) → Platt scaling 5折CV每折仅2个校准样本
  - 12/18 样本 SVM 置信度锁定在 **56.6%** → 概率塌缩
  - 11 个样本决策函数值完全相同 **+0.1484** → 丧失区分能力
  - 6 个 SVM 分歧样本(Z9859/60/61/64/65/67)均为 SVM 将 cattle 误判为 yak

**结论**：RF 独立可信（抗高维, 投票比例概率天然校准）
  → 最终物种鉴定以 RF + BLAST 共识为准
  → SVM 结果仅可于补充材料说明方法学限制，不纳入结论

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 图表清单（最终版）

| 编号 | 内容 | 生成方式 | 输出文件 |
|:----:|------|----------|---------|
| **Fig1** | 覆盖度分层柱状图 | scripts/plot_figures_1to7.py → fig1_qc_coverage() | fig1_qc_coverage |
| **Fig2** | 物种鉴定组成（饼图+BLAST匹配度） | scripts/plot_figures_1to7.py → fig2_species_composition() | fig2_species_composition |
| **Fig3** | Step2c ML鉴定树（cladogram, BS≥70） | scripts/plot_figures_1to7.py → fig3_id_tree() | fig3_id_tree |
| **Fig3b** | Step3 系统发育精细树（69序列） | scripts/plot_step3_tree.py | fig3_step3_tree |
| **Fig4** | 覆盖度 vs N% 散点图 | scripts/plot_figures_1to7.py → fig4_coverage_vs_n() | fig4_quality_scatter |
| **Fig5** | 子遗址mtDNA谱系柱状对比 | scripts/plot_figures_1to7.py → fig5_site_comparison() | fig5_site_comparison |
| **Fig6** | BLAST匹配度分析 | scripts/plot_figures_1to7.py → fig6_blast_identity() | fig6_blast_identity |
| **Fig7** | 形态鉴定×mtDNA谱系矩阵 | scripts/plot_figures_1to7.py → fig7_morphology_vs_mtdna() | fig7_morphology_vs_mtdna |
| **Fig8** | **单倍型网络图（新）** | **scripts/step4_plot_network.py** | fig8_haplotype_network |
| **—** | Step3b去除lowQ验证树（63序列） | scripts/plot_step3_remove_lowQ.py | fig3_step3_remove_lowQ |

配色规范见 → PLOT_GUIDE.md §2.1（禁止自行改动颜色-含义对应）

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## PPT/报告关键科学发现速查

### 1. 主要发现
  - 温江多遗址(4–5世纪) 18个牛族个体，线粒体谱系 **7 条牦牛 + 11 条普通牛**
  - **两个形态×mtDNA不一致候选确认**：
    - Z9863 形态=普通牛 → mtDNA=牦牛支系 **BS=92–100**, 3种方法(BLAST+ML树+k-mer RF)一致
    - Z9867 形态=牦牛 → mtDNA=普通牛支系 **BS=100**, RF+BLAST一致
  - **西南塔 K0914 归属 YakX 谱系(BS=99)**：敏感测试验证通过，n=1，审慎表述

### 2. 子遗址差异（西南塔/N=2 vs 主遗址/N=16）
  - 两个子遗址 **0 个共享单倍型**
  - 西南塔独有 H03(K0909), H04(K0914); 主遗址独有 11 个
  - **两部分均有牦牛+普通牛共存** → 提示共享同一畜牧传统

### 3. 多样性差异（按谱系分组）
  - 牦牛 π=0.21%：中等, 符合预期
  - 普通牛 T3 π=0.058%：极低, 与东亚黄牛瓶颈一致
  - **不得报告混合π=2.56%**（审稿硬伤）

### 4. 方法学验证
  - **k-mer RF vs BLAST：100% (18/18)** — 独立方法验证通过
  - **K0914 YakX 敏感性测试**：排除 lowQ 后仍稳定(BS=99)
  - **3b remove_lowQ 验证**：不一致候选排除低质量样本后位置不变

### 5. 已知局限性（PPT中须主动说明）
  - n=1 的 YakX 谱系归属
  - BLAST 无法识别 YakX（数据库收录不全）
  - 低覆盖度 Z9853/Z9864 排除于群体统计
  - BSP 未执行（样本量小）
  - SVM 分类器不可信（方法学限制）

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 绘图规范 → 09_figures/PLOT_GUIDE.md（必读）

关键规则速查：
  1. 配色以 PLOT_GUIDE.md §2.1 为准，禁止自行改动
  2. 字体：Times New Roman → Georgia → DejaVu Serif
  3. 文件命名：fig{N}_{英文描述}.png + .svg（双格式）
  4. PNG DPI ≥ 300；输出路径固定为 09_figures/
  5. 禁止使用 jet / rainbow / hsv 色图
  6. Fig8 单倍型网络图 §6 专项：节点大小∝频次，两组分色，连线标注步数

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 日志索引（关键日志文件）

| 文件 | 内容 |
|------|------|
| logs/step0_prepare.log | QC覆盖度分层 |
| logs/step2_parse_blast.log | BLAST物种分配明细 |
| logs/step2c_id_tree.log | 鉴定树结果+不一致候选 |
| logs/step3_phylo.log | 精细树结果+模型 |
| logs/step3_remove_lowQ.log | 去除lowQ验证结果 |
| logs/step4_haplotype.log | 单倍型分析+分组多样性 |
| logs/step6_kmer.log | k-mer RF+SVM结果 |
| logs/phylo_inconsistency_check.txt | 不一致候选完整分析 |
| logs/phylo_inconsistency_remove_lowQ.txt | lowQ验证后确认 |
| logs/phylo_sensitivity_test.log | K0914 YakX敏感性测试 |
| logs/best_model.txt | 人工确认的最优模型 |
| logs/sample_metadata.csv | 18样本完整元数据 |
