# CLAUDE.md
# 温江多遗址牛族线粒体分析 · 主控文件 v1.0

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 目录结构

原始数据（只读）：
  ~/LIM/IVPP/Bovini18-information/
  ├── *.fasta              # 18条单独FASTA
  ├── sheet1-西藏牛族mt.csv
  └── 西藏牛线粒体数据-forTraining.xlsx

工作目录：
  ~/LIM/IVPP/analysis/
  ├── 00_raw/              # 合并后序列（由Step0生成，勿手动修改）
  ├── 01_qc/
  ├── 02_ref/core/         # 原始核心参考序列（4条，历史存档，建树已不再直接使用）
  ├── 02_ref/panel_modern/ # 原始现代面板下载位置（历史存档）
  ├── 02_ref/panel_ancient/# 原始古代面板下载位置（含未重命名文献数据）
  ├── 02_ref/renamed/      # 51条统一命名参考面板（建树使用！见README.md）
  ├── 03_blast/
  ├── 04_align/
  ├── 04_tree_identification/  # Step 2c 鉴定树实际输出目录
  ├── 05_phylo/
  ├── 06_haplo/
  ├── 07_beast/
  ├── 08_ml/
  ├── 09_figures/
  │   └── PLOT_GUIDE.md    # ← 绘图前必读
  ├── logs/
  └── scripts/

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## CC行为规则

1. 每步完成输出 DONE 报告（文件数+关键指标）
2. 报错输出 ERROR 报告（错误信息+排查建议），禁止静默跳过
3. 绘图前必须读取 09_figures/PLOT_GUIDE.md，不得自行设定配色
4. PLOT_GUIDE.md 不存在时，停止绘图并报错
5. 禁止修改 Bovini18-information/ 下任何文件
6. 每步操作须在 `logs/` 下写日志文件（如 `logs/step1_qc.log`），记录：时间戳、输入文件、关键参数、输出文件列表、DONE/ERROR 状态。禁止仅依赖对话输出而不落盘
7. 使用 LIM_env 环境完成本项目所有任务，若有依赖缺失可自行安装

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 分析流程索引

详细步骤见 → ANALYSIS_STEPS.md
样本元数据见 → SAMPLE_META.md
绘图规范见   → 09_figures/PLOT_GUIDE.md

## 关键约束速查

- 参考面板权威定义：**02_ref/README.md**（v2.0，51条统一命名序列），禁止直接使用 02_ref/core/ / panel_modern/ / panel_ancient/ 建树
- 外群：REF_Scaf_outgroup_NC020617（非洲水牛，主外群）+ CORE_Bbub_outgroup2_NC006295（家水牛，第二外群）
- 低覆盖度样本（Z9853 / Z9864）：保留在鉴定分析，不纳入BSP/DnaSP群体统计
- 先验阳性对照：Z9863=普通牛，Z9867=牦牛，每步结果须验证与先验一致性
- 最优模型：modeltest-ng结果由人工确认后写入 logs/best_model.txt，禁止自动解析
- BLAST：禁止 -remote，使用NCBI网页版批量提交
- 单倍型识别：须先去除gap列，禁止直接字符串比较

## 禁止事项

1. 修改 Bovini18-information/ 原始数据
2. 使用 blastn -remote
3. 用字符串直接比较单倍型（须去gap列）
4. 硬编码样本ID（统一从 logs/sample_metadata.csv 读取）
5. 用 grep 自动解析 modeltest-ng 最优模型
6. 将 Z9853/Z9864 纳入群体统计
7. 绕过 PLOT_GUIDE.md 自行配色
8. 静默跳过任何报错
