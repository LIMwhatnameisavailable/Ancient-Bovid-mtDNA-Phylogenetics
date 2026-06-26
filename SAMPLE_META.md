# SAMPLE_META.md
# 温江多遗址牛族分析 · 样本元数据速查 v1.0
# 在任何分析步骤中需要样本信息时，直接查阅本文件

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 完整样本表

| ID    | 初始物种 | 遗址       | 标本   | MT覆盖度 | 覆盖度等级 | 先验标注        |
|-------|----------|------------|--------|----------|------------|-----------------|
| K0909 | 牛       | 西南塔遗址 | 牙齿   | 170.4×   | HIGH       | —               |
| K0914 | 牛       | 西南塔遗址 | 骨块   | 149.5×   | HIGH       | —               |
| Z9850 | 牛？     | 主遗址     | 左腕骨 | 163.2×   | HIGH       | UNCERTAIN       |
| Z9851 | 牛       | 主遗址     | 颞骨   | 126.6×   | HIGH       | —               |
| Z9852 | 牛       | 主遗址     | 颞骨   | 133.0×   | HIGH       | —               |
| Z9853 | 牛       | 主遗址     | 牙齿   | 2.1×     | **LOW**    | ⚠️ 排除群体统计 |
| Z9854 | 牛       | 主遗址     | 颞骨   | 143.8×   | HIGH       | —               |
| Z9855 | 牛       | 主遗址     | 牙齿   | 50.3×    | MID        | —               |
| Z9856 | 牛       | 主遗址     | 跟骨   | 80.8×    | MID        | —               |
| Z9857 | 牛？     | 主遗址     | 桡骨   | 84.1×    | MID        | UNCERTAIN       |
| Z9858 | 牛       | 主遗址     | 肩胛骨 | 84.0×    | MID        | —               |
| Z9859 | 牛？     | 主遗址     | 跟骨   | 21.2×    | MID        | UNCERTAIN       |
| Z9860 | 牛       | 主遗址     | 胫骨   | 107.0×   | HIGH       | —               |
| Z9861 | 牛       | 主遗址     | 牙齿   | 22.3×    | MID        | —               |
| Z9863 | 普通牛   | 主遗址     | 趾骨   | 81.5×    | MID        | ✅ PRIOR_CATTLE |
| Z9864 | 牛？     | 主遗址     | 牙齿   | 10.9×    | **LOW**    | ⚠️ 排除群体统计 |
| Z9865 | 牛       | 主遗址     | 牙齿   | 21.9×    | MID        | —               |
| Z9867 | 牦牛     | 主遗址     | 趾骨   | 109.9×   | HIGH       | ✅ PRIOR_YAK    |

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 分组速查

HIGH覆盖度（≥100×，8条）：
  K0909 K0914 Z9850 Z9851 Z9852 Z9854 Z9860 Z9867

MID覆盖度（20–99×，8条）：
  Z9855 Z9856 Z9857 Z9858 Z9859 Z9861 Z9863 Z9865

LOW覆盖度（<20×，2条，排除群体统计）：
  Z9853（2.1×）  Z9864（10.9×）

先验阳性对照：
  Z9863 → 普通牛（Bos taurus），BLAST/ML结果须与此一致
  Z9867 → 牦牛（Bos mutus），BLAST/ML结果须与此一致

形态鉴定不确定（"牛？"，重点关注）：
  Z9850 Z9857 Z9859 Z9864

两个子遗址：
  西南塔遗址：K0909 K0914
  主遗址：Z9850–Z9867（除Z9862外连续编号，注意无Z9862/Z9866）

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## TRAITS分组（PopART/系统发育树着色用 — 以 02_ref/README.md 为准）

### 温江多遗址样本
Ancient_WJD_Main   : Z9850–Z9867（主遗址全部，排除Z9862/Z9866）
Ancient_WJD_Tower  : K0909 K0914（西南塔）

### 现代参考 → MOD（02_ref/renamed/）
Modern_Yak_YakA    : MOD_Bmut_yak_KM233417, KR106993, KY829451, NC025563, MK033130
Modern_Cattle_T3   : MOD_Btau_T3_OQ513044–48（柴达木牛 T3）

### 古代文献参考
#### LIT1 — 古代家牦牛
Ancient_Yak_YakA   : LIT1_Bgru_YakA_BG67（Chen et al. 2023, ~2500 BP）

#### LIT2 — 古代 YakX 谱系（Gilardet & Oppenheimer et al. 2025）
Ancient_YakX       : LIT2_Bgru_YakX_BS791, UP174, UP175, UP193, UP199, UP203, AG036
Ancient_YakX_lowQ  : LIT2_Bgru_YakX_AG018_lowQ, AG028_lowQ（⚠️ 低质量，建树后需检查拓扑稳定性）

#### LIT3 — 东亚 P1a + 欧洲原牛 P（Mannen et al. 2020）
Ancient_P1a_EastAsia   : LIT3_Btau_P1a_LC537308–317, Korea_DQ124389
Ancient_Aurochs_P      : LIT3_Bpri_P_England_GU985279, Poland_JQ437479

#### LIT4 — 中国古代家牛/原牛 C（Cai & Kim et al. 2025）
Ancient_Cattle_T1  : LIT4_Btau_T1_NJ26C
Ancient_Cattle_T3  : LIT4_Btau_T3_CNN13, DSG04C, DSQN4, DSQN19, NJ09C, XH08C, BY02C
Ancient_Aurochs_C  : LIT4_Bpri_C_HH13C, HT07C, ST07C, ST08C, ST09C, YGZ02C
Ancient_Hybrid_C   : LIT4_Btau-Bpri_C_XH10C_hybrid（⚠️ 核基因组=家牛，线粒体=原牛 C）

### 外群
Outgroup_REF       : REF_Scaf_outgroup_NC020617（非洲水牛，主外群）
Outgroup_CORE      : CORE_Bbub_outgroup2_NC006295（家水牛，第二外群）
Core_Bind_T1       : CORE_Bind_T1_NC005971（印度牛 I1，T/I 分支锚定）

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
## 原始数据路径

FASTA文件：~/LIM/IVPP/Bovini18-information/*.fasta
元数据CSV：~/LIM/IVPP/Bovini18-information/sheet1-西藏牛族mt.csv
