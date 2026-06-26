# 02_ref/renamed — 系统发育分析序列集（精细树 Step 2）

## 1. 分析目的

本目录为第二棵精细系统发育树的输入序列集。

- **Step 1（粗略树）**：鉴定样本物种归属（牦牛 vs 家牛 vs 原牛）。
- **Step 2（精细树）**：引入多来源古代与现代参考序列，
  解析单倍群归属、渗入信号及东亚牛母系多样性格局。

---

## 2. 样本构成（共 51 条）

### CORE — 核心参照（2条）

| 文件名 | 登录号 | 物种 | 单倍群 | 用途 |
|--------|--------|------|--------|------|
| `CORE_Bbub_outgroup2_NC006295.fa` | NC_006295 | *Bubalus bubalis* 家水牛 | — | 第二外群，辅助根定位 |
| `CORE_Bind_T1_NC005971.fa` | NC_005971 | *Bos indicus* 印度牛 | **I1** | T/I 单倍群分支锚定 |

- NC_006295 与 REF（非洲水牛 NC_020617）共同构成双外群，
  增强根节点稳定性。
- NC_005971 为 *B. indicus* 标准参考线粒体基因组，
  用于锚定 T 与 I 单倍群的分歧节点。

---

### LIT1 — 古代牦牛参照（1条）

| 文件名 | 登录号 | 物种 | 单倍群 | 来源 |
|--------|--------|------|--------|------|
| `LIT1_Bgru_YakA_BG67.fa` | SRR (PRJNA788155) | *Bos grunniens* | **YakA** | Chen et al. 2023, *Sci Adv* |

- BG-67 为西藏邦噶遗址（~2500 cal BP，海拔 3750 m）最早确认家养牦牛。
- 线粒体单倍群 A，核基因组 PCA 归入现代家牦牛群体。
- 同一遗址同期存在家牛（T3）及牦牛×家牛杂交个体，
  是青藏高原早期混合畜牧系统的直接证据。

---

### LIT2 — 古代 YakX 谱系（9条）

| 文件名 | 登录号 | 遗址/年代 | MIA 深度 | 质量标注 |
|--------|--------|-----------|---------|---------|
| `LIT2_Bgru_YakX_AG018_lowQ.fa` | PX387679 | Denisova EC L12 ~245 ka | 9.68× | ⚠️ lowQ |
| `LIT2_Bgru_YakX_AG028_lowQ.fa` | PX387680 | Denisova EC L13 ~155 ka | 3.92× | ⚠️ lowQ |
| `LIT2_Bgru_YakX_AG036.fa` | PX387681 | Denisova EC L11.3 ~56 ka | 17.20× | — |
| `LIT2_Bgru_YakX_BS791.fa` | PX387676 | Siberia ~44.5 ka (¹⁴C) | 2209× | — |
| `LIT2_Bgru_YakX_UP174.fa` | PX387683 | Siberia ~53 ka | 51.73× | — |
| `LIT2_Bgru_YakX_UP175.fa` | PX387684 | Siberia ~60 ka | 208× | — |
| `LIT2_Bgru_YakX_UP193.fa` | PX387678 | Siberia ~27.3 ka (¹⁴C) | 1590× | — |
| `LIT2_Bgru_YakX_UP199.fa` | PX387674 | Siberia ~75 ka | 280× | — |
| `LIT2_Bgru_YakX_UP203.fa` | PX387675 | Siberia ~39 ka | 196× | — |

- YakX 为 Gilardet & Oppenheimer et al. 2025 (*Genome Biol Evol*,
  DOI: 10.1093/gbe/evaf206) 新发现谱系，与现代牦牛 A–H 单倍群
  互为姐妹群，分歧时间约 675 ka（MIS 17 间冰期）。
- 该谱系在阿尔泰地区存续约 20 万年（~203 ka → ~18 ka），
  末次冰盛期后灭绝，可能对应形态学上的贝加尔牦牛
  (*Bos mutus baikalensis*)。
- AG008（PX387682）因原始覆盖度不足已排除，未纳入本数据集。
- ⚠️ AG018 和 AG028 标注 lowQ：原始 MIA 深度合格，
  已纳入 Gilardet et al. Fig. 2 建树；但低覆盖下采样后
  competitive mapping 不稳定（AG018 downsampled → Bison priscus）。
  建树后需检查二者拓扑位置是否稳定。

---

### LIT3 — 东亚 P1a 家牛 + 欧洲原牛 P（13条）

| 文件名 | 登录号 | 物种/来源 | 单倍群 |
|--------|--------|-----------|--------|
| `LIT3_Btau_P1a_LC537308.fa` | LC537308 | 日本短角牛 JS1 | **P1a** |
| `LIT3_Btau_P1a_LC537309.fa` | LC537309 | 日本短角牛 JS2 | **P1a** |
| `LIT3_Btau_P1a_LC537310.fa` | LC537310 | 日本短角牛 JS3 | **P1a** |
| `LIT3_Btau_P1a_LC537311.fa` | LC537311 | 日本短角牛 JS4 | **P1a** |
| `LIT3_Btau_P1a_LC537312.fa` | LC537312 | 日本短角牛 JS5 | **P1a** |
| `LIT3_Btau_P1a_LC537313.fa` | LC537313 | 日本短角牛 JS6 | **P1a** |
| `LIT3_Btau_P1a_LC537314.fa` | LC537314 | 日本短角牛 JS7 | **P1a** |
| `LIT3_Btau_P1a_LC537315.fa` | LC537315 | 日本短角牛 JS8 | **P1a** |
| `LIT3_Btau_P1a_LC537316.fa` | LC537316 | 日本短角牛 JS9 | **P1a** |
| `LIT3_Btau_P1a_LC537317.fa` | LC537317 | 日本短角牛 JS10 | **P1a** |
| `LIT3_Btau_P1a_Korea_DQ124389.fa` | DQ124389 | 韩国牛 | **P1a** |
| `LIT3_Bpri_P_England_GU985279.fa` | GU985279 | 英格兰原牛（已灭绝） | **P (aurochs)** |
| `LIT3_Bpri_P_Poland_JQ437479.fa` | JQ437479 | 波兰原牛（已灭绝） | **P (aurochs)** |

- P1a 为冰后期欧洲原牛向东亚扩散、约 3700 BP 渗入东北亚
  家牛基因库的谱系（Mannen et al. 2020, *Sci Rep*,
  DOI: 10.1038/s41598-020-78040-8）。
- P1a 在日本短角牛中频率高达 45.9%，在其他东亚品种中极罕见。
- 两条欧洲原牛 P（GU985279、JQ437479）代表 P 单倍群祖先型，
  用于锚定 P1a 与欧洲原牛 P 的分歧节点。

---

### LIT4 — 中国古代家牛 / 原牛 C（15条）

| 文件名 | 遗址 | 省份 | 年代 | 单倍群 | 备注 |
|--------|------|------|------|--------|------|
| `LIT4_Btau_T3_CNN13.fa` | Changning | Qinghai | ~4000 BP | **T3** | |
| `LIT4_Btau_T3_DSG04C.fa` | Dashigou | Shaanxi | ~4000 BP | **T3** | |
| `LIT4_Btau_T3_DSQN4.fa` | Dashanqian | Inner Mongolia | ~3500 BP | **T3** | 已修正（原误标 T1）|
| `LIT4_Btau_T3_DSQN19.fa` | Dashanqian | Inner Mongolia | ~3500 BP | **T3** | 已修正（原误标 T1）|
| `LIT4_Bpri_C_HH13C.fa` | Honghe | Heilongjiang | ~5000 BP | **C** | 东亚原牛 |
| `LIT4_Bpri_C_HT07C.fa` | Houtaomuga | Jilin | ~6000 BP | **C** | 东亚原牛 |
| `LIT4_Btau_T3_NJ09C.fa` | JirentaiGoukou | Xinjiang | ~3300 BP | **T3** | |
| `LIT4_Btau_T1_NJ26C.fa` | JirentaiGoukou | Xinjiang | ~3300 BP | **T1** | |
| `LIT4_Bpri_C_ST07C.fa` | Shuangta | Heilongjiang | ~10000 BP | **C** | 最早东亚原牛 |
| `LIT4_Bpri_C_ST08C.fa` | Shuangta | Heilongjiang | ~10000 BP | **C** | 最早东亚原牛 |
| `LIT4_Bpri_C_ST09C.fa` | Shuangta | Heilongjiang | ~10000 BP | **C** | 最早东亚原牛 |
| `LIT4_Btau_T3_XH08C.fa` | Xiaohe | Xinjiang | ~3500 BP | **T3** | |
| `LIT4_Btau-Bpri_C_XH10C_hybrid.fa` | Xiaohe | Xinjiang | ~3500 BP | **C** | ⚠️ 杂交：核基因组=家牛，线粒体=原牛 C |
| `LIT4_Bpri_C_YGZ02C.fa` | Yangguanzhai | Shaanxi | ~5000 BP | **C** | 东亚原牛 |
| `LIT4_Btau_T3_BY02C.fa` | Shatangbeiyuan | Ningxia | ~3000 BP | **T3** | |

- 数据来源：Cai & Kim et al. 2025 (*Science*,
  DOI: 10.1126/science.adu9904)。
- 东亚原牛（单倍群 C）均属北亚原牛谱系，全新世东亚原牛
  携带约 15% 西方原牛基因流（LGM 后持续东向渗入）。
- XH10C 为迄今最早确认的家牛×原牛杂交个体（Xinjiang ~3500 BP），
  核基因组归属家牛但线粒体为原牛 C，提示雄性原牛渗入家牛群体。
- Shuangta（双塔）三条原牛（~10000 BP）为东亚最早全新世原牛
  基因组，早于家牛引入东亚约 5000 年。

---

### MOD — 现代参照（10条）

| 文件名 | 登录号 | 物种 | 单倍群 |
|--------|--------|------|--------|
| `MOD_Bmut_yak_KM233417.fa` | KM233417 | *Bos mutus* 野牦牛 | YakA |
| `MOD_Bmut_yak_KR106993.fa` | KR106993 | *Bos mutus* 野牦牛 | YakA |
| `MOD_Bmut_yak_KY829451.fa` | KY829451 | *Bos mutus* 野牦牛 | YakA |
| `MOD_Bmut_yak_NC025563.fa` | NC_025563 | *Bos mutus* 野牦牛（参考基因组）| YakA |
| `MOD_Bmut_yak_MK033130.fa` | MK033130 | *Bos mutus* 野牦牛 | YakA |
| `MOD_Btau_T3_OQ513044.fa` | OQ513044 | 柴达木牛 *B. taurus* | **T3** |
| `MOD_Btau_T3_OQ513045.fa` | OQ513045 | 柴达木牛 *B. taurus* | **T3** |
| `MOD_Btau_T3_OQ513046.fa` | OQ513046 | 柴达木牛 *B. taurus* | **T3** |
| `MOD_Btau_T3_OQ513047.fa` | OQ513047 | 柴达木牛 *B. taurus* | **T3** |
| `MOD_Btau_T3_OQ513048.fa` | OQ513048 | 柴达木牛 *B. taurus* | **T3** |

- 野牦牛 5 条覆盖现代 YakA 多样性，用于锚定 YakX 与现代牦牛的
  姐妹群关系。
- 柴达木牛（Qaidam cattle）产自青海柴达木盆地，主要携带 T3
  单倍群，与东亚北方黄牛母系一致（Xia et al. 2021, *Heredity*,
  PMID: 33473210）。
- 柴达木牛单倍群判定依据：Wei X. et al. 2023（NCBI OQ513044–48
  提交记录）及 PMC 综述交叉印证，统一标注为 T3。

---

### REF — 外群（1条）

| 文件名 | 登录号 | 物种 | 用途 |
|--------|--------|------|------|
| `REF_Scaf_outgroup_NC020617.fa` | NC_020617 | *Syncerus caffer* 非洲水牛 | **IQ-TREE2 指定根外群** |

- 非洲水牛为 Bovini 族标准外群，被 Chen et al. 2023、
  Cai et al. 2025 及 Gilardet et al. 2025 一致采用。
- IQ-TREE2 建树时通过 `--outgroup REF_Scaf_outgroup_NC020617`
  指定，与 CORE_Bbub（家水牛）共同提供双外群支持。

---

## 3. 样本数量汇总

| 类别 | 条数 | 内容 |
|------|------|------|
| CORE | 2 | 家水牛 + 印度牛 T1 |
| LIT1 | 1 | 古代家牦牛 YakA |
| LIT2 | 9 | 古代 YakX（含 2 条 lowQ）|
| LIT3 | 13 | P1a 现代牛×11 + 欧洲原牛 P×2 |
| LIT4 | 15 | 古代家牛/原牛 C |
| MOD | 10 | 现代野牦牛×5 + 柴达木牛 T3×5 |
| REF | 1 | 非洲水牛外群 |
| **合计** | **51** | |

---

## 4. 关键文献

| 引用 | 对应样本 | DOI |
|------|---------|-----|
| Gilardet & Oppenheimer et al. 2025, *Genome Biol Evol* | LIT2 (YakX) | 10.1093/gbe/evaf206 |
| Chen et al. 2023, *Sci Adv* | LIT1 (BG67) | 10.1126/sciadv.adi6857 |
| Mannen et al. 2020, *Sci Rep* | LIT3 (P1a) | 10.1038/s41598-020-78040-8 |
| Cai & Kim et al. 2025, *Science* | LIT4 (古代家牛/原牛C) | 10.1126/science.adu9904 |
| Xia et al. 2021, *Heredity* | MOD 柴达木牛 T3 判定依据 | 10.1038/s41437-020-00394-4 |

---

