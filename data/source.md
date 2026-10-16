# 数据来源与当前五病分析口径

当前主分析由 `scripts/01_cohort/reconstruct_cohorts.py`、`scripts/02_models/run_five_disease_models.R` 和 `scripts/05_figures/render_main_figures.R` 顺序重建。输入来源、参数和样本筛选记录在 `output/five_disease_revision/reconstruction_manifest.json` 与 `output/revised/README.md`。

## 覆盖率和抗原成分

唯一覆盖率快照为 `raw/coverage_wuenic_all.csv`，仅保留 WHO/UNICEF WUENIC 国家层面估计。麻疹使用 MCV1、MCV2；风疹使用 RCV1；百日咳和白喉使用 DTP1、DTP3；腮腺炎同时保留严格成分证据层和全 MCV 代理层。DTP3、MCV1 等年度边际覆盖率进入病例模型，MCV2、DTP1 及腮腺炎成分记录用于出生队列重建。

每剂覆盖率按目标月龄映射到 1–14 岁出生队列。两剂流使用边际部分序列状态：零剂、完成第一剂和完成第二剂。个体联合剂次历史未观测，产品成分不确定或历史转换未解析时保留未知状态。腮腺炎严格层的国家证据、时间锚点和来源见 `data/curated/mumps_program_evidence.csv`；全 MCV 层固定使用 12 月和 48 月作为分析情景。

[WHO Immunization Data Portal](https://immunizationdata.who.int/) 提供下载，[WHO/UNICEF WUENIC 方法说明](https://www.who.int/news-room/questions-and-answers/item/who-unicef-estimates-of-national-immunization-coverage) 解释估计过程。本地快照为 `who_immunization_portal_coverage.xlsx` 的整理结果。

## 病例、人口和出生队列权重

病例源表为 `raw/case_data.csv`，质量保留表为 `processed/case_data_quality_preserved.csv`，覆盖麻疹、腮腺炎、百日咳、白喉和风疹的国家年度报告。显式零与空白分别编码为 `reported_zero` 和 `unreported`；主模型保留显式零，空白在病例模型中保持缺失。结局是报告病例数；模型使用 WPP 总人口作为 offset。国家监测能力、病例定义与报告行为会影响计数。

总人口、国家代码和出生人数来自 [UN World Population Prospects 2024](https://population.un.org/wpp/)。出生队列权重为年度出生人数乘 1,000，直接取自 `raw/births_data.csv` 并逐行核对。当前重建没有加入存活或迁移调整，因此年龄汇总表示出生队列规模加权的保护缺口情景，而非当前年龄人口易感比例。

国家范围为 WUENIC、病例、WPP 人口和元数据交集，共 195 个国家/地区；每病种 1,950 个国家–年份组合，观察期为 2014–2023。质量保留病例与完整协变量筛选后覆盖 178 个国家，4,205 个国家–病种–年份记录；病种计数为百日咳 1,378、白喉 1,347、风疹 920、麻疹 306 和腮腺炎 proxy 254。

## 国家协变量

国家年度协变量包括人均 GDP、卫生支出、UHC 服务覆盖指数和城市人口比例，来源及快照信息保存在 `raw/country_covariates.csv`、`country_indicator_metadata.csv` 与 `world_bank/`。模型同时使用前一年报告病例和疾病特异的历史覆盖缺口：百日咳使用 5 年窗口，其余病种使用 10 年窗口。

## 证据和刷新

腮腺炎严格层保存美国、英国、芬兰和中国的历史产品证据；完整 14 岁覆盖的 10 个国家–年份全部来自英国，≥80% 层额外包括中国、芬兰和美国。参数情景见 `data/curated/five_disease_parameters.json`，包括剂次效果、持续时间和标准化接种月龄。

主分析入口使用本地快照，不自动联网刷新。更新快照后，应重新运行队列重建、模型、图形和审计，并同步修改稿件中的样本量与数值。

## 年龄结构传播机制层外部输入

独立的年龄结构 SEIRV 机制层使用 UN WPP 2024 单岁人口 `raw/population_by_single_age.csv`，按 0–4、5–9、…、75+ 汇总。年龄接触矩阵取自 Prem、Cook 和 Jit 的 152 国 projected all-location 矩阵，官方补充材料下载地址为 [PLOS supplementary S2](https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1005697.s002&type=supplementary)，本地原始文件为 `raw/contact_matrices_prem2017/MUestimates_all_locations_1.xlsx` 和 `MUestimates_all_locations_2.xlsx`。模型以各国 2019 年 WPP 年龄人口校正矩阵互易性，校正后的五个国家矩阵及误差检查位于 `output/revised/age_seirv/contact_matrices_reciprocal/` 和 `contact_matrix_audit.csv`。

机制层读取五病出生队列保护缺口；腮腺炎使用 all-MCV proxy。缺失的年龄保护档案使用同一国家–病种的年度覆盖率 fallback，并在内部标记。病例数据不参与该层的 beta 或报告率估计，显式零保留为计数，空白保留为缺失。输出为结构性年度流行潜势，传播阻断系数和成人易感性作为敏感性情景。

## 独立血清学与扩展数据

主血清学层使用 [CDC/NCHS NHANES 2017–March 2020 MMRV surplus-specimen file](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_SSMMRV.htm) 及其人口学文件。分析保留 6–59 岁的完整正权重设计框架，以 SDMVSTRA、SDMVPSU 和 WTSSMMPP 进行 Taylor 线性化；6–14 岁年龄域共 2,541 条记录，其中每种抗原有 1,592 条正权重且检测解释有效的记录。处理脚本为 `scripts/04_validation/validate_external_surveillance.py`，结果位于 `output/revised/external_validation/`。

为独立年代和检测平台敏感性，已保存 [CDC NHANES 2009–2010 MMRV_F](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/MMRV_F.htm) 与 DEMO_F 文件至 `data/external/nhanes_2009_2010/`。该周期使用 Wampole ELISA，2017–March 2020 使用 multiplex bead assay；两者不会直接合并为同一血清学估计。

补充的历史聚合血清学表格包括 WHO ESEN2 的麻疹、腮腺炎和风疹国家–年龄带结果，以及德国 KiGGS MMR 年龄带结果；原始公开表格、转录 CSV 和方法限制保存在 `data/external/esen2_2000s/` 与 `data/external/kiggs_2003_2006/`。这些记录用于年龄梯度背景和敏感性展示，保持与 NHANES 个体层估计和程序性保护缺口分开。
