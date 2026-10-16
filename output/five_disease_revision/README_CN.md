# 五病队列重建（revision）

观察年份为 2014–2023，每个国家—年份保留 1–14 岁完整队列。剂次覆盖率按国家和计划接种年份显式合并。

- 计划年份晚于观察年份：`not_due`，覆盖率为 0。
- 覆盖率缺失或腮腺炎成分无法确认：`unknown`，不填充为 0。
- 两剂覆盖率违反嵌套关系时，`strict` 行标为不一致并排除严格主样本；`projected` 行把第二剂投影至第一剂以内。
- 主样本要求 14 个年龄行、严格剂次状态全部已知；敏感性分析要求至少 80% 年龄行已知。
- 腮腺炎严格成分证据限于 `data/curated/mumps_program_evidence.csv` 的 USA、GBR、FIN、CHN；历史程序来源和出生队列边界均保存在该表。全 MCV 结果另列 `all_mcv_proxy`。
- 腮腺炎初始应答为0.964，指数持续时间均值27.4年；16/27/51年为舍入敏感性网格。白喉部分/完整临床保护指数0.71/0.87，30年为参考衰减情景，检验10/30/100年及不衰减。百日咳5年为产品混合下的指数情景；所有定义见补充参数表。
- 全球疾病程序使用固定标准化接种年龄；腮腺炎严格子集使用证据表中的年龄场景。个人联合剂次、补种和迁移没有观测，出生年份加整年年龄映射写入每剂 `*_target_age_months` 与 `*_year`。
- `cohort_gap_total` 与 `cohort_gap_no_waning` 兼容下游模型；效力只应用一次，时间衰减只应用于效力，年份差不再扣除年龄。

主要表：`country_year_14cohort_primary.csv`、`country_year_14cohort_80pct.csv`、`cohort_level_all_ages.csv`。输入文件和定义见 `reconstruction_manifest.json`。
