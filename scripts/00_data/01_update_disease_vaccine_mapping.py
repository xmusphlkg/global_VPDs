"""
更新疾病-疫苗映射关系
将所有WHO可用的VPD及其对应疫苗整理成映射表
"""

import pandas as pd
import os

# 定义完整的疾病-疫苗映射
# 基于WHO数据中的实际可用数据
mapping_data = [
    # 核心病种 - 已有
    {'disease': 'Measles', 'disease_code': 'MEASLES', 'vaccine': 'MCV1', 'vaccine_code': 'MCV1', 'priority': 1, 'note': '核心病种'},
    {'disease': 'Measles', 'disease_code': 'MEASLES', 'vaccine': 'MCV2', 'vaccine_code': 'MCV2', 'priority': 1, 'note': '核心病种'},
    {'disease': 'Rubella', 'disease_code': 'RUBELLA', 'vaccine': 'RCV1', 'vaccine_code': 'RCV1', 'priority': 1, 'note': '核心病种'},
    {'disease': 'Pertussis', 'disease_code': 'PERTUSSIS', 'vaccine': 'DTPCV3', 'vaccine_code': 'DTPCV3', 'priority': 1, 'note': '核心病种'},
    {'disease': 'Pertussis', 'disease_code': 'PERTUSSIS', 'vaccine': 'DTPCV1', 'vaccine_code': 'DTPCV1', 'priority': 1, 'note': '核心病种'},
    {'disease': 'Diphtheria', 'disease_code': 'DIPHTHERIA', 'vaccine': 'DTPCV3', 'vaccine_code': 'DTPCV3', 'priority': 1, 'note': '核心病种'},
    {'disease': 'Diphtheria', 'disease_code': 'DIPHTHERIA', 'vaccine': 'DTPCV1', 'vaccine_code': 'DTPCV1', 'priority': 1, 'note': '核心病种'},

    # 立即补充 - 优先级1
    {'disease': 'Mumps', 'disease_code': 'MUMPS', 'vaccine': 'MCV1', 'vaccine_code': 'MCV1', 'priority': 1, 'note': '需核实MMR vs MR'},
    {'disease': 'Mumps', 'disease_code': 'MUMPS', 'vaccine': 'MCV2', 'vaccine_code': 'MCV2', 'priority': 1, 'note': '需核实MMR vs MR'},
    {'disease': 'Congenital Rubella Syndrome', 'disease_code': 'CRS', 'vaccine': 'RCV1', 'vaccine_code': 'RCV1', 'priority': 1, 'note': '风疹核心结局'},

    # 破伤风 - 需单独设计
    {'disease': 'Neonatal Tetanus', 'disease_code': 'NTETANUS', 'vaccine': 'TT2+', 'vaccine_code': 'TT2PLUS', 'priority': 2, 'note': '孕妇接种，母体免疫'},
    {'disease': 'Neonatal Tetanus', 'disease_code': 'NTETANUS', 'vaccine': 'PAB', 'vaccine_code': 'PAB', 'priority': 2, 'note': '分娩时保护'},
    {'disease': 'Total Tetanus', 'disease_code': 'TTETANUS', 'vaccine': 'DTPCV3', 'vaccine_code': 'DTPCV3', 'priority': 2, 'note': '婴儿接种'},
    {'disease': 'Total Tetanus', 'disease_code': 'TTETANUS', 'vaccine': 'TT2+', 'vaccine_code': 'TT2PLUS', 'priority': 2, 'note': '其他人群'},

    # 脊髓灰质炎 - 专门分析
    {'disease': 'Poliomyelitis', 'disease_code': 'POLIO', 'vaccine': 'IPV1', 'vaccine_code': 'IPV1', 'priority': 2, 'note': '灭活疫苗'},
    {'disease': 'Poliomyelitis', 'disease_code': 'POLIO', 'vaccine': 'IPVC', 'vaccine_code': 'IPVC', 'priority': 2, 'note': '灭活疫苗完整程序'},

    # 地区性疾病 - 需审计
    {'disease': 'Japanese Encephalitis', 'disease_code': 'JAPENC', 'vaccine': 'JAPENC', 'vaccine_code': 'JAPENC', 'priority': 3, 'note': '亚洲地区'},
    {'disease': 'Yellow Fever', 'disease_code': 'YFEVER', 'vaccine': 'YFV', 'vaccine_code': 'YFV', 'priority': 3, 'note': '非洲/南美地区'},

    # 数据较少的疾病 - 需评估
    {'disease': 'Invasive Meningococcal Disease', 'disease_code': 'INVASIVE_MENING', 'vaccine': 'MenA-Conj', 'vaccine_code': 'MEN_A_CONJ', 'priority': 3, 'note': '脑膜炎带（非洲）'},
    {'disease': 'Invasive Meningococcal Disease', 'disease_code': 'INVASIVE_MENING', 'vaccine': 'MenACYW-Conj', 'vaccine_code': 'MEN_ACYW_CONJ', 'priority': 3, 'note': '多价结合疫苗'},
    {'disease': 'Invasive Meningococcal Disease', 'disease_code': 'INVASIVE_MENING', 'vaccine': 'MenB', 'vaccine_code': 'MEN_B', 'priority': 3, 'note': 'B群疫苗'},
    {'disease': 'Typhoid', 'disease_code': 'TYPHOID', 'vaccine': 'Typhoid', 'vaccine_code': 'TYPHOID', 'priority': 3, 'note': '多糖疫苗'},
    {'disease': 'Typhoid', 'disease_code': 'TYPHOID', 'vaccine': 'Typhoid-Conj', 'vaccine_code': 'TYPHOID_CONJ', 'priority': 3, 'note': '结合疫苗'},
]

# 创建DataFrame
df_mapping = pd.DataFrame(mapping_data)

# 保存
output_file = 'data/processed/vaccine_disease_mapping.csv'
df_mapping.to_csv(output_file, index=False)

print(f"✅ 疾病-疫苗映射已更新")
print(f"   文件: {output_file}")
print(f"   疾病数: {df_mapping['disease'].nunique()}")
print(f"   疾病代码数: {df_mapping['disease_code'].nunique()}")
print(f"   映射关系数: {len(df_mapping)}")
print("\n各优先级疾病数：")
print(df_mapping.groupby('priority')['disease'].nunique())
