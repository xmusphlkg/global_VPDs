"""
从WHO原始数据中提取所有疾病的病例数据
包括所有13种疾病
"""

import pandas as pd
import numpy as np

print("正在读取WHO报告病例数据...")
cases_raw = pd.read_excel('data/raw/who_immunization_portal_reported_cases.xlsx', sheet_name='Data')

print(f"原始数据: {len(cases_raw)}行")
print(f"列名: {list(cases_raw.columns)}")

# 查看疾病列表
diseases = cases_raw['DISEASE'].dropna().unique()
print(f"\n可用疾病数: {len(diseases)}")
for disease in sorted([str(d) for d in diseases]):
    count = len(cases_raw[cases_raw['DISEASE'] == disease])
    print(f"  {disease}: {count}条记录")

# 提取需要的列并重命名
# 根据实际列名调整
print("\n正在处理数据...")

# 标准化列名（先查看实际列名）
print(f"\n实际列名: {cases_raw.columns.tolist()}")

# 选择相关列 (使用实际存在的列名)
cols_to_keep = ['CODE', 'YEAR', 'DISEASE', 'CASES']
available_cols = [col for col in cols_to_keep if col in cases_raw.columns]

if len(available_cols) == len(cols_to_keep):
    df_cases = cases_raw[cols_to_keep].copy()
    df_cases.columns = ['country', 'year', 'disease', 'cases']

    # 清理数据
    df_cases = df_cases.dropna(subset=['country', 'year', 'disease'])
    df_cases['year'] = df_cases['year'].astype(int)
    df_cases['cases'] = pd.to_numeric(df_cases['cases'], errors='coerce').fillna(0)

    # 排序
    df_cases = df_cases.sort_values(['country', 'disease', 'year'])

    # 保存
    output_file = 'data/raw/case_data.csv'
    df_cases.to_csv(output_file, index=False)

    print(f"\n✅ 病例数据已更新")
    print(f"   文件: {output_file}")
    print(f"   总记录数: {len(df_cases)}")
    print(f"   国家数: {df_cases['country'].nunique()}")
    print(f"   疾病数: {df_cases['disease'].nunique()}")
    print(f"   年份范围: {df_cases['year'].min()} - {df_cases['year'].max()}")

    print("\n各疾病记录数:")
    disease_counts = df_cases.groupby('disease').size().sort_values(ascending=False)
    for disease, count in disease_counts.items():
        print(f"  {disease}: {count:,}")

    print("\n数据预览:")
    print(df_cases.head(10))

else:
    print(f"❌ 列名不匹配。需要: {cols_to_keep}, 可用: {available_cols}")
    print("请检查原始数据的列名")
